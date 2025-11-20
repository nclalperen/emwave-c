// ============================================================================
// emwave-c: Dear ImGui Front-End with Simulation Wizard
// ============================================================================

#include "config.h"
#include "types.h"
#include "analysis.h"
#include "app_bootstrap.h"
#include "fdtd_core.h"
#include "ui_render.h"
#include "materials.h"
#include "boundary.h"
#include "sources.h"
#include "expr.h"
#include "core/material_library.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"
#include "implot.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

struct AppState {
    bool basic_mode;
    int selected_source;
    int selected_block;
    bool placing_source;
    bool placing_block;
    bool block_first_set;
    int block_first_i;
    int block_first_j;
    int last_click_i;
    int last_click_j;
    bool show_scope_window;
    bool show_scene_panel;
    bool show_sources_panel;
    bool show_blocks_panel;
    bool show_run_panel;
    bool show_grid_panel;
    bool show_run_settings_panel;
    bool show_probes_panel;
    bool show_log_panel;
    bool show_expression_panel;
    bool show_scenes_panel;
    bool show_grid_overlay;
    /* Simple log buffer */
    char log_lines[128][256];
    int log_count;
    ImVec2 viewport_pos;
    ImVec2 viewport_size;
    bool viewport_valid;
    /* Material browser */
    int selected_material_id;
    int paint_material_id;
    bool material_browser_open;
    char material_search[64];
    bool filter_metals;
    bool filter_dielectrics;
    /* Material visualization */
    bool show_material_legend;
    bool filter_blocks_by_material;
    bool auto_filter_blocks_on_select;
    bool highlight_blocks_by_material;
    int visualization_mode;   // 0=Field,1=Material,2=Overlay
    bool show_material_outlines;
    float material_overlay_alpha;
    std::vector<double> sparam_freq_hz;
    std::vector<double> sparam_freq_ghz;
    std::vector<double> sparam_s21_mag;
    std::vector<double> sparam_s21_db;
    std::vector<double> sparam_vswr;
    std::vector<double> sparam_s11_mag;
    std::vector<double> sparam_s11_db;
    bool sparam_window_open;
    bool sparam_data_loaded;
    double sparam_peak_db;
    double sparam_peak_freq_hz;
    double sparam_f_low_hz;
    double sparam_f_high_hz;
    double sparam_f_center_hz;
    double sparam_bandwidth_hz;
    char sparam_csv_path[260];
    bool smith_chart_open;
    bool smith_show_impedance;
    bool smith_show_vswr_circles;
    double smith_z0;
    int smith_freq_index;
    bool request_rebootstrap;
    char rebootstrap_message[128];
};

struct WizardState {
    SimulationConfig cfg;
    bool open;
    bool advanced;
};

static void apply_wizard_materials_to_sim(const WizardState& wizard,
                                          SimulationBootstrap* bootstrap,
                                          SimulationState* sim);

static void wizard_init_from_config(WizardState& w, const SimulationConfig* cfg) {
    if (cfg) {
        w.cfg = *cfg;
    } else {
        w.cfg = SIM_CONFIG_DEFAULTS;
    }
    w.open = true;
    w.advanced = false;
}

static void ui_log_add(AppState* app, const char* fmt, ...) {
    if (!app) return;
    if (app->log_count >= 128) {
        app->log_count = 0;
    }
    va_list args;
    va_start(args, fmt);
    vsnprintf(app->log_lines[app->log_count], sizeof(app->log_lines[app->log_count]), fmt, args);
    va_end(args);
    app->log_count++;
}

static const char* source_type_label(SourceType t) {
    switch (t) {
        case SRC_CW:           return "CW";
        case SRC_GAUSS_PULSE:  return "Gaussian";
        case SRC_RICKER:       return "Ricker";
        case SRC_EXPR:         return "Expr";
        default:               return "Unknown";
    }
}

static const char* source_field_label(SourceFieldType f) {
    switch (f) {
        case SRC_FIELD_EZ: return "Ez";
        case SRC_FIELD_HX: return "Hx";
        case SRC_FIELD_HY: return "Hy";
        default:           return "Ez";
    }
}

static double clamp01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

static int normalized_to_cell_index(double frac, int n) {
    if (n <= 0) return 0;
    double clamped = clamp01(frac);
    int idx = (int)std::lround(clamped * (double)(n - 1));
    const int pad = 2;
    if (idx < pad) idx = pad;
    if (idx >= n - pad && n > pad) idx = n - pad - 1;
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    return idx;
}

static void free_source_expression(Source* s) {
    if (!s || !s->expr_program) return;
    expr_free((ExprProgram*)s->expr_program);
    s->expr_program = NULL;
}

static void apply_source_spec_to_runtime(SimulationState* sim,
                                         int idx,
                                         const SourceConfigSpec* spec) {
    if (!sim || !spec) return;
    if (idx < 0 || idx >= MAX_SRC) return;
    Source* dst = &sim->sources[idx];
    dst->active = spec->active;
    dst->type = spec->type;
    dst->field = spec->field;
    dst->amp = spec->amp;
    dst->freq = spec->freq;
    dst->sigma2 = (spec->sigma2 > 0.01) ? spec->sigma2 : 0.01;
    dst->ix = normalized_to_cell_index(spec->x, sim->nx);
    dst->iy = normalized_to_cell_index(spec->y, sim->ny);
    std::strncpy(dst->expr_text, spec->expr, SOURCE_EXPR_MAX_LEN);
    dst->expr_text[SOURCE_EXPR_MAX_LEN - 1] = '\0';
    free_source_expression(dst);
    source_reparam(dst);
    if (dst->type == SRC_EXPR && dst->expr_text[0] != '\0') {
        ExprProgram* prog = NULL;
        char errbuf[128];
        if (expr_compile(dst->expr_text, &prog, errbuf, (int)sizeof(errbuf))) {
            dst->expr_program = prog;
        } else {
            dst->expr_program = NULL;
        }
    }
}

static void clear_runtime_source(SimulationState* sim, int idx) {
    if (!sim || idx < 0 || idx >= MAX_SRC) return;
    Source* dst = &sim->sources[idx];
    dst->active = 0;
    free_source_expression(dst);
    dst->expr_text[0] = '\0';
}

static void create_new_source(WizardState* wizard, SimulationState* sim, AppState* app) {
    if (!wizard || !sim) return;
    int count = wizard->cfg.source_count;
    if (count >= MAX_SRC) return;
    SourceConfigSpec* spec = &wizard->cfg.source_configs[count];
    spec->active = 1;
    spec->type = SRC_CW;
    spec->field = SRC_FIELD_EZ;
    spec->x = 0.5;
    spec->y = 0.5;
    spec->amp = 1.0;
    spec->freq = 1.0e9;
    spec->sigma2 = 4.0;
    spec->expr[0] = '\0';
    wizard->cfg.source_count = count + 1;
    apply_source_spec_to_runtime(sim, count, spec);
    if (app) {
        app->selected_source = count;
        ui_log_add(app, "Added source #%d", count);
    }
}

static const Material* guess_material_for_rect(const MaterialRectSpec& r) {
    const Material* best = nullptr;
    double best_score = 1e9;
    int mat_count = material_library_get_count();
    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;
        if (r.tag == 1) {
            if (mat->type == MAT_TYPE_PEC) return mat;
            continue;
        }
        if (r.tag == 2) {
            if (mat->type == MAT_TYPE_PMC) return mat;
            continue;
        }
        if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) continue;
        double diff = std::fabs(mat->epsilon_r - r.epsr);
        if (diff < best_score) {
            best_score = diff;
            best = mat;
        }
    }
    return best;
}

static void apply_material_to_rect(MaterialRectSpec* rect, const Material* mat) {
    if (!rect || !mat) return;
    if (mat->type == MAT_TYPE_PEC) {
        rect->tag = 1;
        rect->epsr = 1.0;
        rect->sigma = mat->conductivity;
    } else if (mat->type == MAT_TYPE_PMC) {
        rect->tag = 2;
        rect->epsr = 1.0;
        rect->sigma = 0.0;
    } else {
        rect->tag = 0;
        rect->epsr = mat->epsilon_r;
        rect->sigma = mat->conductivity;
    }
}

static void remove_block(WizardState* wizard,
                         SimulationBootstrap* bootstrap,
                         SimulationState* sim,
                         AppState* app,
                         int idx) {
    if (!wizard || idx < 0 || idx >= wizard->cfg.material_rect_count) return;
    for (int j = idx; j < wizard->cfg.material_rect_count - 1; ++j) {
        wizard->cfg.material_rects[j] = wizard->cfg.material_rects[j + 1];
    }
    wizard->cfg.material_rect_count--;
    if (bootstrap && sim) {
        apply_wizard_materials_to_sim(*wizard, bootstrap, sim);
    }
    if (app) {
        ui_log_add(app, "Deleted block #%d", idx);
    }
}

static void delete_source_at(WizardState* wizard, SimulationState* sim, AppState* app, int idx) {
    if (!wizard || !sim) return;
    int count = wizard->cfg.source_count;
    if (idx < 0 || idx >= count) return;
    for (int i = idx; i < count - 1; ++i) {
        wizard->cfg.source_configs[i] = wizard->cfg.source_configs[i + 1];
        apply_source_spec_to_runtime(sim, i, &wizard->cfg.source_configs[i]);
    }
    wizard->cfg.source_count = count - 1;
    clear_runtime_source(sim, wizard->cfg.source_count);
    if (app) {
        if (app->selected_source >= wizard->cfg.source_count) {
            app->selected_source = wizard->cfg.source_count - 1;
        }
        ui_log_add(app, "Deleted source #%d", idx);
    }
}

/* Scene overview panel ---------------------------------------------------- */
static void draw_scene_panel(const SimulationState* sim, const WizardState& wizard) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Scene")) return;

    ImGui::Indent();
    if (sim) {
        ImGui::Text("Grid: %d x %d", sim->nx, sim->ny);
        ImGui::Text("Domain: %.3f m x %.3f m", sim->lx, sim->ly);
        ImGui::Text("dt: %.3e s", sim->dt);
        ImGui::Separator();
    } else {
        ImGui::TextUnformatted("Simulation not available");
        ImGui::Unindent();
        return;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Materials");
    ImGui::Indent();
    ImGui::Text("Rectangles: %d", wizard.cfg.material_rect_count);
    ImGui::Unindent();

    ImGui::Separator();
    ImGui::TextUnformatted("Sources");
    ImGui::Indent();
    ImGui::Text("Count: %d (max %d)", wizard.cfg.source_count, MAX_SRC);
    ImGui::Unindent();
}

static void draw_grid_panel(WizardState* wizard, AppState* app) {
    if (!wizard) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Grid & Domain")) return;

    ImGui::Indent();
    int nx = wizard->cfg.nx;
    if (ImGui::InputInt("Cells X (nx)", &nx)) {
        if (nx < 8) nx = 8;
        wizard->cfg.nx = nx;
    }
    int ny = wizard->cfg.ny;
    if (ImGui::InputInt("Cells Y (ny)", &ny)) {
        if (ny < 8) ny = 8;
        wizard->cfg.ny = ny;
    }

    double lx = wizard->cfg.lx;
    if (ImGui::InputDouble("Width (meters)", &lx, 0.01, 0.1, "%.3f")) {
        if (lx < 1e-3) lx = 1e-3;
        wizard->cfg.lx = lx;
    }
    double ly = wizard->cfg.ly;
    if (ImGui::InputDouble("Height (meters)", &ly, 0.01, 0.1, "%.3f")) {
        if (ly < 1e-3) ly = 1e-3;
        wizard->cfg.ly = ly;
    }

    float cfl = (float)wizard->cfg.cfl_safety;
    if (ImGui::SliderFloat("CFL Safety", &cfl, 0.1f, 0.99f, "%.2f")) {
        wizard->cfg.cfl_safety = (double)cfl;
    }

    const char* boundary_items[] = {"CPML", "Mur"};
    int boundary_idx = (wizard->cfg.boundary_mode == SIM_BOUNDARY_MUR) ? 1 : 0;
    if (ImGui::Combo("Boundary", &boundary_idx, boundary_items, IM_ARRAYSIZE(boundary_items))) {
        wizard->cfg.boundary_mode = (boundary_idx == 1) ? SIM_BOUNDARY_MUR : SIM_BOUNDARY_CPML;
    }

    ImGui::Separator();
    if (ImGui::Button("Apply Grid Changes & Restart")) {
        if (app) {
            app->request_rebootstrap = true;
            std::snprintf(app->rebootstrap_message,
                          sizeof(app->rebootstrap_message),
                          "Rebooted after grid/domain update");
        }
    }
    ImGui::Unindent();
}

/* Sources panel ----------------------------------------------------------- */
static void draw_sources_panel(SimulationState* sim, WizardState& wizard, AppState* app) {
    if (!sim || !app) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Sources")) return;

    ImGui::Indent();
    int count = wizard.cfg.source_count;
    if (count < 0) count = 0;
    if (count > MAX_SRC) count = MAX_SRC;
    wizard.cfg.source_count = count;

    if (count == 0) {
        app->selected_source = -1;
    } else {
        if (app->selected_source < 0) app->selected_source = 0;
        if (app->selected_source >= count) app->selected_source = count - 1;
    }

    ImGui::Text("Configured: %d / %d", count, MAX_SRC);
    ImGui::SameLine();
    bool can_add = (count < MAX_SRC);
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ New Source")) {
        create_new_source(&wizard, sim, app);
        count = wizard.cfg.source_count;
    }
    if (!can_add) {
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::TextDisabled("Max reached");
    }
    ImGui::Separator();

    if (count == 0) {
        ImGui::TextDisabled("No sources configured.");
    } else {
        ImGui::TextUnformatted("Source list");
        ImGui::Indent();
        for (int i = 0; i < count; ++i) {
            SourceConfigSpec& spec = wizard.cfg.source_configs[i];
            ImGui::PushID(i);
            char label[64];
            std::snprintf(label,
                          sizeof(label),
                          "#%d %s",
                          i,
                          source_type_label(spec.type));
            bool selected = (app->selected_source == i);
            if (ImGui::Selectable(label, selected)) {
                app->selected_source = i;
            }
            ImGui::SameLine();
            bool active = spec.active != 0;
            if (ImGui::Checkbox("##active", &active)) {
                spec.active = active ? 1 : 0;
                apply_source_spec_to_runtime(sim, i, &spec);
            }
            ImGui::SameLine();
            ImGui::Text("%s", source_field_label(spec.field));
            ImGui::SameLine();
            ImGui::TextDisabled("Pos: (%.2f, %.2f)", spec.x, spec.y);
            ImGui::PopID();
        }
        ImGui::Unindent();
    }

    int idx = app->selected_source;
    if (idx >= 0 && idx < count) {
        SourceConfigSpec& spec = wizard.cfg.source_configs[idx];
        ImGui::Separator();
        ImGui::Text("Source #%d", idx);
        ImGui::Indent();

        const char* types[] = {"CW", "Gaussian", "Ricker", "Expr"};
        int type_idx = (int)spec.type;
        if (type_idx < 0 || type_idx >= (int)IM_ARRAYSIZE(types)) type_idx = 0;
        if (ImGui::Combo("Type", &type_idx, types, IM_ARRAYSIZE(types))) {
            spec.type = (SourceType)type_idx;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }

        const char* field_items[] = {"Ez", "Hx", "Hy"};
        int field_idx = (int)spec.field;
        if (field_idx < 0 || field_idx >= (int)IM_ARRAYSIZE(field_items)) field_idx = 0;
        if (ImGui::Combo("Field", &field_idx, field_items, IM_ARRAYSIZE(field_items))) {
            spec.field = (SourceFieldType)field_idx;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }

        float norm_x = (float)spec.x;
        if (ImGui::SliderFloat("X (0-1)", &norm_x, 0.0f, 1.0f, "%.3f")) {
            spec.x = clamp01(norm_x);
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        float norm_y = (float)spec.y;
        if (ImGui::SliderFloat("Y (0-1)", &norm_y, 0.0f, 1.0f, "%.3f")) {
            spec.y = clamp01(norm_y);
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        const Source& runtime_src = sim->sources[idx];
        ImGui::TextDisabled("Grid cell: (%d, %d)", runtime_src.ix, runtime_src.iy);

        double amp = spec.amp;
        if (ImGui::InputDouble("Amplitude", &amp, 0.1, 1.0, "%.3f")) {
            spec.amp = amp;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        double freq = spec.freq;
        if (ImGui::InputDouble("Frequency (Hz)", &freq, 1e6, 1e8, "%.3e")) {
            if (freq > 0.0) {
                spec.freq = freq;
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
        }

        if (spec.type == SRC_GAUSS_PULSE || spec.type == SRC_RICKER) {
            double sigma = spec.sigma2;
            if (ImGui::InputDouble("Sigma^2 (cells^2)", &sigma, 0.1, 1.0, "%.3f")) {
                if (sigma < 0.1) sigma = 0.1;
                spec.sigma2 = sigma;
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
        } else {
            ImGui::TextDisabled("Sigma^2 applies to Gaussian/Ricker pulses.");
        }

        if (spec.type == SRC_EXPR) {
            if (ImGui::InputTextMultiline("Expression",
                                          spec.expr,
                                          SOURCE_EXPR_MAX_LEN,
                                          ImVec2(-1.0f, ImGui::GetTextLineHeight() * 4))) {
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
            ImGui::TextDisabled("Variables: t (s), amp, freq, pi");

            if (spec.expr[0] != '\0') {
                char errbuf[128];
                ExprProgram* prog = nullptr;
                if (expr_compile(spec.expr, &prog, errbuf, sizeof(errbuf))) {
                    const int N = 128;
                    float values[N];
                    double period = (spec.freq > 0.0) ? (1.0 / spec.freq) : 1e-9;
                    double t_max = (period > 0.0) ? 3.0 * period : 3e-9;
                    double vmin = 0.0;
                    double vmax = 0.0;
                    bool first = true;
                    for (int k = 0; k < N; ++k) {
                        double t = t_max * (double)k / (double)(N - 1);
                        double v = expr_eval(prog, t, spec.amp, spec.freq);
                        values[k] = (float)v;
                        if (first) {
                            vmin = vmax = v;
                            first = false;
                        } else {
                            if (v < vmin) vmin = v;
                            if (v > vmax) vmax = v;
                        }
                    }
                    float ymin = (float)vmin;
                    float ymax = (float)vmax;
                    if (ymin == ymax) {
                        ymin -= 1.0f;
                        ymax += 1.0f;
                    }
                    ImGui::PlotLines("Preview", values, N, 0, nullptr, ymin, ymax,
                                     ImVec2(-1.0f, ImGui::GetTextLineHeight() * 5));
                    expr_free(prog);
                } else {
                    ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.5f, 1.0f),
                                       "Expression error: %s",
                                       errbuf);
                }
            }
        }

        ImGui::Separator();
        if (ImGui::Button("Delete Source")) {
            ImGui::OpenPopup("delete_source_popup");
        }
        if (ImGui::BeginPopupModal("delete_source_popup", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("Delete source #%d?", idx);
            if (ImGui::Button("Delete", ImVec2(120, 0))) {
                delete_source_at(&wizard, sim, app, idx);
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
        ImGui::Unindent();
    }
    ImGui::Unindent();
}

static void draw_run_settings_panel(WizardState* wizard, AppState* app) {
    if (!wizard) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Run Settings")) return;

    ImGui::Indent();
    const char* run_modes[] = {"Fixed steps", "Sweep"};
    int run_mode_idx = (wizard->cfg.run_mode == SIM_RUN_MODE_SWEEP) ? 1 : 0;
    if (ImGui::Combo("Run Mode", &run_mode_idx, run_modes, IM_ARRAYSIZE(run_modes))) {
        wizard->cfg.run_mode = (run_mode_idx == 1) ? SIM_RUN_MODE_SWEEP : SIM_RUN_MODE_FIXED_STEPS;
    }

    if (wizard->cfg.run_mode == SIM_RUN_MODE_FIXED_STEPS) {
        ImGui::InputInt("Run steps", &wizard->cfg.run_steps);
        if (wizard->cfg.run_steps < 1) wizard->cfg.run_steps = 1;
    } else {
        ImGui::InputInt("Sweep points", &wizard->cfg.sweep_points);
        if (wizard->cfg.sweep_points < 1) wizard->cfg.sweep_points = 1;
        ImGui::InputDouble("Sweep start (Hz)", &wizard->cfg.sweep_start_hz, 1e6, 1e8, "%.3e");
        if (wizard->cfg.sweep_start_hz < 1.0) wizard->cfg.sweep_start_hz = 1.0;
        ImGui::InputDouble("Sweep stop (Hz)", &wizard->cfg.sweep_stop_hz, 1e6, 1e8, "%.3e");
        if (wizard->cfg.sweep_stop_hz < wizard->cfg.sweep_start_hz) {
            wizard->cfg.sweep_stop_hz = wizard->cfg.sweep_start_hz;
        }
        ImGui::InputInt("Steps / point", &wizard->cfg.sweep_steps_per_point);
        if (wizard->cfg.sweep_steps_per_point < 1) wizard->cfg.sweep_steps_per_point = 1;
    }

    if (wizard->advanced) {
        ImGui::Separator();
        bool profile = wizard->cfg.enable_profile != 0;
        if (ImGui::Checkbox("Enable profiling", &profile)) {
            wizard->cfg.enable_profile = profile ? 1 : 0;
        }
        bool probe_log = wizard->cfg.enable_probe_log != 0;
        if (ImGui::Checkbox("Enable probe log", &probe_log)) {
            wizard->cfg.enable_probe_log = probe_log ? 1 : 0;
        }
        ImGui::InputText("Probe log path",
                         wizard->cfg.probe_log_path,
                         SIM_PROBE_LOG_PATH_MAX);
    } else {
        ImGui::Separator();
        ImGui::TextDisabled("Switch to Advanced mode for profiling/logging options.");
    }

    ImGui::Separator();
    if (ImGui::Button("Apply Run Settings & Restart")) {
        if (app) {
            app->request_rebootstrap = true;
            std::snprintf(app->rebootstrap_message,
                          sizeof(app->rebootstrap_message),
                          "Rebooted after run settings update");
        }
    }
    ImGui::Unindent();
}

/* Blocks panel ------------------------------------------------------------ */
static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app);

typedef struct {
    int material_id;
    const char* material_name;
    int cell_count;
    float percentage;
} MaterialStats;

static int get_material_statistics(const SimulationState* sim, MaterialStats* stats, int max_stats) {
    (void)sim;
    (void)stats;
    (void)max_stats;
    /* Placeholder until simulation tracks material IDs explicitly. */
    return 0;
}

static int compute_scope_fft(const Scope* scope,
                             double dt,
                             double* freq,
                             double* mag,
                             double* phase,
                             int max_fft) {
    if (!scope || !scope->y || scope->n <= 8 || dt <= 0.0) return 0;
    if (!freq || !mag || !phase || max_fft <= 0) return 0;

    int Nfft = scope->n;
    if (Nfft > max_fft) Nfft = max_fft;
    if (Nfft < 64) return 0;

    std::vector<double> x(Nfft);
    int start = (scope->head - Nfft + scope->n) % scope->n;
    for (int k = 0; k < Nfft; ++k) {
        int idx = (start + k) % scope->n;
        x[k] = scope->y[idx];
    }

    double mean = 0.0;
    for (int k = 0; k < Nfft; ++k) mean += x[k];
    mean /= (double)Nfft;

    const double two_pi = 2.0 * 3.14159265358979323846;
    if (Nfft > 1) {
        for (int k = 0; k < Nfft; ++k) {
            double w = 0.5 * (1.0 - std::cos(two_pi * (double)k / (double)(Nfft - 1)));
            x[k] = (x[k] - mean) * w;
        }
    } else {
        for (int k = 0; k < Nfft; ++k) {
            x[k] = (x[k] - mean);
        }
    }

    int half = Nfft / 2;
    if (half <= 1) return 0;

    for (int k = 0; k < half; ++k) {
        double real = 0.0;
        double imag = 0.0;
        for (int n = 0; n < Nfft; ++n) {
            double angle = -two_pi * (double)k * (double)n / (double)Nfft;
            double c = std::cos(angle);
            double s = std::sin(angle);
            real += x[n] * c;
            imag += x[n] * s;
        }
        real /= (double)Nfft;
        imag /= (double)Nfft;
        mag[k] = std::sqrt(real * real + imag * imag);
        phase[k] = std::atan2(imag, real) * (180.0 / 3.14159265358979323846);
        freq[k] = (double)k / ((double)Nfft * dt);
    }

    return half;
}

static void clear_sparameter_data(AppState* app) {
    if (!app) return;
    app->sparam_freq_hz.clear();
    app->sparam_freq_ghz.clear();
    app->sparam_s21_mag.clear();
    app->sparam_s21_db.clear();
    app->sparam_s11_mag.clear();
    app->sparam_s11_db.clear();
    app->sparam_vswr.clear();
    app->sparam_peak_db = 0.0;
    app->sparam_peak_freq_hz = 0.0;
    app->sparam_f_low_hz = 0.0;
    app->sparam_f_high_hz = 0.0;
    app->sparam_f_center_hz = 0.0;
    app->sparam_bandwidth_hz = 0.0;
    app->sparam_data_loaded = false;
}

static void analyze_sparameter_bandwidth(AppState* app) {
    if (!app || app->sparam_freq_hz.empty()) return;
    app->sparam_peak_db = -1e9;
    app->sparam_peak_freq_hz = 0.0;
    const size_t count = app->sparam_freq_hz.size();
    for (size_t i = 0; i < count; ++i) {
        if (app->sparam_s21_db[i] > app->sparam_peak_db) {
            app->sparam_peak_db = app->sparam_s21_db[i];
            app->sparam_peak_freq_hz = app->sparam_freq_hz[i];
        }
    }
    double threshold = app->sparam_peak_db - 3.0;
    app->sparam_f_low_hz = 0.0;
    app->sparam_f_high_hz = 0.0;
    for (size_t i = 0; i < count; ++i) {
        if (app->sparam_s21_db[i] >= threshold) {
            if (app->sparam_f_low_hz == 0.0) {
                app->sparam_f_low_hz = app->sparam_freq_hz[i];
            }
            app->sparam_f_high_hz = app->sparam_freq_hz[i];
        }
    }
    if (app->sparam_f_high_hz > app->sparam_f_low_hz && app->sparam_f_low_hz > 0.0) {
        app->sparam_f_center_hz = 0.5 * (app->sparam_f_high_hz + app->sparam_f_low_hz);
        app->sparam_bandwidth_hz = app->sparam_f_high_hz - app->sparam_f_low_hz;
    } else {
        app->sparam_f_center_hz = 0.0;
        app->sparam_bandwidth_hz = 0.0;
    }
}

static bool load_sparameter_csv(AppState* app, const char* path) {
    if (!app || !path) return false;
    std::ifstream file(path);
    if (!file.is_open()) {
        return false;
    }
    clear_sparameter_data(app);
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line.find("freq") != std::string::npos) continue;
        std::stringstream ss(line);
        std::string freq_str;
        std::string s21_str;
        if (!std::getline(ss, freq_str, ',')) continue;
        if (!std::getline(ss, s21_str, ',')) continue;
        double freq = std::strtod(freq_str.c_str(), nullptr);
        double s21 = std::strtod(s21_str.c_str(), nullptr);
        if (!std::isfinite(freq) || !std::isfinite(s21)) continue;
        if (freq <= 0.0) continue;
        if (s21 < 0.0) s21 = 0.0;
        app->sparam_freq_hz.push_back(freq);
        app->sparam_s21_mag.push_back(s21);
    }
    if (app->sparam_freq_hz.empty()) {
        app->sparam_data_loaded = false;
        return false;
    }
    const double eps = 1e-12;
    size_t count = app->sparam_freq_hz.size();
    app->sparam_freq_ghz.resize(count);
    app->sparam_s21_db.resize(count);
    app->sparam_vswr.resize(count);
    app->sparam_s11_mag.resize(count);
    app->sparam_s11_db.resize(count);
    for (size_t i = 0; i < count; ++i) {
        double mag = app->sparam_s21_mag[i];
        double mag_clamped = std::min(std::max(mag, 0.0), 0.999999);
        app->sparam_freq_ghz[i] = app->sparam_freq_hz[i] * 1e-9;
        double db = 20.0 * std::log10(mag + eps);
        app->sparam_s21_db[i] = db;
        double s11 = std::sqrt(std::max(0.0, 1.0 - (mag_clamped * mag_clamped)));
        app->sparam_s11_mag[i] = s11;
        app->sparam_s11_db[i] = 20.0 * std::log10(s11 + eps);
        double vswr = (1.0 + s11) / std::max(1e-6, 1.0 - s11);
        app->sparam_vswr[i] = vswr;
    }
    analyze_sparameter_bandwidth(app);
    app->sparam_data_loaded = true;
    std::strncpy(app->sparam_csv_path, path, sizeof(app->sparam_csv_path));
    app->sparam_csv_path[sizeof(app->sparam_csv_path) - 1] = '\0';
    if (app->smith_freq_index >= (int)count) {
        app->smith_freq_index = (int)count - 1;
        if (app->smith_freq_index < 0) app->smith_freq_index = 0;
    }
    return true;
}

static bool string_contains_case_insensitive(const char* haystack, const char* needle) {
    if (!haystack || !needle) return false;
    if (needle[0] == '\0') return true;
    for (const char* h = haystack; *h; ++h) {
        const char* h_it = h;
        const char* n_it = needle;
        while (*h_it && *n_it) {
            char hc = (char)std::tolower((unsigned char)*h_it);
            char nc = (char)std::tolower((unsigned char)*n_it);
            if (hc != nc) break;
            ++h_it;
            ++n_it;
        }
        if (*n_it == '\0') {
            return true;
        }
    }
    return false;
}

static void draw_material_browser(AppState* app, int* active_material_id) {
    if (!app || !app->material_browser_open) return;

    ImGui::SetNextWindowSize(ImVec2(400.0f, 500.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Material Browser", &app->material_browser_open)) {
        ImGui::End();
        return;
    }

    ImGui::TextUnformatted("Search:");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputText("##MaterialSearch", app->material_search, sizeof(app->material_search));

    ImGui::Separator();
    ImGui::Checkbox("Metals", &app->filter_metals);
    ImGui::SameLine();
    ImGui::Checkbox("Dielectrics", &app->filter_dielectrics);

    ImGui::Separator();

    ImGui::BeginChild("MaterialList", ImVec2(0, -170.0f), true);

    int mat_count = material_library_get_count();
    const bool has_query = (app->material_search[0] != '\0');

    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;

        if (mat->category == MAT_CAT_METAL && !app->filter_metals) continue;
        if (mat->category == MAT_CAT_DIELECTRIC && !app->filter_dielectrics) continue;

        if (has_query && !string_contains_case_insensitive(mat->name, app->material_search)) {
            continue;
        }

        ImVec4 color = ImVec4(
            mat->color_r / 255.0f,
            mat->color_g / 255.0f,
            mat->color_b / 255.0f,
            1.0f
        );

        ImGui::PushStyleColor(ImGuiCol_Text, color);
        ImGui::Bullet();
        ImGui::PopStyleColor();
        ImGui::SameLine();

        bool is_selected = (app->selected_material_id == mat->id);
        if (ImGui::Selectable(mat->name, is_selected, ImGuiSelectableFlags_SpanAvailWidth)) {
            app->selected_material_id = mat->id;
        }

        char prop_buf[64];
        if (mat->type == MAT_TYPE_PEC) {
            std::snprintf(prop_buf, sizeof(prop_buf), "sigma=%.2e S/m", mat->conductivity);
        } else {
            std::snprintf(prop_buf, sizeof(prop_buf), "eps=%.2f", mat->epsilon_r);
        }
        ImGui::SameLine(250.0f);
        ImGui::TextDisabled("%s", prop_buf);
    }

    ImGui::EndChild();

    ImGui::Separator();
    ImGui::TextUnformatted("Selected Material:");

    if (app->selected_material_id >= 0) {
        const Material* sel = material_library_get_by_id(app->selected_material_id);
        if (sel) {
            ImGui::Indent();
            ImVec4 color = ImVec4(
                sel->color_r / 255.0f,
                sel->color_g / 255.0f,
                sel->color_b / 255.0f,
                1.0f
            );
            ImGui::PushStyleColor(ImGuiCol_Text, color);
            ImGui::TextUnformatted(sel->name);
            ImGui::PopStyleColor();

            const char* category_label = (sel->category == MAT_CAT_METAL) ? "Metal" : "Dielectric";
            ImGui::TextDisabled("Category: %s", category_label);

            if (sel->type == MAT_TYPE_PEC) {
                ImGui::Text("sigma = %.2e S/m", sel->conductivity);
            } else {
                ImGui::Text("eps = %.2f", sel->epsilon_r);
                if (sel->tan_delta > 0.0) {
                    ImGui::Text("tan(delta) = %.4f", sel->tan_delta);
                }
            }

            ImGui::TextWrapped("%s", sel->description ? sel->description : "");
            ImGui::Unindent();

            if (ImGui::Button("Use for Painting", ImVec2(-1.0f, 0.0f))) {
                if (active_material_id) {
                    *active_material_id = app->selected_material_id;
                }
                app->paint_material_id = app->selected_material_id;
                app->material_browser_open = false;
            }
        }
    } else {
        ImGui::TextDisabled("No material selected");
    }

    ImGui::End();
}

static bool rect_matches_material(const MaterialRectSpec& r, const Material* mat) {
    if (!mat) return false;
    if (r.tag == 1) return mat->type == MAT_TYPE_PEC;
    if (r.tag == 2) return mat->type == MAT_TYPE_PMC;
    if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) return false;
    return std::fabs(r.epsr - mat->epsilon_r) < 0.25;
}

static void draw_material_legend(AppState* app,
                                 SimulationState* sim,
                                 WizardState* wizard,
                                 SimulationBootstrap* bootstrap) {
    (void)sim;
    if (!app || !wizard || !bootstrap) return;
    if (!app->show_material_legend) return;

    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Material Legend")) return;

    ImGui::TextUnformatted("Materials in Scene:");
    ImGui::Separator();
    ImGui::BeginChild("LegendList", ImVec2(0.0f, 160.0f), true);

    int mat_count = material_library_get_count();
    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;

        ImVec4 color = ImVec4(mat->color_r / 255.0f,
                              mat->color_g / 255.0f,
                              mat->color_b / 255.0f,
                              1.0f);
        ImDrawList* dl = ImGui::GetWindowDrawList();
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImVec2 q = ImVec2(p.x + 16.0f, p.y + 16.0f);
        dl->AddRectFilled(p, q, ImGui::ColorConvertFloat4ToU32(color));
        dl->AddRect(p, q, IM_COL32(90, 90, 90, 255));
        ImGui::Dummy(ImVec2(16.0f, 16.0f));
        ImGui::SameLine();

        bool row_selected = (app->selected_material_id == mat->id);
        if (ImGui::Selectable(mat->name, row_selected, ImGuiSelectableFlags_SpanAvailWidth)) {
            app->selected_material_id = mat->id;
            if (app->auto_filter_blocks_on_select) {
                app->filter_blocks_by_material = true;
            }
            for (int b = 0; b < wizard->cfg.material_rect_count; ++b) {
                if (rect_matches_material(wizard->cfg.material_rects[b], mat)) {
                    app->selected_block = b;
                    break;
                }
            }
        }

        ImGui::SameLine(170.0f);
        if (mat->type == MAT_TYPE_PEC) {
            ImGui::TextDisabled("PEC");
        } else if (mat->type == MAT_TYPE_PMC) {
            ImGui::TextDisabled("PMC");
        } else {
            ImGui::TextDisabled("εᵣ=%.1f", mat->epsilon_r);
        }

        ImGui::SameLine();
        bool can_add = (wizard->cfg.material_rect_count < CONFIG_MAX_MATERIAL_RECTS);
        if (!can_add) ImGui::BeginDisabled();
        char btn_label[32];
        std::snprintf(btn_label, sizeof(btn_label), "Add##%d", mat->id);
        if (ImGui::SmallButton(btn_label)) {
            int idx = wizard->cfg.material_rect_count;
            MaterialRectSpec& r = wizard->cfg.material_rects[idx];
            r.x0 = 0.3;
            r.y0 = 0.3;
            r.x1 = 0.7;
            r.y1 = 0.7;
            apply_material_to_rect(&r, mat);
            wizard->cfg.material_rect_count = idx + 1;
            apply_wizard_materials_to_sim(*wizard, bootstrap, sim);
            app->selected_block = idx;
            ui_log_add(app, "Added block #%d from %s", idx, mat->name);
        }
        if (!can_add) ImGui::EndDisabled();
    }

    ImGui::EndChild();
    ImGui::TextDisabled("Legend uses library colors (usage stats coming soon)");
}

static void draw_resistance_circle(double r_norm) {
    if (r_norm < 0.0) return;
    double cx = r_norm / (1.0 + r_norm);
    double cr = 1.0 / (1.0 + r_norm);
    ImPlotPoint center_px = ImPlot::PlotToPixels(cx, 0.0);
    ImPlotPoint edge_px = ImPlot::PlotToPixels(cx + cr, 0.0);
    float radius = (float)std::fabs(edge_px.x - center_px.x);
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    dl->AddCircle(ImVec2((float)center_px.x, (float)center_px.y),
                  radius,
                  IM_COL32(150, 150, 150, 120),
                  96,
                  1.0f);
}

static bool solve_circle_unit_intersection(double cx,
                                           double cy,
                                           double r,
                                           ImPlotPoint* out_point) {
    double d = 0.5 * (cx * cx + cy * cy + 1.0 - r * r);
    double A = cx * cx + cy * cy;
    double B = -2.0 * d * cx;
    double C = d * d - cy * cy;
    double disc = B * B - 4.0 * A * C;
    if (disc <= 0.0) return false;
    double sqrt_disc = std::sqrt(disc);
    double x1 = (-B + sqrt_disc) / (2.0 * A);
    double x2 = (-B - sqrt_disc) / (2.0 * A);
    double chosen = (std::fabs(x1 - 1.0) > 1e-3) ? x1 : x2;
    double y = 0.0;
    if (std::fabs(cy) > 1e-6) {
        y = (d - cx * chosen) / cy;
    } else {
        double val = 1.0 - chosen * chosen;
        if (val < 0.0) val = 0.0;
        y = std::sqrt(val);
        if (cy < 0.0) y = -y;
    }
    if (out_point) {
        *out_point = ImPlotPoint(chosen, y);
    }
    return true;
}

static void draw_reactance_arc(double x_norm) {
    if (std::fabs(x_norm) < 1e-3) return;
    double cy = 1.0 / x_norm;
    double cx = 1.0;
    double radius = 1.0 / std::fabs(x_norm);
    ImPlotPoint other;
    if (!solve_circle_unit_intersection(cx, cy, radius, &other)) {
        return;
    }
    const double pi = 3.14159265358979323846;
    double theta_start = std::atan2(0.0 - cy, 1.0 - cx);
    double theta_end = std::atan2(other.y - cy, other.x - cx);
    if (x_norm > 0.0) {
        if (theta_end < theta_start) theta_end += 2.0 * pi;
    } else {
        if (theta_end > theta_start) theta_end -= 2.0 * pi;
    }
    const int segments = 96;
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImVec2 prev;
    bool has_prev = false;
    for (int i = 0; i <= segments; ++i) {
        double t = (double)i / (double)segments;
        double theta = theta_start + t * (theta_end - theta_start);
        double gx = cx + radius * std::cos(theta);
        double gy = cy + radius * std::sin(theta);
        if (gx * gx + gy * gy > 1.0025) continue;
        ImPlotPoint px = ImPlot::PlotToPixels(gx, gy);
        ImVec2 pos((float)px.x, (float)px.y);
        if (has_prev) {
            dl->AddLine(prev, pos, IM_COL32(150, 150, 150, 120), 1.0f);
        }
        prev = pos;
        has_prev = true;
    }
}

static void draw_smith_grid(const AppState* app) {
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlotPoint center = ImPlot::PlotToPixels(0.0, 0.0);
    ImPlotPoint edge = ImPlot::PlotToPixels(1.0, 0.0);
    float radius = (float)std::fabs(edge.x - center.x);
    dl->AddCircle(ImVec2((float)center.x, (float)center.y),
                  radius,
                  IM_COL32(200, 200, 200, 255),
                  128,
                  2.0f);
    const double r_values[] = {0.2, 0.5, 1.0, 2.0, 5.0};
    for (double r_norm : r_values) {
        draw_resistance_circle(r_norm);
    }
    const double x_values[] = {0.2, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0};
    for (double x_norm : x_values) {
        draw_reactance_arc(x_norm);
    }
    if (app->smith_show_vswr_circles) {
        const double vswr_values[] = {2.0, 3.0, 5.0};
        for (double vs : vswr_values) {
            double gamma = (vs - 1.0) / (vs + 1.0);
            float circle_radius = radius * (float)gamma;
            dl->AddCircle(ImVec2((float)center.x, (float)center.y),
                          circle_radius,
                          IM_COL32(255, 200, 0, 90),
                          96,
                          1.0f);
        }
    }
    dl->AddLine(ImPlot::PlotToPixels(-1.0, 0.0),
                ImPlot::PlotToPixels(1.0, 0.0),
                IM_COL32(255, 255, 255, 80),
                1.0f);
    dl->AddLine(ImPlot::PlotToPixels(0.0, -1.0),
                ImPlot::PlotToPixels(0.0, 1.0),
                IM_COL32(255, 255, 255, 80),
                1.0f);
}

static void plot_smith_reflection_data(AppState* app) {
    const int max_samples = 1024;
    static double gamma_r[max_samples];
    static double gamma_i[max_samples];
    int count = (int)app->sparam_s11_mag.size();
    if (count <= 0) return;
    if (count > max_samples) count = max_samples;
    for (int i = 0; i < count; ++i) {
        double mag = app->sparam_s11_mag[i];
        gamma_r[i] = mag;
        gamma_i[i] = 0.0;
    }
    ImPlot::PlotLine("S11", gamma_r, gamma_i, count);
    if (app->smith_freq_index >= 0 && app->smith_freq_index < count) {
        int idx = app->smith_freq_index;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0, 0, 1));
        ImPlot::PlotScatter("Current", &gamma_r[idx], &gamma_i[idx], 1);
        ImPlot::PopStyleColor();
    }
}

static void draw_smith_chart(AppState* app) {
    if (!app || !app->smith_chart_open) return;
    ImGui::SetNextWindowSize(ImVec2(620.0f, 680.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Smith Chart", &app->smith_chart_open)) {
        ImGui::End();
        return;
    }

    ImGui::Checkbox("Impedance", &app->smith_show_impedance);
    ImGui::SameLine();
    if (!app->smith_show_impedance) {
        ImGui::TextDisabled("(Admittance view)");
    }
    ImGui::Checkbox("VSWR Circles", &app->smith_show_vswr_circles);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Z0 (Ω)", &app->smith_z0, 1.0, 10.0, "%.1f");

    ImGui::Separator();

    if (ImPlot::BeginPlot("##SmithChart",
                          ImVec2(-1, -1),
                          ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxis(ImAxis_X1, nullptr, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.1, 1.1, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.1, 1.1, ImGuiCond_Always);
        draw_smith_grid(app);
        if (!app->sparam_s11_mag.empty()) {
            plot_smith_reflection_data(app);
        } else {
            ImPlot::EndPlot();
            ImGui::TextUnformatted("Load S-parameters to populate chart.");
            ImGui::End();
            return;
        }
        ImPlot::EndPlot();
    }

    int count = (int)app->sparam_freq_hz.size();
    if (count > 0) {
        if (app->smith_freq_index >= count) app->smith_freq_index = count - 1;
        if (app->smith_freq_index < 0) app->smith_freq_index = 0;
        ImGui::Separator();
        ImGui::TextUnformatted("Probe Frequency");
        if (ImGui::SliderInt("##smith_freq_idx", &app->smith_freq_index, 0, count - 1)) {
            /* slider moved */
        }
        int idx = app->smith_freq_index;
        double freq_ghz = app->sparam_freq_hz[idx] * 1e-9;
        double s11_db = app->sparam_s11_db[idx];
        double vswr = app->sparam_vswr[idx];
        ImGui::Text("Frequency: %.3f GHz", freq_ghz);
        ImGui::Text("S11: %.2f dB", s11_db);
        ImGui::Text("VSWR: %.2f", vswr);
    } else {
        ImGui::TextUnformatted("No S-parameter data available.");
    }

    ImGui::End();
}

static void draw_sparameter_window(AppState* app) {
    if (!app || !app->sparam_window_open) return;

    ImGui::SetNextWindowSize(ImVec2(850.0f, 620.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("S-Parameters", &app->sparam_window_open)) {
        ImGui::End();
        return;
    }

    ImGui::InputText("CSV Path", app->sparam_csv_path, sizeof(app->sparam_csv_path));
    ImGui::SameLine();
    if (ImGui::Button("Load CSV")) {
        if (load_sparameter_csv(app, app->sparam_csv_path)) {
            ui_log_add(app, "Loaded S-parameters: %s", app->sparam_csv_path);
        } else {
            ui_log_add(app, "Failed to load: %s", app->sparam_csv_path);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Data")) {
        clear_sparameter_data(app);
        ui_log_add(app, "Cleared S-parameter data");
    }

    if (!app->sparam_data_loaded || app->sparam_freq_hz.empty()) {
        ImGui::Separator();
        ImGui::TextUnformatted("No S-parameter data loaded.");
        ImGui::TextWrapped("Run CLI sweep to generate sweep_results.csv, then load it here.");
        ImGui::End();
        return;
    }

    const int count = (int)app->sparam_freq_hz.size();
    const double* freq_ghz = app->sparam_freq_ghz.data();
    const double* s21_db = app->sparam_s21_db.data();
    const double* s11_db = app->sparam_s11_db.data();
    const double* vswr = app->sparam_vswr.data();

    ImGui::Separator();
    if (ImPlot::BeginPlot("S21 - Transmission", ImVec2(-1, 260))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "S21 (dB)");
        ImPlot::PlotLine("S21", freq_ghz, s21_db, count);

        double threshold = app->sparam_peak_db - 3.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0, 0, 1));
        ImPlot::PlotInfLines("-3 dB", &threshold, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();

        ImPlotRect limits = ImPlot::GetPlotLimits();
        if (app->sparam_f_low_hz > 0.0) {
            double low = app->sparam_f_low_hz * 1e-9;
            double xs[2] = {low, low};
            double ys[2] = {limits.Y.Min, limits.Y.Max};
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0, 1, 0, 1));
            ImPlot::PlotLine("f_low", xs, ys, 2);
            ImPlot::PopStyleColor();
        }
        if (app->sparam_f_high_hz > app->sparam_f_low_hz) {
            double high = app->sparam_f_high_hz * 1e-9;
            double xs[2] = {high, high};
            double ys[2] = {limits.Y.Min, limits.Y.Max};
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0, 1, 0, 1));
            ImPlot::PlotLine("f_high", xs, ys, 2);
            ImPlot::PopStyleColor();
        }

        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("S11 - Reflection", ImVec2(-1, 220))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "S11 (dB)");
        ImPlot::PlotLine("S11", freq_ghz, s11_db, count);
        double match_line = -10.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 1, 0, 1));
        ImPlot::PlotInfLines("-10 dB", &match_line, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("VSWR", ImVec2(-1, 160))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "VSWR");
        ImPlot::SetupAxisLimits(ImAxis_Y1, 1.0, 5.0, ImGuiCond_Always);
        ImPlot::PlotLine("VSWR", freq_ghz, vswr, count);
        double vswr_limit = 2.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0.5f, 0, 1));
        ImPlot::PlotInfLines("VSWR=2", &vswr_limit, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();
        ImPlot::EndPlot();
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Bandwidth Analysis");
    ImGui::Columns(2, "sparam_stats", false);
    ImGui::SetColumnWidth(0, 160.0f);

    ImGui::Text("Peak S21:"); ImGui::NextColumn();
    ImGui::Text("%.2f dB @ %.3f GHz", app->sparam_peak_db, app->sparam_peak_freq_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("Center Freq:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_center_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("Bandwidth:"); ImGui::NextColumn();
    if (app->sparam_bandwidth_hz > 0.0 && app->sparam_f_center_hz > 0.0) {
        ImGui::Text("%.3f GHz (%.1f%%)",
                    app->sparam_bandwidth_hz * 1e-9,
                    (app->sparam_bandwidth_hz / app->sparam_f_center_hz) * 100.0);
    } else {
        ImGui::TextDisabled("N/A");
    }
    ImGui::NextColumn();

    ImGui::Text("f_low:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_low_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("f_high:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_high_hz * 1e-9);
    ImGui::Columns(1);

    ImGui::End();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim);

/* Log panel --------------------------------------------------------------- */
static void draw_log_panel(AppState* app);

static void render_grid_overlay(RenderContext* render, const SimulationState* sim, SDL_Color color) {
    if (!render || !sim) return;
    SDL_Renderer* rr = render->renderer;
    if (!rr) return;
    int width = sim->nx * render->scale;
    int height = sim->ny * render->scale;
    SDL_SetRenderDrawColor(rr, color.r, color.g, color.b, color.a);
    for (int i = 0; i <= sim->nx; ++i) {
        int x = i * render->scale;
        SDL_RenderDrawLine(rr, x, 0, x, height);
    }
    for (int j = 0; j <= sim->ny; ++j) {
        int y = j * render->scale;
        SDL_RenderDrawLine(rr, 0, y, width, y);
    }
}

static SDL_Color material_color_from_epsilon(double eps) {
    SDL_Color color{128, 128, 128, 255};
    if (eps < 1.05) {
        color = SDL_Color{200, 220, 255, 255};
        return color;
    }
    if (eps > 100.0) {
        color = SDL_Color{184, 115, 51, 255};
        return color;
    }
    const Material* best = NULL;
    double best_diff = 1e9;
    int count = material_library_get_count();
    for (int i = 0; i < count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;
        if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) continue;
        double diff = std::fabs(mat->epsilon_r - eps);
        if (diff < best_diff) {
            best_diff = diff;
            best = mat;
        }
    }
    if (best && best_diff < 0.75) {
        color = SDL_Color{best->color_r, best->color_g, best->color_b, 255};
    }
    return color;
}

static void render_material_distribution(SDL_Renderer* renderer,
                                         const SimulationState* sim,
                                         int scale) {
    if (!renderer || !sim || scale <= 0) return;
    for (int j = 0; j < sim->ny; ++j) {
        for (int i = 0; i < sim->nx; ++i) {
            double eps = fdtd_epsilon_at(sim, i, j);
            SDL_Color color = material_color_from_epsilon(eps);
            SDL_Rect rect = {i * scale, j * scale, scale, scale};
            SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
            SDL_RenderFillRect(renderer, &rect);
        }
    }
}

static void render_material_overlay(SDL_Renderer* renderer,
                                    const SimulationState* sim,
                                    int scale,
                                    float alpha) {
    if (!renderer || !sim || scale <= 0) return;
    if (alpha <= 0.0f) return;
    Uint8 a = (Uint8)(std::fmin(std::fmax(alpha, 0.0f), 1.0f) * 255.0f);
    for (int j = 0; j < sim->ny; ++j) {
        for (int i = 0; i < sim->nx; ++i) {
            double eps = fdtd_epsilon_at(sim, i, j);
            SDL_Color base = material_color_from_epsilon(eps);
            SDL_Rect rect = {i * scale, j * scale, scale, scale};
            SDL_SetRenderDrawColor(renderer, base.r, base.g, base.b, a);
            SDL_RenderFillRect(renderer, &rect);
        }
    }
}

static void render_material_outlines(SDL_Renderer* renderer,
                                     const SimulationState* sim,
                                     int scale) {
    if (!renderer || !sim || scale <= 0) return;
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 160);
    for (int j = 0; j < sim->ny - 1; ++j) {
        for (int i = 0; i < sim->nx - 1; ++i) {
            double eps_here = fdtd_epsilon_at(sim, i, j);
            double eps_right = fdtd_epsilon_at(sim, i + 1, j);
            double eps_down = fdtd_epsilon_at(sim, i, j + 1);
            if (std::fabs(eps_here - eps_right) > 0.1) {
                int x = (i + 1) * scale;
                int y0 = j * scale;
                int y1 = (j + 1) * scale;
                SDL_RenderDrawLine(renderer, x, y0, x, y1);
            }
            if (std::fabs(eps_here - eps_down) > 0.1) {
                int x0 = i * scale;
                int x1 = (i + 1) * scale;
                int y = (j + 1) * scale;
                SDL_RenderDrawLine(renderer, x0, y, x1, y);
            }
        }
    }
}

static void apply_paint_brush(SimulationState* sim,
                              int gx,
                              int gy,
                              int brush_radius,
                              int paint_type,
                              double eps) {
    if (!sim) return;
    if (brush_radius < 0) brush_radius = 0;

    int nx = sim->nx;
    int ny = sim->ny;

    int r = brush_radius;
    int r2 = r * r;
    for (int di = -r; di <= r; ++di) {
        int i = gx + di;
        if (i < 0 || i >= nx) continue;
        for (int dj = -r; dj <= r; ++dj) {
            int j = gy + dj;
            if (j < 0 || j >= ny) continue;
            if (di * di + dj * dj <= r2) {
                paint_material_at(sim, i, j, paint_type, eps);
            }
        }
    }
}

static void blocks_paint_pec(SimulationState* sim, int gx, int gy, int brush_radius) {
    apply_paint_brush(sim, gx, gy, brush_radius, 1, 1.0);
}

static void blocks_paint_pmc(SimulationState* sim, int gx, int gy, int brush_radius) {
    apply_paint_brush(sim, gx, gy, brush_radius, 2, 1.0);
}

static void blocks_paint_dielectric(SimulationState* sim,
                                    int gx,
                                    int gy,
                                    int brush_radius,
                                    double epsr) {
    apply_paint_brush(sim, gx, gy, brush_radius, 3, epsr);
}

static void apply_selected_paint(SimulationState* sim,
                                 int gx,
                                 int gy,
                                 int brush_radius,
                                 int paint_material_type,
                                 double paint_epsilon) {
    if (!sim) return;
    switch (paint_material_type) {
        case 0:
            blocks_paint_pec(sim, gx, gy, brush_radius);
            break;
        case 1:
            blocks_paint_pmc(sim, gx, gy, brush_radius);
            break;
        case 2:
        default:
            blocks_paint_dielectric(sim, gx, gy, brush_radius, paint_epsilon);
            break;
    }
}

static void apply_paint_with_material(SimulationState* sim,
                                      int gx,
                                      int gy,
                                      int brush_radius,
                                      int material_id) {
    if (!sim) return;
    const Material* mat = material_library_get_by_id(material_id);
    if (!mat) return;

    int paint_type = 0;
    double epsilon = 1.0;

    switch (mat->type) {
        case MAT_TYPE_PEC:
            paint_type = 1;
            epsilon = 1.0;
            break;
        case MAT_TYPE_PMC:
            paint_type = 2;
            epsilon = 1.0;
            break;
        case MAT_TYPE_DIELECTRIC:
        case MAT_TYPE_LOSSY_DIEL:
        default:
            paint_type = 3;
            epsilon = mat->epsilon_r;
            break;
    }

    apply_paint_brush(sim, gx, gy, brush_radius, paint_type, epsilon);
}

static int find_source_at_position(const SimulationState* sim, int gx, int gy, int tolerance) {
    if (!sim) return -1;
    if (tolerance < 0) tolerance = 0;
    int tol2 = tolerance * tolerance;
    for (int i = 0; i < MAX_SRC; ++i) {
        if (!sim->sources[i].active) continue;
        int sx = sim->sources[i].ix;
        int sy = sim->sources[i].iy;
        int dx = gx - sx;
        int dy = gy - sy;
        if (dx * dx + dy * dy <= tol2) {
            return i;
        }
    }
    return -1;
}

// Theme and visual styling helpers -----------------------------------------

static void apply_theme(int theme_id) {
    if (theme_id == 0) {
        ImGui::StyleColorsDark();
    } else {
        ImGui::StyleColorsLight();
    }

    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.ScrollbarSize = 14.0f;
    style.FramePadding = ImVec2(8.0f, 4.0f);
    style.ItemSpacing = ImVec2(8.0f, 6.0f);
    style.ItemInnerSpacing = ImVec2(6.0f, 4.0f);
    style.IndentSpacing = 20.0f;
    style.WindowPadding = ImVec2(8.0f, 8.0f);
}

static const ImVec4 ACCENT_PALETTES[6][3] = {
    { ImVec4(0.26f, 0.59f, 0.98f, 1.00f),
      ImVec4(0.36f, 0.69f, 1.00f, 1.00f),
      ImVec4(0.16f, 0.49f, 0.88f, 1.00f) },

    { ImVec4(0.20f, 0.70f, 0.50f, 1.00f),
      ImVec4(0.30f, 0.80f, 0.60f, 1.00f),
      ImVec4(0.10f, 0.60f, 0.40f, 1.00f) },

    { ImVec4(0.90f, 0.40f, 0.20f, 1.00f),
      ImVec4(1.00f, 0.50f, 0.30f, 1.00f),
      ImVec4(0.80f, 0.30f, 0.10f, 1.00f) },

    { ImVec4(0.70f, 0.30f, 0.80f, 1.00f),
      ImVec4(0.80f, 0.40f, 0.90f, 1.00f),
      ImVec4(0.60f, 0.20f, 0.70f, 1.00f) },

    { ImVec4(0.90f, 0.20f, 0.30f, 1.00f),
      ImVec4(1.00f, 0.30f, 0.40f, 1.00f),
      ImVec4(0.80f, 0.10f, 0.20f, 1.00f) },

    { ImVec4(0.90f, 0.70f, 0.10f, 1.00f),
      ImVec4(1.00f, 0.80f, 0.20f, 1.00f),
      ImVec4(0.80f, 0.60f, 0.00f, 1.00f) }
};

static const char* ACCENT_NAMES[6] = {
    "Blue", "Teal", "Orange", "Purple", "Red", "Gold"
};

static void apply_accent_palette(int accent_id) {
    if (accent_id < 0 || accent_id >= 6) accent_id = 0;

    ImGuiStyle& style = ImGui::GetStyle();
    const ImVec4* palette = ACCENT_PALETTES[accent_id];

    style.Colors[ImGuiCol_Header] = palette[0];
    style.Colors[ImGuiCol_HeaderHovered] = palette[1];
    style.Colors[ImGuiCol_HeaderActive] = palette[2];

    style.Colors[ImGuiCol_Tab] =
        ImVec4(palette[0].x * 0.7f, palette[0].y * 0.7f, palette[0].z * 0.7f, 0.86f);
    style.Colors[ImGuiCol_TabHovered] = palette[1];
    style.Colors[ImGuiCol_TabActive] = palette[2];

    style.Colors[ImGuiCol_Button] =
        ImVec4(palette[0].x * 0.6f, palette[0].y * 0.6f, palette[0].z * 0.6f, 0.40f);
    style.Colors[ImGuiCol_ButtonHovered] =
        ImVec4(palette[1].x * 0.8f, palette[1].y * 0.8f, palette[1].z * 0.8f, 1.00f);
    style.Colors[ImGuiCol_ButtonActive] = palette[2];

    style.Colors[ImGuiCol_CheckMark] = palette[2];
    style.Colors[ImGuiCol_SliderGrab] = palette[0];
    style.Colors[ImGuiCol_SliderGrabActive] = palette[2];
}

static SDL_Color theme_viewport_clear_color(int theme_id) {
    SDL_Color color;
    if (theme_id == 0) {
        color.r = 8;
        color.g = 8;
        color.b = 10;
        color.a = 255;
    } else {
        color.r = 235;
        color.g = 236;
        color.b = 240;
        color.a = 255;
    }
    return color;
}

static const char* kColormapNames[] = {"Classic", "Viridis", "Plasma"};
static const int kColormapCount = 3;

static bool rebootstrap_simulation(const SimulationConfig* cfg,
                                   SimulationBootstrap* bootstrap,
                                   SimulationState** sim_out,
                                   Scope* scope,
                                   int scale) {
    if (!cfg || !bootstrap || !sim_out || !scope) return false;

    simulation_bootstrap_shutdown(bootstrap);

    if (!simulation_bootstrap_with_config(cfg, bootstrap)) {
        std::fprintf(stderr, "Failed to reinitialise simulation from wizard configuration\n");
        return false;
    }

    *sim_out = bootstrap->sim;

    scope_free(scope);
    if (!scope_init(scope, (*sim_out)->nx * scale)) {
        std::fprintf(stderr, "Failed to initialise oscilloscope after rebootstrap\n");
        simulation_bootstrap_shutdown(bootstrap);
        *sim_out = NULL;
        return false;
    }

    return true;
}

static void apply_wizard_materials_to_sim(const WizardState& wizard,
                                          SimulationBootstrap* bootstrap,
                                          SimulationState* sim) {
    if (!bootstrap || !sim) return;

    int count = wizard.cfg.material_rect_count;
    if (count < 0) count = 0;
    if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;

    bootstrap->config.material_rect_count = count;
    sim->config.material_rect_count = count;
    for (int i = 0; i < count; ++i) {
        bootstrap->config.material_rects[i] = wizard.cfg.material_rects[i];
        sim->config.material_rects[i] = wizard.cfg.material_rects[i];
    }

    materials_init(sim);
}

static void sync_sim_to_wizard_config(SimulationState* sim, WizardState* wizard) {
    if (!sim || !wizard) return;

    int nx = sim->nx;
    int ny = sim->ny;
    double nx1 = (nx > 1) ? (double)(nx - 1) : 1.0;
    double ny1 = (ny > 1) ? (double)(ny - 1) : 1.0;

    int max_src = wizard->cfg.source_count;
    if (max_src < 0) max_src = 0;
    if (max_src > MAX_SRC) max_src = MAX_SRC;

    for (int i = 0; i < max_src; ++i) {
        const Source& s = sim->sources[i];
        SourceConfigSpec* spec = &wizard->cfg.source_configs[i];
        spec->active = s.active ? 1 : 0;
        spec->x = (double)s.ix / nx1;
        spec->y = (double)s.iy / ny1;
        spec->type = s.type;
        spec->field = s.field;
        spec->amp = s.amp;
        spec->freq = s.freq;
        spec->sigma2 = s.sigma2;
        /* Expr text is not round-tripped here; leave cfg.expr as-is. */
    }
}

static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Materials / Blocks")) return;

    ImGui::Indent();
    int count = wizard.cfg.material_rect_count;
    if (count < 0) count = 0;
    if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;
    wizard.cfg.material_rect_count = count;

    ImGui::Text("Blocks: %d / %d", count, CONFIG_MAX_MATERIAL_RECTS);
    ImGui::SameLine();
    bool can_add = (count < CONFIG_MAX_MATERIAL_RECTS);
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ Add Block")) {
        MaterialRectSpec& r = wizard.cfg.material_rects[count];
        r.x0 = 0.25;
        r.y0 = 0.25;
        r.x1 = 0.75;
        r.y1 = 0.75;
        r.epsr = 4.0;
        r.sigma = 0.0;
        r.tag = 0;
        wizard.cfg.material_rect_count = count + 1;
        if (sim && bootstrap) {
            apply_wizard_materials_to_sim(wizard, bootstrap, sim);
        }
        if (app) {
            ui_log_add(app, "Added block #%d (default material)", count);
        }
    }
    ImGui::SameLine();
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ Add Block from Material")) {
        ImGui::OpenPopup("AddBlockMaterialPopup");
    }
    if (!can_add) ImGui::EndDisabled();
    if (ImGui::BeginPopup("AddBlockMaterialPopup")) {
        ImGui::TextUnformatted("Select material:");
        ImGui::Separator();
        int mat_count = material_library_get_count();
        for (int i = 0; i < mat_count; ++i) {
            const Material* mat = material_library_get_by_index(i);
            if (!mat) continue;
            if (ImGui::Selectable(mat->name)) {
                if (wizard.cfg.material_rect_count < CONFIG_MAX_MATERIAL_RECTS) {
                    int idx = wizard.cfg.material_rect_count;
                    MaterialRectSpec& r = wizard.cfg.material_rects[idx];
                    r.x0 = 0.3;
                    r.y0 = 0.3;
                    r.x1 = 0.7;
                    r.y1 = 0.7;
                    apply_material_to_rect(&r, mat);
                    wizard.cfg.material_rect_count = idx + 1;
                    if (sim && bootstrap) {
                        apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                    }
                    if (app) {
                        ui_log_add(app, "Added block #%d (%s)", idx, mat->name);
                    }
                }
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    ImGui::Separator();

    const Material* selected_mat = nullptr;
    if (app && app->selected_material_id >= 0) {
        selected_mat = material_library_get_by_id(app->selected_material_id);
    }

    if (selected_mat) {
        ImGui::Checkbox("Filter by selected material", &app->filter_blocks_by_material);
        ImGui::SameLine();
        ImGui::TextDisabled("Filter: %s", selected_mat->name);
        ImGui::Checkbox("Auto-filter on selection", &app->auto_filter_blocks_on_select);
    } else {
        ImGui::BeginDisabled();
        bool dummy = false;
        ImGui::Checkbox("Filter by selected material", &dummy);
        ImGui::SameLine();
        ImGui::TextDisabled("No material selected");
        ImGui::EndDisabled();
        app->filter_blocks_by_material = false;
    }

    ImGui::Separator();

    if (wizard.cfg.material_rect_count == 0) {
        ImGui::TextDisabled("No material blocks configured.");
    } else {
        for (int i = 0; i < wizard.cfg.material_rect_count; ++i) {
            MaterialRectSpec& r = wizard.cfg.material_rects[i];
            if (app && app->filter_blocks_by_material && selected_mat) {
                if (!rect_matches_material(r, selected_mat)) {
                    ImGui::PopID();
                    continue;
                }
            }
            ImGui::PushID(i);
            bool tree_open = ImGui::CollapsingHeader("", ImGuiTreeNodeFlags_DefaultOpen);
            ImGui::SameLine();
            ImGui::Text("Block %d", i);
            if (tree_open) {
                const Material* current_mat = guess_material_for_rect(r);
                ImVec4 color = current_mat
                                   ? ImVec4(current_mat->color_r / 255.0f,
                                            current_mat->color_g / 255.0f,
                                            current_mat->color_b / 255.0f,
                                            1.0f)
                                   : ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                if (selected_mat && rect_matches_material(r, selected_mat)) {
                    ImGui::PushStyleColor(ImGuiCol_Text, color);
                }
                ImGui::ColorButton("##blkcol",
                                   color,
                                   ImGuiColorEditFlags_NoTooltip |
                                       ImGuiColorEditFlags_NoDragDrop,
                                   ImVec2(16, 16));
                ImGui::SameLine();
                if (selected_mat && rect_matches_material(r, selected_mat)) {
                    ImGui::PopStyleColor();
                }
                const char* combo_label =
                    current_mat ? current_mat->name
                                : (r.tag == 1 ? "PEC"
                                              : (r.tag == 2 ? "PMC" : "Custom"));
                if (ImGui::BeginCombo("Material", combo_label)) {
                    int mat_count = material_library_get_count();
                    for (int m = 0; m < mat_count; ++m) {
                        const Material* mat = material_library_get_by_index(m);
                        if (!mat) continue;
                        bool selected = (current_mat == mat);
                        if (ImGui::Selectable(mat->name, selected)) {
                            apply_material_to_rect(&r, mat);
                            if (sim && bootstrap) {
                                apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                            }
                        }
                        if (selected) ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }

                const char* type_items[] = {"Dielectric", "PEC", "PMC"};
                int type_idx = (r.tag == 1) ? 1 : (r.tag == 2) ? 2 : 0;
                if (ImGui::Combo("Type", &type_idx, type_items, IM_ARRAYSIZE(type_items))) {
                    r.tag = (type_idx == 1) ? 1 : (type_idx == 2) ? 2 : 0;
                    if (sim && bootstrap) apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                }

                bool changed = false;
                float x0 = (float)r.x0;
                float y0 = (float)r.y0;
                float x1 = (float)r.x1;
                float y1 = (float)r.y1;
                changed |= ImGui::SliderFloat("x0", &x0, 0.0f, 0.99f, "%.3f");
                changed |= ImGui::SliderFloat("y0", &y0, 0.0f, 0.99f, "%.3f");
                changed |= ImGui::SliderFloat("x1", &x1, 0.0f, 1.0f, "%.3f");
                changed |= ImGui::SliderFloat("y1", &y1, 0.0f, 1.0f, "%.3f");
                if (changed) {
                    r.x0 = clamp01(x0);
                    r.y0 = clamp01(y0);
                    r.x1 = clamp01(x1);
                    r.y1 = clamp01(y1);
                    if (r.x1 < r.x0) r.x1 = r.x0 + 0.01;
                    if (r.y1 < r.y0) r.y1 = r.y0 + 0.01;
                }

                if (sim) {
                    ImGui::TextDisabled("Meters: (%.3f, %.3f) to (%.3f, %.3f)",
                                        r.x0 * sim->lx,
                                        r.y0 * sim->ly,
                                        r.x1 * sim->lx,
                                        r.y1 * sim->ly);
                }

                if (r.tag == 0) {
                    if (ImGui::InputDouble("epsr", &r.epsr, 0.1, 1.0)) changed = true;
                    if (ImGui::InputDouble("sigma", &r.sigma, 0.0, 0.1)) changed = true;
                } else {
                    ImGui::TextDisabled("PEC/PMC override eps/sigma");
                }

                if (changed && sim && bootstrap) {
                    apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                }

                if (ImGui::Button("Delete Block")) {
                    remove_block(&wizard, bootstrap, sim, app, i);
                    ImGui::PopID();
                    break;
                }
            }
            ImGui::PopID();
        }
    }
    ImGui::Unindent();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Probes")) return;

    ImGui::Indent();
#if EMWAVE_ENABLE_PORTS
    if (!sim) {
        ImGui::Unindent();
        return;
    }
    ImGui::Text("Ports enabled: %s", sim->ports_on ? "yes" : "no");
    ImGui::Separator();
    for (int p = 0; p < MAX_PORTS; ++p) {
        const Port& port = sim->ports[p];
        if (!port.active) continue;
        ImGui::Text("Port %d: x=%d, y0=%d, y1=%d, n=%d", p, port.x, port.y0, port.y1, port.n);
    }
#else
    ImGui::TextUnformatted("Ports are disabled in this build.");
#endif
    ImGui::Unindent();
}

/* Log panel --------------------------------------------------------------- */
static void draw_log_panel(AppState* app) {
    if (!ImGui::Begin("Log")) {
        ImGui::End();
        return;
    }
    if (!app) {
        ImGui::End();
        return;
    }
    for (int i = 0; i < app->log_count; ++i) {
        ImGui::TextUnformatted(app->log_lines[i]);
    }
    ImGui::End();
}

int main(int argc, char** argv) {
    SimulationBootstrap bootstrap;
    int bootstrap_status = simulation_bootstrap_from_args(argc, argv, &bootstrap);
    if (bootstrap_status == 0) {
        return 0;  /* help shown or invalid config */
    }
    if (bootstrap_status < 0) {
        return 1;
    }

    SimulationState* sim = bootstrap.sim;
    if (!sim) {
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }

    material_library_init();

    WizardState wizard;
    wizard_init_from_config(wizard, &bootstrap.config);

    AppState app = {};
    app.basic_mode = true;
    app.selected_source = 0;
    app.selected_block = 0;
    app.placing_source = false;
    app.placing_block = false;
    app.block_first_set = false;
    app.block_first_i = app.block_first_j = -1;
    app.last_click_i = app.last_click_j = -1;
    app.show_scope_window = false;
    app.show_scene_panel = true;
    app.show_sources_panel = true;
    app.show_blocks_panel = true;
    app.show_run_panel = true;
    app.show_grid_panel = true;
    app.show_run_settings_panel = true;
    app.show_probes_panel = true;
    app.show_log_panel = true;
    app.show_expression_panel = true;
    app.show_scenes_panel = true;
    app.show_grid_overlay = true;
    app.log_count = 0;
    app.viewport_pos = ImVec2(0.0f, 0.0f);
    app.viewport_size = ImVec2(0.0f, 0.0f);
    app.viewport_valid = false;
    app.selected_material_id = -1;
    app.paint_material_id = -1;
    app.material_browser_open = false;
    app.material_search[0] = '\0';
    app.filter_metals = true;
    app.filter_dielectrics = true;
    app.show_material_legend = true;
    app.filter_blocks_by_material = false;
    app.auto_filter_blocks_on_select = true;
    app.highlight_blocks_by_material = true;
    app.filter_blocks_by_material = false;
    app.visualization_mode = 0;
    app.show_material_outlines = false;
    app.material_overlay_alpha = 0.5f;
    app.sparam_window_open = false;
    app.sparam_data_loaded = false;
    app.sparam_peak_db = 0.0;
    app.sparam_peak_freq_hz = 0.0;
    app.sparam_f_low_hz = 0.0;
    app.sparam_f_high_hz = 0.0;
    app.sparam_f_center_hz = 0.0;
    app.sparam_bandwidth_hz = 0.0;
    std::strncpy(app.sparam_csv_path, "sweep_results.csv", sizeof(app.sparam_csv_path));
    app.sparam_csv_path[sizeof(app.sparam_csv_path) - 1] = '\0';
    app.smith_chart_open = false;
    app.smith_show_impedance = true;
    app.smith_show_vswr_circles = true;
    app.smith_z0 = 50.0;
    app.smith_freq_index = 0;
    app.request_rebootstrap = false;
    app.rebootstrap_message[0] = '\0';
    ui_log_add(&app, "Simulation started");

    int scale = 2;
    int width = sim->nx * scale;
    int height = sim->ny * scale;
    if (width < 1920) width = 1920;
    if (height < 1080) height = 1080;

    RenderContext* render = render_init("emwave-c (ImGui)", width, height);
    if (!render) {
        std::fprintf(stderr, "Failed to initialise SDL renderer for ImGui front-end\n");
        material_library_shutdown();
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }
    render->scale = scale;

    // Set up Dear ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable docking
    io.FontGlobalScale = 1.0f;
    io.IniFilename = "imgui.ini";  // Enable layout persistence

    ImGui_ImplSDL2_InitForSDLRenderer(render->window, render->renderer);
    ImGui_ImplSDLRenderer2_Init(render->renderer);

    // Oscilloscope for a single probe in the center
    Scope scope = {0};
    if (!scope_init(&scope, sim->nx * scale)) {
        std::fprintf(stderr, "Failed to initialise oscilloscope buffer\n");
        ImGui_ImplSDLRenderer2_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImGui::DestroyContext();
        render_free(render);
        material_library_shutdown();
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }

    bool running = true;
    bool paused = false;
    int steps_per_frame = 5;      // Speed control
    bool paint_mode = false;
    int paint_material_id = -1;   // Material ID for painting (-1 legacy)
    int paint_material_type = 0;  // Legacy: 0=PEC, 1=PMC, 2=dielectric
    double paint_epsilon = 4.0;   // Legacy dielectric epsilon
    int paint_brush_size = 3;     // Radius in grid cells
    bool show_help_overlay = false;
    bool dragging_source = false;
    int dragged_source_idx = -1;
    int current_theme = 0;        // 0=Dark, 1=Light
    int current_colormap = 0;     // 0=Classic, 1=Viridis, 2=Plasma
    int current_accent = 0;       // 0-5 accent palette index
    SDL_Color viewport_clear_color = theme_viewport_clear_color(current_theme);
    bool auto_rescale = true;
    bool hold_color = false;
    bool hold_scope = false;
    double vmax_smooth = sim->step_Ez_absmax;
    if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
    double scope_vmax = 1.0;

    // Apply initial theme, accent, and colormap to both ImGui and SDL renderer
    apply_theme(current_theme);
    apply_accent_palette(current_accent);
    ui_render_set_theme((ThemeMode)current_theme, current_accent);
    ui_render_set_colormap((ColorMapMode)current_colormap);

    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 prev = SDL_GetPerformanceCounter();
    double fps_avg = 0.0;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            ImGui_ImplSDL2_ProcessEvent(&e);
            if (e.type == SDL_QUIT) {
                running = false;
            } else if (e.type == SDL_KEYDOWN) {
                SDL_Keycode key = e.key.keysym.sym;

                // ============================================================
                // GLOBAL SHORTCUTS - Always work, regardless of ImGui focus
                // ============================================================
                bool handled_globally = false;

                switch (key) {
                    case SDLK_ESCAPE:
                    case SDLK_q:
                        running = false;
                        handled_globally = true;
                        break;

                    case SDLK_F1:
                        show_help_overlay = !show_help_overlay;
                        handled_globally = true;
                        break;

                    case SDLK_F2:
                        if (save_screenshot(render, "frame.bmp")) {
                            ui_log_add(&app, "Screenshot saved: frame.bmp");
                        } else {
                            ui_log_add(&app, "Screenshot failed");
                        }
                        handled_globally = true;
                        break;

                    case SDLK_F3:
                        if (dump_scope_fft_csv(&scope, "scope_fft.csv", sim->dt, 1024)) {
                            ui_log_add(&app, "FFT exported: scope_fft.csv");
                        } else {
                            ui_log_add(&app, "FFT export failed");
                        }
                        handled_globally = true;
                        break;

                    case SDLK_F5: {
                        // Load waveguide preset
                        SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                        char errbuf[256];
                        if (config_loader_parse_file("configs/waveguide.json",
                                                     &cfg,
                                                     errbuf,
                                                     sizeof(errbuf))) {
                            if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                wizard_init_from_config(wizard, &cfg);

                                // Resize window to match new grid
                                width = sim->nx * scale;
                                height = sim->ny * scale;
                                if (width < 1920) width = 1920;
                                if (height < 1080) height = 1080;
                                SDL_SetWindowSize(render->window, width, height);

                                // Reset visualization scaling
                                vmax_smooth = sim->step_Ez_absmax;
                                if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                                scope_vmax = 1.0;
                                auto_rescale = true;
                                hold_color = false;
                                hold_scope = false;

                                ui_log_add(&app, "Loaded: waveguide.json");
                            } else {
                                ui_log_add(&app, "Failed to load: waveguide.json");
                            }
                        } else {
                            ui_log_add(&app, "Failed to parse: waveguide.json");
                        }
                        handled_globally = true;
                        break;
                    }

                    case SDLK_F6: {
                        // Load CPW filter preset
                        SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                        char errbuf[256];
                        if (config_loader_parse_file("configs/cpw_filter.json",
                                                     &cfg,
                                                     errbuf,
                                                     sizeof(errbuf))) {
                            if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                wizard_init_from_config(wizard, &cfg);

                                // Resize window to match new grid
                                width = sim->nx * scale;
                                height = sim->ny * scale;
                                if (width < 1920) width = 1920;
                                if (height < 1080) height = 1080;
                                SDL_SetWindowSize(render->window, width, height);

                                // Reset visualization scaling
                                vmax_smooth = sim->step_Ez_absmax;
                                if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                                scope_vmax = 1.0;
                                auto_rescale = true;
                                hold_color = false;
                                hold_scope = false;

                                ui_log_add(&app, "Loaded: cpw_filter.json");
                            } else {
                                ui_log_add(&app, "Failed to load: cpw_filter.json");
                            }
                        } else {
                            ui_log_add(&app, "Failed to parse: cpw_filter.json");
                        }
                        handled_globally = true;
                        break;
                    }

                    default:
                        break;
                }

                // ============================================================
                // Handle Ctrl+S (also global)
                // ============================================================
                if (!handled_globally &&
                    (SDL_GetModState() & KMOD_CTRL) &&
                    key == SDLK_s) {
                    sync_sim_to_wizard_config(sim, &wizard);
                    const char* save_path = "configs/saved_config.json";
                    if (save_config_json(save_path, &wizard.cfg)) {
                        ui_log_add(&app, "Configuration saved: %s", save_path);
                    } else {
                        ui_log_add(&app, "Failed to save configuration");
                    }
                    handled_globally = true;
                }

                // ============================================================
                // CONTEXT SHORTCUTS - Only work when ImGui doesn't want keyboard
                // ============================================================
                if (!handled_globally && !io.WantCaptureKeyboard) {
                    switch (key) {
                        // Simulation controls
                        case SDLK_SPACE:
                            paused = !paused;
                            break;

                        // Source toggles
                        case SDLK_1:
                            if (MAX_SRC > 0) sim->sources[0].active = !sim->sources[0].active;
                            break;
                        case SDLK_2:
                            if (MAX_SRC > 1) sim->sources[1].active = !sim->sources[1].active;
                            break;
                        case SDLK_3:
                            if (MAX_SRC > 2) sim->sources[2].active = !sim->sources[2].active;
                            break;

                        // Source type cycling
                        case SDLK_t:
                        case SDLK_f:
                            sources_cycle_type(sim->sources);
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Cycled source type");
                            break;

                        // Simulation controls
                        case SDLK_r:
                            fdtd_reset(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Simulation reset");
                            break;
                        case SDLK_c:
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Fields cleared");
                            break;

                        // Auto-rescale modes
                        case SDLK_a:
                            auto_rescale = true;
                            hold_color = false;
                            hold_scope = false;
                            ui_log_add(&app, "Rescale mode: Auto (A)");
                            break;
                        case SDLK_h:
                            auto_rescale = false;
                            hold_color = true;
                            hold_scope = false;
                            ui_log_add(&app, "Rescale mode: Hold Color (H)");
                            break;
                        case SDLK_j:
                            auto_rescale = false;
                            hold_color = false;
                            hold_scope = true;
                            ui_log_add(&app, "Rescale mode: Hold Scope (J)");
                            break;
                        case SDLK_l:
                            auto_rescale = false;
                            hold_color = false;
                            hold_scope = false;
                            ui_log_add(&app, "Rescale mode: Manual (L)");
                            break;

                        // Boundary controls
                        case SDLK_y: {
                            BoundaryType type = boundary_get_type(sim);
                            type = (type == BOUNDARY_CPML) ? BOUNDARY_MUR : BOUNDARY_CPML;
                            boundary_set_type(sim, type);
                            if (boundary_is_cpml_enabled(sim)) {
                                cpml_build_coeffs(sim);
                                cpml_zero_psi(sim);
                            }
                            const char* name = (type == BOUNDARY_CPML) ? "CPML" : "Mur";
                            ui_log_add(&app, "Boundary: %s", name);
                            break;
                        }
                        case SDLK_7:
                        case SDLK_8:
                        case SDLK_9: {
                            int idx = (key == SDLK_7) ? 0 : (key == SDLK_8) ? 1 : 2;
                            boundary_set_type(sim, BOUNDARY_CPML);
                            cpml_apply_preset(sim, idx);
                            if (boundary_is_cpml_enabled(sim)) {
                                cpml_zero_psi(sim);
                            }
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "CPML preset %d applied", idx + 1);
                            break;
                        }

                        // S-parameter ports
                        case SDLK_s:
                            sim->ports_on = !sim->ports_on;
                            ui_log_add(&app, "Ports: %s", sim->ports_on ? "ON" : "OFF");
                            break;

                        // Material painting mode
                        case SDLK_m:
                        case SDLK_u:
                            paint_mode = !paint_mode;
                            ui_log_add(&app, "Paint mode: %s", paint_mode ? "ON" : "OFF");
                            break;
                        case SDLK_i:
                            if (paint_mode && paint_material_id < 0) {
                                paint_material_type = (paint_material_type + 1) % 3;
                                const char* mlabel = (paint_material_type == 0)
                                                         ? "PEC"
                                                         : (paint_material_type == 1) ? "PMC"
                                                                                      : "Dielectric";
                                ui_log_add(&app, "Paint type: %s (legacy mode)", mlabel);
                            }
                            break;
                        case SDLK_o:
                            if (paint_mode && paint_material_type == 2 && paint_material_id < 0) {
                                paint_epsilon -= 0.5;
                                if (paint_epsilon < 1.0) paint_epsilon = 1.0;
                                ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            }
                            break;
                        case SDLK_p:
                            if (paint_mode && paint_material_type == 2 && paint_material_id < 0) {
                                paint_epsilon += 0.5;
                                if (paint_epsilon > 20.0) paint_epsilon = 20.0;
                                ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            }
                            break;

                        // Visual controls
                        case SDLK_b:
                            current_theme = (current_theme == 0) ? 1 : 0;
                            apply_theme(current_theme);
                            apply_accent_palette(current_accent);
                            viewport_clear_color = theme_viewport_clear_color(current_theme);
                            ui_render_set_theme((ThemeMode)current_theme, current_accent);
                            ui_log_add(&app, "Theme: %s",
                                       (current_theme == 0) ? "Dark" : "Light");
                            break;
                        case SDLK_k: {
                            int shift = (SDL_GetModState() & KMOD_SHIFT) ? -1 : 1;
                            current_colormap =
                                (current_colormap + shift + kColormapCount) % kColormapCount;
                            ui_render_set_colormap((ColorMapMode)current_colormap);
                            ui_log_add(&app, "Colormap: %s",
                                       kColormapNames[current_colormap]);
                            break;
                        }
                        case SDLK_v: {
                            int shift = (SDL_GetModState() & KMOD_SHIFT) ? -1 : 1;
                            current_accent =
                                (current_accent + shift + 6) % 6;
                            apply_accent_palette(current_accent);
                            ui_render_set_theme((ThemeMode)current_theme, current_accent);
                            ui_log_add(&app, "Accent: %s", ACCENT_NAMES[current_accent]);
                            break;
                        }

                        // Grid overlay
                        case SDLK_g:
                            app.show_grid_overlay = !app.show_grid_overlay;
                            ui_log_add(&app, "Grid overlay: %s",
                                       app.show_grid_overlay ? "ON" : "OFF");
                            break;
                        case SDLK_e: {
                            app.visualization_mode = (app.visualization_mode + 1) % 3;
                            const char* mode_names[] = {"Field", "Material", "Overlay"};
                            int vm = app.visualization_mode;
                            if (vm < 0 || vm > 2) vm = 0;
                            ui_log_add(&app, "View mode: %s", mode_names[vm]);
                            break;
                        }

                        case SDLK_KP_1:
                        case SDLK_KP_2:
                        case SDLK_KP_3:
                        case SDLK_KP_4:
                        case SDLK_KP_5:
                        case SDLK_KP_6:
                        case SDLK_KP_7:
                        case SDLK_KP_8:
                            if (paint_mode) {
                                int idx = key - SDLK_KP_1;
                                const Material* mat = material_library_get_by_index(idx);
                                if (mat) {
                                    paint_material_id = mat->id;
                                    app.paint_material_id = paint_material_id;
                                    ui_log_add(&app, "Paint material: %s", mat->name);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                if (!app.viewport_valid) {
                    continue;
                }

                int mx = e.button.x;
                int my = e.button.y;
                float local_x = (float)mx - app.viewport_pos.x;
                float local_y = (float)my - app.viewport_pos.y;

                if (local_x >= 0.0f && local_y >= 0.0f &&
                    local_x < app.viewport_size.x && local_y < app.viewport_size.y) {
                    int ix = (int)(local_x / (float)scale);
                    int iy = (int)(local_y / (float)scale);
                    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) {
                        continue;
                    }

                    app.last_click_i = ix;
                    app.last_click_j = iy;

                    // PRIORITY 1: Check for source drag (only if NOT in paint mode)
                    if (!paint_mode) {
                        int hit_source = find_source_at_position(sim, ix, iy, 5);
                        if (hit_source >= 0) {
                            dragging_source = true;
                            dragged_source_idx = hit_source;
                            app.selected_source = hit_source;
                            if (hit_source < wizard.cfg.source_count) {
                                double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
                                double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
                                wizard.cfg.source_configs[hit_source].x =
                                    (double)sim->sources[hit_source].ix / nx1;
                                wizard.cfg.source_configs[hit_source].y =
                                    (double)sim->sources[hit_source].iy / ny1;
                            }
                            continue;  // Don't process other click actions
                        }
                    }

                    // PRIORITY 2: Paint mode (if active)
                        if (paint_mode) {
                            if (paint_material_id >= 0) {
                                apply_paint_with_material(sim, ix, iy, paint_brush_size, paint_material_id);
                            } else {
                                apply_selected_paint(sim,
                                                     ix,
                                                     iy,
                                                     paint_brush_size,
                                                     paint_material_type,
                                                     paint_epsilon);
                            }
                        }

                    // PRIORITY 3: Other interactions (block placement, etc.)
                    else if (app.placing_block) {
                        if (!app.block_first_set) {
                            app.block_first_set = true;
                            app.block_first_i = ix;
                            app.block_first_j = iy;
                        } else {
                            int idx = app.selected_block;
                            if (idx < 0) idx = 0;
                            if (idx >= CONFIG_MAX_MATERIAL_RECTS) {
                                idx = CONFIG_MAX_MATERIAL_RECTS - 1;
                            }
                            if (wizard.cfg.material_rect_count <= idx) {
                                wizard.cfg.material_rect_count = idx + 1;
                                if (wizard.cfg.material_rect_count >
                                    CONFIG_MAX_MATERIAL_RECTS) {
                                    wizard.cfg.material_rect_count =
                                        CONFIG_MAX_MATERIAL_RECTS;
                                }
                            }
                            MaterialRectSpec& r = wizard.cfg.material_rects[idx];
                            int i0 = (app.block_first_i < ix) ? app.block_first_i : ix;
                            int i1 = (app.block_first_i > ix) ? app.block_first_i : ix;
                            int j0 = (app.block_first_j < iy) ? app.block_first_j : iy;
                            int j1 = (app.block_first_j > iy) ? app.block_first_j : iy;
                            if (i0 < 0) i0 = 0;
                            if (i1 >= sim->nx) i1 = sim->nx - 1;
                            if (j0 < 0) j0 = 0;
                            if (j1 >= sim->ny) j1 = sim->ny - 1;
                            double nx_d = (double)sim->nx;
                            double ny_d = (double)sim->ny;
                            r.x0 = (double)i0 / nx_d;
                            r.x1 = (double)(i1 + 1) / nx_d;
                            r.y0 = (double)j0 / ny_d;
                            r.y1 = (double)(j1 + 1) / ny_d;
                            app.block_first_set = false;

                            apply_wizard_materials_to_sim(wizard, &bootstrap, sim);
                        }
                    }
                }
            } else if (e.type == SDL_MOUSEMOTION) {
                // Handle source dragging first
                if (dragging_source &&
                    dragged_source_idx >= 0 &&
                    dragged_source_idx < MAX_SRC) {
                    if (!app.viewport_valid) {
                        continue;
                    }

                    int mx = e.motion.x;
                    int my = e.motion.y;
                    float local_x = (float)mx - app.viewport_pos.x;
                    float local_y = (float)my - app.viewport_pos.y;

                    if (local_x < 0.0f || local_y < 0.0f ||
                        local_x >= app.viewport_size.x || local_y >= app.viewport_size.y) {
                        continue;
                    }

                    int ix = (int)(local_x / (float)scale);
                    int iy = (int)(local_y / (float)scale);

                    // Clamp to valid grid range (1 cell margin)
                    if (ix < 1) ix = 1;
                    if (ix >= sim->nx - 1) ix = sim->nx - 2;
                    if (iy < 1) iy = 1;
                    if (iy >= sim->ny - 1) iy = sim->ny - 2;

                    sim->sources[dragged_source_idx].ix = ix;
                    sim->sources[dragged_source_idx].iy = iy;

                    if (dragged_source_idx < wizard.cfg.source_count) {
                        double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
                        double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
                        wizard.cfg.source_configs[dragged_source_idx].x =
                            (double)ix / nx1;
                        wizard.cfg.source_configs[dragged_source_idx].y =
                            (double)iy / ny1;
                    }
                }
                // Then handle paint drag
                else if (paint_mode &&
                         (e.motion.state & SDL_BUTTON_LMASK)) {
                    if (!app.viewport_valid) {
                        continue;
                    }
                    int mx = e.motion.x;
                    int my = e.motion.y;
                    float local_x = (float)mx - app.viewport_pos.x;
                    float local_y = (float)my - app.viewport_pos.y;
                    if (local_x < 0.0f || local_y < 0.0f) {
                        continue;
                    }
                    if (local_x >= app.viewport_size.x || local_y >= app.viewport_size.y) {
                        continue;
                    }
                    int ix = (int)(local_x / (float)scale);
                    int iy = (int)(local_y / (float)scale);
                    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) {
                        continue;
                    }
                    if (paint_material_id >= 0) {
                        apply_paint_with_material(sim, ix, iy, paint_brush_size, paint_material_id);
                    } else {
                        apply_selected_paint(sim,
                                             ix,
                                             iy,
                                             paint_brush_size,
                                             paint_material_type,
                                             paint_epsilon);
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                if (dragging_source) {
                    if (dragged_source_idx >= 0 && dragged_source_idx < MAX_SRC) {
                        ui_log_add(&app,
                                   "Source %d moved to (%d, %d)",
                                   dragged_source_idx,
                                   sim->sources[dragged_source_idx].ix,
                                   sim->sources[dragged_source_idx].iy);
                    }
                    dragging_source = false;
                    dragged_source_idx = -1;
                }
            }
        }

        // Tie wizard advanced flag to app-level basic/advanced mode
        wizard.advanced = !app.basic_mode;

        ImGui_ImplSDL2_NewFrame();
        ImGui_ImplSDLRenderer2_NewFrame();
        ImGui::NewFrame();

        // Main menu bar
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Save Configuration...", "Ctrl+S")) {
                    sync_sim_to_wizard_config(sim, &wizard);
                    const char* save_path = "configs/saved_config.json";
                    if (save_config_json(save_path, &wizard.cfg)) {
                        ui_log_add(&app, "Configuration saved: %s", save_path);
                    } else {
                        ui_log_add(&app, "Failed to save configuration");
                    }
                }

                ImGui::Separator();

                if (ImGui::MenuItem("Load Waveguide Preset", "F5")) {
                    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                    char errbuf[256];
                    if (config_loader_parse_file("configs/waveguide.json", &cfg, errbuf,
                                                 sizeof(errbuf))) {
                        if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                            wizard_init_from_config(wizard, &cfg);

                            // Resize window to match new grid
                            width = sim->nx * scale;
                            height = sim->ny * scale;
                            if (width < 1920) width = 1920;
                            if (height < 1080) height = 1080;
                            SDL_SetWindowSize(render->window, width, height);

                            // Reset visualization scaling
                            vmax_smooth = sim->step_Ez_absmax;
                            if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                            scope_vmax = 1.0;
                            auto_rescale = true;
                            hold_color = false;
                            hold_scope = false;

                            ui_log_add(&app, "Loaded: waveguide.json");
                        } else {
                            ui_log_add(&app, "Failed to load: waveguide.json");
                        }
                    } else {
                        ui_log_add(&app, "Failed to parse: waveguide.json");
                    }
                }

                if (ImGui::MenuItem("Load CPW Filter Preset", "F6")) {
                    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                    char errbuf[256];
                    if (config_loader_parse_file("configs/cpw_filter.json", &cfg, errbuf,
                                                 sizeof(errbuf))) {
                        if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                            wizard_init_from_config(wizard, &cfg);

                            // Resize window to match new grid
                            width = sim->nx * scale;
                            height = sim->ny * scale;
                            if (width < 1920) width = 1920;
                            if (height < 1080) height = 1080;
                            SDL_SetWindowSize(render->window, width, height);

                            // Reset visualization scaling
                            vmax_smooth = sim->step_Ez_absmax;
                            if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                            scope_vmax = 1.0;
                            auto_rescale = true;
                            hold_color = false;
                            hold_scope = false;

                            ui_log_add(&app, "Loaded: cpw_filter.json");
                        } else {
                            ui_log_add(&app, "Failed to load: cpw_filter.json");
                        }
                    } else {
                        ui_log_add(&app, "Failed to parse: cpw_filter.json");
                    }
                }

                ImGui::Separator();

        if (ImGui::MenuItem("Quit", "ESC")) {
            running = false;
        }

        ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Materials")) {
        if (ImGui::MenuItem("Material Browser...", "M")) {
            app.material_browser_open = true;
        }
        ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Analysis")) {
        if (ImGui::MenuItem("S-Parameters...")) {
            app.sparam_window_open = true;
        }
        if (ImGui::MenuItem("Smith Chart...")) {
            app.smith_chart_open = true;
        }
        ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("View")) {
        ImGui::TextUnformatted("Visualization");
        ImGui::Separator();
        if (ImGui::RadioButton("Field (Ez)", app.visualization_mode == 0)) {
            app.visualization_mode = 0;
        }
        if (ImGui::RadioButton("Material Distribution", app.visualization_mode == 1)) {
            app.visualization_mode = 1;
        }
        if (ImGui::RadioButton("Field + Material Overlay", app.visualization_mode == 2)) {
            app.visualization_mode = 2;
        }
        if (app.visualization_mode == 2) {
            ImGui::SliderFloat("Overlay Alpha", &app.material_overlay_alpha, 0.0f, 1.0f);
        }
        ImGui::Checkbox("Material Outlines", &app.show_material_outlines);
        ImGui::Checkbox("Material Legend", &app.show_material_legend);

        ImGui::Separator();
        ImGui::TextUnformatted("Panels");
        ImGui::Separator();
        ImGui::MenuItem("Grid Settings", nullptr, &app.show_grid_panel);
        ImGui::MenuItem("Scene", nullptr, &app.show_scene_panel);
        ImGui::MenuItem("Sources", nullptr, &app.show_sources_panel);
        ImGui::MenuItem("Blocks", nullptr, &app.show_blocks_panel);
        ImGui::MenuItem("Probes", nullptr, &app.show_probes_panel);
        ImGui::MenuItem("Run Controls", nullptr, &app.show_run_panel);
        ImGui::MenuItem("Run Settings", nullptr, &app.show_run_settings_panel);
        ImGui::Checkbox("Highlight blocks for selected material", &app.highlight_blocks_by_material);
        ImGui::Checkbox("Auto-filter blocks on material select", &app.auto_filter_blocks_on_select);
        ImGui::EndMenu();
    }

            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("Keyboard Shortcuts", "F1")) {
                    show_help_overlay = !show_help_overlay;
                }
                ImGui::MenuItem("About", nullptr, false, false);
                ImGui::EndMenu();
            }

            ImGui::EndMainMenuBar();
        }

        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos);
        ImGui::SetNextWindowSize(viewport->Size);
    ImGuiWindowFlags root_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                  ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                                  ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus |
                                  ImGuiWindowFlags_NoBackground;
        ImGui::Begin("RootLayout", nullptr, root_flags);

        ImVec2 full = ImGui::GetContentRegionAvail();
        // Responsive proportions tuned for the 1920x1080 studio layout baseline
        const float target_left_w = 300.0f;
        const float target_right_w = 360.0f;
        const float target_bottom_h = 200.0f;
        const float min_center_w = 640.0f;
        const float min_main_h = 420.0f;

        float width_needed = target_left_w + target_right_w + min_center_w;
        float width_scale = (width_needed > 0.0f && full.x < width_needed)
                                ? (full.x / width_needed)
                                : 1.0f;
        float left_w = target_left_w * width_scale;
        float right_w = target_right_w * width_scale;
        if (left_w + right_w > full.x) {
            float squeeze = (left_w + right_w > 0.0f) ? (full.x / (left_w + right_w)) : 1.0f;
            left_w *= squeeze;
            right_w *= squeeze;
        }
        float center_w = ImMax(0.0f, full.x - left_w - right_w);

        float height_needed = target_bottom_h + min_main_h;
        float height_scale = (height_needed > 0.0f && full.y < height_needed)
                                 ? (full.y / height_needed)
                                 : 1.0f;
        float bottom_h = target_bottom_h * height_scale;
        float main_h = ImMax(0.0f, full.y - bottom_h);

        ImVec2 origin = ImGui::GetCursorPos();
        const ImVec4 column_bg = ImVec4(0.24f, 0.24f, 0.24f, 0.95f);
        const ImVec4 bottom_bg = ImVec4(0.15f, 0.15f, 0.15f, 0.95f);

        // Left column
        ImGui::SetCursorPos(origin);
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.10f, 0.10f, 0.11f, 1.00f));
        ImGui::BeginChild("LeftColumn", ImVec2(left_w, main_h), true);
        bool first_left_panel = true;
        auto left_spacing = [&]() {
            if (!first_left_panel) {
                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
            }
            first_left_panel = false;
        };
        if (app.show_grid_panel) {
            left_spacing();
            draw_grid_panel(&wizard, &app);
        }
        if (app.show_scene_panel) {
            left_spacing();
            draw_scene_panel(sim, wizard);
        }
        if (app.show_sources_panel) {
            left_spacing();
            draw_sources_panel(sim, wizard, &app);
        }
        if (app.show_blocks_panel) {
            left_spacing();
            draw_blocks_panel(wizard, &bootstrap, sim, &app);
        }
        if (app.show_probes_panel) {
            left_spacing();
            draw_probes_panel(sim);
        }
        if (app.show_material_legend) {
            left_spacing();
            draw_material_legend(&app, sim, &wizard, &bootstrap);
        }
        ImGui::EndChild();
        ImGui::PopStyleColor();

        // Center viewport placeholder - keep border but let SDL field stay visible and interactive
        ImGui::SetCursorPos(ImVec2(origin.x + left_w, origin.y));
        ImVec2 center_screen_pos = ImGui::GetCursorScreenPos();
        ImVec2 center_child_size(ImMax(1.0f, center_w), ImMax(1.0f, main_h));
        ImGuiWindowFlags center_flags =
            ImGuiWindowFlags_NoBackground |
            ImGuiWindowFlags_NoScrollbar |
            ImGuiWindowFlags_NoScrollWithMouse |
            ImGuiWindowFlags_NoNavInputs |
            ImGuiWindowFlags_NoNavFocus |
            ImGuiWindowFlags_NoInputs;
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0, 0, 0, 0));
        ImGui::BeginChild("CenterViewport",
                          center_child_size,
                          true,
                          center_flags);
        ImGui::EndChild();
        ImGui::PopStyleColor();
        if (sim) {
            float field_px_w = (float)(sim->nx * scale);
            float field_px_h = (float)(sim->ny * scale);
            float draw_origin_x = center_screen_pos.x;
            float draw_origin_y = center_screen_pos.y;
            float draw_w = center_child_size.x;
            float draw_h = center_child_size.y;
            if (field_px_w <= center_child_size.x) {
                draw_origin_x += (center_child_size.x - field_px_w) * 0.5f;
                draw_w = field_px_w;
            } else {
                draw_w = center_child_size.x;
            }
            if (field_px_h <= center_child_size.y) {
                draw_origin_y += (center_child_size.y - field_px_h) * 0.5f;
                draw_h = field_px_h;
            } else {
                draw_h = center_child_size.y;
            }
            if (draw_w > 4.0f && draw_h > 4.0f) {
                app.viewport_pos = ImVec2(draw_origin_x, draw_origin_y);
                app.viewport_size = ImVec2(draw_w, draw_h);
                app.viewport_valid = true;
            } else {
                app.viewport_valid = false;
            }
        } else {
            app.viewport_valid = false;
        }

        // Right column (Simulation Controls + Wizard)
        ImGui::SetCursorPos(ImVec2(origin.x + left_w + center_w, origin.y));
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.10f, 0.10f, 0.11f, 1.00f));
        ImGui::BeginChild("RightColumn", ImVec2(right_w, main_h), true);

        // Simulation Controls
        if (app.show_run_panel) {
            ImGui::SetNextItemOpen(false, ImGuiCond_Once);
            if (ImGui::CollapsingHeader("Simulation Controls")) {
            // Control buttons
            ImVec2 button_size(right_w * 0.45f, 0.0f);
            if (ImGui::Button(paused ? "Resume" : "Pause", button_size)) {
                paused = !paused;
            }
            ImGui::SameLine();
            ImGui::Text("Space");

            ImGui::Separator();

            // Frequency slider
            ImGui::TextUnformatted("Frequency");
            ImGui::Indent();
            float freq_ghz = (float)(sim->freq * 1e-9);
            if (ImGui::SliderFloat("GHz", &freq_ghz, 0.001f, 5.0f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                double new_freq = (double)freq_ghz * 1e9;
                if (new_freq > 0.0) {
                    fdtd_update_grid_for_freq(sim, new_freq);
                    sources_set_freq(sim->sources, new_freq);
                }
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Adjust simulation frequency (1 MHz - 5 GHz)");
            }
            ImGui::Unindent();

            ImGui::Separator();

            // Steps/frame slider
            ImGui::TextUnformatted("Speed");
            ImGui::Indent();
            if (ImGui::SliderInt("steps/frame", &steps_per_frame, 1, 50)) {
                // Speed is controlled per-frame in the main loop
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Number of simulation steps per rendered frame");
            }
            ImGui::Unindent();

            ImGui::Separator();
            ImGui::TextUnformatted("Paint Mode");
            ImGui::Indent();
            ImGui::Text("Active: %s (M/U)", paint_mode ? "YES" : "NO");
            if (paint_mode) {
                if (paint_material_id >= 0) {
                    const Material* mat = material_library_get_by_id(paint_material_id);
                    if (mat) {
                        ImVec4 color = ImVec4(mat->color_r / 255.0f,
                                              mat->color_g / 255.0f,
                                              mat->color_b / 255.0f,
                                              1.0f);
                        ImGui::PushStyleColor(ImGuiCol_Text, color);
                        ImGui::Text("Material: %s", mat->name);
                        ImGui::PopStyleColor();
                        if (mat->type == MAT_TYPE_PEC) {
                            ImGui::TextDisabled("(Perfect conductor)");
                        } else if (mat->type == MAT_TYPE_PMC) {
                            ImGui::TextDisabled("(Magnetic conductor)");
                        } else {
                            ImGui::Text("εᵣ = %.2f", mat->epsilon_r);
                            if (mat->tan_delta > 0.0) {
                                ImGui::Text("tan(δ) = %.4f", mat->tan_delta);
                            }
                        }
                        if (ImGui::SmallButton("Change Material...")) {
                            app.material_browser_open = true;
                        }
                    }
                } else {
                    const char* types[] = {"PEC", "PMC", "Dielectric"};
                    const int type_idx = (paint_material_type >= 0 && paint_material_type < 3)
                                             ? paint_material_type
                                             : 0;
                    ImGui::Text("Type: %s (I)", types[type_idx]);
                    if (type_idx == 2) {
                        ImGui::Text("Epsilon: %.1f (O/P)", paint_epsilon);
                    }
                    ImGui::TextDisabled("Use Materials menu for library");
                }
                ImGui::Text("Brush radius: %d", paint_brush_size);
            }
            ImGui::Unindent();

            ImGui::Separator();

            // Status section
            ImGui::TextUnformatted("Status");
            ImGui::Indent();
            ImGui::Text("Step: %d", sim->timestep);
            ImGui::Text("FPS: %.1f", fps_avg);
            ImGui::Unindent();

            if (ImGui::TreeNode("Grid Info")) {
                ImGui::Text("Size: %d x %d", sim->nx, sim->ny);
                ImGui::Text("Domain: %.3f x %.3f m", sim->lx, sim->ly);
                ImGui::Text("dt: %.3e s", sim->dt);
                ImGui::Text("freq: %.3f GHz", sim->freq * 1e-9);
                ImGui::TreePop();
            }

            ImGui::Separator();
            ImGui::TextUnformatted("Visualization Scale");
            ImGui::Indent();
            const char* mode_str = "Manual (L)";
            if (auto_rescale && !hold_color && !hold_scope) {
                mode_str = "Auto (A)";
            } else if (hold_color) {
                mode_str = "Hold Color (H)";
            } else if (hold_scope) {
                mode_str = "Hold Scope (J)";
            }
            ImGui::Text("Mode: %s", mode_str);
            ImGui::Text("Field vmax: %.3e", vmax_smooth);
            ImGui::Text("Scope vmax: %.3e", scope_vmax);
            ImGui::Unindent();

            ImGui::Separator();
            ImGui::Checkbox("Basic mode", &app.basic_mode);

            ImGui::Separator();
            if (ImGui::Button(app.show_grid_overlay ? "Hide grid" : "Show grid")) {
                app.show_grid_overlay = !app.show_grid_overlay;
            }
            }
        }
        if (app.show_run_settings_panel) {
            ImGui::Spacing();
            draw_run_settings_panel(&wizard, &app);
        }

        // Interactive Tools section
        ImGui::SetNextItemOpen(false, ImGuiCond_Once);
        if (ImGui::CollapsingHeader("Interactive Tools")) {
            ImGui::TextUnformatted("Source Placement");
            ImGui::Indent();
            int max_src = MAX_SRC;
            if (app.selected_source < 0) app.selected_source = 0;
            if (app.selected_source >= max_src) app.selected_source = max_src - 1;
            char current_src_label[64];
            const Source& ss = sim->sources[app.selected_source];
            std::snprintf(current_src_label, sizeof(current_src_label),
                          "#%d (%s)", app.selected_source,
                          source_type_label(ss.type));
            if (ImGui::BeginCombo("Source", current_src_label)) {
                for (int i = 0; i < max_src; ++i) {
                    const Source& s = sim->sources[i];
                    char label[64];
                    std::snprintf(label, sizeof(label), "#%d (%s, %s)", i,
                                  source_type_label(s.type),
                                  s.active ? "on" : "off");
                    bool selected = (i == app.selected_source);
                    if (ImGui::Selectable(label, selected)) {
                        app.selected_source = i;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            ImGui::Checkbox("Click to move", &app.placing_source);
            ImGui::Unindent();

            ImGui::Separator();
            ImGui::TextUnformatted("Block Drawing");
            ImGui::Indent();
            ImGui::InputInt("Block index", &app.selected_block);
            if (app.selected_block < 0) app.selected_block = 0;
            if (app.selected_block >= CONFIG_MAX_MATERIAL_RECTS) {
                app.selected_block = CONFIG_MAX_MATERIAL_RECTS - 1;
            }
            ImGui::Checkbox("Draw with 2 clicks", &app.placing_block);
            if (app.placing_block) {
                if (!app.block_first_set) {
                    ImGui::TextWrapped("Click first corner, then opposite corner.");
                } else {
                    ImGui::Text("First: (%d,%d)", app.block_first_i, app.block_first_j);
                }
            }
            ImGui::Unindent();

            if (app.last_click_i >= 0 && app.last_click_j >= 0) {
                ImGui::Separator();
                ImGui::Text("Last click: (%d,%d)", app.last_click_i, app.last_click_j);
            }
        }

        // Handle deferred restarts (from panels)
        if (app.request_rebootstrap) {
            const char* restart_msg =
                (app.rebootstrap_message[0] != '\0') ? app.rebootstrap_message
                                                     : "Rebooted simulation";
            if (rebootstrap_simulation(&wizard.cfg, &bootstrap, &sim, &scope, scale)) {
                width = sim->nx * scale;
                height = sim->ny * scale;
                if (width < 1920) width = 1920;
                if (height < 1080) height = 1080;
                SDL_SetWindowSize(render->window, width, height);
                ui_log_add(&app, "%s", restart_msg);
            } else {
                ui_log_add(&app, "Failed to rebootstrap simulation");
            }
            app.request_rebootstrap = false;
            app.rebootstrap_message[0] = '\0';
        }

        ImGui::EndChild();
        ImGui::PopStyleColor();

        // Bottom strip with tabs for Scope and Log
        ImGui::SetCursorPos(ImVec2(origin.x, origin.y + main_h));
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.08f, 0.08f, 0.09f, 1.00f));
        ImGui::BeginChild("BottomStrip", ImVec2(full.x, bottom_h), true);
            if (ImGui::BeginTabBar("BottomTabs", ImGuiTabBarFlags_None)) {
                if (ImGui::BeginTabItem("Scope")) {
                    if (scope.y && scope.n > 0) {
                        const int N = scope.n;
                        static float values[1024];
                    int max_samples = (N <= (int)(sizeof(values) / sizeof(values[0]))) ? N : (int)(sizeof(values) / sizeof(values[0]));
                    double vmin = 0.0, vmax_scope = 0.0;
                    bool first = true;
                    for (int i = 0; i < max_samples; ++i) {
                        int idx = (scope.head + i) % N;
                        double v = scope.y[idx];
                        values[i] = (float)v;
                        if (first) {
                            vmin = vmax_scope = v;
                            first = false;
                        } else {
                            if (v < vmin) vmin = v;
                            if (v > vmax_scope) vmax_scope = v;
                        }
                    }

                    float ymin = 0.0f;
                    float ymax = 0.0f;
                    if (scope_vmax > 0.0) {
                        ymax = (float)scope_vmax;
                        ymin = -ymax;
                    } else {
                        ymin = (float)vmin;
                        ymax = (float)vmax_scope;
                        if (ymin == ymax) {
                            ymin -= 1.0f;
                            ymax += 1.0f;
                        }
                    }

                    ImVec2 plot_size = ImVec2(ImGui::GetContentRegionAvail().x, bottom_h - 60.0f);
                    ImGui::PlotLines("Probe Ez(center)", values, max_samples, 0, nullptr, ymin, ymax, plot_size);
                } else {
                    ImGui::TextUnformatted("Scope: no data available");
                }
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("FFT")) {
                    if (scope.n > 0 && sim) {
                        static bool use_db_scale = true;
                        ImGui::Checkbox("dB Scale", &use_db_scale);
                        static double fft_freq[1024];
                        static double fft_mag[1024];
                        static double fft_mag_db[1024];
                        static double fft_phase[1024];
                        int fft_points =
                            compute_scope_fft(&scope, sim->dt, fft_freq, fft_mag, fft_phase,
                                              IM_ARRAYSIZE(fft_freq));
                        if (fft_points > 1) {
                            int plot_count = fft_points - 1;
                            const double eps = 1e-12;
                            if (use_db_scale) {
                                for (int i = 0; i < fft_points; ++i) {
                                    fft_mag_db[i] = 20.0 * std::log10(fft_mag[i] + eps);
                                }
                            }
                            const double* mag_view = use_db_scale ? fft_mag_db : fft_mag;
                            const char* mag_label = use_db_scale ? "FFT (dB)" : "FFT";
                            const double* freq_view = fft_freq + 1;
                            const double* mag_plot = mag_view + 1;
                            const double* phase_plot = fft_phase + 1;
                            if (plot_count > 0 &&
                                ImPlot::BeginPlot("FFT Magnitude", ImVec2(-1, 300))) {
                                ImPlot::SetupAxis(ImAxis_X1, "Frequency (Hz)");
                                ImPlot::SetupAxis(ImAxis_Y1,
                                                  use_db_scale ? "Magnitude (dB)"
                                                               : "Magnitude");
                                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                                ImPlot::PlotLine(mag_label, freq_view, mag_plot, plot_count);

                                int peak_idx = 1;
                                double peak_val = mag_view[1];
                                for (int i = 2; i < fft_points; ++i) {
                                    if (mag_view[i] > peak_val) {
                                        peak_val = mag_view[i];
                                        peak_idx = i;
                                    }
                                }
                                double peak_freq = fft_freq[peak_idx];
                                ImPlot::Annotation(peak_freq,
                                                   peak_val,
                                                   ImVec4(1, 1, 0, 1),
                                                   ImVec2(10, -10),
                                                   true,
                                                   "Peak: %.2f GHz",
                                                   peak_freq * 1e-9);

                                double sim_freq_marker = sim->freq;
                                ImPlot::DragLineX(1, &sim_freq_marker, ImVec4(1, 0, 0, 1));

                                ImPlot::EndPlot();
                            }

                            if (plot_count > 0 &&
                                ImPlot::BeginPlot("FFT Phase", ImVec2(-1, 200))) {
                                ImPlot::SetupAxis(ImAxis_X1, "Frequency (Hz)");
                                ImPlot::SetupAxis(ImAxis_Y1, "Phase (deg)");
                                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                                ImPlot::PlotLine("Phase", freq_view, phase_plot, plot_count);
                                ImPlot::EndPlot();
                            }

                            if (ImGui::Button("Export FFT to CSV")) {
                                FILE* fp = fopen("fft_export.csv", "w");
                                if (fp) {
                                    fprintf(fp, "# FFT Export\n");
                                    fprintf(fp, "# Center frequency: %.12e Hz\n", sim->freq);
                                    fprintf(fp,
                                            "Frequency(Hz),Magnitude,Magnitude(dB),Phase(deg)\n");
                                    for (int i = 0; i < fft_points; ++i) {
                                        double mag_db = 20.0 * std::log10(fft_mag[i] + eps);
                                        fprintf(fp,
                                                "%.12e,%.12e,%.12e,%.12e\n",
                                                fft_freq[i],
                                                fft_mag[i],
                                                mag_db,
                                                fft_phase[i]);
                                    }
                                    fclose(fp);
                                    ui_log_add(&app, "FFT exported: fft_export.csv");
                                } else {
                                    ui_log_add(&app, "Failed to export FFT");
                                }
                            }
                        } else {
                            ImGui::TextUnformatted("FFT buffer too small");
                        }
                    } else {
                        ImGui::TextUnformatted("No FFT data available");
                    }
                    ImGui::EndTabItem();
                }
            if (ImGui::BeginTabItem("Log")) {
                ImGui::BeginChild("LogScrollRegion", ImVec2(0, 0), false);
                for (int i = 0; i < app.log_count; ++i) {
                    ImGui::TextUnformatted(app.log_lines[i]);
                }
                if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
                    ImGui::SetScrollHereY(1.0f);
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::EndChild();
        ImGui::PopStyleColor();

        ImGui::End();

        // Help overlay (F1)
        if (show_help_overlay) {
            ImGui::SetNextWindowSize(ImVec2(600.0f, 400.0f), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
                                    ImGuiCond_Appearing,
                                    ImVec2(0.5f, 0.5f));

            if (ImGui::Begin("Keyboard Shortcuts & Help",
                             &show_help_overlay,
                             ImGuiWindowFlags_NoCollapse)) {
                if (ImGui::BeginTabBar("HelpTabs")) {
                    if (ImGui::BeginTabItem("Shortcuts")) {
                        ImGui::TextUnformatted("File Operations");
                        ImGui::Separator();
                        ImGui::Columns(2, "file_ops_cols", false);
                        ImGui::SetColumnWidth(0, 140.0f);
                        ImGui::TextUnformatted("Ctrl+S");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Save configuration");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F2");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Save screenshot");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F3");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Export scope FFT");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F5");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Load waveguide preset");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F6");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Load CPW filter preset");
                        ImGui::Columns(1);

                        ImGui::Separator();
                        ImGui::TextWrapped("Workflow:");
                ImGui::BulletText(
                    "Direct manipulation: Paint materials (M/U/I/O/P), drag sources with mouse");
                ImGui::BulletText(
                    "Presets: F5 (waveguide), F6 (CPW filter)");
                ImGui::BulletText(
                    "Save/Load: Ctrl+S saves to configs/saved_config.json");
                ImGui::BulletText(
                    "Materials menu \xE2\x86\x92 Material Browser to select materials");
                ImGui::BulletText(
                    "Numpad 1-8: Quick-select materials");
                ImGui::BulletText(
                    "I/O/P keys only in legacy mode (no material selected)");

                        ImGui::Separator();
                        ImGui::TextUnformatted("Press F1 again to close this window.");
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Materials & Visualization")) {
                        ImGui::TextUnformatted("Material Library:");
                        ImGui::BulletText("Materials menu \xE2\x86\x92 Material Browser");
                        ImGui::BulletText("Select from 11 curated EM materials");
                        ImGui::BulletText("4 Metals: Copper, Gold, Silver, Aluminum");
                        ImGui::BulletText("7 Dielectrics: Air, FR4, Rogers, Teflon, Silicon, Glass, Alumina");
                        ImGui::BulletText("Each material has accurate \xCE\xB5\xE1\xB5\xBD, tan(δ), σ values");
                        ImGui::BulletText("Numpad 1-8: Quick-select materials");
                        ImGui::BulletText("I/O/P: Legacy manual mode (when no material selected)");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("Material Visualization:");
                        ImGui::BulletText("E key: Cycle view modes (Field \xE2\x86\x92 Material \xE2\x86\x92 Overlay)");
                        ImGui::BulletText("View menu \xE2\x86\x92 toggles for overlay alpha, outlines, legend");
                        ImGui::BulletText("Material legend lists library materials with colors");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("View Modes:");
                        ImGui::BulletText("Field: Standard Ez field visualization");
                        ImGui::BulletText("Material: Show \xCE\xB5 distribution with material colors");
                        ImGui::BulletText("Overlay: Field + semi-transparent material layer");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("S-Parameters:");
                        ImGui::BulletText("Analysis menu \xE2\x86\x92 S-Parameters window");
                        ImGui::BulletText("Plots show S21, S11, and VSWR vs frequency");
                        ImGui::BulletText("-3 dB bandwidth markers derived from peak S21");
                        ImGui::BulletText("Export enhanced CSV with VSWR and annotations");
                        ImGui::BulletText("Smith Chart window visualizes reflection coefficients");
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
        }

        app.paint_material_id = paint_material_id;
        draw_material_browser(&app, &paint_material_id);
        app.paint_material_id = paint_material_id;

        // Overlay: show source IDs and wizard blocks on top of the field
        {
            ImDrawList* dl = ImGui::GetForegroundDrawList();
            const ImU32 src_col = IM_COL32(255, 255, 255, 230);
            const ImU32 block_col = IM_COL32(255, 255, 0, 160);
            const ImU32 block_sel_col = IM_COL32(0, 255, 128, 200);

            // Source labels
            if (app.viewport_valid) {
                for (int k = 0; k < MAX_SRC; ++k) {
                    const Source& s = sim->sources[k];
                    if (!s.active) continue;
                    float sx = app.viewport_pos.x + ((float)s.ix + 0.5f) * (float)scale;
                    float sy = app.viewport_pos.y + ((float)s.iy + 0.5f) * (float)scale;
                    char label[16];
                    std::snprintf(label, sizeof(label), "%d", k);
                    dl->AddText(ImVec2(sx + 6.0f, sy - 6.0f), src_col, label);

                    // Highlight if being dragged
                    if (dragging_source && k == dragged_source_idx) {
                        float radius = 10.0f;
                        dl->AddCircle(ImVec2(sx, sy),
                                      radius,
                                      IM_COL32(255, 200, 0, 200),
                                      16,
                                      3.0f);
                    }
                }

                // Block outlines from wizard config
                int rect_count = wizard.cfg.material_rect_count;
                if (rect_count < 0) rect_count = 0;
                if (rect_count > CONFIG_MAX_MATERIAL_RECTS) rect_count = CONFIG_MAX_MATERIAL_RECTS;
                const Material* overlay_sel_mat = nullptr;
                if (app.selected_material_id >= 0) {
                    overlay_sel_mat = material_library_get_by_id(app.selected_material_id);
                }
                for (int i = 0; i < rect_count; ++i) {
                    const MaterialRectSpec& r = wizard.cfg.material_rects[i];
                    float x0 = app.viewport_pos.x + (float)(r.x0 * (double)sim->nx * (double)scale);
                    float x1 = app.viewport_pos.x + (float)(r.x1 * (double)sim->nx * (double)scale);
                    float y0 = app.viewport_pos.y + (float)(r.y0 * (double)sim->ny * (double)scale);
                    float y1 = app.viewport_pos.y + (float)(r.y1 * (double)sim->ny * (double)scale);
                    ImVec2 p0(x0, y0);
                    ImVec2 p1(x1, y1);
                    bool mat_match = overlay_sel_mat && rect_matches_material(r, overlay_sel_mat);
                    ImU32 col = block_col;
                    if (i == app.selected_block) {
                        col = block_sel_col;
                    } else if (mat_match) {
                        col = IM_COL32(0, 200, 255, 200);
                    }
                    dl->AddRect(p0, p1, col, 0.0f, 0, 1.5f);

                    char blabel[16];
                    std::snprintf(blabel, sizeof(blabel), "B%d", i);
                    dl->AddText(ImVec2(x0 + 4.0f, y0 + 4.0f), col, blabel);
                }

                // Paint cursor indicator
                if (paint_mode) {
                    ImVec2 mouse_pos = io.MousePos;
                    float local_x = mouse_pos.x - app.viewport_pos.x;
                    float local_y = mouse_pos.y - app.viewport_pos.y;
                    if (local_x >= 0.0f && local_y >= 0.0f &&
                        local_x < app.viewport_size.x && local_y < app.viewport_size.y) {
                        int ix = (int)(local_x / (float)scale);
                        int iy = (int)(local_y / (float)scale);
                        if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
                            float cx = app.viewport_pos.x + ((float)ix + 0.5f) * (float)scale;
                            float cy = app.viewport_pos.y + ((float)iy + 0.5f) * (float)scale;
                            float radius = (float)paint_brush_size * (float)scale;
                            ImU32 col = IM_COL32(255, 255, 0, 180);
                            if (paint_material_id >= 0) {
                                const Material* mat = material_library_get_by_id(paint_material_id);
                                if (mat) {
                                    col = IM_COL32(mat->color_r, mat->color_g, mat->color_b, 200);
                                }
                            } else {
                                if (paint_material_type == 0) {
                                    col = IM_COL32(255, 80, 80, 200);
                                } else if (paint_material_type == 1) {
                                    col = IM_COL32(80, 160, 255, 200);
                                } else {
                                    col = IM_COL32(120, 255, 120, 200);
                                }
                            }
                            dl->AddCircle(ImVec2(cx, cy), radius, col, 32, 1.5f);
                        }
                    }
                }
            }
        }

        if (!paused) {
            for (int s = 0; s < steps_per_frame; ++s) {
                fdtd_step(sim);
                int px = sim->nx / 2;
                int py = sim->ny / 2;
                double probe_val = fdtd_get_Ez(sim, px, py);
                scope_push(&scope, probe_val);
            }
        }

        Uint64 now = SDL_GetPerformanceCounter();
        double dt = (double)(now - prev) / (double)perf_freq;
        prev = now;
        double fps_inst = (dt > 0.0) ? (1.0 / dt) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        // Update visualization metrics (field and scope scaling)
        {
            // Field amplitude (from last simulation step)
            double field_absmax = sim->step_Ez_absmax;
            if (field_absmax < 1e-12) field_absmax = 1e-12;

            bool field_autop = auto_rescale || hold_scope;
            if (!hold_color && field_autop) {
                // Smooth towards recent maximum
                const double alpha = 0.1;
                if (vmax_smooth <= 0.0) {
                    vmax_smooth = field_absmax;
                } else {
                    vmax_smooth = (1.0 - alpha) * vmax_smooth + alpha * field_absmax;
                }
            } else if (hold_color) {
                // Only allow increases when holding color scale
                if (field_absmax > vmax_smooth) {
                    vmax_smooth = field_absmax;
                }
            }

            // Scope amplitude (use absolute max over buffer)
            double scope_absmax = 0.0;
            if (scope.y && scope.n > 0) {
                bool first = true;
                for (int i = 0; i < scope.n; ++i) {
                    double v = std::fabs(scope.y[i]);
                    if (first || v > scope_absmax) {
                        scope_absmax = v;
                        first = false;
                    }
                }
            }

            bool scope_autop = auto_rescale || hold_color;
            if (!hold_scope && scope_autop && scope_absmax > 0.0) {
                const double alpha_scope = 0.1;
                if (scope_vmax <= 0.0) {
                    scope_vmax = scope_absmax;
                } else {
                    scope_vmax =
                        (1.0 - alpha_scope) * scope_vmax + alpha_scope * scope_absmax;
                }
            } else if (hold_scope && scope_absmax > scope_vmax) {
                scope_vmax = scope_absmax;
            }
        }

        SDL_SetRenderDrawColor(render->renderer,
                               viewport_clear_color.r,
                               viewport_clear_color.g,
                               viewport_clear_color.b,
                               viewport_clear_color.a);
        SDL_RenderClear(render->renderer);

        SDL_Rect sdl_viewport = {0, 0, 0, 0};
        bool viewport_set = false;
        if (app.viewport_valid) {
            sdl_viewport.x = (int)app.viewport_pos.x;
            sdl_viewport.y = (int)app.viewport_pos.y;
            sdl_viewport.w = (int)app.viewport_size.x;
            sdl_viewport.h = (int)app.viewport_size.y;
            if (sdl_viewport.w > 0 && sdl_viewport.h > 0) {
                SDL_RenderSetViewport(render->renderer, &sdl_viewport);
                viewport_set = true;
            }
        }
        if (!viewport_set) {
            SDL_RenderSetViewport(render->renderer, NULL);
        }

        double vmax = vmax_smooth;
        if (vmax <= 0.0) vmax = 1.0;

        switch (app.visualization_mode) {
            case 0:
                render_field_heatmap(render, sim, vmax, 1.0);
                break;
            case 1:
                render_material_distribution(render->renderer, sim, scale);
                break;
            case 2:
            default:
                render_field_heatmap(render, sim, vmax, 1.0);
                SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_BLEND);
                render_material_overlay(render->renderer,
                                        sim,
                                        scale,
                                        app.material_overlay_alpha);
                SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_NONE);
                break;
        }

        if (app.show_material_outlines) {
            SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_BLEND);
            render_material_outlines(render->renderer, sim, scale);
            SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_NONE);
        }

        if (app.show_grid_overlay) {
            SDL_Color grid_col = {40, 40, 50, 255};
            render_grid_overlay(render, sim, grid_col);
        }
        render_sources(render, sim->sources);
        if (viewport_set) {
            SDL_RenderSetViewport(render->renderer, NULL);
        }

        draw_sparameter_window(&app);
        draw_smith_chart(&app);

        ImGui::Render();
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), render->renderer);
        SDL_RenderPresent(render->renderer);

        SDL_Delay(UI_DELAY_MS);
    }

    scope_free(&scope);
    ImGui_ImplSDLRenderer2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    render_free(render);
    material_library_shutdown();
    simulation_bootstrap_shutdown(&bootstrap);

    return 0;
}
