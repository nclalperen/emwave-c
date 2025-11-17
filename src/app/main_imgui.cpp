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

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

#include <cstdio>
#include <cstring>

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
    bool show_wizard_panel;
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
};

struct WizardState {
    SimulationConfig cfg;
    bool open;
    bool advanced;
};

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

/* Scene overview panel ---------------------------------------------------- */
static void draw_scene_panel(const SimulationState* sim, const WizardState& wizard) {
    if (!ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen)) return;

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

/* Sources panel ----------------------------------------------------------- */
static void draw_sources_panel(SimulationState* sim, WizardState& wizard, AppState* app) {
    if (!sim || !app) return;
    if (!ImGui::CollapsingHeader("Sources", ImGuiTreeNodeFlags_DefaultOpen)) return;

    ImGui::Indent();
    int max_src = MAX_SRC;
    if (wizard.cfg.source_count > max_src) {
        wizard.cfg.source_count = max_src;
    }

    ImGui::Text("Configured: %d / %d", wizard.cfg.source_count, max_src);
    ImGui::Separator();

    ImGui::TextUnformatted("Source list");
    ImGui::Indent();
    for (int i = 0; i < max_src; ++i) {
        Source& s = sim->sources[i];
        bool active = (s.active != 0);
        char label[64];
        std::snprintf(label, sizeof(label), "#%d %s", i, source_type_label(s.type));
        bool selected = (app->selected_source == i);
        if (ImGui::Selectable(label, selected)) {
            app->selected_source = i;
        }
        ImGui::SameLine();
        if (ImGui::Checkbox("##on", &active)) {
            s.active = active ? 1 : 0;
            if (i < wizard.cfg.source_count) {
                wizard.cfg.source_configs[i].active = s.active;
            }
        }
        ImGui::SameLine();
        ImGui::Text("%s", source_field_label(s.field));
    }
    ImGui::Unindent();

    int idx = app->selected_source;
    if (idx >= 0 && idx < max_src) {
        ImGui::Separator();
        ImGui::Text("Selected source #%d", idx);
        ImGui::Indent();

        Source& s = sim->sources[idx];

        int type_idx = (int)s.type;
        const char* types[] = { "CW", "Gaussian", "Ricker", "Expr" };
        if (ImGui::Combo("Type", &type_idx, types, IM_ARRAYSIZE(types))) {
            s.type = (SourceType)type_idx;
            source_reparam(&s);
        }

        int field_idx = (int)s.field;
        const char* field_items[] = { "Ez", "Hx", "Hy" };
        if (ImGui::Combo("Field", &field_idx, field_items, IM_ARRAYSIZE(field_items))) {
            s.field = (SourceFieldType)field_idx;
        }

        double amp = s.amp;
        if (ImGui::InputDouble("Amplitude", &amp, 0.1, 1.0, "%.3f")) {
            s.amp = amp;
        }
        double freq = s.freq;
        if (ImGui::InputDouble("Frequency (Hz)", &freq, 1e6, 1e8, "%.3e")) {
            if (freq > 0.0) {
                s.freq = freq;
                source_reparam(&s);
            }
        }

        ImGui::Text("Position: (%d,%d)", s.ix, s.iy);
        ImGui::Unindent();
    }
    ImGui::Unindent();
}

/* Blocks panel ------------------------------------------------------------ */
static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app);

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

static void wizard_grid_tab(WizardState& w) {
    ImGui::TextUnformatted("Domain & Grid");
    ImGui::Separator();

    ImGui::InputInt("nx (cells X)", &w.cfg.nx);
    ImGui::InputInt("ny (cells Y)", &w.cfg.ny);
    ImGui::InputDouble("lx (meters X)", &w.cfg.lx, 0.01, 0.1, "%.3f");
    ImGui::InputDouble("ly (meters Y)", &w.cfg.ly, 0.01, 0.1, "%.3f");

    float cfl = (float)w.cfg.cfl_safety;
    if (ImGui::SliderFloat("CFL safety", &cfl, 0.1f, 0.99f, "%.2f")) {
        w.cfg.cfl_safety = (double)cfl;
    }

    const char* boundary_items[] = { "CPML", "Mur" };
    int bmode = (w.cfg.boundary_mode == SIM_BOUNDARY_MUR) ? 1 : 0;
    if (ImGui::Combo("Boundary", &bmode, boundary_items, IM_ARRAYSIZE(boundary_items))) {
        w.cfg.boundary_mode = (bmode == 1) ? SIM_BOUNDARY_MUR : SIM_BOUNDARY_CPML;
    }

    if (!w.advanced) {
        ImGui::Spacing();
        ImGui::TextWrapped("Tip: Advanced mode exposes sweep and logging options in the Run tab.");
    }
}

static void wizard_materials_tab(WizardState& w) {
    ImGui::TextUnformatted("Materials / Blocks");
    ImGui::Separator();

    int count = w.cfg.material_rect_count;
    ImGui::InputInt("Block count", &count);
    if (count < 0) count = 0;
    if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;
    w.cfg.material_rect_count = count;

    for (int i = 0; i < w.cfg.material_rect_count; ++i) {
        MaterialRectSpec& r = w.cfg.material_rects[i];
        char label[32];
        std::snprintf(label, sizeof(label), "Block %d", i);
        if (ImGui::CollapsingHeader(label, ImGuiTreeNodeFlags_DefaultOpen)) {
            const char* type_items[] = { "Dielectric", "PEC", "PMC" };
            int type_idx = 0;
            if (r.tag == 1) type_idx = 1;
            else if (r.tag == 2) type_idx = 2;
            if (ImGui::Combo("Type", &type_idx, type_items, IM_ARRAYSIZE(type_items))) {
                r.tag = (type_idx == 1) ? 1 : (type_idx == 2) ? 2 : 0;
            }
            ImGui::InputDouble("x0", &r.x0, 0.01, 0.1);
            ImGui::InputDouble("y0", &r.y0, 0.01, 0.1);
            ImGui::InputDouble("x1", &r.x1, 0.01, 0.1);
            ImGui::InputDouble("y1", &r.y1, 0.01, 0.1);
            if (r.tag == 0) {
                ImGui::InputDouble("epsr", &r.epsr, 0.1, 1.0);
                ImGui::InputDouble("sigma", &r.sigma, 0.0, 0.1);
            }
        }
    }
}

static void wizard_sources_tab(WizardState& w) {
    ImGui::TextUnformatted("Sources");
    ImGui::Separator();

    int count = w.cfg.source_count;
    ImGui::InputInt("Source count", &count);
    if (count < 0) count = 0;
    if (count > MAX_SRC) count = MAX_SRC;
    w.cfg.source_count = count;

    static const char* type_items[] = { "CW", "Gaussian", "Ricker", "Expr" };
    static const char* field_items[] = { "Ez", "Hx", "Hy" };

    for (int i = 0; i < w.cfg.source_count; ++i) {
        SourceConfigSpec& s = w.cfg.source_configs[i];
        char label[32];
        std::snprintf(label, sizeof(label), "Source %d", i);
        if (ImGui::CollapsingHeader(label, ImGuiTreeNodeFlags_DefaultOpen)) {
            bool active = (s.active != 0);
            if (ImGui::Checkbox("Active", &active)) {
                s.active = active ? 1 : 0;
            }

            int type_idx = (int)s.type;
            if (type_idx < 0 || type_idx > (int)SRC_EXPR) type_idx = 0;
            if (ImGui::Combo("Type", &type_idx, type_items, IM_ARRAYSIZE(type_items))) {
                s.type = (SourceType)type_idx;
            }

            int field_idx = (int)s.field;
            if (field_idx < 0 || field_idx > (int)SRC_FIELD_HY) field_idx = 0;
            if (ImGui::Combo("Field", &field_idx, field_items, IM_ARRAYSIZE(field_items))) {
                s.field = (SourceFieldType)field_idx;
            }

            ImGui::InputDouble("x (0..1)", &s.x, 0.01, 0.1);
            ImGui::InputDouble("y (0..1)", &s.y, 0.01, 0.1);
            ImGui::InputDouble("amp", &s.amp, 0.1, 1.0);
            ImGui::InputDouble("freq (Hz)", &s.freq, 1e6, 1e9, "%.3e");
            ImGui::InputDouble("sigma2", &s.sigma2, 0.5, 1.0);

            if (s.type == SRC_EXPR) {
                ImGui::InputTextMultiline("expr", s.expr, SOURCE_EXPR_MAX_LEN,
                                          ImVec2(0.0f, ImGui::GetTextLineHeight() * 3));
                ImGui::TextUnformatted("Variables: t (seconds), amp, freq, pi");

                if (s.expr[0] != '\0') {
                    char errbuf[128];
                    ExprProgram* prog = nullptr;
                    if (expr_compile(s.expr, &prog, errbuf, sizeof(errbuf))) {
                        const int N = 128;
                        float values[N];
                        double period = (s.freq > 0.0) ? (1.0 / s.freq) : 1e-9;
                        double t_max = 3.0 * period;
                        if (t_max <= 0.0) t_max = 3e-9;

                        double vmin = 0.0;
                        double vmax = 0.0;
                        bool first = true;
                        for (int k = 0; k < N; ++k) {
                            double t = t_max * (double)k / (double)(N - 1);
                            double v = expr_eval(prog, t, s.amp, s.freq);
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
                        ImGui::PlotLines("expr preview", values, N, 0, nullptr, ymin, ymax,
                                         ImVec2(0.0f, ImGui::GetTextLineHeight() * 4));
                        expr_free(prog);
                    } else {
                        ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.6f, 1.0f),
                                           "Expr error: %s", errbuf);
                    }
                }
            }
        }
    }
}

static void wizard_run_tab(WizardState& w) {
    ImGui::TextUnformatted("Run Settings");
    ImGui::Separator();

    const char* modes[] = { "Fixed steps", "Sweep" };
    int run_mode_idx = (w.cfg.run_mode == SIM_RUN_MODE_SWEEP) ? 1 : 0;
    if (ImGui::Combo("Run mode", &run_mode_idx, modes, IM_ARRAYSIZE(modes))) {
        w.cfg.run_mode = (run_mode_idx == 1) ? SIM_RUN_MODE_SWEEP : SIM_RUN_MODE_FIXED_STEPS;
    }

    if (w.cfg.run_mode == SIM_RUN_MODE_FIXED_STEPS) {
        ImGui::InputInt("Run steps", &w.cfg.run_steps);
    } else {
        ImGui::InputInt("Sweep points", &w.cfg.sweep_points);
        ImGui::InputDouble("Sweep start (Hz)", &w.cfg.sweep_start_hz, 1e6, 1e8, "%.3e");
        ImGui::InputDouble("Sweep stop (Hz)", &w.cfg.sweep_stop_hz, 1e6, 1e8, "%.3e");
        ImGui::InputInt("Steps/point", &w.cfg.sweep_steps_per_point);
    }

    if (w.advanced) {
        ImGui::Separator();
        bool profile = (w.cfg.enable_profile != 0);
        if (ImGui::Checkbox("Enable profiling", &profile)) {
            w.cfg.enable_profile = profile ? 1 : 0;
        }
        bool probe_log = (w.cfg.enable_probe_log != 0);
        if (ImGui::Checkbox("Enable probe log", &probe_log)) {
            w.cfg.enable_probe_log = probe_log ? 1 : 0;
        }
        ImGui::InputText("Probe log path", w.cfg.probe_log_path, SIM_PROBE_LOG_PATH_MAX);
    }
}

static bool wizard_draw(WizardState& w) {
    bool apply = false;

    if (!w.open) return false;
    ImGui::SetNextWindowSize(ImVec2(380.0f, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowPos(ImVec2(12.0f, 50.0f), ImGuiCond_Always);
    if (!ImGui::Begin("Simulation Wizard", &w.open)) {
        ImGui::End();
        return false;
    }

    ImGui::Separator();

    if (ImGui::BeginTabBar("WizardTabs")) {
        if (ImGui::BeginTabItem("Grid")) {
            wizard_grid_tab(w);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Materials")) {
            wizard_materials_tab(w);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Sources")) {
            wizard_sources_tab(w);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Run")) {
            wizard_run_tab(w);
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }

    ImGui::Separator();
    if (ImGui::Button("Apply & Restart")) {
        apply = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Close")) {
        w.open = false;
    }

    ImGui::End();
    return apply;
}

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

static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app) {
    (void)app;
    if (!ImGui::CollapsingHeader("Blocks", ImGuiTreeNodeFlags_DefaultOpen)) return;

    ImGui::Indent();
    int count = wizard.cfg.material_rect_count;
    if (ImGui::InputInt("Count", &count)) {
        if (count < 0) count = 0;
        if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;
        wizard.cfg.material_rect_count = count;
        if (sim && bootstrap) {
            apply_wizard_materials_to_sim(wizard, bootstrap, sim);
        }
    }

    for (int i = 0; i < wizard.cfg.material_rect_count; ++i) {
        MaterialRectSpec& r = wizard.cfg.material_rects[i];
        char label[32];
        std::snprintf(label, sizeof(label), "Block %d", i);
        if (ImGui::TreeNode(label)) {
            const char* type_items[] = { "Dielectric", "PEC", "PMC" };
            int type_idx = 0;
            if (r.tag == 1) type_idx = 1;
            else if (r.tag == 2) type_idx = 2;
            if (ImGui::Combo("Type", &type_idx, type_items, IM_ARRAYSIZE(type_items))) {
                r.tag = (type_idx == 1) ? 1 : (type_idx == 2) ? 2 : 0;
                if (sim && bootstrap) {
                    apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                }
            }
            bool changed = false;
            changed |= ImGui::InputDouble("x0", &r.x0, 0.01, 0.1);
            changed |= ImGui::InputDouble("y0", &r.y0, 0.01, 0.1);
            changed |= ImGui::InputDouble("x1", &r.x1, 0.01, 0.1);
            changed |= ImGui::InputDouble("y1", &r.y1, 0.01, 0.1);
            if (r.tag == 0) {
                changed |= ImGui::InputDouble("epsr", &r.epsr, 0.1, 1.0);
                changed |= ImGui::InputDouble("sigma", &r.sigma, 0.0, 0.1);
            }
            if (changed && sim && bootstrap) {
                apply_wizard_materials_to_sim(wizard, bootstrap, sim);
            }
            ImGui::TreePop();
        }
    }
    ImGui::Unindent();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim) {
    if (!ImGui::CollapsingHeader("Probes", ImGuiTreeNodeFlags_DefaultOpen)) return;

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
    if (!ImGui::CollapsingHeader("Probes", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }
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
    app.show_wizard_panel = true;
    app.show_probes_panel = true;
    app.show_log_panel = true;
    app.show_expression_panel = true;
    app.show_scenes_panel = true;
    app.show_grid_overlay = true;
    app.log_count = 0;
    app.viewport_pos = ImVec2(0.0f, 0.0f);
    app.viewport_size = ImVec2(0.0f, 0.0f);
    app.viewport_valid = false;
    ui_log_add(&app, "Simulation started");

    int scale = 2;
    int width = sim->nx * scale;
    int height = sim->ny * scale;
    if (width < 1920) width = 1920;
    if (height < 1080) height = 1080;

    RenderContext* render = render_init("emwave-c (ImGui)", width, height);
    if (!render) {
        std::fprintf(stderr, "Failed to initialise SDL renderer for ImGui front-end\n");
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }
    render->scale = scale;

    // Set up Dear ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable docking
    io.FontGlobalScale = 1.0f;
    io.IniFilename = "imgui.ini";  // Enable layout persistence

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.ScrollbarSize = 14.0f;
    style.FramePadding = ImVec2(8.0f, 4.0f);
    style.ItemSpacing = ImVec2(8.0f, 6.0f);
    style.ItemInnerSpacing = ImVec2(6.0f, 4.0f);
    style.IndentSpacing = 20.0f;
    style.WindowPadding = ImVec2(8.0f, 8.0f);

    ImVec4* colors = style.Colors;
    colors[ImGuiCol_ChildBg] = ImVec4(0.12f, 0.12f, 0.13f, 1.00f);
    colors[ImGuiCol_Header] = ImVec4(0.26f, 0.59f, 0.98f, 0.31f);
    colors[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
    colors[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);

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
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }

    bool running = true;
    bool paused = false;
    int steps_per_frame = 5;  // Speed control
    bool paint_mode = false;
    int paint_material_type = 0;  // 0=PEC, 1=PMC, 2=dielectric
    double paint_epsilon = 4.0;   // For dielectric
    int paint_brush_size = 3;     // Radius in grid cells

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

                // Handle keyboard shortcuts (only if ImGui doesn't want keyboard input)
                if (!io.WantCaptureKeyboard) {
                    switch (key) {
                        case SDLK_ESCAPE:
                        case SDLK_q:
                            running = false;
                            break;
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
                        case SDLK_i: {
                            paint_material_type = (paint_material_type + 1) % 3;
                            const char* mlabel = (paint_material_type == 0)
                                                     ? "PEC"
                                                     : (paint_material_type == 1) ? "PMC"
                                                                                  : "Dielectric";
                            ui_log_add(&app, "Paint material: %s", mlabel);
                            break;
                        }
                        case SDLK_o:
                            paint_epsilon -= 0.5;
                            if (paint_epsilon < 1.0) paint_epsilon = 1.0;
                            if (paint_epsilon > 20.0) paint_epsilon = 20.0;
                            ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            break;
                        case SDLK_p:
                            paint_epsilon += 0.5;
                            if (paint_epsilon < 1.0) paint_epsilon = 1.0;
                            if (paint_epsilon > 20.0) paint_epsilon = 20.0;
                            ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            break;

                        // Grid overlay
                        case SDLK_g:
                            app.show_grid_overlay = !app.show_grid_overlay;
                            break;

                        // Screenshot
                        case SDLK_F2:
                            if (save_screenshot(render, "frame.bmp")) {
                                ui_log_add(&app, "Screenshot saved: frame.bmp");
                            } else {
                                ui_log_add(&app, "Screenshot failed");
                            }
                            break;

                        // FFT export
                        case SDLK_F3:
                            if (dump_scope_fft_csv(&scope, "scope_fft.csv", sim->dt, 1024)) {
                                ui_log_add(&app, "FFT exported: scope_fft.csv");
                            } else {
                                ui_log_add(&app, "FFT export failed");
                            }
                            break;

                        // Scene presets
                        case SDLK_F5: {
                            // Load waveguide preset
                            SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                            char errbuf[256];
                            if (config_loader_parse_file("configs/waveguide.json", &cfg, errbuf, sizeof(errbuf))) {
                                if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                    wizard_init_from_config(wizard, &cfg);
                                    ui_log_add(&app, "Loaded: waveguide.json");
                                }
                            }
                            break;
                        }
                        case SDLK_F6: {
                            // Load CPW filter preset
                            SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                            char errbuf[256];
                            if (config_loader_parse_file("configs/cpw_filter.json", &cfg, errbuf, sizeof(errbuf))) {
                                if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                    wizard_init_from_config(wizard, &cfg);
                                    ui_log_add(&app, "Loaded: cpw_filter.json");
                                }
                            }
                            break;
                        }

                        default:
                            break;
                    }
                } else {
                    // ImGui wants keyboard, but still handle ESC and Space
                    if (key == SDLK_ESCAPE) {
                        running = false;
                    } else if (key == SDLK_SPACE) {
                        paused = !paused;
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                if (!io.WantCaptureMouse) {
                    if (!app.viewport_valid) {
                        continue;
                    }
                    int mx = e.button.x;
                    int my = e.button.y;
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
                    app.last_click_i = ix;
                    app.last_click_j = iy;
                    if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
                        if (paint_mode) {
                            apply_selected_paint(sim, ix, iy, paint_brush_size, paint_material_type, paint_epsilon);
                        } else if (app.placing_source) {
                            int idx = app.selected_source;
                            if (idx < 0) idx = 0;
                            if (idx >= MAX_SRC) idx = MAX_SRC - 1;
                            int safe_ix = ix;
                            int safe_iy = iy;
                            if (safe_ix < 1) safe_ix = 1;
                            if (safe_ix > sim->nx - 2) safe_ix = sim->nx - 2;
                            if (safe_iy < 1) safe_iy = 1;
                            if (safe_iy > sim->ny - 2) safe_iy = sim->ny - 2;
                            sim->sources[idx].ix = safe_ix;
                            sim->sources[idx].iy = safe_iy;
                            if (idx < wizard.cfg.source_count) {
                                double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
                                double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
                                wizard.cfg.source_configs[idx].x = (double)safe_ix / nx1;
                                wizard.cfg.source_configs[idx].y = (double)safe_iy / ny1;
                            }
                        } else if (app.placing_block) {
                            if (!app.block_first_set) {
                                app.block_first_set = true;
                                app.block_first_i = ix;
                                app.block_first_j = iy;
                            } else {
                                int idx = app.selected_block;
                                if (idx < 0) idx = 0;
                                if (idx >= CONFIG_MAX_MATERIAL_RECTS) idx = CONFIG_MAX_MATERIAL_RECTS - 1;
                                if (wizard.cfg.material_rect_count <= idx) {
                                    wizard.cfg.material_rect_count = idx + 1;
                                    if (wizard.cfg.material_rect_count > CONFIG_MAX_MATERIAL_RECTS) {
                                        wizard.cfg.material_rect_count = CONFIG_MAX_MATERIAL_RECTS;
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
                }
            } else if (e.type == SDL_MOUSEMOTION) {
                if (!io.WantCaptureMouse && paint_mode && (e.motion.state & SDL_BUTTON_LMASK)) {
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
                    apply_selected_paint(sim, ix, iy, paint_brush_size, paint_material_type, paint_epsilon);
                }
            }
        }

        // Tie wizard advanced flag to app-level basic/advanced mode
        wizard.advanced = !app.basic_mode;

        ImGui_ImplSDL2_NewFrame();
        ImGui_ImplSDLRenderer2_NewFrame();
        ImGui::NewFrame();

        // Global menu bar
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("View")) {
                ImGui::MenuItem("Scene", nullptr, &app.show_scene_panel);
                ImGui::MenuItem("Sources", nullptr, &app.show_sources_panel);
                ImGui::MenuItem("Blocks", nullptr, &app.show_blocks_panel);
                ImGui::MenuItem("Probes", nullptr, &app.show_probes_panel);
                ImGui::Separator();
                ImGui::MenuItem("Run Controls", nullptr, &app.show_run_panel);
                ImGui::MenuItem("Wizard", nullptr, &app.show_wizard_panel);
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
        if (app.show_run_panel && ImGui::CollapsingHeader("Simulation Controls", ImGuiTreeNodeFlags_DefaultOpen)) {
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
                const char* types[] = {"PEC", "PMC", "Dielectric"};
                const int type_idx = (paint_material_type >= 0 && paint_material_type < 3)
                                         ? paint_material_type
                                         : 0;
                ImGui::Text("Type: %s (I)", types[type_idx]);
                if (type_idx == 2) {
                    ImGui::Text("Epsilon: %.1f (O/P)", paint_epsilon);
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
            ImGui::Checkbox("Basic mode", &app.basic_mode);

            ImGui::Separator();
            if (ImGui::Button(app.show_grid_overlay ? "Hide grid" : "Show grid")) {
                app.show_grid_overlay = !app.show_grid_overlay;
            }
        }

        // Interactive Tools section
        if (ImGui::CollapsingHeader("Interactive Tools", ImGuiTreeNodeFlags_DefaultOpen)) {
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

        // Embedded Wizard UI
        if (app.show_wizard_panel && ImGui::CollapsingHeader("Wizard", ImGuiTreeNodeFlags_None)) {
            bool apply_wizard = false;
            ImGui::PushID("WizardContent");

            if (ImGui::BeginTabBar("WizardTabs", ImGuiTabBarFlags_None)) {
                if (ImGui::BeginTabItem("Grid")) {
                    wizard_grid_tab(wizard);
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("Materials")) {
                    wizard_materials_tab(wizard);
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("Sources")) {
                    wizard_sources_tab(wizard);
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("Run")) {
                    wizard_run_tab(wizard);
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }

            ImGui::Separator();
            if (ImGui::Button("Apply & Restart")) {
                apply_wizard = true;
            }

            ImGui::PopID();

            // Handle wizard apply (rebootstrap simulation)
            if (apply_wizard) {
                if (rebootstrap_simulation(&wizard.cfg, &bootstrap, &sim, &scope, scale)) {
                    width = sim->nx * scale;
                    height = sim->ny * scale;
                    if (width < 1920) width = 1920;
                    if (height < 1080) height = 1080;
                    SDL_SetWindowSize(render->window, width, height);
                    ui_log_add(&app, "Rebooted simulation from wizard config");
                }
            }
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
                    float ymin = (float)vmin;
                    float ymax = (float)vmax_scope;
                    if (ymin == ymax) {
                        ymin -= 1.0f;
                        ymax += 1.0f;
                    }
                    ImVec2 plot_size = ImVec2(ImGui::GetContentRegionAvail().x, bottom_h - 60.0f);
                    ImGui::PlotLines("Probe Ez(center)", values, max_samples, 0, nullptr, ymin, ymax, plot_size);
                } else {
                    ImGui::TextUnformatted("Scope: no data available");
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
                    float sx = app.viewport_pos.x + (float)(s.ix * scale);
                    float sy = app.viewport_pos.y + (float)(s.iy * scale);
                    char label[16];
                    std::snprintf(label, sizeof(label), "%d", k);
                    dl->AddText(ImVec2(sx + 6.0f, sy - 6.0f), src_col, label);
                }

                // Block outlines from wizard config
                int rect_count = wizard.cfg.material_rect_count;
                if (rect_count < 0) rect_count = 0;
                if (rect_count > CONFIG_MAX_MATERIAL_RECTS) rect_count = CONFIG_MAX_MATERIAL_RECTS;
                for (int i = 0; i < rect_count; ++i) {
                    const MaterialRectSpec& r = wizard.cfg.material_rects[i];
                    float x0 = app.viewport_pos.x + (float)(r.x0 * (double)sim->nx * (double)scale);
                    float x1 = app.viewport_pos.x + (float)(r.x1 * (double)sim->nx * (double)scale);
                    float y0 = app.viewport_pos.y + (float)(r.y0 * (double)sim->ny * (double)scale);
                    float y1 = app.viewport_pos.y + (float)(r.y1 * (double)sim->ny * (double)scale);
                    ImVec2 p0(x0, y0);
                    ImVec2 p1(x1, y1);
                    ImU32 col = (i == app.selected_block) ? block_sel_col : block_col;
                    dl->AddRect(p0, p1, col, 0.0f, 0, 1.0f);

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
                            ImU32 col = IM_COL32(255, 255, 255, 200);
                            if (paint_material_type == 0) {
                                col = IM_COL32(255, 80, 80, 200);
                            } else if (paint_material_type == 1) {
                                col = IM_COL32(80, 160, 255, 200);
                            } else {
                                col = IM_COL32(120, 255, 120, 200);
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

        SDL_SetRenderDrawColor(render->renderer, 0, 0, 0, 255);
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

        double vmax = sim->step_Ez_absmax;
        if (vmax <= 0.0) vmax = 1.0;
        render_field_heatmap(render, sim, vmax, 1.0);
        if (app.show_grid_overlay) {
            SDL_Color grid_col = {40, 40, 50, 255};
            render_grid_overlay(render, sim, grid_col);
        }
        render_sources(render, sim->sources);
        if (viewport_set) {
            SDL_RenderSetViewport(render->renderer, NULL);
        }

        ImGui::Render();
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), render->renderer);
        SDL_RenderPresent(render->renderer);

        SDL_Delay(UI_DELAY_MS);
    }

    scope_free(&scope);
    ImGui_ImplSDLRenderer2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();
    render_free(render);
    simulation_bootstrap_shutdown(&bootstrap);

    return 0;
}
