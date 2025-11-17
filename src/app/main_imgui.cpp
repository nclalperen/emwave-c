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
    /* Simple log buffer */
    char log_lines[128][256];
    int log_count;
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
    if (!ImGui::Begin("Scene")) {
        ImGui::End();
        return;
    }

    if (sim) {
        ImGui::Text("Grid: %d x %d", sim->nx, sim->ny);
        ImGui::Text("Domain: %.3f m x %.3f m", sim->lx, sim->ly);
        ImGui::Text("dt: %.3e s", sim->dt);
        ImGui::Separator();
    }

    ImGui::TextUnformatted("Materials");
    ImGui::Text("Rectangles: %d", wizard.cfg.material_rect_count);
    ImGui::Separator();
    ImGui::TextUnformatted("Sources");
    ImGui::Text("Count: %d (max %d)", wizard.cfg.source_count, MAX_SRC);

    ImGui::End();
}

/* Sources panel ----------------------------------------------------------- */
static void draw_sources_panel(SimulationState* sim, WizardState& wizard, AppState* app) {
    if (!ImGui::Begin("Sources")) {
        ImGui::End();
        return;
    }
    if (!sim || !app) {
        ImGui::End();
        return;
    }

    int max_src = MAX_SRC;
    if (wizard.cfg.source_count > max_src) {
        wizard.cfg.source_count = max_src;
    }

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

    int idx = app->selected_source;
    if (idx >= 0 && idx < max_src) {
        ImGui::Separator();
        Source& s = sim->sources[idx];
        ImGui::Text("Selected source #%d", idx);

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
    }

    ImGui::End();
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
    if (!ImGui::Begin("Blocks")) {
        ImGui::End();
        return;
    }

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

    ImGui::End();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim) {
    if (!ImGui::Begin("Probes")) {
        ImGui::End();
        return;
    }

#if EMWAVE_ENABLE_PORTS
    if (!sim) {
        ImGui::End();
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

    ImGui::End();
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
    app.log_count = 0;
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
    io.ConfigFlags &= ~ImGuiConfigFlags_DockingEnable;
    io.FontGlobalScale = 1.0f;
    io.IniFilename = NULL;

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.ScrollbarSize = 16.0f;
    style.FramePadding = ImVec2(6.0f, 4.0f);

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
                if (e.key.keysym.sym == SDLK_ESCAPE) {
                    running = false;
                } else if (e.key.keysym.sym == SDLK_SPACE) {
                    paused = !paused;
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                if (!io.WantCaptureMouse) {
                    int mx = e.button.x;
                    int my = e.button.y;
                    int ix = mx / scale;
                    int iy = my / scale;
                    app.last_click_i = ix;
                    app.last_click_j = iy;
                    if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
                        if (app.placing_source) {
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
                ImGui::MenuItem("Run", nullptr, &app.show_run_panel);
                ImGui::MenuItem("Scope", nullptr, &app.show_scope_window);
                ImGui::MenuItem("Wizard", nullptr, &app.show_wizard_panel);
                ImGui::MenuItem("Probes", nullptr, &app.show_probes_panel);
                ImGui::MenuItem("Log", nullptr, &app.show_log_panel);
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos);
        ImGui::SetNextWindowSize(viewport->Size);
        ImGuiWindowFlags root_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                      ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                                      ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
        ImGui::Begin("RootLayout", nullptr, root_flags);

        ImVec2 full = ImGui::GetContentRegionAvail();
        float bottom_h = full.y * 0.25f;
        float main_h = full.y - bottom_h;
        float left_w = full.x * 0.22f;
        float right_w = full.x * 0.25f;
        float center_w = full.x - left_w - right_w;

        ImVec2 origin = ImGui::GetCursorPos();

        // Left column
        ImGui::SetCursorPos(origin);
        ImGui::BeginChild("LeftColumn", ImVec2(left_w, main_h), true);
        if (app.show_scene_panel) {
            draw_scene_panel(sim, wizard);
        }
        if (app.show_sources_panel) {
            draw_sources_panel(sim, wizard, &app);
        }
        if (app.show_blocks_panel) {
            draw_blocks_panel(wizard, &bootstrap, sim, &app);
        }
        if (app.show_probes_panel) {
            draw_probes_panel(sim);
        }
        ImGui::EndChild();

        // Center viewport frame (empty ImGui frame; SDL field behind it)
        ImGui::SetCursorPos(ImVec2(origin.x + left_w, origin.y));
        ImGui::BeginChild("CenterViewport", ImVec2(center_w, main_h), true);
        ImGui::EndChild();

        // Right column (Run + Wizard + Log)
        ImGui::SetCursorPos(ImVec2(origin.x + left_w + center_w, origin.y));
        ImGui::BeginChild("RightColumn", ImVec2(right_w, main_h), true);
        if (app.show_run_panel && ImGui::CollapsingHeader("Run", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Basic mode", &app.basic_mode);
            ImGui::Separator();
            ImGui::Text("Grid: %d x %d", sim->nx, sim->ny);
            if (!app.basic_mode) {
                ImGui::Text("dt: %.3e s", sim->dt);
                ImGui::Text("freq: %.3f GHz", sim->freq * 1e-9);
            }
            ImGui::Text("timestep: %d", sim->timestep);
            ImGui::Text("FPS (avg): %.1f", fps_avg);

            if (ImGui::Button(paused ? "Resume (Space)" : "Pause (Space)")) {
                paused = !paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Wizard")) {
                wizard.open = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Scope")) {
                app.show_scope_window = !app.show_scope_window;
            }

            ImGui::Separator();
            ImGui::TextUnformatted("Placement & Tools");
            int max_src = MAX_SRC;
            if (app.selected_source < 0) app.selected_source = 0;
            if (app.selected_source >= max_src) app.selected_source = max_src - 1;
            char current_src_label[64];
            const Source& ss = sim->sources[app.selected_source];
            std::snprintf(current_src_label, sizeof(current_src_label),
                          "#%d (%s, %s)", app.selected_source,
                          source_type_label(ss.type),
                          ss.active ? "on" : "off");
            if (ImGui::BeginCombo("Source to move", current_src_label)) {
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
            ImGui::Checkbox("Click to move source", &app.placing_source);

            ImGui::Separator();
            ImGui::TextUnformatted("Blocks (wizard config)");
            ImGui::InputInt("Block index", &app.selected_block);
            if (app.selected_block < 0) app.selected_block = 0;
            if (app.selected_block >= CONFIG_MAX_MATERIAL_RECTS) {
                app.selected_block = CONFIG_MAX_MATERIAL_RECTS - 1;
            }
            ImGui::Checkbox("Draw block with 2 clicks", &app.placing_block);
            if (app.placing_block) {
                if (!app.block_first_set) {
                    ImGui::TextUnformatted("First click: first corner, second click: opposite corner.");
                } else {
                    ImGui::Text("First corner: (%d,%d)", app.block_first_i, app.block_first_j);
                }
                ImGui::TextUnformatted("Note: geometry changes apply after 'Apply & Restart' in the wizard.");
            }

            ImGui::Separator();
            ImGui::Text("Last field click: (%d,%d)", app.last_click_i, app.last_click_j);
        }
        ImGui::EndChild();

        ImGui::End();

        // Wizard window (may trigger rebootstrap)
        if (app.show_wizard_panel) {
            if (wizard_draw(wizard)) {
                if (rebootstrap_simulation(&wizard.cfg, &bootstrap, &sim, &scope, scale)) {
                    // Resize window to match new grid (respect minimum size)
                    width = sim->nx * scale;
                    height = sim->ny * scale;
                    if (width < 1920) width = 1920;
                    if (height < 1080) height = 1080;
                    SDL_SetWindowSize(render->window, width, height);
                    ui_log_add(&app, "Rebooted simulation from wizard config");
                }
            }
        }

        // Scope window
        if (app.show_scope_window && scope.y && scope.n > 0) {
            if (ImGui::Begin("Scope", &app.show_scope_window)) {
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
                ImGui::PlotLines("Probe Ez(center)", values, max_samples, 0, nullptr, ymin, ymax,
                                   ImVec2(0.0f, ImGui::GetTextLineHeight() * 6));
            }
            ImGui::End();
        }

        // Overlay: show source IDs and wizard blocks on top of the field
        {
            ImDrawList* dl = ImGui::GetForegroundDrawList();
            const ImU32 src_col = IM_COL32(255, 255, 255, 230);
            const ImU32 block_col = IM_COL32(255, 255, 0, 160);
            const ImU32 block_sel_col = IM_COL32(0, 255, 128, 200);

            // Source labels
            for (int k = 0; k < MAX_SRC; ++k) {
                const Source& s = sim->sources[k];
                if (!s.active) continue;
                float sx = (float)(s.ix * scale);
                float sy = (float)(s.iy * scale);
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
                float x0 = (float)(r.x0 * (double)sim->nx * (double)scale);
                float x1 = (float)(r.x1 * (double)sim->nx * (double)scale);
                float y0 = (float)(r.y0 * (double)sim->ny * (double)scale);
                float y1 = (float)(r.y1 * (double)sim->ny * (double)scale);
                ImVec2 p0(x0, y0);
                ImVec2 p1(x1, y1);
                ImU32 col = (i == app.selected_block) ? block_sel_col : block_col;
                dl->AddRect(p0, p1, col, 0.0f, 0, 1.0f);

                char blabel[16];
                std::snprintf(blabel, sizeof(blabel), "B%d", i);
                dl->AddText(ImVec2(x0 + 4.0f, y0 + 4.0f), col, blabel);
            }
        }

        if (!paused) {
            fdtd_step(sim);
            int px = sim->nx / 2;
            int py = sim->ny / 2;
            double probe_val = fdtd_get_Ez(sim, px, py);
            scope_push(&scope, probe_val);
        }

        Uint64 now = SDL_GetPerformanceCounter();
        double dt = (double)(now - prev) / (double)perf_freq;
        prev = now;
        double fps_inst = (dt > 0.0) ? (1.0 / dt) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        SDL_SetRenderDrawColor(render->renderer, 0, 0, 0, 255);
        SDL_RenderClear(render->renderer);

        double vmax = sim->step_Ez_absmax;
        if (vmax <= 0.0) vmax = 1.0;
        render_field_heatmap(render, sim, vmax, 1.0);
        render_sources(render, sim->sources);

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
