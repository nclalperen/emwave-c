// ============================================================================
// emwave-c: Main Entry Point (Modular Version)
// Clean, minimal main loop that coordinates all modules
// ============================================================================

#include "config.h"
#include "types.h"
#include "analysis.h"
#include "app_bootstrap.h"
#include "ui_render.h"
#include "ui_controls.h"
#include "cli_runner.h"
#include "config_loader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct ScenePresetDesc {
    int id;
    const char* name;
    const char* config_path;
} ScenePresetDesc;

static const ScenePresetDesc UI_SCENE_PRESETS[] = {
    {1, "waveguide", "configs/waveguide.json"},
    {2, "cpw_filter", "configs/cpw_filter.json"}
};

static int path_is_absolute(const char* path) {
    if (!path || !path[0]) {
        return 0;
    }
#if defined(_WIN32)
    if (((path[0] >= 'A' && path[0] <= 'Z') || (path[0] >= 'a' && path[0] <= 'z')) &&
        path[1] == ':' &&
        (path[2] == '\\' || path[2] == '/')) {
        return 1;
    }
    if (path[0] == '/' || path[0] == '\\') {
        return 1;
    }
#else
    if (path[0] == '/') {
        return 1;
    }
#endif
    return 0;
}

static int strip_last_component(char* path) {
    if (!path) {
        return 0;
    }
    size_t len = strlen(path);
    if (len == 0) {
        return 0;
    }
    while (len > 0 && (path[len - 1] == '/' || path[len - 1] == '\\')) {
        path[--len] = '\0';
    }
    while (len > 0 && path[len - 1] != '/' && path[len - 1] != '\\') {
        path[--len] = '\0';
    }
    if (len == 0) {
        return 0;
    }
    /* keep the trailing separator */
    path[len] = '\0';
    return 1;
}

static int try_load_scene_from_root(const char* root, const char* relative_path,
                                    SimulationConfig* cfg, char* resolved_path,
                                    size_t resolved_len, char* errbuf, size_t errbuf_len) {
    if (!root || !relative_path) {
        return 0;
    }
    char candidate[PATH_MAX];
    int written = snprintf(candidate, sizeof(candidate), "%s%s", root, relative_path);
    if (written < 0 || written >= (int)sizeof(candidate)) {
        return 0;
    }
    if (config_loader_parse_file(candidate, cfg, errbuf, errbuf_len)) {
        if (resolved_path && resolved_len > 0) {
            snprintf(resolved_path, resolved_len, "%s", candidate);
        }
        return 1;
    }
    return 0;
}

static int load_scene_preset_config(const ScenePresetDesc* preset, SimulationConfig* cfg,
                                    char* resolved_path, size_t resolved_len,
                                    char* errbuf, size_t errbuf_len) {
    if (!preset || !cfg) {
        return 0;
    }
    if (config_loader_parse_file(preset->config_path, cfg, errbuf, errbuf_len)) {
        if (resolved_path && resolved_len > 0) {
            snprintf(resolved_path, resolved_len, "%s", preset->config_path);
        }
        return 1;
    }
    if (path_is_absolute(preset->config_path)) {
        return 0;
    }
    char* base = SDL_GetBasePath();
    if (!base) {
        return 0;
    }
    char search_root[PATH_MAX];
    snprintf(search_root, sizeof(search_root), "%s", base);
    SDL_free(base);
    for (int depth = 0; depth < 5 && search_root[0] != '\0'; ++depth) {
        if (try_load_scene_from_root(search_root, preset->config_path,
                                     cfg, resolved_path, resolved_len,
                                     errbuf, errbuf_len)) {
            return 1;
        }
        if (!strip_last_component(search_root)) {
            break;
        }
    }
    return 0;
}

static const ScenePresetDesc* find_scene_preset(int id) {
    size_t count = sizeof(UI_SCENE_PRESETS) / sizeof(UI_SCENE_PRESETS[0]);
    for (size_t i = 0; i < count; ++i) {
        if (UI_SCENE_PRESETS[i].id == id) {
            return &UI_SCENE_PRESETS[i];
        }
    }
    return NULL;
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

    /* Initialize UI */
    UIState* ui = ui_state_init();
    if (!ui) {
        fprintf(stderr, "Failed to initialize UI state\n");
        fdtd_free(sim);
        return 1;
    }

    int render_width = sim->nx * RENDER_DEFAULT_SCALE
                     + RENDER_DEFAULT_LEFT_PANEL
                     + RENDER_DEFAULT_RIGHT_PANEL;
    int render_height = sim->ny * RENDER_DEFAULT_SCALE
                      + RENDER_DEFAULT_MENU_BAR
                      + RENDER_DEFAULT_TIMELINE_HEIGHT;
    RenderContext* render = render_init("emwave-c (modular)", render_width, render_height);
    if (!render) {
        fprintf(stderr, "Failed to initialize SDL renderer\n");
        ui_state_free(ui);
        fdtd_free(sim);
        return 1;
    }

    ui_state_set_layout(ui, render->scale, render->menu_bar_height, render->timeline_height,
                        render->left_panel_width, render->right_panel_width,
                        sim->nx, sim->ny);
    ui_state_sync_with_sim(ui, sim);

    Scope scope = (Scope){0};
    if (!scope_init(&scope, sim->nx * render->scale)) {
        fprintf(stderr, "Failed to initialize oscilloscope buffer\n");
        render_free(render);
        ui_state_free(ui);
        fdtd_free(sim);
        return 1;
    }

    SimulationRunnerOptions runner_opts;
    simulation_runner_options_from_config(&runner_opts, &bootstrap.config);
    runner_opts.progress_interval = 0;  /* UI handles its own status */
    SimulationRunner runner;
    simulation_runner_init(&runner, &runner_opts);
    simulation_runner_reset_progress(&runner, 0);

    /* Print OpenMP status */
#ifdef _OPENMP
    printf("OpenMP enabled with %d threads\n", omp_get_max_threads());
#else
    printf("OpenMP not available - running sequentially\n");
#endif

    /* FPS tracking */
    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 prev = SDL_GetPerformanceCounter();
    double fps_avg = 0.0;
    const int render_metrics_log_interval = 120;
    double render_metrics_accum = 0.0;
    int render_metrics_samples = 0;

    /* Main loop */
    while (ui->running) {
        Uint64 metrics_start = SDL_GetPerformanceCounter();
        ui_update_metrics(ui, sim, &scope);
        Uint64 metrics_end = SDL_GetPerformanceCounter();
        double metrics_ms = (double)(metrics_end - metrics_start) * 1000.0 / (double)perf_freq;
        if (!ui_handle_events(ui, sim, &scope)) {
            break;
        }

        if (ui->request_scene_reload) {
            int preset_id = ui->requested_scene_preset;
            ui->request_scene_reload = 0;
            const ScenePresetDesc* preset = find_scene_preset(preset_id);
            if (!preset) {
                fprintf(stderr, "Unknown scene preset %d\n", preset_id);
            } else {
                SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                char errbuf[256];
                char resolved_path[PATH_MAX] = {0};
                if (!load_scene_preset_config(preset, &cfg, resolved_path,
                                              sizeof(resolved_path), errbuf, sizeof(errbuf))) {
                    const char* path_display = preset->config_path;
                    if (resolved_path[0]) {
                        path_display = resolved_path;
                    }
                    fprintf(stderr, "Failed to load scene '%s' from %s: %s\n",
                            preset->name, path_display,
                            (errbuf[0] != '\0') ? errbuf : "unknown error");
                } else {
                    simulation_runner_shutdown(&runner);
                    scope_free(&scope);
                    simulation_bootstrap_shutdown(&bootstrap);

                    if (!simulation_bootstrap_with_config(&cfg, &bootstrap)) {
                        fprintf(stderr, "Failed to initialise scene '%s'\n", preset->name);
                        return 1;
                    }

                    sim = bootstrap.sim;

                    if (!scope_init(&scope, sim->nx * render->scale)) {
                        fprintf(stderr, "Failed to initialize oscilloscope buffer for scene '%s'\n", preset->name);
                        simulation_bootstrap_shutdown(&bootstrap);
                        return 1;
                    }

                    ui_state_set_layout(ui, render->scale, render->menu_bar_height, render->timeline_height,
                                        render->left_panel_width, render->right_panel_width,
                                        sim->nx, sim->ny);
                    ui_state_sync_with_sim(ui, sim);
                    if (preset->name) {
                        snprintf(ui->scene_name, sizeof(ui->scene_name), "%s", preset->name);
                    } else {
                        ui->scene_name[0] = '\0';
                    }

                    simulation_runner_options_from_config(&runner_opts, &bootstrap.config);
                    runner_opts.progress_interval = 0;
                    simulation_runner_init(&runner, &runner_opts);
                    simulation_runner_reset_progress(&runner, 0);

                    const char* path_display = resolved_path[0] ? resolved_path : preset->config_path;
                    printf("Scene switched to '%s' (%s)\n", preset->name, path_display);
                }
            }
            /* Skip stepping/rendering on the old state this frame */
            continue;
        }

        if (!ui->paused) {
            for (int s = 0; s < ui->steps_per_frame; ++s) {
                fdtd_step(sim);
                double probe_val = fdtd_get_Ez(sim, ui->probe_x, ui->probe_y);
                scope_push(&scope, probe_val);
                simulation_runner_on_step(&runner, sim, probe_val, ui->log_probe);
            }
        }

        Uint64 now = SDL_GetPerformanceCounter();
        double dt = (double)(now - prev) / (double)perf_freq;
        prev = now;
        double fps_inst = (dt > 0.0) ? (1.0 / dt) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        Uint64 render_start = SDL_GetPerformanceCounter();
        render_frame(render, sim, ui, &scope, fps_avg);
        Uint64 render_end = SDL_GetPerformanceCounter();
        double render_ms = (double)(render_end - render_start) * 1000.0 / (double)perf_freq;
        double render_metrics_ms = metrics_ms + render_ms;
        render_metrics_accum += render_metrics_ms;
        render_metrics_samples++;
        if (render_metrics_samples >= render_metrics_log_interval) {
            double avg_ms = render_metrics_accum / (double)render_metrics_samples;
            printf("[profile] render+metrics avg %.3f ms (metrics %.3f ms, render %.3f ms)\n",
                   avg_ms, metrics_ms, render_ms);
            render_metrics_accum = 0.0;
            render_metrics_samples = 0;
        }

        if (ui->request_screenshot) {
            ui->request_screenshot = 0;
            if (!save_screenshot(render, "frame.bmp")) {
                fprintf(stderr, "Warning: failed to save screenshot to frame.bmp\n");
            } else {
                printf("Saved screenshot to frame.bmp\n");
            }
        }
        SDL_Delay(UI_DELAY_MS);
    }

    simulation_runner_shutdown(&runner);
    scope_free(&scope);
    render_free(render);
    ui_state_free(ui);
    simulation_bootstrap_shutdown(&bootstrap);

    return 0;
}
