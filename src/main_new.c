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

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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

    int render_width = sim->nx * 2 + 260;
    int render_height = sim->ny * 2 + 90;
    RenderContext* render = render_init("emwave-c (modular)", render_width, render_height);
    if (!render) {
        fprintf(stderr, "Failed to initialize SDL renderer\n");
        ui_state_free(ui);
        fdtd_free(sim);
        return 1;
    }

    ui_state_set_layout(ui, render->scale, render->ui_height, render->side_panel_width,
                        sim->nx, sim->ny);
    ui_state_sync_with_sim(ui, sim);

    Scope scope = (Scope){0};
    scope_init(&scope, sim->nx * render->scale);

    FILE* probe_file = probe_open("probe.txt");

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

    /* Main loop */
    while (ui->running) {
        ui_update_metrics(ui, sim, &scope);
        if (!ui_handle_events(ui, sim, &scope, render->scale, render->ui_height, render->side_panel_width)) {
            break;
        }

        if (!ui->paused) {
            for (int s = 0; s < ui->steps_per_frame; ++s) {
                fdtd_step(sim);
                double probe_val = fdtd_get_Ez(sim, ui->probe_x, ui->probe_y);
                scope_push(&scope, probe_val);
                if (ui->log_probe && probe_file) {
                    probe_log(probe_file, sim->timestep, probe_val);
                }
                if (sim->ports_on) {
                    ports_sample(sim, sim->dx, sim->dy);
                }
            }
        }

        Uint64 now = SDL_GetPerformanceCounter();
        double dt = (double)(now - prev) / (double)perf_freq;
        prev = now;
        double fps_inst = (dt > 0.0) ? (1.0 / dt) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        render_frame(render, sim, ui, &scope, fps_avg);
        SDL_Delay(UI_DELAY_MS);
    }

    if (probe_file) fclose(probe_file);
    scope_free(&scope);
    render_free(render);
    ui_state_free(ui);
    simulation_bootstrap_shutdown(&bootstrap);

    return 0;
}
