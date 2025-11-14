// ============================================================================
// emwave-c: Main Entry Point (Modular Version)
// Clean, minimal main loop that coordinates all modules
// ============================================================================

#include "config.h"
#include "types.h"
#include "fdtd_core.h"
#include "boundary.h"
#include "sources.h"
#include "materials.h"
#include "analysis.h"
#include "ui_render.h"
#include "ui_controls.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(void) {
    /* Initialize all subsystems */
    SimulationState* sim = fdtd_init();
    if (!sim) {
        fprintf(stderr, "Failed to initialize simulation\n");
        return 1;
    }

    UIState* ui = ui_state_init();
    if (!ui) {
        fprintf(stderr, "Failed to initialize UI state\n");
        fdtd_free(sim);
        return 1;
    }

    RenderContext* render = render_init("FDTD Electromagnetic Simulator",
                                        NX * 2 + 260, NY * 2 + 80);
    if (!render) {
        fprintf(stderr, "Failed to initialize rendering\n");
        ui_state_free(ui);
        fdtd_free(sim);
        return 1;
    }

    Scope scope = {0};
    scope_init(&scope, NX * 2);

    FILE* probe_log = probe_open("probe.txt");

    /* Print OpenMP status */
#ifdef _OPENMP
    printf("OpenMP enabled with %d threads\n", omp_get_max_threads());
#else
    printf("OpenMP not available - running sequentially\n");
#endif

    /* FPS tracking */
    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 t_prev = SDL_GetPerformanceCounter();
    double fps_avg = 0.0;

    /* Main loop */
    while (ui->running) {
        /* Handle input */
        if (!ui_handle_events(ui, sim, &scope)) {
            break;  /* User requested exit */
        }

        /* Run simulation steps */
        if (!ui->paused) {
            int steps = 2; /* or from UI state */
            for (int s = 0; s < steps; ++s) {
                fdtd_step(sim);

                /* Sample probe */
                double probe_val = fdtd_get_Ez(sim, ui->probe_x, ui->probe_y);
                scope_push(&scope, probe_val);
                if (probe_log) {
                    probe_log(probe_log, sim->timestep, probe_val);
                }

                /* Sample ports if enabled */
                if (sim->ports_on) {
                    ports_sample(sim, sim->dx, sim->dy);
                }
            }
        }

        /* Calculate FPS */
        Uint64 t_now = SDL_GetPerformanceCounter();
        double dt_s = (double)(t_now - t_prev) / perf_freq;
        t_prev = t_now;
        double fps_inst = (dt_s > 0.0) ? (1.0 / dt_s) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        /* Render frame */
        render_frame(render, sim, &scope, fps_avg, 2, ui->paused, ui->show_legend);

        SDL_Delay(UI_DELAY_MS);
    }

    /* Cleanup */
    if (probe_log) fclose(probe_log);
    scope_free(&scope);
    render_free(render);
    ui_state_free(ui);
    fdtd_free(sim);

    return 0;
}
