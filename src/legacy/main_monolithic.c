// ============================================================================
// emwave-c: Modular SDL Front-end
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

#include <SDL2/SDL.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void initialize_boundaries(SimulationState* sim) {
    cpml_on = 1;
    boundary_type = BOUNDARY_CPML;
    cpml_apply_preset(cpml_preset_idx, sim);
    cpml_zero_psi(sim);
}

int main(void) {
    SimulationState* sim = fdtd_init();
    if (!sim) {
        fprintf(stderr, "Failed to allocate SimulationState\n");
        return 1;
    }

    materials_init(sim);
    if (!ports_init(sim->ports, sim->nx, sim->ny)) {
        fprintf(stderr, "Failed to initialize ports\n");
    }
    initialize_boundaries(sim);

    UIState* ui = ui_state_init();
    if (!ui) {
        fprintf(stderr, "Failed to allocate UI state\n");
        fdtd_free(sim);
        return 1;
    }

    RenderContext* render = render_init("emwave-c", sim->nx * 2 + 260, sim->ny * 2 + 90);
    if (!render) {
        fprintf(stderr, "Failed to initialize SDL renderer\n");
        ui_state_free(ui);
        fdtd_free(sim);
        return 1;
    }

    ui_state_set_layout(ui, render->scale, render->ui_height, render->side_panel_width,
                        sim->nx, sim->ny);
    ui_state_sync_with_sim(ui, sim);

    Scope scope = {0};
    scope_init(&scope, sim->nx * render->scale);

    FILE* probe_file = probe_open("probe.txt");

#ifdef _OPENMP
    printf("OpenMP enabled with %d threads\n", omp_get_max_threads());
#else
    printf("OpenMP not available - running sequentially\n");
#endif

    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 prev = SDL_GetPerformanceCounter();
    double fps_avg = 0.0;
    const int render_metrics_log_interval = 120;
    double render_metrics_accum = 0.0;
    int render_metrics_samples = 0;

    while (ui->running) {
        Uint64 metrics_start = SDL_GetPerformanceCounter();
        ui_update_metrics(ui, sim, &scope);
        Uint64 metrics_end = SDL_GetPerformanceCounter();
        double metrics_ms = (double)(metrics_end - metrics_start) * 1000.0 / (double)perf_freq;
        if (!ui_handle_events(ui, sim, &scope, render->scale, render->ui_height, render->side_panel_width)) {
            break;
        }

        if (!ui->paused) {
            for (int s = 0; s < ui->steps_per_frame; ++s) {
                fdtd_step(sim);
                double probe_val = sim->Ez[ui->probe_x][ui->probe_y];
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
        SDL_Delay(UI_DELAY_MS);
    }

    if (probe_file) fclose(probe_file);
    scope_free(&scope);
    render_free(render);
    ui_state_free(ui);
    ports_free(sim->ports);
    fdtd_free(sim);

    return 0;
}
