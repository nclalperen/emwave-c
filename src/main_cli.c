// =============================================================================
// emwave-c: Headless CLI entry point
// Runs the FDTD solver without SDL dependencies so automated builds/tests can
// exercise the runtime configuration and solver loop on machines lacking
// graphics libraries.
// =============================================================================

#include "config.h"
#include "config_loader.h"
#include "fdtd_core.h"
#include "analysis.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int determine_total_steps(const SimulationConfig* cfg) {
    if (!cfg) return 0;
    long steps = (long)cfg->sweep_points * (long)cfg->sweep_steps_per_point;
    if (steps <= 0) {
        steps = 1000;  /* fallback for short demos */
    }
    if (steps > 100000000) {
        steps = 100000000;  /* guard runaway configs */
    }
    return (int)steps;
}

int main(int argc, char** argv) {
    SimulationConfig config;
    if (!config_load_from_args(argc, argv, &config)) {
        return 0;  /* help was shown or validation failed */
    }
    config_print_summary(&config);

#ifdef _OPENMP
    printf("OpenMP enabled (%d threads)\n", omp_get_max_threads());
#else
    printf("OpenMP not available - running sequentially\n");
#endif

    SimulationState* sim = fdtd_init(&config);
    if (!sim) {
        fprintf(stderr, "Failed to initialise simulation state\n");
        return 1;
    }

    FILE* probe_file = probe_open("probe_cli.txt");
    const int total_steps = determine_total_steps(&config);
    const int progress_interval = config.steps_per_frame > 0 ? config.steps_per_frame : 50;

    for (int step = 0; step < total_steps; ++step) {
        fdtd_step(sim);

        /* Sample a probe at the domain centre */
        int px = sim->nx / 2;
        int py = sim->ny / 2;
        double probe_val = fdtd_get_Ez(sim, px, py);
        if (probe_file) {
            probe_log(probe_file, sim->timestep, probe_val);
        }

        if (sim->ports_on) {
            ports_sample(sim, sim->dx, sim->dy);
        }

        if ((step + 1) % progress_interval == 0 || step + 1 == total_steps) {
            double percent = ((double)(step + 1) / (double)total_steps) * 100.0;
            printf("Progress: %d/%d (%.1f%%)\n", step + 1, total_steps, percent);
        }
    }

    if (probe_file) {
        fclose(probe_file);
    }
    fdtd_free(sim);

    printf("Simulation complete. Probe samples saved to probe_cli.txt\n");
    return 0;
}
