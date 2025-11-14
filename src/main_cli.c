// =============================================================================
// emwave-c: Headless CLI entry point
// Runs the FDTD solver without SDL dependencies so automated builds/tests can
// exercise the runtime configuration and solver loop on machines lacking
// graphics libraries.
// =============================================================================

#include "config.h"
#include "app_bootstrap.h"
#include "cli_runner.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    SimulationBootstrap bootstrap;
    int bootstrap_status = simulation_bootstrap_from_args(argc, argv, &bootstrap);
    if (bootstrap_status == 0) {
        return 0;  /* help was shown or validation failed */
    }
    if (bootstrap_status < 0) {
        return 1;
    }

    const SimulationConfig* config = &bootstrap.config;
    SimulationState* sim = bootstrap.sim;

#ifdef _OPENMP
    printf("OpenMP enabled (%d threads)\n", omp_get_max_threads());
#else
    printf("OpenMP not available - running sequentially\n");
#endif

    int ok = cli_runner_execute(config, sim);
    simulation_bootstrap_shutdown(&bootstrap);
    return ok ? 0 : 1;
}
