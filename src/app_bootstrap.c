// =============================================================================
// emwave-c: Bootstrap helpers shared by SDL and CLI entry points
// =============================================================================

#include "app_bootstrap.h"

#include "analysis.h"
#include "boundary.h"
#include "materials.h"
#include "sources.h"

#include <stdio.h>
#include <string.h>

static void configure_boundaries(SimulationState* sim) {
    if (!sim) {
        return;
    }

    cpml_on = 1;
    boundary_type = BOUNDARY_CPML;
    cpml_apply_preset(cpml_preset_idx, sim);
    cpml_zero_psi(sim);
}

static int initialise_state(const SimulationConfig* cfg, SimulationBootstrap* bootstrap) {
    if (!bootstrap) {
        return 0;
    }

    memset(bootstrap, 0, sizeof(*bootstrap));
    bootstrap->config = cfg ? *cfg : SIM_CONFIG_DEFAULTS;

    config_clamp_to_limits(&bootstrap->config);
    config_print_summary(&bootstrap->config);

    bootstrap->sim = fdtd_init(&bootstrap->config);
    if (!bootstrap->sim) {
        fprintf(stderr, "Failed to initialise simulation state\n");
        return 0;
    }

    configure_boundaries(bootstrap->sim);
    return 1;
}

int simulation_bootstrap_with_config(const SimulationConfig* config, SimulationBootstrap* bootstrap) {
    return initialise_state(config, bootstrap);
}

int simulation_bootstrap_from_args(int argc, char** argv, SimulationBootstrap* bootstrap) {
    if (!bootstrap) {
        return -1;
    }

    SimulationConfig parsed_config;
    if (!config_load_from_args(argc, argv, &parsed_config)) {
        return 0; /* help printed or validation failed */
    }

    if (!initialise_state(&parsed_config, bootstrap)) {
        return -1;
    }
    return 1;
}

void simulation_bootstrap_shutdown(SimulationBootstrap* bootstrap) {
    if (!bootstrap) {
        return;
    }
    if (bootstrap->sim) {
        fdtd_free(bootstrap->sim);
        bootstrap->sim = NULL;
    }
}
