// =============================================================================
// emwave-c: Shared bootstrap helpers for SDL and CLI front-ends
// =============================================================================

#ifndef EMWAVE_APP_BOOTSTRAP_H
#define EMWAVE_APP_BOOTSTRAP_H

#include "config_loader.h"
#include "fdtd_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SimulationBootstrap {
    SimulationConfig config;
    SimulationState* sim;
} SimulationBootstrap;

/*
 * Parse command line arguments, validate the resulting configuration, and
 * initialise the SimulationState. Returns 1 on success, 0 if help text was
 * displayed / validation failed, and -1 on fatal errors.
 */
int simulation_bootstrap_from_args(int argc, char** argv, SimulationBootstrap* bootstrap);

/* Initialise the SimulationState directly from a provided configuration. */
int simulation_bootstrap_with_config(const SimulationConfig* config, SimulationBootstrap* bootstrap);

/* Free all resources owned by the bootstrap helper. */
void simulation_bootstrap_shutdown(SimulationBootstrap* bootstrap);

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_APP_BOOTSTRAP_H */
