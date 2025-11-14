// =============================================================================
// emwave-c: Shared simulation runner utilities for CLI and SDL front-ends
// =============================================================================

#ifndef EMWAVE_CLI_RUNNER_H
#define EMWAVE_CLI_RUNNER_H

#include "config.h"
#include "types.h"
#include "analysis.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SimulationRunnerOptions {
    int progress_interval;
    char probe_log_path[SIM_PROBE_LOG_PATH_MAX];
} SimulationRunnerOptions;

typedef struct SimulationRunner {
    SimulationRunnerOptions options;
    ProbeLogWriter probe_writer;
    int steps_completed;
    int total_steps;
    int log_open_failed;
} SimulationRunner;

void simulation_runner_options_from_config(SimulationRunnerOptions* opts, const SimulationConfig* cfg);
void simulation_runner_init(SimulationRunner* runner, const SimulationRunnerOptions* opts);
void simulation_runner_reset_progress(SimulationRunner* runner, int total_steps);
void simulation_runner_on_step(SimulationRunner* runner, SimulationState* sim, double probe_value, int log_probe);
void simulation_runner_flush(SimulationRunner* runner);
void simulation_runner_shutdown(SimulationRunner* runner);
const char* simulation_runner_probe_path(const SimulationRunner* runner);

int cli_runner_execute(const SimulationConfig* config, SimulationState* sim);

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_CLI_RUNNER_H */
