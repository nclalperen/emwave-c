// =============================================================================
// emwave-c: Shared simulation runner implementation
// =============================================================================

#include "cli_runner.h"

#include "analysis.h"
#include "fdtd_core.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static int clamp_progress_interval(int interval) {
    if (interval <= 0) return 50;
    return interval;
}

void simulation_runner_options_from_config(SimulationRunnerOptions* opts, const SimulationConfig* cfg) {
    if (!opts) {
        return;
    }
    memset(opts, 0, sizeof(*opts));
    if (cfg) {
        opts->progress_interval = clamp_progress_interval(cfg->steps_per_frame);
        const char* path = cfg->probe_log_path[0] ? cfg->probe_log_path : "probe.txt";
        strncpy(opts->probe_log_path, path, SIM_PROBE_LOG_PATH_MAX - 1);
        opts->probe_log_path[SIM_PROBE_LOG_PATH_MAX - 1] = '\0';
    } else {
        opts->progress_interval = clamp_progress_interval(STEPS_PER_FRAME);
        strncpy(opts->probe_log_path, "probe.txt", SIM_PROBE_LOG_PATH_MAX - 1);
        opts->probe_log_path[SIM_PROBE_LOG_PATH_MAX - 1] = '\0';
    }
}

void simulation_runner_init(SimulationRunner* runner, const SimulationRunnerOptions* opts) {
    if (!runner) return;
    memset(runner, 0, sizeof(*runner));
    if (opts) {
        runner->options = *opts;
    }
}

void simulation_runner_reset_progress(SimulationRunner* runner, int total_steps) {
    if (!runner) return;
    runner->total_steps = total_steps;
    runner->steps_completed = 0;
}

static void simulation_runner_close_log(SimulationRunner* runner) {
    if (!runner) return;
    if (probe_writer_is_open(&runner->probe_writer)) {
        probe_writer_close(&runner->probe_writer);
    }
    runner->log_open_failed = 0;
}

static void simulation_runner_ensure_log_open(SimulationRunner* runner) {
    if (!runner) return;
    if (probe_writer_is_open(&runner->probe_writer)) {
        return;
    }
    if (!runner->options.probe_log_path[0] || runner->log_open_failed) {
        return;
    }
    if (!probe_writer_open(&runner->probe_writer, runner->options.probe_log_path)) {
        fprintf(stderr, "Failed to open probe log '%s'\n", runner->options.probe_log_path);
        runner->log_open_failed = 1;
    }
}

void simulation_runner_on_step(SimulationRunner* runner, SimulationState* sim, double probe_value, int log_probe) {
    if (!runner || !sim) return;

    if (log_probe) {
        simulation_runner_ensure_log_open(runner);
    } else if (probe_writer_is_open(&runner->probe_writer)) {
        simulation_runner_close_log(runner);
    }

    if (log_probe && probe_writer_is_open(&runner->probe_writer)) {
        probe_writer_append(&runner->probe_writer, sim->timestep, probe_value);
    }

    if (sim->ports_on) {
        ports_sample(sim, sim->dx, sim->dy);
    }

    runner->steps_completed++;
    if (runner->total_steps > 0 && runner->options.progress_interval > 0) {
        if ((runner->steps_completed % runner->options.progress_interval) == 0 ||
            runner->steps_completed == runner->total_steps) {
            double percent = ((double)runner->steps_completed / (double)runner->total_steps) * 100.0;
            printf("Progress: %d/%d (%.1f%%)\n", runner->steps_completed, runner->total_steps, percent);
        }
    }
}

void simulation_runner_flush(SimulationRunner* runner) {
    if (!runner) return;
    probe_writer_flush(&runner->probe_writer);
}

void simulation_runner_shutdown(SimulationRunner* runner) {
    if (!runner) return;
    simulation_runner_close_log(runner);
}

const char* simulation_runner_probe_path(const SimulationRunner* runner) {
    if (!runner) return NULL;
    return runner->options.probe_log_path;
}

static int determine_fixed_steps(const SimulationConfig* cfg) {
    if (!cfg) return 0;
    if (cfg->run_steps > 0) {
        return cfg->run_steps;
    }
    long fallback = (long)cfg->sweep_points * (long)cfg->sweep_steps_per_point;
    if (fallback <= 0) fallback = 1000;
    if (fallback > 100000000) fallback = 100000000;
    return (int)fallback;
}

static void activate_ports(SimulationState* sim) {
    if (!sim) return;
    sim->ports_on = 1;
    for (int p = 0; p < MAX_PORTS; ++p) {
        sim->ports[p].active = 1;
        sim->ports[p].head = 0;
    }
}

static int run_fixed_mode(const SimulationConfig* cfg, SimulationState* sim,
                          SimulationRunner* runner, long long* out_steps) {
    int total_steps = determine_fixed_steps(cfg);
    simulation_runner_reset_progress(runner, total_steps);

    int probe_x = sim->nx / 2;
    int probe_y = sim->ny / 2;

    for (int step = 0; step < total_steps; ++step) {
        fdtd_step(sim);
        double probe_val = fdtd_get_Ez(sim, probe_x, probe_y);
        simulation_runner_on_step(runner, sim, probe_val, cfg->enable_probe_log);
    }

    if (out_steps) {
        *out_steps = (long long)total_steps;
    }

    simulation_runner_flush(runner);
    if (cfg->enable_probe_log && runner->options.probe_log_path[0]) {
        printf("Simulation complete. Probe samples saved to %s\n", runner->options.probe_log_path);
    } else {
        printf("Simulation complete.\n");
    }
    return 1;
}

static double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

static int run_sweep_mode(const SimulationConfig* cfg, SimulationState* sim,
                          SimulationRunner* runner, long long* out_steps) {
    int points = cfg->sweep_points > 0 ? cfg->sweep_points : 1;
    int steps_per_point = cfg->sweep_steps_per_point > 0 ? cfg->sweep_steps_per_point : 200;
    if (points <= 0) return 0;

    activate_ports(sim);
    printf("# freq_Hz,s21_mag\n");

    FILE* sweep_csv = fopen("sweep_s21.csv", "w");
    if (sweep_csv) {
        fprintf(sweep_csv, "# freq_Hz,s21_mag\n");
        fprintf(sweep_csv, "freq_Hz,s21_mag\n");
    }

    for (int idx = 0; idx < points; ++idx) {
        double t = (points > 1) ? ((double)idx / (double)(points - 1)) : 0.0;
        double freq = lerp(cfg->sweep_start_hz, cfg->sweep_stop_hz, t);

        fdtd_reset(sim);
        fdtd_update_grid_for_freq(sim, freq);
        simulation_runner_reset_progress(runner, steps_per_point);

        int probe_x = sim->nx / 2;
        int probe_y = sim->ny / 2;
        for (int step = 0; step < steps_per_point; ++step) {
            fdtd_step(sim);
            double probe_val = fdtd_get_Ez(sim, probe_x, probe_y);
            simulation_runner_on_step(runner, sim, probe_val, cfg->enable_probe_log);
        }

        double s21 = compute_s21(sim->ports, freq, sim->dt);
        printf("%.9e,%.6e\n", freq, s21);
        if (sweep_csv) {
            fprintf(sweep_csv, "%.12e,%.12e\n", freq, s21);
        }
    }

    simulation_runner_flush(runner);
    if (out_steps) {
        *out_steps = (long long)points * (long long)steps_per_point;
    }
    if (sweep_csv) {
        fclose(sweep_csv);
        fprintf(stderr, "Wrote sweep_s21.csv (%d points)\n", points);
    }
    return 1;
}

int cli_runner_execute(const SimulationConfig* config, SimulationState* sim) {
    if (!config || !sim) {
        return 0;
    }

    SimulationRunnerOptions opts;
    simulation_runner_options_from_config(&opts, config);
    SimulationRunner runner;
    simulation_runner_init(&runner, &opts);

    long long steps = 0;
    clock_t t0 = 0;
    if (config->enable_profile) {
        t0 = clock();
    }

    int ok = 0;
    if (config->run_mode == SIM_RUN_MODE_SWEEP) {
        ok = run_sweep_mode(config, sim, &runner, &steps);
    } else {
        ok = run_fixed_mode(config, sim, &runner, &steps);
    }

    clock_t t1 = 0;
    if (config->enable_profile) {
        t1 = clock();
    }

    simulation_runner_shutdown(&runner);

    if (config->enable_profile && ok && steps > 0 && t1 > t0) {
        double elapsed = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
        if (elapsed > 0.0) {
            double rate = (double)steps / elapsed;
            printf("Profile: %lld steps in %.3f s (%.1f steps/s)\n",
                   steps, elapsed, rate);
        }
    }
    return ok;
}
