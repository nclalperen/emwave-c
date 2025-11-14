// =============================================================================
// emwave-c: Runtime configuration loader
// =============================================================================

#include "config_loader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const SimulationConfig SIM_CONFIG_DEFAULTS = {
    .nx = NX_DEFAULT,
    .ny = NY_DEFAULT,
    .lx = LX_DEFAULT,
    .ly = LY_DEFAULT,
    .cfl_safety = CFL_SAFETY_FACTOR,
    .steps_per_frame = STEPS_PER_FRAME,
    .sweep_points = 5,
    .sweep_start_hz = 5e8,
    .sweep_stop_hz = 3e9,
    .sweep_steps_per_point = 200,
};

static int parse_int_arg(const char* value, int* out) {
    if (!value || !out) return 0;
    char* end = NULL;
    long v = strtol(value, &end, 10);
    if (value == end) return 0;
    *out = (int)v;
    return 1;
}

static int parse_double_arg(const char* value, double* out) {
    if (!value || !out) return 0;
    char* end = NULL;
    double v = strtod(value, &end);
    if (value == end) return 0;
    *out = v;
    return 1;
}

void config_clamp_to_limits(SimulationConfig* cfg) {
    if (!cfg) return;
    if (cfg->nx < SIM_MIN_DIM) cfg->nx = SIM_MIN_DIM;
    if (cfg->ny < SIM_MIN_DIM) cfg->ny = SIM_MIN_DIM;
    if (cfg->nx > SIM_MAX_DIM) cfg->nx = SIM_MAX_DIM;
    if (cfg->ny > SIM_MAX_DIM) cfg->ny = SIM_MAX_DIM;

    double cells = (double)cfg->nx * (double)cfg->ny;
    if (cells > (double)SIM_MAX_CELLS) {
        double scale = sqrt((double)SIM_MAX_CELLS / cells);
        cfg->nx = (int)fmax((double)SIM_MIN_DIM, floor(cfg->nx * scale));
        cfg->ny = (int)fmax((double)SIM_MIN_DIM, floor(cfg->ny * scale));
    }

    if (cfg->cfl_safety <= 0.0) cfg->cfl_safety = 0.1;
    if (cfg->cfl_safety >= 1.0) cfg->cfl_safety = 0.99;
    if (cfg->sweep_points < 1) cfg->sweep_points = 1;
    if (cfg->sweep_points > SWEEP_MAX_POINTS) cfg->sweep_points = SWEEP_MAX_POINTS;
    if (cfg->sweep_steps_per_point < 1) cfg->sweep_steps_per_point = 50;
}

int config_validate(const SimulationConfig* cfg, char* errbuf, size_t errbuf_len) {
    if (!cfg) return 0;
    if (cfg->nx < SIM_MIN_DIM || cfg->ny < SIM_MIN_DIM) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Grid too small (min %d)", SIM_MIN_DIM);
        return 0;
    }
    if (cfg->nx > SIM_MAX_DIM || cfg->ny > SIM_MAX_DIM) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Grid exceeds max dim %d", SIM_MAX_DIM);
        return 0;
    }
    if ((double)cfg->nx * (double)cfg->ny > (double)SIM_MAX_CELLS) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Grid cells exceed limit %.0f", (double)SIM_MAX_CELLS);
        return 0;
    }
    if (cfg->cfl_safety <= 0.0 || cfg->cfl_safety >= 1.0) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "CFL safety must be between 0 and 1");
        return 0;
    }
    if (cfg->sweep_points < 1 || cfg->sweep_points > SWEEP_MAX_POINTS) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Sweep points must be 1..%d", SWEEP_MAX_POINTS);
        return 0;
    }
    if (cfg->sweep_start_hz <= 0.0 || cfg->sweep_stop_hz <= cfg->sweep_start_hz) {
        if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Sweep frequency range invalid");
        return 0;
    }
    return 1;
}

static void print_usage(void) {
    printf("Options:\n");
    printf("  --nx=<cells>          Grid cells in X (default %d)\n", NX_DEFAULT);
    printf("  --ny=<cells>          Grid cells in Y (default %d)\n", NY_DEFAULT);
    printf("  --lx=<meters>         Domain length in X (default %.3f m)\n", LX_DEFAULT);
    printf("  --ly=<meters>         Domain length in Y (default %.3f m)\n", LY_DEFAULT);
    printf("  --cfl=<0-1>           CFL safety factor (default %.2f)\n", CFL_SAFETY_FACTOR);
    printf("  --sweep-points=<n>    Number of sweep points (default %d)\n", SIM_CONFIG_DEFAULTS.sweep_points);
    printf("  --sweep-start=<Hz>    Sweep start frequency\n");
    printf("  --sweep-stop=<Hz>     Sweep stop frequency\n");
    printf("  --sweep-steps=<n>     Steps per sweep point\n");
    printf("  --help                Show this help\n");
}

int config_load_from_args(int argc, char** argv, SimulationConfig* out_config) {
    if (!out_config) return 0;
    *out_config = SIM_CONFIG_DEFAULTS;

    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        if (strcmp(arg, "--help") == 0) {
            print_usage();
            return 0;
        } else if (strncmp(arg, "--nx=", 5) == 0) {
            parse_int_arg(arg + 5, &out_config->nx);
        } else if (strncmp(arg, "--ny=", 5) == 0) {
            parse_int_arg(arg + 5, &out_config->ny);
        } else if (strncmp(arg, "--lx=", 5) == 0) {
            parse_double_arg(arg + 5, &out_config->lx);
        } else if (strncmp(arg, "--ly=", 5) == 0) {
            parse_double_arg(arg + 5, &out_config->ly);
        } else if (strncmp(arg, "--cfl=", 6) == 0) {
            parse_double_arg(arg + 6, &out_config->cfl_safety);
        } else if (strncmp(arg, "--sweep-points=", 15) == 0) {
            parse_int_arg(arg + 15, &out_config->sweep_points);
        } else if (strncmp(arg, "--sweep-start=", 14) == 0) {
            parse_double_arg(arg + 14, &out_config->sweep_start_hz);
        } else if (strncmp(arg, "--sweep-stop=", 13) == 0) {
            parse_double_arg(arg + 13, &out_config->sweep_stop_hz);
        } else if (strncmp(arg, "--sweep-steps=", 14) == 0) {
            parse_int_arg(arg + 14, &out_config->sweep_steps_per_point);
        }
    }

    config_clamp_to_limits(out_config);
    char errbuf[128];
    if (!config_validate(out_config, errbuf, sizeof(errbuf))) {
        fprintf(stderr, "Invalid configuration: %s\n", errbuf);
        return 0;
    }

    return 1;
}

void config_print_summary(const SimulationConfig* cfg) {
    if (!cfg) return;
    printf("Grid: %dx%d (%.2fm x %.2fm) CFL=%.2f\n",
           cfg->nx, cfg->ny, cfg->lx, cfg->ly, cfg->cfl_safety);
    printf("Sweep: %d pts, %.3e Hz to %.3e Hz, %d steps/pt\n",
           cfg->sweep_points, cfg->sweep_start_hz, cfg->sweep_stop_hz,
           cfg->sweep_steps_per_point);
}
