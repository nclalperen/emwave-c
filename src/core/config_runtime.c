// =============================================================================
// emwave-c: Simulation configuration runtime helpers
// Centralised defaults, clamping, validation, and summary printing.
// =============================================================================

#include "config_loader.h"

#include <math.h>
#include <stdio.h>

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
    .run_mode = SIM_RUN_MODE_FIXED_STEPS,
    .run_steps = 5 * 200,
    .enable_probe_log = 0,
    .probe_log_path = "probe.txt",
    .material_rect_count = 1,
    .material_rects = {
        { .x0 = 0.4, .y0 = 0.45, .x1 = 0.6, .y1 = 0.55,
          .epsr = EPSR_MAX_SCENE, .sigma = SIGMA_BLOCK, .tag = 0 }
    },
    .source_count = 1,
    .source_configs = {
        { .active = 1, .x = 0.25, .y = 0.5, .type = SRC_CW,
          .amp = 1.0, .freq = 1e9, .sigma2 = 4.0 }
    },
};

static double clamp01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
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
    if (cfg->run_steps < 0) cfg->run_steps = 0;
    if (cfg->material_rect_count < 0) cfg->material_rect_count = 0;
    if (cfg->material_rect_count > CONFIG_MAX_MATERIAL_RECTS) {
        cfg->material_rect_count = CONFIG_MAX_MATERIAL_RECTS;
    }
    for (int i = 0; i < cfg->material_rect_count; i++) {
        cfg->material_rects[i].x0 = clamp01(cfg->material_rects[i].x0);
        cfg->material_rects[i].y0 = clamp01(cfg->material_rects[i].y0);
        cfg->material_rects[i].x1 = clamp01(cfg->material_rects[i].x1);
        cfg->material_rects[i].y1 = clamp01(cfg->material_rects[i].y1);
        if (cfg->material_rects[i].epsr <= 0.0) cfg->material_rects[i].epsr = 1.0;
        if (cfg->material_rects[i].sigma < 0.0) cfg->material_rects[i].sigma = 0.0;
    }

    if (cfg->source_count < 0) cfg->source_count = 0;
    if (cfg->source_count > MAX_SRC) cfg->source_count = MAX_SRC;
    for (int i = 0; i < cfg->source_count; i++) {
        cfg->source_configs[i].x = clamp01(cfg->source_configs[i].x);
        cfg->source_configs[i].y = clamp01(cfg->source_configs[i].y);
        if (cfg->source_configs[i].sigma2 <= 0.0) cfg->source_configs[i].sigma2 = 4.0;
        if (cfg->source_configs[i].freq <= 0.0) cfg->source_configs[i].freq = 1e9;
    }
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
    for (int i = 0; i < cfg->material_rect_count; i++) {
        const MaterialRectSpec* rect = &cfg->material_rects[i];
        if (rect->x1 <= rect->x0 || rect->y1 <= rect->y0) {
            if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Material rect %d has invalid bounds", i);
            return 0;
        }
        if (rect->epsr <= 0.0) {
            if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Material rect %d has invalid epsr", i);
            return 0;
        }
    }
    for (int i = 0; i < cfg->source_count; i++) {
        const SourceConfigSpec* src = &cfg->source_configs[i];
        if (src->freq <= 0.0) {
            if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Source %d has invalid freq", i);
            return 0;
        }
        if (src->sigma2 <= 0.0) {
            if (errbuf && errbuf_len) snprintf(errbuf, errbuf_len, "Source %d has invalid sigma2", i);
            return 0;
        }
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
    const char* mode = (cfg->run_mode == SIM_RUN_MODE_SWEEP) ? "sweep" : "fixed";
    if (cfg->run_mode == SIM_RUN_MODE_FIXED_STEPS) {
        printf("Run mode: %s (%d steps)\n", mode, cfg->run_steps);
    } else {
        printf("Run mode: %s (%d steps/pt)\n", mode, cfg->sweep_steps_per_point);
    }
    if (cfg->enable_probe_log && cfg->probe_log_path[0]) {
        printf("Probe logging: enabled -> %s\n", cfg->probe_log_path);
    } else {
        printf("Probe logging: disabled\n");
    }
    printf("Materials: %d rectangles, Sources: %d\n",
           cfg->material_rect_count, cfg->source_count);
}

