// ============================================================================
// emwave-c: Wave Source Management
// ============================================================================

#include "sources.h"
#include "config.h"
#include "expr.h"

#include <math.h>
#include <string.h>
#include <stdio.h>

static int normalized_to_cell(double frac, int n, int pad) {
    double clamped = frac;
    if (clamped < 0.0) clamped = 0.0;
    if (clamped > 1.0) clamped = 1.0;
    int idx = (int)lround(clamped * (double)(n - 1));
    if (idx < pad) idx = pad;
    if (idx >= n - pad) idx = (n > pad) ? n - pad - 1 : idx;
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    return idx;
}

static void source_free_expr(Source* s) {
    if (!s || !s->expr_program) {
        return;
    }
    expr_free((ExprProgram*)s->expr_program);
    s->expr_program = NULL;
}

static void source_compile_expr(Source* s) {
    source_free_expr(s);
    if (!s) return;
    if (s->type != SRC_EXPR) return;
    if (s->expr_text[0] == '\0') return;

    ExprProgram* prog = NULL;
    char errbuf[128];
    if (!expr_compile(s->expr_text, &prog, errbuf, (int)sizeof(errbuf))) {
        fprintf(stderr, "Warning: failed to compile source expression '%s': %s\n",
                s->expr_text, (errbuf[0] != '\0') ? errbuf : "unknown error");
        s->expr_program = NULL;
        return;
    }
    s->expr_program = (void*)prog;
}

/* Initialize all sources with default parameters */
void sources_init(Source* sources, int nx, int ny, const SimulationConfig* cfg) {
    for (int k = 0; k < MAX_SRC; k++) {
        sources[k].active = (k == 0) ? 1 : 0;  /* Only first source active by default */
        sources[k].type = SRC_CW;
        sources[k].field = SRC_FIELD_EZ;
        sources[k].amp = 1.0;
        sources[k].freq = 1e9;  /* 1 GHz */
        sources[k].sigma2 = 4.0;  /* Spatial footprint */
        sources[k].expr_text[0] = '\0';
        sources[k].expr_program = NULL;

        /* Default positions - spread them out */
        sources[k].ix = nx/2 + (k - MAX_SRC/2) * (nx/8);
        sources[k].iy = ny/2;

        /* Clamp to valid range */
        int pad = 2;
        if (sources[k].ix < pad) sources[k].ix = pad;
        if (sources[k].ix >= nx - pad) sources[k].ix = (nx > pad) ? nx - pad - 1 : 0;
        if (sources[k].iy < pad) sources[k].iy = pad;
        if (sources[k].iy >= ny - pad) sources[k].iy = (ny > pad) ? ny - pad - 1 : 0;

        source_reparam(&sources[k]);
    }

    if (cfg && cfg->source_count > 0) {
        int pad = 2;
        for (int k = 0; k < cfg->source_count && k < MAX_SRC; k++) {
            Source* s = &sources[k];
            const SourceConfigSpec* spec = &cfg->source_configs[k];
            s->active = spec->active;
            s->type = spec->type;
            s->amp = spec->amp;
            s->freq = spec->freq;
            s->sigma2 = spec->sigma2;
            s->field = spec->field;
            s->expr_text[0] = '\0';
            s->expr_program = NULL;
            if (spec->expr[0] != '\0') {
                strncpy(s->expr_text, spec->expr, SOURCE_EXPR_MAX_LEN - 1);
                s->expr_text[SOURCE_EXPR_MAX_LEN - 1] = '\0';
            }
            s->ix = normalized_to_cell(spec->x, nx, pad);
            s->iy = normalized_to_cell(spec->y, ny, pad);
            source_reparam(s);
            source_compile_expr(s);
        }
        for (int k = cfg->source_count; k < MAX_SRC; k++) {
            sources[k].active = 0;
            source_free_expr(&sources[k]);
        }
    }
}

/* Recompute pulse parameters based on frequency */
void source_reparam(Source* s) {
    if (s->type == SRC_GAUSS_PULSE || s->type == SRC_RICKER) {
        double cycles_t0 = 6.0;
        double cycles_tau = 2.0;
        s->t0 = (s->freq > 0) ? (cycles_t0 / s->freq) : 0.0;
        s->tau = (s->freq > 0) ? (cycles_tau / s->freq) : 1e-9;
    } else {
        s->t0 = 0.0;
        s->tau = 1e-9;
    }
}

/* Set frequency for all sources */
void sources_set_freq(Source* sources, double f) {
    for (int k = 0; k < MAX_SRC; k++) {
        sources[k].freq = f;
        source_reparam(&sources[k]);
    }
}

/* Cycle through source types */
void sources_cycle_type(Source* sources) {
    for (int k = 0; k < MAX_SRC; k++) {
        if (sources[k].type == SRC_EXPR) {
            continue; /* expression sources are not cycled from the keyboard */
        }
        int next = ((int)sources[k].type + 1) % (int)SRC_EXPR; /* cycle CW/Gauss/Ricker */
        sources[k].type = (SourceType)next;
        source_reparam(&sources[k]);
        source_compile_expr(&sources[k]);
    }
}

/* Compute source time value */
double source_time_value(const Source* s, int t, double dt) {
    double tt = t * dt;

    switch (s->type) {
        case SRC_CW:
            return s->amp * sin(2.0 * M_PI * s->freq * tt);

        case SRC_GAUSS_PULSE: {
            double x = (tt - s->t0) / s->tau;
            double env = exp(-0.5 * x * x);
            return s->amp * env * sin(2.0 * M_PI * s->freq * (tt - s->t0));
        }

        case SRC_RICKER: {
            double a = M_PI * s->freq * (tt - s->t0);
            double e = exp(-a * a);
            return s->amp * (1.0 - 2.0 * a * a) * e;
        }

        case SRC_EXPR: {
            const ExprProgram* prog = (const ExprProgram*)s->expr_program;
            if (!prog) return 0.0;
            return expr_eval(prog, tt, s->amp, s->freq);
        }

        default:
            return 0.0;
    }
}

/* Inject a single source into the selected field component */
void inject_source_into_Ez(Source* s, SimulationState* state, double dt) {
    if (!s->active) return;
    if (!state) return;

    /* Compute source amplitude at this timestep */
    double A = source_time_value(s, state->timestep, dt);

    double** field = NULL;
    switch (s->field) {
        case SRC_FIELD_HX:
            field = state->Hx;
            break;
        case SRC_FIELD_HY:
            field = state->Hy;
            break;
        case SRC_FIELD_EZ:
        default:
            field = state->Ez;
            break;
    }
    if (!field) return;

    /* Soft source injection with spatial Gaussian footprint */
    for (int di = -2; di <= 2; di++) {
        for (int dj = -2; dj <= 2; dj++) {
            int i = s->ix + di;
            int j = s->iy + dj;

            /* Bounds check */
            if (i > 0 && i < state->nx && j > 0 && j < state->ny) {
                double r2 = (double)(di*di + dj*dj);
                double w = exp(-r2 / s->sigma2);
                double val = A * w;

                /* Saturating injection to reduce reflections */
                field[i][j] += val / (1.0 + fabs(val) * 0.1);
            }
        }
    }
}

void sources_shutdown(Source* sources) {
    if (!sources) return;
    for (int k = 0; k < MAX_SRC; ++k) {
        source_free_expr(&sources[k]);
    }
}

/* Inject all active sources */
void inject_all_sources(SimulationState* state) {
    for (int k = 0; k < MAX_SRC; k++) {
        inject_source_into_Ez(&state->sources[k], state, state->dt);
    }
}

/* Find nearest source to a given position (for UI interaction) */
int find_nearest_source(const Source* sources, int mx, int my, int scale, float max_distance) {
    float bestD2 = max_distance * max_distance;
    int best = -1;

    for (int k = 0; k < MAX_SRC; k++) {
        if (!sources[k].active) continue;

        int sx = sources[k].ix * scale;
        int sy = sources[k].iy * scale;
        int dx = mx - sx;
        int dy = my - sy;
        float d2 = (float)(dx*dx + dy*dy);

        if (d2 < bestD2) {
            bestD2 = d2;
            best = k;
        }
    }

    return best;
}
