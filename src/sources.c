// ============================================================================
// emwave-c: Wave Source Management
// ============================================================================

#include "sources.h"
#include "config.h"
#include <math.h>

/* Initialize all sources with default parameters */
void sources_init(Source* sources, int nx, int ny) {
    for (int k = 0; k < MAX_SRC; k++) {
        sources[k].active = (k == 0) ? 1 : 0;  /* Only first source active by default */
        sources[k].type = SRC_CW;
        sources[k].amp = 1.0;
        sources[k].freq = 1e9;  /* 1 GHz */
        sources[k].sigma2 = 4.0;  /* Spatial footprint */

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
        sources[k].type = (SourceType)((sources[k].type + 1) % 3);
        source_reparam(&sources[k]);
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

        default:
            return 0.0;
    }
}

/* Inject a single source into Ez field */
void inject_source_into_Ez(Source* s, SimulationState* state, double dt) {
    if (!s->active) return;

    /* Compute source amplitude at this timestep */
    double A = source_time_value(s, state->timestep, dt);

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
                state->Ez[i][j] += val / (1.0 + fabs(val) * 0.1);
            }
        }
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
