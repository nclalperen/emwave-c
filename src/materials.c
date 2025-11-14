// ============================================================================
// emwave-c: Material Property Management
// ============================================================================

#include "materials.h"
#include "config.h"
#include <math.h>

static int clamp_index(int v, int min, int max) {
    if (v < min) return min;
    if (v > max) return max;
    return v;
}

static void apply_rect_to_grid(SimulationState* state, const MaterialRectSpec* rect) {
    if (!state || !rect) return;
    int nx = state->nx;
    int ny = state->ny;
    int ix0 = clamp_index((int)floor(fmin(rect->x0, rect->x1) * nx), 0, nx - 1);
    int ix1 = clamp_index((int)ceil(fmax(rect->x0, rect->x1) * nx) - 1, 0, nx - 1);
    int iy0 = clamp_index((int)floor(fmin(rect->y0, rect->y1) * ny), 0, ny - 1);
    int iy1 = clamp_index((int)ceil(fmax(rect->y0, rect->y1) * ny) - 1, 0, ny - 1);

    for (int i = ix0; i <= ix1; i++) {
        for (int j = iy0; j <= iy1; j++) {
            if (rect->tag == 1) {
                state->tag_grid[i][j] = 1;  /* PEC */
                state->epsr[i][j] = 1.0;
                state->sigma_map[i][j] = SIGMA_BG;
            } else if (rect->tag == 2) {
                state->tag_grid[i][j] = 2;  /* PMC */
            } else {
                state->tag_grid[i][j] = 0;
                state->epsr[i][j] = rect->epsr;
                state->sigma_map[i][j] = rect->sigma;
            }
        }
    }
}

/* Initialize material properties to defaults */
void materials_init(SimulationState* state) {
    if (!state) return;
    materials_reset_to_defaults(state);

    const SimulationConfig* cfg = &state->config;
    for (int i = 0; i < cfg->material_rect_count; i++) {
        apply_rect_to_grid(state, &cfg->material_rects[i]);
    }
}

/* Reset all materials to vacuum/background */
void materials_reset_to_defaults(SimulationState* state) {
    for (int i = 0; i < state->nx; i++) {
        for (int j = 0; j < state->ny; j++) {
            state->epsr[i][j] = 1.0;  /* Vacuum */
            state->sigma_map[i][j] = SIGMA_BG;
            state->tag_grid[i][j] = 0;  /* Dielectric */
        }
    }
}

/* Paint material at a specific location (for interactive editing) */
void paint_material_at(SimulationState* state, int i, int j, int paint_type, double paint_eps) {
    /* Bounds check */
    if (i < 0 || i >= state->nx || j < 0 || j >= state->ny) return;

    switch (paint_type) {
        case 1:  /* PEC (Perfect Electric Conductor) */
            if (state->tag_grid[i][j] == 1) {
                /* Toggle off */
                state->tag_grid[i][j] = 0;
            } else {
                /* Toggle on */
                state->tag_grid[i][j] = 1;
            }
            /* Reset material to vacuum when toggling PEC */
            state->epsr[i][j] = 1.0;
            state->sigma_map[i][j] = SIGMA_BG;
            break;

        case 2:  /* PMC (Perfect Magnetic Conductor) */
            if (state->tag_grid[i][j] == 2) {
                /* Toggle off */
                state->tag_grid[i][j] = 0;
            } else {
                /* Toggle on */
                state->tag_grid[i][j] = 2;
            }
            break;

        case 3:  /* Dielectric */
            /* Toggle between vacuum and specified permittivity */
            if (fabs(state->epsr[i][j] - paint_eps) < 1e-6) {
                state->epsr[i][j] = 1.0;
                state->sigma_map[i][j] = SIGMA_BG;
            } else {
                state->epsr[i][j] = paint_eps;
                state->sigma_map[i][j] = SIGMA_BG;
            }
            state->tag_grid[i][j] = 0;  /* Ensure dielectric tag */
            break;

        default:
            break;
    }
}
