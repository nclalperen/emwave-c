// ============================================================================
// emwave-c: Material Property Management
// ============================================================================

#include "materials.h"
#include "config.h"
#include <math.h>

/* Initialize material properties to defaults */
void materials_init(SimulationState* state) {
    materials_reset_to_defaults(state);

    /* Set up central dielectric block (example geometry) */
    int bx0 = state->nx/2 - state->nx/10;
    int bx1 = state->nx/2 + state->nx/10;
    int by0 = state->ny/2 - state->ny/20;
    int by1 = state->ny/2 + state->ny/20;

    for (int i = bx0; i <= bx1; i++) {
        for (int j = by0; j <= by1; j++) {
            if (i >= 0 && i < state->nx && j >= 0 && j < state->ny) {
                state->epsr[i][j] = EPSR_MAX_SCENE;  /* Dielectric block */
                state->sigma_map[i][j] = SIGMA_BLOCK;
            }
        }
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
