// ============================================================================
// emwave-c: FDTD Core Simulation Engine
// Electromagnetic field updates using Yee grid
// ============================================================================

#include "fdtd_core.h"
#include "config.h"
#include "boundary.h"
#include "sources.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External boundary state */
extern int cpml_on;
extern int cpml_N;
extern double kx[], bx[], cx[];
extern double ky[], by[], cy[];

/* CFL timestep calculation */
double fdtd_compute_dt(double dx, double dy) {
    return CFL_SAFETY_FACTOR / (c0 * sqrt(1.0/(dx*dx) + 1.0/(dy*dy)));
}

/* Initialize simulation state */
SimulationState* fdtd_init(void) {
    SimulationState* state = (SimulationState*)malloc(sizeof(SimulationState));
    if (!state) {
        fprintf(stderr, "Failed to allocate SimulationState\n");
        return NULL;
    }

    /* Allocate field arrays */
    state->Ez = malloc(sizeof(double[NX][NY]));
    if (!state->Ez) { fprintf(stderr, "Failed to allocate Ez\n"); free(state); return NULL; }

    state->Hx = malloc(sizeof(double[NX][NY]));
    if (!state->Hx) { free(state->Ez); free(state); return NULL; }

    state->Hy = malloc(sizeof(double[NX][NY]));
    if (!state->Hy) { free(state->Hx); free(state->Ez); free(state); return NULL; }

    state->Ez_old = malloc(sizeof(double[NX][NY]));
    if (!state->Ez_old) { free(state->Hy); free(state->Hx); free(state->Ez); free(state); return NULL; }

    /* Allocate CPML auxiliary fields */
    state->psi_Ezx = malloc(sizeof(double[NX][NY]));
    if (!state->psi_Ezx) {
        free(state->Ez_old); free(state->Hy); free(state->Hx); free(state->Ez); free(state);
        return NULL;
    }

    state->psi_Ezy = malloc(sizeof(double[NX][NY]));
    if (!state->psi_Ezy) {
        free(state->psi_Ezx); free(state->Ez_old); free(state->Hy); free(state->Hx); free(state->Ez); free(state);
        return NULL;
    }

    state->psi_Hyx = malloc(sizeof(double[NX][NY]));
    if (!state->psi_Hyx) {
        free(state->psi_Ezy); free(state->psi_Ezx); free(state->Ez_old);
        free(state->Hy); free(state->Hx); free(state->Ez); free(state);
        return NULL;
    }

    state->psi_Hxy = malloc(sizeof(double[NX][NY]));
    if (!state->psi_Hxy) {
        free(state->psi_Hyx); free(state->psi_Ezy); free(state->psi_Ezx);
        free(state->Ez_old); free(state->Hy); free(state->Hx); free(state->Ez); free(state);
        return NULL;
    }

    /* Initialize grid parameters */
    state->dx = BASE_DX;
    state->dy = BASE_DY;
    state->dt = fdtd_compute_dt(state->dx, state->dy);

    /* Initialize time */
    state->timestep = 0;
    state->freq = 1e9;  /* Default 1 GHz */

    /* Initialize ports */
    state->ports_on = 0;
    for (int p = 0; p < MAX_PORTS; p++) {
        state->ports[p].active = 0;
        state->ports[p].V = NULL;
        state->ports[p].I = NULL;
    }

    /* Clear all fields */
    fdtd_clear_fields(state);

    /* Initialize materials to vacuum */
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            state->epsr[i][j] = 1.0;
            state->sigma_map[i][j] = SIGMA_BG;
            state->tag_grid[i][j] = 0;
        }
    }

    /* Initialize sources */
    sources_init(state->sources);

    return state;
}

/* Free simulation state */
void fdtd_free(SimulationState* state) {
    if (!state) return;

    free(state->Ez);
    free(state->Hx);
    free(state->Hy);
    free(state->Ez_old);
    free(state->psi_Ezx);
    free(state->psi_Ezy);
    free(state->psi_Hyx);
    free(state->psi_Hxy);

    /* Free port buffers */
    for (int p = 0; p < MAX_PORTS; p++) {
        if (state->ports[p].V) free(state->ports[p].V);
        if (state->ports[p].I) free(state->ports[p].I);
    }

    free(state);
}

/* Clear all fields to zero */
void fdtd_clear_fields(SimulationState* state) {
    memset(state->Ez, 0, sizeof(double[NX][NY]));
    memset(state->Hx, 0, sizeof(double[NX][NY]));
    memset(state->Hy, 0, sizeof(double[NX][NY]));
    memset(state->Ez_old, 0, sizeof(double[NX][NY]));
    memset(state->psi_Ezx, 0, sizeof(double[NX][NY]));
    memset(state->psi_Ezy, 0, sizeof(double[NX][NY]));
    memset(state->psi_Hyx, 0, sizeof(double[NX][NY]));
    memset(state->psi_Hxy, 0, sizeof(double[NX][NY]));
}

/* Reset simulation state */
void fdtd_reset(SimulationState* state) {
    fdtd_clear_fields(state);
    state->timestep = 0;

    /* Reset ports */
    for (int p = 0; p < MAX_PORTS; p++) {
        if (state->ports[p].V) memset(state->ports[p].V, 0, sizeof(double) * state->ports[p].n);
        if (state->ports[p].I) memset(state->ports[p].I, 0, sizeof(double) * state->ports[p].n);
        state->ports[p].head = 0;
    }
}

/* Material property access */
double fdtd_epsilon_at(const SimulationState* state, int i, int j) {
    return EPS0 * state->epsr[i][j];
}

double fdtd_sigma_at(const SimulationState* state, int i, int j) {
    return state->sigma_map[i][j];
}

/* Update grid spacing for frequency */
void fdtd_update_grid_for_freq(SimulationState* state, double freq) {
    double lambda = c0 / freq;
    double delta = lambda / (sqrt(EPSR_MAX_SCENE) * TARGET_CPW);

    if (delta <= 0) delta = BASE_DX;

    double max_delta = BASE_DX * MAX_SCALE_FACTOR;
    if (delta > max_delta) delta = max_delta;

    state->dx = delta;
    state->dy = delta;
    state->dt = fdtd_compute_dt(state->dx, state->dy);
    state->freq = freq;
}

/* Main FDTD step - Update electromagnetic fields by one timestep */
void fdtd_step(SimulationState* state) {
    double dx = state->dx;
    double dy = state->dy;
    double dt = state->dt;

    /* Save Ez for Mur-1 boundaries */
    int i, j;
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            state->Ez_old[i][j] = state->Ez[i][j];
        }
    }

    /* Update H fields with optional CPML */
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for (i = 0; i < NX-1; i++) {
        for (j = 0; j < NY-1; j++) {
            /* dEz/dy -> Hx */
            double dEdy = (state->Ez[i][j+1] - state->Ez[i][j]) / dy;
            if (cpml_on && (j < cpml_N || j >= NY - cpml_N)) {
                state->psi_Ezy[i][j] = by[j] * state->psi_Ezy[i][j] + cy[j] * dEdy;
                dEdy = (dEdy / ky[j]) + state->psi_Ezy[i][j];
            }
            state->Hx[i][j] -= (dt / MU0) * dEdy;

            /* dEz/dx -> Hy */
            double dEdx = (state->Ez[i+1][j] - state->Ez[i][j]) / dx;
            if (cpml_on && (i < cpml_N || i >= NX - cpml_N)) {
                state->psi_Ezx[i][j] = bx[i] * state->psi_Ezx[i][j] + cx[i] * dEdx;
                dEdx = (dEdx / kx[i]) + state->psi_Ezx[i][j];
            }
            state->Hy[i][j] += (dt / MU0) * dEdx;

            /* PMC boundary handling */
            if (state->tag_grid[i][j] == 2) {
                state->Hx[i][j] = 0.0;
                state->Hy[i][j] = 0.0;
            }
        }
    }

    /* Update E fields with optional CPML + conductivity */
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for (i = 1; i < NX; i++) {
        for (j = 1; j < NY; j++) {
            /* dHy/dx */
            double dHdx = (state->Hy[i][j] - state->Hy[i-1][j]) / dx;
            if (cpml_on && (i < cpml_N || i >= NX - cpml_N)) {
                state->psi_Hyx[i][j] = bx[i] * state->psi_Hyx[i][j] + cx[i] * dHdx;
                dHdx = (dHdx / kx[i]) + state->psi_Hyx[i][j];
            }

            /* dHx/dy */
            double dHdy = (state->Hx[i][j] - state->Hx[i][j-1]) / dy;
            if (cpml_on && (j < cpml_N || j >= NY - cpml_N)) {
                state->psi_Hxy[i][j] = by[j] * state->psi_Hxy[i][j] + cy[j] * dHdy;
                dHdy = (dHdy / ky[j]) + state->psi_Hxy[i][j];
            }

            double curlH = dHdx - dHdy;
            double epsij = fdtd_epsilon_at(state, i, j);
            double sigma = fdtd_sigma_at(state, i, j);
            double tmp = 0.5 * sigma * dt / (epsij + DIVISION_SAFETY_EPSILON);
            double ceze = (1.0 - tmp) / (1.0 + tmp);
            double cezh = (dt / epsij) / (1.0 + tmp);

            state->Ez[i][j] = ceze * state->Ez[i][j] + cezh * curlH;

            /* Enforce PEC: tag 1 clamps Ez to zero */
            if (state->tag_grid[i][j] == 1) {
                state->Ez[i][j] = 0.0;
            }
        }
    }

    /* Apply Mur boundaries if CPML is off */
    if (!cpml_on) {
        apply_mur_boundaries(state);
    }

    /* Inject all active sources */
    inject_all_sources(state);

    /* Increment timestep counter */
    state->timestep++;
}
