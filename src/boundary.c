// ============================================================================
// emwave-c: Boundary Conditions (CPML and Mur)
// ============================================================================

#include "boundary.h"
#include "config.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Global CPML state */
int cpml_on = 0;
int cpml_N = PML_THICK;
double pml_sigma_max = 1.2;
double pml_kappa_max = 5.0;
double pml_alpha_max = 0.05;
int cpml_preset_idx = 1;
BoundaryType boundary_type = BOUNDARY_MUR;

/* CPML coefficient arrays */
double* kx = NULL;
double* bx = NULL;
double* cx = NULL;
double* ky = NULL;
double* by = NULL;
double* cy = NULL;

static int kx_capacity = 0;
static int ky_capacity = 0;

static int ensure_cpml_capacity(int nx, int ny) {
    if (nx > kx_capacity) {
        double* nkx = (double*)realloc(kx, sizeof(double) * nx);
        double* nbx = (double*)realloc(bx, sizeof(double) * nx);
        double* ncx = (double*)realloc(cx, sizeof(double) * nx);
        if (!nkx || !nbx || !ncx) {
            fprintf(stderr, "Failed to allocate CPML X coefficients\n");
            free(nkx); free(nbx); free(ncx);
            return 0;
        }
        kx = nkx; bx = nbx; cx = ncx;
        kx_capacity = nx;
    }
    if (ny > ky_capacity) {
        double* nky = (double*)realloc(ky, sizeof(double) * ny);
        double* nby = (double*)realloc(by, sizeof(double) * ny);
        double* ncy = (double*)realloc(cy, sizeof(double) * ny);
        if (!nky || !nby || !ncy) {
            fprintf(stderr, "Failed to allocate CPML Y coefficients\n");
            free(nky); free(nby); free(ncy);
            return 0;
        }
        ky = nky; by = nby; cy = ncy;
        ky_capacity = ny;
    }
    return 1;
}

/* CPML presets */
const CpmlPreset CPML_PRESETS[3] = {
    {"Gentle",     1.0, 3.0, 0.03, 10},
    {"Default",    1.2, 5.0, 0.05, 12},
    {"Aggressive", 1.8, 6.0, 0.08, 16},
};

/* Zero all CPML auxiliary fields */
void cpml_zero_psi(SimulationState* state) {
    if (!state) return;
    size_t cells = (size_t)state->nx * (size_t)state->ny;
    if (state->psi_Ezx_data) memset(state->psi_Ezx_data, 0, sizeof(double) * cells);
    if (state->psi_Ezy_data) memset(state->psi_Ezy_data, 0, sizeof(double) * cells);
    if (state->psi_Hyx_data) memset(state->psi_Hyx_data, 0, sizeof(double) * cells);
    if (state->psi_Hxy_data) memset(state->psi_Hxy_data, 0, sizeof(double) * cells);
}

/* Build CPML coefficients for given timestep */
void cpml_build_coeffs(const SimulationState* state) {
    if (!state) return;
    double dt = state->dt;
    int nx = state->nx;
    int ny = state->ny;

    if (!ensure_cpml_capacity(nx, ny)) {
        return;
    }

    const double m = 3.0;  /* cubic grading */

    /* X-direction coefficients */
    for (int i = 0; i < nx; i++) {
        double rx = 0.0;

        /* Compute fractional distance into PML region */
        if (i < cpml_N) {
            rx = (cpml_N - i) / (double)cpml_N;
        } else if (i >= nx - cpml_N) {
            rx = (i - (nx - cpml_N)) / (double)cpml_N;  /* FIXED: was -1 */
        }

        double g = rx > 0 ? pow(rx, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - rx);

        kx[i] = (rx > 0) ? kappa : 1.0;

        if (rx > 0) {
            bx[i] = exp(-(sigma/kappa + alpha) * dt);
            cx[i] = (sigma * (bx[i] - 1.0)) / (kappa * (sigma + kappa*alpha) + DIVISION_SAFETY_EPSILON);
        } else {
            bx[i] = 1.0;
            cx[i] = 0.0;
        }
    }

    /* Y-direction coefficients */
    for (int j = 0; j < ny; j++) {
        double ry = 0.0;

        if (j < cpml_N) {
            ry = (cpml_N - j) / (double)cpml_N;
        } else if (j >= ny - cpml_N) {
            ry = (j - (ny - cpml_N)) / (double)cpml_N;  /* FIXED: was -1 */
        }

        double g = ry > 0 ? pow(ry, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - ry);

        ky[j] = (ry > 0) ? kappa : 1.0;

        if (ry > 0) {
            by[j] = exp(-(sigma/kappa + alpha) * dt);
            cy[j] = (sigma * (by[j] - 1.0)) / (kappa * (sigma + kappa*alpha) + DIVISION_SAFETY_EPSILON);
        } else {
            by[j] = 1.0;
            cy[j] = 0.0;
        }
    }
}

/* Apply a CPML preset */
void cpml_apply_preset(int idx, const SimulationState* state) {
    if (!state) return;
    int n = sizeof(CPML_PRESETS) / sizeof(CPML_PRESETS[0]);
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;

    cpml_preset_idx = idx;
    pml_sigma_max = CPML_PRESETS[idx].smax;
    pml_kappa_max = CPML_PRESETS[idx].kmax;
    pml_alpha_max = CPML_PRESETS[idx].amax;
    cpml_N = CPML_PRESETS[idx].thick;

    /* Validate CPML thickness doesn't exceed grid size */
    int min_dim = (state->nx < state->ny ? state->nx : state->ny);
    int max_thickness = min_dim / 2 - 1;
    if (cpml_N > max_thickness) {
        fprintf(stderr, "Warning: CPML thickness %d exceeds safe limit %d, clamping\n",
                cpml_N, max_thickness);
        cpml_N = max_thickness;
    }
    if (cpml_N < 1) cpml_N = 1;

    cpml_build_coeffs(state);
}

/* Apply Mur first-order absorbing boundaries */
void apply_mur_boundaries(SimulationState* state) {
    double dx = state->dx;
    double dy = state->dy;
    double dt = state->dt;

    /* Top and bottom boundaries (y-direction) */
    int nx = state->nx;
    int ny = state->ny;
    if (nx < 2 || ny < 2) return;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 1; i < nx-1; i++) {
        /* Bottom boundary (j=0) - use local wave speed */
        double epsr_bot = state->epsr[i][1];
        double c_bot = c0 / sqrt(epsr_bot);
        state->Ez[i][0] = state->Ez_old[i][1] +
                          ((c_bot*dt - dy)/(c_bot*dt + dy)) * (state->Ez[i][1] - state->Ez_old[i][0]);

        /* Top boundary (j=NY-1) */
        double epsr_top = state->epsr[i][ny-2];
        double c_top = c0 / sqrt(epsr_top);
        state->Ez[i][ny-1] = state->Ez_old[i][ny-2] +
                             ((c_top*dt - dy)/(c_top*dt + dy)) * (state->Ez[i][ny-2] - state->Ez_old[i][ny-1]);
    }

    /* Left and right boundaries (x-direction) */
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int j = 1; j < ny-1; j++) {
        /* Left boundary (i=0) */
        double epsr_left = state->epsr[1][j];
        double c_left = c0 / sqrt(epsr_left);
        state->Ez[0][j] = state->Ez_old[1][j] +
                          ((c_left*dt - dx)/(c_left*dt + dx)) * (state->Ez[1][j] - state->Ez_old[0][j]);

        /* Right boundary (i=NX-1) */
        double epsr_right = state->epsr[nx-2][j];
        double c_right = c0 / sqrt(epsr_right);
        state->Ez[nx-1][j] = state->Ez_old[nx-2][j] +
                             ((c_right*dt - dx)/(c_right*dt + dx)) * (state->Ez[nx-2][j] - state->Ez_old[nx-1][j]);
    }
}
