// ============================================================================
// emwave-c: Boundary Conditions (CPML and Mur)
// ============================================================================

#include "boundary.h"
#include "config.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

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
double kx[NX], bx[NX], cx[NX];
double ky[NY], by[NY], cy[NY];

/* CPML presets */
const CpmlPreset CPML_PRESETS[3] = {
    {"Gentle",     1.0, 3.0, 0.03, 10},
    {"Default",    1.2, 5.0, 0.05, 12},
    {"Aggressive", 1.8, 6.0, 0.08, 16},
};

/* Zero all CPML auxiliary fields */
void cpml_zero_psi(SimulationState* state) {
    if (state->psi_Ezx) memset(state->psi_Ezx, 0, sizeof(double[NX][NY]));
    if (state->psi_Ezy) memset(state->psi_Ezy, 0, sizeof(double[NX][NY]));
    if (state->psi_Hyx) memset(state->psi_Hyx, 0, sizeof(double[NX][NY]));
    if (state->psi_Hxy) memset(state->psi_Hxy, 0, sizeof(double[NX][NY]));
}

/* Build CPML coefficients for given timestep */
void cpml_build_coeffs(double dt) {
    const double m = 3.0;  /* cubic grading */

    /* X-direction coefficients */
    for (int i = 0; i < NX; i++) {
        double rx = 0.0;

        /* Compute fractional distance into PML region */
        if (i < cpml_N) {
            rx = (cpml_N - i) / (double)cpml_N;
        } else if (i >= NX - cpml_N) {
            rx = (i - (NX - cpml_N)) / (double)cpml_N;  /* FIXED: was -1 */
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
    for (int j = 0; j < NY; j++) {
        double ry = 0.0;

        if (j < cpml_N) {
            ry = (cpml_N - j) / (double)cpml_N;
        } else if (j >= NY - cpml_N) {
            ry = (j - (NY - cpml_N)) / (double)cpml_N;  /* FIXED: was -1 */
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
void cpml_apply_preset(int idx, double dt) {
    int n = sizeof(CPML_PRESETS) / sizeof(CPML_PRESETS[0]);
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;

    cpml_preset_idx = idx;
    pml_sigma_max = CPML_PRESETS[idx].smax;
    pml_kappa_max = CPML_PRESETS[idx].kmax;
    pml_alpha_max = CPML_PRESETS[idx].amax;
    cpml_N = CPML_PRESETS[idx].thick;

    /* Validate CPML thickness doesn't exceed grid size */
    int max_thickness = (NX < NY ? NX : NY) / 2 - 1;
    if (cpml_N > max_thickness) {
        fprintf(stderr, "Warning: CPML thickness %d exceeds safe limit %d, clamping\n",
                cpml_N, max_thickness);
        cpml_N = max_thickness;
    }
    if (cpml_N < 1) cpml_N = 1;

    cpml_build_coeffs(dt);
}

/* Apply Mur first-order absorbing boundaries */
void apply_mur_boundaries(SimulationState* state) {
    double dx = state->dx;
    double dy = state->dy;
    double dt = state->dt;

    /* Top and bottom boundaries (y-direction) */
    int i;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (i = 1; i < NX-1; i++) {
        /* Bottom boundary (j=0) - use local wave speed */
        double epsr_bot = state->epsr[i][1];
        double c_bot = c0 / sqrt(epsr_bot);
        state->Ez[i][0] = state->Ez_old[i][1] +
                          ((c_bot*dt - dy)/(c_bot*dt + dy)) * (state->Ez[i][1] - state->Ez_old[i][0]);

        /* Top boundary (j=NY-1) */
        double epsr_top = state->epsr[i][NY-2];
        double c_top = c0 / sqrt(epsr_top);
        state->Ez[i][NY-1] = state->Ez_old[i][NY-2] +
                             ((c_top*dt - dy)/(c_top*dt + dy)) * (state->Ez[i][NY-2] - state->Ez_old[i][NY-1]);
    }

    /* Left and right boundaries (x-direction) */
    int j;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (j = 1; j < NY-1; j++) {
        /* Left boundary (i=0) */
        double epsr_left = state->epsr[1][j];
        double c_left = c0 / sqrt(epsr_left);
        state->Ez[0][j] = state->Ez_old[1][j] +
                          ((c_left*dt - dx)/(c_left*dt + dx)) * (state->Ez[1][j] - state->Ez_old[0][j]);

        /* Right boundary (i=NX-1) */
        double epsr_right = state->epsr[NX-2][j];
        double c_right = c0 / sqrt(epsr_right);
        state->Ez[NX-1][j] = state->Ez_old[NX-2][j] +
                             ((c_right*dt - dx)/(c_right*dt + dx)) * (state->Ez[NX-2][j] - state->Ez_old[NX-1][j]);
    }
}
