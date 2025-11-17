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

static CpmlState* cpml_state(SimulationState* state) {
    return state ? &state->cpml : NULL;
}

void boundary_init(SimulationState* state) {
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return;

    *cpml = (CpmlState){0};
    cpml->boundary_type = BOUNDARY_MUR;
    cpml->enabled = 0;
    cpml->thickness = PML_THICK;
    cpml->preset_idx = 1;

    const CpmlPreset preset = CPML_PRESETS[cpml->preset_idx];
    cpml->sigma_max = preset.smax;
    cpml->kappa_max = preset.kmax;
    cpml->alpha_max = preset.amax;
}

void boundary_shutdown(SimulationState* state) {
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return;

    free(cpml->kx); cpml->kx = NULL;
    free(cpml->bx); cpml->bx = NULL;
    free(cpml->cx); cpml->cx = NULL;
    free(cpml->ky); cpml->ky = NULL;
    free(cpml->by); cpml->by = NULL;
    free(cpml->cy); cpml->cy = NULL;
    cpml->kx_capacity = 0;
    cpml->ky_capacity = 0;
}

void boundary_set_type(SimulationState* state, BoundaryType type) {
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return;
    cpml->boundary_type = type;
    cpml->enabled = (type == BOUNDARY_CPML);
}

BoundaryType boundary_get_type(const SimulationState* state) {
    const CpmlState* cpml = state ? &state->cpml : NULL;
    return cpml ? cpml->boundary_type : BOUNDARY_MUR;
}

int boundary_is_cpml_enabled(const SimulationState* state) {
    const CpmlState* cpml = state ? &state->cpml : NULL;
    return cpml ? cpml->enabled : 0;
}

int cpml_get_preset_index(const SimulationState* state) {
    const CpmlState* cpml = state ? &state->cpml : NULL;
    return cpml ? cpml->preset_idx : 0;
}

static int ensure_cpml_capacity(SimulationState* state, int nx, int ny) {
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return 0;

    if (nx > cpml->kx_capacity) {
        double* new_ptr = (double*)realloc(cpml->kx, sizeof(double) * nx);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML X coefficients (kx)\n");
            return 0;
        }
        cpml->kx = new_ptr;

        new_ptr = (double*)realloc(cpml->bx, sizeof(double) * nx);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML X coefficients (bx)\n");
            return 0;
        }
        cpml->bx = new_ptr;

        new_ptr = (double*)realloc(cpml->cx, sizeof(double) * nx);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML X coefficients (cx)\n");
            return 0;
        }
        cpml->cx = new_ptr;

        cpml->kx_capacity = nx;
    }
    if (ny > cpml->ky_capacity) {
        double* new_ptr = (double*)realloc(cpml->ky, sizeof(double) * ny);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML Y coefficients (ky)\n");
            return 0;
        }
        cpml->ky = new_ptr;

        new_ptr = (double*)realloc(cpml->by, sizeof(double) * ny);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML Y coefficients (by)\n");
            return 0;
        }
        cpml->by = new_ptr;

        new_ptr = (double*)realloc(cpml->cy, sizeof(double) * ny);
        if (!new_ptr) {
            fprintf(stderr, "Failed to allocate CPML Y coefficients (cy)\n");
            return 0;
        }
        cpml->cy = new_ptr;

        cpml->ky_capacity = ny;
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
void cpml_build_coeffs(SimulationState* state) {
    if (!state) return;
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return;

    double dt = state->dt;
    int nx = state->nx;
    int ny = state->ny;

    if (!ensure_cpml_capacity(state, nx, ny)) {
        return;
    }

    const double m = 3.0;  /* cubic grading */

    /* X-direction coefficients */
    for (int i = 0; i < nx; i++) {
        double rx = 0.0;

        /* Compute fractional distance into PML region */
        if (i < cpml->thickness) {
            rx = (cpml->thickness - i) / (double)cpml->thickness;
        } else if (i >= nx - cpml->thickness) {
            rx = (i - (nx - cpml->thickness)) / (double)cpml->thickness;  /* FIXED: was -1 */
        }

        double g = rx > 0 ? pow(rx, m) : 0.0;
        double sigma = cpml->sigma_max * g;
        double kappa = 1.0 + (cpml->kappa_max - 1.0) * g;
        double alpha = cpml->alpha_max * (1.0 - rx);

        cpml->kx[i] = (rx > 0) ? kappa : 1.0;

        if (rx > 0) {
            cpml->bx[i] = exp(-(sigma/kappa + alpha) * dt);
            cpml->cx[i] = (sigma * (cpml->bx[i] - 1.0)) /
                          (kappa * (sigma + kappa*alpha) + DIVISION_SAFETY_EPSILON);
        } else {
            cpml->bx[i] = 1.0;
            cpml->cx[i] = 0.0;
        }
    }

    /* Y-direction coefficients */
    for (int j = 0; j < ny; j++) {
        double ry = 0.0;

        if (j < cpml->thickness) {
            ry = (cpml->thickness - j) / (double)cpml->thickness;
        } else if (j >= ny - cpml->thickness) {
            ry = (j - (ny - cpml->thickness)) / (double)cpml->thickness;  /* FIXED: was -1 */
        }

        double g = ry > 0 ? pow(ry, m) : 0.0;
        double sigma = cpml->sigma_max * g;
        double kappa = 1.0 + (cpml->kappa_max - 1.0) * g;
        double alpha = cpml->alpha_max * (1.0 - ry);

        cpml->ky[j] = (ry > 0) ? kappa : 1.0;

        if (ry > 0) {
            cpml->by[j] = exp(-(sigma/kappa + alpha) * dt);
            cpml->cy[j] = (sigma * (cpml->by[j] - 1.0)) /
                          (kappa * (sigma + kappa*alpha) + DIVISION_SAFETY_EPSILON);
        } else {
            cpml->by[j] = 1.0;
            cpml->cy[j] = 0.0;
        }
    }
}

/* Apply a CPML preset */
void cpml_apply_preset(SimulationState* state, int idx) {
    if (!state) return;
    CpmlState* cpml = cpml_state(state);
    if (!cpml) return;
    int n = sizeof(CPML_PRESETS) / sizeof(CPML_PRESETS[0]);
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;

    cpml->preset_idx = idx;
    cpml->sigma_max = CPML_PRESETS[idx].smax;
    cpml->kappa_max = CPML_PRESETS[idx].kmax;
    cpml->alpha_max = CPML_PRESETS[idx].amax;
    cpml->thickness = CPML_PRESETS[idx].thick;

    /* Validate CPML thickness doesn't exceed grid size */
    int min_dim = (state->nx < state->ny ? state->nx : state->ny);
    int max_thickness = min_dim / 2 - 1;
    if (cpml->thickness > max_thickness) {
        fprintf(stderr, "Warning: CPML thickness %d exceeds safe limit %d, clamping\n",
                cpml->thickness, max_thickness);
        cpml->thickness = max_thickness;
    }
    if (cpml->thickness < 1) cpml->thickness = 1;

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

    int i;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (i = 1; i < nx - 1; i++) {
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
    int j;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (j = 1; j < ny - 1; j++) {
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
