// ============================================================================
// emwave-c: FDTD Core Simulation Engine
// Electromagnetic field updates using Yee grid
// ============================================================================

#include "fdtd_core.h"
#include "config.h"
#include "boundary.h"
#include "sources.h"
#include "materials.h"
#include "analysis.h"
#include "config_loader.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int alloc_field_double(int nx, int ny, double*** out_rows, double** out_data) {
    size_t cells = (size_t)nx * (size_t)ny;
    double* data = (double*)calloc(cells, sizeof(double));
    if (!data) return 0;
    double** rows = (double**)malloc(sizeof(double*) * nx);
    if (!rows) {
        free(data);
        return 0;
    }
    for (int i = 0; i < nx; i++) {
        rows[i] = data + (size_t)i * ny;
    }
    *out_rows = rows;
    *out_data = data;
    return 1;
}

static int alloc_field_uchar(int nx, int ny, unsigned char*** out_rows, unsigned char** out_data) {
    size_t cells = (size_t)nx * (size_t)ny;
    unsigned char* data = (unsigned char*)calloc(cells, sizeof(unsigned char));
    if (!data) return 0;
    unsigned char** rows = (unsigned char**)malloc(sizeof(unsigned char*) * nx);
    if (!rows) {
        free(data);
        return 0;
    }
    for (int i = 0; i < nx; i++) {
        rows[i] = data + (size_t)i * ny;
    }
    *out_rows = rows;
    *out_data = data;
    return 1;
}

static void free_field_double(double*** rows, double** data) {
    if (*rows) {
        free(*rows);
        *rows = NULL;
    }
    if (*data) {
        free(*data);
        *data = NULL;
    }
}

static void free_field_uchar(unsigned char*** rows, unsigned char** data) {
    if (*rows) {
        free(*rows);
        *rows = NULL;
    }
    if (*data) {
        free(*data);
        *data = NULL;
    }
}

/* CFL timestep calculation */
double fdtd_compute_dt(double dx, double dy, double cfl_safety) {
    double safety = (cfl_safety > 0.0) ? cfl_safety : CFL_SAFETY_FACTOR;
    return safety / (c0 * sqrt(1.0/(dx*dx) + 1.0/(dy*dy)));
}

/* Initialize simulation state */
SimulationState* fdtd_init(const SimulationConfig* cfg) {
    SimulationState* state = (SimulationState*)calloc(1, sizeof(SimulationState));
    if (!state) {
        fprintf(stderr, "Failed to allocate SimulationState\n");
        return NULL;
    }

    SimulationConfig local_cfg = cfg ? *cfg : SIM_CONFIG_DEFAULTS;
    config_clamp_to_limits(&local_cfg);
    char errbuf[128];
    if (!config_validate(&local_cfg, errbuf, sizeof(errbuf))) {
        fprintf(stderr, "Invalid configuration: %s\n", errbuf);
        free(state);
        return NULL;
    }

    state->config = local_cfg;
    state->nx = local_cfg.nx;
    state->ny = local_cfg.ny;
    state->lx = local_cfg.lx;
    state->ly = local_cfg.ly;
    state->cfl_safety = local_cfg.cfl_safety;
    state->sweep_points = local_cfg.sweep_points;
    state->sweep_start_hz = local_cfg.sweep_start_hz;
    state->sweep_stop_hz = local_cfg.sweep_stop_hz;
    state->sweep_steps_per_point = local_cfg.sweep_steps_per_point;

    state->dx = state->lx / (double)state->nx;
    state->dy = state->ly / (double)state->ny;
    state->dt = fdtd_compute_dt(state->dx, state->dy, state->cfl_safety);

    state->timestep = 0;
    state->freq = 1e9;

    boundary_init(state);

    int nx = state->nx;
    int ny = state->ny;

    if (!alloc_field_double(nx, ny, &state->Ez, &state->Ez_data) ||
        !alloc_field_double(nx, ny, &state->Hx, &state->Hx_data) ||
        !alloc_field_double(nx, ny, &state->Hy, &state->Hy_data) ||
        !alloc_field_double(nx, ny, &state->Ez_old, &state->Ez_old_data) ||
        !alloc_field_double(nx, ny, &state->psi_Ezx, &state->psi_Ezx_data) ||
        !alloc_field_double(nx, ny, &state->psi_Ezy, &state->psi_Ezy_data) ||
        !alloc_field_double(nx, ny, &state->psi_Hyx, &state->psi_Hyx_data) ||
        !alloc_field_double(nx, ny, &state->psi_Hxy, &state->psi_Hxy_data) ||
        !alloc_field_double(nx, ny, &state->epsr, &state->epsr_data) ||
        !alloc_field_double(nx, ny, &state->sigma_map, &state->sigma_map_data) ||
        !alloc_field_uchar(nx, ny, &state->tag_grid, &state->tag_grid_data)) {
        fprintf(stderr, "Failed to allocate simulation grids\n");
        fdtd_free(state);
        return NULL;
    }

    state->ports_on = 0;
    for (int p = 0; p < MAX_PORTS; p++) {
        state->ports[p].active = 0;
        state->ports[p].V = NULL;
        state->ports[p].I = NULL;
    }
    ports_init(state->ports, state->nx, state->ny);

    fdtd_clear_fields(state);
    materials_reset_to_defaults(state);
    materials_init(state);
    sources_init(state->sources, state->nx, state->ny, &state->config);
    cpml_apply_preset(state, cpml_get_preset_index(state));
    cpml_zero_psi(state);

    state->step_Ez_absmax = 0.0;

    return state;
}

/* Free simulation state */
void fdtd_free(SimulationState* state) {
    if (!state) return;

    free_field_double(&state->Ez, &state->Ez_data);
    free_field_double(&state->Hx, &state->Hx_data);
    free_field_double(&state->Hy, &state->Hy_data);
    free_field_double(&state->Ez_old, &state->Ez_old_data);
    free_field_double(&state->psi_Ezx, &state->psi_Ezx_data);
    free_field_double(&state->psi_Ezy, &state->psi_Ezy_data);
    free_field_double(&state->psi_Hyx, &state->psi_Hyx_data);
    free_field_double(&state->psi_Hxy, &state->psi_Hxy_data);
    free_field_double(&state->epsr, &state->epsr_data);
    free_field_double(&state->sigma_map, &state->sigma_map_data);
    free_field_uchar(&state->tag_grid, &state->tag_grid_data);

    ports_free(state->ports);
    boundary_shutdown(state);

    free(state);
}

/* Clear all fields to zero */
void fdtd_clear_fields(SimulationState* state) {
    if (!state) return;
    size_t cells = (size_t)state->nx * (size_t)state->ny;
    if (state->Ez_data) memset(state->Ez_data, 0, sizeof(double) * cells);
    if (state->Hx_data) memset(state->Hx_data, 0, sizeof(double) * cells);
    if (state->Hy_data) memset(state->Hy_data, 0, sizeof(double) * cells);
    if (state->Ez_old_data) memset(state->Ez_old_data, 0, sizeof(double) * cells);
    if (state->psi_Ezx_data) memset(state->psi_Ezx_data, 0, sizeof(double) * cells);
    if (state->psi_Ezy_data) memset(state->psi_Ezy_data, 0, sizeof(double) * cells);
    if (state->psi_Hyx_data) memset(state->psi_Hyx_data, 0, sizeof(double) * cells);
    if (state->psi_Hxy_data) memset(state->psi_Hxy_data, 0, sizeof(double) * cells);
    state->step_Ez_absmax = 0.0;
}

/* Reset simulation state */
void fdtd_reset(SimulationState* state) {
    fdtd_clear_fields(state);
    state->timestep = 0;
    cpml_zero_psi(state);
    state->step_Ez_absmax = 0.0;

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
    state->dt = fdtd_compute_dt(state->dx, state->dy, state->cfl_safety);
    state->freq = freq;
    cpml_build_coeffs(state);
    cpml_zero_psi(state);
}

/* Main FDTD step - Update electromagnetic fields by one timestep */
void fdtd_step(SimulationState* state) {
    double dx = state->dx;
    double dy = state->dy;
    double dt = state->dt;
    int nx = state->nx;
    int ny = state->ny;
    if (nx <= 1 || ny <= 1) return;

    CpmlState* cpml = &state->cpml;
    int cpml_active = cpml ? cpml->enabled : 0;
    int cpml_N = cpml ? cpml->thickness : 0;
    double* kx = cpml ? cpml->kx : NULL;
    double* bx = cpml ? cpml->bx : NULL;
    double* cx = cpml ? cpml->cx : NULL;
    double* ky = cpml ? cpml->ky : NULL;
    double* by = cpml ? cpml->by : NULL;
    double* cy = cpml ? cpml->cy : NULL;

    /* Save Ez for Mur-1 boundaries */
    int i, j;
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            state->Ez_old[i][j] = state->Ez[i][j];
        }
    }

    /* Update H fields with optional CPML */
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for (i = 0; i < nx-1; i++) {
        for (j = 0; j < ny-1; j++) {
            /* dEz/dy -> Hx */
            double dEdy = (state->Ez[i][j+1] - state->Ez[i][j]) / dy;
            if (cpml_active && (j < cpml_N || j >= ny - cpml_N)) {
                state->psi_Ezy[i][j] = by[j] * state->psi_Ezy[i][j] + cy[j] * dEdy;
                dEdy = (dEdy / ky[j]) + state->psi_Ezy[i][j];
            }
            state->Hx[i][j] -= (dt / MU0) * dEdy;

            /* dEz/dx -> Hy */
            double dEdx = (state->Ez[i+1][j] - state->Ez[i][j]) / dx;
            if (cpml_active && (i < cpml_N || i >= nx - cpml_N)) {
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
    double ez_absmax = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for private(j) reduction(max:ez_absmax)
    #endif
    for (i = 1; i < nx; i++) {
        for (j = 1; j < ny; j++) {
            /* dHy/dx */
            double dHdx = (state->Hy[i][j] - state->Hy[i-1][j]) / dx;
            if (cpml_active && (i < cpml_N || i >= nx - cpml_N)) {
                state->psi_Hyx[i][j] = bx[i] * state->psi_Hyx[i][j] + cx[i] * dHdx;
                dHdx = (dHdx / kx[i]) + state->psi_Hyx[i][j];
            }

            /* dHx/dy */
            double dHdy = (state->Hx[i][j] - state->Hx[i][j-1]) / dy;
            if (cpml_active && (j < cpml_N || j >= ny - cpml_N)) {
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

            double abs_val = fabs(state->Ez[i][j]);
            if (abs_val > ez_absmax) {
                ez_absmax = abs_val;
            }
        }
    }

    /* Apply Mur boundaries if CPML is off */
    if (!cpml_active) {
        apply_mur_boundaries(state);
    }

    /* Inject all active sources */
    inject_all_sources(state);

    double src_absmax = ez_absmax;
    for (int s = 0; s < MAX_SRC; s++) {
        if (!state->sources[s].active) continue;
        for (int di = -2; di <= 2; di++) {
            for (int dj = -2; dj <= 2; dj++) {
                int i_src = state->sources[s].ix + di;
                int j_src = state->sources[s].iy + dj;
                if (i_src > 0 && i_src < state->nx && j_src > 0 && j_src < state->ny) {
                    double abs_val = fabs(state->Ez[i_src][j_src]);
                    if (abs_val > src_absmax) {
                        src_absmax = abs_val;
                    }
                }
            }
        }
    }

    /* Increment timestep counter */
    state->timestep++;
    state->step_Ez_absmax = src_absmax;
}
