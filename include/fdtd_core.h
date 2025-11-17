// ============================================================================
// emwave-c: FDTD Core Engine Interface
// ============================================================================

#ifndef EMWAVE_FDTD_CORE_H
#define EMWAVE_FDTD_CORE_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Initialization and cleanup */
SimulationState* fdtd_init(const SimulationConfig* cfg);
void fdtd_free(SimulationState* state);

/* Testing hooks */
void fdtd_test_set_alloc_fail_after(int count);

/* Simulation control */
void fdtd_step(SimulationState* state);
void fdtd_reset(SimulationState* state);
void fdtd_clear_fields(SimulationState* state);

/* Grid and timestep management */
double fdtd_compute_dt(double dx, double dy, double cfl_safety);
void fdtd_update_grid_for_freq(SimulationState* state, double freq);

/* Material access */
double fdtd_epsilon_at(const SimulationState* state, int i, int j);
double fdtd_sigma_at(const SimulationState* state, int i, int j);

/* Field access helpers */
#ifdef EMWAVE_BOUNDS_CHECK
static inline double fdtd_get_Ez(const SimulationState* state, int i, int j) {
    if (!state || !state->Ez) return 0.0;
    if (i < 0 || i >= state->nx || j < 0 || j >= state->ny) return 0.0;
    return state->Ez[i][j];
}

static inline double fdtd_get_Hx(const SimulationState* state, int i, int j) {
    if (!state || !state->Hx) return 0.0;
    if (i < 0 || i >= state->nx || j < 0 || j >= state->ny) return 0.0;
    return state->Hx[i][j];
}

static inline double fdtd_get_Hy(const SimulationState* state, int i, int j) {
    if (!state || !state->Hy) return 0.0;
    if (i < 0 || i >= state->nx || j < 0 || j >= state->ny) return 0.0;
    return state->Hy[i][j];
}
#else
static inline double fdtd_get_Ez(const SimulationState* state, int i, int j) {
    return state->Ez[i][j];
}

static inline double fdtd_get_Hx(const SimulationState* state, int i, int j) {
    return state->Hx[i][j];
}

static inline double fdtd_get_Hy(const SimulationState* state, int i, int j) {
    return state->Hy[i][j];
}
#endif

/* Utility functions */
static inline int clampi(int v, int lo, int hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

static inline double clampd(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

static inline float sqrf(float x) {
    return x * x;
}

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_FDTD_CORE_H */
