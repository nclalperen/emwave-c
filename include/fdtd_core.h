// ============================================================================
// emwave-c: FDTD Core Engine Interface
// ============================================================================

#ifndef EMWAVE_FDTD_CORE_H
#define EMWAVE_FDTD_CORE_H

#include "types.h"

/* Initialization and cleanup */
SimulationState* fdtd_init(const SimulationConfig* cfg);
void fdtd_free(SimulationState* state);

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
static inline double fdtd_get_Ez(const SimulationState* state, int i, int j) {
    return state->Ez[i][j];
}

static inline double fdtd_get_Hx(const SimulationState* state, int i, int j) {
    return state->Hx[i][j];
}

static inline double fdtd_get_Hy(const SimulationState* state, int i, int j) {
    return state->Hy[i][j];
}

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

#endif /* EMWAVE_FDTD_CORE_H */
