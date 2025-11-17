// ============================================================================
// emwave-c: Material Property Management
// ============================================================================

#ifndef EMWAVE_MATERIALS_H
#define EMWAVE_MATERIALS_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Material initialization */
void materials_init(SimulationState* state);
void materials_reset_to_defaults(SimulationState* state);

/* Material painting (interactive editing) */
void paint_material_at(SimulationState* state, int i, int j, int paint_type, double paint_eps);

/* Material query helpers (inline for performance) */
static inline int is_pec(const SimulationState* state, int i, int j) {
    return state->tag_grid[i][j] == 1;
}

static inline int is_pmc(const SimulationState* state, int i, int j) {
    return state->tag_grid[i][j] == 2;
}

static inline int is_dielectric(const SimulationState* state, int i, int j) {
    return state->tag_grid[i][j] == 0;
}

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_MATERIALS_H */
