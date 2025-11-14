// ============================================================================
// emwave-c: Boundary Conditions (CPML and Mur)
// ============================================================================

#ifndef EMWAVE_BOUNDARY_H
#define EMWAVE_BOUNDARY_H

#include "types.h"

/* CPML presets */
extern const CpmlPreset CPML_PRESETS[3];

/* Boundary lifecycle */
void boundary_init(SimulationState* state);
void boundary_shutdown(SimulationState* state);

/* Boundary configuration */
void boundary_set_type(SimulationState* state, BoundaryType type);
BoundaryType boundary_get_type(const SimulationState* state);
int boundary_is_cpml_enabled(const SimulationState* state);

/* CPML functions */
void cpml_zero_psi(SimulationState* state);
void cpml_build_coeffs(SimulationState* state);
void cpml_apply_preset(SimulationState* state, int idx);
int cpml_get_preset_index(const SimulationState* state);

/* Mur boundary application */
void apply_mur_boundaries(SimulationState* state);

#endif /* EMWAVE_BOUNDARY_H */
