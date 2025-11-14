// ============================================================================
// emwave-c: Boundary Conditions (CPML and Mur)
// ============================================================================

#ifndef EMWAVE_BOUNDARY_H
#define EMWAVE_BOUNDARY_H

#include "types.h"

/* CPML coefficient arrays (global state) */
extern double* kx;
extern double* bx;
extern double* cx;
extern double* ky;
extern double* by;
extern double* cy;

/* CPML control */
extern int cpml_on;
extern int cpml_N;
extern double pml_sigma_max, pml_kappa_max, pml_alpha_max;
extern int cpml_preset_idx;
extern BoundaryType boundary_type;

/* CPML presets */
extern const CpmlPreset CPML_PRESETS[3];

/* CPML functions */
void cpml_zero_psi(SimulationState* state);
void cpml_build_coeffs(const SimulationState* state);
void cpml_apply_preset(int idx, const SimulationState* state);

/* Mur boundary application */
void apply_mur_boundaries(SimulationState* state);

#endif /* EMWAVE_BOUNDARY_H */
