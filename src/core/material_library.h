#ifndef MATERIAL_LIBRARY_H
#define MATERIAL_LIBRARY_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Material categories
typedef enum {
    MAT_CAT_METAL,       // Perfect conductors (PEC)
    MAT_CAT_DIELECTRIC,  // Insulators (ε_r > 1)
    MAT_CAT_MAGNETIC,    // Magnetic materials (μ_r ≠ 1) [future]
    MAT_CAT_CUSTOM       // User-defined
} MaterialCategory;

// Material type (determines physics model)
typedef enum {
    MAT_TYPE_PEC,        // Perfect electric conductor (σ → ∞)
    MAT_TYPE_PMC,        // Perfect magnetic conductor (boundary condition)
    MAT_TYPE_DIELECTRIC, // Lossless dielectric (ε_r, tan(δ) = 0)
    MAT_TYPE_LOSSY_DIEL  // Lossy dielectric (ε_r, tan(δ) > 0)
} MaterialType;

// Material definition
typedef struct {
    // Identification
    const char* name;         // "FR4", "Copper", etc.
    const char* description;  // "Flame retardant PCB substrate"
    MaterialCategory category;
    MaterialType type;

    // Electromagnetic properties
    double epsilon_r;    // Relative permittivity (ε_r)
    double mu_r;         // Relative permeability (μ_r, default 1.0)
    double tan_delta;    // Loss tangent (tan(δ))
    double conductivity; // σ (S/m) for metals

    // Frequency-dependent flag (for future extensions)
    bool frequency_dependent;

    // Color for visualization (RGB, 0-255)
    unsigned char color_r;
    unsigned char color_g;
    unsigned char color_b;

    // Internal ID
    int id;
} Material;

// Material library API
void material_library_init(void);
void material_library_shutdown(void);

int material_library_get_count(void);
const Material* material_library_get_by_index(int index);
const Material* material_library_get_by_name(const char* name);
const Material* material_library_get_by_id(int id);

// Material queries
int material_library_get_category_count(MaterialCategory cat);
const Material** material_library_get_by_category(MaterialCategory cat, int* count);

// Convert material to simulation parameters
double material_get_epsilon(const Material* mat);
double material_get_mu(const Material* mat);
double material_get_sigma(const Material* mat);

#ifdef __cplusplus
}
#endif

#endif // MATERIAL_LIBRARY_H
