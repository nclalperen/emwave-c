#include "material_library.h"

#include <stddef.h>
#include <string.h>

#define MATERIAL_LIBRARY_MAX_RESULTS 64

// Internal material database
static Material g_materials[] = {
    // === METALS (PEC approximation for RF/microwave) ===
    {
        .name = "Copper",
        .description = "High conductivity metal, PCB traces",
        .category = MAT_CAT_METAL,
        .type = MAT_TYPE_PEC,
        .epsilon_r = 1.0,
        .mu_r = 1.0,
        .tan_delta = 0.0,
        .conductivity = 5.96e7,  // S/m at 20°C
        .frequency_dependent = false,
        .color_r = 184,
        .color_g = 115,
        .color_b = 51,
        .id = 1
    },
    {
        .name = "Gold",
        .description = "Corrosion-resistant contact metal",
        .category = MAT_CAT_METAL,
        .type = MAT_TYPE_PEC,
        .epsilon_r = 1.0,
        .mu_r = 1.0,
        .tan_delta = 0.0,
        .conductivity = 4.1e7,
        .frequency_dependent = false,
        .color_r = 255,
        .color_g = 215,
        .color_b = 0,
        .id = 2
    },
    {
        .name = "Silver",
        .description = "Highest conductivity metal",
        .category = MAT_CAT_METAL,
        .type = MAT_TYPE_PEC,
        .epsilon_r = 1.0,
        .mu_r = 1.0,
        .tan_delta = 0.0,
        .conductivity = 6.3e7,
        .frequency_dependent = false,
        .color_r = 192,
        .color_g = 192,
        .color_b = 192,
        .id = 3
    },
    {
        .name = "Aluminum",
        .description = "Lightweight RF shielding",
        .category = MAT_CAT_METAL,
        .type = MAT_TYPE_PEC,
        .epsilon_r = 1.0,
        .mu_r = 1.0,
        .tan_delta = 0.0,
        .conductivity = 3.5e7,
        .frequency_dependent = false,
        .color_r = 211,
        .color_g = 211,
        .color_b = 211,
        .id = 4
    },

    // === DIELECTRICS ===
    {
        .name = "Air",
        .description = "Vacuum approximation",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_DIELECTRIC,
        .epsilon_r = 1.0,
        .mu_r = 1.0,
        .tan_delta = 0.0,
        .conductivity = 0.0,
        .frequency_dependent = false,
        .color_r = 200,
        .color_g = 220,
        .color_b = 255,
        .id = 10
    },
    {
        .name = "FR4",
        .description = "Standard PCB substrate",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_LOSSY_DIEL,
        .epsilon_r = 4.4,
        .mu_r = 1.0,
        .tan_delta = 0.02,
        .conductivity = 0.0,
        .frequency_dependent = true,
        .color_r = 34,
        .color_g = 139,
        .color_b = 34,
        .id = 11
    },
    {
        .name = "Rogers RO4003",
        .description = "Low-loss RF substrate",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_LOSSY_DIEL,
        .epsilon_r = 3.55,
        .mu_r = 1.0,
        .tan_delta = 0.0027,
        .conductivity = 0.0,
        .frequency_dependent = true,
        .color_r = 100,
        .color_g = 149,
        .color_b = 237,
        .id = 12
    },
    {
        .name = "Teflon (PTFE)",
        .description = "Ultra-low loss insulator",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_DIELECTRIC,
        .epsilon_r = 2.1,
        .mu_r = 1.0,
        .tan_delta = 0.0002,
        .conductivity = 0.0,
        .frequency_dependent = false,
        .color_r = 240,
        .color_g = 255,
        .color_b = 255,
        .id = 13
    },
    {
        .name = "Silicon",
        .description = "Semiconductor substrate",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_LOSSY_DIEL,
        .epsilon_r = 11.68,
        .mu_r = 1.0,
        .tan_delta = 0.001,
        .conductivity = 1e-12,
        .frequency_dependent = false,
        .color_r = 169,
        .color_g = 169,
        .color_b = 169,
        .id = 14
    },
    {
        .name = "Glass",
        .description = "Optical fiber, windows",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_DIELECTRIC,
        .epsilon_r = 6.0,
        .mu_r = 1.0,
        .tan_delta = 0.0001,
        .conductivity = 0.0,
        .frequency_dependent = false,
        .color_r = 173,
        .color_g = 216,
        .color_b = 230,
        .id = 15
    },
    {
        .name = "Alumina (Al2O3)",
        .description = "High-temperature ceramic substrate",
        .category = MAT_CAT_DIELECTRIC,
        .type = MAT_TYPE_DIELECTRIC,
        .epsilon_r = 9.8,
        .mu_r = 1.0,
        .tan_delta = 0.0001,
        .conductivity = 0.0,
        .frequency_dependent = false,
        .color_r = 255,
        .color_g = 228,
        .color_b = 225,
        .id = 16
    }
};

static const int g_material_count = (int)(sizeof(g_materials) / sizeof(g_materials[0]));
static bool g_initialized = false;

static void ensure_initialized(void) {
    if (!g_initialized) {
        material_library_init();
    }
}

void material_library_init(void) {
    g_initialized = true;
}

void material_library_shutdown(void) {
    g_initialized = false;
}

int material_library_get_count(void) {
    ensure_initialized();
    return g_material_count;
}

const Material* material_library_get_by_index(int index) {
    ensure_initialized();
    if (index < 0 || index >= g_material_count) {
        return NULL;
    }
    return &g_materials[index];
}

const Material* material_library_get_by_name(const char* name) {
    ensure_initialized();
    if (!name) return NULL;
    for (int i = 0; i < g_material_count; ++i) {
        if (strcmp(g_materials[i].name, name) == 0) {
            return &g_materials[i];
        }
    }
    return NULL;
}

const Material* material_library_get_by_id(int id) {
    ensure_initialized();
    for (int i = 0; i < g_material_count; ++i) {
        if (g_materials[i].id == id) {
            return &g_materials[i];
        }
    }
    return NULL;
}

int material_library_get_category_count(MaterialCategory cat) {
    ensure_initialized();
    int count = 0;
    for (int i = 0; i < g_material_count; ++i) {
        if (g_materials[i].category == cat) {
            ++count;
        }
    }
    return count;
}

const Material** material_library_get_by_category(MaterialCategory cat, int* out_count) {
    ensure_initialized();
    static const Material* results[MATERIAL_LIBRARY_MAX_RESULTS];
    int count = 0;
    for (int i = 0; i < g_material_count && count < MATERIAL_LIBRARY_MAX_RESULTS; ++i) {
        if (g_materials[i].category == cat) {
            results[count++] = &g_materials[i];
        }
    }
    if (out_count) {
        *out_count = count;
    }
    return results;
}

double material_get_epsilon(const Material* mat) {
    return (mat) ? mat->epsilon_r : 1.0;
}

double material_get_mu(const Material* mat) {
    return (mat) ? mat->mu_r : 1.0;
}

double material_get_sigma(const Material* mat) {
    return (mat) ? mat->conductivity : 0.0;
}

