#include "core/material_library.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

static void test_material_library_init(void) {
    material_library_init();
    assert(material_library_get_count() > 0);
    printf("✓ Material library initialized\n");
}

static void test_get_by_name(void) {
    const Material* mat = material_library_get_by_name("FR4");
    assert(mat != NULL);
    assert(strcmp(mat->name, "FR4") == 0);
    assert(mat->epsilon_r == 4.4);
    assert(mat->category == MAT_CAT_DIELECTRIC);
    printf("✓ Get material by name: FR4\n");
}

static void test_get_by_id(void) {
    const Material* mat = material_library_get_by_id(1);  // Copper
    assert(mat != NULL);
    assert(strcmp(mat->name, "Copper") == 0);
    assert(mat->type == MAT_TYPE_PEC);
    printf("✓ Get material by ID: Copper\n");
}

static void test_category_filtering(void) {
    int count = 0;
    const Material** metals = material_library_get_by_category(MAT_CAT_METAL, &count);
    (void)metals;
    assert(count == 4);
    printf("✓ Category filtering: %d metals found\n", count);
}

int main(void) {
    printf("Running material library tests...\n");
    test_material_library_init();
    test_get_by_name();
    test_get_by_id();
    test_category_filtering();
    material_library_shutdown();
    printf("All tests passed!\n");
    return 0;
}

