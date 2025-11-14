#include "config_loader.h"
#include "config.h"

#include <assert.h>
#include <stdio.h>

static void test_waveguide_config(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/waveguide.json", &cfg);
    assert(ok && "waveguide.json should parse");
    assert(cfg.nx == 512);
    assert(cfg.ny == 256);
    assert(cfg.material_rect_count == 3);
    assert(cfg.source_count == 1);
    assert(cfg.source_configs[0].type == SRC_GAUSS_PULSE);
}

static void test_cpw_config(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/cpw_filter.json", &cfg);
    assert(ok && "cpw_filter.json should parse");
    assert(cfg.nx == 400);
    assert(cfg.ny == 400);
    assert(cfg.material_rect_count == 4);
    assert(cfg.source_count == 2);
    assert(cfg.source_configs[1].active == 0);
}

static void test_invalid_config_file(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/invalid_config.json", &cfg);
    assert(!ok && "invalid_config.json should fail to parse");
    assert(cfg.nx == SIM_CONFIG_DEFAULTS.nx);
    assert(cfg.material_rect_count == SIM_CONFIG_DEFAULTS.material_rect_count);
}

static void test_missing_file(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/does_not_exist.json", &cfg);
    assert(!ok && "Missing file should fail to parse");
}

int main(void) {
    test_waveguide_config();
    test_cpw_config();
    test_invalid_config_file();
    test_missing_file();
    puts("config_loader tests passed");
    return 0;
}
