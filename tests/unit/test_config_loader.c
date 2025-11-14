#include "config_loader.h"
#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void test_waveguide_config(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/waveguide.json", &cfg, NULL, 0);
    assert(ok && "waveguide.json should parse");
    assert(cfg.nx == 512);
    assert(cfg.ny == 256);
    assert(cfg.material_rect_count == 3);
    assert(cfg.source_count == 1);
    assert(cfg.source_configs[0].type == SRC_GAUSS_PULSE);
}

static void test_cpw_config(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/cpw_filter.json", &cfg, NULL, 0);
    assert(ok && "cpw_filter.json should parse");
    assert(cfg.nx == 400);
    assert(cfg.ny == 400);
    assert(cfg.material_rect_count == 4);
    assert(cfg.source_count == 2);
    assert(cfg.source_configs[1].active == 0);
}

static void test_invalid_config_file(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/invalid_config.json", &cfg, NULL, 0);
    assert(!ok && "invalid_config.json should fail to parse");
    assert(cfg.nx == SIM_CONFIG_DEFAULTS.nx);
    assert(cfg.material_rect_count == SIM_CONFIG_DEFAULTS.material_rect_count);
}

static void test_missing_file(void) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(CONFIGS_DIR "/does_not_exist.json", &cfg, errbuf, sizeof(errbuf));
    assert(!ok && "Missing file should fail to parse");
    assert(errbuf[0] != '\0');
}

static void test_oversized_file_guard(void) {
    char path[] = "/tmp/emwave_config_oversizeXXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);

    FILE* f = fdopen(fd, "wb");
    assert(f && "fdopen failed");

    size_t target = (size_t)CONFIG_LOADER_MAX_FILE_BYTES + 1024;
    char block[4096];
    memset(block, 'x', sizeof(block));
    size_t remaining = target;
    while (remaining > 0) {
        size_t chunk = remaining < sizeof(block) ? remaining : sizeof(block);
        size_t written = fwrite(block, 1, chunk, f);
        assert(written == chunk);
        remaining -= chunk;
    }
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    assert(!ok && "Oversized config file should be rejected");
    assert(errbuf[0] != '\0');
    assert(strstr(errbuf, "too large") != NULL);

    unlink(path);
}

int main(void) {
    test_waveguide_config();
    test_cpw_config();
    test_invalid_config_file();
    test_missing_file();
    test_oversized_file_guard();
    puts("config_loader tests passed");
    return 0;
}
