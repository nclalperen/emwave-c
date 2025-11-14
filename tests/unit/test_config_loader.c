#include "config_loader.h"
#include "config.h"

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

START_TEST(test_waveguide_config) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/waveguide.json", &cfg, NULL, 0);
    ck_assert_msg(ok, "waveguide.json should parse");
    ck_assert_int_eq(cfg.nx, 512);
    ck_assert_int_eq(cfg.ny, 256);
    ck_assert_int_eq(cfg.material_rect_count, 3);
    ck_assert_int_eq(cfg.source_count, 1);
    ck_assert_int_eq(cfg.source_configs[0].type, SRC_GAUSS_PULSE);
}
END_TEST

START_TEST(test_cpw_config) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/cpw_filter.json", &cfg, NULL, 0);
    ck_assert_msg(ok, "cpw_filter.json should parse");
    ck_assert_int_eq(cfg.nx, 400);
    ck_assert_int_eq(cfg.ny, 400);
    ck_assert_int_eq(cfg.material_rect_count, 4);
    ck_assert_int_eq(cfg.source_count, 2);
    ck_assert_int_eq(cfg.source_configs[1].active, 0);
}
END_TEST

START_TEST(test_invalid_config_file) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_loader_parse_file(CONFIGS_DIR "/invalid_config.json", &cfg, NULL, 0);
    ck_assert_msg(!ok, "invalid_config.json should fail to parse");
    ck_assert_int_eq(cfg.nx, SIM_CONFIG_DEFAULTS.nx);
    ck_assert_int_eq(cfg.material_rect_count, SIM_CONFIG_DEFAULTS.material_rect_count);
}
END_TEST

START_TEST(test_missing_file) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(CONFIGS_DIR "/does_not_exist.json", &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(!ok, "Missing file should fail to parse");
    ck_assert_msg(errbuf[0] != '\0', "Error buffer should be populated");
}
END_TEST

START_TEST(test_oversized_file_guard) {
    char path[] = "/tmp/emwave_config_oversizeXXXXXX";
    int fd = mkstemp(path);
    ck_assert_msg(fd >= 0, "mkstemp failed");

    FILE* f = fdopen(fd, "wb");
    ck_assert_ptr_nonnull(f);

    size_t target = (size_t)CONFIG_LOADER_MAX_FILE_BYTES + 1024;
    char block[4096];
    memset(block, 'x', sizeof(block));
    size_t remaining = target;
    while (remaining > 0) {
        size_t chunk = remaining < sizeof(block) ? remaining : sizeof(block);
        size_t written = fwrite(block, 1, chunk, f);
        ck_assert_msg(written == chunk, "Failed to write oversized block");
        remaining -= chunk;
    }
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(!ok, "Oversized config file should be rejected");
    ck_assert_msg(errbuf[0] != '\0', "Error buffer should describe failure");
    ck_assert_msg(strstr(errbuf, "too large") != NULL, "Expected error to mention size guard");

    unlink(path);
}
END_TEST

static Suite* config_loader_suite(void) {
    Suite* s = suite_create("config_loader");
    TCase* tc = tcase_create("core");
    tcase_add_test(tc, test_waveguide_config);
    tcase_add_test(tc, test_cpw_config);
    tcase_add_test(tc, test_invalid_config_file);
    tcase_add_test(tc, test_missing_file);
    tcase_add_test(tc, test_oversized_file_guard);
    suite_add_tcase(s, tc);
    return s;
}

int main(void) {
    Suite* s = config_loader_suite();
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_ENV);
    int failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (failed == 0) ? 0 : 1;
}
