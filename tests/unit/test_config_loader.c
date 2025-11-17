#include "config_loader.h"
#include "config.h"

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <unistd.h>
#include <errno.h>
#endif

static int create_temp_config(char* path, size_t path_len, FILE** out_file) {
#if defined(_WIN32)
    char tmp_dir[MAX_PATH];
    DWORD dir_len = GetTempPathA(MAX_PATH, tmp_dir);
    if (dir_len == 0 || dir_len >= MAX_PATH) {
        return 0;
    }
    char tmp_file[MAX_PATH];
    if (GetTempFileNameA(tmp_dir, "ewc", 0, tmp_file) == 0) {
        return 0;
    }
    if (strlen(tmp_file) + 1 > path_len) {
        return 0;
    }
    strncpy(path, tmp_file, path_len);
    path[path_len - 1] = '\0';
    FILE* file = fopen(path, "wb");
    if (!file) {
        return 0;
    }
#else
    char tmpl[512];
#ifdef TEST_TMP_DIR
    int written = snprintf(tmpl, sizeof(tmpl), "%s/emwave_config_XXXXXX", TEST_TMP_DIR);
#else
    int written = snprintf(tmpl, sizeof(tmpl), "/tmp/emwave_config_XXXXXX");
#endif
    if (written <= 0 || (size_t)written >= sizeof(tmpl)) {
        return 0;
    }
    int fd = mkstemp(tmpl);
    if (fd < 0) {
        return 0;
    }
    if (strlen(tmpl) + 1 > path_len) {
        close(fd);
        errno = ENAMETOOLONG;
        return 0;
    }
    strncpy(path, tmpl, path_len);
    path[path_len - 1] = '\0';
    FILE* file = fdopen(fd, "wb");
    if (!file) {
        close(fd);
        return 0;
    }
#endif
    if (out_file) {
        *out_file = file;
    } else {
        fclose(file);
    }
    return 1;
}

static void write_repeated_bytes(FILE* f, size_t total_bytes) {
    char block[4096];
    memset(block, 'x', sizeof(block));
    size_t remaining = total_bytes;
    while (remaining > 0) {
        size_t chunk = remaining < sizeof(block) ? remaining : sizeof(block);
        size_t written = fwrite(block, 1, chunk, f);
        ck_assert_msg(written == chunk, "Failed to write temporary config data");
        remaining -= chunk;
    }
    fflush(f);
}

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

START_TEST(test_directory_path_rejected) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(CONFIGS_DIR, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(!ok, "Directories should fail to parse");
    ck_assert_msg(errbuf[0] != '\0', "Directory failure should populate errbuf");
}
END_TEST

START_TEST(test_oversized_file_guard) {
    char path[512];
    FILE* f = NULL;
    ck_assert_msg(create_temp_config(path, sizeof(path), &f), "Failed to create temp config");

    size_t target = (size_t)CONFIG_LOADER_MAX_FILE_BYTES + 1024;
    write_repeated_bytes(f, target);
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(!ok, "Oversized config file should be rejected");
    ck_assert_msg(errbuf[0] != '\0', "Error buffer should describe failure");
    ck_assert_msg(strstr(errbuf, "too large") != NULL, "Expected error to mention size guard");

    remove(path);
}
END_TEST

START_TEST(test_empty_file_is_rejected) {
    char path[512];
    FILE* f = NULL;
    ck_assert_msg(create_temp_config(path, sizeof(path), &f), "Failed to create temp file");
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(!ok, "Empty config should be rejected");
    ck_assert_msg(errbuf[0] != '\0', "Expected error message for empty config");

    remove(path);
}
END_TEST

START_TEST(test_config_truncates_excess_sources) {
    char path[512];
    FILE* f = NULL;
    ck_assert_msg(create_temp_config(path, sizeof(path), &f), "Failed to create temp config");

    fputs("{\n", f);
    fputs("  \"simulation\": { \"nx\": 32, \"ny\": 32 },\n", f);
    fputs("  \"sources\": [", f);
    for (int i = 0; i < MAX_SRC + 1; ++i) {
        fputs("{ \"type\": \"gaussian\", \"x\": 0.1, \"y\": 0.2 }", f);
        if (i != MAX_SRC) {
            fputc(',', f);
        }
    }
    fputs("]\n}\n", f);
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(ok, "Config with too many sources should still parse");
    ck_assert_int_eq(cfg.source_count, MAX_SRC);
    ck_assert_msg(errbuf[0] == '\0', "No error string expected when truncating sources");

    remove(path);
}
END_TEST

START_TEST(test_ports_array_configures_port_specs) {
    char path[512];
    FILE* f = NULL;
    ck_assert_msg(create_temp_config(path, sizeof(path), &f), "Failed to create temp config");

    fputs("{\n", f);
    fputs("  \"simulation\": { \"nx\": 64, \"ny\": 32 },\n", f);
    fputs("  \"ports\": [\n", f);
    fputs("    { \"active\": true, \"x\": 0.25, \"y0\": 0.1, \"y1\": 0.9 }\n", f);
    fputs("  ]\n", f);
    fputs("}\n", f);
    fclose(f);

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    char errbuf[128];
    int ok = config_loader_parse_file(path, &cfg, errbuf, sizeof(errbuf));
    ck_assert_msg(ok, "Config with ports array should parse");
    ck_assert_int_eq(cfg.port_count, 1);
    ck_assert_int_eq(cfg.port_configs[0].active, 1);
    ck_assert_msg(cfg.port_configs[0].x >= 0.0 && cfg.port_configs[0].x <= 1.0,
                  "Port x must be clamped into [0,1]");
    ck_assert_msg(cfg.port_configs[0].y0 >= 0.0 && cfg.port_configs[0].y0 <= 1.0,
                  "Port y0 must be clamped into [0,1]");
    ck_assert_msg(cfg.port_configs[0].y1 >= 0.0 && cfg.port_configs[0].y1 <= 1.0,
                  "Port y1 must be clamped into [0,1]");

    remove(path);
}
END_TEST

START_TEST(test_cli_overrides_json_boundary_and_mode) {
    const char* argv[] = {
        "emwave_cli",
        "--config=" CONFIGS_DIR "/cpw_filter.json",
        "--boundary=mur",
        "--run-mode=fixed",
        "--run-steps=1234"
    };
    int argc = (int)(sizeof(argv) / sizeof(argv[0]));

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_load_from_args(argc, (char**)argv, &cfg);
    ck_assert_msg(ok, "config_load_from_args should succeed");
    ck_assert_int_eq(cfg.boundary_mode, SIM_BOUNDARY_MUR);
    ck_assert_int_eq(cfg.run_mode, SIM_RUN_MODE_FIXED_STEPS);
    ck_assert_int_eq(cfg.run_steps, 1234);
}
END_TEST

START_TEST(test_cli_profile_flag_sets_enable_profile) {
    const char* argv[] = {
        "emwave_cli",
        "--profile"
    };
    int argc = (int)(sizeof(argv) / sizeof(argv[0]));

    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    int ok = config_load_from_args(argc, (char**)argv, &cfg);
    ck_assert_msg(ok, "config_load_from_args should succeed for --profile only");
    ck_assert_int_eq(cfg.enable_profile, 1);
}
END_TEST

static Suite* config_loader_suite(void) {
    Suite* s = suite_create("config_loader");
    TCase* tc = tcase_create("core");
    tcase_add_test(tc, test_waveguide_config);
    tcase_add_test(tc, test_cpw_config);
    tcase_add_test(tc, test_invalid_config_file);
    tcase_add_test(tc, test_missing_file);
    tcase_add_test(tc, test_directory_path_rejected);
    tcase_add_test(tc, test_oversized_file_guard);
    tcase_add_test(tc, test_empty_file_is_rejected);
    tcase_add_test(tc, test_config_truncates_excess_sources);
    tcase_add_test(tc, test_ports_array_configures_port_specs);
    tcase_add_test(tc, test_cli_overrides_json_boundary_and_mode);
    tcase_add_test(tc, test_cli_profile_flag_sets_enable_profile);
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
