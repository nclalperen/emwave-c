#include "config_loader.h"
#include "config.h"

#include <check.h>
#include <math.h>
#include <string.h>

START_TEST(test_clamp_enforces_min_max_dims_and_cell_limit) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.nx = SIM_MAX_DIM * 4;
    cfg.ny = SIM_MAX_DIM * 4;

    config_clamp_to_limits(&cfg);

    ck_assert_msg(cfg.nx >= SIM_MIN_DIM,
                  "Expected nx >= SIM_MIN_DIM (%d vs %d)", cfg.nx, SIM_MIN_DIM);
    ck_assert_msg(cfg.ny >= SIM_MIN_DIM,
                  "Expected ny >= SIM_MIN_DIM (%d vs %d)", cfg.ny, SIM_MIN_DIM);
    ck_assert_msg(cfg.nx <= SIM_MAX_DIM,
                  "Expected nx <= SIM_MAX_DIM (%d vs %d)", cfg.nx, SIM_MAX_DIM);
    ck_assert_msg(cfg.ny <= SIM_MAX_DIM,
                  "Expected ny <= SIM_MAX_DIM (%d vs %d)", cfg.ny, SIM_MAX_DIM);

    double cells = (double)cfg.nx * (double)cfg.ny;
    ck_assert_msg(cells <= (double)SIM_MAX_CELLS + 1.0,
                  "Expected cells <= SIM_MAX_CELLS (%.0f vs %.0f)",
                  cells, (double)SIM_MAX_CELLS);
}
END_TEST

START_TEST(test_validate_rejects_invalid_cfl) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.cfl_safety = 1.5;
    char errbuf[128];
    int ok = config_validate(&cfg, errbuf, sizeof(errbuf));
    ck_assert_int_eq(ok, 0);
    ck_assert_msg(errbuf[0] != '\0', "Expected error message for invalid cfl_safety");
}
END_TEST

START_TEST(test_validate_rejects_invalid_sweep_range) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.sweep_start_hz = 1e9;
    cfg.sweep_stop_hz = 5e8;
    char errbuf[128];
    int ok = config_validate(&cfg, errbuf, sizeof(errbuf));
    ck_assert_int_eq(ok, 0);
    ck_assert_msg(errbuf[0] != '\0', "Expected error message for invalid sweep range");
}
END_TEST

START_TEST(test_clamp_limits_material_and_source_counts) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.material_rect_count = CONFIG_MAX_MATERIAL_RECTS + 5;
    cfg.source_count = MAX_SRC + 3;

    config_clamp_to_limits(&cfg);

    ck_assert_int_eq(cfg.material_rect_count, CONFIG_MAX_MATERIAL_RECTS);
    ck_assert_int_eq(cfg.source_count, MAX_SRC);
}
END_TEST

START_TEST(test_clamp_normalizes_rect_coords_and_signs) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.material_rect_count = 1;
    cfg.material_rects[0].x0 = -0.5;
    cfg.material_rects[0].y0 = 1.5;
    cfg.material_rects[0].x1 = 2.0;
    cfg.material_rects[0].y1 = -1.0;
    cfg.material_rects[0].epsr = 2.0;
    cfg.material_rects[0].sigma = -10.0;

    cfg.source_count = 1;
    cfg.source_configs[0].x = -0.1;
    cfg.source_configs[0].y = 1.2;
    cfg.source_configs[0].sigma2 = -1.0;
    cfg.source_configs[0].freq = -100.0;

    config_clamp_to_limits(&cfg);

    const MaterialRectSpec* r = &cfg.material_rects[0];
    ck_assert_msg(r->x0 >= 0.0 && r->x0 <= 1.0, "x0 not clamped into [0,1]");
    ck_assert_msg(r->y0 >= 0.0 && r->y0 <= 1.0, "y0 not clamped into [0,1]");
    ck_assert_msg(r->x1 >= 0.0 && r->x1 <= 1.0, "x1 not clamped into [0,1]");
    ck_assert_msg(r->y1 >= 0.0 && r->y1 <= 1.0, "y1 not clamped into [0,1]");
    ck_assert_msg(r->epsr > 0.0, "epsr should remain positive");
    ck_assert_msg(r->sigma >= 0.0, "sigma should be clamped to non-negative");

    const SourceConfigSpec* s = &cfg.source_configs[0];
    ck_assert_msg(s->x >= 0.0 && s->x <= 1.0, "source x not clamped into [0,1]");
    ck_assert_msg(s->y >= 0.0 && s->y <= 1.0, "source y not clamped into [0,1]");
    ck_assert_msg(s->sigma2 > 0.0, "sigma2 should be clamped to positive default");
    ck_assert_msg(s->freq > 0.0, "freq should be clamped to positive default");
}
END_TEST

START_TEST(test_validate_rejects_invalid_material_rect) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.material_rect_count = 1;
    cfg.material_rects[0].x0 = 0.5;
    cfg.material_rects[0].x1 = 0.4; /* x1 <= x0 */
    cfg.material_rects[0].y0 = 0.2;
    cfg.material_rects[0].y1 = 0.8;
    cfg.material_rects[0].epsr = 2.0;

    char errbuf[128];
    int ok = config_validate(&cfg, errbuf, sizeof(errbuf));
    ck_assert_int_eq(ok, 0);
    ck_assert_msg(errbuf[0] != '\0', "Expected error for invalid material rect bounds");
}
END_TEST

START_TEST(test_validate_rejects_invalid_source_sigma2) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.source_count = 1;
    cfg.source_configs[0].freq = 1e9;
    cfg.source_configs[0].sigma2 = -1.0;

    char errbuf[128];
    int ok = config_validate(&cfg, errbuf, sizeof(errbuf));
    ck_assert_int_eq(ok, 0);
    ck_assert_msg(errbuf[0] != '\0', "Expected error for invalid source sigma2");
}
END_TEST

static Suite* config_runtime_suite(void) {
    Suite* s = suite_create("config_runtime");
    TCase* tc = tcase_create("core");
    tcase_add_test(tc, test_clamp_enforces_min_max_dims_and_cell_limit);
    tcase_add_test(tc, test_validate_rejects_invalid_cfl);
    tcase_add_test(tc, test_validate_rejects_invalid_sweep_range);
    tcase_add_test(tc, test_clamp_limits_material_and_source_counts);
    tcase_add_test(tc, test_clamp_normalizes_rect_coords_and_signs);
    tcase_add_test(tc, test_validate_rejects_invalid_material_rect);
    tcase_add_test(tc, test_validate_rejects_invalid_source_sigma2);
    suite_add_tcase(s, tc);
    return s;
}

int main(void) {
    Suite* s = config_runtime_suite();
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_ENV);
    int failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (failed == 0) ? 0 : 1;
}
