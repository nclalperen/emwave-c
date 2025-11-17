#include "fdtd_core.h"
#include "config.h"
#include "analysis.h"

#include <check.h>
#include <math.h>
#include <string.h>

START_TEST(test_dt_uses_requested_safety) {
    double dx = 0.01;
    double dy = 0.02;
    double dt = fdtd_compute_dt(dx, dy, 0.5);
    double denom = c0 * sqrt(1.0/(dx * dx) + 1.0/(dy * dy));
    double expected = 0.5 / denom;
    ck_assert_double_eq_tol(dt, expected, 1e-15);
}
END_TEST

START_TEST(test_dt_falls_back_to_default) {
    double dx = 0.01;
    double dy = 0.02;
    double dt_invalid = fdtd_compute_dt(dx, dy, -1.0);
    double dt_default = fdtd_compute_dt(dx, dy, CFL_SAFETY_FACTOR);
    ck_assert_double_eq_tol(dt_invalid, dt_default, 1e-15);
}
END_TEST

START_TEST(test_fdtd_init_handles_state_alloc_failure) {
    fdtd_test_set_alloc_fail_after(0);
    SimulationState* sim = fdtd_init(NULL);
    ck_assert_ptr_null(sim);
    fdtd_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_fdtd_init_handles_partial_alloc_failure) {
    fdtd_test_set_alloc_fail_after(3);
    SimulationState* sim = fdtd_init(NULL);
    ck_assert_ptr_null(sim);
    fdtd_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_fdtd_init_handles_port_alloc_failure) {
    analysis_test_set_alloc_fail_after(0);
    SimulationState* sim = fdtd_init(NULL);
    ck_assert_ptr_null(sim);
    analysis_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_fdtd_init_clamps_small_grid) {
    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
    cfg.nx = 1;
    cfg.ny = 1;
    SimulationState* sim = fdtd_init(&cfg);
    ck_assert_ptr_nonnull(sim);
    ck_assert_int_eq(sim->nx, SIM_MIN_DIM);
    ck_assert_int_eq(sim->ny, SIM_MIN_DIM);
    fdtd_free(sim);
}
END_TEST

START_TEST(test_fdtd_init_recovers_after_failure) {
    analysis_test_set_alloc_fail_after(0);
    ck_assert_ptr_null(fdtd_init(NULL));
    analysis_test_set_alloc_fail_after(-1);
    SimulationState* sim = fdtd_init(NULL);
    ck_assert_ptr_nonnull(sim);
    fdtd_free(sim);
}
END_TEST

START_TEST(test_fdtd_reset_skips_negative_port_lengths) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 4;
    state.ny = 4;

    double vbuf[2] = {1.0, 2.0};
    double ibuf[2] = {3.0, 4.0};
    state.ports[0].V = vbuf;
    state.ports[0].I = ibuf;
    state.ports[0].n = -10;
    state.ports[0].head = 7;

    fdtd_reset(&state);

    ck_assert_int_eq(state.ports[0].head, 0);
    ck_assert_double_eq_tol(vbuf[0], 1.0, 1e-12);
    ck_assert_double_eq_tol(ibuf[1], 4.0, 1e-12);
}
END_TEST

START_TEST(test_fdtd_reset_normalizes_port_geometry) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 4;
    state.ny = 4;

    double vbuf[4] = {5.0, 6.0, 7.0, 8.0};
    double ibuf[4] = {1.0, 2.0, 3.0, 4.0};
    state.ports[0].V = vbuf;
    state.ports[0].I = ibuf;
    state.ports[0].n = 4;
    state.ports[0].y0 = 5;
    state.ports[0].y1 = 2;
    state.ports[0].len = 99;
    state.ports[0].head = 3;
    state.ports[0].active = 1;

    fdtd_reset(&state);

    ck_assert_int_eq(state.ports[0].len, 0);
    ck_assert_int_eq(state.ports[0].active, 0);
    ck_assert_int_eq(state.ports[0].head, 0);
    for (int i = 0; i < 4; ++i) {
        ck_assert_double_eq_tol(vbuf[i], 0.0, 1e-15);
        ck_assert_double_eq_tol(ibuf[i], 0.0, 1e-15);
    }
}
END_TEST

START_TEST(test_fdtd_reset_reactivates_valid_ports) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 6;
    state.ny = 6;

    double vbuf[4] = {0};
    double ibuf[4] = {0};
    state.ports[0].V = vbuf;
    state.ports[0].I = ibuf;
    state.ports[0].n = 4;
    state.ports[0].y0 = 1;
    state.ports[0].y1 = 3;
    state.ports[0].active = 0;
    state.ports[0].head = 2;

    fdtd_reset(&state);

    ck_assert_int_eq(state.ports[0].active, 1);
    ck_assert_int_eq(state.ports[0].head, 0);
    ck_assert_int_eq(state.ports[0].len, 3);
    for (int i = 0; i < 4; ++i) {
        ck_assert_double_eq_tol(vbuf[i], 0.0, 1e-15);
        ck_assert_double_eq_tol(ibuf[i], 0.0, 1e-15);
    }
}
END_TEST

static Suite* fdtd_suite(void) {
    Suite* s = suite_create("fdtd_core");
    TCase* tc = tcase_create("core");
    tcase_add_test(tc, test_dt_uses_requested_safety);
    tcase_add_test(tc, test_dt_falls_back_to_default);
    tcase_add_test(tc, test_fdtd_init_handles_state_alloc_failure);
    tcase_add_test(tc, test_fdtd_init_handles_partial_alloc_failure);
    tcase_add_test(tc, test_fdtd_init_handles_port_alloc_failure);
    tcase_add_test(tc, test_fdtd_init_clamps_small_grid);
    tcase_add_test(tc, test_fdtd_init_recovers_after_failure);
    tcase_add_test(tc, test_fdtd_reset_skips_negative_port_lengths);
    tcase_add_test(tc, test_fdtd_reset_normalizes_port_geometry);
    tcase_add_test(tc, test_fdtd_reset_reactivates_valid_ports);
    suite_add_tcase(s, tc);
    return s;
}

int main(void) {
    Suite* s = fdtd_suite();
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_ENV);
    int failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (failed == 0) ? 0 : 1;
}
