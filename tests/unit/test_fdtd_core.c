#include "fdtd_core.h"
#include "config.h"

#include <check.h>
#include <math.h>

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

static Suite* fdtd_suite(void) {
    Suite* s = suite_create("fdtd_core");
    TCase* tc = tcase_create("core");
    tcase_add_test(tc, test_dt_uses_requested_safety);
    tcase_add_test(tc, test_dt_falls_back_to_default);
    tcase_add_test(tc, test_fdtd_init_handles_state_alloc_failure);
    tcase_add_test(tc, test_fdtd_init_handles_partial_alloc_failure);
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
