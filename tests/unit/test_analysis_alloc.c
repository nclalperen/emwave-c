#include "analysis.h"

#include <check.h>

START_TEST(test_scope_init_failure_leaves_scope_safe) {
    Scope scope = {0};
    analysis_test_set_alloc_fail_after(0);
    ck_assert_int_eq(scope_init(&scope, 128), 0);
    ck_assert_ptr_null(scope.y);
    ck_assert_int_eq(scope.n, 0);
    ck_assert_int_eq(scope.on, 0);
    scope_free(&scope);
    analysis_test_set_alloc_fail_after(-1);
    ck_assert_int_eq(scope_init(&scope, 128), 1);
    scope_free(&scope);
    analysis_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_ports_init_failure_cleans_up_buffers) {
    Port ports[MAX_PORTS] = {0};
    analysis_test_set_alloc_fail_after(0);
    ck_assert_int_eq(ports_init(ports, 10, 10), 0);
    for (int p = 0; p < MAX_PORTS; ++p) {
        ck_assert_ptr_null(ports[p].V);
        ck_assert_ptr_null(ports[p].I);
        ck_assert_int_eq(ports[p].n, 0);
        ck_assert_int_eq(ports[p].active, 0);
    }
    ports_free(ports);
    analysis_test_set_alloc_fail_after(-1);
    ck_assert_int_eq(ports_init(ports, 10, 10), 1);
    ports_free(ports);
    analysis_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_ports_init_multiple_failures) {
    Port ports[MAX_PORTS] = {0};
    for (int fail_after = 0; fail_after < 4; ++fail_after) {
        analysis_test_set_alloc_fail_after(fail_after);
        ck_assert_int_eq(ports_init(ports, 12, 12), 0);
        ports_free(ports);
    }
    analysis_test_set_alloc_fail_after(-1);
    ck_assert_int_eq(ports_init(ports, 12, 12), 1);
    ports_free(ports);
}
END_TEST

static Suite* analysis_suite(void) {
    Suite* s = suite_create("analysis");
    TCase* tc = tcase_create("alloc");
    tcase_add_test(tc, test_scope_init_failure_leaves_scope_safe);
    tcase_add_test(tc, test_ports_init_failure_cleans_up_buffers);
    tcase_add_test(tc, test_ports_init_multiple_failures);
    suite_add_tcase(s, tc);
    return s;
}

int main(void) {
    Suite* s = analysis_suite();
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_ENV);
    int failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (failed == 0) ? 0 : 1;
}
