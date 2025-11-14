#include "fdtd_core.h"
#include "config.h"

#include <assert.h>
#include <math.h>

static void test_dt_uses_requested_safety(void) {
    double dx = 0.01;
    double dy = 0.02;
    double dt = fdtd_compute_dt(dx, dy, 0.5);
    double denom = c0 * sqrt(1.0/(dx * dx) + 1.0/(dy * dy));
    double expected = 0.5 / denom;
    assert(fabs(dt - expected) < 1e-15);
}

static void test_dt_falls_back_to_default(void) {
    double dx = 0.01;
    double dy = 0.02;
    double dt_invalid = fdtd_compute_dt(dx, dy, -1.0);
    double dt_default = fdtd_compute_dt(dx, dy, CFL_SAFETY_FACTOR);
    assert(fabs(dt_invalid - dt_default) < 1e-15);
}

static void test_fdtd_init_handles_state_alloc_failure(void) {
    fdtd_test_set_alloc_fail_after(0);
    SimulationState* sim = fdtd_init(NULL);
    assert(sim == NULL);
    fdtd_test_set_alloc_fail_after(-1);
}

static void test_fdtd_init_handles_partial_alloc_failure(void) {
    fdtd_test_set_alloc_fail_after(3);
    SimulationState* sim = fdtd_init(NULL);
    assert(sim == NULL);
    fdtd_test_set_alloc_fail_after(-1);
}

int main(void) {
    test_dt_uses_requested_safety();
    test_dt_falls_back_to_default();
    test_fdtd_init_handles_state_alloc_failure();
    test_fdtd_init_handles_partial_alloc_failure();
    return 0;
}
