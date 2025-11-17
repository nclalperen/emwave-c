#include "analysis.h"

#include <check.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double** rows;
    double* data;
} DummyField;

static void dummy_field_init(DummyField* field, int nx, int ny) {
    size_t cells = (size_t)nx * (size_t)ny;
    field->data = (double*)calloc(cells, sizeof(double));
    field->rows = (double**)malloc(sizeof(double*) * (size_t)nx);
    ck_assert_ptr_nonnull(field->data);
    ck_assert_ptr_nonnull(field->rows);
    for (int i = 0; i < nx; ++i) {
        field->rows[i] = field->data + (size_t)i * ny;
    }
}

static void dummy_field_free(DummyField* field) {
    free(field->rows);
    free(field->data);
    field->rows = NULL;
    field->data = NULL;
}

static void alloc_port_buffers(Port* port, int length) {
    port->n = length;
    port->V = (double*)calloc((size_t)length, sizeof(double));
    port->I = (double*)calloc((size_t)length, sizeof(double));
    ck_assert_ptr_nonnull(port->V);
    ck_assert_ptr_nonnull(port->I);
}

static void free_port_buffers(Port* port) {
    free(port->V);
    free(port->I);
    port->V = NULL;
    port->I = NULL;
    port->n = 0;
}

START_TEST(test_scope_init_failure_leaves_scope_safe) {
    Scope scope = {0};
    analysis_test_set_alloc_fail_after(0);
    ck_assert_int_eq(scope_init(&scope, 128), 0);
    ck_assert_ptr_null(scope.y);
    ck_assert_int_eq(scope.n, 0);
    ck_assert_int_eq(scope.on, 0);
    analysis_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_scope_init_enforces_minimum_length) {
    Scope scope = {0};
    ck_assert_int_eq(scope_init(&scope, -5), 1);
    ck_assert_int_eq(scope.n, 64);
    scope_free(&scope);
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
    analysis_test_set_alloc_fail_after(-1);
}
END_TEST

START_TEST(test_ports_init_rejects_too_small_grid) {
    Port ports[MAX_PORTS] = {0};
    ck_assert_int_eq(ports_init(ports, 1, 1), 0);
    for (int p = 0; p < MAX_PORTS; ++p) {
        ck_assert_ptr_null(ports[p].V);
        ck_assert_int_eq(ports[p].active, 0);
    }
}
END_TEST

START_TEST(test_ports_init_enables_ports_on_success) {
    Port ports[MAX_PORTS] = {0};
    ck_assert_int_eq(ports_init(ports, 16, 16), 1);
    for (int p = 0; p < MAX_PORTS; ++p) {
        ck_assert_ptr_nonnull(ports[p].V);
        ck_assert_ptr_nonnull(ports[p].I);
        ck_assert_int_eq(ports[p].active, 1);
        ck_assert_int_eq(ports[p].head, 0);
    }
    ports_free(ports);
}
END_TEST

START_TEST(test_ports_sample_skips_unbacked_ports) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 4;
    state.ny = 4;

    DummyField ez;
    DummyField hy;
    dummy_field_init(&ez, state.nx, state.ny);
    dummy_field_init(&hy, state.nx, state.ny);
    state.Ez = ez.rows;
    state.Ez_data = ez.data;
    state.Hy = hy.rows;
    state.Hy_data = hy.data;

    state.ports[0].active = 1;
    state.ports[0].x = 2;
    state.ports[0].y0 = 1;
    state.ports[0].y1 = 2;
    state.ports[0].n = PORT_SIGNAL_LENGTH;
    state.ports[0].head = 5;

    ports_sample(&state, 1.0, 1.0);

    ck_assert_int_eq(state.ports[0].head, 5);

    dummy_field_free(&ez);
    dummy_field_free(&hy);
}
END_TEST

START_TEST(test_ports_sample_clamps_port_x_inside_grid) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 8;
    state.ny = 6;

    DummyField ez;
    DummyField hy;
    dummy_field_init(&ez, state.nx, state.ny);
    dummy_field_init(&hy, state.nx, state.ny);
    state.Ez = ez.rows;
    state.Ez_data = ez.data;
    state.Hy = hy.rows;
    state.Hy_data = hy.data;

    for (int i = 0; i < state.nx; ++i) {
        for (int j = 0; j < state.ny; ++j) {
            state.Ez[i][j] = 1.0;
            state.Hy[i][j] = 1.0;
        }
    }

    Port* port = &state.ports[0];
    port->active = 1;
    port->x = state.nx + 5;
    port->y0 = 0;
    port->y1 = state.ny - 1;
    alloc_port_buffers(port, 4);

    ports_sample(&state, 1.0, 1.0);

    ck_assert_int_eq(port->x, state.nx - 2);
    ck_assert_int_eq(port->head, 1);
    ck_assert_int_eq(port->len, state.ny);

    free_port_buffers(port);
    dummy_field_free(&ez);
    dummy_field_free(&hy);
}
END_TEST

START_TEST(test_ports_sample_clamps_segment_and_updates_len) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 6;
    state.ny = 4;

    DummyField ez;
    DummyField hy;
    dummy_field_init(&ez, state.nx, state.ny);
    dummy_field_init(&hy, state.nx, state.ny);
    state.Ez = ez.rows;
    state.Ez_data = ez.data;
    state.Hy = hy.rows;
    state.Hy_data = hy.data;

    for (int i = 0; i < state.nx; ++i) {
        for (int j = 0; j < state.ny; ++j) {
            state.Ez[i][j] = 1.0;
            state.Hy[i][j] = 2.0;
        }
    }

    Port* port = &state.ports[0];
    port->active = 1;
    port->x = 3;
    port->y0 = -10;
    port->y1 = 99;
    alloc_port_buffers(port, 8);

    ports_sample(&state, 0.5, 0.25);

    ck_assert_int_eq(port->len, state.ny);
    ck_assert_int_eq(port->head, 1);
    ck_assert_double_eq_tol(port->V[0], state.ny * 0.25, 1e-12);
    ck_assert_double_eq_tol(port->I[0], state.ny * 2.0 * 0.5, 1e-12);

    free_port_buffers(port);
    dummy_field_free(&ez);
    dummy_field_free(&hy);
}
END_TEST

START_TEST(test_ports_sample_requires_minimum_span) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 4;
    state.ny = 1; /* intentionally smaller than kPortMinSegmentLen */

    DummyField ez;
    DummyField hy;
    dummy_field_init(&ez, state.nx, state.ny);
    dummy_field_init(&hy, state.nx, state.ny);
    state.Ez = ez.rows;
    state.Ez_data = ez.data;
    state.Hy = hy.rows;
    state.Hy_data = hy.data;

    Port* port = &state.ports[0];
    port->active = 1;
    port->x = 2;
    port->y0 = 0;
    port->y1 = 0;
    alloc_port_buffers(port, 4);
    port->head = 2;
    port->V[2] = 123.0;
    port->I[2] = 456.0;

    ports_sample(&state, 1.0, 1.0);

    ck_assert_int_eq(port->active, 0);
    ck_assert_int_eq(port->head, 0);
    for (int i = 0; i < port->n; ++i) {
        ck_assert_double_eq_tol(port->V[i], 0.0, 1e-12);
        ck_assert_double_eq_tol(port->I[i], 0.0, 1e-12);
    }

    free_port_buffers(port);
    dummy_field_free(&ez);
    dummy_field_free(&hy);
}
END_TEST

START_TEST(test_ports_sample_deactivates_when_grid_has_no_interior) {
    SimulationState state;
    memset(&state, 0, sizeof(state));
    state.nx = 2;
    state.ny = 4;

    DummyField ez;
    DummyField hy;
    dummy_field_init(&ez, state.nx, state.ny);
    dummy_field_init(&hy, state.nx, state.ny);
    state.Ez = ez.rows;
    state.Ez_data = ez.data;
    state.Hy = hy.rows;
    state.Hy_data = hy.data;

    Port* port = &state.ports[0];
    port->active = 1;
    port->x = 0;
    port->y0 = 0;
    port->y1 = 3;
    alloc_port_buffers(port, 4);
    port->head = 3;

    ports_sample(&state, 1.0, 1.0);

    ck_assert_int_eq(port->active, 0);
    ck_assert_int_eq(port->head, 0);
    ck_assert_int_eq(port->len, 0);

    free_port_buffers(port);
    dummy_field_free(&ez);
    dummy_field_free(&hy);
}
END_TEST

static Suite* analysis_suite(void) {
    Suite* s = suite_create("analysis");
    TCase* tc = tcase_create("alloc");
    tcase_add_test(tc, test_scope_init_failure_leaves_scope_safe);
    tcase_add_test(tc, test_scope_init_enforces_minimum_length);
    tcase_add_test(tc, test_ports_init_failure_cleans_up_buffers);
    tcase_add_test(tc, test_ports_init_rejects_too_small_grid);
    tcase_add_test(tc, test_ports_init_enables_ports_on_success);
    tcase_add_test(tc, test_ports_sample_skips_unbacked_ports);
    tcase_add_test(tc, test_ports_sample_clamps_port_x_inside_grid);
    tcase_add_test(tc, test_ports_sample_clamps_segment_and_updates_len);
    tcase_add_test(tc, test_ports_sample_requires_minimum_span);
    tcase_add_test(tc, test_ports_sample_deactivates_when_grid_has_no_interior);
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
