#include "analysis.h"

#include <assert.h>

static void test_scope_init_failure_leaves_scope_safe(void) {
    Scope scope = {0};
    analysis_test_set_alloc_fail_after(0);
    assert(scope_init(&scope, 128) == 0);
    assert(scope.y == NULL);
    assert(scope.n == 0);
    assert(scope.on == 0);
    analysis_test_set_alloc_fail_after(-1);
}

static void test_ports_init_failure_cleans_up_buffers(void) {
    Port ports[MAX_PORTS] = {0};
    analysis_test_set_alloc_fail_after(0);
    assert(ports_init(ports, 10, 10) == 0);
    for (int p = 0; p < MAX_PORTS; ++p) {
        assert(ports[p].V == NULL);
        assert(ports[p].I == NULL);
        assert(ports[p].n == 0);
        assert(ports[p].active == 0);
    }
    analysis_test_set_alloc_fail_after(-1);
}

int main(void) {
    test_scope_init_failure_leaves_scope_safe();
    test_ports_init_failure_cleans_up_buffers();
    return 0;
}
