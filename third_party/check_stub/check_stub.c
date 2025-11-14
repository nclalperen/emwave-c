#include "check.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct CheckTest {
    TFun func;
    struct CheckTest* next;
};

struct TCase {
    const char* name;
    struct CheckTest* head;
    struct CheckTest* tail;
    struct TCase* next;
};

struct Suite {
    const char* name;
    struct TCase* head;
    struct TCase* tail;
};

struct SRunner {
    Suite* suite;
    int failed;
    int ran;
};

static SRunner* current_runner = NULL;

static void* ck_calloc(size_t count, size_t size) {
    void* ptr = calloc(count, size);
    if (!ptr) {
        fprintf(stderr, "[check_stub] Out of memory\n");
        abort();
    }
    return ptr;
}

Suite* suite_create(const char* name) {
    Suite* s = (Suite*)ck_calloc(1, sizeof(Suite));
    s->name = name ? name : "suite";
    return s;
}

TCase* tcase_create(const char* name) {
    TCase* tc = (TCase*)ck_calloc(1, sizeof(TCase));
    tc->name = name ? name : "tcase";
    return tc;
}

void tcase_add_test(TCase* tc, TFun test) {
    if (!tc || !test) return;
    struct CheckTest* node = (struct CheckTest*)ck_calloc(1, sizeof(struct CheckTest));
    node->func = test;
    if (!tc->head) {
        tc->head = node;
        tc->tail = node;
    } else {
        tc->tail->next = node;
        tc->tail = node;
    }
}

void suite_add_tcase(Suite* s, TCase* tc) {
    if (!s || !tc) return;
    if (!s->head) {
        s->head = tc;
        s->tail = tc;
    } else {
        s->tail->next = tc;
        s->tail = tc;
    }
}

SRunner* srunner_create(Suite* s) {
    SRunner* sr = (SRunner*)ck_calloc(1, sizeof(SRunner));
    sr->suite = s;
    return sr;
}

void srunner_run_all(SRunner* sr, int env) {
    (void)env;
    if (!sr || !sr->suite) return;
    current_runner = sr;
    sr->failed = 0;
    sr->ran = 0;
    for (TCase* tc = sr->suite->head; tc; tc = tc->next) {
        for (struct CheckTest* test = tc->head; test; test = test->next) {
            sr->ran++;
            test->func(0);
        }
    }
    current_runner = NULL;
}

int srunner_ntests_failed(SRunner* sr) {
    return sr ? sr->failed : 0;
}

void srunner_free(SRunner* sr) {
    if (!sr) return;
    if (sr->suite) {
        TCase* tc = sr->suite->head;
        while (tc) {
            struct CheckTest* t = tc->head;
            while (t) {
                struct CheckTest* next_t = t->next;
                free(t);
                t = next_t;
            }
            TCase* next_tc = tc->next;
            free(tc);
            tc = next_tc;
        }
        free(sr->suite);
    }
    free(sr);
}

void ck_assert_failed(const char* file, int line, const char* expr, const char* fmt, ...) {
    if (current_runner) {
        current_runner->failed++;
    }
    fprintf(stderr, "[check] Assertion failed at %s:%d: %s\n", file, line, expr);
    if (fmt) {
        va_list ap;
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        fputc('\n', stderr);
    }
}
