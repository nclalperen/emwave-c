#ifndef CHECK_STUB_H
#define CHECK_STUB_H

#include <math.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*TFun)(int);

typedef struct Suite Suite;
typedef struct TCase TCase;
typedef struct SRunner SRunner;

Suite* suite_create(const char* name);
TCase* tcase_create(const char* name);
void tcase_add_test(TCase* tc, TFun test);
void suite_add_tcase(Suite* s, TCase* tc);
SRunner* srunner_create(Suite* s);
void srunner_run_all(SRunner* sr, int env);
int srunner_ntests_failed(SRunner* sr);
void srunner_free(SRunner* sr);

void ck_assert_failed(const char* file, int line, const char* expr, const char* fmt, ...);

#define CK_ENV 0

#define START_TEST(testname) static void testname(int ck_i)
#define END_TEST /* no-op */

#define ck_assert_msg(expr, ...) \
    do { \
        if (!(expr)) { \
            ck_assert_failed(__FILE__, __LINE__, #expr, __VA_ARGS__); \
            return; \
        } \
    } while (0)

#define ck_assert_int_eq(a, b) \
    ck_assert_msg((long long)(a) == (long long)(b), \
                  "Expected %s == %s (%lld vs %lld)", \
                  #a, #b, (long long)(a), (long long)(b))

#define ck_assert_ptr_null(p) \
    ck_assert_msg((p) == NULL, "Expected %s to be NULL", #p)

#define ck_assert_ptr_nonnull(p) \
    ck_assert_msg((p) != NULL, "Expected %s to be non-NULL", #p)

#define ck_assert_double_eq_tol(a, b, tol) \
    ck_assert_msg(fabs((double)(a) - (double)(b)) <= (double)(tol), \
                  "|%s - %s| (%.12g vs %.12g) exceeded tolerance %s", \
                  #a, #b, (double)(a), (double)(b), #tol)

#ifdef __cplusplus
}
#endif

#endif /* CHECK_STUB_H */
