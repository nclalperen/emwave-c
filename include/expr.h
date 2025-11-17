// ============================================================================
// emwave-c: Lightweight expression engine for time-domain sources
// ============================================================================

#ifndef EMWAVE_EXPR_H
#define EMWAVE_EXPR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ExprProgram ExprProgram;

/* Compile an expression string into an executable program.
 * Supported:
 *   - Variables: t, amp, freq, pi
 *   - Operators: +, -, *, /, ^, unary -
 *   - Functions: sin, cos, exp, sqrt, abs
 * Returns 1 on success, 0 on failure.
 */
int expr_compile(const char* expr_str, ExprProgram** out_prog, char* errbuf, int errbuf_len);

/* Evaluate a compiled program for the given variables. */
double expr_eval(const ExprProgram* prog, double t, double amp, double freq);

/* Free a compiled program. */
void expr_free(ExprProgram* prog);

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_EXPR_H */

