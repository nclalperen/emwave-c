// ============================================================================
// emwave-c: Lightweight expression engine for time-domain sources
// Shunting-yard parser + RPN evaluator with a tiny feature set.
// ============================================================================

#include "expr.h"

#include <ctype.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum {
    EXPR_TOK_NUMBER,
    EXPR_TOK_VAR,       /* t, amp, freq, pi */
    EXPR_TOK_FUNC,      /* sin, cos, exp, sqrt, abs */
    EXPR_TOK_OP,        /* + - * / ^ */
    EXPR_TOK_LPAREN,
    EXPR_TOK_RPAREN,
    EXPR_TOK_COMMA
} ExprTokenType;

typedef enum {
    EXPR_VAR_T,
    EXPR_VAR_AMP,
    EXPR_VAR_FREQ,
    EXPR_VAR_PI
} ExprVar;

typedef enum {
    EXPR_FUNC_SIN,
    EXPR_FUNC_COS,
    EXPR_FUNC_EXP,
    EXPR_FUNC_SQRT,
    EXPR_FUNC_ABS
} ExprFunc;

typedef struct {
    ExprTokenType type;
    union {
        double number;
        ExprVar var;
        ExprFunc func;
        char op;
    } v;
} ExprToken;

typedef enum {
    EXPR_OP_PUSH_CONST,
    EXPR_OP_PUSH_VAR,
    EXPR_OP_FN1,
    EXPR_OP_ADD,
    EXPR_OP_SUB,
    EXPR_OP_MUL,
    EXPR_OP_DIV,
    EXPR_OP_POW
} ExprOpCode;

typedef struct {
    ExprOpCode op;
    union {
        double number;
        ExprVar var;
        ExprFunc func;
    } a;
} ExprOp;

struct ExprProgram {
    ExprOp* ops;
    int count;
};

static int op_precedence(char op) {
    switch (op) {
        case '+':
        case '-':
            return 1;
        case '*':
        case '/':
            return 2;
        case '^':
            return 3;
        default:
            return 0;
    }
}

static int op_is_right_associative(char op) {
    return (op == '^');
}

static int parse_ident(const char* s, int len, ExprToken* out) {
    if (len <= 0 || !out) return 0;
    if (len == 1 && s[0] == 't') {
        out->type = EXPR_TOK_VAR;
        out->v.var = EXPR_VAR_T;
        return 1;
    }
    if (len == 3 && strncmp(s, "amp", 3) == 0) {
        out->type = EXPR_TOK_VAR;
        out->v.var = EXPR_VAR_AMP;
        return 1;
    }
    if (len == 4 && strncmp(s, "freq", 4) == 0) {
        out->type = EXPR_TOK_VAR;
        out->v.var = EXPR_VAR_FREQ;
        return 1;
    }
    if (len == 2 && strncmp(s, "pi", 2) == 0) {
        out->type = EXPR_TOK_VAR;
        out->v.var = EXPR_VAR_PI;
        return 1;
    }
    if (len == 3 && strncmp(s, "sin", 3) == 0) {
        out->type = EXPR_TOK_FUNC;
        out->v.func = EXPR_FUNC_SIN;
        return 1;
    }
    if (len == 3 && strncmp(s, "cos", 3) == 0) {
        out->type = EXPR_TOK_FUNC;
        out->v.func = EXPR_FUNC_COS;
        return 1;
    }
    if (len == 3 && strncmp(s, "exp", 3) == 0) {
        out->type = EXPR_TOK_FUNC;
        out->v.func = EXPR_FUNC_EXP;
        return 1;
    }
    if (len == 4 && strncmp(s, "sqrt", 4) == 0) {
        out->type = EXPR_TOK_FUNC;
        out->v.func = EXPR_FUNC_SQRT;
        return 1;
    }
    if (len == 3 && strncmp(s, "abs", 3) == 0) {
        out->type = EXPR_TOK_FUNC;
        out->v.func = EXPR_FUNC_ABS;
        return 1;
    }
    return 0;
}

static int tokenize(const char* expr, ExprToken** out_tokens, int* out_count, char* errbuf, int errbuf_len) {
    if (errbuf && errbuf_len > 0) errbuf[0] = '\0';
    if (!expr || !out_tokens || !out_count) return 0;

    int cap = 32;
    int count = 0;
    ExprToken* tokens = (ExprToken*)malloc(sizeof(ExprToken) * cap);
    if (!tokens) return 0;

    const char* p = expr;
    while (*p) {
        if (isspace((unsigned char)*p)) {
            p++;
            continue;
        }
        if (count >= cap) {
            cap *= 2;
            ExprToken* tmp = (ExprToken*)realloc(tokens, sizeof(ExprToken) * cap);
            if (!tmp) {
                free(tokens);
                return 0;
            }
            tokens = tmp;
        }
        if (isdigit((unsigned char)*p) || *p == '.') {
            char* endptr = NULL;
            double v = strtod(p, &endptr);
            if (endptr == p) {
                if (errbuf && errbuf_len > 0) {
                    snprintf(errbuf, errbuf_len, "Invalid number near '%.8s'", p);
                }
                free(tokens);
                return 0;
            }
            tokens[count].type = EXPR_TOK_NUMBER;
            tokens[count].v.number = v;
            count++;
            p = endptr;
            continue;
        }
        if (isalpha((unsigned char)*p)) {
            const char* start = p;
            while (isalnum((unsigned char)*p) || *p == '_') p++;
            int len = (int)(p - start);
            ExprToken tok;
            if (!parse_ident(start, len, &tok)) {
                if (errbuf && errbuf_len > 0) {
                    snprintf(errbuf, errbuf_len, "Unknown identifier '%.*s'", len, start);
                }
                free(tokens);
                return 0;
            }
            tokens[count++] = tok;
            continue;
        }
        switch (*p) {
            case '+': case '-': case '*': case '/': case '^':
                tokens[count].type = EXPR_TOK_OP;
                tokens[count].v.op = *p;
                count++;
                p++;
                break;
            case '(':
                tokens[count].type = EXPR_TOK_LPAREN;
                count++;
                p++;
                break;
            case ')':
                tokens[count].type = EXPR_TOK_RPAREN;
                count++;
                p++;
                break;
            case ',':
                tokens[count].type = EXPR_TOK_COMMA;
                count++;
                p++;
                break;
            default:
                if (errbuf && errbuf_len > 0) {
                    snprintf(errbuf, errbuf_len, "Unexpected character '%c'", *p);
                }
                free(tokens);
                return 0;
        }
    }

    *out_tokens = tokens;
    *out_count = count;
    return 1;
}

static int emit_op(ExprOp** ops, int* count, int* cap, ExprOp op) {
    if (*count >= *cap) {
        int new_cap = (*cap == 0) ? 32 : (*cap * 2);
        ExprOp* tmp = (ExprOp*)realloc(*ops, sizeof(ExprOp) * new_cap);
        if (!tmp) return 0;
        *ops = tmp;
        *cap = new_cap;
    }
    (*ops)[*count] = op;
    (*count)++;
    return 1;
}

static int shunting_yard(const ExprToken* tokens, int token_count,
                         ExprOp** out_ops, int* out_count,
                         char* errbuf, int errbuf_len) {
    if (errbuf && errbuf_len > 0) errbuf[0] = '\0';
    ExprToken* op_stack = (ExprToken*)malloc(sizeof(ExprToken) * token_count);
    int op_top = 0;
    ExprOp* ops = NULL;
    int op_count = 0;
    int op_cap = 0;

    int expect_operand = 1;

    for (int i = 0; i < token_count; ++i) {
        ExprToken tok = tokens[i];
        switch (tok.type) {
            case EXPR_TOK_NUMBER: {
                ExprOp op;
                op.op = EXPR_OP_PUSH_CONST;
                op.a.number = tok.v.number;
                if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                expect_operand = 0;
                break;
            }
            case EXPR_TOK_VAR: {
                ExprOp op;
                op.op = EXPR_OP_PUSH_VAR;
                op.a.var = tok.v.var;
                if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                expect_operand = 0;
                break;
            }
            case EXPR_TOK_FUNC:
                op_stack[op_top++] = tok;
                expect_operand = 1;
                break;
            case EXPR_TOK_OP: {
                char op_char = tok.v.op;
                if (expect_operand && op_char == '-') {
                    ExprOp op;
                    op.op = EXPR_OP_PUSH_CONST;
                    op.a.number = 0.0;
                    if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                } else if (expect_operand) {
                    if (errbuf && errbuf_len > 0) {
                        snprintf(errbuf, errbuf_len, "Operator '%c' in unexpected position", op_char);
                    }
                    goto fail;
                }
                while (op_top > 0) {
                    ExprToken top = op_stack[op_top - 1];
                    if (top.type != EXPR_TOK_OP) break;
                    int prec_top = op_precedence(top.v.op);
                    int prec_cur = op_precedence(op_char);
                    if ((prec_top > prec_cur) ||
                        (prec_top == prec_cur && !op_is_right_associative(op_char))) {
                        op_top--;
                        ExprOp op;
                        switch (top.v.op) {
                            case '+': op.op = EXPR_OP_ADD; break;
                            case '-': op.op = EXPR_OP_SUB; break;
                            case '*': op.op = EXPR_OP_MUL; break;
                            case '/': op.op = EXPR_OP_DIV; break;
                            case '^': op.op = EXPR_OP_POW; break;
                            default: goto fail;
                        }
                        if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                    } else {
                        break;
                    }
                }
                op_stack[op_top++] = tok;
                expect_operand = 1;
                break;
            }
            case EXPR_TOK_LPAREN:
                op_stack[op_top++] = tok;
                expect_operand = 1;
                break;
            case EXPR_TOK_RPAREN: {
                int found_lparen = 0;
                while (op_top > 0) {
                    ExprToken top = op_stack[--op_top];
                    if (top.type == EXPR_TOK_LPAREN) {
                        found_lparen = 1;
                        break;
                    }
                    if (top.type == EXPR_TOK_FUNC) {
                        ExprOp op;
                        op.op = EXPR_OP_FN1;
                        op.a.func = top.v.func;
                        if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                    } else if (top.type == EXPR_TOK_OP) {
                        ExprOp op;
                        switch (top.v.op) {
                            case '+': op.op = EXPR_OP_ADD; break;
                            case '-': op.op = EXPR_OP_SUB; break;
                            case '*': op.op = EXPR_OP_MUL; break;
                            case '/': op.op = EXPR_OP_DIV; break;
                            case '^': op.op = EXPR_OP_POW; break;
                            default: goto fail;
                        }
                        if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                    }
                }
                if (!found_lparen) {
                    if (errbuf && errbuf_len > 0) {
                        snprintf(errbuf, errbuf_len, "Mismatched parentheses");
                    }
                    goto fail;
                }
                expect_operand = 0;
                break;
            }
            case EXPR_TOK_COMMA:
                while (op_top > 0 && op_stack[op_top - 1].type != EXPR_TOK_LPAREN) {
                    ExprToken top = op_stack[--op_top];
                    if (top.type == EXPR_TOK_FUNC) {
                        ExprOp op;
                        op.op = EXPR_OP_FN1;
                        op.a.func = top.v.func;
                        if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                    } else if (top.type == EXPR_TOK_OP) {
                        ExprOp op;
                        switch (top.v.op) {
                            case '+': op.op = EXPR_OP_ADD; break;
                            case '-': op.op = EXPR_OP_SUB; break;
                            case '*': op.op = EXPR_OP_MUL; break;
                            case '/': op.op = EXPR_OP_DIV; break;
                            case '^': op.op = EXPR_OP_POW; break;
                            default: goto fail;
                        }
                        if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
                    }
                }
                expect_operand = 1;
                break;
        }
    }

    while (op_top > 0) {
        ExprToken top = op_stack[--op_top];
        if (top.type == EXPR_TOK_LPAREN || top.type == EXPR_TOK_RPAREN) {
            if (errbuf && errbuf_len > 0) {
                snprintf(errbuf, errbuf_len, "Mismatched parentheses at end");
            }
            goto fail;
        }
        if (top.type == EXPR_TOK_FUNC) {
            ExprOp op;
            op.op = EXPR_OP_FN1;
            op.a.func = top.v.func;
            if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
        } else if (top.type == EXPR_TOK_OP) {
            ExprOp op;
            switch (top.v.op) {
                case '+': op.op = EXPR_OP_ADD; break;
                case '-': op.op = EXPR_OP_SUB; break;
                case '*': op.op = EXPR_OP_MUL; break;
                case '/': op.op = EXPR_OP_DIV; break;
                case '^': op.op = EXPR_OP_POW; break;
                default: goto fail;
            }
            if (!emit_op(&ops, &op_count, &op_cap, op)) goto oom;
        }
    }

    free(op_stack);
    *out_ops = ops;
    *out_count = op_count;
    return 1;

oom:
    if (errbuf && errbuf_len > 0) {
        snprintf(errbuf, errbuf_len, "Out of memory while compiling expression");
    }
fail:
    free(op_stack);
    free(ops);
    return 0;
}

int expr_compile(const char* expr_str, ExprProgram** out_prog, char* errbuf, int errbuf_len) {
    if (errbuf && errbuf_len > 0) errbuf[0] = '\0';
    if (!expr_str || !out_prog) return 0;

    ExprToken* tokens = NULL;
    int token_count = 0;
    if (!tokenize(expr_str, &tokens, &token_count, errbuf, errbuf_len)) {
        return 0;
    }

    ExprOp* ops = NULL;
    int op_count = 0;
    if (!shunting_yard(tokens, token_count, &ops, &op_count, errbuf, errbuf_len)) {
        free(tokens);
        return 0;
    }
    free(tokens);

    ExprProgram* prog = (ExprProgram*)malloc(sizeof(ExprProgram));
    if (!prog) {
        free(ops);
        if (errbuf && errbuf_len > 0) {
            snprintf(errbuf, errbuf_len, "Out of memory allocating program");
        }
        return 0;
    }
    prog->ops = ops;
    prog->count = op_count;
    *out_prog = prog;
    return 1;
}

double expr_eval(const ExprProgram* prog, double t, double amp, double freq) {
    if (!prog || prog->count <= 0) return 0.0;
    double stack_local[64];
    double* stack = stack_local;
    int sp = 0;

    for (int i = 0; i < prog->count; ++i) {
        ExprOp op = prog->ops[i];
        switch (op.op) {
            case EXPR_OP_PUSH_CONST:
                if (sp < (int)(sizeof(stack_local) / sizeof(stack_local[0]))) {
                    stack[sp++] = op.a.number;
                }
                break;
            case EXPR_OP_PUSH_VAR:
                if (sp < (int)(sizeof(stack_local) / sizeof(stack_local[0]))) {
                    switch (op.a.var) {
                        case EXPR_VAR_T: stack[sp++] = t; break;
                        case EXPR_VAR_AMP: stack[sp++] = amp; break;
                        case EXPR_VAR_FREQ: stack[sp++] = freq; break;
                        case EXPR_VAR_PI: stack[sp++] = M_PI; break;
                    }
                }
                break;
            case EXPR_OP_FN1:
                if (sp <= 0) break;
                stack[sp - 1] = (op.a.func == EXPR_FUNC_SIN) ? sin(stack[sp - 1]) :
                                (op.a.func == EXPR_FUNC_COS) ? cos(stack[sp - 1]) :
                                (op.a.func == EXPR_FUNC_EXP) ? exp(stack[sp - 1]) :
                                (op.a.func == EXPR_FUNC_SQRT) ? sqrt(fabs(stack[sp - 1])) :
                                fabs(stack[sp - 1]);
                break;
            case EXPR_OP_ADD:
                if (sp >= 2) {
                    stack[sp - 2] = stack[sp - 2] + stack[sp - 1];
                    sp--;
                }
                break;
            case EXPR_OP_SUB:
                if (sp >= 2) {
                    stack[sp - 2] = stack[sp - 2] - stack[sp - 1];
                    sp--;
                }
                break;
            case EXPR_OP_MUL:
                if (sp >= 2) {
                    stack[sp - 2] = stack[sp - 2] * stack[sp - 1];
                    sp--;
                }
                break;
            case EXPR_OP_DIV:
                if (sp >= 2) {
                    double denom = stack[sp - 1];
                    if (denom == 0.0) denom = 1e-18;
                    stack[sp - 2] = stack[sp - 2] / denom;
                    sp--;
                }
                break;
            case EXPR_OP_POW:
                if (sp >= 2) {
                    stack[sp - 2] = pow(stack[sp - 2], stack[sp - 1]);
                    sp--;
                }
                break;
        }
    }

    if (sp <= 0) {
        return 0.0;
    }
    return stack[sp - 1];
}

void expr_free(ExprProgram* prog) {
    if (!prog) return;
    free(prog->ops);
    prog->ops = NULL;
    prog->count = 0;
    free(prog);
}
