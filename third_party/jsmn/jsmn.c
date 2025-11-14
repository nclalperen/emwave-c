#include "jsmn.h"

#include <stdlib.h>

static jsmntok_t* jsmn_alloc_token(jsmn_parser* parser, jsmntok_t* tokens, unsigned int num_tokens) {
    if (!parser || !tokens || parser->toknext >= num_tokens) {
        return NULL;
    }
    jsmntok_t* tok = &tokens[parser->toknext++];
    tok->start = tok->end = -1;
    tok->size = 0;
    tok->parent = -1;
    tok->type = JSMN_UNDEFINED;
    return tok;
}

static void jsmn_fill_token(jsmntok_t* token, jsmntype_t type, int start, int end) {
    if (!token) {
        return;
    }
    token->type = type;
    token->start = start;
    token->end = end;
    token->size = 0;
}

static int jsmn_parse_primitive(jsmn_parser* parser, const char* js, size_t len,
                                jsmntok_t* tokens, unsigned int num_tokens) {
    int start = parser->pos;
    while (parser->pos < len) {
        char c = js[parser->pos];
        switch (c) {
            case '\t':
            case '\r':
            case '\n':
            case ' ':
            case ',':
            case ']':
            case '}':
                goto found;
            default:
                if (c < 32 || c == '"' || c == ':') {
                    return JSMN_ERROR_INVAL;
                }
                parser->pos++;
                break;
        }
    }
found:
    if (!tokens) {
        parser->pos--;
        return 0;
    }
    jsmntok_t* token = jsmn_alloc_token(parser, tokens, num_tokens);
    if (!token) {
        return JSMN_ERROR_NOMEM;
    }
    jsmn_fill_token(token, JSMN_PRIMITIVE, start, parser->pos);
    token->parent = parser->toksuper;
    parser->pos--;
    return 0;
}

static int jsmn_parse_string(jsmn_parser* parser, const char* js, size_t len,
                             jsmntok_t* tokens, unsigned int num_tokens) {
    int start = parser->pos;
    parser->pos++;

    while (parser->pos < len) {
        char c = js[parser->pos];
        if (c == '"') {
            if (!tokens) {
                return 0;
            }
            jsmntok_t* token = jsmn_alloc_token(parser, tokens, num_tokens);
            if (!token) {
                return JSMN_ERROR_NOMEM;
            }
            jsmn_fill_token(token, JSMN_STRING, start + 1, parser->pos);
            token->parent = parser->toksuper;
            return 0;
        }
        if (c == '\\') {
            parser->pos++;
            if (parser->pos == len) {
                return JSMN_ERROR_PART;
            }
            char esc = js[parser->pos];
            switch (esc) {
                case '\"':
                case '\\':
                case '/':
                case 'b':
                case 'f':
                case 'n':
                case 'r':
                case 't':
                    break;
                case 'u':
                    if (parser->pos + 4 >= len) {
                        return JSMN_ERROR_PART;
                    }
                    parser->pos += 4;
                    break;
                default:
                    return JSMN_ERROR_INVAL;
            }
        }
        parser->pos++;
    }
    return JSMN_ERROR_PART;
}

void jsmn_init(jsmn_parser* parser) {
    if (!parser) {
        return;
    }
    parser->pos = 0;
    parser->toknext = 0;
    parser->toksuper = -1;
}

int jsmn_parse(jsmn_parser* parser, const char* js, size_t len,
               jsmntok_t* tokens, unsigned int num_tokens) {
    if (!parser || !js) {
        return JSMN_ERROR_INVAL;
    }

    int r;
    for (; parser->pos < len; parser->pos++) {
        char c = js[parser->pos];
        switch (c) {
            case '{':
            case '[': {
                jsmntok_t* token = jsmn_alloc_token(parser, tokens, num_tokens);
                if (!token) {
                    return JSMN_ERROR_NOMEM;
                }
                if (parser->toksuper != -1) {
                    tokens[parser->toksuper].size++;
                }
                token->type = (c == '{') ? JSMN_OBJECT : JSMN_ARRAY;
                token->start = parser->pos;
                token->parent = parser->toksuper;
                parser->toksuper = parser->toknext - 1;
                break;
            }
            case '}':
            case ']': {
                if (!tokens) {
                    break;
                }
                jsmntype_t type = (c == '}') ? JSMN_OBJECT : JSMN_ARRAY;
                int i = parser->toknext - 1;
                for (; i >= 0; i--) {
                    jsmntok_t* tok = &tokens[i];
                    if (tok->start != -1 && tok->end == -1) {
                        if (tok->type != type) {
                            return JSMN_ERROR_INVAL;
                        }
                        tok->end = parser->pos + 1;
                        parser->toksuper = tok->parent;
                        break;
                    }
                }
                if (i < 0) {
                    return JSMN_ERROR_INVAL;
                }
                break;
            }
            case '"': {
                r = jsmn_parse_string(parser, js, len, tokens, num_tokens);
                if (r < 0) return r;
                if (parser->toksuper != -1 && tokens) {
                    tokens[parser->toksuper].size++;
                }
                break;
            }
            case '\t':
            case '\r':
            case '\n':
            case ' ': {
                break;
            }
            case ':': {
                parser->toksuper = parser->toknext - 1;
                break;
            }
            case ',': {
                if (parser->toksuper != -1 && tokens) {
                    int type = tokens[parser->toksuper].type;
                    if (type != JSMN_ARRAY && type != JSMN_OBJECT) {
                        parser->toksuper = tokens[parser->toksuper].parent;
                    }
                }
                break;
            }
            default: {
                r = jsmn_parse_primitive(parser, js, len, tokens, num_tokens);
                if (r < 0) return r;
                if (parser->toksuper != -1 && tokens) {
                    tokens[parser->toksuper].size++;
                }
                break;
            }
        }
    }

    if (tokens) {
        for (int i = parser->toknext - 1; i >= 0; i--) {
            if (tokens[i].start != -1 && tokens[i].end == -1) {
                return JSMN_ERROR_PART;
            }
        }
    }

    return (int)parser->toknext;
}
