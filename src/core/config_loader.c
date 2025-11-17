
// =============================================================================
// emwave-c: Runtime configuration loader
// JSON/CLI parsing that overlays onto SimulationConfig and delegates
// validation/clamping/summary to config_runtime.
// =============================================================================

#include "config_loader.h"

#include "jsmn.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define JSON_TOKEN_START 256

static void set_error(char* errbuf, size_t errbuf_len, const char* fmt, ...) {
    if (!errbuf || errbuf_len == 0 || !fmt) return;
    va_list args;
    va_start(args, fmt);
    vsnprintf(errbuf, errbuf_len, fmt, args);
    va_end(args);
}

static int parse_int_arg(const char* value, int* out) {
    if (!value || !out) return 0;
    char* end = NULL;
    long v = strtol(value, &end, 10);
    if (value == end) return 0;
    *out = (int)v;
    return 1;
}

static int parse_double_arg(const char* value, double* out) {
    if (!value || !out) return 0;
    char* end = NULL;
    double v = strtod(value, &end);
    if (value == end) return 0;
    *out = v;
    return 1;
}

static int str_ieq(const char* a, const char* b) {
    if (!a || !b) return 0;
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) {
            return 0;
        }
        a++;
        b++;
    }
    return *a == '\0' && *b == '\0';
}

static unsigned char material_tag_from_type(const char* type) {
    if (str_ieq(type, "pec")) return 1;
    if (str_ieq(type, "pmc")) return 2;
    return 0;
}

static SourceType source_type_from_string(const char* type) {
    if (str_ieq(type, "cw") || str_ieq(type, "continuous")) {
        return SRC_CW;
    }
    if (str_ieq(type, "gaussian") || str_ieq(type, "gauss")) {
        return SRC_GAUSS_PULSE;
    }
    if (str_ieq(type, "ricker")) {
        return SRC_RICKER;
    }
    return SRC_CW;
}

static int jsmn_eq(const char* json, const jsmntok_t* tok, const char* s) {
    if (!tok || tok->type != JSMN_STRING) return 0;
    size_t len = (size_t)(tok->end - tok->start);
    return strlen(s) == len && strncmp(json + tok->start, s, len) == 0;
}

static int json_token_to_string(const char* json, const jsmntok_t* tok, char* out, size_t out_len) {
    if (!tok || tok->type != JSMN_STRING || !out || out_len == 0) return 0;
    size_t len = (size_t)(tok->end - tok->start);
    if (len + 1 > out_len) {
        len = out_len - 1;
    }
    memcpy(out, json + tok->start, len);
    out[len] = '\0';
    return 1;
}

static int json_token_to_primitive(const char* json, const jsmntok_t* tok,
                                   char* buffer, size_t buf_len) {
    if (!tok || tok->type != JSMN_PRIMITIVE || !buffer || buf_len == 0) return 0;
    size_t len = (size_t)(tok->end - tok->start);
    if (len + 1 > buf_len) {
        len = buf_len - 1;
    }
    memcpy(buffer, json + tok->start, len);
    buffer[len] = '\0';
    return 1;
}

static int json_token_to_double(const char* json, const jsmntok_t* tok, double* out) {
    char tmp[64];
    if (!json_token_to_primitive(json, tok, tmp, sizeof(tmp))) return 0;
    char* end = NULL;
    double v = strtod(tmp, &end);
    if (end == tmp) return 0;
    if (out) *out = v;
    return 1;
}

static int json_token_to_int(const char* json, const jsmntok_t* tok, int* out) {
    char tmp[32];
    if (!json_token_to_primitive(json, tok, tmp, sizeof(tmp))) return 0;
    char* end = NULL;
    long v = strtol(tmp, &end, 10);
    if (end == tmp) return 0;
    if (out) *out = (int)v;
    return 1;
}

static int json_token_to_bool(const char* json, const jsmntok_t* tok, int* out) {
    char tmp[16];
    if (!json_token_to_primitive(json, tok, tmp, sizeof(tmp))) return 0;
    if (strcmp(tmp, "true") == 0) {
        if (out) *out = 1;
        return 1;
    }
    if (strcmp(tmp, "false") == 0) {
        if (out) *out = 0;
        return 1;
    }
    char* end = NULL;
    long v = strtol(tmp, &end, 10);
    if (end == tmp) return 0;
    if (out) *out = (v != 0);
    return 1;
}

static int json_skip(const jsmntok_t* tokens, int total, int index) {
    if (!tokens || index < 0 || index >= total) return total;
    int next = index + 1;
    const jsmntok_t* tok = &tokens[index];
    if (tok->type == JSMN_ARRAY) {
        for (int i = 0; i < tok->size; i++) {
            next = json_skip(tokens, total, next);
        }
        return next;
    }
    if (tok->type == JSMN_OBJECT) {
        for (int i = 0; i < tok->size; i++) {
            next = json_skip(tokens, total, next);
            next = json_skip(tokens, total, next);
        }
        return next;
    }
    return next;
}

static int json_object_find(const char* json, const jsmntok_t* tokens, int total,
                            int object_index, const char* key) {
    if (!json || !tokens || object_index < 0 || object_index >= total) return -1;
    const jsmntok_t* obj = &tokens[object_index];
    if (obj->type != JSMN_OBJECT) return -1;
    int idx = object_index + 1;
    for (int pair = 0; pair < obj->size; pair++) {
        int key_idx = idx;
        int value_idx = json_skip(tokens, total, key_idx);
        if (value_idx >= total) return -1;
        if (jsmn_eq(json, &tokens[key_idx], key)) {
            return value_idx;
        }
        idx = json_skip(tokens, total, value_idx);
        if (idx > total) return -1;
    }
    return -1;
}

static void json_apply_simulation(const char* json, const jsmntok_t* tokens, int total,
                                  int obj_index, SimulationConfig* cfg) {
    if (!cfg || obj_index < 0) return;
    int idx;
    if ((idx = json_object_find(json, tokens, total, obj_index, "nx")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->nx);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "ny")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->ny);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "lx")) >= 0) {
        json_token_to_double(json, &tokens[idx], &cfg->lx);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "ly")) >= 0) {
        json_token_to_double(json, &tokens[idx], &cfg->ly);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "cfl")) >= 0) {
        json_token_to_double(json, &tokens[idx], &cfg->cfl_safety);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "steps_per_frame")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->steps_per_frame);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "sweep_points")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->sweep_points);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "sweep_start_hz")) >= 0) {
        json_token_to_double(json, &tokens[idx], &cfg->sweep_start_hz);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "sweep_stop_hz")) >= 0) {
        json_token_to_double(json, &tokens[idx], &cfg->sweep_stop_hz);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "sweep_steps_per_point")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->sweep_steps_per_point);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "run_mode")) >= 0 &&
        tokens[idx].type == JSMN_STRING) {
        char mode[16];
        if (json_token_to_string(json, &tokens[idx], mode, sizeof(mode))) {
            if (str_ieq(mode, "sweep")) {
                cfg->run_mode = SIM_RUN_MODE_SWEEP;
            } else {
                cfg->run_mode = SIM_RUN_MODE_FIXED_STEPS;
            }
        }
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "run_steps")) >= 0) {
        json_token_to_int(json, &tokens[idx], &cfg->run_steps);
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "enable_probe_log")) >= 0) {
        int flag = 0;
        if (json_token_to_bool(json, &tokens[idx], &flag)) {
            cfg->enable_probe_log = flag ? 1 : 0;
        }
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "probe_log_path")) >= 0 &&
        tokens[idx].type == JSMN_STRING) {
        char path[SIM_PROBE_LOG_PATH_MAX];
        if (json_token_to_string(json, &tokens[idx], path, sizeof(path))) {
            strncpy(cfg->probe_log_path, path, SIM_PROBE_LOG_PATH_MAX - 1);
            cfg->probe_log_path[SIM_PROBE_LOG_PATH_MAX - 1] = '\0';
        }
    }
    if ((idx = json_object_find(json, tokens, total, obj_index, "boundary")) >= 0 &&
        tokens[idx].type == JSMN_STRING) {
        char bmode[16];
        if (json_token_to_string(json, &tokens[idx], bmode, sizeof(bmode))) {
            if (str_ieq(bmode, "mur")) {
                cfg->boundary_mode = SIM_BOUNDARY_MUR;
            } else {
                cfg->boundary_mode = SIM_BOUNDARY_CPML;
            }
        }
    }
}

static int json_load_material(const char* json, const jsmntok_t* tokens, int total,
                              int obj_index, MaterialRectSpec* out) {
    if (!out || obj_index < 0 || tokens[obj_index].type != JSMN_OBJECT) return 0;
    MaterialRectSpec spec = {0};
    spec.epsr = 1.0;
    spec.sigma = SIGMA_BG;
    spec.tag = 0;
    int idx = obj_index + 1;
    for (int pair = 0; pair < tokens[obj_index].size; pair++) {
        int key_idx = idx;
        int value_idx = json_skip(tokens, total, key_idx);
        char key[32];
        if (!json_token_to_string(json, &tokens[key_idx], key, sizeof(key))) {
            idx = json_skip(tokens, total, value_idx);
            continue;
        }
        if (strcmp(key, "type") == 0) {
            char type[16];
            if (json_token_to_string(json, &tokens[value_idx], type, sizeof(type))) {
                spec.tag = material_tag_from_type(type);
            }
        } else if (strcmp(key, "x0") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.x0);
        } else if (strcmp(key, "y0") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.y0);
        } else if (strcmp(key, "x1") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.x1);
        } else if (strcmp(key, "y1") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.y1);
        } else if (strcmp(key, "epsr") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.epsr);
        } else if (strcmp(key, "sigma") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.sigma);
        }
        idx = json_skip(tokens, total, value_idx);
    }
    *out = spec;
    return 1;
}

static int json_load_source(const char* json, const jsmntok_t* tokens, int total,
                            int obj_index, SourceConfigSpec* out) {
    if (!out || obj_index < 0 || tokens[obj_index].type != JSMN_OBJECT) return 0;
    SourceConfigSpec spec = {0};
    spec.active = 1;
    spec.amp = 1.0;
    spec.freq = 1e9;
    spec.sigma2 = 4.0;
    spec.type = SRC_CW;
    int idx = obj_index + 1;
    for (int pair = 0; pair < tokens[obj_index].size; pair++) {
        int key_idx = idx;
        int value_idx = json_skip(tokens, total, key_idx);
        char key[32];
        if (!json_token_to_string(json, &tokens[key_idx], key, sizeof(key))) {
            idx = json_skip(tokens, total, value_idx);
            continue;
        }
        if (strcmp(key, "type") == 0) {
            char type[16];
            if (json_token_to_string(json, &tokens[value_idx], type, sizeof(type))) {
                spec.type = source_type_from_string(type);
            }
        } else if (strcmp(key, "x") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.x);
        } else if (strcmp(key, "y") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.y);
        } else if (strcmp(key, "amp") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.amp);
        } else if (strcmp(key, "freq") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.freq);
        } else if (strcmp(key, "sigma2") == 0) {
            json_token_to_double(json, &tokens[value_idx], &spec.sigma2);
        } else if (strcmp(key, "active") == 0) {
            json_token_to_bool(json, &tokens[value_idx], &spec.active);
        }
        idx = json_skip(tokens, total, value_idx);
    }
    *out = spec;
    return 1;
}

static void json_apply_materials(const char* json, const jsmntok_t* tokens, int total,
                                 int arr_index, SimulationConfig* cfg) {
    if (!cfg || arr_index < 0) return;
    const jsmntok_t* arr = &tokens[arr_index];
    if (arr->type != JSMN_ARRAY) return;
    cfg->material_rect_count = 0;
    int idx = arr_index + 1;
    for (int i = 0; i < arr->size && cfg->material_rect_count < CONFIG_MAX_MATERIAL_RECTS; i++) {
        MaterialRectSpec spec;
        if (json_load_material(json, tokens, total, idx, &spec)) {
            cfg->material_rects[cfg->material_rect_count++] = spec;
        }
        idx = json_skip(tokens, total, idx);
    }
}

static void json_apply_sources(const char* json, const jsmntok_t* tokens, int total,
                               int arr_index, SimulationConfig* cfg) {
    if (!cfg || arr_index < 0) return;
    const jsmntok_t* arr = &tokens[arr_index];
    if (arr->type != JSMN_ARRAY) return;
    cfg->source_count = 0;
    int idx = arr_index + 1;
    for (int i = 0; i < arr->size && cfg->source_count < MAX_SRC; i++) {
        SourceConfigSpec spec;
        if (json_load_source(json, tokens, total, idx, &spec)) {
            cfg->source_configs[cfg->source_count++] = spec;
        }
        idx = json_skip(tokens, total, idx);
    }
}

static char* read_file(const char* path, size_t* out_len, char* errbuf, size_t errbuf_len) {
    if (!path) {
        set_error(errbuf, errbuf_len, "No config file specified");
        return NULL;
    }
    FILE* f = fopen(path, "rb");
    if (!f) {
        set_error(errbuf, errbuf_len, "Unable to open %s", path);
        return NULL;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        set_error(errbuf, errbuf_len, "Failed to seek in %s", path);
        fclose(f);
        return NULL;
    }
    long len = ftell(f);
    if (len < 0) {
        set_error(errbuf, errbuf_len, "Failed to determine size of %s", path);
        fclose(f);
        return NULL;
    }
    if ((unsigned long)len > (unsigned long)CONFIG_LOADER_MAX_FILE_BYTES) {
        double limit_mb = (double)CONFIG_LOADER_MAX_FILE_BYTES / (1024.0 * 1024.0);
        set_error(errbuf, errbuf_len,
                  "Config file %s is too large (%.2f MB limit)", path, limit_mb);
        fclose(f);
        return NULL;
    }
    if (fseek(f, 0, SEEK_SET) != 0) {
        set_error(errbuf, errbuf_len, "Failed to rewind %s", path);
        fclose(f);
        return NULL;
    }
    size_t alloc = (size_t)len;
    char* buf = (char*)malloc(alloc + 1);
    if (!buf) {
        set_error(errbuf, errbuf_len, "Out of memory reading %s", path);
        fclose(f);
        return NULL;
    }
    size_t read_len = fread(buf, 1, alloc, f);
    fclose(f);
    if (read_len != alloc) {
        set_error(errbuf, errbuf_len, "Failed to read config file %s", path);
        free(buf);
        return NULL;
    }
    buf[len] = '\0';
    if (out_len) *out_len = (size_t)len;
    return buf;
}

static int config_overlay_from_json(const char* path, SimulationConfig* cfg,
                                    char* errbuf, size_t errbuf_len) {
    if (!path || !cfg) {
        set_error(errbuf, errbuf_len, "Invalid arguments");
        return 0;
    }
    size_t len = 0;
    char* data = read_file(path, &len, errbuf, errbuf_len);
    if (!data) {
        if (!errbuf || errbuf_len == 0 || errbuf[0] == '\0') {
            fprintf(stderr, "Failed to read config file %s\n", path);
        }
        return 0;
    }

    size_t tok_capacity = JSON_TOKEN_START;
    jsmntok_t* tokens = NULL;
    int parsed_tokens = 0;
    while (1) {
        jsmn_parser parser;
        jsmn_init(&parser);
        jsmntok_t* tmp = (jsmntok_t*)realloc(tokens, sizeof(jsmntok_t) * tok_capacity);
        if (!tmp) {
            free(tokens);
            free(data);
            return 0;
        }
        tokens = tmp;
        parsed_tokens = jsmn_parse(&parser, data, len, tokens, (unsigned int)tok_capacity);
        if (parsed_tokens == JSMN_ERROR_NOMEM) {
            tok_capacity *= 2;
            continue;
        }
        if (parsed_tokens < 0) {
            if (errbuf && errbuf_len) {
                set_error(errbuf, errbuf_len,
                          "Failed to parse JSON in %s (error %d)", path, parsed_tokens);
            } else {
                fprintf(stderr, "Failed to parse JSON in %s (error %d)\n", path, parsed_tokens);
            }
            free(tokens);
            free(data);
            return 0;
        }
        break;
    }

    if (parsed_tokens == 0 || tokens[0].type != JSMN_OBJECT) {
        if (errbuf && errbuf_len) {
            set_error(errbuf, errbuf_len, "Config file %s must contain a JSON object", path);
        } else {
            fprintf(stderr, "Config file %s must contain a JSON object\n", path);
        }
        free(tokens);
        free(data);
        return 0;
    }

    int sim_idx = json_object_find(data, tokens, parsed_tokens, 0, "simulation");
    if (sim_idx >= 0) {
        json_apply_simulation(data, tokens, parsed_tokens, sim_idx, cfg);
    }
    int mats_idx = json_object_find(data, tokens, parsed_tokens, 0, "materials");
    if (mats_idx >= 0) {
        json_apply_materials(data, tokens, parsed_tokens, mats_idx, cfg);
    }
    int src_idx = json_object_find(data, tokens, parsed_tokens, 0, "sources");
    if (src_idx >= 0) {
        json_apply_sources(data, tokens, parsed_tokens, src_idx, cfg);
    }

    free(tokens);
    free(data);
    return 1;
}

static void print_usage(void) {
    printf("Options:\n");
    printf("  --config=<file>       Load configuration JSON file\n");
    printf("  --nx=<cells>          Grid cells in X (default %d)\n", NX_DEFAULT);
    printf("  --ny=<cells>          Grid cells in Y (default %d)\n", NY_DEFAULT);
    printf("  --lx=<meters>         Domain length in X (default %.3f m)\n", LX_DEFAULT);
    printf("  --ly=<meters>         Domain length in Y (default %.3f m)\n", LY_DEFAULT);
    printf("  --cfl=<0-1>           CFL safety factor (default %.2f)\n", CFL_SAFETY_FACTOR);
    printf("  --sweep-points=<n>    Number of sweep points (default %d)\n", SIM_CONFIG_DEFAULTS.sweep_points);
    printf("  --sweep-start=<Hz>    Sweep start frequency\n");
    printf("  --sweep-stop=<Hz>     Sweep stop frequency\n");
    printf("  --sweep-steps=<n>     Steps per sweep point\n");
    printf("  --run-mode=MODE       fixed or sweep (default fixed)\n");
    printf("  --run-steps=<n>       Steps to run when mode=fixed\n");
    printf("  --boundary=MODE       cpml or mur (default cpml)\n");
    printf("  --probe-log=<path>    Enable probe logging to file\n");
    printf("  --no-probe-log        Disable probe logging\n");
    printf("  --help                Show this help\n");
}

static const char* extract_config_path(int argc, char** argv) {
    const char* path = NULL;
    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        if (strncmp(arg, "--config=", 9) == 0) {
            path = arg + 9;
        } else if (strcmp(arg, "--config") == 0 && i + 1 < argc) {
            path = argv[i + 1];
            i++;
        }
    }
    return path;
}

int config_load_from_args(int argc, char** argv, SimulationConfig* out_config) {
    if (!out_config) return 0;
    *out_config = SIM_CONFIG_DEFAULTS;

    const char* config_path = extract_config_path(argc, argv);
    if (config_path && config_path[0] != '\0') {
        char load_err[256];
        if (!config_loader_parse_file(config_path, out_config, load_err, sizeof(load_err))) {
            fprintf(stderr, "Failed to load config file %s: %s\n", config_path,
                    (load_err[0] != '\0') ? load_err : "unknown error");
            return 0;
        }
    }

    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        if (strcmp(arg, "--help") == 0) {
            print_usage();
            return 0;
        } else if (strncmp(arg, "--config", 8) == 0) {
            if (strcmp(arg, "--config") == 0) {
                i++;
            }
            continue;
        } else if (strncmp(arg, "--nx=", 5) == 0) {
            parse_int_arg(arg + 5, &out_config->nx);
        } else if (strncmp(arg, "--ny=", 5) == 0) {
            parse_int_arg(arg + 5, &out_config->ny);
        } else if (strncmp(arg, "--lx=", 5) == 0) {
            parse_double_arg(arg + 5, &out_config->lx);
        } else if (strncmp(arg, "--ly=", 5) == 0) {
            parse_double_arg(arg + 5, &out_config->ly);
        } else if (strncmp(arg, "--cfl=", 6) == 0) {
            parse_double_arg(arg + 6, &out_config->cfl_safety);
        } else if (strncmp(arg, "--sweep-points=", 15) == 0) {
            parse_int_arg(arg + 15, &out_config->sweep_points);
        } else if (strncmp(arg, "--sweep-start=", 14) == 0) {
            parse_double_arg(arg + 14, &out_config->sweep_start_hz);
        } else if (strncmp(arg, "--sweep-stop=", 13) == 0) {
            parse_double_arg(arg + 13, &out_config->sweep_stop_hz);
        } else if (strncmp(arg, "--sweep-steps=", 14) == 0) {
            parse_int_arg(arg + 14, &out_config->sweep_steps_per_point);
        } else if (strncmp(arg, "--run-mode=", 11) == 0) {
            const char* mode = arg + 11;
            if (strcmp(mode, "sweep") == 0) {
                out_config->run_mode = SIM_RUN_MODE_SWEEP;
            } else {
                out_config->run_mode = SIM_RUN_MODE_FIXED_STEPS;
            }
        } else if (strncmp(arg, "--run-steps=", 12) == 0) {
            parse_int_arg(arg + 12, &out_config->run_steps);
        } else if (strncmp(arg, "--boundary=", 11) == 0) {
            const char* mode = arg + 11;
            if (strcmp(mode, "mur") == 0) {
                out_config->boundary_mode = SIM_BOUNDARY_MUR;
            } else {
                out_config->boundary_mode = SIM_BOUNDARY_CPML;
            }
        } else if (strncmp(arg, "--probe-log=", 12) == 0) {
            const char* path = arg + 12;
            if (*path) {
                strncpy(out_config->probe_log_path, path, SIM_PROBE_LOG_PATH_MAX - 1);
                out_config->probe_log_path[SIM_PROBE_LOG_PATH_MAX - 1] = '\0';
                out_config->enable_probe_log = 1;
            }
        } else if (strcmp(arg, "--no-probe-log") == 0) {
            out_config->enable_probe_log = 0;
        }
    }

    config_clamp_to_limits(out_config);
    char errbuf[128];
    if (!config_validate(out_config, errbuf, sizeof(errbuf))) {
        fprintf(stderr, "Invalid configuration: %s\n", errbuf);
        return 0;
    }

    return 1;
}

int config_loader_parse_file(const char* path, SimulationConfig* cfg,
                             char* errbuf, size_t errbuf_len) {
    if (errbuf && errbuf_len) {
        errbuf[0] = '\0';
    }
    if (!cfg) {
        set_error(errbuf, errbuf_len, "No configuration struct provided");
        return 0;
    }
    return config_overlay_from_json(path, cfg, errbuf, errbuf_len);
}
