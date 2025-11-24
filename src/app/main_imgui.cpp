// ============================================================================
// emwave-c: Dear ImGui Front-End with Simulation Wizard
// ============================================================================

#include "config.h"
#include "types.h"
#include "analysis.h"
#include "app_bootstrap.h"
#include "fdtd_core.h"
#include "ui_render.h"
#include "materials.h"
#include "boundary.h"
#include "sources.h"
#include "expr.h"
#include "core/material_library.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"
#include "implot.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cctype>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <system_error>
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>
#ifdef _WIN32
#include <process.h>
#endif
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <shellapi.h>
#include <direct.h>
#endif

// Minimal PNG write (embedded, trimmed from stb_image_write for PNG-only)
#include <cassert>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#define STBIW_MALLOC(sz)           malloc(sz)
#define STBIW_FREE(p)              free(p)
#define STBIW_MEMMOVE(a,b,sz)      memmove(a,b,sz)
#define STBIW_MEMCPY(a,b,sz)       memcpy(a,b,sz)
#define STBIW_ASSERT(x)            assert(x)
#define STBIW_MIN(a,b)             ((a) < (b) ? (a) : (b))
#define STBIW_MAX(a,b)             ((a) > (b) ? (a) : (b))
#define STBIW_UCHAR(x)             (unsigned char)(x)
#define STBIW_ZLIB_QUALITY 8

static void _stbiw_write_bytes(FILE* f, const void* data, int len) {
    if (!f || !data || len <= 0) return;
    fwrite(data, 1, (size_t)len, f);
}

// CRC
static unsigned int stb_crc32(unsigned char* buffer, int len) {
    static unsigned int crc_table[256];
    static int init = 0;
    if (!init) {
        for (unsigned int i = 0; i < 256; i++) {
            unsigned int c = i;
            for (int j = 0; j < 8; j++) {
                c = (c & 1) ? 0xedb88320 ^ (c >> 1) : c >> 1;
            }
            crc_table[i] = c;
        }
        init = 1;
    }
    unsigned int c = 0xffffffff;
    for (int i = 0; i < len; ++i)
        c = crc_table[(c ^ buffer[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xffffffff;
}

// zlib (very small, not optimized)
static unsigned char* _stbiw__zlib_compress(unsigned char* data, int data_len, int* out_len, int quality) {
    (void)quality;
    // naive: store uncompressed (zlib with no compression)
    int header_size = 2;
    int adler_size = 4;
    int block_overhead = 5; // final block header + len + nlen
    int needed = header_size + block_overhead + data_len + adler_size;
    unsigned char* out = (unsigned char*)STBIW_MALLOC(needed);
    if (!out) return nullptr;
    unsigned char* p = out;
    *p++ = 0x78; // CMF
    *p++ = 0x01; // FLG (no compression)
    // single uncompressed block
    *p++ = 1; // final block, uncompressed
    unsigned short len = (unsigned short)data_len;
    unsigned short nlen = ~len;
    *p++ = (unsigned char)(len & 0xFF);
    *p++ = (unsigned char)(len >> 8);
    *p++ = (unsigned char)(nlen & 0xFF);
    *p++ = (unsigned char)(nlen >> 8);
    STBIW_MEMCPY(p, data, data_len);
    p += data_len;
    // Adler-32
    unsigned int s1 = 1, s2 = 0;
    for (int i = 0; i < data_len; ++i) {
        s1 = (s1 + data[i]) % 65521;
        s2 = (s2 + s1) % 65521;
    }
    unsigned int adler = (s2 << 16) | s1;
    *p++ = (unsigned char)((adler >> 24) & 0xFF);
    *p++ = (unsigned char)((adler >> 16) & 0xFF);
    *p++ = (unsigned char)((adler >> 8) & 0xFF);
    *p++ = (unsigned char)(adler & 0xFF);
    *out_len = (int)(p - out);
    return out;
}

static void _stbiw_write32be(FILE* f, unsigned int v) {
    unsigned char b[4];
    b[0] = (unsigned char)((v >> 24) & 0xFF);
    b[1] = (unsigned char)((v >> 16) & 0xFF);
    b[2] = (unsigned char)((v >> 8) & 0xFF);
    b[3] = (unsigned char)(v & 0xFF);
    _stbiw_write_bytes(f, b, 4);
}

static void _stbiw_write_chunk(FILE* f, const char* tag, unsigned char* data, int len) {
    unsigned char len_bytes[4];
    len_bytes[0] = (unsigned char)((len >> 24) & 0xFF);
    len_bytes[1] = (unsigned char)((len >> 16) & 0xFF);
    len_bytes[2] = (unsigned char)((len >> 8) & 0xFF);
    len_bytes[3] = (unsigned char)(len & 0xFF);
    _stbiw_write_bytes(f, len_bytes, 4);
    _stbiw_write_bytes(f, tag, 4);
    if (data) _stbiw_write_bytes(f, data, len);
    unsigned char crc_src[4096];
    int crc_len = len + 4;
    if (crc_len > (int)sizeof(crc_src)) {
        unsigned char* tmp = (unsigned char*)STBIW_MALLOC(crc_len);
        STBIW_MEMCPY(tmp, tag, 4);
        if (data) STBIW_MEMCPY(tmp + 4, data, len);
        unsigned int crc = stb_crc32(tmp, crc_len);
        STBIW_FREE(tmp);
        unsigned char crc_bytes[4];
        crc_bytes[0] = (unsigned char)((crc >> 24) & 0xFF);
        crc_bytes[1] = (unsigned char)((crc >> 16) & 0xFF);
        crc_bytes[2] = (unsigned char)((crc >> 8) & 0xFF);
        crc_bytes[3] = (unsigned char)(crc & 0xFF);
        _stbiw_write_bytes(f, crc_bytes, 4);
    } else {
        STBIW_MEMCPY(crc_src, tag, 4);
        if (data) STBIW_MEMCPY(crc_src + 4, data, len);
        unsigned int crc = stb_crc32(crc_src, crc_len);
        unsigned char crc_bytes[4];
        crc_bytes[0] = (unsigned char)((crc >> 24) & 0xFF);
        crc_bytes[1] = (unsigned char)((crc >> 16) & 0xFF);
        crc_bytes[2] = (unsigned char)((crc >> 8) & 0xFF);
        crc_bytes[3] = (unsigned char)(crc & 0xFF);
        _stbiw_write_bytes(f, crc_bytes, 4);
    }
}

static int stbi_write_png(const char* filename, int w, int h, int comp, const void* data, int stride_bytes) {
    if (!filename || !data || w <= 0 || h <= 0) return 0;
    if (comp != 3 && comp != 4) return 0;
    FILE* f = fopen(filename, "wb");
    if (!f) return 0;

    static const unsigned char sig[8] = {137,80,78,71,13,10,26,10};
    _stbiw_write_bytes(f, sig, 8);

    unsigned char ihdr[13];
    ihdr[0] = (unsigned char)((w >> 24) & 0xFF);
    ihdr[1] = (unsigned char)((w >> 16) & 0xFF);
    ihdr[2] = (unsigned char)((w >> 8) & 0xFF);
    ihdr[3] = (unsigned char)(w & 0xFF);
    ihdr[4] = (unsigned char)((h >> 24) & 0xFF);
    ihdr[5] = (unsigned char)((h >> 16) & 0xFF);
    ihdr[6] = (unsigned char)((h >> 8) & 0xFF);
    ihdr[7] = (unsigned char)(h & 0xFF);
    ihdr[8] = 8; // bit depth
    ihdr[9] = (unsigned char)(comp == 3 ? 2 : 6); // color type
    ihdr[10] = 0; // compression
    ihdr[11] = 0; // filter
    ihdr[12] = 0; // interlace
    _stbiw_write_chunk(f, "IHDR", ihdr, 13);

    // Pack data with filter byte per row (0)
    int stride = stride_bytes ? stride_bytes : (w * comp);
    int raw_len = (stride + 1) * h;
    unsigned char* raw = (unsigned char*)STBIW_MALLOC(raw_len);
    if (!raw) { fclose(f); return 0; }
    const unsigned char* src = (const unsigned char*)data;
    for (int y = 0; y < h; ++y) {
        unsigned char* row = raw + y * (stride + 1);
        row[0] = 0; // filter type 0
        STBIW_MEMCPY(row + 1, src + y * stride, stride);
    }

    int zlen = 0;
    unsigned char* zdata = _stbiw__zlib_compress(raw, raw_len, &zlen, STBIW_ZLIB_QUALITY);
    STBIW_FREE(raw);
    if (!zdata) { fclose(f); return 0; }
    _stbiw_write_chunk(f, "IDAT", zdata, zlen);
    STBIW_FREE(zdata);

    _stbiw_write_chunk(f, "IEND", nullptr, 0);
    fclose(f);
    return 1;
}

enum ViewportLayout {
    VIEWPORT_SINGLE = 0,
    VIEWPORT_HORIZONTAL = 1,
    VIEWPORT_VERTICAL = 2,
    VIEWPORT_QUAD = 3
};

enum ViewportViz {
    VIEWPORT_VIZ_EZ = 0,
    VIEWPORT_VIZ_EZ_ABS = 1,
    VIEWPORT_VIZ_HX = 2,
    VIEWPORT_VIZ_HY = 3,
    VIEWPORT_VIZ_HMAG = 4,
    VIEWPORT_VIZ_SX = 5,
    VIEWPORT_VIZ_SY = 6,
    VIEWPORT_VIZ_S_MAG = 7,
    VIEWPORT_VIZ_EX = 8,
    VIEWPORT_VIZ_EY = 9,
    VIEWPORT_VIZ_HZ = 10,
    VIEWPORT_VIZ_MATERIAL = 11,
    VIEWPORT_VIZ_OVERLAY = 12
};

struct ViewportInstance {
    ImVec2 pos;
    ImVec2 size;
    float zoom;
    float pan_x;
    float pan_y;
    ViewportViz viz_mode;
    bool active;
    bool valid;
    bool show_grid;
    bool show_sources;
    bool show_vectors;
};

// Forward declarations for helpers defined later in this file (need types above)
static ImVec2 compute_viewport_offset(const ViewportInstance& vp,
                                      const SimulationState* sim,
                                      int scale);
static int compute_grid_step_from_scale(int scale);
static void render_material_distribution(RenderContext* render,
                                         const SimulationState* sim,
                                         int scale);
static void render_material_overlay(RenderContext* render,
                                    const SimulationState* sim,
                                    int scale,
                                    float alpha);
static void render_material_outlines(RenderContext* render,
                                     const SimulationState* sim,
                                     int scale);
static void render_grid_overlay(RenderContext* render,
                                const SimulationState* sim,
                                SDL_Color color,
                                int grid_step);
static int compute_scope_fft(const Scope* scope,
                             double dt,
                             double* freq,
                             double* mag,
                             double* phase,
                             int max_fft);
static const std::string& ffmpeg_command();
static bool has_ffmpeg();
static void log_ffmpeg_attempt(const char* cmd, int rc);
static int run_ffmpeg(const std::vector<std::string>& args);
static std::string join_args(const std::vector<std::string>& args);

enum ComposerItemType {
    COMPOSER_FIELD_VIEW = 0,
    COMPOSER_REGION = 1,
    COMPOSER_SCOPE = 2,
    COMPOSER_FFT = 3,
    COMPOSER_LEGEND = 4,
    COMPOSER_MEAS = 5,
    COMPOSER_SMITH = 6
};

struct ComposerItem {
    int id;
    ComposerItemType type;
    int viewport_idx;   // for field items
    ImVec2 pos;         // page space (px)
    ImVec2 size;        // page space (px)
    bool selected;
    ImVec4 region_norm; // x0,y0,x1,y1 normalized (for region captures)
};

struct ComposerPage {
    char name[32];
    int res_w;
    int res_h;
    ImVec4 bg;
    bool transparent_bg;
    char output_name[64];
    int output_format;  // 0=BMP,1=PNG seq,2=MP4,3=GIF
    int fps;
    int frames;
    int video_kbps;     // target video bitrate for MP4
    std::vector<ComposerItem> items;
};

struct HeadlessComposerOpts {
    bool enabled;
    int fmt;        // 0 bmp,1 png,2 mp4,3 gif
    int fps;
    int frames;
    int res_w;
    int res_h;
    int tmpl_idx;   // -1 = none
    char output_name[64];
    char layout_path[260];
};

struct DistanceMeasurement {
    ImVec2 a;
    ImVec2 b;
    double distance_m;
    double angle_deg;
};

struct AreaMeasurement {
    std::vector<ImVec2> vertices;
    bool closed;
    double area_m2;
    double perimeter_m;
};

struct Annotation {
    ImVec2 grid_pos;
    char text[128];
    ImVec4 color;
    float font_size;
    bool visible;
};

struct MeasurementHistory {
    std::vector<DistanceMeasurement> distances;
    std::vector<AreaMeasurement> areas;
    std::vector<Annotation> annotations;
};

struct AppState {
    bool basic_mode;
    int selected_source;
    int selected_block;
    bool placing_source;
    bool placing_block;
    bool block_first_set;
    int block_first_i;
    int block_first_j;
    int last_click_i;
    int last_click_j;
    bool show_scope_window;
    bool show_scene_panel;
    bool show_sources_panel;
    bool show_blocks_panel;
    bool show_run_panel;
    bool show_grid_panel;
    bool show_run_settings_panel;
    bool show_probes_panel;
    bool show_log_panel;
    bool show_expression_panel;
    bool show_scenes_panel;
    bool show_grid_overlay;
    int theme_preset;
    ImVec4 accent_color;
    /* Simple log buffer */
    char log_lines[128][256];
    int log_count;
    ImVec2 viewport_pos;
    ImVec2 viewport_size;
    bool viewport_valid;
    /* Viewport transform */
    float viewport_zoom;
    float viewport_pan_x;
    float viewport_pan_y;
    bool viewport_panning;
    ImVec2 toolbar_screen_min;
    ImVec2 toolbar_screen_size;
    bool toolbar_valid;
    ImVec2 pan_start_mouse;
    ImVec2 pan_start_offset;
    float hud_zoom_value;
    float hud_zoom_timer;
    ImVec2 hud_pan_value;
    float hud_pan_timer;
    bool show_axis_overlay;
    bool ruler_mode;
    bool ruler_first_point_set;
    ImVec2 ruler_point_a;
    ImVec2 ruler_point_b;
    bool show_context_menu;
    ImVec2 context_menu_pos;
    int context_menu_cell_i;
    int context_menu_cell_j;
    /* Window state */
    int window_width;
    int window_height;
    /* Material browser */
    int selected_material_id;
    int paint_material_id;
    bool material_browser_open;
    char material_search[64];
    bool filter_metals;
    bool filter_dielectrics;
    /* Material visualization */
    bool show_material_legend;
    bool filter_blocks_by_material;
    bool auto_filter_blocks_on_select;
    bool highlight_blocks_by_material;
    int visualization_mode;   // 0=Field,1=Material,2=Overlay
    bool show_material_outlines;
    float material_overlay_alpha;
    std::vector<double> sparam_freq_hz;
    std::vector<double> sparam_freq_ghz;
    std::vector<double> sparam_s21_mag;
    std::vector<double> sparam_s21_db;
    std::vector<double> sparam_vswr;
    std::vector<double> sparam_s11_mag;
    std::vector<double> sparam_s11_db;
    bool sparam_window_open;
    bool sparam_data_loaded;
    double sparam_peak_db;
    double sparam_peak_freq_hz;
    double sparam_f_low_hz;
    double sparam_f_high_hz;
    double sparam_f_center_hz;
    double sparam_bandwidth_hz;
    char sparam_csv_path[260];
    bool smith_chart_open;
    bool smith_show_impedance;
    bool smith_show_vswr_circles;
    double smith_z0;
    int smith_freq_index;
    bool request_rebootstrap;
    char rebootstrap_message[128];

    // Multi-viewport state
    ViewportLayout viewport_layout;
    ViewportInstance viewports[4];
    int active_viewport_idx;
    bool sync_zoom;
    bool sync_pan;

    // Measurements (Prompt #39)
    bool area_mode;
    AreaMeasurement current_area;
    bool annotation_mode;
    Annotation temp_annotation;
    MeasurementHistory measurements;
    bool show_measurement_history;
    bool show_distance_measurements;
    bool show_area_measurements;
    bool show_annotations;

    // Print Composer
    bool show_print_composer;
    std::vector<ComposerPage> composer_pages;
    int composer_active_page;
    int composer_next_item_id;
    char composer_status[128];
    char composer_last_export_path[260];
    SDL_Texture* composer_preview_tex;
    ImVec2 composer_preview_tex_size;
    std::vector<SDL_Texture*> composer_preview_frames;
    int composer_preview_frame_idx;
    int composer_preview_frame_count;
    float composer_preview_fps;
    float composer_preview_time;
    bool composer_preview_playing;
    int composer_preview_max_frames;
    int composer_preview_target_width;
    bool composer_preview_mutate_sim;
    std::vector<ComposerPage> composer_user_templates;
    bool composer_dragging;
    int composer_drag_item;
    ImVec2 composer_drag_offset;
    bool composer_resizing;
    int composer_resize_item;
    int composer_resize_corner;
    ImVec2 composer_resize_start_pos;
    ImVec2 composer_resize_start_size;
    ImVec2 composer_resize_start_mouse;
    bool composer_request_export;
    int composer_request_page;
    bool composer_request_export_all;
    bool composer_request_animation;
    int composer_animation_steps_per_frame;
    bool composer_snap;
    float composer_snap_step;
    bool composer_canvas_grid;
    int composer_canvas_grid_step;
    // Region pick state
    bool composer_pick_region_active;
    bool composer_pick_region_dragging;
    int composer_pick_region_viewport;
    int composer_pick_region_page;
    int composer_pick_region_target_item;
    ImVec2 composer_pick_region_start_norm;
    ImVec2 composer_pick_region_end_norm;
    bool composer_pick_region_snap;
};

static void ui_log_add(AppState* app, const char* fmt, ...);
static void clear_composer_preview(AppState* app);
static SimulationState* clone_simulation_state(const SimulationState* src);
static bool composer_load_layout_json(AppState* app, const char* path);

static double calculate_area_m2(const AreaMeasurement& a, const SimulationState* sim) {
    if (!sim || a.vertices.size() < 3) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        size_t j = (i + 1) % a.vertices.size();
        sum += (double)a.vertices[i].x * (double)a.vertices[j].y;
        sum -= (double)a.vertices[j].x * (double)a.vertices[i].y;
    }
    double area_cells = std::abs(sum) * 0.5;
    double dx = sim->lx / (double)sim->nx;
    double dy = sim->ly / (double)sim->ny;
    return area_cells * dx * dy;
}

static double calculate_perimeter_m(const AreaMeasurement& a, const SimulationState* sim) {
    if (!sim || a.vertices.size() < 2) return 0.0;
    double per_cells = 0.0;
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        size_t j = (i + 1) % a.vertices.size();
        double dx = (double)a.vertices[j].x - (double)a.vertices[i].x;
        double dy = (double)a.vertices[j].y - (double)a.vertices[i].y;
        per_cells += std::sqrt(dx * dx + dy * dy);
    }
    double dx_m = sim->lx / (double)sim->nx;
    double dy_m = sim->ly / (double)sim->ny;
    double avg_cell = 0.5 * (dx_m + dy_m);
    return per_cells * avg_cell;
}

static void close_area_measurement(AppState* app, const SimulationState* sim) {
    if (!app || !sim) return;
    if (app->current_area.vertices.size() < 3) return;
    app->current_area.closed = true;
    app->current_area.area_m2 = calculate_area_m2(app->current_area, sim);
    app->current_area.perimeter_m = calculate_perimeter_m(app->current_area, sim);
    app->measurements.areas.push_back(app->current_area);
    ui_log_add(app, "Area: %.6f m^2 | Perimeter: %.4f m",
               app->current_area.area_m2,
               app->current_area.perimeter_m);
    app->current_area.vertices.clear();
    app->current_area.closed = false;
    app->area_mode = false;
}

static void export_measurements_csv(const AppState* app, const char* path) {
    if (!app || !path) return;
    FILE* f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "Type,ID,Value,Unit,Data\n");
    for (size_t i = 0; i < app->measurements.distances.size(); ++i) {
        const auto& d = app->measurements.distances[i];
        fprintf(f, "Distance,%zu,%.6f,m,\"(%.1f,%.1f)-(%.1f,%.1f) @ %.1f deg\"\n",
                i + 1,
                d.distance_m,
                d.a.x, d.a.y,
                d.b.x, d.b.y,
                d.angle_deg);
    }
    for (size_t i = 0; i < app->measurements.areas.size(); ++i) {
        const auto& a = app->measurements.areas[i];
        fprintf(f, "Area,%zu,%.6f,m2,\"%zu vertices | Perimeter=%.4f m\"\n",
                i + 1,
                a.area_m2,
                a.vertices.size(),
                a.perimeter_m);
    }
    for (size_t i = 0; i < app->measurements.annotations.size(); ++i) {
        const auto& ann = app->measurements.annotations[i];
        fprintf(f, "Annotation,%zu,,,\"(%.1f,%.1f): %s\"\n",
                i + 1,
                ann.grid_pos.x,
                ann.grid_pos.y,
                ann.text);
    }
    fclose(f);
}

static void composer_add_item(AppState* app,
                              ComposerPage& page,
                              ComposerItemType type,
                              int viewport_idx,
                              ImVec2 pos,
                              ImVec2 size) {
    ComposerItem item{};
    item.id = app->composer_next_item_id++;
    item.type = type;
    item.viewport_idx = viewport_idx;
    item.pos = pos;
    item.size = size;
    item.selected = false;
    item.region_norm = ImVec4(0.0f, 0.0f, 1.0f, 1.0f);
    page.items.push_back(item);
}

static void composer_apply_template(AppState* app, ComposerPage& page, int tmpl_idx) {
    if (!app) return;
    page.items.clear();
    float w = (float)page.res_w;
    float h = (float)page.res_h;
    switch (tmpl_idx) {
        case 0: {  // Blank
        } break;
        case 1: {  // Field + Legend
            composer_add_item(app, page, COMPOSER_FIELD_VIEW, 0, ImVec2(w * 0.05f, h * 0.1f), ImVec2(w * 0.65f, h * 0.75f));
            composer_add_item(app, page, COMPOSER_LEGEND, 0, ImVec2(w * 0.72f, h * 0.1f), ImVec2(w * 0.22f, h * 0.3f));
        } break;
        case 2: {  // Field + Scope
            composer_add_item(app, page, COMPOSER_FIELD_VIEW, 0, ImVec2(w * 0.05f, h * 0.05f), ImVec2(w * 0.9f, h * 0.6f));
            composer_add_item(app, page, COMPOSER_SCOPE, 0, ImVec2(w * 0.05f, h * 0.68f), ImVec2(w * 0.42f, h * 0.25f));
            composer_add_item(app, page, COMPOSER_FFT, 0, ImVec2(w * 0.53f, h * 0.68f), ImVec2(w * 0.42f, h * 0.25f));
        } break;
        case 3: {  // Field + FFT + Legend
            composer_add_item(app, page, COMPOSER_FIELD_VIEW, 0, ImVec2(w * 0.05f, h * 0.08f), ImVec2(w * 0.65f, h * 0.55f));
            composer_add_item(app, page, COMPOSER_FFT, 0, ImVec2(w * 0.05f, h * 0.66f), ImVec2(w * 0.65f, h * 0.28f));
            composer_add_item(app, page, COMPOSER_LEGEND, 0, ImVec2(w * 0.73f, h * 0.08f), ImVec2(w * 0.22f, h * 0.3f));
        } break;
        default:
            break;
    }
}

static void composer_select_all(ComposerPage& page, bool select) {
    for (auto& it : page.items) it.selected = select;
}

static void composer_align_selected(ComposerPage& page,
                                    int page_w,
                                    int page_h,
                                    const char* mode,
                                    bool snap,
                                    float snap_step) {
    std::vector<ComposerItem*> sel;
    for (auto& it : page.items) {
        if (it.selected) sel.push_back(&it);
    }
    if (sel.empty()) return;

    auto clamp_to_page = [&](ComposerItem* it) {
        if (!it) return;
        if (it->pos.x < 0.0f) it->pos.x = 0.0f;
        if (it->pos.y < 0.0f) it->pos.y = 0.0f;
        if (it->pos.x + it->size.x > page_w) it->pos.x = page_w - it->size.x;
        if (it->pos.y + it->size.y > page_h) it->pos.y = page_h - it->size.y;
    };

    auto snap_pos = [&](ImVec2 v) {
        if (!snap) return v;
        float step = std::max(1.0f, snap_step);
        v.x = std::round(v.x / step) * step;
        v.y = std::round(v.y / step) * step;
        return v;
    };

    if (strcmp(mode, "left") == 0) {
        float min_x = sel[0]->pos.x;
        for (auto* it : sel) min_x = std::min(min_x, it->pos.x);
        for (auto* it : sel) {
            it->pos.x = min_x;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "right") == 0) {
        float max_r = sel[0]->pos.x + sel[0]->size.x;
        for (auto* it : sel) max_r = std::max(max_r, it->pos.x + it->size.x);
        for (auto* it : sel) {
            it->pos.x = max_r - it->size.x;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "top") == 0) {
        float min_y = sel[0]->pos.y;
        for (auto* it : sel) min_y = std::min(min_y, it->pos.y);
        for (auto* it : sel) {
            it->pos.y = min_y;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "bottom") == 0) {
        float max_b = sel[0]->pos.y + sel[0]->size.y;
        for (auto* it : sel) max_b = std::max(max_b, it->pos.y + it->size.y);
        for (auto* it : sel) {
            it->pos.y = max_b - it->size.y;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "center_x") == 0) {
        float c = 0.0f;
        for (auto* it : sel) c += it->pos.x + it->size.x * 0.5f;
        c /= (float)sel.size();
        for (auto* it : sel) {
            it->pos.x = c - it->size.x * 0.5f;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "center_y") == 0) {
        float c = 0.0f;
        for (auto* it : sel) c += it->pos.y + it->size.y * 0.5f;
        c /= (float)sel.size();
        for (auto* it : sel) {
            it->pos.y = c - it->size.y * 0.5f;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
        }
    } else if (strcmp(mode, "distribute_x") == 0 && sel.size() >= 2) {
        std::sort(sel.begin(), sel.end(), [](const ComposerItem* a, const ComposerItem* b) {
            return a->pos.x < b->pos.x;
        });
        float start = sel.front()->pos.x;
        float end = sel.back()->pos.x + sel.back()->size.x;
        float total_w = 0.0f;
        for (auto* it : sel) total_w += it->size.x;
        float gap = (sel.size() > 1) ? (end - start - total_w) / (float)(sel.size() - 1) : 0.0f;
        float x = start;
        for (auto* it : sel) {
            it->pos.x = x;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
            x += it->size.x + gap;
        }
    } else if (strcmp(mode, "distribute_y") == 0 && sel.size() >= 2) {
        std::sort(sel.begin(), sel.end(), [](const ComposerItem* a, const ComposerItem* b) {
            return a->pos.y < b->pos.y;
        });
        float start = sel.front()->pos.y;
        float end = sel.back()->pos.y + sel.back()->size.y;
        float total_h = 0.0f;
        for (auto* it : sel) total_h += it->size.y;
        float gap = (sel.size() > 1) ? (end - start - total_h) / (float)(sel.size() - 1) : 0.0f;
        float y = start;
        for (auto* it : sel) {
            it->pos.y = y;
            it->pos = snap_pos(it->pos);
            clamp_to_page(it);
            y += it->size.y + gap;
        }
    }
}

static void ensure_composer_initialized(AppState* app) {
    if (!app) return;
    if (!app->composer_pages.empty()) return;
    app->composer_next_item_id = 1;
    ComposerPage page{};
    std::snprintf(page.name, sizeof(page.name), "Page 1");
    page.res_w = 1280;
    page.res_h = 720;
    page.bg = ImVec4(0.08f, 0.08f, 0.10f, 1.0f);
    page.transparent_bg = false;
    std::snprintf(page.output_name, sizeof(page.output_name), "composer_page_1");
    page.output_format = 0;
    page.fps = 30;
    page.frames = 60;
     page.video_kbps = 4000;
    composer_add_item(app, page, COMPOSER_FIELD_VIEW, 0, ImVec2(80, 80), ImVec2(720, 480));
    composer_add_item(app, page, COMPOSER_LEGEND, 0, ImVec2(840, 80), ImVec2(320, 200));
    app->composer_pages.push_back(page);
    app->composer_active_page = 0;
    std::snprintf(app->composer_status, sizeof(app->composer_status), "Composer ready");
    app->composer_last_export_path[0] = '\0';
    app->composer_preview_tex = nullptr;
    app->composer_preview_tex_size = ImVec2(0, 0);
    app->composer_preview_frames.clear();
    app->composer_preview_frame_idx = 0;
    app->composer_preview_frame_count = 0;
    app->composer_preview_fps = 30.0f;
    app->composer_preview_time = 0.0f;
    app->composer_preview_playing = false;
    app->composer_preview_max_frames = 120;
    app->composer_preview_target_width = 480;
    app->composer_preview_mutate_sim = true;
    app->composer_user_templates.clear();
    app->composer_dragging = false;
    app->composer_drag_item = -1;
    app->composer_drag_offset = ImVec2(0, 0);
    app->composer_resizing = false;
    app->composer_resize_item = -1;
    app->composer_resize_corner = -1;
    app->composer_resize_start_pos = ImVec2(0, 0);
    app->composer_resize_start_size = ImVec2(0, 0);
    app->composer_resize_start_mouse = ImVec2(0, 0);
    app->composer_snap = true;
    app->composer_snap_step = 4.0f;
    app->composer_canvas_grid = false;
    app->composer_canvas_grid_step = 32;
    app->composer_request_export_all = false;
    app->composer_pick_region_active = false;
    app->composer_pick_region_dragging = false;
    app->composer_pick_region_viewport = 0;
    app->composer_pick_region_page = 0;
    app->composer_pick_region_target_item = -1;
    app->composer_pick_region_start_norm = ImVec2(0.0f, 0.0f);
    app->composer_pick_region_end_norm = ImVec2(1.0f, 1.0f);
    app->composer_pick_region_snap = true;
    app->composer_request_animation = false;
    app->composer_animation_steps_per_frame = 1;
}

static void draw_measurement_history_panel(AppState* app) {
    if (!app) return;
    if (!ImGui::Begin("Measurement History", &app->show_measurement_history)) {
        ImGui::End();
        return;
    }

    ImGui::Checkbox("Show Distances", &app->show_distance_measurements);
        ImGui::Checkbox("Show Areas", &app->show_area_measurements);
        ImGui::Checkbox("Show Annotations", &app->show_annotations);
        ImGui::Separator();

        if (ImGui::TreeNode("Distances")) {
        for (size_t i = 0; i < app->measurements.distances.size(); ++i) {
            auto& d = app->measurements.distances[i];
            ImGui::PushID((int)i);
            ImGui::Text("%zu: %.6f m @ %.1f deg", i + 1, d.distance_m, d.angle_deg);
            ImGui::SameLine();
            if (ImGui::SmallButton("Delete")) {
                app->measurements.distances.erase(app->measurements.distances.begin() + i);
                ImGui::PopID();
                break;
            }
            ImGui::PopID();
        }
        ImGui::TreePop();
    }

        if (ImGui::TreeNode("Areas")) {
            for (size_t i = 0; i < app->measurements.areas.size(); ++i) {
                auto& a = app->measurements.areas[i];
                ImGui::PushID((int)i);
                ImGui::Text("%zu: Area=%.6f m^2, Perimeter=%.4f m", i + 1, a.area_m2, a.perimeter_m);
            ImGui::SameLine();
            if (ImGui::SmallButton("Delete")) {
                app->measurements.areas.erase(app->measurements.areas.begin() + i);
                ImGui::PopID();
                break;
            }
            ImGui::PopID();
            }
            ImGui::TreePop();
        }
        if (app->measurements.areas.size() > 0) {
            if (ImGui::Button("Clear Areas")) {
                app->measurements.areas.clear();
                ui_log_add(app, "Areas cleared");
            }
        }

        if (ImGui::TreeNode("Annotations")) {
            for (size_t i = 0; i < app->measurements.annotations.size(); ++i) {
                auto& ann = app->measurements.annotations[i];
                ImGui::PushID((int)i);
            ImGui::Checkbox("Visible", &ann.visible);
            ImGui::SameLine();
            ImGui::Text("%zu: \"%s\"", i + 1, ann.text);
            ImGui::SameLine();
            if (ImGui::SmallButton("Delete")) {
                app->measurements.annotations.erase(app->measurements.annotations.begin() + i);
                ImGui::PopID();
                break;
            }
            ImGui::PopID();
        }
        ImGui::TreePop();
    }

    if (ImGui::Button("Export CSV", ImVec2(-1, 0))) {
        export_measurements_csv(app, "measurements.csv");
        ui_log_add(app, "Measurements exported: measurements.csv");
    }
    if (ImGui::Button("Clear All", ImVec2(-1, 0))) {
        app->measurements.distances.clear();
        app->measurements.areas.clear();
        app->measurements.annotations.clear();
        ui_log_add(app, "Measurements cleared");
    }

    ImGui::End();
}

static ImU32 composer_item_color(const ComposerItem& item) {
    switch (item.type) {
        case COMPOSER_FIELD_VIEW: return IM_COL32(70, 150, 240, 200);
        case COMPOSER_REGION:     return IM_COL32(160, 120, 240, 200);
        case COMPOSER_SCOPE:      return IM_COL32(120, 200, 120, 200);
        case COMPOSER_FFT:        return IM_COL32(200, 160, 80, 200);
        case COMPOSER_LEGEND:     return IM_COL32(220, 220, 120, 200);
        case COMPOSER_MEAS:       return IM_COL32(240, 120, 120, 200);
        case COMPOSER_SMITH:      return IM_COL32(120, 200, 200, 200);
        default:                  return IM_COL32(200, 200, 200, 180);
    }
}

static SDL_Color composer_viewport_clear(const AppState* app) {
    int preset = app ? app->theme_preset : 0;
    switch (preset) {
        case 2:  // THEME_PRESET_LIGHT
            return SDL_Color{235, 236, 240, 255};
        case 3:  // THEME_PRESET_HIGH_CONTRAST
            return SDL_Color{4, 4, 4, 255};
        default:  // Dark/Blender
            return SDL_Color{8, 8, 10, 255};
    }
}

static void composer_theme(SDL_Color* bg, SDL_Color* border, SDL_Color* accent) {
    if (bg) *bg = SDL_Color{22, 24, 34, 255};
    if (border) *border = SDL_Color{70, 74, 96, 255};
    if (accent) *accent = SDL_Color{0, 220, 120, 255};
}

static SDL_Surface* composer_render_field_surface(AppState* app,
                                                  RenderContext* render,
                                                  const SimulationState* sim,
                                                  const ViewportInstance& vp_src,
                                                  int target_w,
                                                  int target_h,
                                                  bool use_region,
                                                  ImVec4 region_norm) {
    if (!render || !render->renderer || !sim) return nullptr;
    if (target_w <= 1 || target_h <= 1) return nullptr;

    SDL_Renderer* rr = render->renderer;
    SDL_Texture* tex = SDL_CreateTexture(rr,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         target_w,
                                         target_h);
    if (!tex) return nullptr;

    SDL_Texture* prev_target = SDL_GetRenderTarget(rr);
    SDL_Rect prev_viewport;
    SDL_RenderGetViewport(rr, &prev_viewport);
    int prev_scale = render->scale;
    float prev_off_x = render->offset_x;
    float prev_off_y = render->offset_y;
    SDL_BlendMode prev_blend = SDL_BLENDMODE_NONE;
    SDL_GetRenderDrawBlendMode(rr, &prev_blend);

    SDL_Color clear = composer_viewport_clear(app);
    SDL_SetRenderTarget(rr, tex);
    SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rr, clear.r, clear.g, clear.b, clear.a);
    SDL_RenderClear(rr);

    SDL_Rect full = {0, 0, target_w, target_h};
    SDL_RenderSetViewport(rr, &full);

    ViewportInstance vp = vp_src;
    vp.size = ImVec2((float)target_w, (float)target_h);
    int vp_scale = (int)std::lround(vp.zoom);
    if (vp_scale < 1) vp_scale = 1;
    ImVec2 vp_offset(0.0f, 0.0f);

    if (use_region) {
        float nx = (float)sim->nx;
        float ny = (float)sim->ny;
        float rx0 = std::clamp(region_norm.x, 0.0f, 1.0f);
        float ry0 = std::clamp(region_norm.y, 0.0f, 1.0f);
        float rx1 = std::clamp(region_norm.z, 0.0f, 1.0f);
        float ry1 = std::clamp(region_norm.w, 0.0f, 1.0f);
        if (rx1 < rx0 + 0.001f) rx1 = rx0 + 0.001f;
        if (ry1 < ry0 + 0.001f) ry1 = ry0 + 0.001f;
        float reg_w_cells = (rx1 - rx0) * nx;
        float reg_h_cells = (ry1 - ry0) * ny;
        float scale_x = (float)target_w / std::max(1.0f, reg_w_cells);
        float scale_y = (float)target_h / std::max(1.0f, reg_h_cells);
        float scale_f = std::min(scale_x, scale_y);
        vp_scale = (int)std::lround(scale_f);
        if (vp_scale < 1) vp_scale = 1;
        render->scale = vp_scale;
        float region_px_w = reg_w_cells * (float)vp_scale;
        float region_px_h = reg_h_cells * (float)vp_scale;
        float center_x = ((float)target_w - region_px_w) * 0.5f;
        float center_y = ((float)target_h - region_px_h) * 0.5f;
        vp_offset.x = center_x - rx0 * nx * (float)vp_scale;
        vp_offset.y = center_y - ry0 * ny * (float)vp_scale;
        render->offset_x = vp_offset.x;
        render->offset_y = vp_offset.y;
    } else {
        render->scale = vp_scale;
        vp_offset = compute_viewport_offset(vp, sim, vp_scale);
        render->offset_x = vp_offset.x;
        render->offset_y = vp_offset.y;
    }

    double vmax = sim->step_Ez_absmax;
    if (vmax <= 0.0) vmax = 1.0;
    bool channel_na = false;
    switch (vp.viz_mode) {
        case VIEWPORT_VIZ_EZ:
            render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_EZ_ABS:
            render_field_channel_heatmap(render, sim, FIELD_CH_EZ_ABS, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_HX:
            render_field_channel_heatmap(render, sim, FIELD_CH_HX, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_HY:
            render_field_channel_heatmap(render, sim, FIELD_CH_HY, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_HMAG:
            render_field_channel_heatmap(render, sim, FIELD_CH_H_MAG, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_SX:
            render_field_channel_heatmap(render, sim, FIELD_CH_SX, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_SY:
            render_field_channel_heatmap(render, sim, FIELD_CH_SY, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_S_MAG:
            render_field_channel_heatmap(render, sim, FIELD_CH_S_MAG, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_EX:
            render_field_channel_heatmap(render, sim, FIELD_CH_EX, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_EY:
            render_field_channel_heatmap(render, sim, FIELD_CH_EY, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_HZ:
            render_field_channel_heatmap(render, sim, FIELD_CH_HZ, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_MATERIAL:
            // Skip material-only view for exports to avoid block overlays bleeding
            render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
            break;
        case VIEWPORT_VIZ_OVERLAY:
            render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
            // Skip material overlay/outlines to keep exports clean
            break;
        default:
            render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
            break;
    }

    if (app && app->show_material_outlines) {
        SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);
        render_material_outlines(render, sim, vp_scale);
        SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_NONE);
    }
    // Suppress grid, sources, vectors for clean exports

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_SetRenderTarget(rr, prev_target);
        SDL_RenderSetViewport(rr, &prev_viewport);
        render->scale = prev_scale;
        render->offset_x = prev_off_x;
        render->offset_y = prev_off_y;
        SDL_SetRenderDrawBlendMode(rr, prev_blend);
        SDL_DestroyTexture(tex);
        return nullptr;
    }

    if (SDL_RenderReadPixels(rr, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        surf = nullptr;
    }

    SDL_SetRenderTarget(rr, prev_target);
    SDL_RenderSetViewport(rr, &prev_viewport);
    render->scale = prev_scale;
    render->offset_x = prev_off_x;
    render->offset_y = prev_off_y;
    SDL_SetRenderDrawBlendMode(rr, prev_blend);
    SDL_DestroyTexture(tex);
    return surf;
}

static bool export_composer_page_stub(AppState* app, SDL_Renderer* renderer, int page_idx) {
    if (!app || !renderer) return false;
    ensure_composer_initialized(app);
    if (page_idx < 0 || page_idx >= (int)app->composer_pages.size()) return false;
    const ComposerPage& page = app->composer_pages[page_idx];
    SDL_Texture* old_target = SDL_GetRenderTarget(renderer);
    SDL_Texture* tex = SDL_CreateTexture(renderer,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         page.res_w,
                                         page.res_h);
    if (!tex) return false;
    SDL_SetRenderTarget(renderer, tex);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(renderer,
                           (Uint8)(page.bg.x * 255),
                           (Uint8)(page.bg.y * 255),
                           (Uint8)(page.bg.z * 255),
                           (Uint8)(page.bg.w * 255));
    SDL_RenderClear(renderer);
    for (const auto& item : page.items) {
        SDL_Rect r;
        r.x = (int)std::lround(item.pos.x);
        r.y = (int)std::lround(item.pos.y);
        r.w = (int)std::lround(item.size.x);
        r.h = (int)std::lround(item.size.y);
        ImU32 col = composer_item_color(item);
        SDL_SetRenderDrawColor(renderer,
                               (Uint8)((col >> IM_COL32_R_SHIFT) & 0xFF),
                               (Uint8)((col >> IM_COL32_G_SHIFT) & 0xFF),
                               (Uint8)((col >> IM_COL32_B_SHIFT) & 0xFF),
                               220);
        SDL_RenderFillRect(renderer, &r);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawRect(renderer, &r);
    }
    SDL_SetRenderTarget(renderer, old_target);

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, page.res_w, page.res_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_DestroyTexture(tex);
        return false;
    }
    if (SDL_RenderReadPixels(renderer, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        SDL_DestroyTexture(tex);
        return false;
    }

#ifdef _WIN32
    _mkdir("recordings");
#else
    mkdir("recordings", 0755);
#endif
    char out_path[260];
    std::snprintf(out_path, sizeof(out_path), "recordings/composer_page_%d.png", page_idx + 1);
    SDL_SaveBMP(surf, out_path);
    SDL_FreeSurface(surf);
    SDL_DestroyTexture(tex);
    std::snprintf(app->composer_status, sizeof(app->composer_status), "Exported %s", out_path);
    ui_log_add(app, "Composer exported: %s", out_path);
    return true;
}

static SDL_Surface* capture_screen_surface(SDL_Renderer* renderer, int* out_w, int* out_h) {
    if (!renderer) return nullptr;
    int w = 0, h = 0;
    SDL_GetRendererOutputSize(renderer, &w, &h);
    if (out_w) *out_w = w;
    if (out_h) *out_h = h;
    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) return nullptr;
    if (SDL_RenderReadPixels(renderer, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        return nullptr;
    }
    return surf;
}

static double composer_scope_absmax(const Scope* scope) {
    if (!scope || !scope->y || scope->n <= 0) return 1.0;
    double m = 0.0;
    for (int i = 0; i < scope->n; ++i) {
        double v = std::fabs(scope->y[i]);
        if (v > m) m = v;
    }
    if (m <= 1e-9) m = 1.0;
    return m;
}

static SDL_Surface* composer_render_scope_surface(AppState* app,
                                                  RenderContext* render,
                                                  const Scope* scope,
                                                  int target_w,
                                                  int target_h) {
    if (!render || !render->renderer || !scope) return nullptr;
    if (target_w <= 1 || target_h <= 1) return nullptr;
    SDL_Renderer* rr = render->renderer;
    SDL_Texture* tex = SDL_CreateTexture(rr,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         target_w,
                                         target_h);
    if (!tex) return nullptr;

    SDL_Texture* prev_target = SDL_GetRenderTarget(rr);
    SDL_Rect prev_viewport;
    SDL_RenderGetViewport(rr, &prev_viewport);
    int prev_scale = render->scale;
    float prev_off_x = render->offset_x;
    float prev_off_y = render->offset_y;
    SDL_BlendMode prev_blend = SDL_BLENDMODE_NONE;
    SDL_GetRenderDrawBlendMode(rr, &prev_blend);

    SDL_Rect full = {0, 0, target_w, target_h};
    SDL_SetRenderTarget(rr, tex);
    SDL_RenderSetViewport(rr, &full);
    SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);
    double yscale = composer_scope_absmax(scope);
    render_scope(render, scope, 0, 0, target_w, target_h, yscale);

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_SetRenderTarget(rr, prev_target);
        SDL_RenderSetViewport(rr, &prev_viewport);
        render->scale = prev_scale;
        render->offset_x = prev_off_x;
        render->offset_y = prev_off_y;
        SDL_SetRenderDrawBlendMode(rr, prev_blend);
        SDL_DestroyTexture(tex);
        return nullptr;
    }
    if (SDL_RenderReadPixels(rr, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        surf = nullptr;
    }

    SDL_SetRenderTarget(rr, prev_target);
    SDL_RenderSetViewport(rr, &prev_viewport);
    render->scale = prev_scale;
    render->offset_x = prev_off_x;
    render->offset_y = prev_off_y;
    SDL_SetRenderDrawBlendMode(rr, prev_blend);
    SDL_DestroyTexture(tex);
    return surf;
}

static SDL_Surface* composer_render_fft_surface(AppState* app,
                                                RenderContext* render,
                                                const SimulationState* sim,
                                                const Scope* scope,
                                                int target_w,
                                                int target_h) {
    if (!render || !render->renderer || !sim || !scope) return nullptr;
    if (target_w <= 1 || target_h <= 1) return nullptr;

    static double fft_freq[1024];
    static double fft_mag[1024];
    static double fft_phase[1024];
    static double fft_mag_db[1024];
    int fft_points = compute_scope_fft(scope, sim->dt, fft_freq, fft_mag, fft_phase, IM_ARRAYSIZE(fft_freq));
    if (fft_points < 4) return nullptr;

    const int start_idx = 1;  // skip DC
    int plot_count = fft_points - start_idx;
    double max_db = -1e9;
    const double eps = 1e-9;
    for (int i = 0; i < fft_points; ++i) {
        fft_mag_db[i] = 20.0 * std::log10(fft_mag[i] + eps);
        if (i >= start_idx && fft_mag_db[i] > max_db) max_db = fft_mag_db[i];
    }
    if (max_db < -120.0) max_db = -120.0;
    double min_db = max_db - 80.0;  // 80 dB span

    SDL_Renderer* rr = render->renderer;
    SDL_Texture* tex = SDL_CreateTexture(rr,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         target_w,
                                         target_h);
    if (!tex) return nullptr;

    SDL_Texture* prev_target = SDL_GetRenderTarget(rr);
    SDL_Rect prev_viewport;
    SDL_RenderGetViewport(rr, &prev_viewport);
    SDL_BlendMode prev_blend = SDL_BLENDMODE_NONE;
    SDL_GetRenderDrawBlendMode(rr, &prev_blend);

    SDL_Rect full = {0, 0, target_w, target_h};
    SDL_SetRenderTarget(rr, tex);
    SDL_RenderSetViewport(rr, &full);
    SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);

    SDL_Color bg, border, accent;
    composer_theme(&bg, &border, &accent);
    SDL_SetRenderDrawColor(rr, bg.r, bg.g, bg.b, 255);
    SDL_RenderClear(rr);
    SDL_SetRenderDrawColor(rr, border.r, border.g, border.b, 255);
    SDL_RenderDrawRect(rr, &full);

    SDL_SetRenderDrawColor(rr, accent.r, accent.g, accent.b, accent.a);
    for (int k = start_idx + 1; k < fft_points; ++k) {
        double t0 = (double)(k - 1 - start_idx) / (double)(plot_count - 1);
        double t1 = (double)(k - start_idx) / (double)(plot_count - 1);
        double v0 = fft_mag_db[k - 1];
        double v1 = fft_mag_db[k];
        double y0n = (v0 - min_db) / (max_db - min_db);
        double y1n = (v1 - min_db) / (max_db - min_db);
        if (y0n < 0.0) y0n = 0.0;
        if (y0n > 1.0) y0n = 1.0;
        if (y1n < 0.0) y1n = 0.0;
        if (y1n > 1.0) y1n = 1.0;
        int x0 = (int)std::lround(t0 * (double)(target_w - 1));
        int x1 = (int)std::lround(t1 * (double)(target_w - 1));
        int y0 = target_h - 1 - (int)std::lround(y0n * (double)(target_h - 1));
        int y1 = target_h - 1 - (int)std::lround(y1n * (double)(target_h - 1));
        SDL_RenderDrawLine(rr, x0, y0, x1, y1);
    }

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_SetRenderTarget(rr, prev_target);
        SDL_RenderSetViewport(rr, &prev_viewport);
        SDL_SetRenderDrawBlendMode(rr, prev_blend);
        SDL_DestroyTexture(tex);
        return nullptr;
    }
    if (SDL_RenderReadPixels(rr, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        surf = nullptr;
    }

    SDL_SetRenderTarget(rr, prev_target);
    SDL_RenderSetViewport(rr, &prev_viewport);
    SDL_SetRenderDrawBlendMode(rr, prev_blend);
    SDL_DestroyTexture(tex);
    return surf;
}

static SDL_Surface* composer_render_legend_surface(RenderContext* render,
                                                   int target_w,
                                                   int target_h) {
    if (!render || !render->renderer) return nullptr;
    if (target_w <= 1 || target_h <= 1) return nullptr;
    SDL_Renderer* rr = render->renderer;
    SDL_Texture* tex = SDL_CreateTexture(rr,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         target_w,
                                         target_h);
    if (!tex) return nullptr;

    SDL_Texture* prev_target = SDL_GetRenderTarget(rr);
    SDL_Rect prev_viewport;
    SDL_RenderGetViewport(rr, &prev_viewport);
    SDL_BlendMode prev_blend = SDL_BLENDMODE_NONE;
    SDL_GetRenderDrawBlendMode(rr, &prev_blend);

    SDL_Rect full = {0, 0, target_w, target_h};
    SDL_SetRenderTarget(rr, tex);
    SDL_RenderSetViewport(rr, &full);
    SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);
    SDL_Color bg, border, accent;
    composer_theme(&bg, &border, &accent);
    SDL_SetRenderDrawColor(rr, bg.r, bg.g, bg.b, 255);
    SDL_RenderClear(rr);
    SDL_SetRenderDrawColor(rr, border.r, border.g, border.b, 255);
    SDL_RenderDrawRect(rr, &full);
    render_legend(render, 4, 4, target_w - 8);

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_SetRenderTarget(rr, prev_target);
        SDL_RenderSetViewport(rr, &prev_viewport);
        SDL_SetRenderDrawBlendMode(rr, prev_blend);
        SDL_DestroyTexture(tex);
        return nullptr;
    }
    if (SDL_RenderReadPixels(rr, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        surf = nullptr;
    }

    SDL_SetRenderTarget(rr, prev_target);
    SDL_RenderSetViewport(rr, &prev_viewport);
    SDL_SetRenderDrawBlendMode(rr, prev_blend);
    SDL_DestroyTexture(tex);
    return surf;
}

static SDL_Surface* composer_render_measurements_surface(RenderContext* render,
                                                         const AppState* app,
                                                         int target_w,
                                                         int target_h) {
    if (!render || !render->renderer || !app) return nullptr;
    if (target_w <= 1 || target_h <= 1) return nullptr;
    SDL_Renderer* rr = render->renderer;
    SDL_Texture* tex = SDL_CreateTexture(rr,
                                         SDL_PIXELFORMAT_RGBA8888,
                                         SDL_TEXTUREACCESS_TARGET,
                                         target_w,
                                         target_h);
    if (!tex) return nullptr;

    SDL_Texture* prev_target = SDL_GetRenderTarget(rr);
    SDL_Rect prev_viewport;
    SDL_RenderGetViewport(rr, &prev_viewport);
    SDL_BlendMode prev_blend = SDL_BLENDMODE_NONE;
    SDL_GetRenderDrawBlendMode(rr, &prev_blend);

    SDL_Rect full = {0, 0, target_w, target_h};
    SDL_SetRenderTarget(rr, tex);
    SDL_RenderSetViewport(rr, &full);
    SDL_SetRenderDrawBlendMode(rr, SDL_BLENDMODE_BLEND);
    SDL_Color bg, border, accent;
    composer_theme(&bg, &border, &accent);
    SDL_SetRenderDrawColor(rr, bg.r, bg.g, bg.b, 255);
    SDL_RenderClear(rr);
    SDL_SetRenderDrawColor(rr, border.r, border.g, border.b, 255);
    SDL_RenderDrawRect(rr, &full);

    char line0[128];
    char line1[128];
    char line2[128];
    std::snprintf(line0, sizeof(line0), "Distances: %zu", app->measurements.distances.size());
    std::snprintf(line1, sizeof(line1), "Areas: %zu", app->measurements.areas.size());
    std::snprintf(line2, sizeof(line2), "Annotations: %zu", app->measurements.annotations.size());

    SDL_Color text_col = border;
    int tw, th;
    int cursor_y = 8;
    const int gap = 6;

    auto blit_line = [&](const char* text) {
        SDL_Texture* t = render_text(render, text, text_col, &tw, &th);
        if (!t) return;
        SDL_Rect dst = {8, cursor_y, tw, th};
        SDL_RenderCopy(rr, t, NULL, &dst);
        SDL_DestroyTexture(t);
        cursor_y += th + gap;
    };

    blit_line(line0);
    blit_line(line1);
    blit_line(line2);

    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surf) {
        SDL_SetRenderTarget(rr, prev_target);
        SDL_RenderSetViewport(rr, &prev_viewport);
        SDL_SetRenderDrawBlendMode(rr, prev_blend);
        SDL_DestroyTexture(tex);
        return nullptr;
    }
    if (SDL_RenderReadPixels(rr, NULL, surf->format->format, surf->pixels, surf->pitch) != 0) {
        SDL_FreeSurface(surf);
        surf = nullptr;
    }

    SDL_SetRenderTarget(rr, prev_target);
    SDL_RenderSetViewport(rr, &prev_viewport);
    SDL_SetRenderDrawBlendMode(rr, prev_blend);
    SDL_DestroyTexture(tex);
    return surf;
}

static SDL_Surface* render_composer_page_surface(AppState* app,
                                                  RenderContext* render,
                                                  SimulationState* sim,
                                                  const Scope* scope,
                                                  int page_idx,
                                                  bool* out_any_field,
                                                  bool* out_any_field_blitted) {
    if (out_any_field) *out_any_field = false;
    if (out_any_field_blitted) *out_any_field_blitted = false;
    if (!app || !render || !render->renderer || !sim) return nullptr;
    ensure_composer_initialized(app);
    if (page_idx < 0 || page_idx >= (int)app->composer_pages.size()) return nullptr;
    const ComposerPage& page = app->composer_pages[page_idx];
    std::printf("Composer render page %d items=%zu bg=%.3f,%.3f,%.3f,%.3f transparent=%d\n",
                page_idx,
                page.items.size(),
                page.bg.x, page.bg.y, page.bg.z, page.bg.w,
                page.transparent_bg ? 1 : 0);

    SDL_Surface* out = SDL_CreateRGBSurfaceWithFormat(0, page.res_w, page.res_h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!out) return nullptr;
    Uint8 alpha = page.transparent_bg ? 0 : (Uint8)(page.bg.w * 255);
    std::printf("Composer render fill bg rgba=(%d,%d,%d,%d)\n",
                (int)(page.bg.x * 255),
                (int)(page.bg.y * 255),
                (int)(page.bg.z * 255),
                (int)alpha);
    SDL_FillRect(out, nullptr,
                 SDL_MapRGBA(out->format,
                             (Uint8)(page.bg.x * 255),
                             (Uint8)(page.bg.y * 255),
                             (Uint8)(page.bg.z * 255),
                             alpha));
    bool any_field = false;
    bool any_field_blitted = false;

    for (const auto& item : page.items) {
        SDL_Rect dst;
        dst.x = (int)std::lround(item.pos.x);
        dst.y = (int)std::lround(item.pos.y);
        dst.w = (int)std::lround(item.size.x);
        dst.h = (int)std::lround(item.size.y);
        std::printf("  item id=%d type=%d pos=(%d,%d) size=(%d,%d)\n",
                    item.id,
                    (int)item.type,
                    dst.x, dst.y, dst.w, dst.h);

        bool blitted = false;
        if (item.type == COMPOSER_FIELD_VIEW || item.type == COMPOSER_REGION) {
            if (item.viewport_idx >= 0 && item.viewport_idx < 4) {
                any_field = true;
                const ViewportInstance& vp = app->viewports[item.viewport_idx];
                if (vp.valid) {
                    bool use_region = (item.type == COMPOSER_REGION);
                    bool prev_outlines = app->show_material_outlines;
                    app->show_material_outlines = false;
                    ViewportInstance vp_clean = vp;
                    vp_clean.show_grid = false;
                    vp_clean.show_sources = false;
                    vp_clean.show_vectors = false;
                    SDL_Surface* src = composer_render_field_surface(app,
                                                                     render,
                                                                     sim,
                                                                     vp_clean,
                                                                     dst.w,
                                                                     dst.h,
                                                                     use_region,
                                                                     item.region_norm);
                    app->show_material_outlines = prev_outlines;
                    if (src) {
                        SDL_Rect dest = {dst.x, dst.y, src->w, src->h};
                        SDL_BlitSurface(src, nullptr, out, &dest);
                        SDL_FreeSurface(src);
                        blitted = true;
                        any_field_blitted = true;
                    }
                }
            }
        }

        if (!blitted && item.type == COMPOSER_SCOPE) {
            SDL_Surface* src = composer_render_scope_surface(app, render, scope, dst.w, dst.h);
            if (src) {
                SDL_Rect dest = {dst.x, dst.y, src->w, src->h};
                SDL_BlitSurface(src, nullptr, out, &dest);
                SDL_FreeSurface(src);
                blitted = true;
            }
        }

        if (!blitted && item.type == COMPOSER_FFT) {
            SDL_Surface* src = composer_render_fft_surface(app, render, sim, scope, dst.w, dst.h);
            if (src) {
                SDL_Rect dest = {dst.x, dst.y, src->w, src->h};
                SDL_BlitSurface(src, nullptr, out, &dest);
                SDL_FreeSurface(src);
                blitted = true;
            }
        }

        if (!blitted && item.type == COMPOSER_LEGEND) {
            SDL_Surface* src = composer_render_legend_surface(render, dst.w, dst.h);
            if (src) {
                SDL_Rect dest = {dst.x, dst.y, src->w, src->h};
                SDL_BlitSurface(src, nullptr, out, &dest);
                SDL_FreeSurface(src);
                blitted = true;
            }
        }

        if (!blitted && item.type == COMPOSER_MEAS) {
            SDL_Surface* src = composer_render_measurements_surface(render, app, dst.w, dst.h);
            if (src) {
                SDL_Rect dest = {dst.x, dst.y, src->w, src->h};
                SDL_BlitSurface(src, nullptr, out, &dest);
                SDL_FreeSurface(src);
                blitted = true;
            }
        }

        if (!blitted) {
            Uint32 col = composer_item_color(item);
            SDL_Rect fill = dst;
            SDL_FillRect(out, &fill, SDL_MapRGBA(out->format,
                                                 (col >> IM_COL32_R_SHIFT) & 0xFF,
                                                 (col >> IM_COL32_G_SHIFT) & 0xFF,
                                                 (col >> IM_COL32_B_SHIFT) & 0xFF,
                                                 255));
            std::printf("    filled placeholder color=%08x\n", col);
        }
    }

    // Fallback: if the image is uniform, stamp placeholder rectangles so headless output is never blank.
    auto surface_has_variance = [](SDL_Surface* s) -> bool {
        if (!s) return false;
        Uint8* pixels = (Uint8*)s->pixels;
        int bpp = s->format->BytesPerPixel;
        if (bpp < 3) return false;
        Uint8 r0 = pixels[0], g0 = pixels[1], b0 = pixels[2];
        const int stride = s->pitch;
        for (int y = 0; y < s->h; y += std::max(1, s->h / 16)) {
            Uint8* row = pixels + y * stride;
            for (int x = 0; x < s->w; x += std::max(1, s->w / 16)) {
                Uint8* p = row + x * bpp;
                if (p[0] != r0 || p[1] != g0 || p[2] != b0) return true;
            }
        }
        return false;
    };

    if (!surface_has_variance(out)) {
        std::printf("Composer render: uniform surface, stamping %zu placeholders\n", page.items.size());
        for (const auto& item : page.items) {
            SDL_Rect dst;
            dst.x = (int)std::lround(item.pos.x);
            dst.y = (int)std::lround(item.pos.y);
            dst.w = (int)std::lround(item.size.x);
            dst.h = (int)std::lround(item.size.y);
            Uint32 col = composer_item_color(item);
            std::printf("    fallback fill id=%d type=%d color=%08x\n", item.id, (int)item.type, col);
            SDL_FillRect(out, &dst, SDL_MapRGBA(out->format,
                                                (col >> IM_COL32_R_SHIFT) & 0xFF,
                                                (col >> IM_COL32_G_SHIFT) & 0xFF,
                                                (col >> IM_COL32_B_SHIFT) & 0xFF,
                                                255));
            SDL_Rect border = dst;
            SDL_Rect lines[4] = {
                {border.x, border.y, border.w, 1},
                {border.x, border.y + border.h - 1, border.w, 1},
                {border.x, border.y, 1, border.h},
                {border.x + border.w - 1, border.y, 1, border.h}
            };
            Uint32 edge = SDL_MapRGBA(out->format, 255, 255, 255, 255);
            for (auto& l : lines) {
                SDL_FillRect(out, &l, edge);
            }
        }
        if (out_any_field_blitted) *out_any_field_blitted = false;
    }

    if (out_any_field) *out_any_field = any_field;
    if (out_any_field_blitted) *out_any_field_blitted = any_field_blitted;
    return out;
}

static bool export_composer_page(AppState* app,
                                 RenderContext* render,
                                 SimulationState* sim,
                                 const Scope* scope,
                                 int page_idx) {
    bool any_field = false;
    bool any_field_blitted = false;
    SDL_Surface* out = render_composer_page_surface(app, render, sim, scope, page_idx, &any_field, &any_field_blitted);
    if (!out) return false;
    const ComposerPage& page = app->composer_pages[page_idx];

#ifdef _WIN32
    _mkdir("recordings");
#else
    mkdir("recordings", 0755);
#endif
    char out_path[260];
    const char* base = page.output_name[0] ? page.output_name : "composer_page";
    const char* ext = "bmp";
    bool saved = false;
    std::string ffmpeg_cmd = ffmpeg_command();
    bool ffmpeg_ok = !ffmpeg_cmd.empty();

    if (page.output_format == 0) {
        ext = "bmp";
        std::snprintf(out_path, sizeof(out_path), "recordings/%s.%s", base, ext);
        saved = (SDL_SaveBMP(out, out_path) == 0);
    } else if (page.output_format == 1) {
        ext = "png";
        std::snprintf(out_path, sizeof(out_path), "recordings/%s.%s", base, ext);
        SDL_Surface* conv = SDL_ConvertSurfaceFormat(out, SDL_PIXELFORMAT_RGBA32, 0);
        if (conv) {
            saved = stbi_write_png(out_path,
                                   conv->w,
                                   conv->h,
                                   4,
                                   conv->pixels,
                                   conv->pitch) == 1;
            SDL_FreeSurface(conv);
        } else {
            saved = false;
        }
    } else if (page.output_format == 2 || page.output_format == 3) {
        // MP4 or GIF via ffmpeg; if ffmpeg missing/fails, fall back to PNG
        ext = (page.output_format == 2) ? "mp4" : "gif";
        std::snprintf(out_path, sizeof(out_path), "recordings/%s.%s", base, ext);
        if (ffmpeg_ok) {
            // Write a single-frame raw stream to ffmpeg
            std::filesystem::create_directories("recordings");
            std::vector<std::string> args;
            args.push_back(ffmpeg_cmd);
            args.push_back("-y");
            args.push_back("-f");
            args.push_back("rawvideo");
            args.push_back("-pix_fmt");
            args.push_back("rgba");
            args.push_back("-s");
            char sizebuf[32];
            std::snprintf(sizebuf, sizeof(sizebuf), "%dx%d", page.res_w, page.res_h);
            args.push_back(sizebuf);
            args.push_back("-i");
            args.push_back("-");
            if (page.output_format == 2) {
                int kbps = page.video_kbps > 0 ? page.video_kbps : 4000;
                kbps = std::clamp(kbps, 1000, 50000);
                args.push_back("-c:v");
                args.push_back("mpeg4");
                args.push_back("-b:v");
                args.push_back(std::to_string(kbps * 1000));
                args.push_back("-maxrate");
                args.push_back(std::to_string(kbps * 1000));
                args.push_back("-bufsize");
                args.push_back(std::to_string(kbps * 2000));
                args.push_back("-pix_fmt");
                args.push_back("yuv420p");
            } else {
                args.push_back("-filter_complex");
                args.push_back("fps=15,scale=iw:ih:flags=lanczos");
            }
            args.push_back(out_path);

            // Launch ffmpeg with stdin pipe
#ifdef _WIN32
            // Use _popen for simplicity
            std::string cmdline = "\"" + ffmpeg_cmd + "\"";
            for (size_t i = 1; i < args.size(); ++i) {
                cmdline += " " + args[i];
            }
            FILE* pipe = _popen(cmdline.c_str(), "wb");
#else
            std::string cmdline = ffmpeg_cmd;
            for (size_t i = 1; i < args.size(); ++i) {
                cmdline += " " + args[i];
            }
            FILE* pipe = popen(cmdline.c_str(), "w");
#endif
            if (pipe) {
                SDL_Surface* conv = SDL_ConvertSurfaceFormat(out, SDL_PIXELFORMAT_RGBA32, 0);
                if (conv) {
                    size_t bytes = conv->pitch * conv->h;
                    size_t written = fwrite(conv->pixels, 1, bytes, pipe);
                    saved = (written == bytes);
                    SDL_FreeSurface(conv);
                }
#ifdef _WIN32
                _pclose(pipe);
#else
                pclose(pipe);
#endif
            } else {
                saved = false;
            }
            log_ffmpeg_attempt(cmdline.c_str(), saved ? 0 : -1);
        }
        if (!saved) {
            // fallback to PNG
            std::snprintf(out_path, sizeof(out_path), "recordings/%s.png", base);
            SDL_Surface* conv = SDL_ConvertSurfaceFormat(out, SDL_PIXELFORMAT_RGBA32, 0);
            if (conv) {
                saved = stbi_write_png(out_path,
                                       conv->w,
                                       conv->h,
                                       4,
                                       conv->pixels,
                                       conv->pitch) == 1;
                SDL_FreeSurface(conv);
            }
            if (ffmpeg_ok) {
                ui_log_add(app, "FFmpeg export failed, saved PNG instead: %s", out_path);
            } else {
                ui_log_add(app, "FFmpeg not found, saved PNG instead: %s", out_path);
            }
        }
    }
    SDL_FreeSurface(out);
    if (saved) {
        if (any_field && !any_field_blitted) {
            std::snprintf(app->composer_status, sizeof(app->composer_status), "Exported %s (placeholders used)", out_path);
            ui_log_add(app, "Composer exported with placeholders (viewport not visible): %s", out_path);
        } else {
            std::snprintf(app->composer_status, sizeof(app->composer_status), "Exported %s", out_path);
            ui_log_add(app, "Composer exported: %s", out_path);
        }
        std::filesystem::path p(out_path);
        std::snprintf(app->composer_last_export_path, sizeof(app->composer_last_export_path), "%s", p.parent_path().string().c_str());
    } else {
        std::snprintf(app->composer_status, sizeof(app->composer_status), "Export failed: %s", out_path);
        ui_log_add(app, "Composer export failed: %s", out_path);
    }
    return true;
}

static bool export_composer_animation(AppState* app,
                                      RenderContext* render,
                                      SimulationState* sim,
                                      const Scope* scope,
                                      int page_idx,
                                      bool preview_only) {
    if (!app || !render || !render->renderer || !sim) return false;
    ensure_composer_initialized(app);
    if (page_idx < 0 || page_idx >= (int)app->composer_pages.size()) return false;
    const ComposerPage& page = app->composer_pages[page_idx];
    int frames = std::max(1, page.frames);
    int fps = std::clamp(page.fps, 5, 60);
    bool any_field = false;
    bool any_field_blitted = false;
    int steps_per_frame = std::max(1, app->composer_animation_steps_per_frame);
    // We'll render each frame after stepping sim (mutating sim in-place)
    SDL_Surface* surface = nullptr;

#ifdef _WIN32
    _mkdir("recordings");
#else
    mkdir("recordings", 0755);
#endif
    const char* base = page.output_name[0] ? page.output_name : "composer_page";
    bool saved = false;
    char out_path[260];
    std::string ffmpeg_cmd = ffmpeg_command();
    bool ffmpeg_ok = !ffmpeg_cmd.empty();
    std::error_code ec;

    bool fallback_png = false;
    int effective_format = page.output_format;
    if (!preview_only && (page.output_format == 2 || page.output_format == 3) && !ffmpeg_ok) {
        // Fallback to PNG sequence if ffmpeg unavailable
        effective_format = 1;
        fallback_png = true;
    }

    if (effective_format == 0 || effective_format == 1) {
        // PNG sequence
        std::filesystem::path dir = std::filesystem::path("recordings") / (std::string(base) + "_frames");
        std::filesystem::remove_all(dir, ec);
        std::filesystem::create_directories(dir, ec);
        saved = true;
        for (int i = 0; i < frames; ++i) {
            for (int s = 0; s < steps_per_frame; ++s) {
                fdtd_step(sim);
            }
            SDL_Surface* frame_surface = render_composer_page_surface(app, render, sim, scope, page_idx, &any_field, &any_field_blitted);
            if (!frame_surface) { saved = false; break; }
            SDL_Surface* conv = SDL_ConvertSurfaceFormat(frame_surface, SDL_PIXELFORMAT_RGBA32, 0);
            SDL_FreeSurface(frame_surface);
            if (!conv) { saved = false; break; }
            char frame_path[512];
            std::snprintf(frame_path, sizeof(frame_path), "%s/frame_%04d.png", dir.string().c_str(), i);
            int ok = stbi_write_png(frame_path,
                                    conv->w,
                                    conv->h,
                                    4,
                                    conv->pixels,
                                    conv->pitch);
            SDL_FreeSurface(conv);
            if (ok != 1) {
                saved = false;
                break;
            }
        }
        std::snprintf(out_path, sizeof(out_path), "%s", dir.string().c_str());
    } else {
        // MP4 or GIF via ffmpeg rawvideo pipe (no PNG decode)
        const char* ext = (effective_format == 2) ? "mp4" : "gif";
        std::filesystem::path video_path = std::filesystem::path("recordings") / (std::string(base) + "." + ext);
        std::snprintf(out_path, sizeof(out_path), "%s", video_path.string().c_str());
        saved = false;

        if (ffmpeg_ok && !preview_only) {
            std::vector<std::string> args;
            args.push_back(ffmpeg_cmd);
            args.push_back("-y");
            args.push_back("-f");
            args.push_back("rawvideo");
            args.push_back("-pix_fmt");
            args.push_back("rgba");
            args.push_back("-s");
            char sizebuf[32];
            std::snprintf(sizebuf, sizeof(sizebuf), "%dx%d", page.res_w, page.res_h);
            args.push_back(sizebuf);
            args.push_back("-r");
            args.push_back(std::to_string(fps));
            args.push_back("-i");
            args.push_back("-");
            if (effective_format == 2) {
                int kbps = page.video_kbps > 0 ? page.video_kbps : 4000;
                kbps = std::clamp(kbps, 1000, 50000);
                args.push_back("-c:v");
                args.push_back("mpeg4");
                args.push_back("-b:v");
                args.push_back(std::to_string(kbps * 1000));
                args.push_back("-maxrate");
                args.push_back(std::to_string(kbps * 1000));
                args.push_back("-bufsize");
                args.push_back(std::to_string(kbps * 2000));
                args.push_back("-pix_fmt");
                args.push_back("yuv420p");
            } else {
                args.push_back("-vf");
                args.push_back("fps=" + std::to_string(fps) + ",scale=iw:ih:flags=lanczos");
            }
            args.push_back(video_path.string());

#ifdef _WIN32
            std::string cmdline = "\"" + ffmpeg_cmd + "\"";
            for (size_t i = 1; i < args.size(); ++i) {
                cmdline += " " + args[i];
            }
            FILE* pipe = _popen(cmdline.c_str(), "wb");
#else
            std::string cmdline = join_args(args);
            FILE* pipe = popen(cmdline.c_str(), "w");
#endif
            if (pipe) {
                bool ok = true;
                for (int i = 0; i < frames && ok; ++i) {
                    if (i > 0) {
                        for (int s = 0; s < steps_per_frame; ++s) {
                            fdtd_step(sim);
                        }
                    }
                    SDL_Surface* frame_surface = render_composer_page_surface(app, render, sim, scope, page_idx, &any_field, &any_field_blitted);
                    SDL_Surface* conv = frame_surface ? SDL_ConvertSurfaceFormat(frame_surface, SDL_PIXELFORMAT_RGBA32, 0) : nullptr;
                    if (frame_surface) SDL_FreeSurface(frame_surface);
                    if (!conv) {
                        ok = false;
                        break;
                    }
                    size_t bytes = (size_t)conv->pitch * (size_t)conv->h;
                    size_t written = fwrite(conv->pixels, 1, bytes, pipe);
                    SDL_FreeSurface(conv);
                    if (written != bytes) {
                        ok = false;
                        break;
                    }
                }
#ifdef _WIN32
                int rc = _pclose(pipe);
#else
                int rc = pclose(pipe);
#endif
                saved = ok && (rc == 0);
                log_ffmpeg_attempt(cmdline.c_str(), rc);
            } else {
                log_ffmpeg_attempt(ffmpeg_cmd.c_str(), -1);
                saved = false;
            }
        }
    }

    if (preview_only) {
        std::snprintf(app->composer_status, sizeof(app->composer_status), "Preview done (%d frames @ %d fps)", frames, fps);
        return saved;
    }
    if (saved) {
        std::snprintf(app->composer_status,
                      sizeof(app->composer_status),
                      fallback_png ? "Animation exported PNG sequence: %s" : "Animation exported: %s",
                      out_path);
        std::filesystem::path p(out_path);
        std::snprintf(app->composer_last_export_path, sizeof(app->composer_last_export_path), "%s", p.parent_path().string().c_str());
        ui_log_add(app,
                   fallback_png ? "Composer animation exported PNG sequence: %s (%d frames @ %d fps)"
                                : "Composer animation exported: %s (%d frames @ %d fps)",
                   out_path,
                   frames,
                   fps);
    } else {
        std::snprintf(app->composer_status, sizeof(app->composer_status), "Animation export failed");
        ui_log_add(app, "Composer animation export failed");
    }
    return saved;
}

static bool composer_generate_preview(AppState* app,
                                      RenderContext* render,
                                      SimulationState* sim,
                                      const Scope* scope,
                                      int page_idx) {
    if (!app || !render || !render->renderer || !sim) return false;
    SDL_Surface* frame_surface = render_composer_page_surface(app, render, sim, scope, page_idx, nullptr, nullptr);
    if (!frame_surface) return false;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(render->renderer, frame_surface);
    SDL_FreeSurface(frame_surface);
    if (!tex) return false;
    if (app->composer_preview_tex) {
        SDL_DestroyTexture(app->composer_preview_tex);
    }
    int w = 0, h = 0;
    SDL_QueryTexture(tex, nullptr, nullptr, &w, &h);
    app->composer_preview_tex = tex;
    app->composer_preview_tex_size = ImVec2((float)w, (float)h);
    std::snprintf(app->composer_status, sizeof(app->composer_status), "Preview updated (%dx%d)", w, h);
    return true;
}

static void draw_print_composer(AppState* app, SDL_Renderer* renderer, SimulationState* sim, Scope* scope, RenderContext* render_ctx) {
    if (!app || !app->show_print_composer) return;
    ensure_composer_initialized(app);
    if (!ImGui::Begin("Print Composer", &app->show_print_composer)) {
        ImGui::End();
        return;
    }

    if (app->composer_pages.empty()) {
        ImGui::TextUnformatted("No pages defined.");
        ImGui::End();
        return;
    }
    if (app->composer_active_page >= (int)app->composer_pages.size()) {
        app->composer_active_page = (int)app->composer_pages.size() - 1;
    }

    ImGui::TextUnformatted("Pages:");
    ImGui::SameLine();
    if (ImGui::SmallButton("+ Add Page")) {
        ComposerPage p{};
        int clone_idx = app->composer_active_page;
        if (clone_idx < 0 || clone_idx >= (int)app->composer_pages.size()) clone_idx = 0;
        const ComposerPage& cur = app->composer_pages[clone_idx];
        std::snprintf(p.name, sizeof(p.name), "Page %zu", app->composer_pages.size() + 1);
        p.res_w = cur.res_w;
        p.res_h = cur.res_h;
        p.bg = cur.bg;
        p.transparent_bg = cur.transparent_bg;
        std::snprintf(p.output_name, sizeof(p.output_name), "composer_page_%zu", app->composer_pages.size() + 1);
        p.output_format = cur.output_format;
        p.fps = cur.fps;
        p.frames = cur.frames;
        p.video_kbps = cur.video_kbps;
        app->composer_pages.push_back(p);
        app->composer_active_page = (int)app->composer_pages.size() - 1;
    }

    ImGui::SameLine();
    if (ImGui::SmallButton("Duplicate")) {
        if (!app->composer_pages.empty() && app->composer_active_page >= 0 && app->composer_active_page < (int)app->composer_pages.size()) {
            ComposerPage dup = app->composer_pages[app->composer_active_page];
            std::snprintf(dup.name, sizeof(dup.name), "Page %zu", app->composer_pages.size() + 1);
            std::snprintf(dup.output_name, sizeof(dup.output_name), "composer_page_%zu", app->composer_pages.size() + 1);
            app->composer_pages.push_back(dup);
            app->composer_active_page = (int)app->composer_pages.size() - 1;
        }
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Delete")) {
        if (app->composer_pages.size() > 1 &&
            app->composer_active_page >= 0 &&
            app->composer_active_page < (int)app->composer_pages.size()) {
            app->composer_pages.erase(app->composer_pages.begin() + app->composer_active_page);
            if (app->composer_active_page >= (int)app->composer_pages.size()) {
                app->composer_active_page = (int)app->composer_pages.size() - 1;
            }
        }
    }

    static int tmpl_choice = 0;
    const char* tmpl_labels[] = {"Blank", "Field+Legend", "Field+Scope", "Field+FFT+Legend"};
    ImGui::PushItemWidth(180.0f);
    ImGui::Combo("Template", &tmpl_choice, tmpl_labels, IM_ARRAYSIZE(tmpl_labels));
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::SmallButton("Apply Template")) {
        if (app->composer_active_page >= 0 && app->composer_active_page < (int)app->composer_pages.size()) {
            composer_apply_template(app, app->composer_pages[app->composer_active_page], tmpl_choice);
        }
    }
    static char user_tmpl_name[32] = "MyTemplate";
    ImGui::InputText("Template name", user_tmpl_name, IM_ARRAYSIZE(user_tmpl_name));
    ImGui::SameLine();
    if (ImGui::SmallButton("Save Current as Template")) {
        if (app->composer_active_page >= 0 && app->composer_active_page < (int)app->composer_pages.size()) {
            ComposerPage saved = app->composer_pages[app->composer_active_page];
            std::snprintf(saved.name, sizeof(saved.name), "%s", user_tmpl_name[0] ? user_tmpl_name : "UserTemplate");
            app->composer_user_templates.push_back(saved);
            if (app->composer_user_templates.size() > 8) {
                app->composer_user_templates.erase(app->composer_user_templates.begin());
            }
        }
    }
    if (!app->composer_user_templates.empty()) {
        ImGui::SameLine();
        if (ImGui::BeginCombo("User Templates", app->composer_user_templates.front().name)) {
            for (size_t i = 0; i < app->composer_user_templates.size(); ++i) {
                bool selected = false;
                if (ImGui::Selectable(app->composer_user_templates[i].name, &selected)) {
                    if (app->composer_active_page >= 0 && app->composer_active_page < (int)app->composer_pages.size()) {
                        app->composer_pages[app->composer_active_page] = app->composer_user_templates[i];
                    }
                }
            }
            ImGui::EndCombo();
        }
    }

    ImGui::TextUnformatted("Pages:");
    ImGui::BeginChild("PageStrip", ImVec2(0, 110), true, ImGuiWindowFlags_HorizontalScrollbar);
    for (size_t i = 0; i < app->composer_pages.size(); ++i) {
        ImGui::PushID((int)i);
        ImVec2 thumb = ImVec2(140, 80);
        bool sel = (int)i == app->composer_active_page;
        ImGui::PushStyleColor(ImGuiCol_Button, sel ? ImGui::GetStyle().Colors[ImGuiCol_ButtonActive] : ImGui::GetStyle().Colors[ImGuiCol_Button]);
        if (ImGui::Button(app->composer_pages[i].name, thumb)) {
            app->composer_active_page = (int)i;
        }
        ImGui::PopStyleColor();
        ImGui::SameLine();
        ImGui::PopID();
    }
    ImGui::NewLine();
    ImGui::EndChild();

    if (app->composer_active_page < 0) app->composer_active_page = 0;
    ComposerPage& page = app->composer_pages[app->composer_active_page];
    if (app->composer_drag_item >= (int)page.items.size()) {
        app->composer_drag_item = -1;
        app->composer_dragging = false;
    }
    if (app->composer_resize_item >= (int)page.items.size()) {
        app->composer_resize_item = -1;
        app->composer_resizing = false;
        app->composer_resize_corner = -1;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Sources:");
    if (ImGui::Button("Add Field A")) composer_add_item(app, page, COMPOSER_FIELD_VIEW, 0, ImVec2(80, 80), ImVec2(400, 260));
    ImGui::SameLine();
    if (ImGui::Button("Add Field B")) composer_add_item(app, page, COMPOSER_FIELD_VIEW, 1, ImVec2(520, 80), ImVec2(400, 260));
    ImGui::SameLine();
    if (ImGui::Button("Add Legend")) composer_add_item(app, page, COMPOSER_LEGEND, 0, ImVec2(80, 360), ImVec2(280, 180));
    ImGui::SameLine();
    if (ImGui::Button("Add Scope")) composer_add_item(app, page, COMPOSER_SCOPE, 0, ImVec2(380, 360), ImVec2(420, 220));
    if (ImGui::Button("Add FFT")) composer_add_item(app, page, COMPOSER_FFT, 0, ImVec2(820, 360), ImVec2(420, 220));
    ImGui::SameLine();
    if (ImGui::Button("Add Measurements")) composer_add_item(app, page, COMPOSER_MEAS, 0, ImVec2(840, 80), ImVec2(300, 200));
    ImGui::SameLine();
    if (ImGui::Button("Pick Region")) {
        app->composer_pick_region_active = true;
        app->composer_pick_region_dragging = false;
        app->composer_pick_region_viewport = std::clamp(app->composer_pick_region_viewport, 0, 3);
        app->composer_pick_region_page = app->composer_active_page;
        app->composer_pick_region_target_item = -1;
        std::snprintf(app->composer_status, sizeof(app->composer_status), "Pick a region on viewport %d", app->composer_pick_region_viewport + 1);
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::SliderInt("Viewport##pick", &app->composer_pick_region_viewport, 0, 3, "Viewport %d");
    ImGui::SameLine();
    ImGui::Checkbox("Snap cells", &app->composer_pick_region_snap);
    ImGui::SameLine();
    if (ImGui::Button("Add Region")) composer_add_item(app, page, COMPOSER_REGION, 0, ImVec2(120, 140), ImVec2(360, 240));

    ImGui::Separator();
    ImGui::Text("Page %s", page.name);
    ImGui::InputInt("Resolution W", &page.res_w);
    ImGui::InputInt("Resolution H", &page.res_h);
    page.res_w = std::max(page.res_w, 320);
    page.res_h = std::max(page.res_h, 240);
    ImGui::InputText("Output name", page.output_name, IM_ARRAYSIZE(page.output_name));
    const char* fmt_labels[] = {"BMP", "PNG", "MP4", "GIF"};
    ImGui::Combo("Format", &page.output_format, fmt_labels, IM_ARRAYSIZE(fmt_labels));
    bool ffmpeg_ok = has_ffmpeg();
    ImGui::TextDisabled("FFmpeg: %s", ffmpeg_ok ? "found" : "missing (MP4/GIF will fall back to PNG)");
    ImGui::SliderInt("FPS (for video)", &page.fps, 5, 60);
    ImGui::InputInt("Frames (animation)", &page.frames);
    if (page.frames < 1) page.frames = 1;
    if (page.output_format == 2) {
        ImGui::BeginDisabled(!ffmpeg_ok);
        ImGui::SliderInt("Video bitrate (kbps, MP4)", &page.video_kbps, 1000, 50000);
        ImGui::EndDisabled();
        if (page.video_kbps < 1000) page.video_kbps = 1000;
        if (page.video_kbps > 50000) page.video_kbps = 50000;
    }
    ImGui::Checkbox("Transparent BG (PNG/GIF)", &page.transparent_bg);
    ImGui::ColorEdit4("Background", (float*)&page.bg, ImGuiColorEditFlags_NoInputs);
    ImGui::SliderInt("Preview max frames", &app->composer_preview_max_frames, 10, 240);
    ImGui::SliderInt("Preview width (px)", &app->composer_preview_target_width, 240, 800);
    ImGui::Checkbox("Preview mutates sim", &app->composer_preview_mutate_sim);
    ImGui::Checkbox("Snap to grid", &app->composer_snap);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputFloat("Step", &app->composer_snap_step, 1.0f, 4.0f, "%.1f");
    if (app->composer_snap_step < 1.0f) app->composer_snap_step = 1.0f;
    ImGui::Checkbox("Canvas grid", &app->composer_canvas_grid);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("Grid px", &app->composer_canvas_grid_step);
    if (app->composer_canvas_grid_step < 2) app->composer_canvas_grid_step = 2;
    ImGui::Separator();
    ImGui::InputInt("Anim steps/frame", &app->composer_animation_steps_per_frame);
    if (app->composer_animation_steps_per_frame < 1) app->composer_animation_steps_per_frame = 1;
    ImGui::TextDisabled("Frames: %d @ %d fps -> %.2f sec",
                        page.frames,
                        page.fps,
                        (float)page.frames / (float)page.fps);
    ImGui::Separator();
    int sel_count = 0;
    for (const auto& it : page.items) if (it.selected) sel_count++;
    ImGui::Text("Selected: %d", sel_count);
    ImGui::SameLine();
    if (ImGui::SmallButton("Select All")) composer_select_all(page, true);
    ImGui::SameLine();
    if (ImGui::SmallButton("Clear Sel")) composer_select_all(page, false);
    bool can_align = sel_count >= 1;
    bool can_distribute = sel_count >= 2;
    ImGui::BeginDisabled(!can_align);
    if (ImGui::Button("Align Left")) composer_align_selected(page, page.res_w, page.res_h, "left", app->composer_snap, app->composer_snap_step);
    ImGui::SameLine();
    if (ImGui::Button("Align Right")) composer_align_selected(page, page.res_w, page.res_h, "right", app->composer_snap, app->composer_snap_step);
    ImGui::SameLine();
    if (ImGui::Button("Align Top")) composer_align_selected(page, page.res_w, page.res_h, "top", app->composer_snap, app->composer_snap_step);
    ImGui::SameLine();
    if (ImGui::Button("Align Bottom")) composer_align_selected(page, page.res_w, page.res_h, "bottom", app->composer_snap, app->composer_snap_step);
    if (ImGui::Button("Align Center X")) composer_align_selected(page, page.res_w, page.res_h, "center_x", app->composer_snap, app->composer_snap_step);
    ImGui::SameLine();
    if (ImGui::Button("Align Center Y")) composer_align_selected(page, page.res_w, page.res_h, "center_y", app->composer_snap, app->composer_snap_step);
    ImGui::EndDisabled();
    ImGui::BeginDisabled(!can_distribute);
    if (ImGui::Button("Distribute X")) composer_align_selected(page, page.res_w, page.res_h, "distribute_x", app->composer_snap, app->composer_snap_step);
    ImGui::SameLine();
    if (ImGui::Button("Distribute Y")) composer_align_selected(page, page.res_w, page.res_h, "distribute_y", app->composer_snap, app->composer_snap_step);
    ImGui::EndDisabled();
    if (ImGui::Button("Export Page (BMP)")) {
        app->composer_request_export = true;
        app->composer_request_page = app->composer_active_page;
        std::snprintf(app->composer_status,
                      sizeof(app->composer_status),
                      "Queued export for %s",
                      page.name);
    }
    ImGui::SameLine();
    if (ImGui::Button("Export All Pages (BMP)")) {
        app->composer_request_export_all = true;
        std::snprintf(app->composer_status,
                      sizeof(app->composer_status),
                      "Queued export for all pages (%zu)",
                      app->composer_pages.size());
    }
    ImGui::SameLine();
    if (ImGui::Button("Export Static PNG")) {
        int old_fmt = page.output_format;
        page.output_format = 1;  // PNG
        export_composer_page(app, render_ctx, sim, scope, app->composer_active_page);
        page.output_format = old_fmt;
    }
    ImGui::SameLine();
    if (ImGui::Button("Export All Static PNG")) {
        for (size_t pi = 0; pi < app->composer_pages.size(); ++pi) {
            int old_fmt = app->composer_pages[pi].output_format;
            app->composer_pages[pi].output_format = 1;
            export_composer_page(app, render_ctx, sim, scope, (int)pi);
            app->composer_pages[pi].output_format = old_fmt;
        }
        std::snprintf(app->composer_status, sizeof(app->composer_status), "Exported all pages as PNG");
    }
    ImGui::SameLine();
    if (ImGui::Button("Export Animation (ffmpeg/PNG)")) {
        app->composer_request_export_all = false;
        app->composer_request_animation = true;
        app->composer_request_page = app->composer_active_page;
        std::snprintf(app->composer_status,
                      sizeof(app->composer_status),
                      "Queued animation export for %s",
                      page.name);
    }
    ImGui::SameLine();
    if (ImGui::Button("Preview Animation")) {
        clear_composer_preview(app);
        SimulationState* preview_sim = nullptr;
        if (app->composer_preview_mutate_sim) {
            preview_sim = sim;
        } else {
            preview_sim = clone_simulation_state(sim);
            if (!preview_sim) {
                std::snprintf(app->composer_status, sizeof(app->composer_status), "Preview clone failed");
            }
        }
        if (!preview_sim) {
            // no sim to preview
        } else {
        // Build preview frames (small) without file IO
        const int max_preview_frames = std::min(page.frames, app->composer_preview_max_frames);
        int steps_per_frame = std::max(1, app->composer_animation_steps_per_frame);
        float scale_factor = 1.0f;
        int target_w = page.res_w;
        int target_h = page.res_h;
        int preview_max_w = std::max(240, app->composer_preview_target_width);
        if (target_w > preview_max_w) {
            scale_factor = (float)preview_max_w / (float)target_w;
            target_w = preview_max_w;
            target_h = (int)std::lround(page.res_h * scale_factor);
        }
        bool ok_all = true;
        Uint64 t0 = SDL_GetPerformanceCounter();
        for (int i = 0; i < max_preview_frames; ++i) {
            if (i > 0) {
                for (int s = 0; s < steps_per_frame; ++s) {
                    fdtd_step(preview_sim);
                }
            }
            SDL_Surface* frame_surface = render_composer_page_surface(app, render_ctx, preview_sim, scope, app->composer_active_page, nullptr, nullptr);
            if (!frame_surface) { ok_all = false; break; }
            SDL_Surface* scaled = SDL_CreateRGBSurfaceWithFormat(0, target_w, target_h, 32, SDL_PIXELFORMAT_RGBA32);
            if (scaled) {
                SDL_Rect dst = {0, 0, target_w, target_h};
                SDL_BlitScaled(frame_surface, nullptr, scaled, &dst);
            }
            SDL_FreeSurface(frame_surface);
            if (!scaled) { ok_all = false; break; }
            SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, scaled);
            SDL_FreeSurface(scaled);
            if (!tex) { ok_all = false; break; }
            app->composer_preview_frames.push_back(tex);
        }
        app->composer_preview_frame_count = (int)app->composer_preview_frames.size();
        app->composer_preview_frame_idx = 0;
        app->composer_preview_fps = (float)page.fps;
        app->composer_preview_playing = ok_all && app->composer_preview_frame_count > 0;
        app->composer_preview_time = 0.0f;
        if (app->composer_preview_frame_count > 0) {
            int w = 0, h = 0;
            SDL_QueryTexture(app->composer_preview_frames[0], nullptr, nullptr, &w, &h);
            app->composer_preview_tex_size = ImVec2((float)w, (float)h);
        }
        double dt = 0.0;
        if (SDL_GetPerformanceFrequency() > 0) {
            Uint64 t1 = SDL_GetPerformanceCounter();
            dt = (double)(t1 - t0) / (double)SDL_GetPerformanceFrequency();
        }
        std::snprintf(app->composer_status,
                      sizeof(app->composer_status),
                      ok_all ? "Preview built (%d frames @ %d fps in %.2fs)" : "Preview failed",
                      app->composer_preview_frame_count,
                      page.fps,
                      dt);
        }
        if (!app->composer_preview_mutate_sim && preview_sim && preview_sim != sim) {
            fdtd_free(preview_sim);
        }
    }
    ImGui::SameLine();
    ImGui::TextUnformatted(app->composer_status);
    if (app->composer_last_export_path[0]) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Copy export folder")) {
            SDL_SetClipboardText(app->composer_last_export_path);
        }
    }
    // Preview controls
    if (!app->composer_preview_frames.empty()) {
        ImGui::Separator();
        ImGui::Text("Preview: %d frames @ %.0f fps", app->composer_preview_frame_count, app->composer_preview_fps);
        if (ImGui::Button(app->composer_preview_playing ? "Pause" : "Play")) {
            app->composer_preview_playing = !app->composer_preview_playing;
        }
        ImGui::SameLine();
        if (ImGui::Button("Stop")) {
            app->composer_preview_playing = false;
            app->composer_preview_frame_idx = 0;
            app->composer_preview_time = 0.0f;
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear Preview")) {
            clear_composer_preview(app);
        }
        ImGui::SliderInt("Frame", &app->composer_preview_frame_idx, 0, app->composer_preview_frame_count - 1);
        if (app->composer_preview_playing && app->composer_preview_frame_count > 0) {
            float dt = ImGui::GetIO().DeltaTime;
            app->composer_preview_time += dt;
            float frame_interval = 1.0f / std::max(1.0f, app->composer_preview_fps);
            while (app->composer_preview_time >= frame_interval) {
                app->composer_preview_time -= frame_interval;
                app->composer_preview_frame_idx = (app->composer_preview_frame_idx + 1) % app->composer_preview_frame_count;
            }
        }
        if (app->composer_preview_frame_idx >= 0 && app->composer_preview_frame_idx < app->composer_preview_frame_count) {
            SDL_Texture* tex = app->composer_preview_frames[app->composer_preview_frame_idx];
            if (tex) {
                ImVec2 avail = ImGui::GetContentRegionAvail();
                float scale = 1.0f;
                if (app->composer_preview_tex_size.x > 1.0f && app->composer_preview_tex_size.y > 1.0f) {
                    scale = std::min(avail.x / app->composer_preview_tex_size.x, 240.0f / app->composer_preview_tex_size.y);
                    if (scale > 2.0f) scale = 2.0f;
                    if (scale < 0.1f) scale = 0.1f;
                }
                ImVec2 tex_sz(app->composer_preview_tex_size.x * scale, app->composer_preview_tex_size.y * scale);
                ImGui::Image((ImTextureID)tex, tex_sz);
            }
        }
    }

    ImGui::Separator();
    ImVec2 canvas_size = ImVec2(ImGui::GetContentRegionAvail().x, 480.0f);
    if (canvas_size.x < 200.0f) canvas_size.x = 200.0f;
    if (canvas_size.y < 200.0f) canvas_size.y = 200.0f;
    ImGui::BeginChild("ComposerCanvas", canvas_size, true, ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoScrollbar);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
    float scale = std::min(canvas_size.x / (float)page.res_w, canvas_size.y / (float)page.res_h);
    ImVec2 page_px(page.res_w * scale, page.res_h * scale);
    ImVec2 page_origin = ImVec2(canvas_pos.x + (canvas_size.x - page_px.x) * 0.5f,
                                canvas_pos.y + (canvas_size.y - page_px.y) * 0.5f);
    ImVec2 page_min = page_origin;
    ImVec2 page_max = ImVec2(page_origin.x + page_px.x, page_origin.y + page_px.y);
    dl->AddRectFilled(page_min, page_max, IM_COL32((int)(page.bg.x * 255), (int)(page.bg.y * 255), (int)(page.bg.z * 255), 255));
    dl->AddRect(page_min, page_max, IM_COL32(255, 255, 255, 200), 4.0f, 0, 2.0f);
    if (app->composer_canvas_grid && app->composer_canvas_grid_step > 1) {
        int step = std::max(4, app->composer_canvas_grid_step);
        float step_px = step * scale;
        for (float x = page_min.x; x <= page_max.x + 1.0f; x += step_px) {
            dl->AddLine(ImVec2(x, page_min.y), ImVec2(x, page_max.y), IM_COL32(255, 255, 255, 40));
        }
        for (float y = page_min.y; y <= page_max.y + 1.0f; y += step_px) {
            dl->AddLine(ImVec2(page_min.x, y), ImVec2(page_max.x, y), IM_COL32(255, 255, 255, 40));
        }
    }

    ImGui::InvisibleButton("canvas_btn", canvas_size, ImGuiButtonFlags_MouseButtonLeft);
    bool canvas_hovered = ImGui::IsItemHovered();
    ImVec2 mouse = ImGui::GetIO().MousePos;
    ImVec2 mouse_page = ImVec2((mouse.x - page_origin.x) / scale, (mouse.y - page_origin.y) / scale);

    if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) &&
        ImGui::IsKeyPressed(ImGuiKey_Delete)) {
        page.items.erase(std::remove_if(page.items.begin(), page.items.end(),
                                        [](const ComposerItem& it) { return it.selected; }),
                         page.items.end());
    }

    for (size_t i = 0; i < page.items.size(); ++i) {
        ComposerItem& it = page.items[i];
        ImVec2 minp = ImVec2(page_origin.x + it.pos.x * scale,
                             page_origin.y + it.pos.y * scale);
        ImVec2 maxp = ImVec2(minp.x + it.size.x * scale,
                             minp.y + it.size.y * scale);
        ImU32 col = composer_item_color(it);
        dl->AddRectFilled(minp, maxp, col);
        dl->AddRect(minp, maxp, IM_COL32(255, 255, 255, it.selected ? 255 : 140), 2.0f, 0, it.selected ? 3.0f : 1.5f);
        char label[64];
        const char* type_name = "Item";
        switch (it.type) {
            case COMPOSER_FIELD_VIEW: type_name = "Field"; break;
            case COMPOSER_REGION: type_name = "Region"; break;
            case COMPOSER_SCOPE: type_name = "Scope"; break;
            case COMPOSER_FFT: type_name = "FFT"; break;
            case COMPOSER_LEGEND: type_name = "Legend"; break;
            case COMPOSER_MEAS: type_name = "Measurements"; break;
            case COMPOSER_SMITH: type_name = "Smith"; break;
        }
        std::snprintf(label, sizeof(label), "%s #%d", type_name, it.id);
        dl->AddText(ImVec2(minp.x + 6, minp.y + 6), IM_COL32(20, 20, 20, 255), label);

        if (it.selected) {
            const float handle_sz = 10.0f;
            ImVec2 corners[4] = {minp,
                                 ImVec2(maxp.x, minp.y),
                                 ImVec2(minp.x, maxp.y),
                                 maxp};
            for (int h = 0; h < 4; ++h) {
                ImVec2 hc = corners[h];
                ImVec2 hmin = ImVec2(hc.x - handle_sz * 0.5f, hc.y - handle_sz * 0.5f);
                ImVec2 hmax = ImVec2(hc.x + handle_sz * 0.5f, hc.y + handle_sz * 0.5f);
                bool hhover = mouse.x >= hmin.x && mouse.x <= hmax.x &&
                              mouse.y >= hmin.y && mouse.y <= hmax.y;
                ImU32 hcol = hhover ? IM_COL32(255, 200, 80, 255) : IM_COL32(240, 240, 240, 230);
                if (hhover) {
                    ImGuiMouseCursor cur = (h == 0 || h == 3) ? ImGuiMouseCursor_ResizeNWSE : ImGuiMouseCursor_ResizeNESW;
                    ImGui::SetMouseCursor(cur);
                    if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                        for (auto& other : page.items) other.selected = false;
                        it.selected = true;
                        app->composer_resizing = true;
                        app->composer_resize_item = (int)i;
                        app->composer_resize_corner = h;
                        app->composer_resize_start_pos = it.pos;
                        app->composer_resize_start_size = it.size;
                        app->composer_resize_start_mouse = mouse_page;
                        app->composer_dragging = false;
                        app->composer_drag_item = -1;
                    }
                }
                dl->AddRectFilled(hmin, hmax, hcol, 2.0f);
                dl->AddRect(hmin, hmax, IM_COL32(30, 30, 30, 230), 2.0f);
            }
        }

        bool hovered = canvas_hovered &&
                       mouse.x >= minp.x && mouse.x <= maxp.x &&
                       mouse.y >= minp.y && mouse.y <= maxp.y;
        if (!app->composer_resizing && hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            bool ctrl = ImGui::GetIO().KeyCtrl;
            if (!ctrl) {
                for (auto& other : page.items) other.selected = false;
                it.selected = true;
            } else {
                it.selected = !it.selected;
            }
            app->composer_dragging = true;
            app->composer_drag_item = (int)i;
            app->composer_drag_offset = ImVec2(mouse_page.x - it.pos.x, mouse_page.y - it.pos.y);
        }
    }

    if (app->composer_resizing && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        if (app->composer_resize_item >= 0 && app->composer_resize_item < (int)page.items.size()) {
            ComposerItem& it = page.items[app->composer_resize_item];
            ImVec2 delta = ImVec2(mouse_page.x - app->composer_resize_start_mouse.x,
                                  mouse_page.y - app->composer_resize_start_mouse.y);
            float left = app->composer_resize_start_pos.x;
            float right = left + app->composer_resize_start_size.x;
            float top = app->composer_resize_start_pos.y;
            float bottom = top + app->composer_resize_start_size.y;
            switch (app->composer_resize_corner) {
                case 0: left += delta.x; top += delta.y; break;
                case 1: right += delta.x; top += delta.y; break;
                case 2: left += delta.x; bottom += delta.y; break;
                case 3: right += delta.x; bottom += delta.y; break;
                default: break;
            }
            const float kMin = 24.0f;
            if (right - left < kMin) {
                if (app->composer_resize_corner == 0 || app->composer_resize_corner == 2) {
                    left = right - kMin;
                } else {
                    right = left + kMin;
                }
            }
            if (bottom - top < kMin) {
                if (app->composer_resize_corner == 0 || app->composer_resize_corner == 1) {
                    top = bottom - kMin;
                } else {
                    bottom = top + kMin;
                }
            }
            if (app->composer_snap) {
                float step = std::max(1.0f, app->composer_snap_step);
                left = std::round(left / step) * step;
                right = std::round(right / step) * step;
                top = std::round(top / step) * step;
                bottom = std::round(bottom / step) * step;
            }
            left = std::max(left, 0.0f);
            top = std::max(top, 0.0f);
            right = std::min(right, (float)page.res_w);
            bottom = std::min(bottom, (float)page.res_h);
            if (right - left < kMin) right = left + kMin;
            if (bottom - top < kMin) bottom = top + kMin;
            it.pos = ImVec2(left, top);
            it.size = ImVec2(right - left, bottom - top);
        }
    } else if (app->composer_resizing && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
        app->composer_resizing = false;
        app->composer_resize_item = -1;
        app->composer_resize_corner = -1;
    } else if (app->composer_dragging && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        if (app->composer_drag_item >= 0 && app->composer_drag_item < (int)page.items.size()) {
            ComposerItem& drag = page.items[app->composer_drag_item];
            drag.pos = ImVec2(mouse_page.x - app->composer_drag_offset.x,
                              mouse_page.y - app->composer_drag_offset.y);
            if (app->composer_snap) {
                float step = std::max(1.0f, app->composer_snap_step);
                drag.pos.x = std::round(drag.pos.x / step) * step;
                drag.pos.y = std::round(drag.pos.y / step) * step;
            }
            if (drag.pos.x < 0) drag.pos.x = 0;
            if (drag.pos.y < 0) drag.pos.y = 0;
            if (drag.pos.x + drag.size.x > page.res_w) drag.pos.x = page.res_w - drag.size.x;
            if (drag.pos.y + drag.size.y > page.res_h) drag.pos.y = page.res_h - drag.size.y;
        }
    }
    if ((app->composer_dragging && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) ||
        (!ImGui::IsMouseDown(ImGuiMouseButton_Left) && app->composer_dragging)) {
        app->composer_dragging = false;
        app->composer_drag_item = -1;
    }

    ImGui::EndChild();

    ImGui::Separator();
    // Item properties
    ComposerItem* sel = nullptr;
    for (auto& it : page.items) {
        if (it.selected) { sel = &it; break; }
    }
        if (sel) {
            ImGui::Text("Selected: #%d", sel->id);
            ImGui::TextUnformatted("Type:");
            ImGui::SameLine();
            const char* type_name = "Item";
            switch (sel->type) {
                case COMPOSER_FIELD_VIEW: type_name = "Field"; break;
                case COMPOSER_REGION: type_name = "Region"; break;
                case COMPOSER_SCOPE: type_name = "Scope"; break;
                case COMPOSER_FFT: type_name = "FFT"; break;
                case COMPOSER_LEGEND: type_name = "Legend"; break;
                case COMPOSER_MEAS: type_name = "Measurements"; break;
                case COMPOSER_SMITH: type_name = "Smith"; break;
            }
            ImGui::TextUnformatted(type_name);
            ImGui::InputFloat2("Pos", (float*)&sel->pos);
            ImGui::InputFloat2("Size", (float*)&sel->size);
        if (sel->type == COMPOSER_FIELD_VIEW) {
            ImGui::InputInt("Viewport", &sel->viewport_idx);
            if (sel->viewport_idx < 0) sel->viewport_idx = 0;
            if (sel->viewport_idx > 3) sel->viewport_idx = 3;
        } else if (sel->type == COMPOSER_REGION) {
                ImGui::InputInt("Viewport", &sel->viewport_idx);
                if (sel->viewport_idx < 0) sel->viewport_idx = 0;
                if (sel->viewport_idx > 3) sel->viewport_idx = 3;
                float r[4] = {sel->region_norm.x, sel->region_norm.y, sel->region_norm.z, sel->region_norm.w};
                if (ImGui::InputFloat4("Region (x0,y0,x1,y1)", r)) {
                    sel->region_norm = ImVec4(r[0], r[1], r[2], r[3]);
                    if (sel->region_norm.x < 0.0f) sel->region_norm.x = 0.0f;
                    if (sel->region_norm.y < 0.0f) sel->region_norm.y = 0.0f;
                    if (sel->region_norm.z > 1.0f) sel->region_norm.z = 1.0f;
                    if (sel->region_norm.w > 1.0f) sel->region_norm.w = 1.0f;
                    if (sel->region_norm.z < sel->region_norm.x + 0.001f) sel->region_norm.z = sel->region_norm.x + 0.001f;
                    if (sel->region_norm.w < sel->region_norm.y + 0.001f) sel->region_norm.w = sel->region_norm.y + 0.001f;
                }
                if (ImGui::SmallButton("Re-pick region")) {
                    app->composer_pick_region_active = true;
                    app->composer_pick_region_dragging = false;
                    app->composer_pick_region_viewport = sel->viewport_idx;
                    app->composer_pick_region_page = app->composer_active_page;
                    app->composer_pick_region_target_item = sel->id;
                    app->composer_pick_region_start_norm = ImVec2(sel->region_norm.x, sel->region_norm.y);
                    app->composer_pick_region_end_norm = ImVec2(sel->region_norm.z, sel->region_norm.w);
                    std::snprintf(app->composer_status, sizeof(app->composer_status), "Pick region for item #%d", sel->id);
                }
            }
            if (ImGui::Button("Delete Item")) {
                page.items.erase(std::remove_if(page.items.begin(), page.items.end(),
                                                [&](const ComposerItem& it) { return it.id == sel->id; }),
                                 page.items.end());
            sel = nullptr;
        }
    } else {
        ImGui::TextUnformatted("Select an item to edit properties.");
    }

    ImGui::End();
}

// Lightweight debug stubs (no-op in release builds).
static void debug_logf(const char* fmt, ...) {
    (void)fmt;
}
static void debug_log_close() {}

struct WizardState {
    SimulationConfig cfg;
    bool open;
    bool advanced;
};

static void apply_wizard_materials_to_sim(const WizardState& wizard,
                                          SimulationBootstrap* bootstrap,
                                          SimulationState* sim);

static void wizard_init_from_config(WizardState& w, const SimulationConfig* cfg) {
    if (cfg) {
        w.cfg = *cfg;
    } else {
        w.cfg = SIM_CONFIG_DEFAULTS;
    }
    w.open = true;
    w.advanced = false;
}

static void ui_log_add(AppState* app, const char* fmt, ...) {
    if (!app) return;
    if (app->log_count >= 128) {
        app->log_count = 0;
    }
    va_list args;
    va_start(args, fmt);
    vsnprintf(app->log_lines[app->log_count], sizeof(app->log_lines[app->log_count]), fmt, args);
    va_end(args);
    app->log_count++;
}

static const char* source_type_label(SourceType t) {
    switch (t) {
        case SRC_CW:           return "CW";
        case SRC_GAUSS_PULSE:  return "Gaussian";
        case SRC_RICKER:       return "Ricker";
        case SRC_EXPR:         return "Expr";
        default:               return "Unknown";
    }
}

static const char* source_field_label(SourceFieldType f) {
    switch (f) {
        case SRC_FIELD_EZ: return "Ez";
        case SRC_FIELD_HX: return "Hx";
        case SRC_FIELD_HY: return "Hy";
        default:           return "Ez";
    }
}

static double clamp01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

static int normalized_to_cell_index(double frac, int n) {
    if (n <= 0) return 0;
    double clamped = clamp01(frac);
    int idx = (int)std::lround(clamped * (double)(n - 1));
    const int pad = 2;
    if (idx < pad) idx = pad;
    if (idx >= n - pad && n > pad) idx = n - pad - 1;
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    return idx;
}

static ImVec2 compute_center_offset(const ViewportInstance& vp,
                                    const SimulationState* sim,
                                    int scale) {
    ImVec2 offset(0.0f, 0.0f);
    if (!sim || scale <= 0) return offset;
    float field_px_w = (float)(sim->nx * scale);
    float field_px_h = (float)(sim->ny * scale);
    if (vp.size.x > field_px_w) {
        offset.x = (vp.size.x - field_px_w) * 0.5f;
    }
    if (vp.size.y > field_px_h) {
        offset.y = (vp.size.y - field_px_h) * 0.5f;
    }
    return offset;
}

static ImVec2 compute_viewport_offset(const ViewportInstance& vp,
                                      const SimulationState* sim,
                                      int scale) {
    ImVec2 centered = compute_center_offset(vp, sim, scale);
    return ImVec2(vp.pan_x + centered.x, vp.pan_y + centered.y);
}

static int compute_grid_step_from_scale(int scale) {
    if (scale < 2) return 10;
    if (scale < 8) return 5;
    return 1;
}

static void compute_viewport_layout(AppState* app,
                                    ImVec2 container_size) {
    if (!app) return;
    const float gap = 4.0f;
    const float min_dim = 200.0f;
    for (int i = 0; i < 4; ++i) {
        app->viewports[i].valid = false;
        app->viewports[i].active = false;
        app->viewports[i].pos = ImVec2(0.0f, 0.0f);
        app->viewports[i].size = ImVec2(0.0f, 0.0f);
    }

    switch (app->viewport_layout) {
        case VIEWPORT_SINGLE: {
            app->viewports[0].pos = ImVec2(0.0f, 0.0f);
            app->viewports[0].size = container_size;
            app->viewports[0].valid = (container_size.x > 4.0f && container_size.y > 4.0f);
        } break;
        case VIEWPORT_HORIZONTAL: {
            float h = (container_size.y - gap) * 0.5f;
            if (h < min_dim) {
                app->viewport_layout = VIEWPORT_SINGLE;
                compute_viewport_layout(app, container_size);
                return;
            }
            app->viewports[0].pos = ImVec2(0.0f, 0.0f);
            app->viewports[0].size = ImVec2(container_size.x, h);
            app->viewports[0].valid = true;
            app->viewports[1].pos = ImVec2(0.0f, h + gap);
            app->viewports[1].size = ImVec2(container_size.x, h);
            app->viewports[1].valid = true;
        } break;
        case VIEWPORT_VERTICAL: {
            float w = (container_size.x - gap) * 0.5f;
            if (w < min_dim) {
                app->viewport_layout = VIEWPORT_SINGLE;
                compute_viewport_layout(app, container_size);
                return;
            }
            app->viewports[0].pos = ImVec2(0.0f, 0.0f);
            app->viewports[0].size = ImVec2(w, container_size.y);
            app->viewports[0].valid = true;
            app->viewports[1].pos = ImVec2(w + gap, 0.0f);
            app->viewports[1].size = ImVec2(w, container_size.y);
            app->viewports[1].valid = true;
        } break;
        case VIEWPORT_QUAD: {
            float w = (container_size.x - gap) * 0.5f;
            float h = (container_size.y - gap) * 0.5f;
            if (w < min_dim || h < min_dim) {
                app->viewport_layout = VIEWPORT_SINGLE;
                compute_viewport_layout(app, container_size);
                return;
            }
            app->viewports[0].pos = ImVec2(0.0f, 0.0f);
            app->viewports[0].size = ImVec2(w, h);
            app->viewports[0].valid = true;
            app->viewports[1].pos = ImVec2(w + gap, 0.0f);
            app->viewports[1].size = ImVec2(w, h);
            app->viewports[1].valid = true;
            app->viewports[2].pos = ImVec2(0.0f, h + gap);
            app->viewports[2].size = ImVec2(w, h);
            app->viewports[2].valid = true;
            app->viewports[3].pos = ImVec2(w + gap, h + gap);
            app->viewports[3].size = ImVec2(w, h);
            app->viewports[3].valid = true;
        } break;
        default:
            break;
    }
}

static int get_viewport_at_mouse(const AppState& app, int mx, int my) {
    float local_x = (float)mx - app.viewport_pos.x;
    float local_y = (float)my - app.viewport_pos.y;
    for (int i = 0; i < 4; ++i) {
        const ViewportInstance& vp = app.viewports[i];
        if (!vp.valid) continue;
        if (local_x >= vp.pos.x && local_x < vp.pos.x + vp.size.x &&
            local_y >= vp.pos.y && local_y < vp.pos.y + vp.size.y) {
            return i;
        }
    }
    return -1;
}

enum LayoutPreset {
    LAYOUT_BEGINNER = 0,
    LAYOUT_POWER_USER = 1,
    LAYOUT_ANALYSIS = 2,
    LAYOUT_CANVAS_FIRST = 3
};

static constexpr float kDefaultViewportZoom = 2.0f;

static void apply_layout_preset(ImGuiID dockspace_id, LayoutPreset preset) {
    if (dockspace_id == 0) return;
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);

    ImGuiID dock_center = dockspace_id;
    ImGuiID dock_bottom = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Down, 0.25f, nullptr, &dock_center);
    ImGuiID dock_left = dock_center;
    ImGuiID dock_right = dock_center;

    switch (preset) {
        case LAYOUT_BEGINNER: {
            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Run Controls", dock_right);
            ImGui::DockBuilderDockWindow("Sources", dock_right);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("Log", dock_bottom);
            break;
        }
        case LAYOUT_POWER_USER: {
            dock_left = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Left, 0.2f, nullptr, &dock_center);
            dock_right = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Right, 0.22f, nullptr, &dock_center);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Sources", dock_left);
            ImGui::DockBuilderDockWindow("Materials / Blocks", dock_left);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("Log", dock_bottom);
            ImGui::DockBuilderDockWindow("Material Legend", dock_right);
            ImGui::DockBuilderDockWindow("Run Controls", dock_right);
            break;
        }
        case LAYOUT_ANALYSIS: {
            dock_right = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Right, 0.4f, nullptr, &dock_center);
            ImGuiID dock_right_bottom = ImGui::DockBuilderSplitNode(dock_right, ImGuiDir_Down, 0.5f, nullptr, &dock_right);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("S-Parameters", dock_right);
            ImGui::DockBuilderDockWindow("Smith Chart", dock_right_bottom);
            ImGui::DockBuilderDockWindow("Log", dock_bottom);
            break;
        }
        case LAYOUT_CANVAS_FIRST: {
            dock_left = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Left, 0.16f, nullptr, &dock_center);
            dock_right = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Right, 0.18f, nullptr, &dock_center);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Sources", dock_left);
            ImGui::DockBuilderDockWindow("Materials / Blocks", dock_left);
            ImGui::DockBuilderDockWindow("Grid & Domain", dock_left);
            ImGui::DockBuilderDockWindow("Run Controls", dock_right);
            ImGui::DockBuilderDockWindow("Material Legend", dock_right);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("Log", dock_bottom);
            break;
        }
        default:
            break;
    }

    ImGui::DockBuilderFinish(dockspace_id);
}

static void clamp_pan_offset(ViewportInstance* vp,
                             const SimulationState* sim,
                             int scale) {
    if (!vp || !sim || scale <= 0) return;
    ImVec2 center_offset = compute_center_offset(*vp, sim, scale);
    float field_w = (float)(sim->nx * scale);
    float field_h = (float)(sim->ny * scale);
    float viewport_w = vp->size.x;
    float viewport_h = vp->size.y;

    float min_pan_x = -(field_w * 0.9f);
    float max_pan_x = viewport_w * 0.9f;
    float min_pan_y = -(field_h * 0.9f);
    float max_pan_y = viewport_h * 0.9f;

    float effective_pan_x = vp->pan_x + center_offset.x;
    float effective_pan_y = vp->pan_y + center_offset.y;

    effective_pan_x = std::clamp(effective_pan_x, min_pan_x, max_pan_x);
    effective_pan_y = std::clamp(effective_pan_y, min_pan_y, max_pan_y);

    vp->pan_x = std::round(effective_pan_x - center_offset.x);
    vp->pan_y = std::round(effective_pan_y - center_offset.y);
}

static void apply_zoom_at_point(ViewportInstance* vp,
                                RenderContext* render,
                                const SimulationState* sim,
                                int* scale,
                                float focus_x,
                                float focus_y,
                                float zoom_multiplier,
                                bool smooth) {
    if (!vp || !render || !sim || !scale) return;
    if (*scale <= 0) return;

    ImVec2 center_before = compute_center_offset(*vp, sim, *scale);
    float field_x = focus_x - (vp->pan_x + center_before.x);
    float field_y = focus_y - (vp->pan_y + center_before.y);
    float grid_x = field_x / (float)(*scale);
    float grid_y = field_y / (float)(*scale);

    const float target_zoom = std::clamp(vp->zoom * zoom_multiplier, 0.5f, 32.0f);
    float new_zoom = target_zoom;
    if (smooth) {
        const float zoom_alpha = 0.25f;  // ease wheel only; buttons stay snappy
        new_zoom = ImLerp(vp->zoom, target_zoom, zoom_alpha);
    }
    vp->zoom = new_zoom;

    *scale = (int)std::lround(new_zoom);
    if (*scale < 1) *scale = 1;
    render->scale = *scale;

    ImVec2 center_after = compute_center_offset(*vp, sim, *scale);
    float effective_pan_x = focus_x - grid_x * (float)(*scale);
    float effective_pan_y = focus_y - grid_y * (float)(*scale);
    vp->pan_x = effective_pan_x - center_after.x;
    vp->pan_y = effective_pan_y - center_after.y;

    clamp_pan_offset(vp, sim, *scale);
}

static void reset_view_transform(ViewportInstance* vp,
                                 RenderContext* render,
                                 const SimulationState* sim,
                                 int* scale,
                                 float default_zoom) {
    if (!vp || !render || !sim || !scale) return;
    vp->zoom = default_zoom;
    *scale = (int)std::lround(default_zoom);
    if (*scale < 1) *scale = 1;
    render->scale = *scale;
    vp->pan_x = 0.0f;
    vp->pan_y = 0.0f;
}

static SDL_Surface* capture_frame(SDL_Renderer* renderer,
                                  int width,
                                  int height,
                                  float scale_factor) {
    if (!renderer) return nullptr;
    int out_w = (int)std::lround(width * scale_factor);
    int out_h = (int)std::lround(height * scale_factor);
    if (out_w < 2) out_w = 2;
    if (out_h < 2) out_h = 2;
    out_w = (out_w / 2) * 2;
    out_h = (out_h / 2) * 2;

    SDL_Surface* surface = SDL_CreateRGBSurface(0, out_w, out_h, 24,
                                                0x00FF0000,
                                                0x0000FF00,
                                                0x000000FF,
                                                0);
    if (!surface) return nullptr;

    SDL_Rect src = {0, 0, width, height};
    if (SDL_RenderReadPixels(renderer, &src, SDL_PIXELFORMAT_RGB24, surface->pixels, surface->pitch) != 0) {
        SDL_FreeSurface(surface);
        return nullptr;
    }

    if (out_w != width || out_h != height) {
        SDL_Surface* scaled = SDL_CreateRGBSurface(0, out_w, out_h, 24,
                                                   0x00FF0000,
                                                   0x0000FF00,
                                                   0x000000FF,
                                                   0);
        if (!scaled) {
            SDL_FreeSurface(surface);
            return nullptr;
        }
        SDL_BlitScaled(surface, nullptr, scaled, nullptr);
        SDL_FreeSurface(surface);
        surface = scaled;
    }
    return surface;
}

static std::string ffmpeg_base_dir() {
    char* base = SDL_GetBasePath();
    if (!base) return std::string();
    std::string path(base);
    SDL_free(base);
    return path;
}

static const std::string& ffmpeg_command() {
    static bool initialized = false;
    static std::string cmd;
    if (initialized) return cmd;
    initialized = true;

    const char* env_path = std::getenv("FFMPEG_PATH");
    std::vector<std::filesystem::path> candidates;
    if (env_path && env_path[0] != '\0') {
        candidates.emplace_back(env_path);
    }
#ifdef _WIN32
    candidates.emplace_back("ffmpeg\\bin\\ffmpeg.exe");
    candidates.emplace_back("ffmpeg\\ffmpeg.exe");
    candidates.emplace_back(".\\ffmpeg.exe");
    candidates.emplace_back("..\\ffmpeg\\bin\\ffmpeg.exe");
    candidates.emplace_back("..\\..\\ffmpeg\\bin\\ffmpeg.exe");
#else
    candidates.emplace_back("ffmpeg/bin/ffmpeg");
    candidates.emplace_back("ffmpeg/ffmpeg");
    candidates.emplace_back("./ffmpeg");
    candidates.emplace_back("../ffmpeg/bin/ffmpeg");
    candidates.emplace_back("../../ffmpeg/bin/ffmpeg");
#endif

    std::string base_dir = ffmpeg_base_dir();
    if (!base_dir.empty()) {
        std::filesystem::path base_path(base_dir);
#ifdef _WIN32
        candidates.emplace_back(base_path / "ffmpeg" / "bin" / "ffmpeg.exe");
#else
        candidates.emplace_back(base_path / "ffmpeg" / "bin" / "ffmpeg");
#endif
    }

    for (const auto& c : candidates) {
        if (std::filesystem::exists(c)) {
            std::filesystem::path abs = std::filesystem::absolute(c);
            cmd = abs.string();
            return cmd;
        }
    }

#ifdef _WIN32
    if (std::system("ffmpeg -version >NUL 2>&1") == 0) {
        cmd = "ffmpeg";
    }
#else
    if (std::system("ffmpeg -version >/dev/null 2>&1") == 0) {
        cmd = "ffmpeg";
    }
#endif
    return cmd;
}

static bool has_ffmpeg() {
    return !ffmpeg_command().empty();
}

static bool save_surface_png(const std::filesystem::path& path, SDL_Surface* surf) {
    if (!surf) return false;
    SDL_Surface* conv = SDL_ConvertSurfaceFormat(surf, SDL_PIXELFORMAT_RGBA32, 0);
    if (!conv) return false;
    int ok = stbi_write_png(path.string().c_str(),
                            conv->w,
                            conv->h,
                            4,
                            conv->pixels,
                            conv->pitch);
    SDL_FreeSurface(conv);
    return ok == 1;
}

static void clear_composer_preview(AppState* app) {
    if (!app) return;
    if (app->composer_preview_tex) {
        SDL_DestroyTexture(app->composer_preview_tex);
        app->composer_preview_tex = nullptr;
    }
    for (SDL_Texture* t : app->composer_preview_frames) {
        if (t) SDL_DestroyTexture(t);
    }
    app->composer_preview_frames.clear();
    app->composer_preview_frame_idx = 0;
    app->composer_preview_frame_count = 0;
    app->composer_preview_playing = false;
    app->composer_preview_time = 0.0f;
    app->composer_preview_tex_size = ImVec2(0, 0);
    app->composer_preview_fps = 30.0f;
}

static SimulationState* clone_simulation_state(const SimulationState* src) {
    if (!src) return nullptr;
    SimulationState* dst = fdtd_init(&src->config);
    if (!dst) return nullptr;
    int nxy = src->nx * src->ny;
    auto copy_buf = [](void* dstb, const void* srcb, size_t bytes) {
        if (dstb && srcb && bytes > 0) std::memcpy(dstb, srcb, bytes);
    };
    copy_buf(dst->Ez_data, src->Ez_data, sizeof(double) * nxy);
    copy_buf(dst->Hx_data, src->Hx_data, sizeof(double) * nxy);
    copy_buf(dst->Hy_data, src->Hy_data, sizeof(double) * nxy);
    copy_buf(dst->Ez_old_data, src->Ez_old_data, sizeof(double) * nxy);
    copy_buf(dst->psi_Ezx_data, src->psi_Ezx_data, sizeof(double) * nxy);
    copy_buf(dst->psi_Ezy_data, src->psi_Ezy_data, sizeof(double) * nxy);
    copy_buf(dst->psi_Hyx_data, src->psi_Hyx_data, sizeof(double) * nxy);
    copy_buf(dst->psi_Hxy_data, src->psi_Hxy_data, sizeof(double) * nxy);
    copy_buf(dst->epsr_data, src->epsr_data, sizeof(double) * nxy);
    copy_buf(dst->sigma_map_data, src->sigma_map_data, sizeof(double) * nxy);
    copy_buf(dst->tag_grid_data, src->tag_grid_data, sizeof(unsigned char) * nxy);
    dst->timestep = src->timestep;
    dst->freq = src->freq;
    dst->step_Ez_absmax = src->step_Ez_absmax;
    dst->ports_on = src->ports_on;
    std::memcpy(dst->sources, src->sources, sizeof(dst->sources));
    for (int i = 0; i < MAX_PORTS; ++i) {
        dst->ports[i].active = src->ports[i].active;
        dst->ports[i].head = src->ports[i].head;
        dst->ports[i].n = src->ports[i].n;
        dst->ports[i].len = src->ports[i].len;
        dst->ports[i].x = src->ports[i].x;
        dst->ports[i].y0 = src->ports[i].y0;
        dst->ports[i].y1 = src->ports[i].y1;
        copy_buf(dst->ports[i].V, src->ports[i].V, sizeof(double) * src->ports[i].n);
        copy_buf(dst->ports[i].I, src->ports[i].I, sizeof(double) * src->ports[i].n);
    }
    // CPML buffers
    copy_buf(dst->cpml.kx, src->cpml.kx, sizeof(double) * dst->cpml.kx_capacity);
    copy_buf(dst->cpml.bx, src->cpml.bx, sizeof(double) * dst->cpml.kx_capacity);
    copy_buf(dst->cpml.cx, src->cpml.cx, sizeof(double) * dst->cpml.kx_capacity);
    copy_buf(dst->cpml.ky, src->cpml.ky, sizeof(double) * dst->cpml.ky_capacity);
    copy_buf(dst->cpml.by, src->cpml.by, sizeof(double) * dst->cpml.ky_capacity);
    copy_buf(dst->cpml.cy, src->cpml.cy, sizeof(double) * dst->cpml.ky_capacity);
    dst->cpml.enabled = src->cpml.enabled;
    dst->cpml.thickness = src->cpml.thickness;
    dst->cpml.preset_idx = src->cpml.preset_idx;
    dst->cpml.sigma_max = src->cpml.sigma_max;
    dst->cpml.kappa_max = src->cpml.kappa_max;
    dst->cpml.alpha_max = src->cpml.alpha_max;
    dst->cpml.boundary_type = src->cpml.boundary_type;
    return dst;
}

// Minimal JSON parser for composer layouts (expects array "pages" with fields we use)
static bool composer_load_layout_json(AppState* app, const char* path) {
    if (!app || !path) return false;
    std::error_code ec;
    std::string data;
    if (!std::filesystem::exists(path, ec)) return false;
    {
        std::ifstream f(path, std::ios::in | std::ios::binary);
        if (!f) return false;
        std::ostringstream oss;
        oss << f.rdbuf();
        data = oss.str();
    }
    if (data.empty()) return false;
    // Very naive JSON-ish reader for our saved format. This is intentionally permissive:
    // it looks for keys and brackets without needing strict JSON.
    auto read_int = [&](const std::string& key, const std::string& src, size_t from, int defv) -> int {
        std::string pat = "\"" + key + "\"";
        size_t pos = src.find(pat, from);
        if (pos == std::string::npos) return defv;
        pos = src.find(":", pos);
        if (pos == std::string::npos) return defv;
        size_t start = src.find_first_of("-0123456789", pos + 1);
        if (start == std::string::npos) return defv;
        size_t end = src.find_first_not_of("0123456789", start);
        return std::atoi(src.substr(start, end - start).c_str());
    };
    auto read_float = [&](const std::string& key, const std::string& src, size_t from, float defv) -> float {
        std::string pat = "\"" + key + "\"";
        size_t pos = src.find(pat, from);
        if (pos == std::string::npos) return defv;
        pos = src.find(":", pos);
        if (pos == std::string::npos) return defv;
        size_t start = src.find_first_of("-0123456789.", pos + 1);
        if (start == std::string::npos) return defv;
        size_t end = src.find_first_not_of("0123456789.", start);
        return std::atof(src.substr(start, end - start).c_str());
    };
    auto read_string = [&](const std::string& key, const std::string& src, size_t from) -> std::string {
        std::string pat = "\"" + key + "\"";
        size_t pos = src.find(pat, from);
        if (pos == std::string::npos) return std::string();
        pos = src.find(":", pos);
        if (pos == std::string::npos) return std::string();
        size_t start = src.find("\"", pos);
        if (start == std::string::npos) return std::string();
        size_t end = src.find("\"", start + 1);
        if (end == std::string::npos) return std::string();
        return src.substr(start + 1, end - start - 1);
    };
    auto read_vec2 = [&](const std::string& key, const std::string& src, size_t from, ImVec2 defv) -> ImVec2 {
        std::string pat = "\"" + key + "\"";
        size_t pos = src.find(pat, from);
        if (pos == std::string::npos) return defv;
        size_t lb = src.find("[", pos);
        size_t rb = (lb != std::string::npos) ? src.find("]", lb) : std::string::npos;
        if (lb == std::string::npos || rb == std::string::npos) return defv;
        float a = defv.x, b = defv.y;
        std::sscanf(src.substr(lb + 1, rb - lb - 1).c_str(), "%f , %f", &a, &b);
        return ImVec2(a, b);
    };
    auto read_vec4 = [&](const std::string& key, const std::string& src, size_t from, ImVec4 defv) -> ImVec4 {
        std::string pat = "\"" + key + "\"";
        size_t pos = src.find(pat, from);
        if (pos == std::string::npos) return defv;
        size_t lb = src.find("[", pos);
        size_t rb = (lb != std::string::npos) ? src.find("]", lb) : std::string::npos;
        if (lb == std::string::npos || rb == std::string::npos) return defv;
        float a = defv.x, b = defv.y, c = defv.z, d = defv.w;
        std::sscanf(src.substr(lb + 1, rb - lb - 1).c_str(), "%f , %f , %f , %f", &a, &b, &c, &d);
        return ImVec4(a, b, c, d);
    };
    auto read_bool = [&](const std::string& key, const std::string& src, size_t from, bool defv) -> bool {
        return read_int(key, src, from, defv ? 1 : 0) != 0;
    };

    size_t pages_pos = data.find("\"pages\"");
    if (pages_pos == std::string::npos) return false;
    size_t arr_pos = data.find("[", pages_pos);
    if (arr_pos == std::string::npos) return false;

    std::vector<ComposerPage> loaded;
    int max_item_id = 1;
    size_t pos = arr_pos;
    while (true) {
        size_t obj_start = data.find("{", pos);
        if (obj_start == std::string::npos) break;
        int brace = 0;
        size_t obj_end = std::string::npos;
        for (size_t i = obj_start; i < data.size(); ++i) {
            if (data[i] == '{') brace++;
            else if (data[i] == '}') {
                brace--;
                if (brace == 0) {
                    obj_end = i;
                    break;
                }
            }
        }
        if (obj_end == std::string::npos) break;
        std::string obj = data.substr(obj_start, obj_end - obj_start + 1);

        ComposerPage p{};
        std::snprintf(p.name, sizeof(p.name), "Page %zu", loaded.size() + 1);
        std::string pname = read_string("name", obj, 0);
        if (!pname.empty()) std::snprintf(p.name, sizeof(p.name), "%s", pname.c_str());
        p.res_w = read_int("res_w", obj, 0, 1280);
        p.res_h = read_int("res_h", obj, 0, 720);
        std::string oname = read_string("output_name", obj, 0);
        if (!oname.empty()) std::snprintf(p.output_name, sizeof(p.output_name), "%s", oname.c_str());
        p.output_format = read_int("output_format", obj, 0, 0);
        p.fps = read_int("fps", obj, 0, 30);
        p.frames = read_int("frames", obj, 0, 60);
        p.video_kbps = read_int("video_kbps", obj, 0, 4000);
        p.transparent_bg = read_bool("transparent_bg", obj, 0, false);
        ImVec4 bg_from_vec = read_vec4("bg", obj, 0, ImVec4(-1, -1, -1, -1));
        if (bg_from_vec.x >= 0.0f) {
            p.bg = bg_from_vec;
        } else {
            float bgx = read_float("bg_r", obj, 0, 0.08f);
            float bgy = read_float("bg_g", obj, 0, 0.08f);
            float bgz = read_float("bg_b", obj, 0, 0.10f);
            float bga = read_float("bg_a", obj, 0, 1.0f);
            p.bg = ImVec4(bgx, bgy, bgz, bga);
        }

        // Items
        size_t items_pos = obj.find("\"items\"");
        if (items_pos != std::string::npos) {
            size_t arr_start = obj.find("[", items_pos);
            if (arr_start != std::string::npos) {
                size_t arr_end = std::string::npos;
                int arr_depth = 0;
                for (size_t i = arr_start; i < obj.size(); ++i) {
                    if (obj[i] == '[') arr_depth++;
                    else if (obj[i] == ']') {
                        arr_depth--;
                        if (arr_depth == 0) {
                            arr_end = i;
                            break;
                        }
                    }
                }
                size_t cursor = arr_start + 1;
                while (arr_end != std::string::npos && cursor < arr_end) {
                    size_t it_start = obj.find("{", cursor);
                    if (it_start == std::string::npos || it_start > arr_end) break;
                    int it_depth = 0;
                    size_t it_end = std::string::npos;
                    for (size_t i = it_start; i <= arr_end; ++i) {
                        if (obj[i] == '{') it_depth++;
                        else if (obj[i] == '}') {
                            it_depth--;
                            if (it_depth == 0) {
                                it_end = i;
                                break;
                            }
                        }
                    }
                    if (it_end == std::string::npos) break;
                    std::string item_str = obj.substr(it_start, it_end - it_start + 1);
                    ComposerItem it{};
                    it.id = read_int("id", item_str, 0, max_item_id++);
                    max_item_id = std::max(max_item_id, it.id + 1);
                    std::string tname = read_string("type", item_str, 0);
                    if (!tname.empty()) {
                        std::string lower = tname;
                        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return (char)std::tolower(c); });
                        if (lower.find("field") != std::string::npos) it.type = COMPOSER_FIELD_VIEW;
                        else if (lower.find("legend") != std::string::npos) it.type = COMPOSER_LEGEND;
                        else if (lower.find("scope") != std::string::npos) it.type = COMPOSER_SCOPE;
                        else if (lower.find("fft") != std::string::npos) it.type = COMPOSER_FFT;
                        else if (lower.find("region") != std::string::npos) it.type = COMPOSER_REGION;
                        else if (lower.find("meas") != std::string::npos) it.type = COMPOSER_MEAS;
                        else if (lower.find("smith") != std::string::npos) it.type = COMPOSER_SMITH;
                    } else {
                        it.type = (ComposerItemType)read_int("type", item_str, 0, (int)COMPOSER_FIELD_VIEW);
                    }
                    if (it.type < COMPOSER_FIELD_VIEW || it.type > COMPOSER_SMITH) {
                        it.type = COMPOSER_FIELD_VIEW;
                    }
                    it.viewport_idx = read_int("viewport_idx", item_str, 0, 0);
                    it.viewport_idx = read_int("viewport", item_str, 0, it.viewport_idx);
                    it.viewport_idx = read_int("vp", item_str, 0, it.viewport_idx);
                    it.pos = read_vec2("pos", item_str, 0, ImVec2(80, 80));
                    it.size = read_vec2("size", item_str, 0, ImVec2(320, 240));
                    if (it.size.x <= 0.0f) it.size.x = read_float("w", item_str, 0, 320.0f);
                    if (it.size.y <= 0.0f) it.size.y = read_float("h", item_str, 0, 240.0f);
                    ImVec4 reg = read_vec4("region", item_str, 0, ImVec4(0, 0, 1, 1));
                    ImVec4 reg_norm = read_vec4("region_norm", item_str, 0, reg);
                    it.region_norm = reg_norm;
                    it.selected = false;
                    p.items.push_back(it);
                    cursor = it_end + 1;
                }
            }
        }

        loaded.push_back(p);
        pos = obj_end + 1;
        size_t next_obj = data.find("{", pos);
        size_t closing = data.find("]", pos);
        if (closing != std::string::npos && (next_obj == std::string::npos || closing < next_obj)) break;
    }
    if (!loaded.empty()) {
        app->composer_pages = loaded;
        app->composer_active_page = 0;
        app->composer_next_item_id = max_item_id;
        return true;
    }
    return false;
}

static std::string join_args(const std::vector<std::string>& args) {
    std::string out;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) out.push_back(' ');
        out += args[i];
    }
    return out;
}

static int run_ffmpeg(const std::vector<std::string>& args) {
    if (args.empty()) return -1;
#ifdef _WIN32
    std::vector<const char*> cargs;
    cargs.reserve(args.size() + 1);
    for (const auto& s : args) cargs.push_back(s.c_str());
    cargs.push_back(nullptr);
    return _spawnv(_P_WAIT, args[0].c_str(), cargs.data());
#else
    std::string cmd = join_args(args);
    return std::system(cmd.c_str());
#endif
}

static void log_ffmpeg_attempt(const char* cmd, int rc) {
    if (!cmd) return;
    std::error_code ec;
    std::filesystem::create_directories("recordings", ec);
    std::ofstream f("recordings/ffmpeg_last.txt", std::ios::out | std::ios::trunc);
    if (!f) return;
    f << "rc=" << rc << "\n";
    f << "cmd:\n" << cmd << "\n";
}

static void free_source_expression(Source* s) {
    if (!s || !s->expr_program) return;
    expr_free((ExprProgram*)s->expr_program);
    s->expr_program = NULL;
}

static void apply_source_spec_to_runtime(SimulationState* sim,
                                         int idx,
                                         const SourceConfigSpec* spec) {
    if (!sim || !spec) return;
    if (idx < 0 || idx >= MAX_SRC) return;
    Source* dst = &sim->sources[idx];
    dst->active = spec->active;
    dst->type = spec->type;
    dst->field = spec->field;
    dst->amp = spec->amp;
    dst->freq = spec->freq;
    dst->sigma2 = (spec->sigma2 > 0.01) ? spec->sigma2 : 0.01;
    dst->ix = normalized_to_cell_index(spec->x, sim->nx);
    dst->iy = normalized_to_cell_index(spec->y, sim->ny);
    std::strncpy(dst->expr_text, spec->expr, SOURCE_EXPR_MAX_LEN);
    dst->expr_text[SOURCE_EXPR_MAX_LEN - 1] = '\0';
    free_source_expression(dst);
    source_reparam(dst);
    if (dst->type == SRC_EXPR && dst->expr_text[0] != '\0') {
        ExprProgram* prog = NULL;
        char errbuf[128];
        if (expr_compile(dst->expr_text, &prog, errbuf, (int)sizeof(errbuf))) {
            dst->expr_program = prog;
        } else {
            dst->expr_program = NULL;
        }
    }
}

static void clear_runtime_source(SimulationState* sim, int idx) {
    if (!sim || idx < 0 || idx >= MAX_SRC) return;
    Source* dst = &sim->sources[idx];
    dst->active = 0;
    free_source_expression(dst);
    dst->expr_text[0] = '\0';
}

static void create_new_source(WizardState* wizard, SimulationState* sim, AppState* app) {
    if (!wizard || !sim) return;
    int count = wizard->cfg.source_count;
    if (count >= MAX_SRC) return;
    SourceConfigSpec* spec = &wizard->cfg.source_configs[count];
    spec->active = 1;
    spec->type = SRC_CW;
    spec->field = SRC_FIELD_EZ;
    spec->x = 0.5;
    spec->y = 0.5;
    spec->amp = 1.0;
    spec->freq = 1.0e9;
    spec->sigma2 = 4.0;
    spec->expr[0] = '\0';
    wizard->cfg.source_count = count + 1;
    apply_source_spec_to_runtime(sim, count, spec);
    if (app) {
        app->selected_source = count;
        ui_log_add(app, "Added source #%d", count);
    }
}

static void create_new_source_at(WizardState* wizard,
                                 SimulationState* sim,
                                 AppState* app,
                                 int ix,
                                 int iy) {
    if (!wizard || !sim) return;
    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) return;
    int count = wizard->cfg.source_count;
    if (count >= MAX_SRC) return;
    SourceConfigSpec* spec = &wizard->cfg.source_configs[count];
    spec->active = 1;
    spec->type = SRC_CW;
    spec->field = SRC_FIELD_EZ;
    spec->amp = 1.0;
    spec->freq = sim->freq > 0.0 ? sim->freq : 1.0e9;
    spec->sigma2 = 4.0;
    spec->expr[0] = '\0';

    double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
    double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
    spec->x = (double)ix / nx1;
    spec->y = (double)iy / ny1;

    wizard->cfg.source_count = count + 1;
    apply_source_spec_to_runtime(sim, count, spec);
    if (app) {
        app->selected_source = count;
        ui_log_add(app, "Added source #%d at (%d, %d)", count, ix, iy);
    }
}

static void create_block_at(WizardState* wizard,
                            SimulationBootstrap* bootstrap,
                            SimulationState* sim,
                            AppState* app,
                            int ix,
                            int iy) {
    if (!wizard || !sim || !bootstrap) return;
    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) return;
    if (wizard->cfg.material_rect_count >= CONFIG_MAX_MATERIAL_RECTS) return;

    int idx = wizard->cfg.material_rect_count;
    MaterialRectSpec& r = wizard->cfg.material_rects[idx];
    double nx_d = (double)sim->nx;
    double ny_d = (double)sim->ny;
    r.x0 = (double)ix / nx_d;
    r.x1 = (double)(ix + 1) / nx_d;
    r.y0 = (double)iy / ny_d;
    r.y1 = (double)(iy + 1) / ny_d;
    r.epsr = 4.0;
    r.sigma = 0.0;
    r.tag = 0;
    wizard->cfg.material_rect_count = idx + 1;
    apply_wizard_materials_to_sim(*wizard, bootstrap, sim);
    if (app) {
        app->selected_block = idx;
        ui_log_add(app, "Added block at (%d, %d)", ix, iy);
    }
}

static const Material* guess_material_for_rect(const MaterialRectSpec& r) {
    const Material* best = nullptr;
    double best_score = 1e9;
    int mat_count = material_library_get_count();
    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;
        if (r.tag == 1) {
            if (mat->type == MAT_TYPE_PEC) return mat;
            continue;
        }
        if (r.tag == 2) {
            if (mat->type == MAT_TYPE_PMC) return mat;
            continue;
        }
        if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) continue;
        double diff = std::fabs(mat->epsilon_r - r.epsr);
        if (diff < best_score) {
            best_score = diff;
            best = mat;
        }
    }
    return best;
}

static void apply_material_to_rect(MaterialRectSpec* rect, const Material* mat) {
    if (!rect || !mat) return;
    if (mat->type == MAT_TYPE_PEC) {
        rect->tag = 1;
        rect->epsr = 1.0;
        rect->sigma = mat->conductivity;
    } else if (mat->type == MAT_TYPE_PMC) {
        rect->tag = 2;
        rect->epsr = 1.0;
        rect->sigma = 0.0;
    } else {
        rect->tag = 0;
        rect->epsr = mat->epsilon_r;
        rect->sigma = mat->conductivity;
    }
}

static void remove_block(WizardState* wizard,
                         SimulationBootstrap* bootstrap,
                         SimulationState* sim,
                         AppState* app,
                         int idx) {
    if (!wizard || idx < 0 || idx >= wizard->cfg.material_rect_count) return;
    for (int j = idx; j < wizard->cfg.material_rect_count - 1; ++j) {
        wizard->cfg.material_rects[j] = wizard->cfg.material_rects[j + 1];
    }
    wizard->cfg.material_rect_count--;
    if (bootstrap && sim) {
        apply_wizard_materials_to_sim(*wizard, bootstrap, sim);
    }
    if (app) {
        ui_log_add(app, "Deleted block #%d", idx);
    }
}

static void delete_source_at(WizardState* wizard, SimulationState* sim, AppState* app, int idx) {
    if (!wizard || !sim) return;
    int count = wizard->cfg.source_count;
    if (idx < 0 || idx >= count) return;
    for (int i = idx; i < count - 1; ++i) {
        wizard->cfg.source_configs[i] = wizard->cfg.source_configs[i + 1];
        apply_source_spec_to_runtime(sim, i, &wizard->cfg.source_configs[i]);
    }
    wizard->cfg.source_count = count - 1;
    clear_runtime_source(sim, wizard->cfg.source_count);
    if (app) {
        if (app->selected_source >= wizard->cfg.source_count) {
            app->selected_source = wizard->cfg.source_count - 1;
        }
        ui_log_add(app, "Deleted source #%d", idx);
    }
}

/* Scene overview panel ---------------------------------------------------- */
static void draw_scene_panel(const SimulationState* sim, const WizardState& wizard) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Scene")) return;

    ImGui::Indent();
    if (sim) {
        ImGui::Text("Grid: %d x %d", sim->nx, sim->ny);
        ImGui::Text("Domain: %.3f m x %.3f m", sim->lx, sim->ly);
        ImGui::Text("dt: %.3e s", sim->dt);
        ImGui::Separator();
    } else {
        ImGui::TextUnformatted("Simulation not available");
        ImGui::Unindent();
        return;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Materials");
    ImGui::Indent();
    ImGui::Text("Rectangles: %d", wizard.cfg.material_rect_count);
    ImGui::Unindent();

    ImGui::Separator();
    ImGui::TextUnformatted("Sources");
    ImGui::Indent();
    ImGui::Text("Count: %d (max %d)", wizard.cfg.source_count, MAX_SRC);
    ImGui::Unindent();
}

static void draw_grid_panel(WizardState* wizard, AppState* app) {
    if (!wizard) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Grid & Domain")) return;

    ImGui::Indent();
    int nx = wizard->cfg.nx;
    if (ImGui::InputInt("Cells X (nx)", &nx)) {
        if (nx < 8) nx = 8;
        wizard->cfg.nx = nx;
    }
    int ny = wizard->cfg.ny;
    if (ImGui::InputInt("Cells Y (ny)", &ny)) {
        if (ny < 8) ny = 8;
        wizard->cfg.ny = ny;
    }

    double lx = wizard->cfg.lx;
    if (ImGui::InputDouble("Width (meters)", &lx, 0.01, 0.1, "%.3f")) {
        if (lx < 1e-3) lx = 1e-3;
        wizard->cfg.lx = lx;
    }
    double ly = wizard->cfg.ly;
    if (ImGui::InputDouble("Height (meters)", &ly, 0.01, 0.1, "%.3f")) {
        if (ly < 1e-3) ly = 1e-3;
        wizard->cfg.ly = ly;
    }

    float cfl = (float)wizard->cfg.cfl_safety;
    if (ImGui::SliderFloat("CFL Safety", &cfl, 0.1f, 0.99f, "%.2f")) {
        wizard->cfg.cfl_safety = (double)cfl;
    }

    const char* boundary_items[] = {"CPML", "Mur"};
    int boundary_idx = (wizard->cfg.boundary_mode == SIM_BOUNDARY_MUR) ? 1 : 0;
    if (ImGui::Combo("Boundary", &boundary_idx, boundary_items, IM_ARRAYSIZE(boundary_items))) {
        wizard->cfg.boundary_mode = (boundary_idx == 1) ? SIM_BOUNDARY_MUR : SIM_BOUNDARY_CPML;
    }

    ImGui::Separator();
    if (ImGui::Button("Apply Grid Changes & Restart")) {
        if (app) {
            app->request_rebootstrap = true;
            std::snprintf(app->rebootstrap_message,
                          sizeof(app->rebootstrap_message),
                          "Rebooted after grid/domain update");
        }
    }
    ImGui::Unindent();
}

/* Sources panel ----------------------------------------------------------- */
static void draw_sources_panel(SimulationState* sim, WizardState& wizard, AppState* app) {
    if (!sim || !app) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Sources")) return;

    ImGui::Indent();
    int count = wizard.cfg.source_count;
    if (count < 0) count = 0;
    if (count > MAX_SRC) count = MAX_SRC;
    wizard.cfg.source_count = count;

    if (count == 0) {
        app->selected_source = -1;
    } else {
        if (app->selected_source < 0) app->selected_source = 0;
        if (app->selected_source >= count) app->selected_source = count - 1;
    }

    ImGui::Text("Configured: %d / %d", count, MAX_SRC);
    ImGui::SameLine();
    bool can_add = (count < MAX_SRC);
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ New Source")) {
        create_new_source(&wizard, sim, app);
        count = wizard.cfg.source_count;
    }
    if (!can_add) {
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::TextDisabled("Max reached");
    }
    ImGui::Separator();

    if (count == 0) {
        ImGui::TextDisabled("No sources configured.");
    } else {
        ImGui::TextUnformatted("Source list");
        ImGui::Indent();
        for (int i = 0; i < count; ++i) {
            SourceConfigSpec& spec = wizard.cfg.source_configs[i];
            ImGui::PushID(i);
            char label[64];
            std::snprintf(label,
                          sizeof(label),
                          "#%d %s",
                          i,
                          source_type_label(spec.type));
            bool selected = (app->selected_source == i);
            if (ImGui::Selectable(label, selected)) {
                app->selected_source = i;
            }
            ImGui::SameLine();
            bool active = spec.active != 0;
            if (ImGui::Checkbox("##active", &active)) {
                spec.active = active ? 1 : 0;
                apply_source_spec_to_runtime(sim, i, &spec);
            }
            ImGui::SameLine();
            ImGui::Text("%s", source_field_label(spec.field));
            ImGui::SameLine();
            ImGui::TextDisabled("Pos: (%.2f, %.2f)", spec.x, spec.y);
            ImGui::PopID();
        }
        ImGui::Unindent();
    }

    int idx = app->selected_source;
    if (idx >= 0 && idx < count) {
        SourceConfigSpec& spec = wizard.cfg.source_configs[idx];
        ImGui::Separator();
        ImGui::Text("Source #%d", idx);
        ImGui::Indent();

        const char* types[] = {"CW", "Gaussian", "Ricker", "Expr"};
        int type_idx = (int)spec.type;
        if (type_idx < 0 || type_idx >= (int)IM_ARRAYSIZE(types)) type_idx = 0;
        if (ImGui::Combo("Type", &type_idx, types, IM_ARRAYSIZE(types))) {
            spec.type = (SourceType)type_idx;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }

        const char* field_items[] = {"Ez", "Hx", "Hy"};
        int field_idx = (int)spec.field;
        if (field_idx < 0 || field_idx >= (int)IM_ARRAYSIZE(field_items)) field_idx = 0;
        if (ImGui::Combo("Field", &field_idx, field_items, IM_ARRAYSIZE(field_items))) {
            spec.field = (SourceFieldType)field_idx;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }

        float norm_x = (float)spec.x;
        if (ImGui::SliderFloat("X (0-1)", &norm_x, 0.0f, 1.0f, "%.3f")) {
            spec.x = clamp01(norm_x);
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        float norm_y = (float)spec.y;
        if (ImGui::SliderFloat("Y (0-1)", &norm_y, 0.0f, 1.0f, "%.3f")) {
            spec.y = clamp01(norm_y);
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        const Source& runtime_src = sim->sources[idx];
        ImGui::TextDisabled("Grid cell: (%d, %d)", runtime_src.ix, runtime_src.iy);

        double amp = spec.amp;
        if (ImGui::InputDouble("Amplitude", &amp, 0.1, 1.0, "%.3f")) {
            spec.amp = amp;
            apply_source_spec_to_runtime(sim, idx, &spec);
        }
        double freq = spec.freq;
        if (ImGui::InputDouble("Frequency (Hz)", &freq, 1e6, 1e8, "%.3e")) {
            if (freq > 0.0) {
                spec.freq = freq;
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
        }

        if (spec.type == SRC_GAUSS_PULSE || spec.type == SRC_RICKER) {
            double sigma = spec.sigma2;
            if (ImGui::InputDouble("Sigma^2 (cells^2)", &sigma, 0.1, 1.0, "%.3f")) {
                if (sigma < 0.1) sigma = 0.1;
                spec.sigma2 = sigma;
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
        } else {
            ImGui::TextDisabled("Sigma^2 applies to Gaussian/Ricker pulses.");
        }

        if (spec.type == SRC_EXPR) {
            if (ImGui::InputTextMultiline("Expression",
                                          spec.expr,
                                          SOURCE_EXPR_MAX_LEN,
                                          ImVec2(-1.0f, ImGui::GetTextLineHeight() * 4))) {
                apply_source_spec_to_runtime(sim, idx, &spec);
            }
            ImGui::TextDisabled("Variables: t (s), amp, freq, pi");

            if (spec.expr[0] != '\0') {
                char errbuf[128];
                ExprProgram* prog = nullptr;
                if (expr_compile(spec.expr, &prog, errbuf, sizeof(errbuf))) {
                    const int N = 128;
                    float values[N];
                    double period = (spec.freq > 0.0) ? (1.0 / spec.freq) : 1e-9;
                    double t_max = (period > 0.0) ? 3.0 * period : 3e-9;
                    double vmin = 0.0;
                    double vmax = 0.0;
                    bool first = true;
                    for (int k = 0; k < N; ++k) {
                        double t = t_max * (double)k / (double)(N - 1);
                        double v = expr_eval(prog, t, spec.amp, spec.freq);
                        values[k] = (float)v;
                        if (first) {
                            vmin = vmax = v;
                            first = false;
                        } else {
                            if (v < vmin) vmin = v;
                            if (v > vmax) vmax = v;
                        }
                    }
                    float ymin = (float)vmin;
                    float ymax = (float)vmax;
                    if (ymin == ymax) {
                        ymin -= 1.0f;
                        ymax += 1.0f;
                    }
                    ImGui::PlotLines("Preview", values, N, 0, nullptr, ymin, ymax,
                                     ImVec2(-1.0f, ImGui::GetTextLineHeight() * 5));
                    expr_free(prog);
                } else {
                    ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.5f, 1.0f),
                                       "Expression error: %s",
                                       errbuf);
                }
            }
        }

        ImGui::Separator();
        if (ImGui::Button("Delete Source")) {
            ImGui::OpenPopup("delete_source_popup");
        }
        if (ImGui::BeginPopupModal("delete_source_popup", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("Delete source #%d?", idx);
            if (ImGui::Button("Delete", ImVec2(120, 0))) {
                delete_source_at(&wizard, sim, app, idx);
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
        ImGui::Unindent();
    }
    ImGui::Unindent();
}

static void draw_run_settings_panel(WizardState* wizard, AppState* app) {
    if (!wizard) return;
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Run Settings")) return;

    ImGui::Indent();
    const char* run_modes[] = {"Fixed steps", "Sweep"};
    int run_mode_idx = (wizard->cfg.run_mode == SIM_RUN_MODE_SWEEP) ? 1 : 0;
    if (ImGui::Combo("Run Mode", &run_mode_idx, run_modes, IM_ARRAYSIZE(run_modes))) {
        wizard->cfg.run_mode = (run_mode_idx == 1) ? SIM_RUN_MODE_SWEEP : SIM_RUN_MODE_FIXED_STEPS;
    }

    if (wizard->cfg.run_mode == SIM_RUN_MODE_FIXED_STEPS) {
        ImGui::InputInt("Run steps", &wizard->cfg.run_steps);
        if (wizard->cfg.run_steps < 1) wizard->cfg.run_steps = 1;
    } else {
        ImGui::InputInt("Sweep points", &wizard->cfg.sweep_points);
        if (wizard->cfg.sweep_points < 1) wizard->cfg.sweep_points = 1;
        ImGui::InputDouble("Sweep start (Hz)", &wizard->cfg.sweep_start_hz, 1e6, 1e8, "%.3e");
        if (wizard->cfg.sweep_start_hz < 1.0) wizard->cfg.sweep_start_hz = 1.0;
        ImGui::InputDouble("Sweep stop (Hz)", &wizard->cfg.sweep_stop_hz, 1e6, 1e8, "%.3e");
        if (wizard->cfg.sweep_stop_hz < wizard->cfg.sweep_start_hz) {
            wizard->cfg.sweep_stop_hz = wizard->cfg.sweep_start_hz;
        }
        ImGui::InputInt("Steps / point", &wizard->cfg.sweep_steps_per_point);
        if (wizard->cfg.sweep_steps_per_point < 1) wizard->cfg.sweep_steps_per_point = 1;
    }

    if (wizard->advanced) {
        ImGui::Separator();
        bool profile = wizard->cfg.enable_profile != 0;
        if (ImGui::Checkbox("Enable profiling", &profile)) {
            wizard->cfg.enable_profile = profile ? 1 : 0;
        }
        bool probe_log = wizard->cfg.enable_probe_log != 0;
        if (ImGui::Checkbox("Enable probe log", &probe_log)) {
            wizard->cfg.enable_probe_log = probe_log ? 1 : 0;
        }
        ImGui::InputText("Probe log path",
                         wizard->cfg.probe_log_path,
                         SIM_PROBE_LOG_PATH_MAX);
    } else {
        ImGui::Separator();
        ImGui::TextDisabled("Switch to Advanced mode for profiling/logging options.");
    }

    ImGui::Separator();
    if (ImGui::Button("Apply Run Settings & Restart")) {
        if (app) {
            app->request_rebootstrap = true;
            std::snprintf(app->rebootstrap_message,
                          sizeof(app->rebootstrap_message),
                          "Rebooted after run settings update");
        }
    }
    ImGui::Unindent();
}

/* Blocks panel ------------------------------------------------------------ */
static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app);

typedef struct {
    int material_id;
    const char* material_name;
    int cell_count;
    float percentage;
} MaterialStats;

static int get_material_statistics(const SimulationState* sim, MaterialStats* stats, int max_stats) {
    (void)sim;
    (void)stats;
    (void)max_stats;
    /* Placeholder until simulation tracks material IDs explicitly. */
    return 0;
}

static int compute_scope_fft(const Scope* scope,
                             double dt,
                             double* freq,
                             double* mag,
                             double* phase,
                             int max_fft) {
    if (!scope || !scope->y || scope->n <= 8 || dt <= 0.0) return 0;
    if (!freq || !mag || !phase || max_fft <= 0) return 0;

    int Nfft = scope->n;
    if (Nfft > max_fft) Nfft = max_fft;
    if (Nfft < 64) return 0;

    std::vector<double> x(Nfft);
    int start = (scope->head - Nfft + scope->n) % scope->n;
    for (int k = 0; k < Nfft; ++k) {
        int idx = (start + k) % scope->n;
        x[k] = scope->y[idx];
    }

    double mean = 0.0;
    for (int k = 0; k < Nfft; ++k) mean += x[k];
    mean /= (double)Nfft;

    const double two_pi = 2.0 * 3.14159265358979323846;
    if (Nfft > 1) {
        for (int k = 0; k < Nfft; ++k) {
            double w = 0.5 * (1.0 - std::cos(two_pi * (double)k / (double)(Nfft - 1)));
            x[k] = (x[k] - mean) * w;
        }
    } else {
        for (int k = 0; k < Nfft; ++k) {
            x[k] = (x[k] - mean);
        }
    }

    int half = Nfft / 2;
    if (half <= 1) return 0;

    for (int k = 0; k < half; ++k) {
        double real = 0.0;
        double imag = 0.0;
        for (int n = 0; n < Nfft; ++n) {
            double angle = -two_pi * (double)k * (double)n / (double)Nfft;
            double c = std::cos(angle);
            double s = std::sin(angle);
            real += x[n] * c;
            imag += x[n] * s;
        }
        real /= (double)Nfft;
        imag /= (double)Nfft;
        mag[k] = std::sqrt(real * real + imag * imag);
        phase[k] = std::atan2(imag, real) * (180.0 / 3.14159265358979323846);
        freq[k] = (double)k / ((double)Nfft * dt);
    }

    return half;
}

static void clear_sparameter_data(AppState* app) {
    if (!app) return;
    app->sparam_freq_hz.clear();
    app->sparam_freq_ghz.clear();
    app->sparam_s21_mag.clear();
    app->sparam_s21_db.clear();
    app->sparam_s11_mag.clear();
    app->sparam_s11_db.clear();
    app->sparam_vswr.clear();
    app->sparam_peak_db = 0.0;
    app->sparam_peak_freq_hz = 0.0;
    app->sparam_f_low_hz = 0.0;
    app->sparam_f_high_hz = 0.0;
    app->sparam_f_center_hz = 0.0;
    app->sparam_bandwidth_hz = 0.0;
    app->sparam_data_loaded = false;
}

static void analyze_sparameter_bandwidth(AppState* app) {
    if (!app || app->sparam_freq_hz.empty()) return;
    app->sparam_peak_db = -1e9;
    app->sparam_peak_freq_hz = 0.0;
    const size_t count = app->sparam_freq_hz.size();
    for (size_t i = 0; i < count; ++i) {
        if (app->sparam_s21_db[i] > app->sparam_peak_db) {
            app->sparam_peak_db = app->sparam_s21_db[i];
            app->sparam_peak_freq_hz = app->sparam_freq_hz[i];
        }
    }
    double threshold = app->sparam_peak_db - 3.0;
    app->sparam_f_low_hz = 0.0;
    app->sparam_f_high_hz = 0.0;
    for (size_t i = 0; i < count; ++i) {
        if (app->sparam_s21_db[i] >= threshold) {
            if (app->sparam_f_low_hz == 0.0) {
                app->sparam_f_low_hz = app->sparam_freq_hz[i];
            }
            app->sparam_f_high_hz = app->sparam_freq_hz[i];
        }
    }
    if (app->sparam_f_high_hz > app->sparam_f_low_hz && app->sparam_f_low_hz > 0.0) {
        app->sparam_f_center_hz = 0.5 * (app->sparam_f_high_hz + app->sparam_f_low_hz);
        app->sparam_bandwidth_hz = app->sparam_f_high_hz - app->sparam_f_low_hz;
    } else {
        app->sparam_f_center_hz = 0.0;
        app->sparam_bandwidth_hz = 0.0;
    }
}

static bool load_sparameter_csv(AppState* app, const char* path) {
    if (!app || !path) return false;
    std::ifstream file(path);
    if (!file.is_open()) {
        return false;
    }
    clear_sparameter_data(app);
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line.find("freq") != std::string::npos) continue;
        std::stringstream ss(line);
        std::string freq_str;
        std::string s21_str;
        if (!std::getline(ss, freq_str, ',')) continue;
        if (!std::getline(ss, s21_str, ',')) continue;
        double freq = std::strtod(freq_str.c_str(), nullptr);
        double s21 = std::strtod(s21_str.c_str(), nullptr);
        if (!std::isfinite(freq) || !std::isfinite(s21)) continue;
        if (freq <= 0.0) continue;
        if (s21 < 0.0) s21 = 0.0;
        app->sparam_freq_hz.push_back(freq);
        app->sparam_s21_mag.push_back(s21);
    }
    if (app->sparam_freq_hz.empty()) {
        app->sparam_data_loaded = false;
        return false;
    }
    const double eps = 1e-12;
    size_t count = app->sparam_freq_hz.size();
    app->sparam_freq_ghz.resize(count);
    app->sparam_s21_db.resize(count);
    app->sparam_vswr.resize(count);
    app->sparam_s11_mag.resize(count);
    app->sparam_s11_db.resize(count);
    for (size_t i = 0; i < count; ++i) {
        double mag = app->sparam_s21_mag[i];
        double mag_clamped = std::min(std::max(mag, 0.0), 0.999999);
        app->sparam_freq_ghz[i] = app->sparam_freq_hz[i] * 1e-9;
        double db = 20.0 * std::log10(mag + eps);
        app->sparam_s21_db[i] = db;
        double s11 = std::sqrt(std::max(0.0, 1.0 - (mag_clamped * mag_clamped)));
        app->sparam_s11_mag[i] = s11;
        app->sparam_s11_db[i] = 20.0 * std::log10(s11 + eps);
        double vswr = (1.0 + s11) / std::max(1e-6, 1.0 - s11);
        app->sparam_vswr[i] = vswr;
    }
    analyze_sparameter_bandwidth(app);
    app->sparam_data_loaded = true;
    std::strncpy(app->sparam_csv_path, path, sizeof(app->sparam_csv_path));
    app->sparam_csv_path[sizeof(app->sparam_csv_path) - 1] = '\0';
    if (app->smith_freq_index >= (int)count) {
        app->smith_freq_index = (int)count - 1;
        if (app->smith_freq_index < 0) app->smith_freq_index = 0;
    }
    return true;
}

static bool string_contains_case_insensitive(const char* haystack, const char* needle) {
    if (!haystack || !needle) return false;
    if (needle[0] == '\0') return true;
    for (const char* h = haystack; *h; ++h) {
        const char* h_it = h;
        const char* n_it = needle;
        while (*h_it && *n_it) {
            char hc = (char)std::tolower((unsigned char)*h_it);
            char nc = (char)std::tolower((unsigned char)*n_it);
            if (hc != nc) break;
            ++h_it;
            ++n_it;
        }
        if (*n_it == '\0') {
            return true;
        }
    }
    return false;
}

static void draw_material_browser(AppState* app, int* active_material_id) {
    if (!app || !app->material_browser_open) return;

    ImGui::SetNextWindowSize(ImVec2(400.0f, 500.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Material Browser", &app->material_browser_open)) {
        ImGui::End();
        return;
    }

    ImGui::TextUnformatted("Search:");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputText("##MaterialSearch", app->material_search, sizeof(app->material_search));

    ImGui::Separator();
    ImGui::Checkbox("Metals", &app->filter_metals);
    ImGui::SameLine();
    ImGui::Checkbox("Dielectrics", &app->filter_dielectrics);

    ImGui::Separator();

    ImGui::BeginChild("MaterialList", ImVec2(0, -170.0f), true);

    int mat_count = material_library_get_count();
    const bool has_query = (app->material_search[0] != '\0');

    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;

        if (mat->category == MAT_CAT_METAL && !app->filter_metals) continue;
        if (mat->category == MAT_CAT_DIELECTRIC && !app->filter_dielectrics) continue;

        if (has_query && !string_contains_case_insensitive(mat->name, app->material_search)) {
            continue;
        }

        ImVec4 color = ImVec4(
            mat->color_r / 255.0f,
            mat->color_g / 255.0f,
            mat->color_b / 255.0f,
            1.0f
        );

        ImGui::PushStyleColor(ImGuiCol_Text, color);
        ImGui::Bullet();
        ImGui::PopStyleColor();
        ImGui::SameLine();

        bool is_selected = (app->selected_material_id == mat->id);
        if (ImGui::Selectable(mat->name, is_selected, ImGuiSelectableFlags_SpanAvailWidth)) {
            app->selected_material_id = mat->id;
        }

        char prop_buf[64];
        if (mat->type == MAT_TYPE_PEC) {
            std::snprintf(prop_buf, sizeof(prop_buf), "sigma=%.2e S/m", mat->conductivity);
        } else {
            std::snprintf(prop_buf, sizeof(prop_buf), "eps=%.2f", mat->epsilon_r);
        }
        ImGui::SameLine(250.0f);
        ImGui::TextDisabled("%s", prop_buf);
    }

    ImGui::EndChild();

    ImGui::Separator();
    ImGui::TextUnformatted("Selected Material:");

    if (app->selected_material_id >= 0) {
        const Material* sel = material_library_get_by_id(app->selected_material_id);
        if (sel) {
            ImGui::Indent();
            ImVec4 color = ImVec4(
                sel->color_r / 255.0f,
                sel->color_g / 255.0f,
                sel->color_b / 255.0f,
                1.0f
            );
            ImGui::PushStyleColor(ImGuiCol_Text, color);
            ImGui::TextUnformatted(sel->name);
            ImGui::PopStyleColor();

            const char* category_label = (sel->category == MAT_CAT_METAL) ? "Metal" : "Dielectric";
            ImGui::TextDisabled("Category: %s", category_label);

            if (sel->type == MAT_TYPE_PEC) {
                ImGui::Text("sigma = %.2e S/m", sel->conductivity);
            } else {
                ImGui::Text("eps = %.2f", sel->epsilon_r);
                if (sel->tan_delta > 0.0) {
                    ImGui::Text("tan(delta) = %.4f", sel->tan_delta);
                }
            }

            ImGui::TextWrapped("%s", sel->description ? sel->description : "");
            ImGui::Unindent();

            if (ImGui::Button("Use for Painting", ImVec2(-1.0f, 0.0f))) {
                if (active_material_id) {
                    *active_material_id = app->selected_material_id;
                }
                app->paint_material_id = app->selected_material_id;
                app->material_browser_open = false;
            }
        }
    } else {
        ImGui::TextDisabled("No material selected");
    }

    ImGui::End();
}

static bool rect_matches_material(const MaterialRectSpec& r, const Material* mat) {
    if (!mat) return false;
    if (r.tag == 1) return mat->type == MAT_TYPE_PEC;
    if (r.tag == 2) return mat->type == MAT_TYPE_PMC;
    if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) return false;
    return std::fabs(r.epsr - mat->epsilon_r) < 0.25;
}

static void draw_material_legend(AppState* app,
                                 SimulationState* sim,
                                 WizardState* wizard,
                                 SimulationBootstrap* bootstrap) {
    (void)sim;
    if (!app || !wizard || !bootstrap) return;
    if (!app->show_material_legend) return;

    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Material Legend")) return;

    ImGui::TextUnformatted("Materials in Scene:");
    ImGui::Separator();
    ImGui::BeginChild("LegendList", ImVec2(0.0f, 160.0f), true);

    int mat_count = material_library_get_count();
    for (int i = 0; i < mat_count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;

        ImVec4 color = ImVec4(mat->color_r / 255.0f,
                              mat->color_g / 255.0f,
                              mat->color_b / 255.0f,
                              1.0f);
        ImDrawList* dl = ImGui::GetWindowDrawList();
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImVec2 q = ImVec2(p.x + 16.0f, p.y + 16.0f);
        dl->AddRectFilled(p, q, ImGui::ColorConvertFloat4ToU32(color));
        dl->AddRect(p, q, IM_COL32(90, 90, 90, 255));
        ImGui::Dummy(ImVec2(16.0f, 16.0f));
        ImGui::SameLine();

        bool row_selected = (app->selected_material_id == mat->id);
        if (ImGui::Selectable(mat->name, row_selected, ImGuiSelectableFlags_SpanAvailWidth)) {
            app->selected_material_id = mat->id;
            if (app->auto_filter_blocks_on_select) {
                app->filter_blocks_by_material = true;
            }
            for (int b = 0; b < wizard->cfg.material_rect_count; ++b) {
                if (rect_matches_material(wizard->cfg.material_rects[b], mat)) {
                    app->selected_block = b;
                    break;
                }
            }
        }

        ImGui::SameLine(170.0f);
        if (mat->type == MAT_TYPE_PEC) {
            ImGui::TextDisabled("PEC");
        } else if (mat->type == MAT_TYPE_PMC) {
            ImGui::TextDisabled("PMC");
        } else {
            ImGui::TextDisabled("εᵣ=%.1f", mat->epsilon_r);
        }

        ImGui::SameLine();
        bool can_add = (wizard->cfg.material_rect_count < CONFIG_MAX_MATERIAL_RECTS);
        if (!can_add) ImGui::BeginDisabled();
        char btn_label[32];
        std::snprintf(btn_label, sizeof(btn_label), "Add##%d", mat->id);
        if (ImGui::SmallButton(btn_label)) {
            int idx = wizard->cfg.material_rect_count;
            MaterialRectSpec& r = wizard->cfg.material_rects[idx];
            r.x0 = 0.3;
            r.y0 = 0.3;
            r.x1 = 0.7;
            r.y1 = 0.7;
            apply_material_to_rect(&r, mat);
            wizard->cfg.material_rect_count = idx + 1;
            apply_wizard_materials_to_sim(*wizard, bootstrap, sim);
            app->selected_block = idx;
            ui_log_add(app, "Added block #%d from %s", idx, mat->name);
        }
        if (!can_add) ImGui::EndDisabled();
    }

    ImGui::EndChild();
    ImGui::TextDisabled("Legend uses library colors (usage stats coming soon)");
}

static void draw_resistance_circle(double r_norm) {
    if (r_norm < 0.0) return;
    double cx = r_norm / (1.0 + r_norm);
    double cr = 1.0 / (1.0 + r_norm);
    ImPlotPoint center_px = ImPlot::PlotToPixels(cx, 0.0);
    ImPlotPoint edge_px = ImPlot::PlotToPixels(cx + cr, 0.0);
    float radius = (float)std::fabs(edge_px.x - center_px.x);
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    dl->AddCircle(ImVec2((float)center_px.x, (float)center_px.y),
                  radius,
                  IM_COL32(150, 150, 150, 120),
                  96,
                  1.0f);
}

static bool solve_circle_unit_intersection(double cx,
                                           double cy,
                                           double r,
                                           ImPlotPoint* out_point) {
    double d = 0.5 * (cx * cx + cy * cy + 1.0 - r * r);
    double A = cx * cx + cy * cy;
    double B = -2.0 * d * cx;
    double C = d * d - cy * cy;
    double disc = B * B - 4.0 * A * C;
    if (disc <= 0.0) return false;
    double sqrt_disc = std::sqrt(disc);
    double x1 = (-B + sqrt_disc) / (2.0 * A);
    double x2 = (-B - sqrt_disc) / (2.0 * A);
    double chosen = (std::fabs(x1 - 1.0) > 1e-3) ? x1 : x2;
    double y = 0.0;
    if (std::fabs(cy) > 1e-6) {
        y = (d - cx * chosen) / cy;
    } else {
        double val = 1.0 - chosen * chosen;
        if (val < 0.0) val = 0.0;
        y = std::sqrt(val);
        if (cy < 0.0) y = -y;
    }
    if (out_point) {
        *out_point = ImPlotPoint(chosen, y);
    }
    return true;
}

static void draw_reactance_arc(double x_norm) {
    if (std::fabs(x_norm) < 1e-3) return;
    double cy = 1.0 / x_norm;
    double cx = 1.0;
    double radius = 1.0 / std::fabs(x_norm);
    ImPlotPoint other;
    if (!solve_circle_unit_intersection(cx, cy, radius, &other)) {
        return;
    }
    const double pi = 3.14159265358979323846;
    double theta_start = std::atan2(0.0 - cy, 1.0 - cx);
    double theta_end = std::atan2(other.y - cy, other.x - cx);
    if (x_norm > 0.0) {
        if (theta_end < theta_start) theta_end += 2.0 * pi;
    } else {
        if (theta_end > theta_start) theta_end -= 2.0 * pi;
    }
    const int segments = 96;
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImVec2 prev;
    bool has_prev = false;
    for (int i = 0; i <= segments; ++i) {
        double t = (double)i / (double)segments;
        double theta = theta_start + t * (theta_end - theta_start);
        double gx = cx + radius * std::cos(theta);
        double gy = cy + radius * std::sin(theta);
        if (gx * gx + gy * gy > 1.0025) continue;
        ImPlotPoint px = ImPlot::PlotToPixels(gx, gy);
        ImVec2 pos((float)px.x, (float)px.y);
        if (has_prev) {
            dl->AddLine(prev, pos, IM_COL32(150, 150, 150, 120), 1.0f);
        }
        prev = pos;
        has_prev = true;
    }
}

static void draw_smith_grid(const AppState* app) {
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlotPoint center = ImPlot::PlotToPixels(0.0, 0.0);
    ImPlotPoint edge = ImPlot::PlotToPixels(1.0, 0.0);
    float radius = (float)std::fabs(edge.x - center.x);
    dl->AddCircle(ImVec2((float)center.x, (float)center.y),
                  radius,
                  IM_COL32(200, 200, 200, 255),
                  128,
                  2.0f);
    const double r_values[] = {0.2, 0.5, 1.0, 2.0, 5.0};
    for (double r_norm : r_values) {
        draw_resistance_circle(r_norm);
    }
    const double x_values[] = {0.2, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0};
    for (double x_norm : x_values) {
        draw_reactance_arc(x_norm);
    }
    if (app->smith_show_vswr_circles) {
        const double vswr_values[] = {2.0, 3.0, 5.0};
        for (double vs : vswr_values) {
            double gamma = (vs - 1.0) / (vs + 1.0);
            float circle_radius = radius * (float)gamma;
            dl->AddCircle(ImVec2((float)center.x, (float)center.y),
                          circle_radius,
                          IM_COL32(255, 200, 0, 90),
                          96,
                          1.0f);
        }
    }
    dl->AddLine(ImPlot::PlotToPixels(-1.0, 0.0),
                ImPlot::PlotToPixels(1.0, 0.0),
                IM_COL32(255, 255, 255, 80),
                1.0f);
    dl->AddLine(ImPlot::PlotToPixels(0.0, -1.0),
                ImPlot::PlotToPixels(0.0, 1.0),
                IM_COL32(255, 255, 255, 80),
                1.0f);
}

static void plot_smith_reflection_data(AppState* app) {
    const int max_samples = 1024;
    static double gamma_r[max_samples];
    static double gamma_i[max_samples];
    int count = (int)app->sparam_s11_mag.size();
    if (count <= 0) return;
    if (count > max_samples) count = max_samples;
    for (int i = 0; i < count; ++i) {
        double mag = app->sparam_s11_mag[i];
        gamma_r[i] = mag;
        gamma_i[i] = 0.0;
    }
    ImPlot::PlotLine("S11", gamma_r, gamma_i, count);
    if (app->smith_freq_index >= 0 && app->smith_freq_index < count) {
        int idx = app->smith_freq_index;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0, 0, 1));
        ImPlot::PlotScatter("Current", &gamma_r[idx], &gamma_i[idx], 1);
        ImPlot::PopStyleColor();
    }
}

static void draw_smith_chart(AppState* app) {
    if (!app || !app->smith_chart_open) return;
    ImGui::SetNextWindowSize(ImVec2(620.0f, 680.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Smith Chart", &app->smith_chart_open)) {
        ImGui::End();
        return;
    }

    ImGui::Checkbox("Impedance", &app->smith_show_impedance);
    ImGui::SameLine();
    if (!app->smith_show_impedance) {
        ImGui::TextDisabled("(Admittance view)");
    }
    ImGui::Checkbox("VSWR Circles", &app->smith_show_vswr_circles);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Z0 (Ω)", &app->smith_z0, 1.0, 10.0, "%.1f");

    ImGui::Separator();

    if (ImPlot::BeginPlot("##SmithChart",
                          ImVec2(-1, -1),
                          ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxis(ImAxis_X1, nullptr, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.1, 1.1, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.1, 1.1, ImGuiCond_Always);
        draw_smith_grid(app);
        if (!app->sparam_s11_mag.empty()) {
            plot_smith_reflection_data(app);
        } else {
            ImPlot::EndPlot();
            ImGui::TextUnformatted("Load S-parameters to populate chart.");
            ImGui::End();
            return;
        }
        ImPlot::EndPlot();
    }

    int count = (int)app->sparam_freq_hz.size();
    if (count > 0) {
        if (app->smith_freq_index >= count) app->smith_freq_index = count - 1;
        if (app->smith_freq_index < 0) app->smith_freq_index = 0;
        ImGui::Separator();
        ImGui::TextUnformatted("Probe Frequency");
        if (ImGui::SliderInt("##smith_freq_idx", &app->smith_freq_index, 0, count - 1)) {
            /* slider moved */
        }
        int idx = app->smith_freq_index;
        double freq_ghz = app->sparam_freq_hz[idx] * 1e-9;
        double s11_db = app->sparam_s11_db[idx];
        double vswr = app->sparam_vswr[idx];
        ImGui::Text("Frequency: %.3f GHz", freq_ghz);
        ImGui::Text("S11: %.2f dB", s11_db);
        ImGui::Text("VSWR: %.2f", vswr);
    } else {
        ImGui::TextUnformatted("No S-parameter data available.");
    }

    ImGui::End();
}

static void draw_sparameter_window(AppState* app) {
    if (!app || !app->sparam_window_open) return;

    ImGui::SetNextWindowSize(ImVec2(850.0f, 620.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("S-Parameters", &app->sparam_window_open)) {
        ImGui::End();
        return;
    }

    ImGui::InputText("CSV Path", app->sparam_csv_path, sizeof(app->sparam_csv_path));
    ImGui::SameLine();
    if (ImGui::Button("Load CSV")) {
        if (load_sparameter_csv(app, app->sparam_csv_path)) {
            ui_log_add(app, "Loaded S-parameters: %s", app->sparam_csv_path);
        } else {
            ui_log_add(app, "Failed to load: %s", app->sparam_csv_path);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Data")) {
        clear_sparameter_data(app);
        ui_log_add(app, "Cleared S-parameter data");
    }

    if (!app->sparam_data_loaded || app->sparam_freq_hz.empty()) {
        ImGui::Separator();
        ImGui::TextUnformatted("No S-parameter data loaded.");
        ImGui::TextWrapped("Run CLI sweep to generate sweep_results.csv, then load it here.");
        ImGui::End();
        return;
    }

    const int count = (int)app->sparam_freq_hz.size();
    const double* freq_ghz = app->sparam_freq_ghz.data();
    const double* s21_db = app->sparam_s21_db.data();
    const double* s11_db = app->sparam_s11_db.data();
    const double* vswr = app->sparam_vswr.data();

    ImGui::Separator();
    if (ImPlot::BeginPlot("S21 - Transmission", ImVec2(-1, 260))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "S21 (dB)");
        ImPlot::PlotLine("S21", freq_ghz, s21_db, count);

        double threshold = app->sparam_peak_db - 3.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0, 0, 1));
        ImPlot::PlotInfLines("-3 dB", &threshold, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();

        ImPlotRect limits = ImPlot::GetPlotLimits();
        if (app->sparam_f_low_hz > 0.0) {
            double low = app->sparam_f_low_hz * 1e-9;
            double xs[2] = {low, low};
            double ys[2] = {limits.Y.Min, limits.Y.Max};
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0, 1, 0, 1));
            ImPlot::PlotLine("f_low", xs, ys, 2);
            ImPlot::PopStyleColor();
        }
        if (app->sparam_f_high_hz > app->sparam_f_low_hz) {
            double high = app->sparam_f_high_hz * 1e-9;
            double xs[2] = {high, high};
            double ys[2] = {limits.Y.Min, limits.Y.Max};
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0, 1, 0, 1));
            ImPlot::PlotLine("f_high", xs, ys, 2);
            ImPlot::PopStyleColor();
        }

        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("S11 - Reflection", ImVec2(-1, 220))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "S11 (dB)");
        ImPlot::PlotLine("S11", freq_ghz, s11_db, count);
        double match_line = -10.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 1, 0, 1));
        ImPlot::PlotInfLines("-10 dB", &match_line, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("VSWR", ImVec2(-1, 160))) {
        ImPlot::SetupAxis(ImAxis_X1, "Frequency (GHz)");
        ImPlot::SetupAxis(ImAxis_Y1, "VSWR");
        ImPlot::SetupAxisLimits(ImAxis_Y1, 1.0, 5.0, ImGuiCond_Always);
        ImPlot::PlotLine("VSWR", freq_ghz, vswr, count);
        double vswr_limit = 2.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1, 0.5f, 0, 1));
        ImPlot::PlotInfLines("VSWR=2", &vswr_limit, 1, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PopStyleColor();
        ImPlot::EndPlot();
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Bandwidth Analysis");
    ImGui::Columns(2, "sparam_stats", false);
    ImGui::SetColumnWidth(0, 160.0f);

    ImGui::Text("Peak S21:"); ImGui::NextColumn();
    ImGui::Text("%.2f dB @ %.3f GHz", app->sparam_peak_db, app->sparam_peak_freq_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("Center Freq:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_center_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("Bandwidth:"); ImGui::NextColumn();
    if (app->sparam_bandwidth_hz > 0.0 && app->sparam_f_center_hz > 0.0) {
        ImGui::Text("%.3f GHz (%.1f%%)",
                    app->sparam_bandwidth_hz * 1e-9,
                    (app->sparam_bandwidth_hz / app->sparam_f_center_hz) * 100.0);
    } else {
        ImGui::TextDisabled("N/A");
    }
    ImGui::NextColumn();

    ImGui::Text("f_low:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_low_hz * 1e-9);
    ImGui::NextColumn();

    ImGui::Text("f_high:"); ImGui::NextColumn();
    ImGui::Text("%.3f GHz", app->sparam_f_high_hz * 1e-9);
    ImGui::Columns(1);

    ImGui::End();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim);

/* Log panel --------------------------------------------------------------- */
static void draw_log_panel(AppState* app);

static void render_grid_overlay(RenderContext* render, const SimulationState* sim, SDL_Color color, int grid_step) {
    if (!render || !sim) return;
    SDL_Renderer* rr = render->renderer;
    if (!rr) return;
    if (grid_step < 1) grid_step = 1;
    int width = sim->nx * render->scale;
    int height = sim->ny * render->scale;
    int offset_x = (int)std::lround(render->offset_x);
    int offset_y = (int)std::lround(render->offset_y);
    SDL_SetRenderDrawColor(rr, color.r, color.g, color.b, color.a);
    for (int i = 0; i <= sim->nx; i += grid_step) {
        int x = i * render->scale + offset_x;
        SDL_RenderDrawLine(rr, x, offset_y, x, height + offset_y);
    }
    for (int j = 0; j <= sim->ny; j += grid_step) {
        int y = j * render->scale + offset_y;
        SDL_RenderDrawLine(rr, offset_x, y, width + offset_x, y);
    }
}

static SDL_Color material_color_from_epsilon(double eps) {
    SDL_Color color{128, 128, 128, 255};
    if (eps < 1.05) {
        color = SDL_Color{200, 220, 255, 255};
        return color;
    }
    if (eps > 100.0) {
        color = SDL_Color{184, 115, 51, 255};
        return color;
    }
    const Material* best = NULL;
    double best_diff = 1e9;
    int count = material_library_get_count();
    for (int i = 0; i < count; ++i) {
        const Material* mat = material_library_get_by_index(i);
        if (!mat) continue;
        if (mat->type == MAT_TYPE_PEC || mat->type == MAT_TYPE_PMC) continue;
        double diff = std::fabs(mat->epsilon_r - eps);
        if (diff < best_diff) {
            best_diff = diff;
            best = mat;
        }
    }
    if (best && best_diff < 0.75) {
        color = SDL_Color{best->color_r, best->color_g, best->color_b, 255};
    }
    return color;
}

static void render_material_distribution(RenderContext* render,
                                         const SimulationState* sim,
                                         int scale) {
    if (!render || !render->renderer || !sim || scale <= 0) return;
    int offset_x = (int)std::lround(render->offset_x);
    int offset_y = (int)std::lround(render->offset_y);
    for (int j = 0; j < sim->ny; ++j) {
        for (int i = 0; i < sim->nx; ++i) {
            double eps = fdtd_epsilon_at(sim, i, j);
            SDL_Color color = material_color_from_epsilon(eps);
            SDL_Rect rect = {i * scale + offset_x, j * scale + offset_y, scale, scale};
            SDL_SetRenderDrawColor(render->renderer, color.r, color.g, color.b, color.a);
            SDL_RenderFillRect(render->renderer, &rect);
        }
    }
}

static void render_material_overlay(RenderContext* render,
                                    const SimulationState* sim,
                                    int scale,
                                    float alpha) {
    if (!render || !render->renderer || !sim || scale <= 0) return;
    if (alpha <= 0.0f) return;
    int offset_x = (int)std::lround(render->offset_x);
    int offset_y = (int)std::lround(render->offset_y);
    Uint8 a = (Uint8)(std::fmin(std::fmax(alpha, 0.0f), 1.0f) * 255.0f);
    for (int j = 0; j < sim->ny; ++j) {
        for (int i = 0; i < sim->nx; ++i) {
            double eps = fdtd_epsilon_at(sim, i, j);
            SDL_Color base = material_color_from_epsilon(eps);
            SDL_Rect rect = {i * scale + offset_x, j * scale + offset_y, scale, scale};
            SDL_SetRenderDrawColor(render->renderer, base.r, base.g, base.b, a);
            SDL_RenderFillRect(render->renderer, &rect);
        }
    }
}

static void render_material_outlines(RenderContext* render,
                                     const SimulationState* sim,
                                     int scale) {
    if (!render || !render->renderer || !sim || scale <= 0) return;
    int offset_x = (int)std::lround(render->offset_x);
    int offset_y = (int)std::lround(render->offset_y);
    SDL_SetRenderDrawColor(render->renderer, 255, 255, 255, 160);
    for (int j = 0; j < sim->ny - 1; ++j) {
        for (int i = 0; i < sim->nx - 1; ++i) {
            double eps_here = fdtd_epsilon_at(sim, i, j);
            double eps_right = fdtd_epsilon_at(sim, i + 1, j);
            double eps_down = fdtd_epsilon_at(sim, i, j + 1);
            if (std::fabs(eps_here - eps_right) > 0.1) {
                int x = (i + 1) * scale + offset_x;
                int y0 = j * scale + offset_y;
                int y1 = (j + 1) * scale + offset_y;
                SDL_RenderDrawLine(render->renderer, x, y0, x, y1);
            }
            if (std::fabs(eps_here - eps_down) > 0.1) {
                int x0 = i * scale + offset_x;
                int x1 = (i + 1) * scale + offset_x;
                int y = (j + 1) * scale + offset_y;
                SDL_RenderDrawLine(render->renderer, x0, y, x1, y);
            }
        }
    }
}

static void apply_paint_brush(SimulationState* sim,
                              int gx,
                              int gy,
                              int brush_radius,
                              int paint_type,
                              double eps) {
    if (!sim) return;
    if (brush_radius < 0) brush_radius = 0;

    int nx = sim->nx;
    int ny = sim->ny;

    int r = brush_radius;
    int r2 = r * r;
    for (int di = -r; di <= r; ++di) {
        int i = gx + di;
        if (i < 0 || i >= nx) continue;
        for (int dj = -r; dj <= r; ++dj) {
            int j = gy + dj;
            if (j < 0 || j >= ny) continue;
            if (di * di + dj * dj <= r2) {
                paint_material_at(sim, i, j, paint_type, eps);
            }
        }
    }
}

static void blocks_paint_pec(SimulationState* sim, int gx, int gy, int brush_radius) {
    apply_paint_brush(sim, gx, gy, brush_radius, 1, 1.0);
}

static void blocks_paint_pmc(SimulationState* sim, int gx, int gy, int brush_radius) {
    apply_paint_brush(sim, gx, gy, brush_radius, 2, 1.0);
}

static void blocks_paint_dielectric(SimulationState* sim,
                                    int gx,
                                    int gy,
                                    int brush_radius,
                                    double epsr) {
    apply_paint_brush(sim, gx, gy, brush_radius, 3, epsr);
}

static void apply_selected_paint(SimulationState* sim,
                                 int gx,
                                 int gy,
                                 int brush_radius,
                                 int paint_material_type,
                                 double paint_epsilon) {
    if (!sim) return;
    switch (paint_material_type) {
        case 0:
            blocks_paint_pec(sim, gx, gy, brush_radius);
            break;
        case 1:
            blocks_paint_pmc(sim, gx, gy, brush_radius);
            break;
        case 2:
        default:
            blocks_paint_dielectric(sim, gx, gy, brush_radius, paint_epsilon);
            break;
    }
}

static void apply_paint_with_material(SimulationState* sim,
                                      int gx,
                                      int gy,
                                      int brush_radius,
                                      int material_id) {
    if (!sim) return;
    const Material* mat = material_library_get_by_id(material_id);
    if (!mat) return;

    int paint_type = 0;
    double epsilon = 1.0;

    switch (mat->type) {
        case MAT_TYPE_PEC:
            paint_type = 1;
            epsilon = 1.0;
            break;
        case MAT_TYPE_PMC:
            paint_type = 2;
            epsilon = 1.0;
            break;
        case MAT_TYPE_DIELECTRIC:
        case MAT_TYPE_LOSSY_DIEL:
        default:
            paint_type = 3;
            epsilon = mat->epsilon_r;
            break;
    }

    apply_paint_brush(sim, gx, gy, brush_radius, paint_type, epsilon);
}

static int find_source_at_position(const SimulationState* sim, int gx, int gy, int tolerance) {
    if (!sim) return -1;
    if (tolerance < 0) tolerance = 0;
    int tol2 = tolerance * tolerance;
    for (int i = 0; i < MAX_SRC; ++i) {
        if (!sim->sources[i].active) continue;
        int sx = sim->sources[i].ix;
        int sy = sim->sources[i].iy;
        int dx = gx - sx;
        int dy = gy - sy;
        if (dx * dx + dy * dy <= tol2) {
            return i;
        }
    }
    return -1;
}

// Theme and visual styling helpers -----------------------------------------

enum ThemePreset {
    THEME_PRESET_DARK_PRO = 0,
    THEME_PRESET_BLENDER = 1,
    THEME_PRESET_LIGHT = 2,
    THEME_PRESET_HIGH_CONTRAST = 3
};

static const char* theme_preset_name(ThemePreset preset) {
    switch (preset) {
        case THEME_PRESET_DARK_PRO: return "Dark Professional";
        case THEME_PRESET_BLENDER: return "Blender Inspired";
        case THEME_PRESET_LIGHT: return "Light Mode";
        case THEME_PRESET_HIGH_CONTRAST: return "High Contrast";
        default: return "Unknown";
    }
}

static const ImVec4 ACCENT_PALETTES[6] = {
    ImVec4(0.26f, 0.59f, 0.98f, 1.00f),
    ImVec4(0.20f, 0.70f, 0.50f, 1.00f),
    ImVec4(0.90f, 0.40f, 0.20f, 1.00f),
    ImVec4(0.70f, 0.30f, 0.80f, 1.00f),
    ImVec4(0.90f, 0.20f, 0.30f, 1.00f),
    ImVec4(0.90f, 0.70f, 0.10f, 1.00f)
};

static const char* ACCENT_NAMES[6] = {
    "Blue", "Teal", "Orange", "Purple", "Red", "Gold"
};

static ImVec4 accent_from_palette(int accent_id) {
    if (accent_id < 0 || accent_id >= 6) accent_id = 0;
    return ACCENT_PALETTES[accent_id];
}

static void apply_accent_to_style(ImGuiStyle& style, const ImVec4& accent) {
    ImVec4 hover = ImVec4(std::min(accent.x + 0.08f, 1.0f),
                          std::min(accent.y + 0.08f, 1.0f),
                          std::min(accent.z + 0.08f, 1.0f),
                          accent.w);
    ImVec4 active = ImVec4(std::min(accent.x + 0.12f, 1.0f),
                           std::min(accent.y + 0.12f, 1.0f),
                           std::min(accent.z + 0.12f, 1.0f),
                           accent.w);

    style.Colors[ImGuiCol_Header] = accent;
    style.Colors[ImGuiCol_HeaderHovered] = hover;
    style.Colors[ImGuiCol_HeaderActive] = active;
    style.Colors[ImGuiCol_Button] =
        ImVec4(accent.x * 0.5f, accent.y * 0.5f, accent.z * 0.5f, 0.70f);
    style.Colors[ImGuiCol_ButtonHovered] = hover;
    style.Colors[ImGuiCol_ButtonActive] = active;
    style.Colors[ImGuiCol_CheckMark] = active;
    style.Colors[ImGuiCol_SliderGrab] = accent;
    style.Colors[ImGuiCol_SliderGrabActive] = active;
    style.Colors[ImGuiCol_Tab] =
        ImVec4(accent.x * 0.55f, accent.y * 0.55f, accent.z * 0.55f, 0.80f);
    style.Colors[ImGuiCol_TabHovered] = hover;
    style.Colors[ImGuiCol_TabActive] = active;
}

static void apply_theme(ThemePreset preset, const ImVec4& accent) {
    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark();

    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 3.0f;
    style.TabRounding = 4.0f;
    style.ScrollbarSize = 14.0f;
    style.FramePadding = ImVec2(8.0f, 4.0f);
    style.ItemSpacing = ImVec2(8.0f, 6.0f);
    style.ItemInnerSpacing = ImVec2(6.0f, 4.0f);
    style.IndentSpacing = 20.0f;
    style.WindowPadding = ImVec2(8.0f, 8.0f);

    ImVec4* colors = style.Colors;
    switch (preset) {
        case THEME_PRESET_DARK_PRO: {
            colors[ImGuiCol_WindowBg] = ImVec4(0.10f, 0.10f, 0.11f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.12f, 0.12f, 0.13f, 1.00f);
            colors[ImGuiCol_FrameBg] = ImVec4(0.18f, 0.18f, 0.20f, 1.00f);
            colors[ImGuiCol_FrameBgHovered] = ImVec4(0.25f, 0.25f, 0.28f, 1.00f);
            colors[ImGuiCol_FrameBgActive] = ImVec4(0.30f, 0.30f, 0.35f, 1.00f);
            colors[ImGuiCol_TitleBg] = ImVec4(0.08f, 0.08f, 0.09f, 1.00f);
            colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.12f, 0.15f, 1.00f);
            break;
        }
        case THEME_PRESET_BLENDER: {
            colors[ImGuiCol_WindowBg] = ImVec4(0.16f, 0.16f, 0.16f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.18f, 0.18f, 0.18f, 1.00f);
            colors[ImGuiCol_FrameBg] = ImVec4(0.22f, 0.22f, 0.24f, 1.00f);
            colors[ImGuiCol_TitleBg] = ImVec4(0.12f, 0.12f, 0.12f, 1.00f);
            colors[ImGuiCol_TitleBgActive] = ImVec4(0.26f, 0.46f, 0.76f, 1.00f);
            break;
        }
        case THEME_PRESET_LIGHT: {
            ImGui::StyleColorsLight();
            style = ImGui::GetStyle();
            style.WindowRounding = 6.0f;
            style.FrameRounding = 4.0f;
            style.TabRounding = 4.0f;
            style.FramePadding = ImVec2(8.0f, 5.0f);
            colors = style.Colors;
            colors[ImGuiCol_WindowBg] = ImVec4(0.95f, 0.95f, 0.96f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.97f, 0.97f, 0.98f, 1.00f);
            colors[ImGuiCol_FrameBg] = ImVec4(0.90f, 0.90f, 0.92f, 1.00f);
            colors[ImGuiCol_FrameBgHovered] = ImVec4(0.86f, 0.90f, 0.96f, 1.00f);
            colors[ImGuiCol_TitleBg] = ImVec4(0.88f, 0.88f, 0.90f, 1.00f);
            colors[ImGuiCol_TitleBgActive] = ImVec4(0.82f, 0.84f, 0.88f, 1.00f);
            break;
        }
        case THEME_PRESET_HIGH_CONTRAST: {
            ImGui::StyleColorsDark();
            style = ImGui::GetStyle();
            style.WindowRounding = 0.0f;
            style.FrameRounding = 0.0f;
            style.TabRounding = 0.0f;
            colors = style.Colors;
            colors[ImGuiCol_WindowBg] = ImVec4(0.02f, 0.02f, 0.02f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.03f, 0.03f, 0.03f, 1.00f);
            colors[ImGuiCol_Text] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
            colors[ImGuiCol_FrameBg] = ImVec4(0.10f, 0.10f, 0.10f, 1.00f);
            colors[ImGuiCol_FrameBgHovered] = ImVec4(0.15f, 0.15f, 0.15f, 1.00f);
            colors[ImGuiCol_FrameBgActive] = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
            break;
        }
        default:
            break;
    }

    apply_accent_to_style(style, accent);
}

static SDL_Color theme_viewport_clear_color(ThemePreset preset) {
    SDL_Color color;
    switch (preset) {
        case THEME_PRESET_LIGHT:
            color = SDL_Color{235, 236, 240, 255};
            break;
        case THEME_PRESET_HIGH_CONTRAST:
            color = SDL_Color{4, 4, 4, 255};
            break;
        case THEME_PRESET_BLENDER:
        case THEME_PRESET_DARK_PRO:
        default:
            color = SDL_Color{8, 8, 10, 255};
            break;
    }
    return color;
}

static const char* kColormapNames[] = {"Classic", "Viridis", "Plasma"};
static const int kColormapCount = 3;

static bool rebootstrap_simulation(const SimulationConfig* cfg,
                                   SimulationBootstrap* bootstrap,
                                   SimulationState** sim_out,
                                   Scope* scope,
                                   int scale) {
    if (!cfg || !bootstrap || !sim_out || !scope) return false;

    simulation_bootstrap_shutdown(bootstrap);

    if (!simulation_bootstrap_with_config(cfg, bootstrap)) {
        std::fprintf(stderr, "Failed to reinitialise simulation from wizard configuration\n");
        return false;
    }

    *sim_out = bootstrap->sim;

    scope_free(scope);
    if (!scope_init(scope, (*sim_out)->nx * scale)) {
        std::fprintf(stderr, "Failed to initialise oscilloscope after rebootstrap\n");
        simulation_bootstrap_shutdown(bootstrap);
        *sim_out = NULL;
        return false;
    }

    return true;
}

static void apply_wizard_materials_to_sim(const WizardState& wizard,
                                          SimulationBootstrap* bootstrap,
                                          SimulationState* sim) {
    if (!bootstrap || !sim) return;

    int count = wizard.cfg.material_rect_count;
    if (count < 0) count = 0;
    if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;

    bootstrap->config.material_rect_count = count;
    sim->config.material_rect_count = count;
    for (int i = 0; i < count; ++i) {
        bootstrap->config.material_rects[i] = wizard.cfg.material_rects[i];
        sim->config.material_rects[i] = wizard.cfg.material_rects[i];
    }

    materials_init(sim);
}

static void sync_sim_to_wizard_config(SimulationState* sim, WizardState* wizard) {
    if (!sim || !wizard) return;

    int nx = sim->nx;
    int ny = sim->ny;
    double nx1 = (nx > 1) ? (double)(nx - 1) : 1.0;
    double ny1 = (ny > 1) ? (double)(ny - 1) : 1.0;

    int max_src = wizard->cfg.source_count;
    if (max_src < 0) max_src = 0;
    if (max_src > MAX_SRC) max_src = MAX_SRC;

    for (int i = 0; i < max_src; ++i) {
        const Source& s = sim->sources[i];
        SourceConfigSpec* spec = &wizard->cfg.source_configs[i];
        spec->active = s.active ? 1 : 0;
        spec->x = (double)s.ix / nx1;
        spec->y = (double)s.iy / ny1;
        spec->type = s.type;
        spec->field = s.field;
        spec->amp = s.amp;
        spec->freq = s.freq;
        spec->sigma2 = s.sigma2;
        /* Expr text is not round-tripped here; leave cfg.expr as-is. */
    }
}

static void draw_blocks_panel(WizardState& wizard,
                              SimulationBootstrap* bootstrap,
                              SimulationState* sim,
                              AppState* app) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Materials / Blocks")) return;

    ImGui::Indent();
    int count = wizard.cfg.material_rect_count;
    if (count < 0) count = 0;
    if (count > CONFIG_MAX_MATERIAL_RECTS) count = CONFIG_MAX_MATERIAL_RECTS;
    wizard.cfg.material_rect_count = count;

    ImGui::Text("Blocks: %d / %d", count, CONFIG_MAX_MATERIAL_RECTS);
    ImGui::SameLine();
    bool can_add = (count < CONFIG_MAX_MATERIAL_RECTS);
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ Add Block")) {
        MaterialRectSpec& r = wizard.cfg.material_rects[count];
        r.x0 = 0.25;
        r.y0 = 0.25;
        r.x1 = 0.75;
        r.y1 = 0.75;
        r.epsr = 4.0;
        r.sigma = 0.0;
        r.tag = 0;
        wizard.cfg.material_rect_count = count + 1;
        if (sim && bootstrap) {
            apply_wizard_materials_to_sim(wizard, bootstrap, sim);
        }
        if (app) {
            ui_log_add(app, "Added block #%d (default material)", count);
        }
    }
    ImGui::SameLine();
    if (!can_add) ImGui::BeginDisabled();
    if (ImGui::Button("+ Add Block from Material")) {
        ImGui::OpenPopup("AddBlockMaterialPopup");
    }
    if (!can_add) ImGui::EndDisabled();
    if (ImGui::BeginPopup("AddBlockMaterialPopup")) {
        ImGui::TextUnformatted("Select material:");
        ImGui::Separator();
        int mat_count = material_library_get_count();
        for (int i = 0; i < mat_count; ++i) {
            const Material* mat = material_library_get_by_index(i);
            if (!mat) continue;
            if (ImGui::Selectable(mat->name)) {
                if (wizard.cfg.material_rect_count < CONFIG_MAX_MATERIAL_RECTS) {
                    int idx = wizard.cfg.material_rect_count;
                    MaterialRectSpec& r = wizard.cfg.material_rects[idx];
                    r.x0 = 0.3;
                    r.y0 = 0.3;
                    r.x1 = 0.7;
                    r.y1 = 0.7;
                    apply_material_to_rect(&r, mat);
                    wizard.cfg.material_rect_count = idx + 1;
                    if (sim && bootstrap) {
                        apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                    }
                    if (app) {
                        ui_log_add(app, "Added block #%d (%s)", idx, mat->name);
                    }
                }
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    ImGui::Separator();

    const Material* selected_mat = nullptr;
    if (app && app->selected_material_id >= 0) {
        selected_mat = material_library_get_by_id(app->selected_material_id);
    }

    if (selected_mat) {
        ImGui::Checkbox("Filter by selected material", &app->filter_blocks_by_material);
        ImGui::SameLine();
        ImGui::TextDisabled("Filter: %s", selected_mat->name);
        ImGui::Checkbox("Auto-filter on selection", &app->auto_filter_blocks_on_select);
    } else {
        ImGui::BeginDisabled();
        bool dummy = false;
        ImGui::Checkbox("Filter by selected material", &dummy);
        ImGui::SameLine();
        ImGui::TextDisabled("No material selected");
        ImGui::EndDisabled();
        app->filter_blocks_by_material = false;
    }

    ImGui::Separator();

    if (wizard.cfg.material_rect_count == 0) {
        ImGui::TextDisabled("No material blocks configured.");
    } else {
        for (int i = 0; i < wizard.cfg.material_rect_count; ++i) {
            MaterialRectSpec& r = wizard.cfg.material_rects[i];
            if (app && app->filter_blocks_by_material && selected_mat) {
                if (!rect_matches_material(r, selected_mat)) {
                    ImGui::PopID();
                    continue;
                }
            }
            ImGui::PushID(i);
            bool tree_open = ImGui::CollapsingHeader("", ImGuiTreeNodeFlags_DefaultOpen);
            ImGui::SameLine();
            ImGui::Text("Block %d", i);
            if (tree_open) {
                const Material* current_mat = guess_material_for_rect(r);
                ImVec4 color = current_mat
                                   ? ImVec4(current_mat->color_r / 255.0f,
                                            current_mat->color_g / 255.0f,
                                            current_mat->color_b / 255.0f,
                                            1.0f)
                                   : ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                if (selected_mat && rect_matches_material(r, selected_mat)) {
                    ImGui::PushStyleColor(ImGuiCol_Text, color);
                }
                ImGui::ColorButton("##blkcol",
                                   color,
                                   ImGuiColorEditFlags_NoTooltip |
                                       ImGuiColorEditFlags_NoDragDrop,
                                   ImVec2(16, 16));
                ImGui::SameLine();
                if (selected_mat && rect_matches_material(r, selected_mat)) {
                    ImGui::PopStyleColor();
                }
                const char* combo_label =
                    current_mat ? current_mat->name
                                : (r.tag == 1 ? "PEC"
                                              : (r.tag == 2 ? "PMC" : "Custom"));
                if (ImGui::BeginCombo("Material", combo_label)) {
                    int mat_count = material_library_get_count();
                    for (int m = 0; m < mat_count; ++m) {
                        const Material* mat = material_library_get_by_index(m);
                        if (!mat) continue;
                        bool selected = (current_mat == mat);
                        if (ImGui::Selectable(mat->name, selected)) {
                            apply_material_to_rect(&r, mat);
                            if (sim && bootstrap) {
                                apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                            }
                        }
                        if (selected) ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }

                const char* type_items[] = {"Dielectric", "PEC", "PMC"};
                int type_idx = (r.tag == 1) ? 1 : (r.tag == 2) ? 2 : 0;
                if (ImGui::Combo("Type", &type_idx, type_items, IM_ARRAYSIZE(type_items))) {
                    r.tag = (type_idx == 1) ? 1 : (type_idx == 2) ? 2 : 0;
                    if (sim && bootstrap) apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                }

                bool changed = false;
                float x0 = (float)r.x0;
                float y0 = (float)r.y0;
                float x1 = (float)r.x1;
                float y1 = (float)r.y1;
                changed |= ImGui::SliderFloat("x0", &x0, 0.0f, 0.99f, "%.3f");
                changed |= ImGui::SliderFloat("y0", &y0, 0.0f, 0.99f, "%.3f");
                changed |= ImGui::SliderFloat("x1", &x1, 0.0f, 1.0f, "%.3f");
                changed |= ImGui::SliderFloat("y1", &y1, 0.0f, 1.0f, "%.3f");
                if (changed) {
                    r.x0 = clamp01(x0);
                    r.y0 = clamp01(y0);
                    r.x1 = clamp01(x1);
                    r.y1 = clamp01(y1);
                    if (r.x1 < r.x0) r.x1 = r.x0 + 0.01;
                    if (r.y1 < r.y0) r.y1 = r.y0 + 0.01;
                }

                if (sim) {
                    ImGui::TextDisabled("Meters: (%.3f, %.3f) to (%.3f, %.3f)",
                                        r.x0 * sim->lx,
                                        r.y0 * sim->ly,
                                        r.x1 * sim->lx,
                                        r.y1 * sim->ly);
                }

                if (r.tag == 0) {
                    if (ImGui::InputDouble("epsr", &r.epsr, 0.1, 1.0)) changed = true;
                    if (ImGui::InputDouble("sigma", &r.sigma, 0.0, 0.1)) changed = true;
                } else {
                    ImGui::TextDisabled("PEC/PMC override eps/sigma");
                }

                if (changed && sim && bootstrap) {
                    apply_wizard_materials_to_sim(wizard, bootstrap, sim);
                }

                if (ImGui::Button("Delete Block")) {
                    remove_block(&wizard, bootstrap, sim, app, i);
                    ImGui::PopID();
                    break;
                }
            }
            ImGui::PopID();
        }
    }
    ImGui::Unindent();
}

/* Probes panel ------------------------------------------------------------ */
static void draw_probes_panel(const SimulationState* sim) {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Probes")) return;

    ImGui::Indent();
#if EMWAVE_ENABLE_PORTS
    if (!sim) {
        ImGui::Unindent();
        return;
    }
    ImGui::Text("Ports enabled: %s", sim->ports_on ? "yes" : "no");
    ImGui::Separator();
    for (int p = 0; p < MAX_PORTS; ++p) {
        const Port& port = sim->ports[p];
        if (!port.active) continue;
        ImGui::Text("Port %d: x=%d, y0=%d, y1=%d, n=%d", p, port.x, port.y0, port.y1, port.n);
    }
#else
    ImGui::TextUnformatted("Ports are disabled in this build.");
#endif
    ImGui::Unindent();
}

/* Log panel --------------------------------------------------------------- */
static void draw_log_panel(AppState* app) {
    if (!ImGui::Begin("Log")) {
        ImGui::End();
        return;
    }
    if (!app) {
        ImGui::End();
        return;
    }
    for (int i = 0; i < app->log_count; ++i) {
        ImGui::TextUnformatted(app->log_lines[i]);
    }
    ImGui::End();
}

int main(int argc, char** argv) {
    HeadlessComposerOpts headless{};
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (std::strcmp(a, "--composer-headless") == 0) {
            headless.enabled = true;
        } else if (std::strncmp(a, "--ch-format=", 12) == 0) {
            const char* v = a + 12;
            if (_stricmp(v, "bmp") == 0) headless.fmt = 0;
            else if (_stricmp(v, "png") == 0) headless.fmt = 1;
            else if (_stricmp(v, "mp4") == 0) headless.fmt = 2;
            else if (_stricmp(v, "gif") == 0) headless.fmt = 3;
        } else if (std::strncmp(a, "--ch-fps=", 9) == 0) {
            headless.fps = std::atoi(a + 9);
        } else if (std::strncmp(a, "--ch-frames=", 12) == 0) {
            headless.frames = std::atoi(a + 12);
        } else if (std::strncmp(a, "--ch-output=", 12) == 0) {
            std::snprintf(headless.output_name, sizeof(headless.output_name), "%s", a + 12);
        } else if (std::strncmp(a, "--ch-width=", 11) == 0) {
            headless.res_w = std::atoi(a + 11);
        } else if (std::strncmp(a, "--ch-height=", 12) == 0) {
            headless.res_h = std::atoi(a + 12);
        } else if (std::strncmp(a, "--ch-template=", 14) == 0) {
            const char* v = a + 14;
            if (_stricmp(v, "blank") == 0) headless.tmpl_idx = 0;
            else if (_stricmp(v, "field+legend") == 0) headless.tmpl_idx = 1;
            else if (_stricmp(v, "field+scope") == 0) headless.tmpl_idx = 2;
            else if (_stricmp(v, "field+fft+legend") == 0) headless.tmpl_idx = 3;
        } else if (std::strncmp(a, "--ch-layout=", 12) == 0) {
            std::snprintf(headless.layout_path, sizeof(headless.layout_path), "%s", a + 12);
            headless.enabled = true;
            headless.tmpl_idx = -2; // signal external layout load
        }
    }
    if (headless.fmt < 0 || headless.fmt > 3) headless.fmt = 1;
    if (headless.fps <= 0) headless.fps = 30;
    if (headless.frames <= 0) headless.frames = 60;
    if (headless.res_w <= 0) headless.res_w = 1280;
    if (headless.res_h <= 0) headless.res_h = 720;
    if (headless.output_name[0] == '\0') {
        std::snprintf(headless.output_name, sizeof(headless.output_name), "composer_page_1");
    }
    SimulationBootstrap bootstrap;
    int bootstrap_status = simulation_bootstrap_from_args(argc, argv, &bootstrap);
    if (bootstrap_status == 0) {
        return 0;  /* help shown or invalid config */
    }
    if (bootstrap_status < 0) {
        return 1;
    }

    SimulationState* sim = bootstrap.sim;
    if (!sim) {
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }

    material_library_init();

    WizardState wizard;
    wizard_init_from_config(wizard, &bootstrap.config);

    AppState app = {};
    app.basic_mode = true;
    app.selected_source = 0;
    app.selected_block = 0;
    app.placing_source = false;
    app.placing_block = false;
    app.block_first_set = false;
    app.block_first_i = app.block_first_j = -1;
    app.last_click_i = app.last_click_j = -1;
    app.show_scope_window = true;
    app.show_scene_panel = true;
    app.show_sources_panel = true;
    app.show_blocks_panel = true;
    app.show_run_panel = true;
    app.show_grid_panel = true;
    app.show_run_settings_panel = true;
    app.show_probes_panel = true;
    app.show_log_panel = true;
    app.show_expression_panel = true;
    app.show_scenes_panel = true;
    app.show_grid_overlay = true;
    app.theme_preset = 0;
    app.accent_color = ImVec4(0.26f, 0.59f, 0.98f, 1.0f);
    app.log_count = 0;
    app.viewport_pos = ImVec2(0.0f, 0.0f);
    app.viewport_size = ImVec2(0.0f, 0.0f);
    app.viewport_valid = false;
    app.viewport_zoom = 1.0f;
    app.viewport_pan_x = 0.0f;
    app.viewport_pan_y = 0.0f;
    app.viewport_panning = false;
    app.toolbar_screen_min = ImVec2(0.0f, 0.0f);
    app.toolbar_screen_size = ImVec2(0.0f, 0.0f);
    app.toolbar_valid = false;
    app.pan_start_mouse = ImVec2(0.0f, 0.0f);
    app.pan_start_offset = ImVec2(0.0f, 0.0f);
    app.hud_zoom_value = kDefaultViewportZoom;
    app.hud_zoom_timer = 0.0f;
    app.hud_pan_value = ImVec2(0.0f, 0.0f);
    app.hud_pan_timer = 0.0f;
    app.show_axis_overlay = true;
    app.ruler_mode = false;
    app.ruler_first_point_set = false;
    app.ruler_point_a = ImVec2(0.0f, 0.0f);
    app.ruler_point_b = ImVec2(0.0f, 0.0f);
    app.show_context_menu = false;
    app.context_menu_pos = ImVec2(0.0f, 0.0f);
    app.context_menu_cell_i = 0;
    app.context_menu_cell_j = 0;
    app.window_width = 0;
    app.window_height = 0;
    app.selected_material_id = -1;
    app.paint_material_id = -1;
    app.material_browser_open = false;
    app.material_search[0] = '\0';
    app.filter_metals = true;
    app.filter_dielectrics = true;
    app.show_material_legend = true;
    app.filter_blocks_by_material = false;
    app.auto_filter_blocks_on_select = true;
    app.highlight_blocks_by_material = true;
    app.filter_blocks_by_material = false;
    app.visualization_mode = 0;
    app.show_material_outlines = false;
    app.material_overlay_alpha = 0.5f;
    app.sparam_window_open = false;
    app.sparam_data_loaded = false;
    app.sparam_peak_db = 0.0;
    app.sparam_peak_freq_hz = 0.0;
    app.sparam_f_low_hz = 0.0;
    app.sparam_f_high_hz = 0.0;
    app.sparam_f_center_hz = 0.0;
    app.sparam_bandwidth_hz = 0.0;
    std::strncpy(app.sparam_csv_path, "sweep_results.csv", sizeof(app.sparam_csv_path));
    app.sparam_csv_path[sizeof(app.sparam_csv_path) - 1] = '\0';
    app.smith_chart_open = false;
    app.smith_show_impedance = true;
    app.smith_show_vswr_circles = true;
    app.smith_z0 = 50.0;
    app.smith_freq_index = 0;
    app.request_rebootstrap = false;
    app.rebootstrap_message[0] = '\0';
    app.viewport_layout = VIEWPORT_SINGLE;
    for (int i = 0; i < 4; ++i) {
        app.viewports[i].pos = ImVec2(0.0f, 0.0f);
        app.viewports[i].size = ImVec2(0.0f, 0.0f);
        app.viewports[i].zoom = 1.0f;
        app.viewports[i].pan_x = 0.0f;
        app.viewports[i].pan_y = 0.0f;
        app.viewports[i].viz_mode = VIEWPORT_VIZ_EZ;
        app.viewports[i].active = (i == 0);
        app.viewports[i].valid = (i == 0);
        app.viewports[i].show_grid = true;
        app.viewports[i].show_sources = true;
        app.viewports[i].show_vectors = false;
    }
    app.active_viewport_idx = 0;
    app.sync_zoom = false;
    app.sync_pan = false;
    app.area_mode = false;
    app.annotation_mode = false;
    app.temp_annotation.text[0] = '\0';
    app.temp_annotation.color = ImVec4(1.0f, 1.0f, 0.2f, 1.0f);
    app.temp_annotation.font_size = 14.0f;
    app.temp_annotation.visible = false;
    app.current_area.vertices.clear();
    app.current_area.closed = false;
    app.current_area.area_m2 = 0.0;
    app.current_area.perimeter_m = 0.0;
    app.measurements.distances.clear();
    app.measurements.areas.clear();
    app.measurements.annotations.clear();
    app.show_measurement_history = false;
    app.show_distance_measurements = true;
    app.show_area_measurements = true;
    app.show_annotations = true;
    app.show_print_composer = false;
    app.composer_pages.clear();
    app.composer_active_page = 0;
    app.composer_next_item_id = 1;
    app.composer_status[0] = '\0';
    app.composer_dragging = false;
    app.composer_drag_item = -1;
    app.composer_drag_offset = ImVec2(0, 0);
    app.composer_resizing = false;
    app.composer_resize_item = -1;
    app.composer_resize_corner = -1;
    app.composer_resize_start_pos = ImVec2(0, 0);
    app.composer_resize_start_size = ImVec2(0, 0);
    app.composer_resize_start_mouse = ImVec2(0, 0);
    app.composer_request_export = false;
    app.composer_request_page = 0;
    app.composer_request_export_all = false;
    app.composer_request_animation = false;
    debug_logf("init: sim grid %dx%d", sim->nx, sim->ny);
    ui_log_add(&app, "Simulation started");
    app.viewports[0].viz_mode = VIEWPORT_VIZ_EZ;
    app.viewports[1].viz_mode = VIEWPORT_VIZ_EZ;
    app.viewports[2].viz_mode = VIEWPORT_VIZ_EZ;
    app.viewports[3].viz_mode = VIEWPORT_VIZ_EZ;

    const int min_window_width = 1280;
    const int min_window_height = 720;
    int scale = (int)std::lround(kDefaultViewportZoom);
    int width = sim->nx * scale;
    int height = sim->ny * scale;
    if (width < min_window_width) width = min_window_width;
    if (height < min_window_height) height = min_window_height;
    app.viewport_zoom = kDefaultViewportZoom;
    app.window_width = width;
    app.window_height = height;

    RenderContext* render = render_init("emwave-c (ImGui)", width, height);
    if (!render) {
        std::fprintf(stderr, "Failed to initialise SDL renderer for ImGui front-end\n");
        material_library_shutdown();
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }
    debug_logf("init: render_init ok window %dx%d scale %d", app.window_width, app.window_height, scale);
    SDL_GetWindowSize(render->window, &app.window_width, &app.window_height);
    render->scale = scale;

    // Set up Dear ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable docking
    io.FontGlobalScale = 1.0f;
    io.IniFilename = "imgui.ini";  // Enable layout persistence
    float dpi_scale = 1.0f;
    int display_index = SDL_GetWindowDisplayIndex(render->window);
    float ddpi = 0.0f, hdpi = 0.0f, vdpi = 0.0f;
    if (SDL_GetDisplayDPI(display_index, &ddpi, &hdpi, &vdpi) == 0 && ddpi > 1.0f) {
        dpi_scale = ddpi / 96.0f;
    }
    auto find_font_path = []() -> std::string {
        const char* candidates[] = {
            "third_party/fonts/DejaVuSans.ttf",
            "../third_party/fonts/DejaVuSans.ttf",
            "../../third_party/fonts/DejaVuSans.ttf",
            "assets/fonts/DejaVuSans.ttf"
        };
        for (const char* c : candidates) {
            if (std::filesystem::exists(c)) {
                return std::string(c);
            }
        }
        return std::string();
    };
    std::string font_path = find_font_path();
    float base_font_px = (dpi_scale > 1.4f) ? 18.0f * dpi_scale : 16.0f * dpi_scale;
    io.Fonts->Clear();
    if (!font_path.empty()) {
        io.Fonts->AddFontFromFileTTF(font_path.c_str(), base_font_px);
    } else {
        io.Fonts->AddFontDefault();
    }
    if (dpi_scale > 1.1f) {
        ImGui::GetStyle().ScaleAllSizes(dpi_scale * 0.9f);
    }
    bool layout_initialized = false;
    bool has_saved_layout = (io.IniFilename && std::filesystem::exists(io.IniFilename));

    ImGui_ImplSDL2_InitForSDLRenderer(render->window, render->renderer);
    ImGui_ImplSDLRenderer2_Init(render->renderer);

    // Oscilloscope for a single probe in the center
    Scope scope = {0};
    if (!scope_init(&scope, sim->nx * scale)) {
        std::fprintf(stderr, "Failed to initialise oscilloscope buffer\n");
        ImGui_ImplSDLRenderer2_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImGui::DestroyContext();
        render_free(render);
        material_library_shutdown();
        simulation_bootstrap_shutdown(&bootstrap);
        return 1;
    }

    if (headless.enabled) {
        ensure_composer_initialized(&app);
        if (headless.tmpl_idx == -2 && headless.layout_path[0]) {
            if (!composer_load_layout_json(&app, headless.layout_path)) {
                ui_log_add(&app, "Headless: failed to load layout '%s'", headless.layout_path);
            }
        } else if (headless.tmpl_idx >= 0) {
            composer_apply_template(&app, app.composer_pages[0], headless.tmpl_idx);
        }
        ComposerPage& hp = app.composer_pages[0];
        hp.output_format = headless.fmt;
        hp.fps = std::clamp(headless.fps, 5, 60);
        hp.frames = std::max(1, headless.frames);
        hp.res_w = headless.res_w;
        hp.res_h = headless.res_h;
        std::snprintf(hp.output_name, sizeof(hp.output_name), "%s", headless.output_name);
        // Log summary for headless runs
        int total_items = 0;
        for (const auto& p : app.composer_pages) total_items += (int)p.items.size();
        std::printf("Headless composer: pages=%zu items=%d fmt=%d fps=%d frames=%d res=%dx%d\n",
                    app.composer_pages.size(),
                    total_items,
                    hp.output_format,
                    hp.fps,
                    hp.frames,
                    hp.res_w,
                    hp.res_h);
        bool ok = true;
        if (hp.frames > 1) {
            ok = export_composer_animation(&app, render, sim, &scope, 0, false);
        } else {
            ok = export_composer_page(&app, render, sim, &scope, 0);
        }
        scope_free(&scope);
        ImGui_ImplSDLRenderer2_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        render_free(render);
        material_library_shutdown();
        simulation_bootstrap_shutdown(&bootstrap);
        return ok ? 0 : 1;
    }

    bool running = true;
    bool paused = false;
    int steps_per_frame = 5;      // Speed control
    bool paint_mode = false;
    int paint_material_id = -1;   // Material ID for painting (-1 legacy)
    int paint_material_type = 0;  // Legacy: 0=PEC, 1=PMC, 2=dielectric
    double paint_epsilon = 4.0;   // Legacy dielectric epsilon
    int paint_brush_size = 3;     // Radius in grid cells
    bool show_help_overlay = false;
    bool dragging_source = false;
    int dragged_source_idx = -1;
    bool shift_right_panning = false;
    int panning_viewport_idx = -1;
    LayoutPreset current_layout = LAYOUT_POWER_USER;
    bool layout_dirty = false;
    ImGuiID last_dockspace_id = 0;
    ThemePreset current_theme = THEME_PRESET_DARK_PRO;
    int current_colormap = 0;     // 0=Classic, 1=Viridis, 2=Plasma
    int current_accent = 0;       // 0-5 accent palette index
    SDL_Color viewport_clear_color = theme_viewport_clear_color(current_theme);
    bool auto_rescale = true;
    bool hold_color = false;
    bool hold_scope = false;
    double vmax_smooth = sim->step_Ez_absmax;
    if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
    double scope_vmax = 1.0;
    ViewportLayout last_layout = app.viewport_layout;

    // Apply initial theme, accent, and colormap to both ImGui and SDL renderer
    app.accent_color = accent_from_palette(current_accent);
    apply_theme(current_theme, app.accent_color);
    ThemeMode render_theme = (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
    ui_render_set_theme(render_theme, current_accent);
    ui_render_set_colormap((ColorMapMode)current_colormap);
    app.theme_preset = (int)current_theme;

    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 prev = SDL_GetPerformanceCounter();
    double fps_avg = 0.0;
    int frame_counter = 0;

    while (running) {
        frame_counter++;
        if (frame_counter <= 120) {
            debug_logf("frame %d: viewport_valid=%d size=%.1fx%.1f layout=%d",
                       frame_counter,
                       app.viewport_valid ? 1 : 0,
                       app.viewport_size.x,
                       app.viewport_size.y,
                       (int)app.viewport_layout);
        }
        compute_viewport_layout(&app, app.viewport_size);
        bool layout_switched = (app.viewport_layout != last_layout);
        if (layout_switched && app.viewport_layout == VIEWPORT_QUAD) {
            const ViewportViz quad_defaults[4] = {
                VIEWPORT_VIZ_EZ,     // A
                VIEWPORT_VIZ_HX,     // B
                VIEWPORT_VIZ_HY,     // C
                VIEWPORT_VIZ_S_MAG    // D
            };
            for (int i = 0; i < 4; ++i) {
                if (app.viewports[i].valid) {
                    app.viewports[i].viz_mode = quad_defaults[i];
                }
            }
            app.active_viewport_idx = 0;
        }
        last_layout = app.viewport_layout;
        auto ensure_active_viewport = [&]() -> ViewportInstance* {
            if (app.active_viewport_idx >= 0 &&
                app.active_viewport_idx < 4 &&
                app.viewports[app.active_viewport_idx].valid) {
                return &app.viewports[app.active_viewport_idx];
            }
            for (int i = 0; i < 4; ++i) {
                if (app.viewports[i].valid) {
                    app.active_viewport_idx = i;
                    return &app.viewports[i];
                }
            }
            return nullptr;
        };

        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            ImGui_ImplSDL2_ProcessEvent(&e);
            if (e.type == SDL_QUIT) {
                running = false;
            } else if (e.type == SDL_WINDOWEVENT) {
                if (e.window.event == SDL_WINDOWEVENT_RESIZED ||
                    e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                    int new_w = e.window.data1;
                    int new_h = e.window.data2;
                    if (new_w < min_window_width) new_w = min_window_width;
                    if (new_h < min_window_height) new_h = min_window_height;
                    app.window_width = new_w;
                    app.window_height = new_h;
                    width = new_w;
                    height = new_h;
                } else if (e.window.event == SDL_WINDOWEVENT_MAXIMIZED ||
                           e.window.event == SDL_WINDOWEVENT_RESTORED) {
                    SDL_GetWindowSize(render->window, &width, &height);
                    app.window_width = width;
                    app.window_height = height;
                }
            } else if (e.type == SDL_MOUSEWHEEL) {
                // Zoom with wheel on the hovered viewport (per-viewport).
                int mx, my;
                SDL_GetMouseState(&mx, &my);
                bool over_toolbar = false;
                if (app.toolbar_valid) {
                    float tx0 = app.toolbar_screen_min.x;
                    float ty0 = app.toolbar_screen_min.y;
                    float tx1 = tx0 + app.toolbar_screen_size.x;
                    float ty1 = ty0 + app.toolbar_screen_size.y;
                    over_toolbar = ((float)mx >= tx0 && (float)mx <= tx1 &&
                                    (float)my >= ty0 && (float)my <= ty1);
                }
                int hit = (app.viewport_valid && sim) ? get_viewport_at_mouse(app, mx, my) : -1;
                if (hit >= 0 && !over_toolbar) {
                    app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (vp && vp->valid) {
                        float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                        float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                        if (local_x >= 0.0f && local_y >= 0.0f &&
                            local_x < vp->size.x && local_y < vp->size.y) {
                            int active_scale = (int)std::lround(vp->zoom);
                            if (active_scale < 1) active_scale = 1;
                            float steps = (float)e.wheel.y;
                            if (steps == 0.0f && e.wheel.preciseY != 0.0f) {
                                steps = (float)e.wheel.preciseY;
                            }
                            if (steps != 0.0f) {
                                float base = (steps > 0.0f) ? 1.1f : 0.9f;
                                float mult = std::pow(base, std::fabs(steps));
                                apply_zoom_at_point(vp,
                                                    render,
                                                    sim,
                                                    &active_scale,
                                                    local_x,
                                                    local_y,
                                                    mult,
                                                    true);
                                if (app.sync_zoom) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        app.viewports[i].zoom = vp->zoom;
                                    }
                                }
                                scale = active_scale;
                                app.hud_zoom_value = vp->zoom;
                                app.hud_zoom_timer = 1.0f;
                            }
                        }
                    }
                    continue;
                }
                if (io.WantCaptureMouse) continue;
            } else if (e.type == SDL_KEYDOWN) {
                SDL_Keycode key = e.key.keysym.sym;

                // ============================================================
                // GLOBAL SHORTCUTS - Always work, regardless of ImGui focus
                // ============================================================
                bool handled_globally = false;

                switch (key) {
                    case SDLK_ESCAPE:
                        if (app.area_mode) {
                            app.area_mode = false;
                            app.current_area.vertices.clear();
                            app.current_area.closed = false;
                            ui_log_add(&app, "Area tool: cancelled");
                            handled_globally = true;
                        } else {
                            running = false;
                            handled_globally = true;
                        }
                        break;
                    case SDLK_q:
                        running = false;
                        handled_globally = true;
                        break;

                    case SDLK_F1:
                        show_help_overlay = !show_help_overlay;
                        handled_globally = true;
                        break;

                    case SDLK_F2:
                        current_layout = LAYOUT_BEGINNER;
                        layout_dirty = true;
                        ui_log_add(&app, "Layout: Beginner");
                        handled_globally = true;
                        break;

                    case SDLK_F3:
                        current_layout = LAYOUT_POWER_USER;
                        layout_dirty = true;
                        ui_log_add(&app, "Layout: Power User");
                        handled_globally = true;
                        break;

                    case SDLK_F4:
                        current_layout = LAYOUT_ANALYSIS;
                        layout_dirty = true;
                        ui_log_add(&app, "Layout: Analysis");
                        handled_globally = true;
                        break;

                    case SDLK_F7:
                        current_layout = LAYOUT_CANVAS_FIRST;
                        layout_dirty = true;
                        ui_log_add(&app, "Layout: Canvas First");
                        handled_globally = true;
                        break;

                    case SDLK_1:
                    case SDLK_KP_1:
                        if (SDL_GetModState() & KMOD_ALT) {
                            app.viewport_layout = VIEWPORT_SINGLE;
                            ui_log_add(&app, "Layout: Single viewport");
                            handled_globally = true;
                        }
                        break;
                    case SDLK_2:
                    case SDLK_KP_2:
                        if (SDL_GetModState() & KMOD_ALT) {
                            app.viewport_layout = VIEWPORT_HORIZONTAL;
                            ui_log_add(&app, "Layout: Horizontal split");
                            handled_globally = true;
                        }
                        break;
                    case SDLK_3:
                    case SDLK_KP_3:
                        if (SDL_GetModState() & KMOD_ALT) {
                            app.viewport_layout = VIEWPORT_VERTICAL;
                            ui_log_add(&app, "Layout: Vertical split");
                            handled_globally = true;
                        }
                        break;
                    case SDLK_4:
                    case SDLK_KP_4:
                        if (SDL_GetModState() & KMOD_ALT) {
                            app.viewport_layout = VIEWPORT_QUAD;
                            ui_log_add(&app, "Layout: Quad view");
                            handled_globally = true;
                        }
                        break;

                    case SDLK_F5: {
                        // Load waveguide preset
                        SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                        char errbuf[256];
                        if (config_loader_parse_file("configs/waveguide.json",
                                                     &cfg,
                                                     errbuf,
                                                     sizeof(errbuf))) {
                            if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                wizard_init_from_config(wizard, &cfg);
                                for (int i = 0; i < 4; ++i) {
                                    app.viewports[i].zoom = kDefaultViewportZoom;
                                    app.viewports[i].pan_x = 0.0f;
                                    app.viewports[i].pan_y = 0.0f;
                                }
                                scale = (int)std::lround(kDefaultViewportZoom);
                                if (scale < 1) scale = 1;
                                render->scale = scale;

                                // Resize window to match new grid
                                width = sim->nx * scale;
                                height = sim->ny * scale;
                                if (width < min_window_width) width = min_window_width;
                                if (height < min_window_height) height = min_window_height;
                                SDL_SetWindowSize(render->window, width, height);
                                app.window_width = width;
                                app.window_height = height;

                                // Reset visualization scaling
                                vmax_smooth = sim->step_Ez_absmax;
                                if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                                scope_vmax = 1.0;
                                auto_rescale = true;
                                hold_color = false;
                                hold_scope = false;

                                ui_log_add(&app, "Loaded: waveguide.json");
                            } else {
                                ui_log_add(&app, "Failed to load: waveguide.json");
                            }
                        } else {
                            ui_log_add(&app, "Failed to parse: waveguide.json");
                        }
                        handled_globally = true;
                        break;
                    }

                    case SDLK_F6: {
                        // Load CPW filter preset
                        SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                        char errbuf[256];
                        if (config_loader_parse_file("configs/cpw_filter.json",
                                                     &cfg,
                                                     errbuf,
                                                     sizeof(errbuf))) {
                            if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                                wizard_init_from_config(wizard, &cfg);
                                for (int i = 0; i < 4; ++i) {
                                    app.viewports[i].zoom = kDefaultViewportZoom;
                                    app.viewports[i].pan_x = 0.0f;
                                    app.viewports[i].pan_y = 0.0f;
                                }
                                scale = (int)std::lround(kDefaultViewportZoom);
                                if (scale < 1) scale = 1;
                                render->scale = scale;

                                // Resize window to match new grid
                                width = sim->nx * scale;
                                height = sim->ny * scale;
                                if (width < min_window_width) width = min_window_width;
                                if (height < min_window_height) height = min_window_height;
                                SDL_SetWindowSize(render->window, width, height);
                                app.window_width = width;
                                app.window_height = height;

                                // Reset visualization scaling
                                vmax_smooth = sim->step_Ez_absmax;
                                if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                                scope_vmax = 1.0;
                                auto_rescale = true;
                                hold_color = false;
                                hold_scope = false;

                                ui_log_add(&app, "Loaded: cpw_filter.json");
                            } else {
                                ui_log_add(&app, "Failed to load: cpw_filter.json");
                            }
                        } else {
                            ui_log_add(&app, "Failed to parse: cpw_filter.json");
                        }
                        handled_globally = true;
                        break;
                    }

                    default:
                        break;
                }

                // ============================================================
                // Handle Ctrl+S (also global)
                // ============================================================
                if (!handled_globally &&
                    (SDL_GetModState() & KMOD_CTRL) &&
                    key == SDLK_s) {
                    sync_sim_to_wizard_config(sim, &wizard);
                    const char* save_path = "configs/saved_config.json";
                    if (save_config_json(save_path, &wizard.cfg)) {
                        ui_log_add(&app, "Configuration saved: %s", save_path);
                    } else {
                        ui_log_add(&app, "Failed to save configuration");
                    }
                    handled_globally = true;
                }

                // ============================================================
                // CONTEXT SHORTCUTS - Only work when ImGui doesn't want keyboard
                // ============================================================
                if (!handled_globally && !io.WantCaptureKeyboard) {
                    SDL_Keymod mod = SDL_GetModState();
                    ViewportInstance* active_vp = ensure_active_viewport();
                    int active_scale = active_vp ? (int)std::lround(active_vp->zoom) : scale;
                    if (active_scale < 1) active_scale = 1;
                        switch (key) {
                            case SDLK_EQUALS:
                            case SDLK_KP_PLUS:
                                if (mod & KMOD_CTRL) {
                                    if (app.viewport_valid && sim && active_vp) {
                                        float focus_x = active_vp->size.x * 0.5f;
                                        float focus_y = active_vp->size.y * 0.5f;
                                        apply_zoom_at_point(active_vp,
                                                            render,
                                                            sim,
                                                            &active_scale,
                                                            focus_x,
                                                            focus_y,
                                                            1.2f,
                                                            false);
                                        if (app.sync_zoom) {
                                            for (int i = 0; i < 4; ++i) {
                                                if (!app.viewports[i].valid) continue;
                                                app.viewports[i].zoom = active_vp->zoom;
                                            }
                                        }
                                        scale = active_scale;
                                        app.hud_zoom_value = active_vp->zoom;
                                        app.hud_zoom_timer = 1.0f;
                                    }
                                    break;
                                }
                            break;

                        case SDLK_MINUS:
                        case SDLK_KP_MINUS:
                                if (mod & KMOD_CTRL) {
                                    if (app.viewport_valid && sim && active_vp) {
                                        float focus_x = active_vp->size.x * 0.5f;
                                        float focus_y = active_vp->size.y * 0.5f;
                                        apply_zoom_at_point(active_vp,
                                                            render,
                                                            sim,
                                                            &active_scale,
                                                            focus_x,
                                                            focus_y,
                                                            0.8f,
                                                            false);
                                        if (app.sync_zoom) {
                                            for (int i = 0; i < 4; ++i) {
                                                if (!app.viewports[i].valid) continue;
                                                app.viewports[i].zoom = active_vp->zoom;
                                            }
                                        }
                                        scale = active_scale;
                                        app.hud_zoom_value = active_vp->zoom;
                                        app.hud_zoom_timer = 1.0f;
                                    }
                                    break;
                                }
                            break;

                        case SDLK_LEFT:
                            if (app.viewport_valid && sim && active_vp) {
                                float pan_step = (mod & KMOD_SHIFT) ? 50.0f : 10.0f;
                                active_vp->pan_x += pan_step;
                                clamp_pan_offset(active_vp, sim, active_scale);
                                if (app.sync_pan) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        app.viewports[i].pan_x = active_vp->pan_x;
                                        app.viewports[i].pan_y = active_vp->pan_y;
                                    }
                                }
                                app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                                app.hud_pan_timer = 1.0f;
                            }
                            break;
                        case SDLK_RIGHT:
                            if (app.viewport_valid && sim && active_vp) {
                                float pan_step = (mod & KMOD_SHIFT) ? 50.0f : 10.0f;
                                active_vp->pan_x -= pan_step;
                                clamp_pan_offset(active_vp, sim, active_scale);
                                if (app.sync_pan) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        app.viewports[i].pan_x = active_vp->pan_x;
                                        app.viewports[i].pan_y = active_vp->pan_y;
                                    }
                                }
                                app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                                app.hud_pan_timer = 1.0f;
                            }
                            break;
                        case SDLK_UP:
                            if (app.viewport_valid && sim && active_vp) {
                                float pan_step = (mod & KMOD_SHIFT) ? 50.0f : 10.0f;
                                active_vp->pan_y += pan_step;
                                clamp_pan_offset(active_vp, sim, active_scale);
                                if (app.sync_pan) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        app.viewports[i].pan_x = active_vp->pan_x;
                                        app.viewports[i].pan_y = active_vp->pan_y;
                                    }
                                }
                                app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                                app.hud_pan_timer = 1.0f;
                            }
                            break;
                        case SDLK_DOWN:
                            if (app.viewport_valid && sim && active_vp) {
                                float pan_step = (mod & KMOD_SHIFT) ? 50.0f : 10.0f;
                                active_vp->pan_y -= pan_step;
                                clamp_pan_offset(active_vp, sim, active_scale);
                                if (app.sync_pan) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        app.viewports[i].pan_x = active_vp->pan_x;
                                        app.viewports[i].pan_y = active_vp->pan_y;
                                    }
                                }
                                app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                                app.hud_pan_timer = 1.0f;
                            }
                            break;

                        case SDLK_0:
                        case SDLK_HOME:
                            if ((mod & KMOD_CTRL) || key == SDLK_HOME) {
                                if (app.viewport_valid && sim && active_vp) {
                                    for (int i = 0; i < 4; ++i) {
                                        if (!app.viewports[i].valid) continue;
                                        int tmp_scale = (int)std::lround(kDefaultViewportZoom);
                                        if (tmp_scale < 1) tmp_scale = 1;
                                        reset_view_transform(&app.viewports[i],
                                                             render,
                                                             sim,
                                                             &tmp_scale,
                                                             kDefaultViewportZoom);
                                    }
                                    active_scale = (int)std::lround(kDefaultViewportZoom);
                                    if (active_scale < 1) active_scale = 1;
                                    scale = active_scale;
                                    app.hud_zoom_value = kDefaultViewportZoom;
                                    app.hud_pan_value = ImVec2(0.0f, 0.0f);
                                    app.hud_zoom_timer = 1.0f;
                                    app.hud_pan_timer = 1.0f;
                                }
                                break;
                            }
                            break;

                        // Simulation controls
                        case SDLK_SPACE:
                            paused = !paused;
                            break;

                        // Source toggles
                        case SDLK_1:
                            if (MAX_SRC > 0) sim->sources[0].active = !sim->sources[0].active;
                            break;
                        case SDLK_2:
                            if (MAX_SRC > 1) sim->sources[1].active = !sim->sources[1].active;
                            break;
                        case SDLK_3:
                            if (MAX_SRC > 2) sim->sources[2].active = !sim->sources[2].active;
                            break;

                        // Source type cycling
                        case SDLK_t:
                            if (mod & KMOD_SHIFT) {
                                if (app.viewport_valid) {
                                    int mx, my;
                                    SDL_GetMouseState(&mx, &my);
                                    int hit = get_viewport_at_mouse(app, mx, my);
                                    if (hit >= 0) app.active_viewport_idx = hit;
                                    ViewportInstance* vp = ensure_active_viewport();
                                    if (vp) {
                                        int scale_local = (int)std::lround(vp->zoom);
                                        if (scale_local < 1) scale_local = 1;
                                        float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                                        float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                                        ImVec2 offset = compute_viewport_offset(*vp, sim, scale_local);
                                        float fx = local_x - offset.x;
                                        float fy = local_y - offset.y;
                                        int ix = (int)(fx / (float)scale_local);
                                        int iy = (int)(fy / (float)scale_local);
                                        if (ix >= 0 && iy >= 0 && ix < sim->nx && iy < sim->ny) {
                                            app.annotation_mode = true;
                                            std::snprintf(app.temp_annotation.text,
                                                          sizeof(app.temp_annotation.text),
                                                          "Note (%d,%d)",
                                                          ix,
                                                          iy);
                                            app.temp_annotation.text[sizeof(app.temp_annotation.text) - 1] = '\0';
                                            app.temp_annotation.grid_pos = ImVec2((float)ix, (float)iy);
                                            app.temp_annotation.visible = true;
                                            app.temp_annotation.font_size = 14.0f;
                                            ui_log_add(&app,
                                                       "Annotation: set at (%d, %d)",
                                                       ix,
                                                       iy);
                                        }
                                    }
                                }
                                break;
                            }
                            sources_cycle_type(sim->sources);
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Cycled source type");
                            break;
                        case SDLK_f:
                            sources_cycle_type(sim->sources);
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Cycled source type");
                            break;

                        // Simulation controls
                        case SDLK_r:
                            if (mod & KMOD_CTRL) {
                                fdtd_reset(sim);
                                scope_clear(&scope);
                                ui_log_add(&app, "Simulation reset");
                            } else {
                                app.ruler_mode = !app.ruler_mode;
                                app.ruler_first_point_set = false;
                                ui_log_add(&app, "Ruler tool: %s",
                                           app.ruler_mode ? "ON (click two points)" : "OFF");
                            }
                            break;
                        case SDLK_c:
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "Fields cleared");
                            break;

                        // Auto-rescale modes
                    case SDLK_a:
                            if (mod & KMOD_SHIFT) {
                                app.area_mode = !app.area_mode;
                                app.current_area.vertices.clear();
                                app.current_area.closed = false;
                                if (app.area_mode) {
                                    ui_log_add(&app, "Area tool: click vertices, right-click to close");
                                } else {
                                    ui_log_add(&app, "Area tool: OFF");
                                }
                            } else {
                                auto_rescale = true;
                                hold_color = false;
                                hold_scope = false;
                                ui_log_add(&app, "Rescale mode: Auto (A)");
                            }
                            break;
                        case SDLK_h:
                            auto_rescale = false;
                            hold_color = true;
                            hold_scope = false;
                            ui_log_add(&app, "Rescale mode: Hold Color (H)");
                            break;
                        case SDLK_j:
                            auto_rescale = false;
                            hold_color = false;
                            hold_scope = true;
                            ui_log_add(&app, "Rescale mode: Hold Scope (J)");
                            break;
                        case SDLK_l:
                            auto_rescale = false;
                            hold_color = false;
                            hold_scope = false;
                            ui_log_add(&app, "Rescale mode: Manual (L)");
                            break;

                        // Boundary controls
                        case SDLK_y: {
                            BoundaryType type = boundary_get_type(sim);
                            type = (type == BOUNDARY_CPML) ? BOUNDARY_MUR : BOUNDARY_CPML;
                            boundary_set_type(sim, type);
                            if (boundary_is_cpml_enabled(sim)) {
                                cpml_build_coeffs(sim);
                                cpml_zero_psi(sim);
                            }
                            const char* name = (type == BOUNDARY_CPML) ? "CPML" : "Mur";
                            ui_log_add(&app, "Boundary: %s", name);
                            break;
                        }
                        case SDLK_7:
                        case SDLK_8:
                        case SDLK_9: {
                            int idx = (key == SDLK_7) ? 0 : (key == SDLK_8) ? 1 : 2;
                            boundary_set_type(sim, BOUNDARY_CPML);
                            cpml_apply_preset(sim, idx);
                            if (boundary_is_cpml_enabled(sim)) {
                                cpml_zero_psi(sim);
                            }
                            fdtd_clear_fields(sim);
                            scope_clear(&scope);
                            ui_log_add(&app, "CPML preset %d applied", idx + 1);
                            break;
                        }

                        // S-parameter ports
                        case SDLK_s:
                            sim->ports_on = !sim->ports_on;
                            ui_log_add(&app, "Ports: %s", sim->ports_on ? "ON" : "OFF");
                            break;

                        // Material painting mode
                        case SDLK_m:
                        case SDLK_u:
                            paint_mode = !paint_mode;
                            ui_log_add(&app, "Paint mode: %s", paint_mode ? "ON" : "OFF");
                            break;
                        case SDLK_i:
                            if (paint_mode && paint_material_id < 0) {
                                paint_material_type = (paint_material_type + 1) % 3;
                                const char* mlabel = (paint_material_type == 0)
                                                         ? "PEC"
                                                         : (paint_material_type == 1) ? "PMC"
                                                                                      : "Dielectric";
                                ui_log_add(&app, "Paint type: %s (legacy mode)", mlabel);
                            }
                            break;
                        case SDLK_o:
                            if (paint_mode && paint_material_type == 2 && paint_material_id < 0) {
                                paint_epsilon -= 0.5;
                                if (paint_epsilon < 1.0) paint_epsilon = 1.0;
                                ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            }
                            break;
                        case SDLK_p:
                            if (paint_mode && paint_material_type == 2 && paint_material_id < 0) {
                                paint_epsilon += 0.5;
                                if (paint_epsilon > 20.0) paint_epsilon = 20.0;
                                ui_log_add(&app, "Paint epsilon: %.1f", paint_epsilon);
                            }
                            break;

                        // Visual controls
                        case SDLK_b:
                            current_theme = (ThemePreset)(((int)current_theme + 1) % 4);
                            apply_theme(current_theme, app.accent_color);
                            viewport_clear_color = theme_viewport_clear_color(current_theme);
                            render_theme = (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
                            ui_render_set_theme(render_theme, current_accent);
                            app.theme_preset = (int)current_theme;
                            ui_log_add(&app, "Theme: %s", theme_preset_name(current_theme));
                            break;
                        case SDLK_k: {
                            int shift = (SDL_GetModState() & KMOD_SHIFT) ? -1 : 1;
                            current_colormap =
                                (current_colormap + shift + kColormapCount) % kColormapCount;
                            ui_render_set_colormap((ColorMapMode)current_colormap);
                            ui_log_add(&app, "Colormap: %s",
                                       kColormapNames[current_colormap]);
                            break;
                        }
                        case SDLK_v: {
                            int shift = (SDL_GetModState() & KMOD_SHIFT) ? -1 : 1;
                            current_accent =
                                (current_accent + shift + 6) % 6;
                            app.accent_color = accent_from_palette(current_accent);
                            apply_theme(current_theme, app.accent_color);
                            render_theme = (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
                            ui_render_set_theme(render_theme, current_accent);
                            ui_log_add(&app, "Accent: %s", ACCENT_NAMES[current_accent]);
                            break;
                        }

                        // Grid overlay
                        case SDLK_g:
                            app.show_grid_overlay = !app.show_grid_overlay;
                            ui_log_add(&app, "Grid overlay: %s",
                                       app.show_grid_overlay ? "ON" : "OFF");
                            break;
                        case SDLK_e: {
                            app.visualization_mode = (app.visualization_mode + 1) % 3;
                            const char* mode_names[] = {"Field", "Material", "Overlay"};
                            int vm = app.visualization_mode;
                            if (vm < 0 || vm > 2) vm = 0;
                            ui_log_add(&app, "View mode: %s", mode_names[vm]);
                            break;
                        }

                        case SDLK_KP_1:
                        case SDLK_KP_2:
                        case SDLK_KP_3:
                        case SDLK_KP_4:
                        case SDLK_KP_5:
                        case SDLK_KP_6:
                        case SDLK_KP_7:
                        case SDLK_KP_8:
                            if (paint_mode) {
                                int idx = key - SDLK_KP_1;
                                const Material* mat = material_library_get_by_index(idx);
                                if (mat) {
                                    paint_material_id = mat->id;
                                    app.paint_material_id = paint_material_id;
                                    ui_log_add(&app, "Paint material: %s", mat->name);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_MIDDLE) {
                if (!app.viewport_valid || !sim || io.WantCaptureMouse) {
                    continue;
                }
                int mx = e.button.x;
                int my = e.button.y;
                if (app.toolbar_valid) {
                    float tx0 = app.toolbar_screen_min.x;
                    float ty0 = app.toolbar_screen_min.y;
                    float tx1 = tx0 + app.toolbar_screen_size.x;
                    float ty1 = ty0 + app.toolbar_screen_size.y;
                    if ((float)mx >= tx0 && (float)mx <= tx1 &&
                        (float)my >= ty0 && (float)my <= ty1) {
                        continue;
                    }
                }
                int hit = get_viewport_at_mouse(app, mx, my);
                if (hit < 0) continue;
                app.active_viewport_idx = hit;
                ViewportInstance* vp = ensure_active_viewport();
                if (!vp) continue;
                float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                if (local_x >= 0.0f && local_y >= 0.0f &&
                    local_x < vp->size.x && local_y < vp->size.y) {
                    app.viewport_panning = true;
                    panning_viewport_idx = app.active_viewport_idx;
                    app.pan_start_mouse = ImVec2((float)mx, (float)my);
                    app.pan_start_offset = ImVec2(vp->pan_x, vp->pan_y);
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT) {
                SDL_Keymod btn_mod = SDL_GetModState();
                if (!app.viewport_valid || !sim) {
                    continue;
                }
                if (app.toolbar_valid) {
                    int mx = e.button.x;
                    int my = e.button.y;
                    float tx0 = app.toolbar_screen_min.x;
                    float ty0 = app.toolbar_screen_min.y;
                    float tx1 = tx0 + app.toolbar_screen_size.x;
                    float ty1 = ty0 + app.toolbar_screen_size.y;
                    if ((float)mx >= tx0 && (float)mx <= tx1 &&
                        (float)my >= ty0 && (float)my <= ty1) {
                        continue;  // Let ImGui handle toolbar interactions
                    }
                }
                if (btn_mod & KMOD_SHIFT) {
                    int mx = e.button.x;
                    int my = e.button.y;
                    int hit = get_viewport_at_mouse(app, mx, my);
                    if (hit >= 0) app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (!vp) continue;
                    float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                    if (local_x >= 0.0f && local_y >= 0.0f &&
                        local_x < vp->size.x && local_y < vp->size.y) {
                        shift_right_panning = true;
                        app.viewport_panning = true;
                        panning_viewport_idx = app.active_viewport_idx;
                        app.pan_start_mouse = ImVec2((float)mx, (float)my);
                        app.pan_start_offset = ImVec2(vp->pan_x, vp->pan_y);
                    }
                } else {
                    int mx = e.button.x;
                    int my = e.button.y;
                    int hit = get_viewport_at_mouse(app, mx, my);
                    if (hit >= 0) app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (!vp) continue;
                    int scale_local = (int)std::lround(vp->zoom);
                    if (scale_local < 1) scale_local = 1;
                    float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                    if (local_x >= 0.0f && local_y >= 0.0f &&
                        local_x < vp->size.x && local_y < vp->size.y) {
                        ImVec2 viewport_offset = compute_viewport_offset(*vp, sim, scale_local);
                        float field_x = local_x - viewport_offset.x;
                        float field_y = local_y - viewport_offset.y;
                        int ix = (int)(field_x / (float)scale_local);
                        int iy = (int)(field_y / (float)scale_local);
                        if (ix >= 0 && iy >= 0 && ix < sim->nx && iy < sim->ny) {
                            if (app.area_mode) {
                                if (app.current_area.vertices.size() >= 3) {
                                    close_area_measurement(&app, sim);
                                } else {
                                    ui_log_add(&app, "Area tool: need at least 3 vertices to finish");
                                }
                                continue;
                            }
                            app.show_context_menu = true;
                            app.context_menu_pos = ImVec2((float)mx, (float)my);
                            app.context_menu_cell_i = ix;
                            app.context_menu_cell_j = iy;
                        }
                    }
        }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                if (!app.viewport_valid) {
                    continue;
                }
                if (SDL_GetModState() & KMOD_ALT) {
                    int mx = e.button.x;
                    int my = e.button.y;
                    int hit = get_viewport_at_mouse(app, mx, my);
                    if (hit < 0) continue;
                    app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (!vp) continue;
                    float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                    if (local_x >= 0.0f && local_y >= 0.0f &&
                        local_x < vp->size.x && local_y < vp->size.y) {
                        app.viewport_panning = true;
                        panning_viewport_idx = app.active_viewport_idx;
                        app.pan_start_mouse = ImVec2((float)mx, (float)my);
                        app.pan_start_offset = ImVec2(vp->pan_x, vp->pan_y);
                        continue;
                    }
                }

                int mx = e.button.x;
                int my = e.button.y;
            if (app.toolbar_valid) {
                float tx0 = app.toolbar_screen_min.x;
                float ty0 = app.toolbar_screen_min.y;
                float tx1 = tx0 + app.toolbar_screen_size.x;
                float ty1 = ty0 + app.toolbar_screen_size.y;
                if ((float)mx >= tx0 && (float)mx <= tx1 &&
                    (float)my >= ty0 && (float)my <= ty1) {
                    continue;  // Ignore clicks on the viewport toolbar
                }
            }
            int hit = get_viewport_at_mouse(app, mx, my);
            if (hit < 0) continue;
                app.active_viewport_idx = hit;
                ViewportInstance* vp = ensure_active_viewport();
                if (!vp) continue;
                int scale_local = (int)std::lround(vp->zoom);
                if (scale_local < 1) scale_local = 1;
                float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);

                if (local_x >= 0.0f && local_y >= 0.0f &&
                    local_x < vp->size.x && local_y < vp->size.y) {
                    ImVec2 viewport_offset = compute_viewport_offset(*vp, sim, scale_local);
                    float field_x = local_x - viewport_offset.x;
                    float field_y = local_y - viewport_offset.y;
                    int ix = (int)(field_x / (float)scale_local);
                    int iy = (int)(field_y / (float)scale_local);
                    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) {
                        continue;
                    }
                    if (app.area_mode) {
                        app.current_area.vertices.push_back(ImVec2((float)ix, (float)iy));
                        continue;
                    }
                    if (app.ruler_mode) {
                        ImVec2 grid_pt((float)ix, (float)iy);
                        if (!app.ruler_first_point_set) {
                            app.ruler_point_a = grid_pt;
                            app.ruler_first_point_set = true;
                            ui_log_add(&app,
                                       "Ruler: first point (%d, %d)",
                                       (int)grid_pt.x,
                                       (int)grid_pt.y);
                        } else {
                            app.ruler_point_b = grid_pt;
                            double dx_m = ((double)app.ruler_point_b.x - (double)app.ruler_point_a.x) *
                                          (sim->lx / (double)sim->nx);
                            double dy_m = ((double)app.ruler_point_b.y - (double)app.ruler_point_a.y) *
                                          (sim->ly / (double)sim->ny);
                            double distance_m = std::sqrt(dx_m * dx_m + dy_m * dy_m);
                            double angle_deg = std::atan2(dy_m, dx_m) * 180.0 / IM_PI;
                            ui_log_add(&app,
                                       "Ruler: %.3f m @ %.1f deg",
                                       distance_m,
                                       angle_deg);
                            if (app.show_distance_measurements) {
                                DistanceMeasurement dm;
                                dm.a = app.ruler_point_a;
                                dm.b = app.ruler_point_b;
                                dm.distance_m = distance_m;
                                dm.angle_deg = angle_deg;
                                app.measurements.distances.push_back(dm);
                            }
                            app.ruler_first_point_set = false;
                        }
                        continue;
                    }

                    app.last_click_i = ix;
                    app.last_click_j = iy;

                    // PRIORITY 1: Check for source drag (only if NOT in paint mode)
                    if (!paint_mode) {
                        int hit_source = find_source_at_position(sim, ix, iy, 5);
                        if (hit_source >= 0) {
                            dragging_source = true;
                            dragged_source_idx = hit_source;
                            app.selected_source = hit_source;
                            if (hit_source < wizard.cfg.source_count) {
                                double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
                                double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
                                wizard.cfg.source_configs[hit_source].x =
                                    (double)sim->sources[hit_source].ix / nx1;
                                wizard.cfg.source_configs[hit_source].y =
                                    (double)sim->sources[hit_source].iy / ny1;
                            }
                            continue;  // Don't process other click actions
                        }
                    }

                    // PRIORITY 2: Paint mode (if active)
                    if (paint_mode) {
                        if (paint_material_id >= 0) {
                            apply_paint_with_material(sim, ix, iy, paint_brush_size, paint_material_id);
                        } else {
                            apply_selected_paint(sim,
                                                 ix,
                                                 iy,
                                                 paint_brush_size,
                                                 paint_material_type,
                                                 paint_epsilon);
                        }
                    }

                    // PRIORITY 3: Other interactions (block placement, etc.)
                    else if (app.placing_block) {
                        if (!app.block_first_set) {
                            app.block_first_set = true;
                            app.block_first_i = ix;
                            app.block_first_j = iy;
                        } else {
                            int idx = app.selected_block;
                            if (idx < 0) idx = 0;
                            if (idx >= CONFIG_MAX_MATERIAL_RECTS) {
                                idx = CONFIG_MAX_MATERIAL_RECTS - 1;
                            }
                            if (wizard.cfg.material_rect_count <= idx) {
                                wizard.cfg.material_rect_count = idx + 1;
                                if (wizard.cfg.material_rect_count >
                                    CONFIG_MAX_MATERIAL_RECTS) {
                                    wizard.cfg.material_rect_count =
                                        CONFIG_MAX_MATERIAL_RECTS;
                                }
                            }
                            MaterialRectSpec& r = wizard.cfg.material_rects[idx];
                            int i0 = (app.block_first_i < ix) ? app.block_first_i : ix;
                            int i1 = (app.block_first_i > ix) ? app.block_first_i : ix;
                            int j0 = (app.block_first_j < iy) ? app.block_first_j : iy;
                            int j1 = (app.block_first_j > iy) ? app.block_first_j : iy;
                            if (i0 < 0) i0 = 0;
                            if (i1 >= sim->nx) i1 = sim->nx - 1;
                            if (j0 < 0) j0 = 0;
                            if (j1 >= sim->ny) j1 = sim->ny - 1;
                            double nx_d = (double)sim->nx;
                            double ny_d = (double)sim->ny;
                            r.x0 = (double)i0 / nx_d;
                            r.x1 = (double)(i1 + 1) / nx_d;
                            r.y0 = (double)j0 / ny_d;
                            r.y1 = (double)(j1 + 1) / ny_d;
                            app.block_first_set = false;

                            apply_wizard_materials_to_sim(wizard, &bootstrap, sim);
                        }
                    }
                }
            } else if (e.type == SDL_MOUSEMOTION) {
                if (app.viewport_panning) {
                    if (!app.viewport_valid || !sim) {
                        continue;
                    }
                    ViewportInstance* vp = nullptr;
                    if (panning_viewport_idx >= 0 &&
                        panning_viewport_idx < 4 &&
                        app.viewports[panning_viewport_idx].valid) {
                        vp = &app.viewports[panning_viewport_idx];
                    } else {
                        vp = ensure_active_viewport();
                    }
                    if (!vp) continue;
                    int scale_local = (int)std::lround(vp->zoom);
                    if (scale_local < 1) scale_local = 1;
                    float delta_x = (float)e.motion.x - app.pan_start_mouse.x;
                    float delta_y = (float)e.motion.y - app.pan_start_mouse.y;
                    vp->pan_x = app.pan_start_offset.x + delta_x;
                    vp->pan_y = app.pan_start_offset.y + delta_y;
                    clamp_pan_offset(vp, sim, scale_local);
                    if (app.sync_pan) {
                        for (int i = 0; i < 4; ++i) {
                            if (!app.viewports[i].valid) continue;
                            if (i == panning_viewport_idx) continue;
                            app.viewports[i].pan_x = vp->pan_x;
                            app.viewports[i].pan_y = vp->pan_y;
                        }
                    }
                    app.hud_pan_value = ImVec2(vp->pan_x, vp->pan_y);
                    app.hud_pan_timer = 1.0f;
                    continue;
                }
                // Handle source dragging first
                if (dragging_source &&
                    dragged_source_idx >= 0 &&
                    dragged_source_idx < MAX_SRC) {
                    if (!app.viewport_valid) {
                        continue;
                    }

                    int mx = e.motion.x;
                    int my = e.motion.y;
                    int hit = get_viewport_at_mouse(app, mx, my);
                    if (hit >= 0) app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (!vp) continue;
                    int scale_local = (int)std::lround(vp->zoom);
                    if (scale_local < 1) scale_local = 1;
                    float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);

                    if (local_x < 0.0f || local_y < 0.0f ||
                        local_x >= vp->size.x || local_y >= vp->size.y) {
                        continue;
                    }

                    ImVec2 viewport_offset = compute_viewport_offset(*vp, sim, scale_local);
                    float field_x = local_x - viewport_offset.x;
                    float field_y = local_y - viewport_offset.y;
                    int ix = (int)(field_x / (float)scale_local);
                    int iy = (int)(field_y / (float)scale_local);

                    // Clamp to valid grid range (1 cell margin)
                    if (ix < 1) ix = 1;
                    if (ix >= sim->nx - 1) ix = sim->nx - 2;
                    if (iy < 1) iy = 1;
                    if (iy >= sim->ny - 1) iy = sim->ny - 2;

                    sim->sources[dragged_source_idx].ix = ix;
                    sim->sources[dragged_source_idx].iy = iy;

                    if (dragged_source_idx < wizard.cfg.source_count) {
                        double nx1 = (sim->nx > 1) ? (double)(sim->nx - 1) : 1.0;
                        double ny1 = (sim->ny > 1) ? (double)(sim->ny - 1) : 1.0;
                        wizard.cfg.source_configs[dragged_source_idx].x =
                            (double)ix / nx1;
                        wizard.cfg.source_configs[dragged_source_idx].y =
                            (double)iy / ny1;
                    }
                }
                // Then handle paint drag
                else if (paint_mode &&
                         (e.motion.state & SDL_BUTTON_LMASK)) {
                    if (!app.viewport_valid) {
                        continue;
                    }
                    int mx = e.motion.x;
                    int my = e.motion.y;
                    int hit = get_viewport_at_mouse(app, mx, my);
                    if (hit >= 0) app.active_viewport_idx = hit;
                    ViewportInstance* vp = ensure_active_viewport();
                    if (!vp) continue;
                    int scale_local = (int)std::lround(vp->zoom);
                    if (scale_local < 1) scale_local = 1;
                    float local_x = (float)mx - (app.viewport_pos.x + vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + vp->pos.y);
                    if (local_x < 0.0f || local_y < 0.0f) {
                        continue;
                    }
                    if (local_x >= vp->size.x || local_y >= vp->size.y) {
                        continue;
                    }
                    ImVec2 viewport_offset = compute_viewport_offset(*vp, sim, scale_local);
                    float field_x = local_x - viewport_offset.x;
                    float field_y = local_y - viewport_offset.y;
                    int ix = (int)(field_x / (float)scale_local);
                    int iy = (int)(field_y / (float)scale_local);
                    if (ix < 0 || iy < 0 || ix >= sim->nx || iy >= sim->ny) {
                        continue;
                    }
                    if (paint_material_id >= 0) {
                        apply_paint_with_material(sim, ix, iy, paint_brush_size, paint_material_id);
                    } else {
                        apply_selected_paint(sim,
                                             ix,
                                             iy,
                                             paint_brush_size,
                                             paint_material_type,
                                             paint_epsilon);
                    }
                }
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_MIDDLE) {
                if (app.viewport_panning) {
                    app.viewport_panning = false;
                    shift_right_panning = false;
                    panning_viewport_idx = -1;
                }
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_RIGHT) {
                if (shift_right_panning || app.viewport_panning) {
                    app.viewport_panning = false;
                    shift_right_panning = false;
                    panning_viewport_idx = -1;
                }
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                if (app.viewport_panning) {
                    app.viewport_panning = false;
                    shift_right_panning = false;
                    panning_viewport_idx = -1;
                }

                if (dragging_source) {
                    if (dragged_source_idx >= 0 && dragged_source_idx < MAX_SRC) {
                        ui_log_add(&app,
                                   "Source %d moved to (%d, %d)",
                                   dragged_source_idx,
                                   sim->sources[dragged_source_idx].ix,
                                   sim->sources[dragged_source_idx].iy);
                    }
                    dragging_source = false;
                    dragged_source_idx = -1;
                }
            }
        }

        // Tie wizard advanced flag to app-level basic/advanced mode
        wizard.advanced = !app.basic_mode;

        ImGui_ImplSDL2_NewFrame();
        ImGui_ImplSDLRenderer2_NewFrame();
        ImGui::NewFrame();
        if (frame_counter <= 120) {
            debug_logf("frame %d: after ImGui NewFrame", frame_counter);
        }

        ImGuiViewport* main_viewport = ImGui::GetMainViewport();
        const float status_bar_h = 26.0f;
        ImVec2 host_pos = main_viewport->Pos;
        ImVec2 host_size = ImVec2(main_viewport->Size.x, main_viewport->Size.y - status_bar_h);
        if (host_size.y < 0.0f) host_size.y = 0.0f;
        ImGui::SetNextWindowPos(host_pos);
        ImGui::SetNextWindowSize(host_size);
        ImGui::SetNextWindowViewport(main_viewport->ID);
        ImGuiWindowFlags host_flags = ImGuiWindowFlags_NoDocking |
                                      ImGuiWindowFlags_NoTitleBar |
                                      ImGuiWindowFlags_NoCollapse |
                                      ImGuiWindowFlags_NoResize |
                                      ImGuiWindowFlags_NoMove |
                                      ImGuiWindowFlags_NoBringToFrontOnFocus |
                                      ImGuiWindowFlags_NoNavFocus |
                                      ImGuiWindowFlags_NoBackground |
                                      ImGuiWindowFlags_MenuBar;
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("DockSpaceHost", nullptr, host_flags);
        ImGui::PopStyleVar(3);

        ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
        last_dockspace_id = dockspace_id;

        if (!layout_initialized && !has_saved_layout) {
            layout_dirty = true;
        }
        if (layout_dirty) {
            apply_layout_preset(dockspace_id, current_layout);
            layout_dirty = false;
            layout_initialized = true;
        }

        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_None);

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Save Configuration...", "Ctrl+S")) {
                    sync_sim_to_wizard_config(sim, &wizard);
                    const char* save_path = "configs/saved_config.json";
                    if (save_config_json(save_path, &wizard.cfg)) {
                        ui_log_add(&app, "Configuration saved: %s", save_path);
                    } else {
                        ui_log_add(&app, "Failed to save configuration");
                    }
                }

                ImGui::Separator();

                if (ImGui::MenuItem("Load Waveguide Preset", "F5")) {
                    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                    char errbuf[256];
                    if (config_loader_parse_file("configs/waveguide.json", &cfg, errbuf,
                                                 sizeof(errbuf))) {
                        if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                            wizard_init_from_config(wizard, &cfg);
                            for (int i = 0; i < 4; ++i) {
                                int tmp_scale = (int)std::lround(kDefaultViewportZoom);
                                if (tmp_scale < 1) tmp_scale = 1;
                                reset_view_transform(&app.viewports[i],
                                                     render,
                                                     sim,
                                                     &tmp_scale,
                                                     kDefaultViewportZoom);
                            }
                            scale = (int)std::lround(kDefaultViewportZoom);
                            if (scale < 1) scale = 1;

                            width = sim->nx * scale;
                            height = sim->ny * scale;
                            if (width < min_window_width) width = min_window_width;
                            if (height < min_window_height) height = min_window_height;
                            SDL_SetWindowSize(render->window, width, height);
                            app.window_width = width;
                            app.window_height = height;

                            vmax_smooth = sim->step_Ez_absmax;
                            if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                            scope_vmax = 1.0;
                            auto_rescale = true;
                            hold_color = false;
                            hold_scope = false;

                            ui_log_add(&app, "Loaded: waveguide.json");
                        } else {
                            ui_log_add(&app, "Failed to load: waveguide.json");
                        }
                    } else {
                        ui_log_add(&app, "Failed to parse: waveguide.json");
                    }
                }

                if (ImGui::MenuItem("Load CPW Filter Preset", "F6")) {
                    SimulationConfig cfg = SIM_CONFIG_DEFAULTS;
                    char errbuf[256];
                    if (config_loader_parse_file("configs/cpw_filter.json", &cfg, errbuf,
                                                 sizeof(errbuf))) {
                        if (rebootstrap_simulation(&cfg, &bootstrap, &sim, &scope, scale)) {
                            wizard_init_from_config(wizard, &cfg);
                            for (int i = 0; i < 4; ++i) {
                                int tmp_scale = (int)std::lround(kDefaultViewportZoom);
                                if (tmp_scale < 1) tmp_scale = 1;
                                reset_view_transform(&app.viewports[i],
                                                     render,
                                                     sim,
                                                     &tmp_scale,
                                                     kDefaultViewportZoom);
                            }
                            scale = (int)std::lround(kDefaultViewportZoom);
                            if (scale < 1) scale = 1;

                            width = sim->nx * scale;
                            height = sim->ny * scale;
                            if (width < min_window_width) width = min_window_width;
                            if (height < min_window_height) height = min_window_height;
                            SDL_SetWindowSize(render->window, width, height);
                            app.window_width = width;
                            app.window_height = height;

                            vmax_smooth = sim->step_Ez_absmax;
                            if (vmax_smooth <= 0.0) vmax_smooth = 1.0;
                            scope_vmax = 1.0;
                            auto_rescale = true;
                            hold_color = false;
                            hold_scope = false;

                            ui_log_add(&app, "Loaded: cpw_filter.json");
                        } else {
                            ui_log_add(&app, "Failed to load: cpw_filter.json");
                        }
                    } else {
                        ui_log_add(&app, "Failed to parse: cpw_filter.json");
                    }
                }

                ImGui::Separator();

                if (ImGui::MenuItem("Quit", "ESC")) {
                    running = false;
                }

                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("View")) {
                ImGui::MenuItem("Grid & Domain", nullptr, &app.show_grid_panel);
                ImGui::MenuItem("Scenes", nullptr, &app.show_scene_panel);
                ImGui::MenuItem("Sources", nullptr, &app.show_sources_panel);
                ImGui::MenuItem("Blocks", nullptr, &app.show_blocks_panel);
                ImGui::MenuItem("Material Legend", nullptr, &app.show_material_legend);
                ImGui::MenuItem("Run Controls", nullptr, &app.show_run_panel);
                ImGui::MenuItem("Run Settings", nullptr, &app.show_run_settings_panel);
                ImGui::MenuItem("Probes", nullptr, &app.show_probes_panel);
                ImGui::MenuItem("Scope / FFT", nullptr, &app.show_scope_window);
                ImGui::MenuItem("Log", nullptr, &app.show_log_panel);
                ImGui::Separator();
                ImGui::MenuItem("Grid Overlay", nullptr, &app.show_grid_overlay);
                ImGui::MenuItem("Axis Overlay", nullptr, &app.show_axis_overlay);
                if (ImGui::BeginMenu("Theme")) {
                    auto set_theme = [&](ThemePreset preset) {
                        current_theme = preset;
                        app.theme_preset = (int)preset;
                        apply_theme(current_theme, app.accent_color);
                        viewport_clear_color = theme_viewport_clear_color(current_theme);
                        render_theme = (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
                        ui_render_set_theme(render_theme, current_accent);
                    };
                    if (ImGui::MenuItem("Dark Professional", nullptr, current_theme == THEME_PRESET_DARK_PRO)) {
                        set_theme(THEME_PRESET_DARK_PRO);
                    }
                    if (ImGui::MenuItem("Blender Inspired", nullptr, current_theme == THEME_PRESET_BLENDER)) {
                        set_theme(THEME_PRESET_BLENDER);
                    }
                    if (ImGui::MenuItem("Light Mode", nullptr, current_theme == THEME_PRESET_LIGHT)) {
                        set_theme(THEME_PRESET_LIGHT);
                    }
                    if (ImGui::MenuItem("High Contrast", nullptr, current_theme == THEME_PRESET_HIGH_CONTRAST)) {
                        set_theme(THEME_PRESET_HIGH_CONTRAST);
                    }
                    ImGui::Separator();
                    if (ImGui::BeginMenu("Accent Palette")) {
                        for (int i = 0; i < 6; ++i) {
                            if (ImGui::MenuItem(ACCENT_NAMES[i], nullptr, current_accent == i)) {
                                current_accent = i;
                                app.accent_color = accent_from_palette(current_accent);
                                apply_theme(current_theme, app.accent_color);
                                render_theme =
                                    (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
                                ui_render_set_theme(render_theme, current_accent);
                            }
                        }
                        ImGui::EndMenu();
                    }
                    if (ImGui::ColorEdit3("Accent Color",
                                          (float*)&app.accent_color,
                                          ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_DisplayRGB)) {
                        apply_theme(current_theme, app.accent_color);
                        render_theme = (current_theme == THEME_PRESET_LIGHT) ? THEME_LIGHT : THEME_DARK;
                        ui_render_set_theme(render_theme, current_accent);
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Layout")) {
                if (ImGui::MenuItem("Beginner", "F2")) {
                    current_layout = LAYOUT_BEGINNER;
                    layout_dirty = true;
                }
                if (ImGui::MenuItem("Power User", "F3")) {
                    current_layout = LAYOUT_POWER_USER;
                    layout_dirty = true;
                }
                if (ImGui::MenuItem("Analysis", "F4")) {
                    current_layout = LAYOUT_ANALYSIS;
                    layout_dirty = true;
                }
                if (ImGui::MenuItem("Canvas First", "F7")) {
                    current_layout = LAYOUT_CANVAS_FIRST;
                    layout_dirty = true;
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Reset Layout")) {
                    layout_dirty = true;
                }
                if (ImGui::MenuItem("Save Layout...")) {
                    ImGui::SaveIniSettingsToDisk("custom_layout.ini");
                    ui_log_add(&app, "Layout saved: custom_layout.ini");
                }
                if (ImGui::MenuItem("Load Layout...")) {
                    ImGui::LoadIniSettingsFromDisk("custom_layout.ini");
                    ui_log_add(&app, "Layout loaded: custom_layout.ini");
                    layout_dirty = true;
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Materials")) {
                if (ImGui::MenuItem("Material Browser...", "M")) {
                    app.material_browser_open = true;
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Analysis")) {
                if (ImGui::MenuItem("S-Parameters...")) {
                    app.sparam_window_open = true;
                }
                if (ImGui::MenuItem("Smith Chart...")) {
                    app.smith_chart_open = true;
                }
                if (ImGui::MenuItem("Measurement History...", nullptr, &app.show_measurement_history)) {
                }
                if (paused) {
                    ImGui::Separator();
                    if (ImGui::MenuItem("Save Scope Snapshot (PNG)")) {
                        std::filesystem::create_directories("recordings");
                        SDL_Surface* surf = composer_render_scope_surface(&app, render, &scope, 800, 480);
                        if (surf) {
                            time_t now = time(nullptr);
                            struct tm* tm_info = localtime(&now);
                            char fname[128];
                            std::snprintf(fname,
                                          sizeof(fname),
                                          "scope_%04d%02d%02d_%02d%02d%02d.png",
                                          tm_info->tm_year + 1900,
                                          tm_info->tm_mon + 1,
                                          tm_info->tm_mday,
                                          tm_info->tm_hour,
                                          tm_info->tm_min,
                                          tm_info->tm_sec);
                            std::filesystem::path outp = std::filesystem::path("recordings") / fname;
                            if (save_surface_png(outp, surf)) {
                                std::snprintf(app.composer_status, sizeof(app.composer_status), "Saved scope: %s", outp.string().c_str());
                                std::snprintf(app.composer_last_export_path, sizeof(app.composer_last_export_path), "%s", outp.parent_path().string().c_str());
                                ui_log_add(&app, "Saved scope snapshot: %s", outp.string().c_str());
                            }
                            SDL_FreeSurface(surf);
                        }
                    }
                    if (ImGui::MenuItem("Save FFT Snapshot (PNG)")) {
                        std::filesystem::create_directories("recordings");
                        SDL_Surface* surf = composer_render_fft_surface(&app, render, sim, &scope, 800, 480);
                        if (surf) {
                            time_t now = time(nullptr);
                            struct tm* tm_info = localtime(&now);
                            char fname[128];
                            std::snprintf(fname,
                                          sizeof(fname),
                                          "fft_%04d%02d%02d_%02d%02d%02d.png",
                                          tm_info->tm_year + 1900,
                                          tm_info->tm_mon + 1,
                                          tm_info->tm_mday,
                                          tm_info->tm_hour,
                                          tm_info->tm_min,
                                          tm_info->tm_sec);
                            std::filesystem::path outp = std::filesystem::path("recordings") / fname;
                            if (save_surface_png(outp, surf)) {
                                std::snprintf(app.composer_status, sizeof(app.composer_status), "Saved FFT: %s", outp.string().c_str());
                                std::snprintf(app.composer_last_export_path, sizeof(app.composer_last_export_path), "%s", outp.parent_path().string().c_str());
                                ui_log_add(&app, "Saved FFT snapshot: %s", outp.string().c_str());
                            }
                            SDL_FreeSurface(surf);
                        }
                    }
                    if (ImGui::MenuItem("Save Legend Snapshot (PNG)")) {
                        std::filesystem::create_directories("recordings");
                        SDL_Surface* surf = composer_render_legend_surface(render, 400, 300);
                        if (surf) {
                            time_t now = time(nullptr);
                            struct tm* tm_info = localtime(&now);
                            char fname[128];
                            std::snprintf(fname,
                                          sizeof(fname),
                                          "legend_%04d%02d%02d_%02d%02d%02d.png",
                                          tm_info->tm_year + 1900,
                                          tm_info->tm_mon + 1,
                                          tm_info->tm_mday,
                                          tm_info->tm_hour,
                                          tm_info->tm_min,
                                          tm_info->tm_sec);
                            std::filesystem::path outp = std::filesystem::path("recordings") / fname;
                            if (save_surface_png(outp, surf)) {
                                std::snprintf(app.composer_status, sizeof(app.composer_status), "Saved legend: %s", outp.string().c_str());
                                std::snprintf(app.composer_last_export_path, sizeof(app.composer_last_export_path), "%s", outp.parent_path().string().c_str());
                                ui_log_add(&app, "Saved legend snapshot: %s", outp.string().c_str());
                            }
                            SDL_FreeSurface(surf);
                        }
                    }
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Tools")) {
                ImGui::MenuItem("Print Composer...", nullptr, &app.show_print_composer);
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        ImGui::End();

        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGuiWindowFlags viewport_flags = ImGuiWindowFlags_NoScrollbar |
                                          ImGuiWindowFlags_NoScrollWithMouse |
                                          ImGuiWindowFlags_NoCollapse |
                                          ImGuiWindowFlags_NoSavedSettings |
                                          ImGuiWindowFlags_NoBackground;
        ImGui::SetNextWindowDockID(dockspace_id, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Viewport", nullptr, viewport_flags)) {
            ImVec2 win_pos = ImGui::GetWindowPos();
            ImVec2 content_min = ImGui::GetWindowContentRegionMin();
            ImVec2 content_max = ImGui::GetWindowContentRegionMax();
            ImVec2 content_pos = ImVec2(win_pos.x + content_min.x, win_pos.y + content_min.y);
            ImVec2 avail = ImVec2(content_max.x - content_min.x, content_max.y - content_min.y);
            app.viewport_pos = content_pos;
            app.viewport_size = avail;
            app.viewport_valid = (avail.x > 4.0f && avail.y > 4.0f);
        } else {
            app.viewport_valid = false;
        }
        ImGui::End();
        ImGui::PopStyleVar();

        // If docking collapses the viewport window to nearly zero, fall back to host area
        if (!app.viewport_valid) {
            float menu_h = ImGui::GetFrameHeight();
            app.viewport_pos = ImVec2(host_pos.x, host_pos.y + menu_h);
            app.viewport_size = ImVec2(host_size.x, host_size.y - menu_h);
            app.viewport_valid = (app.viewport_size.x > 4.0f && app.viewport_size.y > 4.0f);
        }
        compute_viewport_layout(&app, app.viewport_size);
        if (frame_counter <= 120) {
            debug_logf("frame %d: viewport window pos %.1f,%.1f size %.1f,%.1f valid=%d",
                       frame_counter,
                       app.viewport_pos.x,
                       app.viewport_pos.y,
                       app.viewport_size.x,
                       app.viewport_size.y,
                       app.viewport_valid ? 1 : 0);
        }

        if (app.show_grid_panel) {
            if (ImGui::Begin("Grid & Domain", &app.show_grid_panel)) {
                draw_grid_panel(&wizard, &app);
            }
            ImGui::End();
        }
        if (app.show_scene_panel) {
            if (ImGui::Begin("Scenes", &app.show_scene_panel)) {
                draw_scene_panel(sim, wizard);
            }
            ImGui::End();
        }
        if (app.show_sources_panel) {
            if (ImGui::Begin("Sources", &app.show_sources_panel)) {
                draw_sources_panel(sim, wizard, &app);
            }
            ImGui::End();
        }
        if (app.show_blocks_panel) {
            if (ImGui::Begin("Materials / Blocks", &app.show_blocks_panel)) {
                draw_blocks_panel(wizard, &bootstrap, sim, &app);
            }
            ImGui::End();
        }
        if (app.show_material_legend) {
            if (ImGui::Begin("Material Legend", &app.show_material_legend)) {
                draw_material_legend(&app, sim, &wizard, &bootstrap);
            }
            ImGui::End();
        }
        if (app.show_probes_panel) {
            if (ImGui::Begin("Probes", &app.show_probes_panel)) {
                draw_probes_panel(sim);
            }
            ImGui::End();
        }
        if (app.show_run_panel) {
            if (ImGui::Begin("Run Controls", &app.show_run_panel)) {
                ImGui::SetNextItemOpen(false, ImGuiCond_Once);
                if (ImGui::CollapsingHeader("Simulation Controls")) {
                    float right_w = ImGui::GetContentRegionAvail().x;
                    ImVec2 button_size(ImMax(0.0f, right_w * 0.45f), 0.0f);
                    if (ImGui::Button(paused ? "Resume" : "Pause", button_size)) {
                        paused = !paused;
                    }
                    ImGui::SameLine();
                    ImGui::Text("Space");

                    ImGui::Separator();

                    ImGui::TextUnformatted("Frequency");
                    ImGui::Indent();
                    float freq_ghz = (float)(sim->freq * 1e-9);
                    if (ImGui::SliderFloat("GHz", &freq_ghz, 0.001f, 5.0f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                        double new_freq = (double)freq_ghz * 1e9;
                        if (new_freq > 0.0) {
                            fdtd_update_grid_for_freq(sim, new_freq);
                            sources_set_freq(sim->sources, new_freq);
                        }
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Adjust simulation frequency (1 MHz - 5 GHz)");
                    }
                    ImGui::Unindent();

                    ImGui::Separator();

                    ImGui::TextUnformatted("Speed");
                    ImGui::Indent();
                    if (ImGui::SliderInt("steps/frame", &steps_per_frame, 1, 50)) {
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Number of simulation steps per rendered frame");
                    }
                    ImGui::Unindent();

                    ImGui::Separator();
                    ImGui::TextUnformatted("Paint Mode");
                    ImGui::Indent();
                    ImGui::Text("Active: %s (M/U)", paint_mode ? "YES" : "NO");
                    if (paint_mode) {
                        if (paint_material_id >= 0) {
                            const Material* mat = material_library_get_by_id(paint_material_id);
                            if (mat) {
                                ImVec4 color = ImVec4(mat->color_r / 255.0f,
                                                      mat->color_g / 255.0f,
                                                      mat->color_b / 255.0f,
                                                      1.0f);
                                ImGui::PushStyleColor(ImGuiCol_Text, color);
                                ImGui::Text("Material: %s", mat->name);
                                ImGui::PopStyleColor();
                                if (mat->type == MAT_TYPE_PEC) {
                                    ImGui::TextDisabled("(Perfect conductor)");
                                } else if (mat->type == MAT_TYPE_PMC) {
                                    ImGui::TextDisabled("(Magnetic conductor)");
                                } else {
                                    ImGui::Text("eps_r = %.2f", mat->epsilon_r);
                                    if (mat->tan_delta > 0.0) {
                                        ImGui::Text("tan(delta) = %.4f", mat->tan_delta);
                                    }
                                }
                            }
                        } else {
                            const char* type_label = (paint_material_type == 0)
                                                         ? "PEC"
                                                         : (paint_material_type == 1) ? "PMC"
                                                                                       : "Dielectric";
                            ImGui::Text("Type: %s", type_label);
                            if (paint_material_type == 2) {
                                ImGui::Text("eps_r = %.2f", paint_epsilon);
                            }
                        }
                        ImGui::Text("Brush radius: %d", paint_brush_size);
                    }
                    ImGui::Unindent();

                    ImGui::Separator();

                    ImGui::TextUnformatted("Status");
                    ImGui::Indent();
                    ImGui::Text("Step: %d", sim->timestep);
                    ImGui::Text("FPS: %.1f", fps_avg);
                    ImGui::Unindent();

                    if (ImGui::TreeNode("Grid Info")) {
                        ImGui::Text("Size: %d x %d", sim->nx, sim->ny);
                        ImGui::Text("Domain: %.3f x %.3f m", sim->lx, sim->ly);
                        ImGui::Text("dt: %.3e s", sim->dt);
                        ImGui::Text("freq: %.3f GHz", sim->freq * 1e-9);
                        ImGui::TreePop();
                    }

                    ImGui::Separator();
                    ImGui::TextUnformatted("Visualization Scale");
                    ImGui::Indent();
                    const char* mode_str = "Manual (L)";
                    if (auto_rescale && !hold_color && !hold_scope) {
                        mode_str = "Auto (A)";
                    } else if (hold_color) {
                        mode_str = "Hold Color (H)";
                    } else if (hold_scope) {
                        mode_str = "Hold Scope (J)";
                    }
                    ImGui::Text("Mode: %s", mode_str);
                    ImGui::Text("Field vmax: %.3e", vmax_smooth);
                    ImGui::Text("Scope vmax: %.3e", scope_vmax);
                    ImGui::Unindent();

                    ImGui::Separator();
                    ImGui::Checkbox("Basic mode", &app.basic_mode);

                    ImGui::Separator();
                    if (ImGui::Button(app.show_grid_overlay ? "Hide grid" : "Show grid")) {
                        app.show_grid_overlay = !app.show_grid_overlay;
                    }
                }
            }
            ImGui::End();
        }
        if (app.show_run_settings_panel) {
            if (ImGui::Begin("Run Settings", &app.show_run_settings_panel)) {
                draw_run_settings_panel(&wizard, &app);
            }
            ImGui::End();
        }

        if (app.show_scope_window) {
            if (ImGui::Begin("Scope", &app.show_scope_window)) {
                if (ImGui::BeginTabBar("ScopeTabs", ImGuiTabBarFlags_None)) {
                    if (ImGui::BeginTabItem("Scope")) {
                        if (scope.y && scope.n > 0) {
                            const int N = scope.n;
                            static float values[1024];
                            int max_samples = (N <= (int)(sizeof(values) / sizeof(values[0]))) ? N : (int)(sizeof(values) / sizeof(values[0]));
                            double vmin = 0.0, vmax_scope = 0.0;
                            bool first = true;
                            for (int i = 0; i < max_samples; ++i) {
                                int idx = (scope.head + i) % N;
                                double v = scope.y[idx];
                                values[i] = (float)v;
                                if (first) {
                                    vmin = vmax_scope = v;
                                    first = false;
                                } else {
                                    if (v < vmin) vmin = v;
                                    if (v > vmax_scope) vmax_scope = v;
                                }
                            }

                            float ymin = 0.0f;
                            float ymax = 0.0f;
                            if (scope_vmax > 0.0) {
                                ymax = (float)scope_vmax;
                                ymin = -ymax;
                            } else {
                                ymin = (float)vmin;
                                ymax = (float)vmax_scope;
                                if (ymin == ymax) {
                                    ymin -= 1.0f;
                                    ymax += 1.0f;
                                }
                            }

                            ImVec2 plot_size = ImGui::GetContentRegionAvail();
                            if (plot_size.y < 120.0f) plot_size.y = 120.0f;
                            ImGui::PlotLines("Probe Ez(center)", values, max_samples, 0, nullptr, ymin, ymax, plot_size);
                        } else {
                            ImGui::TextUnformatted("Scope: no data available");
                        }
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("FFT")) {
                        if (scope.n > 0 && sim) {
                            static bool use_db_scale = true;
                            ImGui::Checkbox("dB Scale", &use_db_scale);
                            static double fft_freq[1024];
                            static double fft_mag[1024];
                            static double fft_mag_db[1024];
                            static double fft_phase[1024];
                            int fft_points =
                                compute_scope_fft(&scope, sim->dt, fft_freq, fft_mag, fft_phase,
                                                  IM_ARRAYSIZE(fft_freq));
                            if (fft_points > 1) {
                                int plot_count = fft_points - 1;
                                const double eps = 1e-12;
                                if (use_db_scale) {
                                    for (int i = 0; i < fft_points; ++i) {
                                        fft_mag_db[i] = 20.0 * std::log10(fft_mag[i] + eps);
                                    }
                                }
                                const double* mag_view = use_db_scale ? fft_mag_db : fft_mag;
                                const char* mag_label = use_db_scale ? "FFT (dB)" : "FFT";
                                const double* freq_view = fft_freq + 1;
                                const double* mag_plot = mag_view + 1;
                                const double* phase_plot = fft_phase + 1;
                                if (plot_count > 0 &&
                                    ImPlot::BeginPlot("FFT Magnitude", ImVec2(-1, 240))) {
                                    ImPlot::SetupAxis(ImAxis_X1, "Frequency (Hz)");
                                    ImPlot::SetupAxis(ImAxis_Y1,
                                                      use_db_scale ? "Magnitude (dB)"
                                                                   : "Magnitude");
                                    ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                                    ImPlot::PlotLine(mag_label, freq_view, mag_plot, plot_count);

                                    int peak_idx = 1;
                                    double peak_val = mag_view[1];
                                    for (int i = 2; i < fft_points; ++i) {
                                        if (mag_view[i] > peak_val) {
                                            peak_val = mag_view[i];
                                            peak_idx = i;
                                        }
                                    }
                                    double peak_freq = fft_freq[peak_idx];
                                    ImPlot::Annotation(peak_freq,
                                                       peak_val,
                                                       ImVec4(1, 1, 0, 1),
                                                       ImVec2(10, -10),
                                                       true,
                                                       "Peak: %.2f GHz",
                                                       peak_freq * 1e-9);

                                    double sim_freq_marker = sim->freq;
                                    ImPlot::DragLineX(1, &sim_freq_marker, ImVec4(1, 0, 0, 1));

                                    ImPlot::EndPlot();
                                }

                                if (plot_count > 0 &&
                                    ImPlot::BeginPlot("FFT Phase", ImVec2(-1, 180))) {
                                    ImPlot::SetupAxis(ImAxis_X1, "Frequency (Hz)");
                                    ImPlot::SetupAxis(ImAxis_Y1, "Phase (deg)");
                                    ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                                    ImPlot::PlotLine("Phase", freq_view, phase_plot, plot_count);
                                    ImPlot::EndPlot();
                                }

                                if (ImGui::Button("Export FFT to CSV")) {
                                    FILE* fp = fopen("fft_export.csv", "w");
                                    if (fp) {
                                        fprintf(fp, "# FFT Export\n");
                                        fprintf(fp, "# Center frequency: %.12e Hz\n", sim->freq);
                                        fprintf(fp,
                                                "Frequency(Hz),Magnitude,Magnitude(dB),Phase(deg)\n");
                                        for (int i = 0; i < fft_points; ++i) {
                                            double mag_db = 20.0 * std::log10(fft_mag[i] + eps);
                                            fprintf(fp,
                                                    "%.12e,%.12e,%.12e,%.12e\n",
                                                    fft_freq[i],
                                                    fft_mag[i],
                                                    mag_db,
                                                    fft_phase[i]);
                                        }
                                        fclose(fp);
                                        ui_log_add(&app, "FFT exported: fft_export.csv");
                                    } else {
                                        ui_log_add(&app, "Failed to export FFT");
                                    }
                                }
                            } else {
                                ImGui::TextUnformatted("FFT buffer too small");
                            }
                        } else {
                            ImGui::TextUnformatted("No FFT data available");
                        }
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
        }

        if (app.show_log_panel) {
        draw_log_panel(&app);
    }

    if (app.show_measurement_history) {
        draw_measurement_history_panel(&app);
    }

    if (app.show_print_composer) {
        draw_print_composer(&app, render->renderer, sim, &scope, render);
    }

        if (app.annotation_mode) {
            ImGui::OpenPopup("AddAnnotation");
            app.annotation_mode = false;
        }
        if (ImGui::BeginPopup("AddAnnotation")) {
            ImGui::InputText("Text", app.temp_annotation.text, IM_ARRAYSIZE(app.temp_annotation.text));
            ImGui::ColorEdit3("Color", (float*)&app.temp_annotation.color);
            ImGui::SliderFloat("Size", &app.temp_annotation.font_size, 8.0f, 28.0f, "%.0f px");
            if (ImGui::Button("Add", ImVec2(120, 0))) {
                app.measurements.annotations.push_back(app.temp_annotation);
                ImGui::CloseCurrentPopup();
                app.annotation_mode = false;
                ui_log_add(&app, "Annotation added");
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
                app.annotation_mode = false;
            }
            ImGui::EndPopup();
        }

        if (app.request_rebootstrap) {
            const char* restart_msg =
                (app.rebootstrap_message[0] != '\0') ? app.rebootstrap_message
                                                     : "Rebooted simulation";
            if (rebootstrap_simulation(&wizard.cfg, &bootstrap, &sim, &scope, scale)) {
                for (int i = 0; i < 4; ++i) {
                    int tmp_scale = (int)std::lround(kDefaultViewportZoom);
                    if (tmp_scale < 1) tmp_scale = 1;
                    reset_view_transform(&app.viewports[i], render, sim, &tmp_scale, kDefaultViewportZoom);
                }
                scale = (int)std::lround(kDefaultViewportZoom);
                if (scale < 1) scale = 1;
                width = sim->nx * scale;
                height = sim->ny * scale;
                if (width < min_window_width) width = min_window_width;
                if (height < min_window_height) height = min_window_height;
                SDL_SetWindowSize(render->window, width, height);
                app.window_width = width;
                app.window_height = height;
                ui_log_add(&app, "%s", restart_msg);
            } else {
                ui_log_add(&app, "Failed to rebootstrap simulation");
            }
            app.request_rebootstrap = false;
            app.rebootstrap_message[0] = '\0';
        }

        ImGui::SetNextWindowPos(ImVec2(main_viewport->Pos.x,
                                      main_viewport->Pos.y + main_viewport->Size.y - status_bar_h));
        ImGui::SetNextWindowSize(ImVec2(main_viewport->Size.x, status_bar_h));
        ImGuiWindowFlags status_flags = ImGuiWindowFlags_NoTitleBar |
                                        ImGuiWindowFlags_NoResize |
                                        ImGuiWindowFlags_NoMove |
                                        ImGuiWindowFlags_NoScrollbar |
                                        ImGuiWindowFlags_NoDocking |
                                        ImGuiWindowFlags_NoSavedSettings |
                                        ImGuiWindowFlags_NoCollapse;
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8.0f, 4.0f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.08f, 0.08f, 0.09f, 1.0f));
        if (ImGui::Begin("StatusBarOverlay", nullptr, status_flags)) {
            double sim_time_ns = sim ? (sim->timestep * sim->dt * 1e9) : 0.0;
            ImGui::Text("Step: %d | Time: %.2f ns | FPS: %.0f",
                        sim ? sim->timestep : 0,
                        sim_time_ns,
                        fps_avg);
            ImGui::SameLine();
            ImGui::SetCursorPosX(ImGui::GetContentRegionAvail().x + ImGui::GetStyle().WindowPadding.x - 260.0f);
            ViewportInstance* active_vp = ensure_active_viewport();
            float status_zoom = active_vp ? active_vp->zoom : app.viewport_zoom;
            float status_pan_x = active_vp ? active_vp->pan_x : app.viewport_pan_x;
            float status_pan_y = active_vp ? active_vp->pan_y : app.viewport_pan_y;
            ImGui::Text("Zoom: %.1fx | Pan: (%.0f, %.0f)",
                        status_zoom,
                        status_pan_x,
                        status_pan_y);
        }
        ImGui::End();
        ImGui::PopStyleColor();
        ImGui::PopStyleVar();

        // Help overlay (F1)
    if (show_help_overlay) {
        ImGui::SetNextWindowSize(ImVec2(600.0f, 400.0f), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
                                ImGuiCond_Appearing,
                                ImVec2(0.5f, 0.5f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.08f, 0.08f, 0.10f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.08f, 0.08f, 0.10f, 1.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 6.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 1.0f);
        bool help_open = ImGui::Begin("Keyboard Shortcuts & Help",
                                      &show_help_overlay,
                                      ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoDocking);
        if (help_open) {
            if (ImGui::BeginTabBar("HelpTabs")) {
                if (ImGui::BeginTabItem("Shortcuts")) {
                    ImGui::TextUnformatted("File Operations");
                        ImGui::Separator();
                        ImGui::Columns(2, "file_ops_cols", false);
                        ImGui::SetColumnWidth(0, 140.0f);
                        ImGui::TextUnformatted("Ctrl+S");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Save configuration");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F2");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Save screenshot");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F3");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Export scope FFT");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F5");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Load waveguide preset");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("F6");
                        ImGui::NextColumn();
                        ImGui::TextUnformatted("Load CPW filter preset");
                        ImGui::Columns(1);

                        ImGui::Separator();
                        ImGui::TextWrapped("Workflow:");
                ImGui::BulletText(
                    "Direct manipulation: Paint materials (M/U/I/O/P), drag sources with mouse");
                ImGui::BulletText(
                    "Presets: F5 (waveguide), F6 (CPW filter)");
                ImGui::BulletText(
                    "Save/Load: Ctrl+S saves to configs/saved_config.json");
                ImGui::BulletText(
                    "Materials menu \xE2\x86\x92 Material Browser to select materials");
                ImGui::BulletText(
                    "Numpad 1-8: Quick-select materials");
                ImGui::BulletText(
                    "I/O/P keys only in legacy mode (no material selected)");
                ImGui::BulletText("Viewports: Alt+1/2/3/4 for Single/H/V/Quad; hover sets active");
                ImGui::BulletText("Viz mode per viewport via toolbar; sync zoom/pan toggles available");
                ImGui::BulletText("Zoom: Mouse wheel, Ctrl+= / Ctrl+-, Reset: Ctrl+0 or Home");
                ImGui::BulletText("Pan: Middle drag or Shift+Right drag, Arrow keys (Shift = fast)");

                        ImGui::Separator();
                        ImGui::TextUnformatted("Press F1 again to close this window.");
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Materials & Visualization")) {
                        ImGui::TextUnformatted("Material Library:");
                        ImGui::BulletText("Materials menu \xE2\x86\x92 Material Browser");
                        ImGui::BulletText("Select from 11 curated EM materials");
                        ImGui::BulletText("4 Metals: Copper, Gold, Silver, Aluminum");
                        ImGui::BulletText("7 Dielectrics: Air, FR4, Rogers, Teflon, Silicon, Glass, Alumina");
                        ImGui::BulletText("Each material has accurate \xCE\xB5\xE1\xB5\xBD, tan(δ), σ values");
                        ImGui::BulletText("Numpad 1-8: Quick-select materials");
                        ImGui::BulletText("I/O/P: Legacy manual mode (when no material selected)");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("Material Visualization:");
                        ImGui::BulletText("E key: Cycle view modes (Field \xE2\x86\x92 Material \xE2\x86\x92 Overlay)");
                        ImGui::BulletText("View menu \xE2\x86\x92 toggles for overlay alpha, outlines, legend");
                        ImGui::BulletText("Material legend lists library materials with colors");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("View Modes:");
                        ImGui::BulletText("Field: Standard Ez field visualization");
                        ImGui::BulletText("Material: Show eps distribution with material colors");
                        ImGui::BulletText("Overlay: Field + semi-transparent material layer");
                        ImGui::BulletText("Per-viewport viz: toolbar selector (Field/Material/Overlay/Magnitude)");
                        ImGui::BulletText("Layouts: Alt+1/2/3/4 (Single/Horz/Vert/Quad), hover sets active pane");
                        ImGui::BulletText("Sync: optional sync zoom/pan toggles in toolbar");
                        ImGui::Spacing();
                        ImGui::TextUnformatted("S-Parameters:");
                        ImGui::BulletText("Analysis menu \xE2\x86\x92 S-Parameters window");
                        ImGui::BulletText("Plots show S21, S11, and VSWR vs frequency");
                        ImGui::BulletText("-3 dB bandwidth markers derived from peak S21");
                        ImGui::BulletText("Export enhanced CSV with VSWR and annotations");
                        ImGui::BulletText("Smith Chart window visualizes reflection coefficients");
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }
        }
        ImGui::End();
        ImGui::PopStyleVar(2);
        ImGui::PopStyleColor(2);
    }
        SDL_Texture* viewport_texture = NULL;
        static SDL_Texture* field_texture = NULL;
        static int field_tex_w = 0;
        static int field_tex_h = 0;

        if (frame_counter <= 120) {
            debug_logf("frame %d: before step paused=%d", frame_counter, paused ? 1 : 0);
        }

        if (!paused) {
            for (int s = 0; s < steps_per_frame; ++s) {
                fdtd_step(sim);
                int px = sim->nx / 2;
                int py = sim->ny / 2;
                double probe_val = fdtd_get_Ez(sim, px, py);
                scope_push(&scope, probe_val);
            }
        }

        if (frame_counter <= 120) {
            debug_logf("frame %d: after step scope_head=%d", frame_counter, scope.head);
        }

        Uint64 now = SDL_GetPerformanceCounter();
        double dt = (double)(now - prev) / (double)perf_freq;
        prev = now;
        double fps_inst = (dt > 0.0) ? (1.0 / dt) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;
        if (app.hud_zoom_timer > 0.0f) {
            app.hud_zoom_timer = std::max(0.0f, app.hud_zoom_timer - (float)dt);
        }
        if (app.hud_pan_timer > 0.0f) {
            app.hud_pan_timer = std::max(0.0f, app.hud_pan_timer - (float)dt);
        }

        // Update visualization metrics (field and scope scaling)
        {
            // Field amplitude (from last simulation step)
            double field_absmax = sim->step_Ez_absmax;
            if (field_absmax < 1e-12) field_absmax = 1e-12;

            bool field_autop = auto_rescale || hold_scope;
            if (!hold_color && field_autop) {
                // Smooth towards recent maximum
                const double alpha = 0.1;
                if (vmax_smooth <= 0.0) {
                    vmax_smooth = field_absmax;
                } else {
                    vmax_smooth = (1.0 - alpha) * vmax_smooth + alpha * field_absmax;
                }
            } else if (hold_color) {
                // Only allow increases when holding color scale
                if (field_absmax > vmax_smooth) {
                    vmax_smooth = field_absmax;
                }
            }

            // Scope amplitude (use absolute max over buffer)
            double scope_absmax = 0.0;
            if (scope.y && scope.n > 0) {
                bool first = true;
                for (int i = 0; i < scope.n; ++i) {
                    double v = std::fabs(scope.y[i]);
                    if (first || v > scope_absmax) {
                        scope_absmax = v;
                        first = false;
                    }
                }
            }

            bool scope_autop = auto_rescale || hold_color;
            if (!hold_scope && scope_autop && scope_absmax > 0.0) {
                const double alpha_scope = 0.1;
                if (scope_vmax <= 0.0) {
                    scope_vmax = scope_absmax;
                } else {
                    scope_vmax =
                        (1.0 - alpha_scope) * scope_vmax + alpha_scope * scope_absmax;
                }
            } else if (hold_scope && scope_absmax > scope_vmax) {
                scope_vmax = scope_absmax;
            }
        }

        // Render field into a texture that we display in ImGui
        if (app.viewport_valid) {
            int tex_w = (int)app.viewport_size.x;
            int tex_h = (int)app.viewport_size.y;
            if (tex_w < 1) tex_w = 1;
            if (tex_h < 1) tex_h = 1;
            if (!field_texture || tex_w != field_tex_w || tex_h != field_tex_h) {
                if (field_texture) {
                    SDL_DestroyTexture(field_texture);
                }
                field_texture =
                    SDL_CreateTexture(render->renderer,
                                      SDL_PIXELFORMAT_RGBA8888,
                                      SDL_TEXTUREACCESS_TARGET,
                                      tex_w,
                                      tex_h);
                if (frame_counter <= 120) {
                    debug_logf("frame %d: recreate field_texture %dx%d %s",
                               frame_counter,
                               tex_w,
                               tex_h,
                               field_texture ? "ok" : SDL_GetError());
                }
                field_tex_w = tex_w;
                field_tex_h = tex_h;
            }
            viewport_texture = field_texture;
            SDL_SetRenderTarget(render->renderer, viewport_texture);
            SDL_SetRenderDrawColor(render->renderer,
                                   viewport_clear_color.r,
                                   viewport_clear_color.g,
                                   viewport_clear_color.b,
                                   viewport_clear_color.a);
            SDL_RenderClear(render->renderer);

            bool channel_na_flags[4] = {false, false, false, false};
            for (int vp_idx = 0; vp_idx < 4; ++vp_idx) {
                ViewportInstance& vp = app.viewports[vp_idx];
                if (!vp.valid) continue;
                SDL_Rect vp_rect = {(int)vp.pos.x, (int)vp.pos.y, (int)vp.size.x, (int)vp.size.y};
                SDL_RenderSetViewport(render->renderer, &vp_rect);

                int vp_scale = (int)std::lround(vp.zoom);
                if (vp_scale < 1) vp_scale = 1;
                render->scale = vp_scale;
                ImVec2 vp_offset = compute_viewport_offset(vp, sim, vp_scale);
                render->offset_x = vp_offset.x;
                render->offset_y = vp_offset.y;

                double vmax = vmax_smooth;
                if (vmax <= 0.0) vmax = 1.0;

                bool channel_na = false;
                switch (vp.viz_mode) {
                    case VIEWPORT_VIZ_EZ:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_EZ_ABS:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EZ_ABS, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_HX:
                        render_field_channel_heatmap(render, sim, FIELD_CH_HX, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_HY:
                        render_field_channel_heatmap(render, sim, FIELD_CH_HY, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_HMAG:
                        render_field_channel_heatmap(render, sim, FIELD_CH_H_MAG, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_SX:
                        render_field_channel_heatmap(render, sim, FIELD_CH_SX, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_SY:
                        render_field_channel_heatmap(render, sim, FIELD_CH_SY, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_S_MAG:
                        render_field_channel_heatmap(render, sim, FIELD_CH_S_MAG, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_EX:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EX, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_EY:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EY, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_HZ:
                        render_field_channel_heatmap(render, sim, FIELD_CH_HZ, vmax, 1.0, &channel_na);
                        break;
                    case VIEWPORT_VIZ_MATERIAL:
                        render_material_distribution(render, sim, vp_scale);
                        break;
                    case VIEWPORT_VIZ_OVERLAY:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
                        SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_BLEND);
                        render_material_overlay(render,
                                                sim,
                                                vp_scale,
                                                app.material_overlay_alpha);
                        SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_NONE);
                        break;
                    default:
                        render_field_channel_heatmap(render, sim, FIELD_CH_EZ, vmax, 1.0, &channel_na);
                        break;
                }
                channel_na_flags[vp_idx] = channel_na;

                if (app.show_material_outlines) {
                    SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_BLEND);
                    render_material_outlines(render, sim, vp_scale);
                    SDL_SetRenderDrawBlendMode(render->renderer, SDL_BLENDMODE_NONE);
                }

                if (vp.show_grid) {
                    SDL_Color grid_col = {40, 40, 50, 255};
                    int grid_step = compute_grid_step_from_scale(vp_scale);
                    render_grid_overlay(render, sim, grid_col, grid_step);
                }
                if (vp.show_sources) {
                    render_sources(render, sim->sources);
                }
                if (vp.show_vectors) {
                    bool use_h = (vp.viz_mode == VIEWPORT_VIZ_HX ||
                                  vp.viz_mode == VIEWPORT_VIZ_HY ||
                                  vp.viz_mode == VIEWPORT_VIZ_HMAG);
                    bool use_s = (vp.viz_mode == VIEWPORT_VIZ_SX ||
                                  vp.viz_mode == VIEWPORT_VIZ_SY ||
                                  vp.viz_mode == VIEWPORT_VIZ_S_MAG);
                    if (use_h || use_s) {
                        SDL_SetRenderDrawColor(render->renderer, 80, 220, 255, 180);
                        int stride = (int)std::max(4.0f, 24.0f / (float)vp_scale);
                        float arrow_scale = std::min(18.0f, (float)vp_scale * 3.0f);
                        for (int i = 0; i < sim->nx; i += stride) {
                            for (int j = 0; j < sim->ny; j += stride) {
                                float vx = 0.0f, vy = 0.0f;
                                if (use_h) {
                                    vx = (float)sim->Hy[i][j];
                                    vy = -(float)sim->Hx[i][j];
                                } else if (use_s) {
                                    vx = (float)(sim->Ez[i][j] * sim->Hy[i][j]);
                                    vy = (float)(-sim->Ez[i][j] * sim->Hx[i][j]);
                                }
                                float mag = std::sqrt(vx * vx + vy * vy);
                                if (mag < 1e-6f) continue;
                                vx /= mag;
                                vy /= mag;
                                float len = arrow_scale;
                                float cx = vp_offset.x + ((float)i + 0.5f) * (float)vp_scale;
                                float cy = vp_offset.y + ((float)j + 0.5f) * (float)vp_scale;
                                SDL_RenderDrawLine(render->renderer,
                                                   (int)cx,
                                                   (int)cy,
                                                   (int)(cx + vx * len),
                                                   (int)(cy + vy * len));
                            }
                        }
                    }
                }
            }

            render->offset_x = 0.0f;
            render->offset_y = 0.0f;
            SDL_RenderSetViewport(render->renderer, NULL);
            SDL_SetRenderTarget(render->renderer, NULL);
        }

        // Clear backbuffer for ImGui
        if (frame_counter <= 120) {
            debug_logf("frame %d: after field render valid=%d", frame_counter, app.viewport_valid ? 1 : 0);
        }
        SDL_SetRenderDrawColor(render->renderer,
                               viewport_clear_color.r,
                               viewport_clear_color.g,
                               viewport_clear_color.b,
                               viewport_clear_color.a);
        SDL_RenderClear(render->renderer);

        // Viewport content (image + overlay toolbar)
        bool channel_na_flags[4] = {false, false, false, false};
        app.toolbar_valid = false;
        if (app.viewport_valid && viewport_texture) {
            if (frame_counter <= 120) {
                debug_logf("frame %d: entering viewport overlay window", frame_counter);
            }
            ImGui::SetNextWindowPos(app.viewport_pos);
            ImGui::SetNextWindowSize(app.viewport_size);
            ImGui::SetNextWindowViewport(ImGui::GetMainViewport()->ID);
            ImGuiWindowFlags overlay_flags = ImGuiWindowFlags_NoTitleBar |
                                             ImGuiWindowFlags_NoResize |
                                             ImGuiWindowFlags_NoMove |
                                             ImGuiWindowFlags_NoScrollbar |
                                             ImGuiWindowFlags_NoSavedSettings |
                                             ImGuiWindowFlags_NoDocking |
                                             ImGuiWindowFlags_NoCollapse |
                                             ImGuiWindowFlags_NoBackground;
            if (ImGui::Begin("ViewportOverlayWindow", nullptr, overlay_flags)) {
                ImGui::SetCursorPos(ImVec2(0, 0));
                ImGui::Image((ImTextureID)viewport_texture, app.viewport_size);

                ViewportInstance* active_vp = ensure_active_viewport();
                int active_scale = active_vp ? (int)std::lround(active_vp->zoom) : 1;
                if (active_scale < 1) active_scale = 1;
                if (frame_counter <= 120) {
                    debug_logf("frame %d: overlay start active_idx=%d active_scale=%d",
                               frame_counter,
                               app.active_viewport_idx,
                               active_scale);
                }

                ImDrawList* overlay_dl = ImGui::GetWindowDrawList();
                ImVec2 win_origin = ImGui::GetWindowPos();
                const char* labels[] = {"A", "B", "C", "D"};
                const char* viz_labels[] = {"Ez", "|Ez|", "Hx", "Hy", "|H|", "Sx", "Sy", "|S|", "Ex", "Ey", "Hz", "Material", "Overlay"};
                int hovered_idx = -1;
                const ImGuiIO& io = ImGui::GetIO();
                ImVec2 mp = ImGui::GetIO().MousePos;
                bool over_toolbar = false;
                if (app.toolbar_valid) {
                    float tx0 = app.toolbar_screen_min.x;
                    float ty0 = app.toolbar_screen_min.y;
                    float tx1 = tx0 + app.toolbar_screen_size.x;
                    float ty1 = ty0 + app.toolbar_screen_size.y;
                    over_toolbar = (mp.x >= tx0 && mp.x <= tx1 && mp.y >= ty0 && mp.y <= ty1);
                }
                for (int vp_idx = 0; vp_idx < 4; ++vp_idx) {
                    const ViewportInstance& vp = app.viewports[vp_idx];
                    if (!vp.valid) continue;
                    ImVec2 p0 = ImVec2(win_origin.x + vp.pos.x, win_origin.y + vp.pos.y);
                    ImVec2 p1 = ImVec2(p0.x + vp.size.x, p0.y + vp.size.y);
                    if (mp.x >= p0.x && mp.x < p1.x && mp.y >= p0.y && mp.y < p1.y) {
                        hovered_idx = vp_idx;
                    }
                    ImU32 col = (vp_idx == app.active_viewport_idx)
                                    ? IM_COL32(80, 160, 255, 255)
                                    : IM_COL32(80, 80, 80, 180);
                    overlay_dl->AddRect(p0, p1, col, 0.0f, 0, 2.0f);

                    char label_buf[64];
                    const char* viz = viz_labels[(int)vp.viz_mode];
                    std::snprintf(label_buf, sizeof(label_buf), "%s: %s", labels[vp_idx], viz);
                    ImVec2 pad(6.0f, 3.0f);
                    ImVec2 text_sz = ImGui::CalcTextSize(label_buf);
                    ImVec2 l0 = ImVec2(p0.x + 8.0f, p0.y + 8.0f);
                    ImVec2 l1 = ImVec2(l0.x + text_sz.x + pad.x * 2.0f, l0.y + text_sz.y + pad.y * 2.0f);
                    overlay_dl->AddRectFilled(l0, l1, IM_COL32(0, 0, 0, 180), 3.0f);
                    overlay_dl->AddText(ImVec2(l0.x + pad.x, l0.y + pad.y),
                                        IM_COL32(230, 230, 230, 255),
                                        label_buf);
                    if (channel_na_flags[vp_idx]) {
                        ImVec2 warn_pos = ImVec2(p0.x + 10.0f, p0.y + 30.0f);
                        ImVec2 b0 = ImVec2(warn_pos.x - 4.0f, warn_pos.y - 4.0f);
                        ImVec2 b1 = ImVec2(warn_pos.x + 200.0f, warn_pos.y + 22.0f);
                        overlay_dl->AddRectFilled(b0, b1, IM_COL32(30, 18, 18, 200), 4.0f);
                        overlay_dl->AddText(warn_pos, IM_COL32(255, 120, 120, 255), "Channel unavailable (2D model)");
                    }

                    // Constrain overlays to the viewport rect to avoid UI bleed
                    overlay_dl->PushClipRect(p0, p1, true);
                    int vp_scale = (int)std::lround(vp.zoom);
                    if (vp_scale < 1) vp_scale = 1;
                    ImVec2 viewport_offset = compute_viewport_offset(vp, sim, vp_scale);

                    if (app.composer_pick_region_active && ImGui::IsKeyPressed(ImGuiKey_Escape)) {
                        app.composer_pick_region_active = false;
                        app.composer_pick_region_dragging = false;
                        std::snprintf(app.composer_status, sizeof(app.composer_status), "Region pick canceled");
                    }

                    if (!show_help_overlay) {
                        // Region pick overlay (composer)
                        bool picking_this_vp = app.composer_pick_region_active &&
                                               (app.composer_pick_region_viewport == vp_idx);
                        if (picking_this_vp) {
                            float local_x = io.MousePos.x - p0.x;
                            float local_y = io.MousePos.y - p0.y;
                            bool inside = (local_x >= 0.0f && local_y >= 0.0f &&
                                           local_x < vp.size.x && local_y < vp.size.y);
                            auto to_norm = [&](float lx, float ly) {
                                float fx = lx - viewport_offset.x;
                                float fy = ly - viewport_offset.y;
                                int ix = (int)std::floor(fx / (float)vp_scale);
                                int iy = (int)std::floor(fy / (float)vp_scale);
                                ix = std::clamp(ix, 0, sim->nx - 1);
                                iy = std::clamp(iy, 0, sim->ny - 1);
                                return ImVec2((float)ix / (float)sim->nx,
                                              (float)iy / (float)sim->ny);
                            };

                            if (!app.composer_pick_region_dragging && inside && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                                app.composer_pick_region_start_norm = to_norm(local_x, local_y);
                                app.composer_pick_region_end_norm = app.composer_pick_region_start_norm;
                                app.composer_pick_region_dragging = true;
                            } else if (app.composer_pick_region_dragging && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                                app.composer_pick_region_end_norm = to_norm(local_x, local_y);
                            } else if (app.composer_pick_region_dragging && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                                app.composer_pick_region_dragging = false;
                                ImVec2 a = app.composer_pick_region_start_norm;
                                ImVec2 b = app.composer_pick_region_end_norm;
                                float x0 = std::clamp(std::min(a.x, b.x), 0.0f, 1.0f);
                                float y0 = std::clamp(std::min(a.y, b.y), 0.0f, 1.0f);
                                float x1 = std::clamp(std::max(a.x, b.x), 0.0f, 1.0f);
                                float y1 = std::clamp(std::max(a.y, b.y), 0.0f, 1.0f);
                                const float min_delta = 1.0f / std::max(4.0f, (float)std::min(sim->nx, sim->ny));
                                if (x1 - x0 < min_delta) x1 = x0 + min_delta;
                                if (y1 - y0 < min_delta) y1 = y0 + min_delta;
                                if (x1 > 1.0f) x1 = 1.0f;
                                if (y1 > 1.0f) y1 = 1.0f;

                                bool updated = false;
                                if (app.composer_pick_region_page >= 0 &&
                                    app.composer_pick_region_page < (int)app.composer_pages.size()) {
                                    ComposerPage& pg = app.composer_pages[app.composer_pick_region_page];
                                    if (app.composer_pick_region_target_item >= 0) {
                                        for (auto& it : pg.items) {
                                            if (it.id == app.composer_pick_region_target_item &&
                                                it.type == COMPOSER_REGION) {
                                                it.region_norm = ImVec4(x0, y0, x1, y1);
                                                it.viewport_idx = app.composer_pick_region_viewport;
                                                updated = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!updated) {
                                        float region_w = x1 - x0;
                                        float region_h = y1 - y0;
                                        float aspect = region_h > 1e-6f ? (region_w / region_h) : 1.0f;
                                        float target_w = pg.res_w * 0.45f;
                                        float target_h = target_w / aspect;
                                        if (target_h > pg.res_h * 0.8f) {
                                            target_h = pg.res_h * 0.8f;
                                            target_w = target_h * aspect;
                                        }
                                        ImVec2 pos(pg.res_w * 0.05f, pg.res_h * 0.1f);
                                        composer_add_item(&app,
                                                          pg,
                                                          COMPOSER_REGION,
                                                          app.composer_pick_region_viewport,
                                                          pos,
                                                          ImVec2(target_w, target_h));
                                        pg.items.back().region_norm = ImVec4(x0, y0, x1, y1);
                                        updated = true;
                                    }
                                }
                                app.composer_pick_region_active = false;
                                std::snprintf(app.composer_status,
                                              sizeof(app.composer_status),
                                              updated ? "Region captured (viewport %d)" : "Region capture failed",
                                              app.composer_pick_region_viewport + 1);
                            }

                            if (app.composer_pick_region_dragging || picking_this_vp) {
                                ImVec2 a = app.composer_pick_region_dragging ? app.composer_pick_region_start_norm
                                                                             : app.composer_pick_region_start_norm;
                                ImVec2 b = app.composer_pick_region_dragging ? app.composer_pick_region_end_norm
                                                                             : app.composer_pick_region_end_norm;
                                float x0_px = p0.x + viewport_offset.x + a.x * (float)sim->nx * (float)vp_scale;
                                float y0_px = p0.y + viewport_offset.y + a.y * (float)sim->ny * (float)vp_scale;
                                float x1_px = p0.x + viewport_offset.x + b.x * (float)sim->nx * (float)vp_scale;
                                float y1_px = p0.y + viewport_offset.y + b.y * (float)sim->ny * (float)vp_scale;
                                ImVec2 rp0(std::min(x0_px, x1_px), std::min(y0_px, y1_px));
                                ImVec2 rp1(std::max(x0_px, x1_px), std::max(y0_px, y1_px));
                                overlay_dl->AddRect(rp0, rp1, IM_COL32(0, 200, 255, 255), 0.0f, 0, 2.0f);
                                overlay_dl->AddRectFilled(rp0, rp1, IM_COL32(0, 200, 255, 50));
                            }
                        }

                        const ImU32 src_col = IM_COL32(255, 255, 255, 230);
                        for (int k = 0; k < MAX_SRC; ++k) {
                            const Source& s = sim->sources[k];
                            if (!s.active) continue;
                            float sx = p0.x + viewport_offset.x + ((float)s.ix + 0.5f) * (float)vp_scale;
                            float sy = p0.y + viewport_offset.y + ((float)s.iy + 0.5f) * (float)vp_scale;
                            char label[16];
                            std::snprintf(label, sizeof(label), "%d", k);
                            overlay_dl->AddText(ImVec2(sx + 6.0f, sy - 6.0f), src_col, label);
                            if (dragging_source && k == dragged_source_idx) {
                                float radius = 10.0f;
                                overlay_dl->AddCircle(ImVec2(sx, sy),
                                                      radius,
                                                      IM_COL32(255, 200, 0, 200),
                                                      16,
                                                      3.0f);
                            }
                        }

                        const ImU32 block_col = IM_COL32(255, 255, 0, 160);
                        const ImU32 block_sel_col = IM_COL32(0, 255, 128, 200);
                        int rect_count = wizard.cfg.material_rect_count;
                        if (rect_count < 0) rect_count = 0;
                        if (rect_count > CONFIG_MAX_MATERIAL_RECTS) rect_count = CONFIG_MAX_MATERIAL_RECTS;
                        const Material* overlay_sel_mat = nullptr;
                        if (app.selected_material_id >= 0) {
                            overlay_sel_mat = material_library_get_by_id(app.selected_material_id);
                        }
                        for (int i = 0; i < rect_count; ++i) {
                            const MaterialRectSpec& r = wizard.cfg.material_rects[i];
                            float x0 = p0.x + viewport_offset.x +
                                       (float)(r.x0 * (double)sim->nx * (double)vp_scale);
                            float x1 = p0.x + viewport_offset.x +
                                       (float)(r.x1 * (double)sim->nx * (double)vp_scale);
                            float y0 = p0.y + viewport_offset.y +
                                       (float)(r.y0 * (double)sim->ny * (double)vp_scale);
                            float y1 = p0.y + viewport_offset.y +
                                       (float)(r.y1 * (double)sim->ny * (double)vp_scale);
                            ImVec2 b0(x0, y0);
                            ImVec2 b1(x1, y1);
                            bool mat_match = overlay_sel_mat && rect_matches_material(r, overlay_sel_mat);
                            ImU32 bcol = block_col;
                            if (i == app.selected_block) {
                                bcol = block_sel_col;
                            } else if (mat_match) {
                                bcol = IM_COL32(0, 200, 255, 200);
                            }
                            overlay_dl->AddRect(b0, b1, bcol, 0.0f, 0, 1.5f);
                            char blabel[16];
                            std::snprintf(blabel, sizeof(blabel), "B%d", i);
                            overlay_dl->AddText(ImVec2(x0 + 4.0f, y0 + 4.0f), bcol, blabel);
                        }

                        if (paint_mode && app.active_viewport_idx == vp_idx) {
                            float local_x = io.MousePos.x - p0.x;
                            float local_y = io.MousePos.y - p0.y;
                            if (local_x >= 0.0f && local_y >= 0.0f &&
                                local_x < vp.size.x && local_y < vp.size.y) {
                                float field_x = local_x - viewport_offset.x;
                                float field_y = local_y - viewport_offset.y;
                                int ix = (int)(field_x / (float)vp_scale);
                                int iy = (int)(field_y / (float)vp_scale);
                                if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
                                    float cx = p0.x + viewport_offset.x +
                                               ((float)ix + 0.5f) * (float)vp_scale;
                                    float cy = p0.y + viewport_offset.y +
                                               ((float)iy + 0.5f) * (float)vp_scale;
                                    float radius = (float)paint_brush_size * (float)vp_scale;
                                    ImU32 col_p = IM_COL32(255, 255, 0, 180);
                                    if (paint_material_id >= 0) {
                                        const Material* mat = material_library_get_by_id(paint_material_id);
                                        if (mat) {
                                            col_p = IM_COL32(mat->color_r, mat->color_g, mat->color_b, 200);
                                        }
                                    } else {
                                        if (paint_material_type == 0) {
                                            col_p = IM_COL32(255, 80, 80, 200);
                                        } else if (paint_material_type == 1) {
                                            col_p = IM_COL32(80, 160, 255, 200);
                                        } else {
                                            col_p = IM_COL32(120, 255, 120, 200);
                                        }
                                    }
                                    overlay_dl->AddCircle(ImVec2(cx, cy), radius, col_p, 32, 1.5f);
                                }
                            }
                        }
                    }
                    overlay_dl->PopClipRect();
                }

                if (!over_toolbar &&
                    hovered_idx >= 0 &&
                    ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows)) {
                    app.active_viewport_idx = hovered_idx;
                    active_vp = ensure_active_viewport();
                    active_scale = active_vp ? (int)std::lround(active_vp->zoom) : active_scale;
                    if (active_scale < 1) active_scale = 1;
                }

                if (frame_counter <= 120) {
                    debug_logf("frame %d: overlay after labels", frame_counter);
                }

                const float toolbar_w = 260.0f;
                const float toolbar_h = 90.0f;
                bool viewport_hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows |
                                                               ImGuiHoveredFlags_AllowWhenBlockedByPopup |
                                                               ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
                float toolbar_alpha = viewport_hovered ? 1.0f : 0.35f;
                ImVec2 toolbar_pos;
                if (active_vp) {
                    toolbar_pos = ImVec2(active_vp->pos.x + active_vp->size.x - toolbar_w - 12.0f,
                                         active_vp->pos.y + 12.0f);
                } else {
                    toolbar_pos = ImVec2(app.viewport_size.x - toolbar_w - 12.0f,
                                         12.0f);
                }
                if (toolbar_pos.x < 0.0f) toolbar_pos.x = 0.0f;
                if (toolbar_pos.y < 0.0f) toolbar_pos.y = 0.0f;
                app.toolbar_screen_min = ImVec2(app.viewport_pos.x + toolbar_pos.x,
                                                app.viewport_pos.y + toolbar_pos.y);
                app.toolbar_screen_size = ImVec2(toolbar_w, toolbar_h);
                app.toolbar_valid = true;
                ImGui::SetCursorPos(toolbar_pos);
                ImVec4 toolbar_bg(0.08f, 0.08f, 0.10f, 0.9f * toolbar_alpha);
                ImGui::PushStyleColor(ImGuiCol_ChildBg, toolbar_bg);
                ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(10.0f, 8.0f));
                        ImGui::PushStyleVar(ImGuiStyleVar_Alpha, 0.75f + 0.25f * toolbar_alpha);
                            if (ImGui::BeginChild("ViewportToolbar",
                                                  ImVec2(toolbar_w, toolbar_h),
                                                  true,
                                                  ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings)) {
                                if (frame_counter <= 120) {
                                    debug_logf("frame %d: overlay toolbar begin", frame_counter);
                                }
                                if (ImGui::Button("Cfg", ImVec2(36.0f, 0.0f))) {
                                    ImGui::OpenPopup("ViewportCfgPopup");
                                }
                                if (ImGui::BeginPopup("ViewportCfgPopup")) {
                                    if (active_vp) {
                                        const char* viz_options[] = {"Ez", "|Ez|", "Hx", "Hy", "|H|", "Sx", "Sy", "|S|", "Ex (n/a)", "Ey (n/a)", "Hz (n/a)", "Material", "Overlay"};
                                        int v = (int)active_vp->viz_mode;
                                        if (ImGui::Combo("Channel", &v, viz_options, IM_ARRAYSIZE(viz_options))) {
                                            active_vp->viz_mode = (ViewportViz)v;
                                        }
                                        ImGui::Separator();
                                        ImGui::Checkbox("Grid", &active_vp->show_grid);
                                        ImGui::Checkbox("Sources", &active_vp->show_sources);
                                        ImGui::Checkbox("Vectors", &active_vp->show_vectors);
                                    }
                                    ImGui::Separator();
                                    ImGui::Checkbox("Sync Zoom", &app.sync_zoom);
                                    ImGui::Checkbox("Sync Pan", &app.sync_pan);
                                    ImGui::EndPopup();
                            }
                            ImVec2 focus = active_vp ? ImVec2(active_vp->size.x * 0.5f, active_vp->size.y * 0.5f)
                                                     : ImVec2(app.viewport_size.x * 0.5f, app.viewport_size.y * 0.5f);
                            ImGui::TextUnformatted("View");
                    ImGui::PushButtonRepeat(true);
                    if (ImGui::Button("-", ImVec2(26.0f, 0.0f))) {
                        if (active_vp) {
                            apply_zoom_at_point(active_vp, render, sim, &active_scale, focus.x, focus.y, 0.9f, false);
                            if (app.sync_zoom) {
                                for (int i = 0; i < 4; ++i) {
                                    if (!app.viewports[i].valid) continue;
                                    app.viewports[i].zoom = active_vp->zoom;
                                }
                            }
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_zoom_timer = 1.0f;
                        }
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Fit", ImVec2(50.0f, 0.0f))) {
                        if (app.viewport_valid && sim && active_vp) {
                            float fit_scale_x = active_vp->size.x / (float)sim->nx;
                            float fit_scale_y = active_vp->size.y / (float)sim->ny;
                            float fit_zoom = std::max(0.5f, std::min(fit_scale_x, fit_scale_y));
                            active_vp->zoom = fit_zoom;
                            active_scale = (int)std::lround(fit_zoom);
                            if (active_scale < 1) active_scale = 1;
                            render->scale = active_scale;
                            active_vp->pan_x = 0.0f;
                            active_vp->pan_y = 0.0f;
                            clamp_pan_offset(active_vp, sim, active_scale);
                            if (app.sync_zoom) {
                                for (int i = 0; i < 4; ++i) {
                                    if (!app.viewports[i].valid) continue;
                                    app.viewports[i].zoom = active_vp->zoom;
                                    app.viewports[i].pan_x = active_vp->pan_x;
                                    app.viewports[i].pan_y = active_vp->pan_y;
                                }
                            }
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                            app.hud_zoom_timer = 1.0f;
                            app.hud_pan_timer = 1.0f;
                        }
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("1:1", ImVec2(40.0f, 0.0f))) {
                        if (app.viewport_valid && sim && active_vp) {
                            active_vp->zoom = 1.0f;
                            active_scale = 1;
                            render->scale = active_scale;
                            active_vp->pan_x = 0.0f;
                            active_vp->pan_y = 0.0f;
                            clamp_pan_offset(active_vp, sim, active_scale);
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                            app.hud_zoom_timer = 1.0f;
                            app.hud_pan_timer = 1.0f;
                        }
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Reset", ImVec2(60.0f, 0.0f))) {
                        if (active_vp) {
                            reset_view_transform(active_vp, render, sim, &active_scale, kDefaultViewportZoom);
                            clamp_pan_offset(active_vp, sim, active_scale);
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                            app.hud_zoom_timer = 1.0f;
                            app.hud_pan_timer = 1.0f;
                        }
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("+", ImVec2(26.0f, 0.0f))) {
                        if (active_vp) {
                            apply_zoom_at_point(active_vp, render, sim, &active_scale, focus.x, focus.y, 1.1f, false);
                            if (app.sync_zoom) {
                                for (int i = 0; i < 4; ++i) {
                                    if (!app.viewports[i].valid) continue;
                                    app.viewports[i].zoom = active_vp->zoom;
                                }
                            }
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_zoom_timer = 1.0f;
                        }
                    }
                    ImGui::PopButtonRepeat();
                    if (frame_counter <= 120) {
                        debug_logf("frame %d: toolbar after zoom buttons", frame_counter);
                    }

                    ImGui::Spacing();
                    ImGui::Separator();
                    ImGui::Spacing();

                    auto toggle_button = [&](const char* label, bool* value, const ImVec4& on_col) {
                        bool is_on = value ? *value : false;
                        if (!value) {
                            ImGui::BeginDisabled();
                        }
                        ImVec4 off_col(0.16f, 0.16f, 0.18f, 0.7f);
                        ImVec4 off_hover(0.20f, 0.20f, 0.24f, 0.8f);
                        ImVec4 off_active(0.24f, 0.24f, 0.28f, 0.9f);
                        auto clamp_color = [](const ImVec4& c, float delta) {
                            return ImVec4(std::clamp(c.x + delta, 0.0f, 1.0f),
                                          std::clamp(c.y + delta, 0.0f, 1.0f),
                                          std::clamp(c.z + delta, 0.0f, 1.0f),
                                          c.w);
                        };
                        ImVec4 on_hover = clamp_color(on_col, 0.05f);
                        ImVec4 on_active = clamp_color(on_col, 0.08f);
                        ImGui::PushStyleColor(ImGuiCol_Button, is_on ? on_col : off_col);
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, is_on ? on_hover : off_hover);
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, is_on ? on_active : off_active);
                        bool clicked = ImGui::Button(label, ImVec2(70.0f, 0.0f));
                        ImGui::PopStyleColor(3);
                        if (clicked && value) {
                            *value = !*value;
                            is_on = *value;
                        }
                        if (!value) {
                            ImGui::EndDisabled();
                        }
                    };

                    toggle_button("Grid",
                                  active_vp ? &active_vp->show_grid : nullptr,
                                  ImVec4(0.16f, 0.35f, 0.22f, 0.9f));
                    ImGui::SameLine();
                    toggle_button("Axis", &app.show_axis_overlay, ImVec4(0.20f, 0.32f, 0.50f, 0.9f));
                    ImGui::SameLine();
                    if (ImGui::Button("Center", ImVec2(70.0f, 0.0f))) {
                        if (app.viewport_valid && sim && active_vp) {
                            active_vp->pan_x = 0.0f;
                            active_vp->pan_y = 0.0f;
                            clamp_pan_offset(active_vp, sim, active_scale);
                            if (app.sync_pan) {
                                for (int i = 0; i < 4; ++i) {
                                    if (!app.viewports[i].valid) continue;
                                    app.viewports[i].pan_x = active_vp->pan_x;
                                    app.viewports[i].pan_y = active_vp->pan_y;
                                }
                            }
                            app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                            app.hud_pan_timer = 1.0f;
                        }
                    }
                    ImGui::SameLine();
                    ImGui::TextDisabled("|");
                    ImGui::SameLine();

                    ImVec2 offset = active_vp ? compute_viewport_offset(*active_vp, sim, active_scale)
                                              : ImVec2(0.0f, 0.0f);
                    int mx = 0;
                    int my = 0;
                    SDL_GetMouseState(&mx, &my);
                    float local_x = active_vp ? (float)mx - (app.viewport_pos.x + active_vp->pos.x) : -1.0f;
                    float local_y = active_vp ? (float)my - (app.viewport_pos.y + active_vp->pos.y) : -1.0f;
                    if (active_vp &&
                        local_x >= 0.0f && local_x < active_vp->size.x &&
                        local_y >= 0.0f && local_y < active_vp->size.y) {
                        float field_x = local_x - offset.x;
                        float field_y = local_y - offset.y;
                        int ix = (int)(field_x / (float)active_scale);
                        int iy = (int)(field_y / (float)active_scale);
                        if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
                            ImGui::Text("Cell: (%d, %d)", ix, iy);
                        } else {
                            ImGui::TextUnformatted("Cell: --");
                        }
                    } else {
                        ImGui::TextUnformatted("Cell: --");
                    }

                    ImGui::Separator();
                    if (frame_counter <= 120) {
                        debug_logf("frame %d: toolbar after cell block", frame_counter);
                    }
                    if (active_vp) {
                        const char* viz_options[] = {"Ez", "|Ez|", "Hx", "Hy", "|H|", "Sx", "Sy", "|S|", "Ex (n/a)", "Ey (n/a)", "Hz (n/a)", "Material", "Overlay"};
                        int v = (int)active_vp->viz_mode;
                        if (ImGui::Combo("Visualization", &v, viz_options, IM_ARRAYSIZE(viz_options))) {
                            active_vp->viz_mode = (ViewportViz)v;
                        }
                        ImGui::Checkbox("Sync Zoom", &app.sync_zoom);
                        ImGui::SameLine();
                        ImGui::Checkbox("Sync Pan", &app.sync_pan);
                    }
                    if (frame_counter <= 120) {
                        debug_logf("frame %d: toolbar end block", frame_counter);
                    }
                }
                ImGui::EndChild();
                ImGui::PopStyleVar(2);
                ImGui::PopStyleColor();

                // HUD badges
                ImDrawList* dl = ImGui::GetWindowDrawList();
                ImVec2 hud_pos = active_vp ? ImVec2(app.viewport_pos.x + active_vp->pos.x + 12.0f,
                                                 app.viewport_pos.y + active_vp->pos.y + 12.0f)
                                           : ImVec2(app.viewport_pos.x + 12.0f, app.viewport_pos.y + 12.0f);
                if (app.hud_zoom_timer > 0.0f) {
                    char buf[64];
                    std::snprintf(buf, sizeof(buf), "Zoom %.2fx", app.hud_zoom_value);
                    ImU32 bg = IM_COL32(20, 20, 30, (int)(200 * app.hud_zoom_timer));
                    ImU32 fg = IM_COL32(230, 230, 240, (int)(255 * app.hud_zoom_timer));
                    ImVec2 text_sz = ImGui::CalcTextSize(buf);
                    ImVec2 pad(8.0f, 4.0f);
                    ImVec2 box_min = hud_pos;
                    ImVec2 box_max = ImVec2(hud_pos.x + text_sz.x + pad.x * 2.0f,
                                            hud_pos.y + text_sz.y + pad.y * 2.0f);
                    dl->AddRectFilled(box_min, box_max, bg, 4.0f);
                    dl->AddText(ImVec2(box_min.x + pad.x, box_min.y + pad.y), fg, buf);
                    hud_pos.y += text_sz.y + pad.y * 2.0f + 6.0f;
                }
                if (app.hud_pan_timer > 0.0f) {
                    char buf[64];
                    std::snprintf(buf, sizeof(buf), "Pan (%.0f, %.0f)", app.hud_pan_value.x, app.hud_pan_value.y);
                    ImU32 bg = IM_COL32(20, 20, 30, (int)(200 * app.hud_pan_timer));
                    ImU32 fg = IM_COL32(230, 230, 240, (int)(255 * app.hud_pan_timer));
                    ImVec2 text_sz = ImGui::CalcTextSize(buf);
                    ImVec2 pad(8.0f, 4.0f);
                    ImVec2 box_min = hud_pos;
                    ImVec2 box_max = ImVec2(hud_pos.x + text_sz.x + pad.x * 2.0f,
                                            hud_pos.y + text_sz.y + pad.y * 2.0f);
                    dl->AddRectFilled(box_min, box_max, bg, 4.0f);
                    dl->AddText(ImVec2(box_min.x + pad.x, box_min.y + pad.y), fg, buf);
                }

                // Optional axis crosshair at center
                if (app.show_axis_overlay && active_vp) {
                    ImVec2 vp_origin = ImVec2(app.viewport_pos.x + active_vp->pos.x,
                                              app.viewport_pos.y + active_vp->pos.y);
                    ImVec2 center = ImVec2(vp_origin.x + active_vp->size.x * 0.5f,
                                           vp_origin.y + active_vp->size.y * 0.5f);
                    ImU32 axis_col = IM_COL32(120, 180, 255, 100);
                    dl->AddLine(ImVec2(vp_origin.x, center.y),
                                ImVec2(vp_origin.x + active_vp->size.x, center.y),
                                axis_col, 1.0f);
                    dl->AddLine(ImVec2(center.x, vp_origin.y),
                                ImVec2(center.x, vp_origin.y + active_vp->size.y),
                                axis_col, 1.0f);
                }

                // Rulers, labels, and origin marker
                if (sim && active_vp && active_vp->show_grid) {
                    int grid_step = compute_grid_step_from_scale(active_scale);
                    if (grid_step < 1) grid_step = 1;
                    ImVec2 offset = compute_viewport_offset(*active_vp, sim, active_scale);
                    ImVec2 vp_origin = ImVec2(app.viewport_pos.x + active_vp->pos.x,
                                              app.viewport_pos.y + active_vp->pos.y);
                    float start_x_px = vp_origin.x + offset.x;
                    float start_y_px = vp_origin.y + offset.y;
                    float end_x_px = vp_origin.x + active_vp->size.x;
                    float end_y_px = vp_origin.y + active_vp->size.y;

                    if (start_x_px >= vp_origin.x - 12.0f &&
                        start_x_px <= end_x_px + 12.0f &&
                        start_y_px >= vp_origin.y - 12.0f &&
                        start_y_px <= end_y_px + 12.0f) {
                        ImU32 origin_col = IM_COL32(255, 80, 80, 220);
                        dl->AddLine(ImVec2(start_x_px - 10.0f, start_y_px),
                                    ImVec2(start_x_px + 10.0f, start_y_px),
                                    origin_col, 2.0f);
                        dl->AddLine(ImVec2(start_x_px, start_y_px - 10.0f),
                                    ImVec2(start_x_px, start_y_px + 10.0f),
                                    origin_col, 2.0f);
                    }

                    const float tick_short = 6.0f;
                    const float tick_long = 10.0f;
                    const float pad = 4.0f;
                    double cell_w_m = (sim->nx > 0) ? (sim->lx / (double)sim->nx) : 0.0;
                    double cell_h_m = (sim->ny > 0) ? (sim->ly / (double)sim->ny) : 0.0;

                    auto format_dist = [](double meters, char* out, size_t sz) {
                        if (meters >= 1.0) {
                            std::snprintf(out, sz, "%.2f m", meters);
                        } else if (meters >= 0.01) {
                            std::snprintf(out, sz, "%.1f cm", meters * 100.0);
                        } else {
                            std::snprintf(out, sz, "%.0f mm", meters * 1000.0);
                        }
                    };

                    int first_cell_x = (int)std::floor((-offset.x) / (float)active_scale);
                    int last_cell_x = (int)std::ceil((active_vp->size.x - offset.x) / (float)active_scale);
                    first_cell_x = std::max(0, first_cell_x);
                    last_cell_x = std::min(sim->nx, last_cell_x);

                    int first_cell_y = (int)std::floor((-offset.y) / (float)active_scale);
                    int last_cell_y = (int)std::ceil((active_vp->size.y - offset.y) / (float)active_scale);
                    first_cell_y = std::max(0, first_cell_y);
                    last_cell_y = std::min(sim->ny, last_cell_y);

                    for (int cx = first_cell_x - (first_cell_x % grid_step); cx <= last_cell_x; cx += grid_step) {
                        float x = vp_origin.x + offset.x + (float)cx * (float)active_scale;
                        if (x < vp_origin.x || x > end_x_px) continue;
                        bool major = ((cx / grid_step) % 5) == 0;
                        float tick = major ? tick_long : tick_short;
                        dl->AddLine(ImVec2(x, vp_origin.y),
                                    ImVec2(x, vp_origin.y + tick),
                                    IM_COL32(180, 180, 190, 190),
                                    1.0f);
                        if (major && cell_w_m > 0.0) {
                            char buf[32];
                            format_dist(cell_w_m * (double)cx, buf, sizeof(buf));
                            ImVec2 ts = ImGui::CalcTextSize(buf);
                            ImVec2 pos = ImVec2(x - ts.x * 0.5f, vp_origin.y + tick + pad);
                            dl->AddText(pos, IM_COL32(210, 210, 220, 200), buf);
                        }
                    }

                    for (int cy = first_cell_y - (first_cell_y % grid_step); cy <= last_cell_y; cy += grid_step) {
                        float y = vp_origin.y + offset.y + (float)cy * (float)active_scale;
                        if (y < vp_origin.y || y > end_y_px) continue;
                        bool major = ((cy / grid_step) % 5) == 0;
                        float tick = major ? tick_long : tick_short;
                        dl->AddLine(ImVec2(vp_origin.x, y),
                                    ImVec2(vp_origin.x + tick, y),
                                    IM_COL32(180, 180, 190, 190),
                                    1.0f);
                        if (major && cell_h_m > 0.0) {
                            char buf[32];
                            format_dist(cell_h_m * (double)cy, buf, sizeof(buf));
                            ImVec2 ts = ImGui::CalcTextSize(buf);
                            ImVec2 pos = ImVec2(vp_origin.x + tick + pad,
                                                y - ts.y * 0.5f);
                            dl->AddText(pos, IM_COL32(210, 210, 220, 200), buf);
                        }
                    }
                }

                // Area/annotation overlays
                if (sim && app.viewport_valid && active_vp) {
                    ImVec2 offset = compute_viewport_offset(*active_vp, sim, active_scale);
                    int mx = 0, my = 0;
                    SDL_GetMouseState(&mx, &my);
                    float local_x = (float)mx - (app.viewport_pos.x + active_vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + active_vp->pos.y);
                    bool mouse_in_view =
                        (local_x >= 0.0f && local_y >= 0.0f &&
                         local_x < active_vp->size.x && local_y < active_vp->size.y);

                    ImVec2 vp_origin = ImVec2(app.viewport_pos.x + active_vp->pos.x,
                                              app.viewport_pos.y + active_vp->pos.y);
                    auto grid_to_screen = [&](const ImVec2& g) {
                        return ImVec2(vp_origin.x + offset.x + g.x * (float)active_scale,
                                      vp_origin.y + offset.y + g.y * (float)active_scale);
                    };

                    if (app.show_area_measurements) {
                        ImU32 area_fill = IM_COL32(80, 180, 80, 80);
                        ImU32 area_line = IM_COL32(100, 255, 100, 200);
                        ImU32 vert_col = IM_COL32(255, 255, 120, 255);
                        for (const auto& a : app.measurements.areas) {
                            if (a.vertices.size() < 3) continue;
                            std::vector<ImVec2> pts;
                            pts.reserve(a.vertices.size());
                            for (const auto& v : a.vertices) pts.push_back(grid_to_screen(v));
                            dl->AddConvexPolyFilled(pts.data(), (int)pts.size(), area_fill);
                            for (size_t i = 0; i < pts.size(); ++i) {
                                size_t j = (i + 1) % pts.size();
                                dl->AddLine(pts[i], pts[j], area_line, 2.0f);
                            }
                            for (const auto& p : pts) dl->AddCircleFilled(p, 4.0f, vert_col);
                            ImVec2 centroid(0, 0);
                            for (const auto& p : pts) { centroid.x += p.x; centroid.y += p.y; }
                            centroid.x /= (float)pts.size();
                            centroid.y /= (float)pts.size();
                            char buf[96];
                            std::snprintf(buf, sizeof(buf), "%.6f m^2\n%.4f m", a.area_m2, a.perimeter_m);
                            ImVec2 ts = ImGui::CalcTextSize(buf);
                            ImVec2 pad(4.0f, 3.0f);
                            ImVec2 b0 = ImVec2(centroid.x - ts.x * 0.5f - pad.x,
                                               centroid.y - ts.y * 0.5f - pad.y);
                            ImVec2 b1 = ImVec2(centroid.x + ts.x * 0.5f + pad.x,
                                               centroid.y + ts.y * 0.5f + pad.y);
                            dl->AddRectFilled(b0, b1, IM_COL32(0, 0, 0, 200), 3.0f);
                            dl->AddText(ImVec2(b0.x + pad.x, b0.y + pad.y), IM_COL32(120, 255, 120, 255), buf);
                        }

                        if (app.area_mode) {
                            ImVec2 mouse_pos = io.MousePos;
                            float lx = mouse_pos.x - (app.viewport_pos.x + active_vp->pos.x);
                            float ly = mouse_pos.y - (app.viewport_pos.y + active_vp->pos.y);
                            ImVec2 preview = ImVec2(0, 0);
                            bool preview_valid = false;
                            if (lx >= 0 && ly >= 0 && lx < active_vp->size.x && ly < active_vp->size.y) {
                                float fx = lx - offset.x;
                                float fy = ly - offset.y;
                                int ix = (int)(fx / (float)active_scale);
                                int iy = (int)(fy / (float)active_scale);
                                if (ix >= 0 && iy >= 0 && ix < sim->nx && iy < sim->ny) {
                                    preview = ImVec2((float)ix, (float)iy);
                                    preview_valid = true;
                                }
                            }
                            std::vector<ImVec2> pts;
                            pts.reserve(app.current_area.vertices.size() + (preview_valid ? 1 : 0));
                            for (const auto& v : app.current_area.vertices) pts.push_back(grid_to_screen(v));
                            if (preview_valid) pts.push_back(grid_to_screen(preview));
                            for (size_t i = 0; i + 1 < pts.size(); ++i) {
                                dl->AddLine(pts[i], pts[i + 1], area_line, 2.0f);
                            }
                            for (const auto& p : pts) dl->AddCircleFilled(p, 5.0f, vert_col);

                            // Instruction badge
                            const char* msg = "Area: click vertices, right-click to finish, Esc to cancel";
                            ImVec2 ts = ImGui::CalcTextSize(msg);
                            ImVec2 pad(8.0f, 4.0f);
                            ImVec2 pos = ImVec2(vp_origin.x + 12.0f, vp_origin.y + 12.0f);
                            dl->AddRectFilled(pos,
                                              ImVec2(pos.x + ts.x + pad.x * 2.0f,
                                                     pos.y + ts.y + pad.y * 2.0f),
                                              IM_COL32(10, 10, 14, 200),
                                              4.0f);
                            dl->AddText(ImVec2(pos.x + pad.x, pos.y + pad.y),
                                        IM_COL32(230, 230, 240, 255),
                                        msg);
                        }
                    }

                    if (app.show_annotations) {
                        for (const auto& ann : app.measurements.annotations) {
                            if (!ann.visible) continue;
                            ImVec2 sp = grid_to_screen(ann.grid_pos);
                            ImVec2 ts = ImGui::CalcTextSize(ann.text);
                            ImVec2 pad(6.0f, 4.0f);
                            ImVec2 b0 = ImVec2(sp.x - pad.x, sp.y - pad.y);
                            ImVec2 b1 = ImVec2(sp.x + ts.x + pad.x, sp.y + ts.y + pad.y);
                            dl->AddRectFilled(b0, b1, IM_COL32(0, 0, 0, 200), 4.0f);
                            dl->AddRect(b0, b1, ImColor(ann.color), 4.0f, 0, 2.0f);
                            dl->AddText(sp, ImColor(ann.color), ann.text);
                        }
                    }
                }

                // Ruler overlays & cursor readout
                if (sim && app.viewport_valid && active_vp) {
                    ImVec2 offset = compute_viewport_offset(*active_vp, sim, active_scale);
                    int mx = 0, my = 0;
                    SDL_GetMouseState(&mx, &my);
                    float local_x = (float)mx - (app.viewport_pos.x + active_vp->pos.x);
                    float local_y = (float)my - (app.viewport_pos.y + active_vp->pos.y);
                    bool mouse_in_view =
                        (local_x >= 0.0f && local_y >= 0.0f &&
                         local_x < active_vp->size.x && local_y < active_vp->size.y);

                    ImVec2 vp_origin = ImVec2(app.viewport_pos.x + active_vp->pos.x,
                                              app.viewport_pos.y + active_vp->pos.y);
                    auto grid_to_screen = [&](const ImVec2& g) {
                        return ImVec2(vp_origin.x + offset.x + g.x * (float)active_scale,
                                      vp_origin.y + offset.y + g.y * (float)active_scale);
                    };

                    ImVec2 live_cursor_grid(0.0f, 0.0f);
                    bool cursor_on_grid = false;
                    if (mouse_in_view) {
                        float field_x = local_x - offset.x;
                        float field_y = local_y - offset.y;
                        int ix = (int)(field_x / (float)active_scale);
                        int iy = (int)(field_y / (float)active_scale);
                        if (ix >= 0 && iy >= 0 && ix < sim->nx && iy < sim->ny) {
                            live_cursor_grid = ImVec2((float)ix, (float)iy);
                            cursor_on_grid = true;
                        }
                    }

                    if (app.ruler_mode) {
                        ImU32 ruler_col = IM_COL32(255, 220, 80, 255);
                        ImU32 ruler_text = IM_COL32(255, 230, 140, 255);
                        if (app.ruler_first_point_set) {
                            ImVec2 a_grid = app.ruler_point_a;
                            ImVec2 b_grid = cursor_on_grid ? live_cursor_grid : app.ruler_point_b;
                            ImVec2 a_px = grid_to_screen(a_grid);
                            ImVec2 b_px = grid_to_screen(b_grid);
                            dl->AddLine(a_px, b_px, ruler_col, 2.0f);

                            double dx_m = ((double)b_grid.x - (double)a_grid.x) *
                                          (sim->lx / (double)sim->nx);
                            double dy_m = ((double)b_grid.y - (double)a_grid.y) *
                                          (sim->ly / (double)sim->ny);
                            double distance_m = std::sqrt(dx_m * dx_m + dy_m * dy_m);
                            double angle_deg = std::atan2(dy_m, dx_m) * 180.0 / IM_PI;
                            char buf[64];
                            std::snprintf(buf, sizeof(buf), "%.3f m @ %.1f deg", distance_m, angle_deg);
                            ImVec2 mid = ImVec2((a_px.x + b_px.x) * 0.5f, (a_px.y + b_px.y) * 0.5f);
                            ImVec2 ts = ImGui::CalcTextSize(buf);
                            ImVec2 box_min = ImVec2(mid.x - ts.x * 0.5f - 6.0f,
                                                    mid.y - ts.y * 0.5f - 4.0f);
                            ImVec2 box_max = ImVec2(mid.x + ts.x * 0.5f + 6.0f,
                                                    mid.y + ts.y * 0.5f + 4.0f);
                            dl->AddRectFilled(box_min, box_max, IM_COL32(20, 20, 20, 200), 4.0f);
                            dl->AddText(ImVec2(box_min.x + 6.0f, box_min.y + 4.0f), ruler_text, buf);
                        } else {
                            if (cursor_on_grid) {
                                ImVec2 hint_pos = grid_to_screen(live_cursor_grid);
                                dl->AddCircleFilled(hint_pos, 4.0f, ruler_col, 16);
                                dl->AddText(ImVec2(hint_pos.x + 10.0f, hint_pos.y - 6.0f),
                                            ruler_text,
                                            "Ruler: click first point (R to toggle)");
                            }
                        }
                    }

                    if (cursor_on_grid) {
                        double x_m = (double)live_cursor_grid.x * (sim->lx / (double)sim->nx);
                        double y_m = (double)live_cursor_grid.y * (sim->ly / (double)sim->ny);
                        char info[96];
                        std::snprintf(info,
                                      sizeof(info),
                                      "Cell: (%d, %d)  |  %.4f m, %.4f m",
                                      (int)live_cursor_grid.x,
                                      (int)live_cursor_grid.y,
                                      x_m,
                                      y_m);
                        ImVec2 info_sz = ImGui::CalcTextSize(info);
                        ImVec2 info_pad(8.0f, 5.0f);
                        ImVec2 rect_min = ImVec2(vp_origin.x + 12.0f,
                                                 vp_origin.y + active_vp->size.y - info_sz.y - info_pad.y * 2.0f - 12.0f);
                        ImVec2 rect_max = ImVec2(rect_min.x + info_sz.x + info_pad.x * 2.0f,
                                                 rect_min.y + info_sz.y + info_pad.y * 2.0f);
                        dl->AddRectFilled(rect_min, rect_max, IM_COL32(10, 10, 14, 200), 4.0f);
                        dl->AddText(ImVec2(rect_min.x + info_pad.x, rect_min.y + info_pad.y),
                                    IM_COL32(230, 230, 240, 255),
                                    info);
                    }
                }

                if (app.show_context_menu) {
                    ImGui::OpenPopup("ViewportContextMenu");
                    ImGui::SetNextWindowPos(app.context_menu_pos, ImGuiCond_Always);
                    app.show_context_menu = false; // consume request
                }
                if (ImGui::BeginPopup("ViewportContextMenu")) {
                    int cx = app.context_menu_cell_i;
                    int cy = app.context_menu_cell_j;
                    ImGui::Text("Cell (%d, %d)", cx, cy);
                    ImGui::Separator();
                    if (paused && app.viewport_valid) {
                        if (ImGui::MenuItem("Save viewport snapshot (PNG)")) {
                            std::filesystem::create_directories("recordings");
                            time_t now = time(nullptr);
                            struct tm* tm_info = localtime(&now);
                            char fname[128];
                            std::snprintf(fname,
                                          sizeof(fname),
                                          "viewport_%04d%02d%02d_%02d%02d%02d.png",
                                          tm_info->tm_year + 1900,
                                          tm_info->tm_mon + 1,
                                          tm_info->tm_mday,
                                          tm_info->tm_hour,
                                          tm_info->tm_min,
                                          tm_info->tm_sec);
                            std::filesystem::path outp = std::filesystem::path("recordings") / fname;
                            SDL_Surface* snap = capture_frame(render->renderer,
                                                              (int)app.viewport_size.x,
                                                              (int)app.viewport_size.y,
                                                              1.0f);
                            if (snap) {
                                SDL_Surface* conv = SDL_ConvertSurfaceFormat(snap, SDL_PIXELFORMAT_RGBA32, 0);
                                SDL_FreeSurface(snap);
                                if (conv) {
                                    int ok = stbi_write_png(outp.string().c_str(),
                                                            conv->w,
                                                            conv->h,
                                                            4,
                                                            conv->pixels,
                                                            conv->pitch);
                                    SDL_FreeSurface(conv);
                                    if (ok == 1) {
                                        std::snprintf(app.composer_status,
                                                      sizeof(app.composer_status),
                                                      "Saved viewport snapshot: %s",
                                                      outp.string().c_str());
                                        ui_log_add(&app, "Saved viewport snapshot: %s", outp.string().c_str());
                                        std::snprintf(app.composer_last_export_path,
                                                      sizeof(app.composer_last_export_path),
                                                      "%s",
                                                      outp.parent_path().string().c_str());
                                    } else {
                                        ui_log_add(&app, "Viewport snapshot failed: %s", outp.string().c_str());
                                    }
                                }
                            }
                        }
                        ImGui::Separator();
                    }
                    if (paused && ImGui::MenuItem("Save viewport snapshot (PNG)")) {
                        std::filesystem::create_directories("recordings");
                        time_t now = time(nullptr);
                        struct tm* tm_info = localtime(&now);
                        char fname[128];
                        std::snprintf(fname,
                                      sizeof(fname),
                                      "viewport_%04d%02d%02d_%02d%02d%02d.png",
                                      tm_info->tm_year + 1900,
                                      tm_info->tm_mon + 1,
                                      tm_info->tm_mday,
                                      tm_info->tm_hour,
                                      tm_info->tm_min,
                                      tm_info->tm_sec);
                        std::filesystem::path outp = std::filesystem::path("recordings") / fname;
                        SDL_Surface* snap = capture_frame(render->renderer,
                                                          (int)app.viewport_size.x,
                                                          (int)app.viewport_size.y,
                                                          1.0f);
                        if (snap) {
                            if (save_surface_png(outp, snap)) {
                                std::snprintf(app.composer_status,
                                              sizeof(app.composer_status),
                                              "Saved viewport snapshot: %s",
                                              outp.string().c_str());
                                ui_log_add(&app, "Saved viewport snapshot: %s", outp.string().c_str());
                                std::snprintf(app.composer_last_export_path,
                                              sizeof(app.composer_last_export_path),
                                              "%s",
                                              outp.parent_path().string().c_str());
                            } else {
                                ui_log_add(&app, "Viewport snapshot failed: %s", outp.string().c_str());
                            }
                            SDL_FreeSurface(snap);
                        }
                    }
                    if (ImGui::MenuItem("Add Source Here")) {
                        create_new_source_at(&wizard, sim, &app, cx, cy);
                        app.show_context_menu = false;
                    }
                    if (ImGui::MenuItem("Add Material Block Here")) {
                        create_block_at(&wizard, &bootstrap, sim, &app, cx, cy);
                    }
                    ImGui::Separator();
                    if (ImGui::MenuItem("Measure Distance (R)")) {
                        app.ruler_mode = true;
                        app.ruler_first_point_set = false;
                    }
                    if (ImGui::MenuItem("Zoom to Fit")) {
                        if (app.viewport_valid && sim && active_vp) {
                            float fit_scale_x = active_vp->size.x / (float)sim->nx;
                            float fit_scale_y = active_vp->size.y / (float)sim->ny;
                            float fit_zoom = std::max(0.5f, std::min(fit_scale_x, fit_scale_y));
                            active_vp->zoom = fit_zoom;
                            active_scale = (int)std::lround(fit_zoom);
                            if (active_scale < 1) active_scale = 1;
                            render->scale = active_scale;
                            active_vp->pan_x = 0.0f;
                            active_vp->pan_y = 0.0f;
                            clamp_pan_offset(active_vp, sim, active_scale);
                            scale = active_scale;
                            app.hud_zoom_value = active_vp->zoom;
                            app.hud_pan_value = ImVec2(active_vp->pan_x, active_vp->pan_y);
                            app.hud_zoom_timer = 1.0f;
                            app.hud_pan_timer = 1.0f;
                        }
                    }
                    if (app.area_mode) {
                        if (ImGui::MenuItem("Finish Area")) {
                            if (app.current_area.vertices.size() >= 3) {
                                close_area_measurement(&app, sim);
                            } else {
                                ui_log_add(&app, "Area tool: need at least 3 vertices");
                            }
                        }
                        if (ImGui::MenuItem("Cancel Area")) {
                            app.area_mode = false;
                            app.current_area.vertices.clear();
                            app.current_area.closed = false;
                            ui_log_add(&app, "Area tool: cancelled");
                        }
                    }
                    if (ImGui::MenuItem("Copy Coordinates")) {
                        char buf[64];
                        std::snprintf(buf, sizeof(buf), "(%d, %d)", cx, cy);
                        SDL_SetClipboardText(buf);
                    }
                    ImGui::EndPopup();
                } else {
                    app.show_context_menu = false;
                }
            }
            ImGui::End();
            if (frame_counter <= 120) {
                debug_logf("frame %d: overlay window end", frame_counter);
            }
        }

        app.paint_material_id = paint_material_id;
        draw_material_browser(&app, &paint_material_id);
        app.paint_material_id = paint_material_id;

        draw_sparameter_window(&app);
        draw_smith_chart(&app);

        if (frame_counter <= 120) {
            debug_logf("frame %d: before imgui render", frame_counter);
        }
        ImGui::Render();
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), render->renderer);
        if (app.composer_request_animation) {
            bool ok = export_composer_animation(&app, render, sim, &scope, app.composer_request_page, false);
            if (!ok) {
                std::snprintf(app.composer_status,
                              sizeof(app.composer_status),
                              "Animation export failed (page %d)",
                              app.composer_request_page + 1);
                ui_log_add(&app, "Composer animation export failed for page %d", app.composer_request_page + 1);
            }
            app.composer_request_animation = false;
        } else if (app.composer_request_export) {
            bool ok = export_composer_page(&app, render, sim, &scope, app.composer_request_page);
            if (!ok) {
                std::snprintf(app.composer_status,
                              sizeof(app.composer_status),
                              "Export failed (page %d)",
                              app.composer_request_page + 1);
                ui_log_add(&app, "Composer export failed for page %d", app.composer_request_page + 1);
            }
            app.composer_request_export = false;
        }
        if (app.composer_request_export_all) {
            for (int pi = 0; pi < (int)app.composer_pages.size(); ++pi) {
                bool ok = export_composer_page(&app, render, sim, &scope, pi);
                if (!ok) {
                    ui_log_add(&app, "Composer export failed for page %d", pi + 1);
                }
            }
            std::snprintf(app.composer_status,
                          sizeof(app.composer_status),
                          "Exported all pages (%zu)",
                          app.composer_pages.size());
            app.composer_request_export_all = false;
        }
        SDL_RenderPresent(render->renderer);

        if (frame_counter <= 120) {
            debug_logf("frame %d: after present", frame_counter);
        }
        SDL_Delay(UI_DELAY_MS);
    }

    scope_free(&scope);
    ImGui_ImplSDLRenderer2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    render_free(render);
    material_library_shutdown();
    simulation_bootstrap_shutdown(&bootstrap);
    if (app.composer_preview_tex) {
        SDL_DestroyTexture(app.composer_preview_tex);
        app.composer_preview_tex = nullptr;
    }
    debug_log_close();

    return 0;
}
