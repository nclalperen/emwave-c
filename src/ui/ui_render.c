// ============================================================================
// emwave-c: SDL Rendering Helpers
// ============================================================================

#include "ui_render.h"

#include "config.h"
#include "materials.h"
#include "analysis.h"
#include "sources.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct {
    SDL_Color window_bg;
    SDL_Color menu_bg;
    SDL_Color toolbox_bg;
    SDL_Color properties_bg;
    SDL_Color timeline_bg;
    SDL_Color viewport_bg;
    SDL_Color viewport_border;
    SDL_Color timeline_border;
    SDL_Color slider_track;
    SDL_Color slider_knob;
    SDL_Color info_label;
    SDL_Color legend_text;
    SDL_Color text_muted;
    SDL_Color scope_bg;
    SDL_Color scope_border;
} ThemeColors;

typedef struct {
    SDL_Color primary;
    SDL_Color secondary;
} AccentPalette;

static const ThemeColors THEME_TABLE[THEME_COUNT] = {
    [THEME_DARK] = {
        .window_bg = {6, 8, 12, 255},
        .menu_bg = {20, 22, 34, 255},
        .toolbox_bg = {16, 20, 30, 255},
        .properties_bg = {20, 24, 36, 255},
        .timeline_bg = {12, 14, 24, 255},
        .viewport_bg = {10, 12, 20, 255},
        .viewport_border = {42, 46, 60, 255},
        .timeline_border = {48, 52, 70, 255},
        .slider_track = {70, 74, 94, 255},
        .slider_knob = {230, 230, 238, 255},
        .info_label = {194, 196, 210, 255},
        .legend_text = {182, 186, 202, 255},
        .text_muted = {170, 174, 190, 255},
        .scope_bg = {22, 24, 34, 255},
        .scope_border = {70, 74, 96, 255},
    },
    [THEME_LIGHT] = {
        .window_bg = {240, 242, 247, 255},
        .menu_bg = {224, 226, 234, 255},
        .toolbox_bg = {232, 234, 241, 255},
        .properties_bg = {236, 238, 246, 255},
        .timeline_bg = {220, 222, 232, 255},
        .viewport_bg = {245, 247, 252, 255},
        .viewport_border = {180, 184, 196, 255},
        .timeline_border = {168, 170, 185, 255},
        .slider_track = {190, 194, 210, 255},
        .slider_knob = {90, 92, 100, 255},
        .info_label = {70, 72, 80, 255},
        .legend_text = {70, 74, 86, 255},
        .text_muted = {94, 98, 110, 255},
        .scope_bg = {232, 234, 240, 255},
        .scope_border = {170, 174, 186, 255},
    }
};

static const AccentPalette ACCENT_TABLE[UI_ACCENT_PRESET_COUNT] = {
    {{0, 220, 120, 255}, {0, 169, 255, 255}},
    {{253, 126, 20, 255}, {255, 200, 67, 255}},
    {{122, 162, 247, 255}, {56, 189, 248, 255}},
    {{255, 92, 146, 255}, {255, 147, 206, 255}},
    {{180, 160, 255, 255}, {255, 255, 255, 255}},
    {{90, 220, 255, 255}, {90, 180, 255, 255}}
};

static const ThemeColors* g_active_theme = NULL;
static const AccentPalette* g_active_palette = NULL;

static void theme_set_active(ThemeMode mode, int accent_index) {
    int theme_idx = (mode >= 0 && mode < THEME_COUNT) ? mode : THEME_DARK;
    int palette_idx = accent_index;
#if UI_ACCENT_PRESET_COUNT > 0
    palette_idx %= UI_ACCENT_PRESET_COUNT;
    if (palette_idx < 0) palette_idx += UI_ACCENT_PRESET_COUNT;
#else
    palette_idx = 0;
#endif
    g_active_theme = &THEME_TABLE[theme_idx];
    g_active_palette = &ACCENT_TABLE[palette_idx];
}

static const ThemeColors* theme_colors(void) {
    if (!g_active_theme) {
        g_active_theme = &THEME_TABLE[THEME_DARK];
    }
    return g_active_theme;
}

static const AccentPalette* theme_palette(void) {
    if (!g_active_palette) {
        g_active_palette = &ACCENT_TABLE[0];
    }
    return g_active_palette;
}

static double clampd(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static SDL_Color colormap_heat(double v, double vmax) {
    double t = 0.5;
    if (vmax > 0.0) {
        t = 0.5 + 0.5 * clampd(v / vmax, -1.0, 1.0);
    }
    double r = clampd(4.0 * (t - 0.25), 0.0, 1.0);
    double g = clampd(4.0 * fabs(t - 0.5) - 1.0, 0.0, 1.0);
    double b = clampd(4.0 * (0.75 - t), 0.0, 1.0);
    SDL_Color c = { (Uint8)(r * 255.0), (Uint8)(g * 255.0), (Uint8)(b * 255.0), 255 };
    return c;
}

static void fill_panel(SDL_Renderer* renderer, const IntRect* rect,
                       Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    if (!renderer || !rect || rect->w <= 0 || rect->h <= 0) return;
    SDL_Rect sdl_rect = { rect->x, rect->y, rect->w, rect->h };
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
    SDL_RenderFillRect(renderer, &sdl_rect);
}

static SDL_Rect to_sdl_rect(const IntRect* rect) {
    SDL_Rect r = {0, 0, 0, 0};
    if (rect) {
        r.x = rect->x;
        r.y = rect->y;
        r.w = rect->w;
        r.h = rect->h;
    }
    return r;
}

RenderContext* render_init(const char* title, int width, int height) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return NULL;
    }
    if (TTF_Init() < 0) {
        fprintf(stderr, "TTF_Init failed: %s\n", TTF_GetError());
        SDL_Quit();
        return NULL;
    }

    RenderContext* ctx = (RenderContext*)calloc(1, sizeof(RenderContext));
    if (!ctx) {
        TTF_Quit();
        SDL_Quit();
        return NULL;
    }

    ctx->window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                   width, height, SDL_WINDOW_SHOWN);
    if (!ctx->window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        free(ctx);
        TTF_Quit();
        SDL_Quit();
        return NULL;
    }

    ctx->renderer = SDL_CreateRenderer(ctx->window, -1, SDL_RENDERER_ACCELERATED);
    if (!ctx->renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(ctx->window);
        free(ctx);
        TTF_Quit();
        SDL_Quit();
        return NULL;
    }

    const char* bundled_relative_font = "assets/fonts/DejaVuSans.ttf";
    char* base_path = SDL_GetBasePath();
    char* bundled_font_path = NULL;
    if (base_path) {
        size_t len = strlen(base_path) + strlen(bundled_relative_font) + 1;
        bundled_font_path = (char*)malloc(len);
        if (bundled_font_path) {
            snprintf(bundled_font_path, len, "%s%s", base_path, bundled_relative_font);
        }
        SDL_free(base_path);
    } else {
        size_t len = strlen(bundled_relative_font) + 1;
        bundled_font_path = (char*)malloc(len);
        if (bundled_font_path) {
            memcpy(bundled_font_path, bundled_relative_font, len);
        }
    }

    const char* font_paths[5] = {0};
    int font_count = 0;
    if (bundled_font_path) {
        font_paths[font_count++] = bundled_font_path;
    }
    font_paths[font_count++] = bundled_relative_font;
    font_paths[font_count++] = "DejaVuSans.ttf";
    font_paths[font_count++] = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf";

    for (int i = 0; i < font_count; ++i) {
        ctx->font = TTF_OpenFont(font_paths[i], 14);
        if (ctx->font) break;
    }

    free(bundled_font_path);

    if (!ctx->font) {
        fprintf(stderr, "TTF_OpenFont failed: %s\n", TTF_GetError());
        SDL_DestroyRenderer(ctx->renderer);
        SDL_DestroyWindow(ctx->window);
        free(ctx);
        TTF_Quit();
        SDL_Quit();
        return NULL;
    }

    ctx->scale = RENDER_DEFAULT_SCALE;
    ctx->logical_w = width;
    ctx->logical_h = height;
    ctx->menu_bar_height = RENDER_DEFAULT_MENU_BAR;
    ctx->timeline_height = RENDER_DEFAULT_TIMELINE_HEIGHT;
    ctx->left_panel_width = RENDER_DEFAULT_LEFT_PANEL;
    ctx->right_panel_width = RENDER_DEFAULT_RIGHT_PANEL;

    return ctx;
}

void render_free(RenderContext* ctx) {
    if (!ctx) return;
    if (ctx->font) TTF_CloseFont(ctx->font);
    if (ctx->renderer) SDL_DestroyRenderer(ctx->renderer);
    if (ctx->window) SDL_DestroyWindow(ctx->window);
    free(ctx);
    TTF_Quit();
    SDL_Quit();
}

SDL_Texture* render_text(RenderContext* ctx, const char* text, SDL_Color color, int* w, int* h) {
    if (!ctx || !ctx->font) return NULL;
    SDL_Surface* surf = TTF_RenderUTF8_Blended(ctx->font, text, color);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ctx->renderer, surf);
    if (tex && w && h) {
        *w = surf->w;
        *h = surf->h;
    }
    SDL_FreeSurface(surf);
    return tex;
}

SDL_Texture* render_text_wrapped(RenderContext* ctx, const char* text, SDL_Color color,
                                 unsigned int wrap_width, int* outw, int* outh) {
    if (!ctx || !ctx->font) return NULL;
    SDL_Surface* surf = TTF_RenderUTF8_Blended_Wrapped(ctx->font, text, color, wrap_width);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ctx->renderer, surf);
    if (tex && outw && outh) {
        *outw = surf->w;
        *outh = surf->h;
    }
    SDL_FreeSurface(surf);
    return tex;
}

static double draw_grid(RenderContext* ctx, const SimulationState* state, double vmax) {
    if (!ctx || !state) return 0.0;
    SDL_Rect pixel = {0, 0, ctx->scale, ctx->scale};
    double field_max = 0.0;

    for (int i = 0; i < state->nx; ++i) {
        for (int j = 0; j < state->ny; ++j) {
            double v = state->Ez[i][j];
            double abs_v = fabs(v);
            if (abs_v > field_max) {
                field_max = abs_v;
            }
            SDL_Color c = colormap_heat(v, vmax);
            if (state->tag_grid[i][j] == 1) {
                c = (SDL_Color){200, 200, 200, 255};
            } else if (state->tag_grid[i][j] == 2) {
                c = (SDL_Color){120, 180, 200, 255};
            }
            SDL_SetRenderDrawColor(ctx->renderer, c.r, c.g, c.b, c.a);
            pixel.x = i * ctx->scale;
            pixel.y = j * ctx->scale;
            SDL_RenderFillRect(ctx->renderer, &pixel);
        }
    }

    return field_max;
}

double render_field_heatmap(RenderContext* ctx, const SimulationState* state,
                            double vmax, double color_scale) {
    (void)color_scale;
    return draw_grid(ctx, state, vmax);
}

void render_sources(RenderContext* ctx, const Source* sources) {
    if (!ctx || !sources) return;
    SDL_Color accent = theme_palette()->primary;
    SDL_SetRenderDrawColor(ctx->renderer, accent.r, accent.g, accent.b, accent.a);
    for (int k = 0; k < MAX_SRC; ++k) {
        if (!sources[k].active) continue;
        int x = sources[k].ix * ctx->scale;
        int y = sources[k].iy * ctx->scale;
        SDL_RenderDrawLine(ctx->renderer, x - 4, y, x + 4, y);
        SDL_RenderDrawLine(ctx->renderer, x, y - 4, x, y + 4);
    }
}

void render_block_outline(RenderContext* ctx, const RenderLayout* layout) {
    if (!ctx || !layout) return;
    SDL_Rect r = { layout->block_outline.x, layout->block_outline.y,
                   layout->block_outline.w, layout->block_outline.h };
    if (r.w <= 0 || r.h <= 0) return;
    SDL_Color accent = theme_palette()->secondary;
    SDL_SetRenderDrawColor(ctx->renderer, accent.r, accent.g, accent.b, accent.a);
    SDL_RenderDrawRect(ctx->renderer, &r);
}

void render_colorbar(RenderContext* ctx, const RenderLayout* layout, double vmin, double vmax) {
    if (!ctx || !layout) return;
    SDL_Rect bar = { layout->colorbar.x, layout->colorbar.y, layout->colorbar.w, layout->colorbar.h };
    if (bar.h <= 1 || bar.w <= 0) return;
    SDL_Color bg = theme_colors()->properties_bg;
    SDL_SetRenderDrawColor(ctx->renderer, bg.r, bg.g, bg.b, 255);
    SDL_RenderFillRect(ctx->renderer, &bar);
    for (int y = 0; y < bar.h; ++y) {
        double t = 1.0 - (double)y / (double)(bar.h - 1);
        double v = vmin + t * (vmax - vmin);
        SDL_Color c = colormap_heat(v, vmax);
        SDL_SetRenderDrawColor(ctx->renderer, c.r, c.g, c.b, c.a);
        SDL_RenderDrawLine(ctx->renderer, bar.x, bar.y + y, bar.x + bar.w, bar.y + y);
    }
}

double render_scope(RenderContext* ctx, const Scope* scope, int x, int y, int w, int h, double yscale) {
    if (!scope || !scope->y || scope->n <= 0) return 0.0;
    SDL_Rect frame = { x, y, w, h };
    const ThemeColors* theme = theme_colors();
    SDL_SetRenderDrawColor(ctx->renderer, theme->scope_bg.r, theme->scope_bg.g, theme->scope_bg.b, 255);
    SDL_RenderFillRect(ctx->renderer, &frame);
    SDL_SetRenderDrawColor(ctx->renderer, theme->scope_border.r, theme->scope_border.g, theme->scope_border.b, 255);
    SDL_RenderDrawRect(ctx->renderer, &frame);

    int samples = scope->n;
    double vmax = (yscale > 0.0) ? yscale : 1.0;
    double scope_max = 0.0;

    for (int k = 1; k < samples; ++k) {
        int idx0 = (scope->head + k - 1) % samples;
        int idx1 = (scope->head + k) % samples;
        double v0 = scope->y[idx0];
        double v1 = scope->y[idx1];
        double abs_v0 = fabs(v0);
        double abs_v1 = fabs(v1);
        if (abs_v0 > scope_max) scope_max = abs_v0;
        if (abs_v1 > scope_max) scope_max = abs_v1;
        int px0 = x + (k - 1) * w / samples;
        int px1 = x + k * w / samples;
        int py0 = y + h / 2 - (int)((v0 / vmax) * (h / 2));
        int py1 = y + h / 2 - (int)((v1 / vmax) * (h / 2));
        SDL_Color accent = theme_palette()->primary;
        SDL_SetRenderDrawColor(ctx->renderer, accent.r, accent.g, accent.b, accent.a);
        SDL_RenderDrawLine(ctx->renderer, px0, py0, px1, py1);
    }

    return scope_max;
}

void render_info_panel(RenderContext* ctx, const SimulationState* state, const UIState* ui,
                       double fps_avg, const RenderLayout* layout) {
    if (!ctx || !state || !ui || !layout) return;
    const ThemeColors* theme = theme_colors();
    SDL_Color accent = theme_palette()->primary;
    const char* labels[] = {"freq", "steps", "dx", "dt", "fps"};
    char values[5][32];
    snprintf(values[0], sizeof(values[0]), "%.3f GHz", state->freq * 1e-9);
    snprintf(values[1], sizeof(values[1]), "%d", ui->steps_per_frame);
    snprintf(values[2], sizeof(values[2]), "%.3f mm", state->dx * 1e3);
    snprintf(values[3], sizeof(values[3]), "%.3e s", state->dt);
    snprintf(values[4], sizeof(values[4]), "%.1f", fps_avg);

    typedef struct {
        SDL_Texture* tex;
        int w;
        int h;
    } TextSprite;

    TextSprite label_tex[5] = {0};
    TextSprite value_tex[5] = {0};
    int max_label_w = 0;
    for (int i = 0; i < 5; ++i) {
        label_tex[i].tex = render_text(ctx, labels[i], theme->info_label, &label_tex[i].w, &label_tex[i].h);
        if (label_tex[i].w > max_label_w) max_label_w = label_tex[i].w;
        value_tex[i].tex = render_text(ctx, values[i], accent, &value_tex[i].w, &value_tex[i].h);
    }

    int base_x = layout->properties_panel.x + 16;
    int base_y = layout->properties_panel.y + 16;
    int cursor_y = base_y;
    const int line_gap = 6;

    for (int i = 0; i < 5; ++i) {
        int row_height = label_tex[i].h;
        if (value_tex[i].h > row_height) row_height = value_tex[i].h;

        SDL_Rect label_dst = { base_x, cursor_y, label_tex[i].w, label_tex[i].h };
        if (label_tex[i].tex) SDL_RenderCopy(ctx->renderer, label_tex[i].tex, NULL, &label_dst);

        SDL_Rect value_dst = { base_x + max_label_w + 12, cursor_y, value_tex[i].w, value_tex[i].h };
        if (value_tex[i].tex) SDL_RenderCopy(ctx->renderer, value_tex[i].tex, NULL, &value_dst);

        cursor_y += row_height + line_gap;
    }

    for (int i = 0; i < 5; ++i) {
        if (label_tex[i].tex) SDL_DestroyTexture(label_tex[i].tex);
        if (value_tex[i].tex) SDL_DestroyTexture(value_tex[i].tex);
    }
}

void render_legend(RenderContext* ctx, int x, int y, int max_width) {
    const char* legend =
        "Space: pause/resume   ESC/Q: quit\n"
        "L: legend on/off      M/U: paint mode\n"
        "I: paint type  O/P: paint eps -/+\n"
        "1/2/3: toggle sources  T: cycle src type\n"
        "C: clear fields       Y: CPML/Mur\n"
        "H/J: hold color/scope S: ports  G: log probe\n"
        "Bottom sliders: frequency and steps per frame";
    SDL_Color c = theme_colors()->legend_text;
    int tw, th;
    if (max_width <= 0) max_width = 200;
    SDL_Texture* tex = render_text_wrapped(ctx, legend, c, (unsigned int)max_width, &tw, &th);
    if (!tex) return;
    SDL_Rect dst = { x, y, tw, th };
    SDL_RenderCopy(ctx->renderer, tex, NULL, &dst);
    SDL_DestroyTexture(tex);
}

void slider_draw(RenderContext* ctx, const Slider* slider) {
    if (!ctx || !slider) return;
    const ThemeColors* theme = theme_colors();
    SDL_Color accent = theme_palette()->primary;
    SDL_Rect track = { slider->x, slider->y + slider->h / 2 - 2, slider->w, 4 };
    SDL_SetRenderDrawColor(ctx->renderer, theme->slider_track.r, theme->slider_track.g, theme->slider_track.b, 255);
    SDL_RenderFillRect(ctx->renderer, &track);
    double t = (slider->value - slider->minv) / (slider->maxv - slider->minv);
    int kx = slider->x + (int)(t * slider->w);
    SDL_Rect knob = { kx - 6, slider->y, 12, slider->h };
    SDL_SetRenderDrawColor(ctx->renderer, accent.r, accent.g, accent.b, accent.a);
    SDL_RenderFillRect(ctx->renderer, &knob);
}

void render_frame(RenderContext* ctx, const SimulationState* state, UIState* ui,
                  const Scope* scope, double fps_avg) {
    if (!ctx || !state || !ui) return;

    RenderLayout layout;
    render_layout_compute(&layout, state->nx, state->ny, ctx->scale,
                          ui->left_panel_width, ui->right_panel_width,
                          ui->timeline_height, ui->menu_bar_height);

    theme_set_active(ui->theme_mode, ui->accent_index);
    const ThemeColors* colors = theme_colors();

    if (layout.window_w != ctx->logical_w || layout.window_h != ctx->logical_h) {
        SDL_SetWindowSize(ctx->window, layout.window_w, layout.window_h);
        ctx->logical_w = layout.window_w;
        ctx->logical_h = layout.window_h;
    }
    SDL_RenderSetLogicalSize(ctx->renderer, layout.window_w, layout.window_h);

    SDL_SetRenderDrawColor(ctx->renderer, colors->window_bg.r, colors->window_bg.g, colors->window_bg.b, 255);
    SDL_RenderClear(ctx->renderer);

    fill_panel(ctx->renderer, &layout.menu_bar, colors->menu_bg.r, colors->menu_bg.g, colors->menu_bg.b, 255);
    fill_panel(ctx->renderer, &layout.toolbox_panel, colors->toolbox_bg.r, colors->toolbox_bg.g, colors->toolbox_bg.b, 255);
    fill_panel(ctx->renderer, &layout.properties_panel, colors->properties_bg.r, colors->properties_bg.g, colors->properties_bg.b, 255);
    fill_panel(ctx->renderer, &layout.timeline_panel, colors->timeline_bg.r, colors->timeline_bg.g, colors->timeline_bg.b, 255);
    SDL_SetRenderDrawColor(ctx->renderer, colors->timeline_border.r, colors->timeline_border.g, colors->timeline_border.b, 255);
    SDL_RenderDrawLine(ctx->renderer,
                       layout.timeline_panel.x,
                       layout.timeline_panel.y,
                       layout.timeline_panel.x + layout.timeline_panel.w,
                       layout.timeline_panel.y);

    SDL_SetRenderDrawColor(ctx->renderer, colors->viewport_border.r, colors->viewport_border.g, colors->viewport_border.b, 255);
    SDL_Rect view_border = to_sdl_rect(&layout.viewport);
    SDL_RenderDrawRect(ctx->renderer, &view_border);

    const char* menu = "File   Edit   View   Simulation   Help";
    int menu_w = 0, menu_h = 0;
    SDL_Texture* menu_tex = render_text(ctx, menu, colors->text_muted, &menu_w, &menu_h);
    if (menu_tex) {
        int menu_x = layout.menu_bar.x + 18;
        int menu_y = layout.menu_bar.y + (layout.menu_bar.h - menu_h) / 2;
        SDL_Rect dst = { menu_x, menu_y, menu_w, menu_h };
        SDL_RenderCopy(ctx->renderer, menu_tex, NULL, &dst);
        SDL_DestroyTexture(menu_tex);
    }

    double vmax = ui->hold_color ? ui->held_vmax : ui->vmax_smooth;
    if (vmax < 1e-6) vmax = 1e-6;

    SDL_Rect viewport_rect = to_sdl_rect(&layout.viewport);
    SDL_RenderSetViewport(ctx->renderer, &viewport_rect);
    SDL_SetRenderDrawColor(ctx->renderer, colors->viewport_bg.r, colors->viewport_bg.g, colors->viewport_bg.b, 255);
    SDL_RenderClear(ctx->renderer);
    double field_vmax = render_field_heatmap(ctx, state, vmax, 1.0);
    render_sources(ctx, state->sources);
    render_block_outline(ctx, &layout);
    SDL_RenderSetViewport(ctx->renderer, NULL);

    render_colorbar(ctx, &layout, -vmax, vmax);
    render_info_panel(ctx, state, ui, fps_avg, &layout);

    double scope_vmax = 0.0;
    const IntRect* dock_rect = (ui->scope_dock == SCOPE_DOCK_TIMELINE)
                                   ? &layout.timeline_scope
                                   : &layout.properties_scope;
    if (scope && dock_rect && dock_rect->w > 40 && dock_rect->h > 40) {
        scope_vmax = render_scope(ctx, scope,
                                  dock_rect->x, dock_rect->y,
                                  dock_rect->w, dock_rect->h,
                                  ui->hold_scope ? ui->held_scope_vmax : ui->scope_vmax_smooth);
    }

    slider_draw(ctx, &ui->freq_slider);
    slider_draw(ctx, &ui->speed_slider);

    if (ui->show_legend) {
        int legend_w = layout.toolbox_panel.w - 24;
        if (legend_w <= 0) legend_w = 200;
        int legend_x = layout.toolbox_panel.x + 12;
        int legend_y = layout.toolbox_panel.y + layout.toolbox_panel.h - 150;
        if (legend_y < layout.toolbox_panel.y + 12) {
            legend_y = layout.toolbox_panel.y + 12;
        }
        render_legend(ctx, legend_x, legend_y, legend_w);
    }

    SDL_RenderPresent(ctx->renderer);

    if (ui) {
        ui->render_field_vmax = field_vmax;
        ui->render_scope_vmax = scope_vmax;
        ui->render_metrics_valid = 1;
    }
}

int save_screenshot(RenderContext* ctx, const char* filename) {
    if (!ctx || !ctx->renderer) return 0;
    int w, h;
    SDL_GetRendererOutputSize(ctx->renderer, &w, &h);
    SDL_Surface* surface = SDL_CreateRGBSurfaceWithFormat(0, w, h, 24, SDL_PIXELFORMAT_RGB24);
    if (!surface) return 0;
    if (SDL_RenderReadPixels(ctx->renderer, NULL, surface->format->format, surface->pixels, surface->pitch) != 0) {
        SDL_FreeSurface(surface);
        return 0;
    }
    int ok = SDL_SaveBMP(surface, filename);
    SDL_FreeSurface(surface);
    return (ok == 0);
}
