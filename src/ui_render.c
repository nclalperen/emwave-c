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

    const char* font_paths[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "DejaVuSans.ttf",
        NULL
    };

    for (int i = 0; font_paths[i]; ++i) {
        ctx->font = TTF_OpenFont(font_paths[i], 14);
        if (ctx->font) break;
    }

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
    ctx->ui_height = RENDER_DEFAULT_UI_HEIGHT;
    ctx->side_panel_width = RENDER_DEFAULT_SIDE_PANEL;

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

static void draw_grid(RenderContext* ctx, const SimulationState* state, double vmax) {
    if (!ctx || !state) return;
    SDL_Rect pixel = {0, 0, ctx->scale, ctx->scale};

    for (int i = 0; i < state->nx; ++i) {
        for (int j = 0; j < state->ny; ++j) {
            double v = state->Ez[i][j];
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
}

void render_field_heatmap(RenderContext* ctx, const SimulationState* state,
                          double vmax, double color_scale) {
    (void)color_scale;
    draw_grid(ctx, state, vmax);
}

void render_sources(RenderContext* ctx, const Source* sources) {
    if (!ctx || !sources) return;
    SDL_SetRenderDrawColor(ctx->renderer, 255, 255, 255, 255);
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
    SDL_SetRenderDrawColor(ctx->renderer, 220, 220, 220, 255);
    SDL_RenderDrawRect(ctx->renderer, &r);
}

void render_colorbar(RenderContext* ctx, const RenderLayout* layout, double vmin, double vmax) {
    if (!ctx || !layout) return;
    SDL_Rect bar = { layout->colorbar.x, layout->colorbar.y, layout->colorbar.w, layout->colorbar.h };
    if (bar.h <= 1 || bar.w <= 0) return;
    for (int y = 0; y < bar.h; ++y) {
        double t = 1.0 - (double)y / (double)(bar.h - 1);
        double v = vmin + t * (vmax - vmin);
        SDL_Color c = colormap_heat(v, vmax);
        SDL_SetRenderDrawColor(ctx->renderer, c.r, c.g, c.b, c.a);
        SDL_RenderDrawLine(ctx->renderer, bar.x, bar.y + y, bar.x + bar.w, bar.y + y);
    }
}

void render_scope(RenderContext* ctx, const Scope* scope, int x, int y, int w, int h, double yscale) {
    if (!scope || !scope->y || scope->n <= 0) return;
    SDL_Rect frame = { x, y, w, h };
    SDL_SetRenderDrawColor(ctx->renderer, 30, 30, 30, 255);
    SDL_RenderFillRect(ctx->renderer, &frame);
    SDL_SetRenderDrawColor(ctx->renderer, 80, 80, 80, 255);
    SDL_RenderDrawRect(ctx->renderer, &frame);

    int samples = scope->n;
    double vmax = (yscale > 0.0) ? yscale : 1.0;

    for (int k = 1; k < samples; ++k) {
        int idx0 = (scope->head + k - 1) % samples;
        int idx1 = (scope->head + k) % samples;
        double v0 = scope->y[idx0];
        double v1 = scope->y[idx1];
        int px0 = x + (k - 1) * w / samples;
        int px1 = x + k * w / samples;
        int py0 = y + h / 2 - (int)((v0 / vmax) * (h / 2));
        int py1 = y + h / 2 - (int)((v1 / vmax) * (h / 2));
        SDL_SetRenderDrawColor(ctx->renderer, 0, 220, 120, 255);
        SDL_RenderDrawLine(ctx->renderer, px0, py0, px1, py1);
    }
}

void render_info_panel(RenderContext* ctx, const SimulationState* state, const UIState* ui,
                       double fps_avg, const RenderLayout* layout) {
    if (!ctx || !state || !ui || !layout) return;
    char buffer[256];
    snprintf(buffer, sizeof(buffer),
             "freq: %.3f GHz\nsteps: %d\nΔx=%.3f mm\nΔt=%.3e s\nFPS: %.1f",
             state->freq * 1e-9, ui->steps_per_frame,
             state->dx * 1e3, state->dt, fps_avg);
    SDL_Color c = { 230, 230, 230, 255 };
    int tw, th;
    SDL_Texture* tex = render_text_wrapped(ctx, buffer, c, ctx->side_panel_width - 30, &tw, &th);
    if (!tex) return;
    SDL_Rect dst = { layout->canvas_w + 30, 40, tw, th };
    SDL_RenderCopy(ctx->renderer, tex, NULL, &dst);
    SDL_DestroyTexture(tex);
}

void render_legend(RenderContext* ctx, int x, int y) {
    const char* legend =
        "Space: pause/resume   ESC/Q: quit\n"
        "L: legend on/off      M/U: paint mode\n"
        "I: paint type  O/P: paint eps -/+\n"
        "1/2/3: toggle sources  T: cycle src type\n"
        "C: clear fields       Y: CPML/Mur\n"
        "H/J: hold color/scope S: ports  G: log probe\n"
        "Bottom sliders: frequency and steps per frame";
    SDL_Color c = { 180, 180, 180, 255 };
    int tw, th;
    SDL_Texture* tex = render_text_wrapped(ctx, legend, c, ctx->side_panel_width - 30, &tw, &th);
    if (!tex) return;
    SDL_Rect dst = { x, y, tw, th };
    SDL_RenderCopy(ctx->renderer, tex, NULL, &dst);
    SDL_DestroyTexture(tex);
}

void slider_draw(RenderContext* ctx, const Slider* slider) {
    if (!ctx || !slider) return;
    SDL_Rect track = { slider->x, slider->y + slider->h / 2 - 2, slider->w, 4 };
    SDL_SetRenderDrawColor(ctx->renderer, 80, 80, 80, 255);
    SDL_RenderFillRect(ctx->renderer, &track);
    double t = (slider->value - slider->minv) / (slider->maxv - slider->minv);
    int kx = slider->x + (int)(t * slider->w);
    SDL_Rect knob = { kx - 6, slider->y, 12, slider->h };
    SDL_SetRenderDrawColor(ctx->renderer, 200, 200, 200, 255);
    SDL_RenderFillRect(ctx->renderer, &knob);
}

void render_frame(RenderContext* ctx, const SimulationState* state, const UIState* ui,
                  const Scope* scope, double fps_avg) {
    if (!ctx || !state || !ui) return;

    RenderLayout layout;
    render_layout_compute(&layout, state->nx, state->ny, ctx->scale,
                          ctx->side_panel_width, ctx->ui_height);

    SDL_RenderSetLogicalSize(ctx->renderer, layout.window_w, layout.window_h);

    SDL_SetRenderDrawColor(ctx->renderer, 8, 8, 12, 255);
    SDL_RenderClear(ctx->renderer);

    double vmax = ui->hold_color ? ui->held_vmax : ui->vmax_smooth;
    if (vmax < 1e-6) vmax = 1e-6;

    render_field_heatmap(ctx, state, vmax, 1.0);
    render_sources(ctx, state->sources);
    render_block_outline(ctx, &layout);

    SDL_Rect panel = { layout.canvas_w, 0, ctx->side_panel_width, layout.window_h };
    SDL_SetRenderDrawColor(ctx->renderer, 20, 20, 30, 255);
    SDL_RenderFillRect(ctx->renderer, &panel);

    render_colorbar(ctx, &layout, -vmax, vmax);
    render_info_panel(ctx, state, ui, fps_avg, &layout);

    int scope_x = layout.canvas_w + 30;
    int scope_y = layout.canvas_h - 160;
    if (scope_y < 20) scope_y = 20;
    render_scope(ctx, scope, scope_x, scope_y, ctx->side_panel_width - 60, 120,
                 ui->hold_scope ? ui->held_scope_vmax : ui->scope_vmax_smooth);

    slider_draw(ctx, &ui->freq_slider);
    slider_draw(ctx, &ui->speed_slider);

    if (ui->show_legend) {
        render_legend(ctx, layout.canvas_w + 30, scope_y + 140);
    }

    SDL_RenderPresent(ctx->renderer);
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
