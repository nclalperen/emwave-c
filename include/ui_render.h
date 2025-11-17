// ============================================================================
// emwave-c: UI Rendering Interface
// This module contains ALL SDL2 rendering code
// Perfect isolation for future UI/UX overhaul
// ============================================================================

#ifndef EMWAVE_UI_RENDER_H
#define EMWAVE_UI_RENDER_H

#include "types.h"
#include "ui_controls.h"
#include "ui_layout.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Rendering context structure */
typedef struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    TTF_Font* font;
    int scale;
    int logical_w;
    int logical_h;
    int menu_bar_height;
    int timeline_height;
    int left_panel_width;
    int right_panel_width;
} RenderContext;

/* Initialization */
RenderContext* render_init(const char* title, int width, int height);
void render_free(RenderContext* ctx);

/* Main rendering function */
void render_frame(RenderContext* ctx, const SimulationState* state, UIState* ui,
                  const Scope* scope, double fps_avg);

/* Individual rendering components */
double render_field_heatmap(RenderContext* ctx, const SimulationState* state,
                            double vmax, double color_scale);
void render_sources(RenderContext* ctx, const Source* sources);
void render_block_outline(RenderContext* ctx, const RenderLayout* layout);
void render_colorbar(RenderContext* ctx, const RenderLayout* layout,
                     ColorMapMode mode);
double render_scope(RenderContext* ctx, const Scope* scope, int x, int y, int w, int h, double yscale);
void render_info_panel(RenderContext* ctx, const SimulationState* state, const UIState* ui,
                       double fps_avg, const RenderLayout* layout);
void render_legend(RenderContext* ctx, int x, int y, int max_width);

/* Text rendering helpers */
SDL_Texture* render_text(RenderContext* ctx, const char* text, SDL_Color color, int* w, int* h);
SDL_Texture* render_text_wrapped(RenderContext* ctx, const char* text, SDL_Color color,
                                  unsigned int wrap_width, int* outw, int* outh);

/* Slider rendering */
void slider_draw(RenderContext* ctx, const Slider* slider);

/* Screenshot */
int save_screenshot(RenderContext* ctx, const char* filename);

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_UI_RENDER_H */
