// =============================================================================
// emwave-c: UI Layout Helpers
// =============================================================================

#ifndef EMWAVE_UI_LAYOUT_H
#define EMWAVE_UI_LAYOUT_H

#define RENDER_DEFAULT_SCALE 2
#define RENDER_DEFAULT_UI_HEIGHT 90
#define RENDER_DEFAULT_SIDE_PANEL 240

typedef struct {
    int x;
    int y;
    int w;
    int h;
} IntRect;

typedef struct {
    int canvas_w;
    int canvas_h;
    int window_w;
    int window_h;
    IntRect colorbar;
    IntRect block_outline;
} RenderLayout;

void render_layout_compute(RenderLayout* layout,
                           int nx,
                           int ny,
                           int scale,
                           int side_panel_width,
                           int ui_height);

#endif /* EMWAVE_UI_LAYOUT_H */
