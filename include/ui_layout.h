// =============================================================================
// emwave-c: UI Layout Helpers
// =============================================================================

#ifndef EMWAVE_UI_LAYOUT_H
#define EMWAVE_UI_LAYOUT_H

#define RENDER_DEFAULT_SCALE 2
#define RENDER_DEFAULT_MENU_BAR 28
#define RENDER_DEFAULT_TIMELINE_HEIGHT 120
#define RENDER_DEFAULT_LEFT_PANEL 160
#define RENDER_DEFAULT_RIGHT_PANEL 240

#define RENDER_MIN_PANEL_WIDTH 96
#define RENDER_MAX_PANEL_WIDTH 512
#define RENDER_MIN_TIMELINE_HEIGHT 72
#define RENDER_MAX_TIMELINE_HEIGHT 300

/* Backwards-compat aliases for older codepaths */
#define RENDER_DEFAULT_UI_HEIGHT RENDER_DEFAULT_TIMELINE_HEIGHT
#define RENDER_DEFAULT_SIDE_PANEL RENDER_DEFAULT_RIGHT_PANEL

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
    IntRect viewport;
    IntRect toolbox_panel;
    IntRect properties_panel;
    IntRect timeline_panel;
    IntRect menu_bar;
    IntRect properties_scope;
    IntRect timeline_scope;
    IntRect colorbar;
    IntRect block_outline;
} RenderLayout;

void render_layout_compute(RenderLayout* layout,
                           int nx,
                           int ny,
                           int scale,
                           int left_panel_width,
                           int right_panel_width,
                           int timeline_height,
                           int menu_bar_height);

#endif /* EMWAVE_UI_LAYOUT_H */
