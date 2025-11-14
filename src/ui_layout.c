// =============================================================================
// emwave-c: UI Layout Helpers Implementation
// =============================================================================

#include "ui_layout.h"

static int clamp_int(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

void render_layout_compute(RenderLayout* layout,
                           int nx,
                           int ny,
                           int scale,
                           int side_panel_width,
                           int ui_height) {
    if (!layout) return;

    if (nx < 0) nx = 0;
    if (ny < 0) ny = 0;
    if (scale < 1) scale = 1;
    if (side_panel_width < 0) side_panel_width = 0;
    if (ui_height < 0) ui_height = 0;

    layout->canvas_w = nx * scale;
    layout->canvas_h = ny * scale;
    layout->window_w = layout->canvas_w + side_panel_width;
    layout->window_h = layout->canvas_h + ui_height;

    layout->colorbar.x = layout->canvas_w + 20;
    layout->colorbar.y = 20;
    layout->colorbar.w = 20;
    int cb_h = layout->canvas_h - 40;
    if (cb_h < 0) cb_h = layout->canvas_h;
    if (cb_h < 0) cb_h = 0;
    layout->colorbar.h = cb_h;
    if (layout->colorbar.y + layout->colorbar.h > layout->window_h) {
        layout->colorbar.h = clamp_int(layout->window_h - layout->colorbar.y, 0, layout->window_h);
    }
    if (layout->colorbar.x + layout->colorbar.w > layout->window_w) {
        layout->colorbar.w = clamp_int(layout->window_w - layout->colorbar.x, 0, layout->colorbar.w);
    }

    int half_block_w = nx / 10;
    int half_block_h = ny / 20;
    int center_x = nx / 2;
    int center_y = ny / 2;

    int bx0 = clamp_int(center_x - half_block_w, 0, (nx > 0) ? nx - 1 : 0);
    int bx1 = clamp_int(center_x + half_block_w, bx0 + 1, nx);
    if (bx1 <= bx0) bx1 = clamp_int(bx0 + 1, bx0 + 1, (nx > 0) ? nx : bx0 + 1);

    int by0 = clamp_int(center_y - half_block_h, 0, (ny > 0) ? ny - 1 : 0);
    int by1 = clamp_int(center_y + half_block_h, by0 + 1, ny);
    if (by1 <= by0) by1 = clamp_int(by0 + 1, by0 + 1, (ny > 0) ? ny : by0 + 1);

    layout->block_outline.x = bx0 * scale;
    layout->block_outline.y = by0 * scale;
    layout->block_outline.w = (bx1 - bx0) * scale;
    layout->block_outline.h = (by1 - by0) * scale;

    if (layout->block_outline.x + layout->block_outline.w > layout->canvas_w) {
        int excess = layout->block_outline.x + layout->block_outline.w - layout->canvas_w;
        layout->block_outline.w = clamp_int(layout->block_outline.w - excess, 0, layout->block_outline.w);
    }
    if (layout->block_outline.y + layout->block_outline.h > layout->canvas_h) {
        int excess = layout->block_outline.y + layout->block_outline.h - layout->canvas_h;
        layout->block_outline.h = clamp_int(layout->block_outline.h - excess, 0, layout->block_outline.h);
    }
}
