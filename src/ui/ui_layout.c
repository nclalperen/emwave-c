// =============================================================================
// emwave-c: UI Layout Helpers Implementation
// =============================================================================

#include "ui_layout.h"
#include "util.h"

static int clamp_panel(int v) {
    return util_clamp_int(v, RENDER_MIN_PANEL_WIDTH, RENDER_MAX_PANEL_WIDTH);
}

static int clamp_timeline(int v) {
    return util_clamp_int(v, RENDER_MIN_TIMELINE_HEIGHT, RENDER_MAX_TIMELINE_HEIGHT);
}

void render_layout_compute(RenderLayout* layout,
                           int nx,
                           int ny,
                           int scale,
                           int left_panel_width,
                           int right_panel_width,
                           int timeline_height,
                           int menu_bar_height) {
    if (!layout) return;

    if (nx < 0) nx = 0;
    if (ny < 0) ny = 0;
    if (scale < 1) scale = 1;
    left_panel_width = clamp_panel(left_panel_width);
    right_panel_width = clamp_panel(right_panel_width);
    timeline_height = clamp_timeline(timeline_height);
    if (menu_bar_height < 0) menu_bar_height = 0;

    layout->canvas_w = nx * scale;
    layout->canvas_h = ny * scale;

    int viewport_w = layout->canvas_w;
    int viewport_h = layout->canvas_h;
    int total_width = left_panel_width + viewport_w + right_panel_width;
    int total_height = menu_bar_height + viewport_h + timeline_height;

    layout->viewport.x = left_panel_width;
    layout->viewport.y = menu_bar_height;
    layout->viewport.w = viewport_w;
    layout->viewport.h = viewport_h;

    layout->toolbox_panel.x = 0;
    layout->toolbox_panel.y = menu_bar_height;
    layout->toolbox_panel.w = left_panel_width;
    layout->toolbox_panel.h = viewport_h;

    layout->properties_panel.x = left_panel_width + viewport_w;
    layout->properties_panel.y = menu_bar_height;
    layout->properties_panel.w = right_panel_width;
    layout->properties_panel.h = viewport_h;

    layout->timeline_panel.x = left_panel_width;
    layout->timeline_panel.y = menu_bar_height + viewport_h;
    layout->timeline_panel.w = viewport_w;
    layout->timeline_panel.h = timeline_height;

    layout->menu_bar.x = 0;
    layout->menu_bar.y = 0;
    layout->menu_bar.w = total_width;
    layout->menu_bar.h = menu_bar_height;

    layout->window_w = total_width;
    layout->window_h = total_height;

    const int bar_margin = 18;
    const int bar_width = 18;
    layout->colorbar.w = (right_panel_width > bar_width + 8) ? bar_width : right_panel_width;
    layout->colorbar.x = layout->properties_panel.x + layout->properties_panel.w - layout->colorbar.w - bar_margin;
    layout->colorbar.y = layout->properties_panel.y + bar_margin;
    layout->colorbar.h = layout->properties_panel.h - (bar_margin * 2);
    layout->colorbar.h = util_clamp_int(layout->colorbar.h, 0, layout->properties_panel.h);
    if (layout->colorbar.x < layout->properties_panel.x + 4) {
        layout->colorbar.x = layout->properties_panel.x + 4;
    }
    int colorbar_right_limit = layout->properties_panel.x + layout->properties_panel.w - 4;
    if (layout->colorbar.x + layout->colorbar.w > colorbar_right_limit) {
        layout->colorbar.x = colorbar_right_limit - layout->colorbar.w;
    }
    if (layout->colorbar.y < layout->properties_panel.y + 4) {
        layout->colorbar.y = layout->properties_panel.y + 4;
    }
    int colorbar_bottom_limit = layout->properties_panel.y + layout->properties_panel.h - 4;
    if (layout->colorbar.y + layout->colorbar.h > colorbar_bottom_limit) {
        layout->colorbar.h = util_clamp_int(colorbar_bottom_limit - layout->colorbar.y, 0, layout->properties_panel.h);
    }

    int max_scope_height = util_clamp_int(layout->properties_panel.h - 32, 0, layout->properties_panel.h);
    if (max_scope_height < 0) {
        max_scope_height = 0;
    }
    int desired_scope_height = layout->properties_panel.h / 3;
    int min_scope_height = 72;
    if (max_scope_height < min_scope_height) {
        min_scope_height = max_scope_height;
    }
    if (min_scope_height < 0) min_scope_height = 0;
    int properties_scope_height = util_clamp_int(desired_scope_height, min_scope_height, max_scope_height);
    layout->properties_scope.x = layout->properties_panel.x + 16;
    layout->properties_scope.w = util_clamp_int(layout->properties_panel.w - 32, 0, layout->properties_panel.w);
    layout->properties_scope.h = properties_scope_height;
    layout->properties_scope.y =
        layout->properties_panel.y + layout->properties_panel.h - properties_scope_height - 16;
    if (layout->properties_scope.y < layout->properties_panel.y + 16) {
        layout->properties_scope.y = layout->properties_panel.y + 16;
    }

    layout->timeline_scope.x = layout->viewport.x + 16;
    layout->timeline_scope.y = layout->timeline_panel.y + 16;
    layout->timeline_scope.w = util_clamp_int(layout->viewport.w - 32, 0, layout->viewport.w);
    const int timeline_controls_reserved = 72;
    int timeline_scope_height = layout->timeline_panel.h - (timeline_controls_reserved + 24);
    layout->timeline_scope.h = util_clamp_int(timeline_scope_height, 0, layout->timeline_panel.h);

    int half_block_w = nx / 10;
    int half_block_h = ny / 20;
    int center_x = nx / 2;
    int center_y = ny / 2;

    int bx0 = util_clamp_int(center_x - half_block_w, 0, (nx > 0) ? nx - 1 : 0);
    int bx1 = util_clamp_int(center_x + half_block_w, bx0 + 1, nx);
    if (bx1 <= bx0) bx1 = util_clamp_int(bx0 + 1, bx0 + 1, (nx > 0) ? nx : bx0 + 1);

    int by0 = util_clamp_int(center_y - half_block_h, 0, (ny > 0) ? ny - 1 : 0);
    int by1 = util_clamp_int(center_y + half_block_h, by0 + 1, ny);
    if (by1 <= by0) by1 = util_clamp_int(by0 + 1, by0 + 1, (ny > 0) ? ny : by0 + 1);

    layout->block_outline.x = bx0 * scale;
    layout->block_outline.y = by0 * scale;
    layout->block_outline.w = (bx1 - bx0) * scale;
    layout->block_outline.h = (by1 - by0) * scale;

    if (layout->block_outline.x + layout->block_outline.w > layout->canvas_w) {
        int excess = layout->block_outline.x + layout->block_outline.w - layout->canvas_w;
        layout->block_outline.w = util_clamp_int(layout->block_outline.w - excess, 0, layout->block_outline.w);
    }
    if (layout->block_outline.y + layout->block_outline.h > layout->canvas_h) {
        int excess = layout->block_outline.y + layout->block_outline.h - layout->canvas_h;
        layout->block_outline.h = util_clamp_int(layout->block_outline.h - excess, 0, layout->block_outline.h);
    }
}
