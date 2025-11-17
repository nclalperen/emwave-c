// =============================================================================
// emwave-c: Render Layout Regression Test
// =============================================================================

#include "ui_layout.h"

#include <stdio.h>

static int validate_layout(int nx, int ny) {
    RenderLayout layout;
    render_layout_compute(&layout, nx, ny,
                          RENDER_DEFAULT_SCALE,
                          RENDER_DEFAULT_LEFT_PANEL,
                          RENDER_DEFAULT_RIGHT_PANEL,
                          RENDER_DEFAULT_TIMELINE_HEIGHT,
                          RENDER_DEFAULT_MENU_BAR);

    int ok = 1;
    int cb_right = layout.colorbar.x + layout.colorbar.w;
    int cb_bottom = layout.colorbar.y + layout.colorbar.h;
    int props_right = layout.properties_panel.x + layout.properties_panel.w;
    int props_bottom = layout.properties_panel.y + layout.properties_panel.h;
    if (layout.colorbar.x < layout.properties_panel.x ||
        cb_right > props_right ||
        layout.colorbar.y < layout.properties_panel.y ||
        cb_bottom > props_bottom) {
        printf("Colorbar out of bounds for %dx%d (rect %d,%d,%d,%d within props %d,%d,%d,%d)\n",
               nx, ny,
               layout.colorbar.x, layout.colorbar.y,
               layout.colorbar.w, layout.colorbar.h,
               layout.properties_panel.x, layout.properties_panel.y,
               layout.properties_panel.w, layout.properties_panel.h);
        ok = 0;
    }

    int block_right = layout.block_outline.x + layout.block_outline.w;
    int block_bottom = layout.block_outline.y + layout.block_outline.h;
    if (layout.block_outline.x < 0 || layout.block_outline.y < 0 ||
        block_right > layout.canvas_w || block_bottom > layout.canvas_h) {
        printf("Block outline out of bounds for %dx%d (rect %d,%d,%d,%d)\n",
               nx, ny,
               layout.block_outline.x, layout.block_outline.y,
               layout.block_outline.w, layout.block_outline.h);
        ok = 0;
    }

    int expected_width = layout.toolbox_panel.w + layout.viewport.w + layout.properties_panel.w;
    if (layout.window_w != expected_width) {
        printf("Window width mismatch for %dx%d (got %d, expected %d)\n",
               nx, ny, layout.window_w, expected_width);
        ok = 0;
    }

    int expected_height = layout.menu_bar.h + layout.viewport.h + layout.timeline_panel.h;
    if (layout.window_h != expected_height) {
        printf("Window height mismatch for %dx%d (got %d, expected %d)\n",
               nx, ny, layout.window_h, expected_height);
        ok = 0;
    }

    if (layout.timeline_panel.y != layout.menu_bar.h + layout.viewport.h) {
        printf("Timeline offset mismatch for %dx%d (got %d)\n",
               nx, ny, layout.timeline_panel.y);
        ok = 0;
    }

    return ok;
}

int main(void) {
    int ok = 1;
    ok &= validate_layout(128, 96);
    ok &= validate_layout(512, 384);
    if (!ok) {
        fprintf(stderr, "Render layout validation failed.\n");
        return 1;
    }
    printf("Render layout validation passed.\n");
    return 0;
}
