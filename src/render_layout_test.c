// =============================================================================
// emwave-c: Render Layout Regression Test
// =============================================================================

#include "ui_layout.h"

#include <stdio.h>

static int validate_layout(int nx, int ny) {
    RenderLayout layout;
    render_layout_compute(&layout, nx, ny,
                          RENDER_DEFAULT_SCALE,
                          RENDER_DEFAULT_SIDE_PANEL,
                          RENDER_DEFAULT_UI_HEIGHT);

    int ok = 1;
    int cb_right = layout.colorbar.x + layout.colorbar.w;
    int cb_bottom = layout.colorbar.y + layout.colorbar.h;
    if (layout.colorbar.x < layout.canvas_w ||
        cb_right > layout.window_w ||
        cb_bottom > layout.window_h) {
        printf("Colorbar out of bounds for %dx%d (rect %d,%d,%d,%d)\n",
               nx, ny,
               layout.colorbar.x, layout.colorbar.y,
               layout.colorbar.w, layout.colorbar.h);
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

    if (layout.window_w != layout.canvas_w + RENDER_DEFAULT_SIDE_PANEL) {
        printf("Window width mismatch for %dx%d (got %d, expected %d)\n",
               nx, ny, layout.window_w,
               layout.canvas_w + RENDER_DEFAULT_SIDE_PANEL);
        ok = 0;
    }
    if (layout.window_h != layout.canvas_h + RENDER_DEFAULT_UI_HEIGHT) {
        printf("Window height mismatch for %dx%d (got %d, expected %d)\n",
               nx, ny, layout.window_h,
               layout.canvas_h + RENDER_DEFAULT_UI_HEIGHT);
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
