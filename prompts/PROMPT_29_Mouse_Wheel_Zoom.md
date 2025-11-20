# Prompt #29: Mouse Wheel Zoom (Zoom-to-Cursor)

**Phase:** 2.75A (Essential Modernization - Week 1)
**Effort:** ~4 hours
**Priority:** HIGH - Core viewport interaction
**Status:** Ready for implementation
**Dependencies:** Prompt #28 (Window Resizability)

---

## Objective

Implement smooth mouse wheel zoom with **zoom-to-cursor** behavior (like Blender). This is the most requested feature for viewport interaction and dramatically improves usability for inspecting fine details and viewing large domains.

**Strategic Value:**
- **Blender-like UX:** Industry-standard viewport interaction
- **Inspection:** Zoom in to see single-cell details (32× magnification)
- **Overview:** Zoom out to see entire large domains (0.5× scale)
- **Precision:** Zoom to cursor keeps mouse position steady (not center)

---

## Background

**Current State:**
- Fixed `scale = 2` (1 simulation cell = 2 screen pixels)
- No zoom controls
- Field always rendered at 2× scale
- Viewport centers field if smaller than available space

**Target State:**
- `scale` variable dynamically adjusted (0.5× to 32×)
- Mouse wheel scrolls zoom in/out
- **Zoom to cursor:** Mouse position stays at same simulation coordinate
- Smooth, responsive interaction

---

## Technical Approach

### Zoom-to-Cursor Algorithm

**Problem:** When zooming, we want the point under the mouse cursor to stay fixed.

**Solution:**
1. Get mouse position in viewport coordinates
2. Convert mouse → simulation grid coordinates (using current zoom)
3. Change zoom level
4. Adjust pan offset so same grid point is under mouse
5. Re-render field with new zoom

**Math:**
```
Before zoom:
  mouse_x, mouse_y = mouse position in viewport
  grid_x = (mouse_x - viewport_offset_x) / old_scale
  grid_y = (mouse_y - viewport_offset_y) / old_scale

After zoom:
  new_scale = old_scale * zoom_factor
  viewport_offset_x = mouse_x - grid_x * new_scale
  viewport_offset_y = mouse_y - grid_y * new_scale
```

---

## Requirements

### 1. Handle Mouse Wheel Events

**File:** `src/app/main_imgui.cpp`
**Location:** Event loop (lines ~2150-2800)

**Implementation:**

```cpp
// In SDL event loop, add mouse wheel handler
if (event.type == SDL_MOUSEWHEEL) {
    // Check if mouse is over viewport
    int mx, my;
    SDL_GetMouseState(&mx, &my);

    if (app.viewport_valid) {
        float local_x = (float)mx - app.viewport_pos.x;
        float local_y = (float)my - app.viewport_pos.y;

        // Check if mouse is inside viewport bounds
        if (local_x >= 0.0f && local_x < app.viewport_size.x &&
            local_y >= 0.0f && local_y < app.viewport_size.y) {

            // Zoom factor: 10% per wheel tick
            float zoom_delta = event.wheel.y * 0.1f;  // +1 or -1 typically

            // Calculate grid coordinates under mouse (before zoom)
            float field_px_w = (float)(sim->nx * scale);
            float field_px_h = (float)(sim->ny * scale);

            // Mouse position relative to field origin
            float field_local_x = local_x;
            float field_local_y = local_y;

            // Grid coordinates (0 to nx, 0 to ny)
            float grid_x = (field_local_x - app.viewport_pan_x) / (float)scale;
            float grid_y = (field_local_y - app.viewport_pan_y) / (float)scale;

            // Apply zoom
            float old_scale = (float)scale;
            float new_scale = old_scale * (1.0f + zoom_delta);

            // Clamp zoom range: 0.5× to 32×
            if (new_scale < 0.5f) new_scale = 0.5f;
            if (new_scale > 32.0f) new_scale = 32.0f;

            // Update scale (cast to int for SDL rendering)
            scale = (int)(new_scale + 0.5f);  // Round to nearest int
            if (scale < 1) scale = 1;

            // Adjust pan to keep grid_x, grid_y under mouse
            app.viewport_pan_x = field_local_x - grid_x * (float)scale;
            app.viewport_pan_y = field_local_y - grid_y * (float)scale;

            // Update render scale
            render->scale = scale;

            // Log zoom level (optional)
            ui_log_add(&app, "Zoom: %.1f× (scale=%d)", new_scale, scale);
        }
    }
}
```

**Key Points:**
- `event.wheel.y > 0` = scroll up = zoom in
- `event.wheel.y < 0` = scroll down = zoom out
- `zoom_delta = 0.1` means 10% zoom per tick
- Only zoom when mouse over viewport (don't zoom in panels)

---

### 2. Add Keyboard Zoom Shortcuts

**Additional zoom controls for accessibility:**

```cpp
// In keyboard event handler
if (event.type == SDL_KEYDOWN) {
    SDL_Keymod mod = SDL_GetModState();

    // Ctrl++ : Zoom in (centered)
    if ((mod & KMOD_CTRL) && event.key.keysym.sym == SDLK_EQUALS) {
        float new_scale = (float)scale * 1.2f;  // 20% zoom in
        if (new_scale > 32.0f) new_scale = 32.0f;
        scale = (int)(new_scale + 0.5f);
        render->scale = scale;
        ui_log_add(&app, "Zoom in: %dx", scale);
    }

    // Ctrl+- : Zoom out (centered)
    if ((mod & KMOD_CTRL) && event.key.keysym.sym == SDLK_MINUS) {
        float new_scale = (float)scale * 0.8f;  // 20% zoom out
        if (new_scale < 0.5f) new_scale = 0.5f;
        scale = (int)(new_scale + 0.5f);
        if (scale < 1) scale = 1;
        render->scale = scale;
        ui_log_add(&app, "Zoom out: %dx", scale);
    }

    // Ctrl+0 or Home : Reset zoom to 1:1 (and center view)
    if (((mod & KMOD_CTRL) && event.key.keysym.sym == SDLK_0) ||
        event.key.keysym.sym == SDLK_HOME) {
        scale = 2;  // Reset to default 2x
        render->scale = scale;
        app.viewport_pan_x = 0.0f;
        app.viewport_pan_y = 0.0f;
        app.viewport_zoom = 1.0f;
        ui_log_add(&app, "Reset zoom to %dx", scale);
    }
}
```

---

### 3. Update Rendering with Pan Offset

**File:** `src/app/main_imgui.cpp`
**Location:** SDL rendering calls (lines ~1400-1500)

**Current rendering:**
```cpp
// Current: render at (i * scale, j * scale)
SDL_Rect rect = {i * scale, j * scale, scale, scale};
SDL_RenderFillRect(renderer, &rect);
```

**Updated rendering with pan:**
```cpp
// Add pan offset to rendering
int draw_x = (int)(i * scale + app.viewport_pan_x);
int draw_y = (int)(j * scale + app.viewport_pan_y);
SDL_Rect rect = {draw_x, draw_y, scale, scale};
SDL_RenderFillRect(renderer, &rect);
```

**Apply to all rendering functions:**
1. `render_field_to_texture()` or equivalent
2. Grid overlay rendering
3. Material overlay rendering
4. Source markers
5. Any other viewport overlays

---

### 4. Handle Float→Int Scale Conversion

**Challenge:** SDL requires integer pixel coordinates, but smooth zoom needs float precision.

**Solution:** Use two variables:
```cpp
// In AppState (from Prompt #28)
float viewport_zoom;    // Float: 0.5 to 32.0 (smooth)
int scale;              // Int: Quantized for SDL rendering

// During zoom update:
app.viewport_zoom *= (1.0f + zoom_delta);
app.viewport_zoom = clamp(app.viewport_zoom, 0.5f, 32.0f);
scale = (int)(app.viewport_zoom + 0.5f);  // Round
if (scale < 1) scale = 1;
render->scale = scale;
```

**Benefit:** `viewport_zoom` accumulates smooth changes, `scale` is quantized for rendering.

---

### 5. Clamp Pan Offset to Reasonable Bounds

**Prevent viewport from panning infinitely off-screen:**

```cpp
// After adjusting pan offset
// Limit pan so at least 10% of field is visible

float field_w = (float)(sim->nx * scale);
float field_h = (float)(sim->ny * scale);
float viewport_w = app.viewport_size.x;
float viewport_h = app.viewport_size.y;

// Allow panning, but keep some field visible
float min_pan_x = -(field_w * 0.9f);  // Can pan 90% off-screen
float max_pan_x = viewport_w * 0.9f;
float min_pan_y = -(field_h * 0.9f);
float max_pan_y = viewport_h * 0.9f;

app.viewport_pan_x = clamp(app.viewport_pan_x, min_pan_x, max_pan_x);
app.viewport_pan_y = clamp(app.viewport_pan_y, min_pan_y, max_pan_y);
```

---

## Files to Modify

### 1. `src/app/main_imgui.cpp`
**Changes:**
1. **Event loop** (lines ~2150-2800)
   - Add `SDL_MOUSEWHEEL` handler
   - Add keyboard shortcuts (Ctrl+±, Ctrl+0, Home)
2. **Rendering functions** (lines ~1400-1500, ~1446-1500)
   - Add `app.viewport_pan_x/y` to all draw positions
3. **AppState usage**
   - Update `app.viewport_zoom` float value
   - Update `scale` and `render->scale` integer values
4. **Main loop**
   - Apply pan offsets to all viewport rendering

---

## Testing Checklist

### Basic Zoom Functionality:
- [ ] Launch application with default 2× zoom
- [ ] Scroll mouse wheel up (over viewport) - field zooms in
- [ ] Scroll mouse wheel down - field zooms out
- [ ] Zoom in to 32× (maximum) - wheel stops zooming in
- [ ] Zoom out to 0.5× (minimum) - wheel stops zooming out
- [ ] Zoom with mouse in different viewport positions - zooms to cursor
- [ ] Zoom with mouse outside viewport - no zoom (panels unaffected)

### Zoom-to-Cursor Behavior:
- [ ] Place mouse over specific cell (e.g., material block)
- [ ] Scroll wheel up - that cell stays under mouse (not center)
- [ ] Scroll wheel down - cell still under mouse
- [ ] Move mouse to different location, zoom again - new point fixed
- [ ] Compare to center-zoom (wrong behavior) - should NOT center-zoom

### Keyboard Shortcuts:
- [ ] Press Ctrl++ - zooms in (centered)
- [ ] Press Ctrl+- - zooms out (centered)
- [ ] Press Ctrl+0 - resets to 2× zoom, centered view
- [ ] Press Home - resets zoom and pan
- [ ] Shortcuts work regardless of mouse position

### Zoom Levels:
- [ ] 0.5× zoom - see large domain (1000×1000 cells visible)
- [ ] 1× zoom - 1:1 pixel mapping
- [ ] 2× zoom - default (current behavior)
- [ ] 4× zoom - see finer details
- [ ] 8× zoom - inspect individual cells
- [ ] 16× zoom - very close inspection
- [ ] 32× zoom - maximum magnification (single cell fills screen)

### Rendering Quality:
- [ ] No visual glitches at any zoom level
- [ ] Grid overlay scales correctly
- [ ] Material colors render correctly
- [ ] Source markers scale appropriately
- [ ] Text overlays remain readable
- [ ] No jitter or flickering during zoom

### Performance:
- [ ] Smooth zoom (no lag or stutter)
- [ ] Simulation continues running while zooming
- [ ] FPS remains stable (60 FPS target)
- [ ] Zoom works on large simulations (1000×1000+)
- [ ] Memory usage stable (no leaks)

### Integration:
- [ ] Mouse interactions still work (source placement, painting)
- [ ] Keyboard shortcuts still work (pause, step, etc.)
- [ ] Panels still functional after zooming
- [ ] Config save/load preserves zoom level (optional)
- [ ] Log shows zoom level changes

---

## Success Criteria

**Phase 2.75A Step 2 complete when:**
1. ✅ Mouse wheel zooms in/out over viewport
2. ✅ Zoom-to-cursor behavior works (point under mouse stays fixed)
3. ✅ Zoom range: 0.5× to 32× enforced
4. ✅ Keyboard shortcuts (Ctrl+±, Ctrl+0, Home) work
5. ✅ Smooth, responsive zoom (no lag)
6. ✅ Rendering updated with pan offsets
7. ✅ No visual glitches or performance issues

---

## Implementation Notes

### Zoom Range Rationale:
- **0.5× (half):** See 4× more area (useful for 1000×1000 domains)
- **1× (1:1):** One simulation cell = one screen pixel
- **2× (default):** Current behavior (comfortable viewing)
- **4×-8×:** Standard working range
- **16×-32×:** Extreme close-up (single cell inspection)

### Performance Considerations:
- Zoom only changes rendering scale, not simulation grid
- SDL hardware acceleration handles scaling efficiently
- Clipping outside viewport prevents overdraw
- Integer scale avoids sub-pixel rendering (fast)

### Alternative: Sub-Pixel Rendering (Future)
- Use `SDL_RenderSetScale()` for smooth float zoom
- Requires OpenGL/DirectX backend
- More complex but smoother visuals
- Consider for Prompt #35 (Polish)

---

## Known Issues & Workarounds

### Issue 1: Integer Scale Quantization
**Problem:** Scale jumps between 1×, 2×, 3× (not smooth)
**Impact:** Low (zoom accumulation in `viewport_zoom` helps)
**Workaround:** Use `viewport_zoom` float for accumulation, quantize to `scale` int
**Future Fix:** Migrate to ImGui::Image with float scaling (Prompt #36+)

### Issue 2: High Zoom Aliasing
**Problem:** At 16×-32×, individual pixels very large (blocky)
**Impact:** Low (expected for pixel art / simulation data)
**Workaround:** None needed (feature, not bug)
**Future Enhancement:** Bilinear interpolation (optional polish)

### Issue 3: Pan Drift at Extreme Zoom
**Problem:** Float→Int conversion may cause 1-pixel drift
**Impact:** Very low (< 1 pixel error)
**Workaround:** Clamp pan offsets to grid alignment
**Solution:** `app.viewport_pan_x = round(app.viewport_pan_x)`

---

## Documentation Updates

After implementation, add to user guide:

**Section: "Viewport Navigation - Zoom"**

Content:
1. **Mouse Wheel Zoom:**
   - Scroll up: Zoom in (toward mouse)
   - Scroll down: Zoom out
   - Zoom range: 0.5× to 32×
   - Zoom to cursor (not center)

2. **Keyboard Zoom:**
   - Ctrl++: Zoom in (centered)
   - Ctrl+-: Zoom out (centered)
   - Ctrl+0: Reset zoom to 2×
   - Home: Reset zoom and pan

3. **Zoom Levels Guide:**
   - 0.5×-1×: Large domain overview
   - 2×-4×: Normal working range
   - 8×-16×: Detailed inspection
   - 32×: Maximum magnification

4. **Tips:**
   - Hold Ctrl for faster zoom (2× per tick)
   - Use Home to reset if lost
   - Zoom to cursor keeps reference point

---

## Expected Outcome

**Before:**
```
❌ Fixed 2× scale (no zoom)
❌ Can't inspect fine details
❌ Can't see large domains
❌ No zoom controls
```

**After:**
```
✅ Smooth mouse wheel zoom (0.5× to 32×)
✅ Zoom to cursor (Blender-like)
✅ Keyboard shortcuts (Ctrl+±, Home)
✅ Inspect single cells at 32×
✅ See entire large domains at 0.5×
✅ Professional viewport interaction
```

**Visual Example:**
```
At 2× (default):     At 8× (zoomed in):    At 0.5× (zoomed out):
┌─────────────┐      ┌─────────────┐        ┌─────────────┐
│ ░░▓▓░░      │      │ ░░░░░░      │        │░▓░▓░▓░▓░▓░▓░│
│ ░░▓▓░░      │  →   │ ░░░░░░      │   ←    │░▓░▓░▓░▓░▓░▓░│
│ ░░▓▓░░      │      │ ░░░░░░      │        │░▓░▓░▓░▓░▓░▓░│
│ ░░▓▓░░      │      │ ░░░░░░      │        │░▓░▓░▓░▓░▓░▓░│
└─────────────┘      └─────────────┘        └─────────────┘
   100×100            25×25 visible         400×400 visible
```

---

## Next Steps

**After this prompt:**
1. ✅ Mouse wheel zoom working
2. ✅ Zoom-to-cursor behavior implemented
3. ➡️ **Prompt #30:** Implement middle mouse pan
4. ➡️ **Prompt #31:** Add viewport toolbar with zoom UI

**Testing before moving on:**
- Zoom at various levels (0.5× to 32×)
- Verify zoom-to-cursor works correctly
- Check performance with large simulations
- Test keyboard shortcuts

---

**Ready for implementation!** 🚀

Hand this prompt to your implementation agent and proceed with Phase 2.75A Step 2.

---

*Created: 2025-11-20*
*Status: Ready for Implementation*
*Priority: HIGH - Core Feature*
