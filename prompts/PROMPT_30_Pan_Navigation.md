# Prompt #30: Middle Mouse Pan & Keyboard Navigation

**Phase:** 2.75A (Essential Modernization - Week 1)
**Effort:** ~3 hours
**Priority:** HIGH - Complete viewport interaction trinity (resize/zoom/pan)
**Status:** Ready for implementation
**Dependencies:** Prompts #28 (Resizability), #29 (Zoom)

---

## Objective

Implement viewport panning with middle mouse drag (Blender standard) and keyboard arrow keys. Combined with zoom (#29), this completes professional viewport navigation.

**Strategic Value:**
- **Blender UX:** Middle mouse drag is industry standard
- **Accessibility:** Keyboard arrows + Shift+Right for laptops
- **Large domains:** Navigate 1000×1000+ grids easily
- **Precision:** Combined with zoom, enables precise positioning

---

## Background

**Current State:**
- Field always centered in viewport
- No pan controls
- Can't navigate off-center regions
- `viewport_pan_x/y` initialized but unused

**Target State:**
- Middle mouse drag pans viewport
- Shift+Right drag as alternative (laptop users)
- Arrow keys pan in 10px increments
- Shift+Arrow pans 50px (fast)
- Pan offset persists across frames

---

## Requirements

### 1. Middle Mouse Drag Pan

**File:** `src/app/main_imgui.cpp`
**Location:** Event loop (lines ~2150-2800)

**Implementation:**

```cpp
// Handle middle mouse button down
if (event.type == SDL_MOUSEBUTTONDOWN) {
    if (event.button.button == SDL_BUTTON_MIDDLE) {
        // Check if mouse is over viewport
        int mx = event.button.x;
        int my = event.button.y;

        if (app.viewport_valid) {
            float local_x = (float)mx - app.viewport_pos.x;
            float local_y = (float)my - app.viewport_pos.y;

            if (local_x >= 0.0f && local_x < app.viewport_size.x &&
                local_y >= 0.0f && local_y < app.viewport_size.y) {
                // Start panning
                app.viewport_panning = true;
                app.pan_start_mouse = ImVec2((float)mx, (float)my);
                app.pan_start_offset = ImVec2(app.viewport_pan_x, app.viewport_pan_y);

                // Optional: Change cursor to "move" icon
                SDL_SetCursor(SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_SIZEALL));
            }
        }
    }
}

// Handle mouse motion (pan drag)
if (event.type == SDL_MOUSEMOTION) {
    if (app.viewport_panning) {
        int mx = event.motion.x;
        int my = event.motion.y;

        // Calculate delta from pan start
        float delta_x = (float)mx - app.pan_start_mouse.x;
        float delta_y = (float)my - app.pan_start_mouse.y;

        // Update pan offset
        app.viewport_pan_x = app.pan_start_offset.x + delta_x;
        app.viewport_pan_y = app.pan_start_offset.y + delta_y;

        // Optional: Clamp pan to reasonable bounds (see Section 4)
    }
}

// Handle middle mouse button up
if (event.type == SDL_MOUSEBUTTONUP) {
    if (event.button.button == SDL_BUTTON_MIDDLE) {
        if (app.viewport_panning) {
            app.viewport_panning = false;

            // Restore cursor
            SDL_SetCursor(SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_ARROW));

            ui_log_add(&app, "Pan offset: (%.0f, %.0f)",
                       app.viewport_pan_x, app.viewport_pan_y);
        }
    }
}
```

---

### 2. Alternative Pan Control (Shift+Right Drag)

**For laptop users without middle mouse button:**

```cpp
// Track Shift+Right drag as alternative pan
static bool shift_right_panning = false;

if (event.type == SDL_MOUSEBUTTONDOWN) {
    if (event.button.button == SDL_BUTTON_RIGHT) {
        SDL_Keymod mod = SDL_GetModState();
        if (mod & KMOD_SHIFT) {
            // Same as middle mouse logic
            shift_right_panning = true;
            // ... (copy middle mouse start logic)
        }
    }
}

if (event.type == SDL_MOUSEMOTION) {
    if (shift_right_panning) {
        // Same as middle mouse motion logic
        // ... (copy pan delta calculation)
    }
}

if (event.type == SDL_MOUSEBUTTONUP) {
    if (event.button.button == SDL_BUTTON_RIGHT) {
        if (shift_right_panning) {
            shift_right_panning = false;
            // ... (copy pan end logic)
        }
    }
}
```

---

### 3. Keyboard Arrow Key Pan

**Accessibility and precision control:**

```cpp
// In keyboard event handler
if (event.type == SDL_KEYDOWN) {
    SDL_Keymod mod = SDL_GetModState();
    float pan_step = (mod & KMOD_SHIFT) ? 50.0f : 10.0f;  // Shift = fast

    switch (event.key.keysym.sym) {
        case SDLK_LEFT:
            app.viewport_pan_x += pan_step;
            ui_log_add(&app, "Pan left");
            break;

        case SDLK_RIGHT:
            app.viewport_pan_x -= pan_step;
            ui_log_add(&app, "Pan right");
            break;

        case SDLK_UP:
            app.viewport_pan_y += pan_step;
            ui_log_add(&app, "Pan up");
            break;

        case SDLK_DOWN:
            app.viewport_pan_y -= pan_step;
            ui_log_add(&app, "Pan down");
            break;

        case SDLK_HOME:
            // Reset pan and zoom (already in Prompt #29)
            app.viewport_pan_x = 0.0f;
            app.viewport_pan_y = 0.0f;
            app.viewport_zoom = 1.0f;
            scale = 2;
            render->scale = 2;
            ui_log_add(&app, "Reset view");
            break;
    }

    // Optional: Clamp pan after arrow key adjustment
}
```

---

### 4. Clamp Pan Offset to Reasonable Bounds

**Prevent viewport from panning infinitely off-screen:**

```cpp
// Helper function to clamp pan offsets
static void clamp_pan_offset(AppState* app, const SimulationState* sim, int scale) {
    if (!sim || scale <= 0) return;

    float field_w = (float)(sim->nx * scale);
    float field_h = (float)(sim->ny * scale);
    float viewport_w = app->viewport_size.x;
    float viewport_h = app->viewport_size.y;

    // Allow panning until 90% of field is off-screen
    // (keeps at least 10% visible for orientation)
    float min_pan_x = -(field_w * 0.9f);
    float max_pan_x = viewport_w * 0.9f;
    float min_pan_y = -(field_h * 0.9f);
    float max_pan_y = viewport_h * 0.9f;

    app->viewport_pan_x = fmaxf(min_pan_x, fminf(max_pan_x, app->viewport_pan_x));
    app->viewport_pan_y = fmaxf(min_pan_y, fminf(max_pan_y, app->viewport_pan_y));
}

// Call after any pan operation
clamp_pan_offset(&app, sim, scale);
```

---

### 5. Update Mouse→Grid Coordinate Conversion

**Account for pan offset when converting mouse clicks to grid cells:**

**Current code (lines ~2525-2533):**
```cpp
float local_x = (float)mx - app.viewport_pos.x;
float local_y = (float)my - app.viewport_pos.y;
int ix = (int)(local_x / (float)scale);
int iy = (int)(local_y / (float)scale);
```

**Updated code with pan:**
```cpp
float local_x = (float)mx - app.viewport_pos.x;
float local_y = (float)my - app.viewport_pos.y;

// Subtract pan offset to get field-relative coordinates
float field_x = local_x - app.viewport_pan_x;
float field_y = local_y - app.viewport_pan_y;

// Convert to grid indices
int ix = (int)(field_x / (float)scale);
int iy = (int)(field_y / (float)scale);

// Bounds check
if (ix < 0 || ix >= sim->nx || iy < 0 || iy >= sim->ny) {
    // Mouse outside field bounds
    return;
}
```

**Apply this update to:**
- Source placement (lines ~2618-2633)
- Material painting (lines ~2656-2670)
- Block placement (similar areas)
- Any other mouse interactions

---

### 6. Visual Feedback for Panning

**Optional but helpful:**

```cpp
// While panning, show drag cursor and visual feedback
if (app.viewport_panning) {
    // Draw subtle overlay indicating pan mode
    ImGui::SetCursorPos(ImVec2(app.viewport_pos.x + 10, app.viewport_pos.y + 10));
    ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 0.7f),
                       "Panning... (%.0f, %.0f)",
                       app.viewport_pan_x, app.viewport_pan_y);
}
```

---

## Files to Modify

### 1. `src/app/main_imgui.cpp`
**Changes:**
1. **Event loop** (lines ~2150-2800)
   - Add middle mouse button handlers (down/motion/up)
   - Add Shift+Right drag handlers
   - Add arrow key pan handlers
2. **Mouse→Grid conversion** (multiple locations)
   - Update ~lines 2525-2533 (source placement)
   - Update ~lines 2618-2633 (material painting)
   - Update ~lines 2656-2670 (block placement)
   - Add pan offset subtraction to all conversions
3. **Helper function**
   - Add `clamp_pan_offset()` function
4. **Pan clamping**
   - Call `clamp_pan_offset()` after pan operations

---

## Testing Checklist

### Middle Mouse Pan:
- [ ] Click and hold middle mouse button over viewport
- [ ] Drag mouse - field pans with cursor
- [ ] Release middle mouse - pan stops, offset retained
- [ ] Pan in all directions (up/down/left/right)
- [ ] Pan while zoomed in - smooth, responsive
- [ ] Pan while zoomed out - works correctly
- [ ] Cursor changes to "move" icon during pan (optional)

### Shift+Right Pan (Alternative):
- [ ] Hold Shift, click and hold right mouse
- [ ] Drag mouse - field pans (same as middle mouse)
- [ ] Release - pan stops
- [ ] Works identically to middle mouse
- [ ] Useful on laptop trackpad

### Keyboard Arrow Pan:
- [ ] Press Left arrow - field pans right (10px)
- [ ] Press Right arrow - field pans left (10px)
- [ ] Press Up arrow - field pans down (10px)
- [ ] Press Down arrow - field pans up (10px)
- [ ] Hold Shift+Arrow - pans faster (50px per press)
- [ ] Press Home - resets pan and zoom

### Pan Clamping:
- [ ] Pan until field is 90% off-screen - stops at limit
- [ ] Can't pan completely off-screen
- [ ] At least 10% of field always visible
- [ ] Clamping smooth (no jitter)
- [ ] Works at all zoom levels

### Mouse Interactions with Pan:
- [ ] Click to place source - correct grid cell even with pan
- [ ] Paint material - brush follows mouse correctly
- [ ] Place material block - coordinates correct
- [ ] Hover over cell - shows correct coordinates
- [ ] All mouse interactions account for pan offset

### Integration with Zoom:
- [ ] Zoom in, then pan - both work together
- [ ] Pan, then zoom - zoom-to-cursor still works
- [ ] Pan at 0.5× zoom - entire domain navigable
- [ ] Pan at 32× zoom - fine positioning works
- [ ] Combo: Zoom in, pan to region, zoom more - smooth

### Edge Cases:
- [ ] Pan while simulation running - no lag
- [ ] Pan with window resize - correct behavior
- [ ] Pan with multiple sources/materials - rendering correct
- [ ] Pan limits at different zoom levels - appropriate
- [ ] Reset view (Home) from any pan position - works

---

## Success Criteria

**Phase 2.75A Step 3 complete when:**
1. ✅ Middle mouse drag pans viewport
2. ✅ Shift+Right drag pans viewport (alternative)
3. ✅ Arrow keys pan in 10px increments
4. ✅ Shift+Arrow pans in 50px increments
5. ✅ Home key resets pan and zoom
6. ✅ Pan clamping prevents infinite off-screen
7. ✅ Mouse interactions account for pan offset
8. ✅ Smooth, responsive panning (no lag)

---

## Implementation Notes

### Pan Direction Convention:
```
Middle mouse drag RIGHT → Field appears to move RIGHT (pan offset +X)
Arrow key LEFT → Field appears to move LEFT (pan offset -X)
```
This matches Blender behavior (viewport follows drag, not camera).

### Performance:
- Pan only updates offset values (cheap)
- Rendering already supports offset (from Prompt #29)
- No texture recreation needed
- 60 FPS maintained during pan

### Cursor Icons:
```cpp
SDL_SYSTEM_CURSOR_ARROW      // Normal
SDL_SYSTEM_CURSOR_SIZEALL    // Pan (4-way arrows)
SDL_SYSTEM_CURSOR_HAND       // Alternative
```

---

## Known Issues & Workarounds

### Issue 1: Pan Drift with Integer Rendering
**Problem:** Float offset quantized to int may cause 1px drift
**Impact:** Negligible (< 1 pixel)
**Workaround:** Round pan offsets: `app.viewport_pan_x = roundf(app.viewport_pan_x)`

### Issue 2: Trackpad Scroll vs. Middle Click
**Problem:** Laptops may not have middle mouse button
**Solution:** Shift+Right drag provided as alternative
**Alternative:** Two-finger trackpad drag (platform-dependent)

### Issue 3: Pan Limits Feel Restrictive
**Problem:** 90% off-screen limit may feel limiting
**Solution:** Adjustable - change `0.9f` to `0.95f` or `1.0f` (unlimited)
**Recommendation:** Start with 90%, gather user feedback

---

## Documentation Updates

After implementation, add to user guide:

**Section: "Viewport Navigation - Pan"**

Content:
1. **Middle Mouse Pan:**
   - Click and drag middle mouse button
   - Viewport follows cursor
   - Industry standard (Blender, Maya, etc.)

2. **Alternative Pan (Laptops):**
   - Hold Shift + drag right mouse
   - Same behavior as middle mouse
   - For trackpads without middle button

3. **Keyboard Pan:**
   - Arrow keys: Pan 10 pixels
   - Shift+Arrow: Pan 50 pixels (fast)
   - Home: Reset view (zoom + pan)

4. **Pan Tips:**
   - Combine with zoom for precise navigation
   - Pan limits keep 10% field visible
   - Use Home if lost

5. **Workflow Example:**
   ```
   1. Zoom out (Ctrl+-) to see entire domain
   2. Zoom in (wheel) on region of interest
   3. Pan (middle drag) to fine-tune position
   4. Inspect details at high zoom
   ```

---

## Expected Outcome

**Before:**
```
❌ Field always centered
❌ Can't navigate large domains
❌ Can't position viewport precisely
❌ No pan controls
```

**After:**
```
✅ Middle mouse drag pan (Blender standard)
✅ Shift+Right drag (laptop alternative)
✅ Arrow key pan (keyboard precision)
✅ Navigate 1000×1000+ domains easily
✅ Combined zoom+pan = professional UX
✅ Complete viewport navigation trinity
```

**Visual Example:**
```
Pan Left:               Pan Right:              Pan Up/Down:
┌─────────────┐        ┌─────────────┐         ┌─────────────┐
│    ░░▓▓░░   │   →    │  ░░▓▓░░     │    ↕    │ ░░▓▓░░      │
│    ░░▓▓░░   │        │  ░░▓▓░░     │         │ ░░▓▓░░      │
│    ░░▓▓░░   │        │  ░░▓▓░░     │         │             │
└─────────────┘        └─────────────┘         └─────────────┘
```

---

## Next Steps

**After this prompt:**
1. ✅ Window resizable (Prompt #28)
2. ✅ Mouse wheel zoom (Prompt #29)
3. ✅ Middle mouse pan (Prompt #30)
4. ➡️ **Prompt #31:** Add viewport toolbar with zoom/pan UI
5. ➡️ **Phase 2.75A COMPLETE!**

**Testing before moving on:**
- Test all three pan methods (middle, Shift+Right, arrows)
- Verify pan works with zoom at all levels
- Check mouse interactions account for pan
- Test on laptop trackpad

---

**Ready for implementation!** 🚀

Hand this prompt to your implementation agent and proceed with Phase 2.75A Step 3.

---

*Created: 2025-11-20*
*Status: Ready for Implementation*
*Priority: HIGH - Completes Viewport Trinity*
