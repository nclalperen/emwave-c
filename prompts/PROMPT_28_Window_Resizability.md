# Prompt #28: Window Resizability & Base Infrastructure

**Phase:** 2.75A (Essential Modernization - Week 1)
**Effort:** ~2 hours
**Priority:** CRITICAL - Foundation for all viewport improvements
**Status:** Ready for implementation
**Dependencies:** None

---

## Objective

Make the application window resizable and establish the infrastructure for dynamic viewport zoom/pan. This is the foundation for Blender-like viewport interaction.

**Strategic Value:**
- **User expectation:** Modern apps must be resizable
- **Multi-monitor support:** Essential for professional workflows
- **Foundation:** Enables zoom/pan in subsequent prompts
- **Quick win:** Simple change, huge impact

---

## Background

**Current State:**
- Window created with `SDL_WINDOW_SHOWN` only ([src/ui/ui_render.c:422-423](../src/ui/ui_render.c#L422))
- Fixed minimum size: 1920×1080
- ImGui layout calculated with fixed panel widths
- No resize event handling

**Target State:**
- Window resizable with `SDL_WINDOW_RESIZABLE` flag
- Handle `SDL_WINDOWEVENT_RESIZED` events
- Dynamic ImGui layout reflow
- Maintain minimum size (e.g., 1280×720)

---

## Requirements

### 1. Enable Window Resizability

**File:** `src/ui/ui_render.c`
**Location:** Lines 422-423

**Current Code:**
```c
ctx->window = SDL_CreateWindow(title,
                               SDL_WINDOWPOS_CENTERED,
                               SDL_WINDOWPOS_CENTERED,
                               width, height,
                               SDL_WINDOW_SHOWN);
```

**Target Code:**
```c
ctx->window = SDL_CreateWindow(title,
                               SDL_WINDOWPOS_CENTERED,
                               SDL_WINDOWPOS_CENTERED,
                               width, height,
                               SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
```

**Additional Features (Optional but Recommended):**
```c
// Set minimum window size to prevent unusably small windows
SDL_SetWindowMinimumSize(ctx->window, 1280, 720);

// Optionally set maximum size (e.g., for 4K displays)
// SDL_SetWindowMaximumSize(ctx->window, 3840, 2160);
```

---

### 2. Handle Window Resize Events

**File:** `src/app/main_imgui.cpp`
**Location:** Event loop (approximately lines 2150-2800)

**Implementation:**

```cpp
// In the SDL event loop
if (event.type == SDL_WINDOWEVENT) {
    if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
        int new_width = event.window.data1;
        int new_height = event.window.data2;

        // Update internal tracking
        width = new_width;
        height = new_height;

        // Optional: Log resize event
        ui_log_add(&app, "Window resized to %dx%d", new_width, new_height);

        // Note: ImGui automatically handles viewport size changes
        // SDL renderer automatically handles texture resizing
        // No manual framebuffer management needed
    }
}
```

**Event Types to Handle:**
- `SDL_WINDOWEVENT_RESIZED` - User resized window
- `SDL_WINDOWEVENT_SIZE_CHANGED` - Window size changed (programmatically or by user)
- `SDL_WINDOWEVENT_MAXIMIZED` - Window maximized
- `SDL_WINDOWEVENT_RESTORED` - Window restored from minimize/maximize

---

### 3. Dynamic Layout Reflow

**Current Behavior:**
- Panel widths calculated once at startup (lines ~2856-2890)
- Fixed column layout: `left_w`, `center_w`, `right_w`

**Target Behavior:**
- Recalculate panel widths every frame (already happens in current code!)
- ImGui naturally reflows with window size changes
- Viewport center area automatically adjusts

**Validation:**
Check that existing layout code (lines 2856-2890) recalculates dimensions:
```cpp
const ImGuiViewport* viewport = ImGui::GetMainViewport();
ImGui::SetNextWindowPos(viewport->Pos);
ImGui::SetNextWindowSize(viewport->Size);
// ... layout calculations happen every frame
```

**Expected:** This code already runs per-frame, so resizing should work automatically once SDL_WINDOW_RESIZABLE is enabled!

---

### 4. Add Viewport State to AppState

**File:** `src/app/main_imgui.cpp`
**Location:** AppState struct definition (lines ~35-102)

**Add fields for future zoom/pan support:**
```cpp
struct AppState {
    // ... existing fields ...

    // Viewport transformation (for Prompts #29, #30)
    float viewport_zoom;        // Zoom level (1.0 = normal, 2.0 = 2x, 0.5 = half)
    float viewport_pan_x;       // Pan offset in pixels (for future use)
    float viewport_pan_y;
    bool viewport_panning;      // Middle mouse pan active flag
    ImVec2 pan_start_mouse;     // Mouse position at pan start
    ImVec2 pan_start_offset;    // Offset at pan start

    // Window state
    int window_width;           // Current window width
    int window_height;          // Current window height
};
```

**Initialize in main():**
```cpp
app.viewport_zoom = 1.0f;
app.viewport_pan_x = 0.0f;
app.viewport_pan_y = 0.0f;
app.viewport_panning = false;
app.pan_start_mouse = ImVec2(0.0f, 0.0f);
app.pan_start_offset = ImVec2(0.0f, 0.0f);
app.window_width = width;
app.window_height = height;
```

---

### 5. Update Window Size Tracking

**In resize event handler:**
```cpp
if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
    app.window_width = event.window.data1;
    app.window_height = event.window.data2;
}
```

---

## Files to Modify

### 1. `src/ui/ui_render.c`
**Function:** `render_init()`
**Lines:** ~422-423
**Change:** Add `SDL_WINDOW_RESIZABLE` flag
**Change:** Add `SDL_SetWindowMinimumSize()` call

### 2. `src/app/main_imgui.cpp`
**Changes:**
1. **AppState struct** (lines ~35-102)
   - Add viewport zoom/pan fields
2. **main() initialization** (lines ~2050-2106)
   - Initialize new AppState fields
3. **Event loop** (lines ~2150-2800)
   - Add `SDL_WINDOWEVENT_RESIZED` handler
4. **Layout code** (lines ~2856-2890)
   - Verify dynamic recalculation (should already work!)

---

## Testing Checklist

### Basic Functionality:
- [ ] Launch application
- [ ] Grab window edge/corner with mouse
- [ ] Resize window smaller - UI reflows smoothly
- [ ] Resize window larger - UI reflows smoothly
- [ ] Panels maintain relative proportions
- [ ] Viewport center area resizes correctly
- [ ] No visual glitches or rendering artifacts

### Edge Cases:
- [ ] Resize to minimum size (1280×720) - window stops shrinking
- [ ] Try to resize smaller than minimum - prevented by SDL
- [ ] Maximize window - fills screen correctly
- [ ] Restore from maximize - returns to previous size
- [ ] Drag window between monitors (different DPI) - scaling correct
- [ ] Resize while simulation running - no crashes or freezes

### Layout Verification:
- [ ] Left column (Scene/Sources/Blocks) - resizes proportionally
- [ ] Center viewport - maintains aspect ratio of field
- [ ] Right column (Run Controls) - resizes proportionally
- [ ] Bottom panels (Scope/FFT) - adjusts to new width
- [ ] All text remains readable (no cutoff)
- [ ] Buttons remain clickable (no overlap)

### Cross-Platform:
- [ ] Windows: Resize from edges and corners works
- [ ] Windows: Maximize/restore works
- [ ] Linux: Resize from edges and corners works
- [ ] Linux: Maximize/restore works
- [ ] macOS: Resize from bottom-right corner works
- [ ] macOS: Green maximize button works

### Integration:
- [ ] Simulation continues running during resize
- [ ] Mouse interactions still work after resize
- [ ] Keyboard shortcuts still work after resize
- [ ] Material painting still works after resize
- [ ] Source placement still works after resize

---

## Success Criteria

**Phase 2.75A Step 1 complete when:**
1. ✅ Window can be resized by dragging edges/corners
2. ✅ Minimum size enforced (1280×720)
3. ✅ ImGui layout reflows smoothly
4. ✅ No crashes or visual glitches during resize
5. ✅ AppState prepared for zoom/pan (fields added)
6. ✅ Window size tracked in AppState
7. ✅ Resize events logged (optional but helpful)

---

## Implementation Notes

### Performance Considerations:
- ImGui layout calculation is cheap (<1ms per frame)
- SDL automatically handles framebuffer resizing
- No need to recreate textures or shaders
- Existing code already recalculates layout per-frame

### Platform-Specific Notes:

**Windows:**
- Resize from any edge or corner
- Maximize button in title bar
- Snap to screen edges (Win+Left/Right)

**Linux (X11/Wayland):**
- Window manager controls resize behavior
- May need compositor for smooth resize
- Test on GNOME, KDE, XFCE

**macOS:**
- Resize from bottom-right corner only (by default)
- Green maximize button (full screen toggle)
- Retina displays: Handle DPI scaling

### Debugging Tips:
- Print window size in resize handler to verify events fire
- Check `ImGui::GetMainViewport()->Size` matches SDL window size
- Use `ImGui::ShowMetricsWindow()` for layout debugging

---

## Known Issues & Workarounds

### Issue 1: Aspect Ratio
**Problem:** Field viewport may stretch if window aspect ratio changes drastically
**Solution:** Already handled! Existing code centers field in viewport (lines 2956-2983)

### Issue 2: High DPI (4K/Retina)
**Problem:** UI may appear too small on high-DPI displays
**Solution:** SDL2 handles DPI automatically on Windows/macOS. Linux may need `SDL_HINT_VIDEO_X11_FORCE_EGL` or `SDL_HINT_VIDEO_WAYLAND_ALLOW_LIBDECOR`.

### Issue 3: Minimum Size Too Small
**Problem:** 1280×720 might still be cramped on some layouts
**Solution:** Adjust `SDL_SetWindowMinimumSize()` to 1600×900 if needed (test with users)

---

## Documentation Updates

After implementation, add to user guide:

**Section: "Window Management"**

Content:
1. **Resizing the Window:**
   - Drag window edges or corners
   - Maximize button (platform-specific)
   - Minimum size: 1280×720
   - UI automatically reflows

2. **Multi-Monitor Setup:**
   - Drag window to secondary monitor
   - DPI scaling handled automatically
   - Panels can be undocked (in future prompts)

3. **Keyboard Shortcuts:**
   - Alt+Enter: Toggle fullscreen (future)
   - F11: Toggle fullscreen (future)

---

## Expected Outcome

**Before:**
```
❌ Fixed 1920×1080 window
❌ Can't resize for different monitors
❌ Wasted space on smaller screens
❌ Cramped on laptop displays
```

**After:**
```
✅ Resize to any size (minimum 1280×720)
✅ Maximize to full screen
✅ UI reflows smoothly
✅ Works on any monitor size
✅ Foundation for zoom/pan (Prompts #29-30)
```

---

## Next Steps

**After this prompt:**
1. ✅ Window is resizable
2. ✅ AppState has zoom/pan infrastructure
3. ➡️ **Prompt #29:** Implement mouse wheel zoom
4. ➡️ **Prompt #30:** Implement middle mouse pan
5. ➡️ **Prompt #31:** Add viewport toolbar & zoom UI

**Testing before moving on:**
- Verify resize works on target platforms (Windows/Linux/macOS)
- Confirm no performance degradation
- Check that existing features still work

---

**Ready for implementation!** 🚀

Hand this prompt to your implementation agent and proceed with Phase 2.75A Step 1.

---

*Created: 2025-11-20*
*Status: Ready for Implementation*
*Priority: CRITICAL - Foundation*
