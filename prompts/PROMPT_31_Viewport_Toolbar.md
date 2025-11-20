# Prompt #31: Viewport Toolbar & Zoom Controls UI

**Phase:** 2.75A (Essential Modernization - Week 1 - Final)
**Effort:** ~3 hours
**Priority:** MEDIUM - UI polish for zoom/pan controls
**Status:** Ready for implementation
**Dependencies:** Prompts #28-30 (Resize/Zoom/Pan)

---

## Objective

Add a floating toolbar overlay to the viewport with visual zoom controls, reset button, and status display. Provides discoverability for new users and quick access to zoom/pan features.

---

## Implementation

### 1. Floating Toolbar (Bottom-Right Corner)

```cpp
// In viewport rendering section (after field draw, before end of viewport)
if (app.viewport_valid && sim) {
    ImGui::SetCursorPos(ImVec2(app.viewport_pos.x, app.viewport_pos.y));
    ImGui::BeginChild("ViewportOverlay", app.viewport_size,
                      false,  // No border
                      ImGuiWindowFlags_NoScrollbar |
                      ImGuiWindowFlags_NoScrollWithMouse |
                      ImGuiWindowFlags_NoBackground);

    // Position toolbar at bottom-right
    float toolbar_w = 200.0f;
    float toolbar_h = 80.0f;
    ImGui::SetCursorPos(ImVec2(app.viewport_size.x - toolbar_w - 10.0f,
                                app.viewport_size.y - toolbar_h - 10.0f));

    // Semi-transparent background
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.1f, 0.1f, 0.1f, 0.7f));
    ImGui::BeginChild("ViewportToolbar", ImVec2(toolbar_w, toolbar_h),
                      true, ImGuiWindowFlags_NoScrollbar);

    // Zoom controls
    ImGui::Text("Zoom: %.1f×", app.viewport_zoom);

    ImGui::PushButtonRepeat(true);  // Allow holding button
    if (ImGui::Button("-", ImVec2(30, 0))) {
        // Zoom out 10%
        app.viewport_zoom *= 0.9f;
        if (app.viewport_zoom < 0.5f) app.viewport_zoom = 0.5f;
        scale = (int)(app.viewport_zoom + 0.5f);
        if (scale < 1) scale = 1;
        render->scale = scale;
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset", ImVec2(50, 0))) {
        app.viewport_zoom = 1.0f;
        scale = 2;
        render->scale = 2;
        app.viewport_pan_x = 0.0f;
        app.viewport_pan_y = 0.0f;
    }
    ImGui::SameLine();
    if (ImGui::Button("+", ImVec2(30, 0))) {
        // Zoom in 10%
        app.viewport_zoom *= 1.1f;
        if (app.viewport_zoom > 32.0f) app.viewport_zoom = 32.0f;
        scale = (int)(app.viewport_zoom + 0.5f);
        render->scale = scale;
    }
    ImGui::PopButtonRepeat();

    // Pan indicator
    ImGui::Text("Pan: (%.0f, %.0f)", app.viewport_pan_x, app.viewport_pan_y);

    // Optional: Grid cell under mouse
    int mx, my;
    SDL_GetMouseState(&mx, &my);
    float local_x = (float)mx - app.viewport_pos.x;
    float local_y = (float)my - app.viewport_pos.y;
    if (local_x >= 0 && local_x < app.viewport_size.x &&
        local_y >= 0 && local_y < app.viewport_size.y) {
        float field_x = local_x - app.viewport_pan_x;
        float field_y = local_y - app.viewport_pan_y;
        int ix = (int)(field_x / (float)scale);
        int iy = (int)(field_y / (float)scale);
        if (ix >= 0 && ix < sim->nx && iy >= 0 && iy < sim->ny) {
            ImGui::Text("Cell: (%d, %d)", ix, iy);
        }
    }

    ImGui::EndChild();
    ImGui::PopStyleColor();

    ImGui::EndChild();  // ViewportOverlay
}
```

---

### 2. Status Bar (Bottom of Window)

```cpp
// After main layout, before ImGui::End()
ImGui::SetCursorPos(ImVec2(origin.x, origin.y + main_h));
ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.08f, 0.08f, 0.09f, 1.00f));
ImGui::BeginChild("StatusBar", ImVec2(full.x, 25.0f), true,
                  ImGuiWindowFlags_NoScrollbar);

// Left side: Simulation info
ImGui::Text("Step: %d | Time: %.2f ns | FPS: %.0f",
            sim ? sim->step : 0,
            sim ? (sim->step * sim->dt * 1e9) : 0.0,
            ImGui::GetIO().Framerate);

// Right side: Viewport info
ImGui::SameLine(full.x - 300.0f);
ImGui::Text("Zoom: %.1f× | Pan: (%.0f, %.0f)",
            app.viewport_zoom,
            app.viewport_pan_x,
            app.viewport_pan_y);

ImGui::EndChild();
ImGui::PopStyleColor();
```

---

### 3. Keyboard Shortcut Help (Optional)

```cpp
// Toggle help overlay with F1 or H key
static bool show_help_overlay = false;

if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_F1) {
    show_help_overlay = !show_help_overlay;
}

if (show_help_overlay && app.viewport_valid) {
    ImGui::SetCursorPos(ImVec2(app.viewport_pos.x + 10, app.viewport_pos.y + 10));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.0f, 0.0f, 0.0f, 0.8f));
    ImGui::BeginChild("HelpOverlay", ImVec2(300, 200), true);

    ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Viewport Controls");
    ImGui::Separator();
    ImGui::Text("Mouse Wheel: Zoom in/out");
    ImGui::Text("Middle Drag: Pan view");
    ImGui::Text("Shift+Right Drag: Pan (alt)");
    ImGui::Text("Arrow Keys: Pan 10px");
    ImGui::Text("Shift+Arrows: Pan 50px");
    ImGui::Text("Ctrl+±: Zoom in/out");
    ImGui::Text("Ctrl+0 / Home: Reset view");
    ImGui::Separator();
    ImGui::Text("Press F1 to close");

    ImGui::EndChild();
    ImGui::PopStyleColor();
}
```

---

## Testing Checklist

- [ ] Toolbar appears in bottom-right of viewport
- [ ] Semi-transparent background (readable over field)
- [ ] "-" button zooms out
- [ ] "+" button zooms in
- [ ] "Reset" button resets zoom and pan
- [ ] Zoom level displayed correctly (0.5× to 32×)
- [ ] Pan offset displayed correctly
- [ ] Cell coordinates show under mouse (optional)
- [ ] Status bar shows simulation step, time, FPS
- [ ] Status bar shows zoom and pan on right side
- [ ] F1 toggles help overlay
- [ ] Help overlay shows all shortcuts
- [ ] Toolbar doesn't block important viewport area
- [ ] Works at all window sizes (resize test)

---

## Success Criteria

1. ✅ Floating toolbar visible in viewport
2. ✅ Zoom controls (+/-/Reset) functional
3. ✅ Zoom and pan status displayed
4. ✅ Status bar at bottom of window
5. ✅ Help overlay toggleable (F1)
6. ✅ Professional, polished appearance

---

**Phase 2.75A COMPLETE after this prompt!**

Next: Phase 2.75B (Docking UI)

---

*Created: 2025-11-20*
