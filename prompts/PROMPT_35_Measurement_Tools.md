# Prompt #35: Measurement Tools & Context Menu

**Phase:** 2.75C (Polish - Week 3)
**Effort:** ~6 hours
**Priority:** MEDIUM - Research feature
**Status:** Ready for implementation
**Dependencies:** Phase 2.75A-B complete

---

## Objective

Add ruler tool for distance/angle measurement and right-click context menu for quick actions. Essential for research and documentation.

---

## Implementation

### 1. Ruler Tool (R Key)

```cpp
// AppState additions
struct AppState {
    bool ruler_mode;
    bool ruler_first_point_set;
    ImVec2 ruler_point_a;
    ImVec2 ruler_point_b;
};

// Keyboard toggle
if (event.key.keysym.sym == SDLK_r) {
    app.ruler_mode = !app.ruler_mode;
    app.ruler_first_point_set = false;
    ui_log_add(&app, ruler_mode ? "Ruler tool ON (click 2 points)" :
                                   "Ruler tool OFF");
}

// Click handling
if (app.ruler_mode && event.type == SDL_MOUSEBUTTONDOWN) {
    if (event.button.button == SDL_BUTTON_LEFT) {
        if (!app.ruler_first_point_set) {
            app.ruler_point_a = ImVec2(mx, my);
            app.ruler_first_point_set = true;
        } else {
            app.ruler_point_b = ImVec2(mx, my);

            // Calculate distance and angle
            float dx = ruler_point_b.x - ruler_point_a.x;
            float dy = ruler_point_b.y - ruler_point_a.y;
            float distance_px = sqrtf(dx*dx + dy*dy);
            float distance_m = distance_px / (float)scale * (sim->lx / sim->nx);
            float angle_deg = atan2f(dy, dx) * 180.0f / M_PI;

            ui_log_add(&app, "Ruler: %.3f m, angle %.1f°",
                       distance_m, angle_deg);

            // Reset for next measurement
            app.ruler_first_point_set = false;
        }
    }
}

// Draw ruler line
if (app.ruler_mode && app.ruler_first_point_set) {
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    draw_list->AddLine(ruler_point_a, ruler_point_b,
                       IM_COL32(255, 255, 0, 255), 2.0f);

    // Draw label at midpoint
    ImVec2 mid = ImVec2((ruler_point_a.x + ruler_point_b.x) * 0.5f,
                        (ruler_point_a.y + ruler_point_b.y) * 0.5f);
    draw_list->AddText(mid, IM_COL32(255, 255, 0, 255), "Click second point");
}
```

---

### 2. Right-Click Context Menu

```cpp
// In viewport mouse handling
static bool show_context_menu = false;
static ImVec2 context_menu_pos;

if (event.type == SDL_MOUSEBUTTONDOWN &&
    event.button.button == SDL_BUTTON_RIGHT) {

    // Don't show if Shift+Right (pan mode)
    if (!(SDL_GetModState() & KMOD_SHIFT)) {
        show_context_menu = true;
        context_menu_pos = ImVec2((float)mx, (float)my);
        ImGui::OpenPopup("ViewportContextMenu");
    }
}

// Draw context menu
if (ImGui::BeginPopup("ViewportContextMenu")) {
    int ix = ...; // Grid cell at mouse
    int iy = ...;

    if (ImGui::MenuItem("Add Source Here")) {
        create_new_source_at(wizard, sim, &app, ix, iy);
        show_context_menu = false;
    }

    if (ImGui::MenuItem("Add Material Block Here")) {
        create_block_at(wizard, &bootstrap, sim, &app, ix, iy);
        show_context_menu = false;
    }

    ImGui::Separator();

    if (ImGui::MenuItem("Measure Distance (R)")) {
        app.ruler_mode = true;
        show_context_menu = false;
    }

    if (ImGui::MenuItem("Zoom to Fit")) {
        // Reset zoom + pan
        app.viewport_zoom = 1.0f;
        scale = 2;
        app.viewport_pan_x = 0.0f;
        app.viewport_pan_y = 0.0f;
        show_context_menu = false;
    }

    ImGui::Separator();

    if (ImGui::MenuItem("Copy Coordinates")) {
        char buf[64];
        snprintf(buf, sizeof(buf), "(%d, %d)", ix, iy);
        SDL_SetClipboardText(buf);
        show_context_menu = false;
    }

    ImGui::EndPopup();
}
```

---

### 3. Coordinate Display (Cursor Position)

```cpp
// Bottom-left corner overlay
if (app.viewport_valid && mouse_in_viewport) {
    ImGui::SetCursorPos(ImVec2(app.viewport_pos.x + 10,
                               app.viewport_pos.y + app.viewport_size.y - 30));

    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.0f, 0.0f, 0.0f, 0.6f));
    ImGui::BeginChild("CursorInfo", ImVec2(200, 25), true);

    int ix = ...; // Grid cell under mouse
    int iy = ...;
    float x_m = ix * (sim->lx / sim->nx);
    float y_m = iy * (sim->ly / sim->ny);

    ImGui::Text("Cell: (%d, %d) | %.3f, %.3f m", ix, iy, x_m, y_m);

    ImGui::EndChild();
    ImGui::PopStyleColor();
}
```

---

## Testing Checklist

- [ ] Press R - ruler tool activates
- [ ] Click first point - line draws from point to cursor
- [ ] Click second point - distance and angle displayed
- [ ] Measurements accurate (verify with known distances)
- [ ] Right-click in viewport - context menu appears
- [ ] "Add Source Here" - creates source at click location
- [ ] "Add Material Block Here" - creates block
- [ ] "Measure Distance" - activates ruler tool
- [ ] "Zoom to Fit" - resets view
- [ ] "Copy Coordinates" - copies to clipboard
- [ ] Cursor position display shows correct cell and meters
- [ ] Works at all zoom levels

---

## Success Criteria

1. ✅ Ruler tool measures distance and angle
2. ✅ Right-click context menu functional
3. ✅ Quick actions accessible
4. ✅ Coordinate display under cursor
5. ✅ Professional research tool

---

*Created: 2025-11-20*
