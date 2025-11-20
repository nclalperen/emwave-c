# Prompt #37: Multi-Viewport Split View

**Phase:** 2.75D (Week 1)  
**Effort:** ~8 hours  
**Priority:** High  
**Status:** Ready to implement  
**Dependencies:** Phase 2.75A-C navigation/ruler/grid completed

---

## Objective

Deliver Blender-style viewport layouts with four modes (Single, Horizontal, Vertical, Quad), per-viewport transforms, and quick layout hotkeys.

---

## Visual Layouts (Ascii)

```
Single (Alt+1)                Horizontal (Alt+2)
+--------------------+        +--------------------+
|                    |        |        A           |
|     Full view      |        +--------------------+
|                    |        |        B           |
+--------------------+        +--------------------+

Vertical (Alt+3)               Quad (Alt+4)
+-----------+-----------+     +-----------+-----------+
|     A     |     B     |     |     A     |     B     |
|           |           |     +-----------+-----------+
+-----------+-----------+     |     C     |     D     |
                               +-----------+-----------+
```

---

## Data Structures

```cpp
enum ViewportLayout { VIEWPORT_SINGLE, VIEWPORT_HORIZONTAL, VIEWPORT_VERTICAL, VIEWPORT_QUAD };
enum ViewportViz { VIEWPORT_VIZ_FIELD, VIEWPORT_VIZ_MATERIAL, VIEWPORT_VIZ_OVERLAY, VIEWPORT_VIZ_MAG };

struct ViewportInstance {
    ImVec2 pos, size;
    float zoom;
    float pan_x, pan_y;
    ViewportViz viz_mode;
    bool active;
    bool valid;
    bool show_grid;
    bool show_sources;
};

struct AppState {
    ViewportLayout viewport_layout;
    ViewportInstance viewports[4];
    int active_viewport_idx;   // 0-3 or -1
    bool sync_zoom;
    bool sync_pan;
};
```

---

## Layout Manager

```cpp
static void compute_viewport_layout(AppState* app,
                                    ImVec2 pos,
                                    ImVec2 size) {
    const float gap = 4.0f;
    const float min_dim = 200.0f;
    for (int i = 0; i < 4; ++i) app->viewports[i].valid = false;

    switch (app->viewport_layout) {
        case VIEWPORT_SINGLE: {
            app->viewports[0].pos = pos;
            app->viewports[0].size = size;
            app->viewports[0].valid = true;
        } break;
        case VIEWPORT_HORIZONTAL: {
            float h = (size.y - gap) * 0.5f;
            if (h < min_dim) { app->viewport_layout = VIEWPORT_SINGLE; compute_viewport_layout(app, pos, size); return; }
            app->viewports[0] = {ImVec2(pos.x, pos.y), ImVec2(size.x, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
            app->viewports[1] = {ImVec2(pos.x, pos.y + h + gap), ImVec2(size.x, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
        } break;
        case VIEWPORT_VERTICAL: {
            float w = (size.x - gap) * 0.5f;
            if (w < min_dim) { app->viewport_layout = VIEWPORT_SINGLE; compute_viewport_layout(app, pos, size); return; }
            app->viewports[0] = {ImVec2(pos.x, pos.y), ImVec2(w, size.y), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
            app->viewports[1] = {ImVec2(pos.x + w + gap, pos.y), ImVec2(w, size.y), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
        } break;
        case VIEWPORT_QUAD: {
            float w = (size.x - gap) * 0.5f;
            float h = (size.y - gap) * 0.5f;
            if (w < min_dim || h < min_dim) { app->viewport_layout = VIEWPORT_SINGLE; compute_viewport_layout(app, pos, size); return; }
            app->viewports[0] = {ImVec2(pos.x, pos.y), ImVec2(w, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
            app->viewports[1] = {ImVec2(pos.x + w + gap, pos.y), ImVec2(w, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
            app->viewports[2] = {ImVec2(pos.x, pos.y + h + gap), ImVec2(w, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
            app->viewports[3] = {ImVec2(pos.x + w + gap, pos.y + h + gap), ImVec2(w, h), 0,0,0, VIEWPORT_VIZ_FIELD, false, true};
        } break;
    }
}
```

---

## Rendering Loop (Pseudo-code)

```cpp
for (int i = 0; i < 4; ++i) {
    ViewportInstance& vp = app.viewports[i];
    if (!vp.valid) continue;

    SDL_Rect rect = {(int)vp.pos.x, (int)vp.pos.y, (int)vp.size.x, (int)vp.size.y};
    SDL_RenderSetViewport(renderer, &rect);

    // Apply per-viewport transform
    render->scale = (int)std::lround(std::max(1.0f, vp.zoom));
    render->offset_x = vp.pan_x;
    render->offset_y = vp.pan_y;

    // Render content by viz mode
    switch (vp.viz_mode) {
        case VIEWPORT_VIZ_FIELD:      render_field(renderer, sim); break;
        case VIEWPORT_VIZ_MATERIAL:   render_material(renderer, sim); break;
        case VIEWPORT_VIZ_OVERLAY:    render_field(renderer, sim); render_material_overlay(renderer, sim, 0.3f); break;
        case VIEWPORT_VIZ_MAG:        render_field_magnitude(renderer, sim); break;
    }

    // Optional overlays
    if (vp.show_grid) render_grid(renderer, sim, render->scale);
    if (vp.show_sources) render_sources(renderer, sim->sources);

    SDL_RenderSetViewport(renderer, nullptr); // reset for ImGui

    // ImGui overlays (border + label)
    ImDrawList* dl = ImGui::GetForegroundDrawList();
    ImVec2 min = vp.pos;
    ImVec2 max = ImVec2(vp.pos.x + vp.size.x, vp.pos.y + vp.size.y);
    ImU32 color = (app.active_viewport_idx == i) ? IM_COL32(80,160,255,255) : IM_COL32(80,80,80,180);
    dl->AddRect(min, max, color, 0.0f, 0, 2.0f);
    const char* labels[] = {"A","B","C","D"};
    char buf[64]; snprintf(buf, sizeof(buf), "%s: %s", labels[i], (vp.viz_mode == VIEWPORT_VIZ_FIELD ? "Field" :
                                                                   vp.viz_mode == VIEWPORT_VIZ_MATERIAL ? "Material" :
                                                                   vp.viz_mode == VIEWPORT_VIZ_OVERLAY ? "Overlay" : "Magnitude"));
    dl->AddText(ImVec2(min.x + 8, min.y + 8), IM_COL32(230,230,230,255), buf);
}
```

---

## Interaction and Hotkeys

```cpp
// Active viewport follows mouse
int hit = get_viewport_at_mouse(app, mouse_x, mouse_y);
if (hit >= 0) app.active_viewport_idx = hit;

// Mouse wheel zoom on active viewport (with optional sync)
if (wheel && app.active_viewport_idx >= 0) {
    float factor = 1.0f + wheel * 0.1f;
    ViewportInstance& active = app.viewports[app.active_viewport_idx];
    active.zoom = std::clamp(active.zoom * factor, 0.5f, 32.0f);
    if (app.sync_zoom) for (int i = 0; i < 4; ++i) if (app.viewports[i].valid) app.viewports[i].zoom = active.zoom;
}

// Keyboard hotkeys (Alt+1/2/3/4)
if (alt && key == SDLK_1) app.viewport_layout = VIEWPORT_SINGLE;
if (alt && key == SDLK_2) app.viewport_layout = VIEWPORT_HORIZONTAL;
if (alt && key == SDLK_3) app.viewport_layout = VIEWPORT_VERTICAL;
if (alt && key == SDLK_4) app.viewport_layout = VIEWPORT_QUAD;
```

---

## Viewport Settings Panel (ImGui)

- Layout radio buttons: Single / Horizontal / Vertical / Quad  
- Checkboxes: Sync Zoom, Sync Pan  
- Active viewport readout: A-D  
- Per-viewport controls for viz mode, show grid, show sources, zoom slider  

---

## Testing Checklist (16 cases)

- [ ] Alt+1 switches to Single; Alt+2 Horizontal; Alt+3 Vertical; Alt+4 Quad  
- [ ] Hovering a viewport sets it active; border highlights active  
- [ ] Per-viewport zoom changes do not affect others unless sync enabled  
- [ ] Per-viewport pan changes do not affect others unless sync enabled  
- [ ] Visualization mode can differ per viewport (Field/Material/Overlay/Mag)  
- [ ] Viewport labels show correct name and viz mode  
- [ ] Layout falls back to Single when space is too small  
- [ ] Grid toggle works per viewport  
- [ ] Source toggle works per viewport  
- [ ] Ruler/measuring tools operate in active viewport without cross-talk  
- [ ] Docking layout persists after switching layouts and reopening  
- [ ] Minimal performance regression when all 4 viewports are active  
- [ ] Mouse wheel zoom centers on cursor within active viewport  
- [ ] Keyboard pan continues to work in active viewport  
- [ ] No artifacts when resetting viewport to full window  
- [ ] Cleanup: `SDL_RenderSetViewport` reset before ImGui overlays

---

## Success Criteria

1. Four layout modes available via hotkeys and UI controls.  
2. Independent zoom/pan per viewport with optional sync.  
3. Per-viewport visualization modes selectable and persisted during session.  
4. Active viewport clearly indicated on hover with label and outline.  
5. Performance remains acceptable with four viewports enabled.  
6. Fits cleanly into existing Phase 2.75A-C input and rendering patterns.

*Estimated effort: 8 hours*  
