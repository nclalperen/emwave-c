# Prompt #33: Layout Presets & Save/Load

**Phase:** 2.75B (Docking UI - Week 2 - Final)
**Effort:** ~3 hours
**Priority:** MEDIUM - User convenience
**Status:** Ready for implementation
**Dependencies:** Prompt #32 (DockSpace)

---

## Objective

Add 3 preset layouts (Beginner/Power User/Analysis) and save/load custom layouts to files. Users can quickly switch between workflow-optimized configurations.

---

## Implementation

### 1. Layout Presets

```cpp
enum LayoutPreset {
    LAYOUT_BEGINNER = 0,
    LAYOUT_POWER_USER = 1,
    LAYOUT_ANALYSIS = 2,
    LAYOUT_CUSTOM = 3
};

static void apply_layout_preset(ImGuiID dockspace_id, LayoutPreset preset) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);

    ImGuiID dock_left, dock_right, dock_bottom, dock_center;

    switch (preset) {
        case LAYOUT_BEGINNER:
            // Simple: Viewport | Controls
            dock_right = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Right,
                                                     0.25f, nullptr, &dock_center);
            dock_bottom = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Down,
                                                      0.25f, nullptr, &dock_center);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Run Controls", dock_right);
            ImGui::DockBuilderDockWindow("Sources", dock_right);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            break;

        case LAYOUT_POWER_USER:
            // 3-column: Tools | Viewport | Controls
            dock_left = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left,
                                                    0.2f, nullptr, &dock_center);
            dock_right = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Right,
                                                     0.2f, nullptr, &dock_center);
            dock_bottom = ImGui::DockBuilderSplitNode(dock_left, ImGuiDir_Down,
                                                      0.5f, nullptr, &dock_left);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Sources", dock_left);
            ImGui::DockBuilderDockWindow("Materials / Blocks", dock_left);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("Material Legend", dock_right);
            ImGui::DockBuilderDockWindow("Run Controls", dock_right);
            break;

        case LAYOUT_ANALYSIS:
            // Analysis: Viewport + Plots side-by-side
            dock_right = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Right,
                                                     0.5f, nullptr, &dock_center);
            dock_bottom = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Down,
                                                      0.5f, nullptr, &dock_center);
            ImGuiID dock_right_bottom = ImGui::DockBuilderSplitNode(dock_right,
                                                                    ImGuiDir_Down,
                                                                    0.5f, nullptr,
                                                                    &dock_right);

            ImGui::DockBuilderDockWindow("Viewport", dock_center);
            ImGui::DockBuilderDockWindow("Scope", dock_bottom);
            ImGui::DockBuilderDockWindow("S-Parameters", dock_right);
            ImGui::DockBuilderDockWindow("Smith Chart", dock_right_bottom);
            break;
    }

    ImGui::DockBuilderFinish(dockspace_id);
}
```

---

### 2. Layout Menu

```cpp
// In menu bar
if (ImGui::BeginMenu("Layout")) {
    if (ImGui::MenuItem("Beginner")) {
        apply_layout_preset(dockspace_id, LAYOUT_BEGINNER);
    }
    if (ImGui::MenuItem("Power User")) {
        apply_layout_preset(dockspace_id, LAYOUT_POWER_USER);
    }
    if (ImGui::MenuItem("Analysis")) {
        apply_layout_preset(dockspace_id, LAYOUT_ANALYSIS);
    }
    ImGui::Separator();
    if (ImGui::MenuItem("Save Layout...")) {
        // Save imgui.ini to custom file
        ImGui::SaveIniSettingsToDisk("layouts/my_layout.ini");
    }
    if (ImGui::MenuItem("Load Layout...")) {
        // Load from custom file
        ImGui::LoadIniSettingsFromDisk("layouts/my_layout.ini");
    }
    ImGui::EndMenu();
}
```

---

### 3. Layout Hotkeys

```cpp
// F2, F3, F4 for quick layout switching
if (event.type == SDL_KEYDOWN) {
    if (event.key.keysym.sym == SDLK_F2) {
        apply_layout_preset(dockspace_id, LAYOUT_BEGINNER);
        ui_log_add(&app, "Layout: Beginner");
    }
    else if (event.key.keysym.sym == SDLK_F3) {
        apply_layout_preset(dockspace_id, LAYOUT_POWER_USER);
        ui_log_add(&app, "Layout: Power User");
    }
    else if (event.key.keysym.sym == SDLK_F4) {
        apply_layout_preset(dockspace_id, LAYOUT_ANALYSIS);
        ui_log_add(&app, "Layout: Analysis");
    }
}
```

---

## Testing Checklist

- [ ] Layout → Beginner - simple 2-column layout
- [ ] Layout → Power User - 3-column layout
- [ ] Layout → Analysis - viewport + plots
- [ ] F2/F3/F4 hotkeys switch layouts
- [ ] Save Layout saves to file
- [ ] Load Layout restores from file
- [ ] Layout persists across restarts (imgui.ini)
- [ ] Custom layouts saveable
- [ ] Multiple layout files supported

---

## Success Criteria

1. ✅ 3 preset layouts implemented
2. ✅ Layout menu functional
3. ✅ F2/F3/F4 hotkeys work
4. ✅ Save/Load custom layouts
5. ✅ Professional workflow support

---

**Phase 2.75B COMPLETE!**
**Next: Phase 2.75C (Polish)**

---

*Created: 2025-11-20*
