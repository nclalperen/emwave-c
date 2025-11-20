# Prompt #32: Full DockSpace Layout Conversion

**Phase:** 2.75B (Docking UI - Week 2)
**Effort:** ~6 hours
**Priority:** HIGH - Unlock existing docking feature
**Status:** Ready for implementation
**Dependencies:** Phase 2.75A complete

---

## Objective

Convert fixed-layout UI to full DockSpace, enabling users to drag panels anywhere, create tabs, split viewports, and customize layout like Blender/VS Code.

**Key Discovery:** Docking is ALREADY ENABLED (`ImGuiConfigFlags_DockingEnable` line 2128) but unused!

---

## Implementation

### 1. Replace Fixed Layout with DockSpace

**File:** `src/app/main_imgui.cpp`
**Location:** Main window section (lines ~2856-2988)

**Remove:**
```cpp
// DELETE fixed column calculations (lines 2856-2890)
float left_w = target_left_w * width_scale;
float right_w = target_right_w * width_scale;
// ... all fixed positioning
```

**Replace with:**
```cpp
// Use DockSpace over main viewport
ImGuiViewport* viewport = ImGui::GetMainViewport();
ImGui::SetNextWindowPos(viewport->Pos);
ImGui::SetNextWindowSize(viewport->Size);
ImGui::SetNextWindowViewport(viewport->ID);

ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar |
                                ImGuiWindowFlags_NoDocking |
                                ImGuiWindowFlags_NoTitleBar |
                                ImGuiWindowFlags_NoCollapse |
                                ImGuiWindowFlags_NoResize |
                                ImGuiWindowFlags_NoMove |
                                ImGuiWindowFlags_NoBringToFrontOnFocus |
                                ImGuiWindowFlags_NoNavFocus;

ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
ImGui::Begin("DockSpaceWindow", nullptr, window_flags);
ImGui::PopStyleVar(3);

ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
ImGui::DockSpace(dockspace_id,
                 ImVec2(0.0f, 0.0f),
                 ImGuiDockNodeFlags_PassthruCentralNode);

ImGui::End();  // DockSpaceWindow
```

---

### 2. Convert Child Windows to Dockable Windows

**Before (fixed child):**
```cpp
ImGui::BeginChild("LeftColumn", ImVec2(left_w, main_h), true);
draw_sources_panel(&wizard, &bootstrap, sim, &app);
ImGui::EndChild();
```

**After (dockable window):**
```cpp
if (app.show_sources_panel) {
    ImGui::Begin("Sources", &app.show_sources_panel);
    draw_sources_panel_content(&wizard, &bootstrap, sim, &app);
    ImGui::End();
}
```

**Convert all panels:**
- Scene panel
- Sources panel
- Blocks panel
- Material Legend
- Run Controls
- Grid & Domain
- Run Settings
- Probes panel
- Scope window
- FFT window (if separate)
- Material Browser
- S-Parameter window
- Smith Chart window

---

### 3. Central Viewport Node

**Keep field viewport in center (non-dockable):**

```cpp
// After DockSpace creation
ImGui::SetNextWindowDockID(dockspace_id, ImGuiCond_FirstUseEver);
ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
ImGui::Begin("Viewport", nullptr,
             ImGuiWindowFlags_NoScrollbar |
             ImGuiWindowFlags_NoScrollWithMouse);

// SDL field rendering happens here (existing code)
// Viewport overlay toolbar (from Prompt #31)

ImGui::End();
ImGui::PopStyleVar();
```

---

### 4. Default Layout Setup

**First run initialization:**

```cpp
static bool first_layout = true;
if (first_layout) {
    first_layout = false;

    // Split dockspace into regions
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, viewport->Size);

    ImGuiID dock_left, dock_right, dock_bottom, dock_center;

    // Split: Left 20% | Center 60% | Right 20%
    dock_left = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left,
                                            0.2f, nullptr, &dock_center);
    dock_right = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Right,
                                             0.25f, nullptr, &dock_center);

    // Split bottom 30% from center
    dock_bottom = ImGui::DockBuilderSplitNode(dock_center, ImGuiDir_Down,
                                              0.3f, nullptr, &dock_center);

    // Assign windows to docks
    ImGui::DockBuilderDockWindow("Viewport", dock_center);
    ImGui::DockBuilderDockWindow("Sources", dock_left);
    ImGui::DockBuilderDockWindow("Materials / Blocks", dock_left);
    ImGui::DockBuilderDockWindow("Material Legend", dock_left);
    ImGui::DockBuilderDockWindow("Run Controls", dock_right);
    ImGui::DockBuilderDockWindow("Grid & Domain", dock_right);
    ImGui::DockBuilderDockWindow("Run Settings", dock_right);
    ImGui::DockBuilderDockWindow("Scope", dock_bottom);
    ImGui::DockBuilderDockWindow("Log", dock_bottom);

    ImGui::DockBuilderFinish(dockspace_id);
}
```

---

### 5. Menu Bar with View Menu

```cpp
if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("View")) {
        ImGui::MenuItem("Sources", nullptr, &app.show_sources_panel);
        ImGui::MenuItem("Blocks", nullptr, &app.show_blocks_panel);
        ImGui::MenuItem("Material Legend", nullptr, &app.show_material_legend);
        ImGui::MenuItem("Run Controls", nullptr, &app.show_run_panel);
        ImGui::MenuItem("Grid & Domain", nullptr, &app.show_grid_panel);
        ImGui::MenuItem("Run Settings", nullptr, &app.show_run_settings_panel);
        ImGui::MenuItem("Scope", nullptr, &app.show_scope_window);
        ImGui::MenuItem("Log", nullptr, &app.show_log_panel);
        ImGui::Separator();
        if (ImGui::MenuItem("Reset Layout")) {
            first_layout = true;  // Trigger layout rebuild
        }
        ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
}
```

---

## Testing Checklist

- [ ] All panels appear on first launch
- [ ] Default layout: Left panels | Center viewport | Right controls | Bottom scope
- [ ] Can drag panel tabs to reposition
- [ ] Can create tabbed panel groups
- [ ] Can split panels horizontally
- [ ] Can split panels vertically
- [ ] Can undock panels to floating windows
- [ ] Can close panels (X button)
- [ ] View menu shows/hides panels
- [ ] "Reset Layout" restores default
- [ ] Layout persists across sessions (imgui.ini)
- [ ] Viewport always visible in center
- [ ] SDL field rendering still works
- [ ] Mouse interactions still work (source/block placement)

---

## Success Criteria

1. ✅ Full DockSpace implemented
2. ✅ All panels dockable
3. ✅ Default layout sensible
4. ✅ View menu for panel visibility
5. ✅ Reset Layout function
6. ✅ Layout saved in imgui.ini
7. ✅ Blender-like customizability

---

**Next: Prompt #33 (Layout Presets)**

---

*Created: 2025-11-20*
