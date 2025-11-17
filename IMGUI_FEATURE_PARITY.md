# ImGui Feature Parity Checklist

This document tracks feature parity between the SDL2 UI ([main_new.c](src/app/main_new.c)) and the ImGui UI ([main_imgui.cpp](src/app/main_imgui.cpp)).

## ✅ Features Already Implemented in ImGui

| Feature | SDL2 | ImGui | Notes |
|---------|------|-------|-------|
| **Basic Simulation** |
| Pause/Resume | ✅ Space | ✅ Space | Working |
| Quit | ✅ ESC/Q | ✅ ESC | Working |
| Step counter display | ✅ | ✅ | In status panel |
| FPS display | ✅ | ✅ | In status panel |
| **Source Management** |
| Source placement by click | ✅ | ✅ | Interactive tools panel |
| Source selection | ✅ | ✅ | Combo box |
| Source list display | ✅ | ✅ | Sources panel |
| Source position display | ✅ | ✅ | Sources panel |
| Source type selection | ✅ | ✅ | Dropdown in wizard |
| Source amplitude editing | ✅ | ✅ | Wizard/Sources panel |
| Source frequency editing | ✅ | ✅ | Wizard/Sources panel |
| **Material/Block Management** |
| Block drawing (2-click) | ✅ | ✅ | Interactive tools panel |
| Block list display | ✅ | ✅ | Blocks panel |
| Block property editing | ✅ | ✅ | Blocks panel |
| Material type selection | ✅ | ✅ | Dropdown (PEC/PMC/dielectric) |
| **Visualization** |
| Field heatmap | ✅ | ✅ | Viewport |
| Source indicators | ✅ | ✅ | Overlay on viewport |
| Block outlines | ✅ | ✅ | Overlay on viewport |
| Grid overlay | ✅ | ✅ | Toggle button |
| Oscilloscope | ✅ | ✅ | Scope tab (basic) |
| **Configuration** |
| Wizard/Setup UI | ❌ | ✅ | ImGui advantage! |
| Grid configuration | ❌ | ✅ | Wizard tabs |
| Boundary type selection | ✅ | ✅ | Wizard |
| CFL safety setting | ❌ | ✅ | Wizard |
| Run mode (fixed/sweep) | ❌ | ✅ | Wizard |
| **Panels** |
| Scene info panel | ✅ | ✅ | Left column |
| Sources panel | ✅ | ✅ | Left column |
| Blocks panel | ✅ | ✅ | Left column |
| Simulation controls | ✅ | ✅ | Right column |
| Log panel | ❌ | ✅ | Bottom tabs |

## ❌ Features Missing in ImGui (Need Implementation)

### Critical (Phase 1 - Required for Feature Parity)

| Feature | SDL2 | ImGui | Priority | Effort |
|---------|------|-------|----------|--------|
| **Interactive Controls** |
| Frequency slider | ✅ | ❌ | **HIGH** | 2h |
| Steps/frame slider | ✅ | ❌ | **HIGH** | 1h |
| Source dragging | ✅ | ❌ | **HIGH** | 2h |
| Material painting/brush | ✅ | ❌ | **HIGH** | 4h |
| Probe positioning | ✅ | ❌ | **MED** | 1h |
| **Keyboard Shortcuts** |
| 1-3: Toggle sources | ✅ | ❌ | **HIGH** | 1h |
| T: Cycle source type | ✅ | ❌ | **MED** | 1h |
| F: Source type (alt) | ✅ | ❌ | **LOW** | 0.5h |
| R: Reset simulation | ✅ | ❌ | **HIGH** | 0.5h |
| C: Clear fields | ✅ | ❌ | **HIGH** | 0.5h |
| Y: Toggle CPML/Mur | ✅ | ❌ | **MED** | 1h |
| 7-9: CPML presets | ✅ | ❌ | **MED** | 1h |
| S: Toggle S-param ports | ✅ | ❌ | **MED** | 0.5h |
| **Visual Controls** |
| B: Theme switch (dark/light) | ✅ | ❌ | **HIGH** | 2h |
| V/Shift+V: Accent palette | ✅ | ❌ | **MED** | 1h |
| K/Shift+K: Colormap | ✅ | ❌ | **HIGH** | 2h |
| H: Hold color scale | ✅ | ❌ | **MED** | 1h |
| J: Hold scope Y-axis | ✅ | ❌ | **MED** | 1h |
| A: Auto-rescale toggle | ✅ | ❌ | **LOW** | 0.5h |
| L: Toggle legend | ✅ | ❌ | **LOW** | 0.5h |
| **Help & Documentation** |
| F1: Help overlay | ✅ | ❌ | **HIGH** | 3h |
| Legend panel | ✅ | ❌ | **MED** | 2h |
| Status badges | ✅ | ❌ | **MED** | 2h |
| **File Operations** |
| F2: Screenshot | ✅ | ❌ | **MED** | 1h |
| F3: FFT export | ✅ | ❌ | **LOW** | 2h |
| F5/F6: Scene presets | ✅ | ❌ | **MED** | 2h |
| **Layout** |
| [/]: Resize left panel | ✅ | ❌ | **LOW** | n/a (ImGui docking) |
| ,/.: Resize right panel | ✅ | ❌ | **LOW** | n/a (ImGui docking) |
| -/=: Resize timeline | ✅ | ❌ | **LOW** | n/a (ImGui docking) |
| \\: Scope docking | ✅ | ❌ | **MED** | 2h |
| **Paint Mode** |
| M/U: Toggle paint | ✅ | ❌ | **HIGH** | 1h |
| I: Cycle paint type | ✅ | ❌ | **HIGH** | 1h |
| O/P: Paint eps -/+ | ✅ | ❌ | **MED** | 1h |
| Mouse paint dragging | ✅ | ❌ | **HIGH** | 2h |
| **Advanced** |
| G: Probe logging toggle | ✅ | ❌ | **LOW** | 0.5h |

**Total Critical Items:** 35 features
**Estimated Effort:** ~40 hours

### Nice-to-Have (Phase 2 - Enhancements Beyond SDL2)

| Feature | SDL2 | ImGui | Priority | Effort |
|---------|------|-------|----------|--------|
| **ImGui Advantages** |
| Material browser with presets | ❌ | ❌ TODO | **HIGH** | 8h |
| ImPlot oscilloscope | ❌ | ❌ TODO | **HIGH** | 6h |
| ImPlot FFT panel | ❌ | ❌ TODO | **HIGH** | 6h |
| Smith chart (S-params) | ❌ | ❌ TODO | **MED** | 8h |
| Scene gallery with thumbnails | ❌ | ❌ TODO | **MED** | 6h |
| File dialogs (save/load) | ❌ | ❌ TODO | **MED** | 4h |
| Tooltips everywhere | ❌ | ❌ TODO | **HIGH** | 4h |
| Docking layout persistence | ❌ | ❌ TODO | **MED** | 3h |
| Keyboard shortcut panel | ❌ | ❌ TODO | **MED** | 3h |
| Welcome screen/tutorial | ❌ | ❌ TODO | **LOW** | 6h |
| Expression editor (syntax HL) | ❌ | ❌ TODO | **LOW** | 8h |
| Material property graphs | ❌ | ❌ TODO | **LOW** | 4h |
| Animation timeline controls | ❌ | ❌ TODO | **LOW** | 6h |

**Total Enhancement Items:** 13 features
**Estimated Effort:** ~72 hours

## 📊 Feature Parity Summary

| Category | SDL2 Features | ImGui Has | Missing | Parity % |
|----------|---------------|-----------|---------|----------|
| **Basic Simulation** | 4 | 4 | 0 | 100% |
| **Source Management** | 7 | 7 | 0 | 100% |
| **Material/Block Management** | 5 | 5 | 0 | 100% |
| **Visualization** | 5 | 5 | 0 | 100% |
| **Configuration** | 5 | 5 | 0 | 100% |
| **Interactive Controls** | 5 | 0 | 5 | 0% |
| **Keyboard Shortcuts** | 12 | 2 | 10 | 17% |
| **Visual Controls** | 7 | 0 | 7 | 0% |
| **Help & Documentation** | 3 | 0 | 3 | 0% |
| **File Operations** | 3 | 0 | 3 | 0% |
| **Layout** | 4 | 0 | 4 | 0% |
| **Paint Mode** | 4 | 0 | 4 | 0% |
| **Advanced** | 1 | 0 | 1 | 0% |
| **TOTAL** | **60** | **28** | **37** | **47%** |

## 🎯 Implementation Roadmap

### Week 1: Core Interactions (16h)
1. ✅ Feature audit (DONE)
2. Frequency & speed sliders (3h)
3. All keyboard shortcuts (8h)
4. Source dragging (2h)
5. Theme & colormap switching (3h)

### Week 2: Material System (16h)
1. Material painting/brush mode (6h)
2. Paint keyboard shortcuts (M/U/I/O/P) (2h)
3. Material preset library (4h)
4. Material browser UI (4h)

### Week 3: Help & Visualization (16h)
1. Help overlay (F1) (4h)
2. Legend panel (2h)
3. Status badges (2h)
4. Screenshot & FFT export (3h)
5. ImPlot oscilloscope (5h)

### Week 4: Advanced Features (16h)
1. ImPlot FFT panel (6h)
2. Scene gallery (6h)
3. Tooltips (4h)

### Week 5: Polish & Testing (16h)
1. Docking layout persistence (3h)
2. Welcome screen (6h)
3. Performance testing (3h)
4. Cross-platform builds (4h)

**Total: 80 hours (5 weeks full-time, or 10 weeks part-time)**

## 📝 Notes

### Why ImGui is Better (Long-Term)
- **Docking**: Users can arrange panels freely
- **Tabs**: Multiple views in one area (Scope/FFT/Log)
- **Widgets**: Rich controls (color pickers, combos, tables)
- **Tooltips**: Contextual help everywhere
- **ImPlot**: Professional plots with zoom/pan
- **Styling**: Consistent theming via ImGui::StyleColors*

### Migration Strategy
1. **Phase 1**: Achieve feature parity (40h)
2. **Phase 2**: Add ImGui enhancements (72h)
3. **Phase 3**: Deprecate SDL2 UI (mark as legacy)
4. **Phase 4**: Update all documentation

### Testing Checklist
- [ ] All keyboard shortcuts work
- [ ] Material painting works
- [ ] Source dragging works
- [ ] Sliders respond correctly
- [ ] Theme switching applies immediately
- [ ] Colormap cycling updates heatmap
- [ ] Screenshot saves correctly
- [ ] FFT export generates CSV
- [ ] Scene presets load correctly
- [ ] Wizard Apply & Restart works
- [ ] ImPlot scope shows waveform
- [ ] ImPlot FFT shows spectrum
- [ ] Tooltips appear on hover
- [ ] Docking layout persists
- [ ] Welcome screen shows on first run
- [ ] Performance matches SDL2 (60+ FPS)
- [ ] Windows build works
- [ ] Linux build works
- [ ] macOS build works

## 🐛 Known Issues
- Docking is disabled (line 668: `io.ConfigFlags &= ~ImGuiConfigFlags_DockingEnable`)
- No IniFilename (line 670: `io.IniFilename = NULL`) - layout not persisted
- Fixed window size (1920x1080 minimum) - should be resizable
- Scope is simple plot, not ImPlot
- No tooltips anywhere
- No keyboard shortcut documentation visible

## 🚀 Quick Wins (Do These First)
1. Enable docking: `io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;`
2. Enable layout persistence: `io.IniFilename = "imgui.ini";`
3. Add frequency slider in "Simulation Controls"
4. Add steps/frame slider in "Simulation Controls"
5. Add keyboard shortcut handler in main loop (before ImGui::NewFrame())
6. Add tooltips with `ImGui::SetItemTooltip("...")` after each widget
7. Replace hardcoded window size with dynamic calculation
8. Add theme toggle button
9. Add colormap cycling buttons
10. Add "Help" menu with F1 overlay
