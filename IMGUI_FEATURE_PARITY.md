# ImGui Feature Parity Checklist

Status (2025-11-20): SDL2 feature parity achieved in `main_imgui.cpp`. The items that were flagged as missing (paint mode, visual controls, help overlay, frequency/speed controls, keyboard shortcuts) are implemented and wired to the solver/UI.

## Parity Highlights
- Interactive controls: frequency slider and steps-per-frame slider in the Simulation Controls panel; source dragging works in the viewport; block placement uses two-click flow; probes panel handles logging toggles.
- Paint mode: `M`/`U` toggle paint mode; `I` cycles legacy paint type when no material is selected; `O`/`P` adjust dielectric epsilon; mouse drag paints; Numpad 1-8 picks materials from the library; material-browser selection sets `paint_material_id`.
- Visual controls: `B` switches theme (dark/light), `K`/`Shift+K` cycles colormaps, `V`/`Shift+V` cycles accent palettes, `A/H/J/L` set auto/hold modes, `G` toggles grid overlay, `E` cycles Field/Material/Overlay view.
- Help system: `F1` opens the full-screen help/shortcuts overlay; also reachable via Help -> Keyboard Shortcuts.
- File operations: `F2` screenshot, `F3` FFT export, `F5/F6` load presets, `Ctrl+S` saves config; layout persists via `imgui.ini`.
- Paint + materials coexist: paint uses the active material or legacy types; material blocks remain separate and filter/auto-apply correctly.

## Code References (main_imgui.cpp, UTF-8)
- SDL shortcut handling: lines ~2189-2513 (paint mode, visual controls, presets, save, boundary toggles, source toggles).
- Paint brush application: lines ~1503-1596 (brush helpers) and ~2558-2684 (click/drag painting paths).
- Help overlay window: lines ~3359-3447.
- Material legend toggle: lines ~981-1015 and menu hook ~2828-2930.
- Theme/colormap/accent application: initialization around lines ~2159-2175 and updates in the shortcut handler.

## Remaining Gaps vs SDL2 list
- None blocking parity. Probe logging is controlled from the Probes panel instead of a dedicated key; layout resize hotkeys are replaced by ImGui docking/layout persistence.

## Validation Notes
- Release build and `ctest -C Release --output-on-failure` pass (see `validation_log.txt` for latest run).
- Manual UI interactions (paint, visual controls, help overlay) should be exercised during interactive validation to confirm shortcuts with/without ImGui focus.
