# Current Status & Findings (ImGui Front-End)

## Summary
- Wizard UI fully decomposed into dedicated panels (Grid, Run Settings, Sources, Materials/Blocks).
- Material legend integrates with block creation and selection; viewport outlines highlight blocks that match the selected material.
- Build verified: `cmake --build build-imgui --config Debug`.

## Issues / Gaps
- **Highlight toggle ignored**: `app.highlight_blocks_by_material` is never checked in the viewport overlay; matching blocks are always highlighted when a material is selected. Outcome: highlight cannot be turned off despite View-menu toggle.
- **Material filter auto-enable**: Selecting a material forces `filter_blocks_by_material = true` even if the user had disabled filtering. Consider respecting the toggle when `auto_filter_blocks_on_select` is off (currently default is true, but manual toggle may be overridden on selection).
- **Block filtering continue path**: When filtering, the code calls `ImGui::PopID()` inside the skip branch (good), but still leaves the filter toggle enabled even when no material is selected; we already reset `filter_blocks_by_material = false` when no selection, so this is OK—just note the behavior.
- **Material property mapping**: `apply_material_to_rect` sets dielectric `sigma` from the library’s conductivity field (most dielectrics have 0 except Silicon). Loss tangent isn’t propagated, so lossy dielectrics from the library won’t reflect tan δ in blocks. If tan δ-based sigma is desired, extra mapping is needed.

## Suggested Fixes
- Gate viewport block highlighting with `app.highlight_blocks_by_material` so users can disable it.
- When selecting a material, only enable block filtering if `auto_filter_blocks_on_select` is true; otherwise preserve the user’s filter setting.
- Consider mapping tan δ to an effective conductivity when applying dielectric materials to blocks (if desired for simulation fidelity).

## Recent Additions (for context)
- Legend row “Add” inserts a block with the chosen material and logs it.
- Blocks panel filter toggle and auto-filter option in View menu.
- Blocks list shows material swatches; matching legend selection highlights blocks in the viewport.
