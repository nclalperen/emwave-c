# Phase 2.5: Wizard Function Integration Plan

**Date:** 2025-11-19
**Status:** COMPLETE (Phase 2.5 delivered)
**Objective:** Reintegrate wizard functions into modern Phase 2 UI without recreating wizard window

---

## Strategic Context

After completing Phase 2 (Professional Transformation), we identified that valuable wizard functions exist but are hidden. Rather than revive the monolithic wizard window, we'll **redistribute functionality** into existing panels where users expect to find them.

**Key Principle:** Modern UI paradigm - put controls where they're needed, not in a separate configuration wizard.

**Completion Summary (2025-11-20):**
- Sources, Materials/Blocks, and Simulation Setup are integrated into the ImGui panels with create/edit/delete flows and live apply paths.
- Material blocks pull directly from the Phase 2 material library, support default creation and from-library creation, allow filtering/swatching, and auto-apply to the solver.
- Grid/domain/boundary controls are exposed in the setup panel with restart handling; documentation/validation now tracked for sign-off.

---

## Integration Roadmap

### **Week 1: Sources Panel Enhancement** (~6h)
**Prompt #24**

**Objective:** Transform Sources panel from basic display into full source editor with expression support

**Current State:**
- Sources panel shows existing sources (read-only mostly)
- Mouse placement works but limited
- No create/delete UI
- No expression editor (despite existing code!)

**Target State:**
- Create/delete sources with buttons
- Full property editor (type, field, position, amplitude, frequency)
- **Expression editor with live preview** (extract from wizard_sources_tab)
- Active/inactive toggle per source
- Drag-to-reposition (if not already working)

**Files to Modify:**
- `src/app/main_imgui.cpp` - draw_sources_panel() function
- Extract expression editor logic from wizard_sources_tab() (lines 1371-1415)

**Key Features:**
1. **Source List with Controls**
   ```
   Source 0: CW @ 2.4 GHz [✓ Active] [Edit] [Delete]
   Source 1: Gaussian [  Inactive] [Edit] [Delete]
   [+ New Source]
   ```

2. **Property Editor** (when "Edit" clicked)
   ```
   Type: [CW ▼] [Gaussian] [Ricker] [Expression]
   Field Component: [Ez ▼] [Hx] [Hy]
   Position (normalized): X: 0.25  Y: 0.50
   Amplitude: 1.0
   Frequency: 2.4e9 Hz (2.40 GHz)
   Sigma²: 1.0 (for Gaussian/Ricker)
   ```

3. **Expression Editor** (for SRC_EXPR type)
   ```
   Expression:
   ┌────────────────────────────────────┐
   │ amp * sin(2*pi*freq*t)             │
   │   * exp(-(t-1e-9)^2/sigma2)        │
   └────────────────────────────────────┘
   Variables: t (time), amp, freq, pi

   Preview:
   [Live waveform plot - 128 samples]
   ```

**Significance:**
- **Unique capability** - expression-based sources don't exist in other EM tools
- **Educational value** - students visualize signal shapes before running
- **Research capability** - custom excitations for S-parameter characterization
- **UX improvement** - no more wizard, edit sources in context

**Estimated Effort:** 6 hours
- Extract expression editor: 2h
- Source CRUD UI: 2h
- Property editor layout: 1h
- Testing & polish: 1h

---

### **Week 2: Material Blocks Integration** (~4h)
**Prompt #25**

**Objective:** Integrate structured material rectangles with Phase 2 material library

**Current State:**
- Material library exists (11 materials)
- Material browser works
- Paint mode works (freehand)
- Blocks panel exists but basic
- Wizard materials tab hidden (precise rectangles)

**Target State:**
- "Material Blocks" section in Blocks panel
- Create rectangles using library materials
- Edit coordinates (x0, y0, x1, y1) with numeric precision
- Material dropdown populated from library
- Apply changes updates simulation

**Files to Modify:**
- `src/app/main_imgui.cpp` - draw_blocks_panel() function
- Link to material_library API

**Key Features:**
1. **Block List with Material Names**
   ```
   ┏━━━━━━━━━━━━━━━━━━━━━━━┓
   ┃ Material Blocks       ┃
   ┗━━━━━━━━━━━━━━━━━━━━━━━┛

   Block 0: Copper Ground Plane
     Material: [Copper ▼] (σ = 5.96e7 S/m)
     Region: (0.00, 0.00) → (0.10, 0.01) m
     [Edit Coords] [Delete]

   Block 1: FR4 Substrate
     Material: [FR4 ▼] (εᵣ = 4.4, tan δ = 0.02)
     Region: (0.00, 0.01) → (0.10, 0.03) m
     [Edit Coords] [Delete]

   [+ Add Block from Library]
   ```

2. **Coordinate Editor**
   ```
   Edit Block 0 Coordinates:
   x0: 0.000 m    x1: 0.100 m
   y0: 0.000 m    y1: 0.010 m
   [Apply] [Cancel]
   ```

3. **Material Picker** (dropdown or browser button)
   ```
   Material: [Copper ▼]
   ─────────────────
   Metals:
     • Copper (σ = 5.96e7 S/m)
     • Gold
     • Silver
     • Aluminum
   Dielectrics:
     • Air (εᵣ = 1.0)
     • FR4 (εᵣ = 4.4)
     • Rogers RO4003
     ...
   ```

**Integration Strategy:**
- Use existing `MaterialRectSpec` from wizard.cfg
- Populate material properties from library by ID/name
- Store both rect coords AND material library reference
- On "Apply", call existing `apply_wizard_materials_to_sim()`

**Significance:**
- **Precision** - rectangles are exact, paint is approximate
- **Reproducibility** - save/load exact geometries to JSON
- **Professional** - microstrip/waveguide design needs numeric accuracy
- **Synergy** - material library (properties) + rectangles (geometry) = complete system

**Estimated Effort:** 4 hours
- UI layout: 1.5h
- Material library integration: 1.5h
- Coordinate editor: 0.5h
- Testing: 0.5h

---

### **Week 3: Simulation Setup Panel** (~3h)
**Prompt #26**

**Objective:** Expose grid/domain/boundary configuration without full wizard

**Current State:**
- Grid size fixed at startup
- Boundary toggle via 'Y' key only
- CFL not exposed
- No way to change domain size in UI

**Target State:**
- "Simulation Setup" collapsing header in right column
- Grid size (nx, ny) with warning on change
- Physical domain (lx, ly) in meters
- CFL safety slider (0.1-0.99)
- Boundary mode toggle (CPML ↔ Mur)
- "Apply & Restart" button (clears fields)

**Files to Modify:**
- `src/app/main_imgui.cpp` - new draw_simulation_setup_panel()
- Extract from wizard_grid_tab() (lines 1270-1294)

**Key Features:**
1. **Setup Panel Layout**
   ```
   ┏━━━━━━━━━━━━━━━━━━━━━━━┓
   ┃ Simulation Setup      ┃  [Collapse ▼]
   ┗━━━━━━━━━━━━━━━━━━━━━━━┛

   Grid Size:
     nx: [256   ] cells  (X-direction)
     ny: [256   ] cells  (Y-direction)

   Physical Domain:
     lx: [0.100 ] m  (width)
     ly: [0.100 ] m  (height)

   Cell size: dx = 0.39 mm, dy = 0.39 mm

   CFL Safety: [████████████░░] 0.95

   Boundary Conditions:
     ( ) Mur (1st-order ABC)
     (•) CPML (PML absorbing layers)

   ⚠ Changes require restart
   [Apply & Restart Simulation]
   ```

2. **Validation & Feedback**
   - Calculate and display cell size (dx, dy)
   - Warn if grid too large (memory estimate)
   - Show CFL limit (c·dt < dx/√2)
   - Confirm before restart (clear fields warning)

**Integration Strategy:**
- Use existing WizardState.cfg for config storage
- On "Apply", call bootstrap_sim_from_config()
- Preserve sources and materials where possible
- Log to UI: "Restarted with 512×512 grid, CPML boundaries"

**Significance:**
- **Experimentation** - try different grid resolutions
- **Educational** - students learn CFL stability condition
- **Flexibility** - change problem size without editing JSON
- **Professional** - quick convergence studies

**Estimated Effort:** 3 hours
- UI layout: 1h
- Validation logic: 0.5h
- Restart integration: 1h
- Testing: 0.5h

---

### **Phase 3: Analysis & Automation** (~8h)
**Prompt #27** (Future)

**Objective:** In-app frequency sweep launcher and result viewer

**Deferred to Phase 3 because:**
- Requires process management (launch CLI as subprocess)
- Progress tracking UI needed
- Most complex integration
- Current workaround (manual CLI sweep) is acceptable

**Planned Features:**
1. Sweep configuration UI (start/stop freq, points, steps/point)
2. "Launch Sweep" button (spawns `emwave_cli --mode sweep`)
3. Progress bar with real-time updates
4. Auto-load results into S-parameter window when complete
5. Batch sweep queue (multiple runs)

**Will integrate:** wizard_run_tab() logic (lines 1420-1451)

---

## Implementation Schedule

| Week | Prompt | Feature | Effort | Priority |
|------|--------|---------|--------|----------|
| **Week 1** | #24 | Sources Panel Enhancement | 6h | 🔴 HIGH |
| **Week 2** | #25 | Material Blocks Integration | 4h | 🔴 HIGH |
| **Week 3** | #26 | Simulation Setup Panel | 3h | 🟡 MEDIUM |
| **Phase 3** | #27 | Analysis & Automation | 8h | 🟢 FUTURE |

**Total Phase 2.5 Effort:** 13 hours (Weeks 1-3)
**Phase 3 Addition:** +8 hours (Week 4+)

---

## Success Criteria

### Week 1 Complete When:
- [x] Sources can be created/deleted via UI
- [x] Expression editor works with live preview
- [x] Source properties editable (type, freq, amp, position)
- [x] Active/inactive toggle works per source
- [x] No need to use wizard for source configuration

### Week 2 Complete When:
- [x] Material blocks list shows material library names
- [x] Blocks can be created with library material selection
- [x] Coordinates editable with numeric precision
- [x] "Apply" updates simulation with correct materials
- [x] Rectangles and paint mode coexist

### Week 3 Complete When:
- [x] Grid size changeable via UI
- [x] Domain dimensions editable
- [x] Boundary mode toggleable (CPML ↔ Mur)
- [x] CFL slider functional
- [x] "Apply & Restart" rebuilds simulation correctly
- [x] Cell size displayed and validated

### Phase 2.5 Complete When:
- [x] All wizard tab functions accessible in modern UI
- [x] No need for wizard window
- [x] User workflow smoother than before
- [x] Documentation updated

---

## Testing Plan

### Integration Testing (after each week):
1. **Sources Week:**
   - Create CW source at 2.4 GHz via UI
   - Create expression source: `amp*sin(2*pi*freq*t)*exp(-t/1e-9)`
   - Verify preview plot shows decaying sinusoid
   - Toggle source inactive, verify simulation stops excitation
   - Delete source, verify removed from list

2. **Materials Week:**
   - Create Copper rectangle (0, 0) → (0.1, 0.01)
   - Create FR4 rectangle (0, 0.01) → (0.1, 0.03)
   - Verify visualization shows correct colors
   - Edit coordinates, verify changes apply
   - Change material to Silver, verify conductivity updates

3. **Setup Week:**
   - Change grid from 256×256 → 512×512
   - Verify cell size recalculated
   - Change boundary CPML → Mur, verify absorption changes
   - Adjust CFL to 0.5, verify timestep changes
   - "Apply & Restart" without crash

### Regression Testing:
- [ ] Phase 2 features still work (material library, S-params, Smith chart)
- [ ] Paint mode still works alongside rectangles
- [ ] Keyboard shortcuts unchanged
- [ ] Performance maintained (60+ FPS)

---

## Risk Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Expression editor breaks existing sources | Low | High | Test with CW/Gaussian first, expr last |
| Material blocks conflict with paint | Medium | Medium | Keep separate arrays, paint → epsilon_map, blocks → material_rects |
| Grid restart loses user work | High | High | Save current config before restart, offer "Save & Restart" |
| UI clutter from added panels | Medium | Low | Use collapsing headers, default collapsed |

---

## Code Architecture Notes

### Existing Functions to Reuse:
```cpp
// Sources (lines 1330-1418)
wizard_sources_tab(WizardState& w)
  - Expression compiler integration (lines 1371-1415)
  - Live preview plot logic

// Materials (lines 1296-1328)
wizard_materials_tab(WizardState& w)
  - MaterialRectSpec editing
  - Type selector (Dielectric/PEC/PMC)

// Grid (lines 1270-1294)
wizard_grid_tab(WizardState& w)
  - Grid/domain inputs
  - CFL slider
  - Boundary combo

// Application functions (already exist)
apply_wizard_materials_to_sim() // line 1526
sync_sim_to_wizard_config()     // line 1545
bootstrap_sim_from_config()     // in app_bootstrap
```

### New Functions to Create:
```cpp
// Week 1
draw_source_editor(SourceConfigSpec& src, int index)
draw_expression_editor(SourceConfigSpec& src)
create_new_source(WizardState& w)
delete_source(WizardState& w, int index)

// Week 2
draw_material_block_editor(MaterialRectSpec& rect, int index)
create_material_block_from_library(WizardState& w, int material_id)
update_block_material(MaterialRectSpec& rect, int material_id)

// Week 3
draw_simulation_setup_panel(WizardState& w, AppState& app)
validate_grid_parameters(int nx, int ny)
calculate_memory_estimate(int nx, int ny)
restart_simulation_with_config(WizardState& w)
```

---

## Documentation Updates Required

After Phase 2.5 completion:

1. **User Guide Additions:**
   - "Creating Custom Sources with Expressions"
   - "Precise Material Placement with Blocks"
   - "Changing Simulation Grid Size"

2. **Tutorial Updates:**
   - Replace wizard screenshots with new panel screenshots
   - Update quickstart guide

3. **Progress Document Updates:**
   - Mark STRATEGIC_VISION_SUMMARY.md Phase 2.5 complete
   - Update IMGUI_IMPLEMENTATION_PROGRESS.md

4. **API Documentation:**
   - Document source expression syntax
   - List available expression variables (t, amp, freq, pi)
   - Example expressions library

---

## Benefits Summary

### For Users:
- ✅ **More discoverable** - controls in logical panels, not hidden wizard
- ✅ **More powerful** - expression sources, precise blocks, grid experimentation
- ✅ **More efficient** - no need to open separate wizard window
- ✅ **More professional** - modern UI paradigm

### For Development:
- ✅ **Code reuse** - 80% of logic already exists
- ✅ **Maintainability** - distributed UI easier to extend
- ✅ **Consistency** - matches Phase 2 UI patterns
- ✅ **Foundation** - enables Phase 3 automation

### For Education:
- ✅ **Expression editor** - unique teaching tool (visualize waveforms)
- ✅ **Grid experimentation** - learn convergence/stability
- ✅ **Material precision** - understand boundary effects

---

*Created: 2025-11-19*
*Status: Ready for Prompt #24 (Week 1)*
