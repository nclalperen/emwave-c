# Prompt #24: Sources Panel Enhancement with Expression Editor

**Phase:** 2.5 (Wizard Integration - Week 1)
**Effort:** ~6 hours
**Priority:** HIGH
**Status:** Ready for implementation

---

## Objective

Transform the Sources panel from a basic display into a **full source editor** with create/delete/edit capabilities and an **expression-based source editor with live preview**. Extract the expression editor logic from the hidden wizard and integrate it into the modern UI.

---

## Background

**Current State:**
- Sources panel exists in [src/app/main_imgui.cpp](../src/app/main_imgui.cpp) (function: `draw_sources_panel()`)
- Shows list of existing sources (mostly read-only)
- Mouse placement works but limited
- **Expression editor exists but hidden** in `wizard_sources_tab()` (lines 1330-1418)
- No UI for create/delete sources
- No property editor for type, frequency, amplitude, field component

**Why This Matters:**
- **Expression-based sources are unique** - no other EM simulator has this capability
- **Educational value** - students can visualize custom waveforms before running
- **Research capability** - precise excitation for S-parameter characterization
- **UX improvement** - edit sources in context, not in separate wizard

---

## Requirements

### 1. Source List with Controls

**Current behavior:** Static list showing source positions
**Target behavior:** Interactive list with per-source controls

**UI Layout:**
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Sources                              ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

Source 0: CW @ 2.400 GHz
  Position: (0.25, 0.50)
  [✓ Active] [✎ Edit] [🗑 Delete]

Source 1: Gaussian Pulse
  Position: (0.75, 0.50)
  [  Inactive] [✎ Edit] [🗑 Delete]

[+ New Source]
```

**Implementation:**
- Iterate through `sim->source_count` sources
- For each source, show:
  - Source type label (CW/Gaussian/Ricker/Expr)
  - Frequency (for CW) or "Pulse" label
  - Position in normalized coordinates
  - **Active checkbox** - toggle `sources[i].active`
  - **Edit button** - opens property editor (below)
  - **Delete button** - confirms, then removes source

**New Source Button:**
- Adds new source to `sim->sources[]` array
- Increments `sim->source_count`
- Initializes with defaults:
  ```c
  sources[i].type = SRC_CW;
  sources[i].field = SRC_FIELD_EZ;
  sources[i].active = 1;
  sources[i].x = 0.5;
  sources[i].y = 0.5;
  sources[i].amp = 1.0;
  sources[i].freq = 2.4e9;  // 2.4 GHz
  sources[i].sigma2 = 1.0;
  sources[i].expr[0] = '\0';
  ```

---

### 2. Source Property Editor

**Behavior:** When "Edit" button clicked, expand editor inline OR open in popup window

**UI Layout:**
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Edit Source 0                        ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

Type: [CW ▼]
Options: CW | Gaussian | Ricker | Expression

Field Component: [Ez ▼]
Options: Ez | Hx | Hy

Position (normalized 0-1):
  X: [0.25    ]  (0.025 m)
  Y: [0.50    ]  (0.050 m)

Amplitude: [1.0     ]

Frequency: [2.400e9 ] Hz  (2.400 GHz)

Sigma² (pulse width): [1.0     ]
  (for Gaussian/Ricker types)

[Apply Changes] [Cancel]
```

**Implementation Details:**

1. **Type Selector:**
   ```cpp
   const char* type_items[] = {"CW", "Gaussian", "Ricker", "Expression"};
   int type_idx = (int)sources[i].type;
   if (ImGui::Combo("Type", &type_idx, type_items, 4)) {
       sources[i].type = (SourceType)type_idx;
   }
   ```

2. **Field Component:**
   ```cpp
   const char* field_items[] = {"Ez", "Hx", "Hy"};
   int field_idx = (int)sources[i].field;
   if (ImGui::Combo("Field", &field_idx, field_items, 3)) {
       sources[i].field = (SourceFieldType)field_idx;
   }
   ```

3. **Position Inputs:**
   ```cpp
   ImGui::InputDouble("X (0..1)", &sources[i].x, 0.01, 0.1);
   ImGui::InputDouble("Y (0..1)", &sources[i].y, 0.01, 0.1);

   // Show physical coordinates as hint
   double phys_x = sources[i].x * sim->lx;
   double phys_y = sources[i].y * sim->ly;
   ImGui::SameLine();
   ImGui::TextDisabled("(%.3f, %.3f m)", phys_x, phys_y);
   ```

4. **Frequency Input:**
   ```cpp
   double freq_ghz = sources[i].freq * 1e-9;
   if (ImGui::InputDouble("Frequency (GHz)", &freq_ghz, 0.1, 1.0, "%.3f")) {
       sources[i].freq = freq_ghz * 1e9;
   }
   ImGui::SameLine();
   ImGui::TextDisabled("(%.3e Hz)", sources[i].freq);
   ```

5. **Conditional Sigma² Input:**
   - Only show for Gaussian/Ricker types:
   ```cpp
   if (sources[i].type == SRC_GAUSS_PULSE || sources[i].type == SRC_RICKER) {
       ImGui::InputDouble("Sigma²", &sources[i].sigma2, 0.1, 1.0);
   }
   ```

---

### 3. Expression Editor (CRITICAL FEATURE)

**Behavior:** When source type = "Expression", show multi-line text editor with **live preview plot**

**Extract from:** `wizard_sources_tab()` lines 1371-1415 in [src/app/main_imgui.cpp](../src/app/main_imgui.cpp)

**UI Layout:**
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Expression Editor                    ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

Expression (C-like syntax):
┌────────────────────────────────────────┐
│ amp * sin(2*pi*freq*t)                 │
│   * exp(-(t - 1e-9)^2 / sigma2)        │
│                                        │
└────────────────────────────────────────┘

Available Variables:
  t      - Time (seconds)
  amp    - Amplitude
  freq   - Frequency (Hz)
  pi     - π (3.14159...)

Example Expressions:
  CW: amp * sin(2*pi*freq*t)
  Gaussian: amp * exp(-(t-1e-9)^2/sigma2)
  Modulated: amp * sin(2*pi*freq*t) * (1 + 0.5*sin(2*pi*1e8*t))

Preview (3 periods):
┌────────────────────────────────────────┐
│  1.0┤     ╱╲    ╱╲    ╱╲               │
│  0.0┤────╱──╲──╱──╲──╱──╲──            │
│ -1.0┤        ╲╱    ╲╱    ╲╱             │
│     └────────────────────> Time (ns)   │
└────────────────────────────────────────┘

✓ Expression compiled successfully

[Apply Expression]
```

**Implementation (extract from wizard code):**

```cpp
// Multi-line text input for expression
ImGui::InputTextMultiline("##expr",
                          sources[i].expr,
                          SOURCE_EXPR_MAX_LEN,
                          ImVec2(-1.0f, ImGui::GetTextLineHeight() * 4));

ImGui::TextDisabled("Variables: t (seconds), amp, freq, pi");

// Live preview compilation and plotting
if (sources[i].expr[0] != '\0') {
    char errbuf[128];
    ExprProgram* prog = nullptr;

    // Compile expression
    if (expr_compile(sources[i].expr, &prog, errbuf, sizeof(errbuf))) {
        const int N = 128;
        float values[N];

        // Calculate time range (3 periods)
        double period = (sources[i].freq > 0.0) ? (1.0 / sources[i].freq) : 1e-9;
        double t_max = 3.0 * period;
        if (t_max <= 0.0) t_max = 3e-9;

        // Evaluate expression over time range
        double vmin = 0.0, vmax = 0.0;
        bool first = true;
        for (int k = 0; k < N; ++k) {
            double t = t_max * (double)k / (double)(N - 1);
            double v = expr_eval(prog, t, sources[i].amp, sources[i].freq);
            values[k] = (float)v;

            if (first) {
                vmin = vmax = v;
                first = false;
            } else {
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
        }

        // Plot waveform
        float ymin = (float)vmin;
        float ymax = (float)vmax;
        if (ymin == ymax) {
            ymin -= 1.0f;
            ymax += 1.0f;
        }

        ImGui::PlotLines("Preview", values, N, 0, nullptr, ymin, ymax,
                        ImVec2(-1.0f, ImGui::GetTextLineHeight() * 6));

        ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1.0f),
                          "✓ Expression compiled successfully");

        expr_free(prog);
    } else {
        // Show compilation error
        ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.6f, 1.0f),
                          "Expression error: %s", errbuf);
    }
}
```

**Example Expressions to Document:**

1. **CW (continuous wave):**
   ```
   amp * sin(2*pi*freq*t)
   ```

2. **Gaussian pulse:**
   ```
   amp * exp(-(t - 1e-9)^2 / sigma2)
   ```

3. **Gaussian-modulated sine:**
   ```
   amp * sin(2*pi*freq*t) * exp(-(t - 1e-9)^2 / sigma2)
   ```

4. **Amplitude modulation:**
   ```
   amp * sin(2*pi*freq*t) * (1 + 0.5*sin(2*pi*1e8*t))
   ```

5. **Chirp (frequency sweep):**
   ```
   amp * sin(2*pi*(freq + 1e9*t)*t)
   ```

---

### 4. Source Deletion

**Behavior:** Confirm deletion, then remove source from array

**Implementation:**
```cpp
if (ImGui::Button("Delete")) {
    ImGui::OpenPopup("delete_confirm");
}

if (ImGui::BeginPopupModal("delete_confirm", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::Text("Delete Source %d?", i);
    ImGui::Separator();

    if (ImGui::Button("Yes, Delete", ImVec2(120, 0))) {
        // Shift sources down
        for (int j = i; j < sim->source_count - 1; ++j) {
            sim->sources[j] = sim->sources[j + 1];
        }
        sim->source_count--;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(120, 0))) {
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}
```

---

## Files to Modify

### Primary File: `src/app/main_imgui.cpp`

**Functions to Modify:**
1. **`draw_sources_panel()`** - Main sources panel (currently lines ~175-250)
   - Replace basic list with interactive list
   - Add New/Edit/Delete buttons
   - Add active/inactive toggles

**New Functions to Add:**
```cpp
// Source property editor
static void draw_source_editor(SimulationState* sim, int source_index, bool* open);

// Expression editor with preview
static void draw_expression_editor(SimulationState* sim, int source_index);

// Source creation
static void create_new_source(SimulationState* sim);

// Source deletion
static void delete_source(SimulationState* sim, int source_index);
```

**Code to Extract from Wizard:**
- Lines 1371-1415 in `wizard_sources_tab()` → expression editor logic
- Adapt to work with `sim->sources[i]` instead of `wizard.cfg.source_configs[i]`

---

## Testing Checklist

### Basic Functionality:
- [ ] Open Sources panel
- [ ] Click "New Source" - source created with defaults
- [ ] Source appears in list with "Active" checked
- [ ] Click "Edit" - property editor opens
- [ ] Change type from CW → Gaussian - sigma² field appears
- [ ] Change frequency to 5.0 GHz - updates correctly
- [ ] Uncheck "Active" - source stops emitting in simulation
- [ ] Click "Delete" - confirmation popup appears
- [ ] Confirm delete - source removed from list

### Expression Editor:
- [ ] Create new source, set type to "Expression"
- [ ] Expression editor appears with text box
- [ ] Enter: `amp * sin(2*pi*freq*t)`
- [ ] Preview plot shows sinusoid (3 periods)
- [ ] Change to: `amp * exp(-t/1e-9)` - plot shows exponential decay
- [ ] Enter invalid expression: `amp * (` - error message shows in red
- [ ] Enter valid expression again - "✓ Expression compiled successfully"
- [ ] Run simulation - custom waveform excites field

### Edge Cases:
- [ ] Create 10 sources - all editable independently
- [ ] Delete middle source (index 5) - array shifts correctly
- [ ] Delete all sources - simulation runs without excitation
- [ ] Create source with freq = 0 - preview shows 3 ns window (default)
- [ ] Expression with division by zero - error caught and displayed

### Integration:
- [ ] Sources created via UI visible in Scene panel
- [ ] Mouse-placed sources editable via panel
- [ ] Keyboard shortcuts (1-3 toggle) still work
- [ ] F5/F6 preset loads preserve new sources

---

## Success Criteria

**Phase 2.5 Week 1 complete when:**
1. ✅ Sources can be created/deleted via UI
2. ✅ All source properties editable (type, field, position, freq, amp, sigma²)
3. ✅ Expression editor functional with live preview plot
4. ✅ Active/inactive toggle works per source
5. ✅ No need to use wizard for source configuration
6. ✅ Expression compilation errors displayed clearly
7. ✅ Example expressions documented in help text

---

## Implementation Notes

### AppState Additions (if needed):
```cpp
struct AppState {
    // ... existing fields ...

    // Source editor state
    int editing_source_index;  // -1 = none, 0+ = editing source i
    bool show_source_editor;   // Editor popup open flag
    bool show_expr_help;       // Expression help overlay
};
```

### Expression Syntax Reference:
**Operators:** `+`, `-`, `*`, `/`, `^` (power)
**Functions:** `sin()`, `cos()`, `exp()`, `sqrt()`, `abs()`
**Variables:** `t`, `amp`, `freq`, `pi`
**Constants:** Numeric literals (e.g., `1.0`, `2e-9`, `3.14`)

### Performance Considerations:
- Expression compilation only on text change, not every frame
- Preview plot limited to 128 samples (fast evaluation)
- Cache compiled ExprProgram* for active expressions

### Error Handling:
- Invalid expression → show error in red, don't crash
- Division by zero → clamp to ±1e30, show warning
- Empty expression → treat as `0.0` (no excitation)

---

## Documentation Updates

After implementation, add to user guide:

**Section: "Creating Custom Source Waveforms"**

Content:
1. How to create a new source
2. Source type descriptions (CW/Gaussian/Ricker/Expression)
3. Expression syntax guide
4. Example expressions with screenshots of preview plots
5. Tips for debugging expression errors

---

## Expected Outcome

**Before:** Users had to use hidden wizard or edit JSON to configure sources
**After:** Full source editor in Sources panel with unique expression capability

**Screenshot Goal:**
```
╔═══════════════════════════════════════════════╗
║ Sources Panel                                 ║
╠═══════════════════════════════════════════════╣
║ Source 0: Custom Chirp                        ║
║   [✓ Active] [Edit ▼]                         ║
║                                               ║
║   Type: Expression                            ║
║   Position: (0.25, 0.50)                      ║
║   Amplitude: 1.0                              ║
║   Frequency: 2.4 GHz                          ║
║                                               ║
║   Expression:                                 ║
║   ┌─────────────────────────────────────────┐ ║
║   │ amp*sin(2*pi*(freq+5e8*t)*t)            │ ║
║   └─────────────────────────────────────────┘ ║
║                                               ║
║   Preview:                                    ║
║   [Chirp waveform plot with increasing freq] ║
║   ✓ Expression compiled successfully          ║
║                                               ║
║   [Apply] [Cancel] [Delete]                   ║
║                                               ║
║ [+ New Source]                                ║
╚═══════════════════════════════════════════════╝
```

---

**Ready for implementation!** 🚀

Hand this prompt to your implementation agent and proceed with Week 1.
