# Phase 2.75D Quick Reference Card

**Status:** ✅ Implementation verified in code (2025-11-21)
**Ready for:** Smoke testing and deployment

---

## Verified Hotkeys

### Multi-Viewport Layouts (Prompt #37)
| Hotkey | Layout | Code Location | Log Message |
|--------|--------|---------------|-------------|
| **Alt+1** | Single full viewport | main_imgui.cpp:3527 | "Layout: Single viewport" |
| **Alt+2** | Horizontal split (top/bottom) | main_imgui.cpp:3535 | "Layout: Horizontal split" |
| **Alt+3** | Vertical split (left/right) | main_imgui.cpp:3543 | "Layout: Vertical split" |
| **Alt+4** | Quad view (2×2 grid) | main_imgui.cpp:3551 | "Layout: Quad view" |

### Advanced Measurements (Prompt #39)
| Hotkey | Function | Code Location | Popup/UI |
|--------|----------|---------------|----------|
| **Shift+A** | Toggle area measurement mode | main_imgui.cpp:3907 | Log: "Area mode: ON/OFF" |
| **Shift+T** | Add text annotation | main_imgui.cpp:3863 | Opens "AddAnnotation" popup |
| **R** | Ruler tool (Phase 2.75C) | (existing) | Distance/angle readout |

### General Navigation (Phase 2.75A-C)
| Input | Function | Notes |
|-------|----------|-------|
| **Mouse wheel** | Zoom in/out | 0.5× to 32× range, zoom-to-cursor |
| **Middle mouse drag** | Pan viewport | Alternative: Shift+Right drag |
| **Shift+Arrow keys** | Keyboard pan | Step size adapts to zoom |

---

## Verified Menu Items

### Tools Menu
Located in main menu bar:

```
Tools
├── Measurement History...    (main_imgui.cpp:4712)
│   └── Shows: distances, areas, annotations
│       with delete/visibility controls
│
└── Animation Recording...    (main_imgui.cpp:4718)
    └── Controls: duration, fps, format (GIF/MP4/PNG)
        resolution scale, auto-play
```

### View Menu
```
View
└── Sync Zoom                 (for multi-viewport)
└── Sync Pan                  (for multi-viewport)
```

---

## Verified UI Panels

### 1. Animation Recording Panel
**Window Title:** "Animation Recording"
**Code:** main_imgui.cpp:416-500
**Open via:** Tools → Animation Recording...

**Controls:**
- Duration slider: 1-60 seconds
- Framerate slider: 10-60 fps
- Resolution scale: 0.5× to 2.0×
- Format combo: GIF / MP4 / PNG Sequence
- Auto-play checkbox
- Start/Stop buttons
- Progress bar + status message

**Output Location:** `recordings/emwave_YYYYMMDD_HHMMSS.{gif|mp4|png}`

### 2. Measurement History Panel
**Window Title:** "Measurement History"
**Code:** main_imgui.cpp:338-400
**Open via:** Tools → Measurement History...

**Sections:**
- Distances (from ruler tool)
- Areas (from Shift+A polygon tool)
- Annotations (from Shift+T text tool)

**Controls:**
- Show/hide toggles for each section
- Delete buttons per measurement
- Export to CSV button
- Clear All button

---

## Multi-Viewport Rendering

### Viewport Layout Manager
**Function:** `compute_viewport_layout()`
**Code:** main_imgui.cpp:560-625
**Minimum size:** 200×200 px per pane (auto-falls back to Single if too small)
**Gap between panes:** 4px

### Per-Viewport State
Each of 4 viewports (A, B, C, D) has independent:
- Zoom level (0.5× to 32×)
- Pan offset (x, y)
- Visualization mode: Field / Material / Overlay / Magnitude
- Grid visibility toggle
- Source overlay toggle

**Sync modes:** When enabled, zoom/pan changes apply to all viewports

### Viewport Labels
**Location:** Top-left corner of each viewport
**Format:** `"A: Field"`, `"B: Material"`, etc.
**Active highlight:** Blue outline (IM_COL32(80,160,255,255))
**Inactive:** Gray outline (IM_COL32(80,80,80,180))

---

## Animation Recording Details

### Supported Formats

**GIF Export**
- Library: jo_gif.h (single-header)
- Delay: 1000/framerate ms
- Function: `export_gif()` at main_imgui.cpp:895

**MP4 Export**
- Pipeline: PNG temp → FFmpeg subprocess → cleanup
- FFmpeg command: `-framerate <fps> -i temp/frame_%04d.png -c:v libx264 -pix_fmt yuv420p`
- Fallback: GIF/PNG if FFmpeg not found
- Function: `export_mp4()` at main_imgui.cpp:952

**PNG Sequence**
- Output: `<path>_frames/frame_0001.png`, `frame_0002.png`, ...
- Function: `export_png_sequence()` at main_imgui.cpp:838

### Capture Process
1. `start_recording()` → reserves frame buffer, sets target count
2. Render loop → `capture_frame_if_recording()` → `SDL_RenderReadPixels`
3. Progress updates: `rec->progress = frame_count / target_frames`
4. On completion → export to selected format → free frames

---

## Advanced Measurement Details

### Area Tool (Shoelace Formula)
**Activation:** Shift+A
**Workflow:**
1. Left-click to add vertices (minimum 3)
2. Right-click to close polygon (or click near first vertex)
3. Area computed: `|Σ(x_i * y_{i+1} - x_{i+1} * y_i)| / 2`
4. Perimeter: Sum of edge lengths

**Rendering:**
- Fill: Green semi-transparent (IM_COL32(80,180,80,80))
- Outline: Bright green (IM_COL32(100,255,100,200))
- Label at centroid: `"%.6f m²\n%.4f m"`

**Code:** `close_area_measurement()` at main_imgui.cpp:260-300

### Text Annotations
**Activation:** Shift+T
**Popup:** "AddAnnotation" (main_imgui.cpp:5090-5108)
**Fields:**
- Text input (128 chars max)
- Color picker (default: yellow 1,1,0.2,1)
- Font size slider (8-28 px, default 14)

**Rendering:**
- Black background with rounded corners
- Colored border matching text color
- Positioned at grid coordinates

### CSV Export
**Function:** `export_measurements_csv()` at main_imgui.cpp:303-330
**Format:**
```csv
Type,ID,Value,Unit,Data
Distance,1,0.125000,m,"(10.0,20.0)-(30.0,40.0) @ 45.0 deg"
Area,1,0.003125,m2,"4 vertices | Perimeter=0.2500 m"
Annotation,1,,,"(50.0,60.0): Test Note"
```

---

## Code Statistics

**Total implementation:** 1,134 new lines in main_imgui.cpp
**Data structures:** 8 new structs (ViewportInstance, AnimationRecorder, etc.)
**Functions added:** 12+ (layout manager, recording, export, measurements)
**UI panels:** 2 new (Recording, Measurement History)
**Hotkeys:** 4 layout + 2 measurement = 6 new hotkeys
**Menu items:** 2 new (Tools menu)

---

## Testing Status

**Build verification:** ✅ emwave_imgui.exe (932 KB, Nov 21 01:30)
**Launch test:** ✅ App starts without crash
**Code validation:** ✅ All hotkeys/menus/functions verified in source
**Next step:** 5-minute smoke test (see PHASE_2.75D_SMOKE_TEST.md)

---

## Dependencies

**Required (already present):**
- SDL2, SDL2_ttf
- ImGui 1.89+

**Optional (for MP4 export):**
- FFmpeg (runtime dependency, graceful fallback if absent)

**Single-header libraries (to bundle):**
- jo_gif.h (for GIF export)

---

*Verified: 2025-11-21*
*Author: Code Validation Agent*
