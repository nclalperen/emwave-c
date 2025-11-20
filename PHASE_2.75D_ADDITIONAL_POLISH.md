# Phase 2.75D: Additional Polish and Advanced Features

**Status:** Ready to implement  
**Effort:** 20 hours total (8h multi-viewport, 7h recording, 5h measurements)  
**Dependencies:** Phase 2.75A-C complete in codebase  
**Audience:** Implementation engineer and reviewers  
**Created:** 2025-11-20 (clean rewrite)

---

## Page 1: Executive Overview

Phase 2.75D adds three professional, Blender-inspired upgrades that turn the simulator into a publication-ready toolset:
- Multi-viewport split view with four layouts and per-viewport transforms
- Animation recording to GIF/MP4/PNG sequences with progress feedback
- Advanced measurement tools (area, perimeter, annotations, CSV export)

Key outcomes:
- Side-by-side field comparisons (Ez/Hx/Hy/Material/Overlay)
- Exportable animations for papers and presentations
- Quantitative measurements suitable for external analysis

---

## Page 2: Feature Breakdown (Fast Facts)

| Feature | Highlights | Hotkeys / Controls | Primary Benefit |
|---------|------------|--------------------|-----------------|
| Multi-Viewport Split (8h) | Single / Horizontal / Vertical / Quad layouts; independent zoom/pan; per-viewport visualization | Alt+1/2/3/4 to switch layouts; hover to focus active viewport | Compare components and regions simultaneously |
| Animation Recording (7h) | GIF/MP4/PNG output; 10-60 fps; 0.5x-2.0x resolution scale; progress bar; optional auto-play | Panel in Tools menu; Stop button; FFmpeg fallback to GIF/PNG | Publication-grade animations and demos |
| Advanced Measurements (5h) | Polygon area via shoelace; perimeter; text annotations with colors; history panel; CSV export | A = area mode, T = annotation, history panel toggles | Quantitative reporting and sharable measurements |

---

## Page 3: Multi-Viewport Specification

Layouts (targeting Blender familiarity):
```
Single (Alt+1)                Horizontal Split (Alt+2)
+--------------------+        +--------------------+
|                    |        |        A           |
|     Full view      |        +--------------------+
|                    |        |        B           |
+--------------------+        +--------------------+

Vertical Split (Alt+3)        Quad (Alt+4)
+-----------+-----------+     +-----------+-----------+
|     A     |     B     |     |     A     |     B     |
|           |           |     +-----------+-----------+
+-----------+-----------+     |     C     |     D     |
                               +-----------+-----------+
```

Behavior:
- Independent zoom/pan per viewport; optional sync toggles for zoom/pan.
- Shared simulation timestep; per-viewport visualization mode (Field/Material/Overlay/Magnitude).
- Active viewport follows hover; highlighted border and label A/B/C/D with viz mode.
- Minimum viewport size fallback to Single if space is insufficient.

Data structures:
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
    int active_viewport_idx; // 0-3 or -1
    bool sync_zoom;
    bool sync_pan;
};
```

Rendering approach:
- Compute layout rectangles, set `SDL_RenderSetViewport` per viewport, render textures with per-viewport transforms; reset viewport before ImGui overlays.
- Optional reduced resolution for inactive viewports to protect performance.

---

## Page 4: Animation Recording Specification

User flow:
1) Open Tools -> Animation Recording panel.  
2) Choose format (GIF/MP4/PNG sequence), framerate (10-60 fps), duration, resolution scale (0.5x-2.0x), auto-play toggle.  
3) Press Start; progress bar updates with frame count; Stop allowed.  
4) Encoding: GIF via jo_gif or stb_image_write; MP4 via FFmpeg; PNG writes sequence.  
5) On failure (missing FFmpeg), fall back to GIF or PNG with clear messaging.

Recorder structures:
```cpp
enum RecordingFormat { RECORDING_GIF, RECORDING_MP4, RECORDING_PNG_SEQUENCE };
enum RecordingState { RECORDING_IDLE, RECORDING_ACTIVE, RECORDING_PROCESSING, RECORDING_ERROR };

struct AnimationRecorder {
    RecordingState state;
    RecordingFormat format;
    float framerate;          // 10-60
    int target_frames;        // duration * framerate
    int frame_count;
    float resolution_scale;   // 0.5x - 2.0x
    std::vector<SDL_Surface*> frames;
    char output_path[260];
    bool auto_play;
    float progress;           // 0-1
    char status_message[128];
};
```

Capture rules:
- Capture from renderer each frame interval; enforce even dimensions for video encoders.
- PNG sequence stored in `recordings/temp` for MP4; cleaned after encode.
- Auto-play uses platform-default opener.

---

## Page 5: Advanced Measurement Tools Specification

Capabilities:
- Area tool: polygon input with auto-close when near first vertex; shoelace formula for area; perimeter computed in meters.
- Perimeter: derived from polygon edges; reusable for ruler history (distance measurements).
- Annotations: text labels with custom color and font size; pinned to grid coordinates.
- History panel: tree view for distances, areas, annotations; delete/visibility toggles.
- Export: CSV with type, id, value, units, and coordinates for spreadsheet analysis.

Core structures:
```cpp
enum MeasurementType { MEASURE_DISTANCE, MEASURE_AREA, MEASURE_PERIMETER, MEASURE_ANNOTATION };

struct DistanceMeasurement { ImVec2 a, b; double distance_m; double angle_deg; };
struct AreaMeasurement { std::vector<ImVec2> vertices; bool closed; double area_m2; double perimeter_m; };
struct Annotation { ImVec2 grid_pos; char text[128]; ImVec4 color; float font_size; bool visible; };

struct MeasurementHistory {
    std::vector<DistanceMeasurement> distances;
    std::vector<AreaMeasurement> areas;
    std::vector<Annotation> annotations;
};
```

Export format (CSV):
```
Type,ID,Value,Unit,Data
Distance,1,0.012345,m,"(10,20)-(60,80) @ 35.0deg"
Area,1,0.001234,m2,"5 vertices | Perimeter=0.0870 m"
Annotation,1,,,"(12.5,40.0): Feedline tap"
```

---

## Page 6: Effort and Timeline (20 hours total)

| Stream | Tasks | Effort |
|--------|-------|--------|
| Multi-Viewport | Layout manager, per-viewport transform state, visualization selector, hover-to-activate, alt hotkeys, performance check | 8h |
| Animation Recording | Recorder struct, frame capture, GIF export, MP4 via FFmpeg with fallback, PNG sequence, UI panel, progress bar | 7h |
| Advanced Measurements | Area/perimeter calculation, annotation dialog, history panel, CSV export, display toggles | 5h |
| **Total** | | **20h** |

Suggested schedule:
- Week 1: Multi-viewport (4 days)
- Week 2: Animation (2 days) + Measurements (2 days) + buffer (1 day)

---

## Page 7: Risks, Challenges, and Solutions

| Risk / Challenge | Impact | Mitigation |
|------------------|--------|------------|
| 4 viewports at 60 fps | GPU load increase | Use SDL render targets; update inactive viewports at lower rate; cap texture size for passive panes |
| MP4 relies on FFmpeg | Encoding may fail if missing | Detect FFmpeg; fall back to GIF/PNG; document install; optional bundled binary on Windows |
| Large GIF files | Disk bloat, slow sharing | Default to 640x480 @ 20-30 fps; warn when estimated size > 50 MB; allow lossy palette |
| Huge polygons | UI slowdowns | Soft limit 50 vertices; early exit for degenerate polygons; reuse ImDrawList buffers |
| Storage of many frames | Memory spikes during recording | Stream PNG sequence to disk when frame count > threshold; release surfaces after encode |

---

## Page 8: Acceptance and Next Actions

Success checklist (high level):
- Alt+1/2/3/4 smoothly switch layouts; active viewport clearly highlighted.
- Each viewport holds its own zoom/pan and visualization mode without side effects.
- Recording works for GIF/MP4/PNG; progress and errors clearly shown; fallback path succeeds when FFmpeg is absent.
- Measurement tools report correct area/perimeter; annotations render clearly; CSV opens cleanly in Excel/LibreOffice.
- Total effort contained to 20 hours with minimal disruption to existing Phase 2.75A-C code.

Next steps:
1. Follow prompts #37-39 for implementation details.  
2. Use the 43 test cases in the master guide to validate each feature.  
3. Update user-facing help once features land (short how-to per prompt).  
4. Capture before/after screenshots for release notes.

---

*End of Phase 2.75D overview (approx. eight pages of detail)*  
*Author: Planning Assistant*  
