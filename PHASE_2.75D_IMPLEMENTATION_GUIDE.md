# Phase 2.75D Implementation Guide

**Scope:** Multi-Viewport Split, Animation Recording, Advanced Measurement Tools  
**Effort:** 20 hours total | 3 prompts | 43 test cases  
**Status:** Ready for implementation (Phase 2.75A-C already in codebase)  
**Created:** 2025-11-20 (clean rewrite)

---

## Prompt Summary

| Prompt | Feature | Effort | Priority | Test Cases | Status |
|--------|---------|--------|----------|------------|--------|
| #37 | Multi-Viewport Split View | 8h | High | 16 | Ready |
| #38 | Animation Recording | 7h | Medium | 12 | Ready |
| #39 | Advanced Measurement Tools | 5h | Medium | 15 | Ready |
| **Total** | | **20h** | | **43** | |

Notes:
- All three prompts assume the Phase 2.75A-C navigation, grid, and measurement baselines are present.
- Prompts live in `prompts/PROMPT_37_Multi_Viewport_Split.md`, `prompts/PROMPT_38_Animation_Recording.md`, `prompts/PROMPT_39_Advanced_Measurement_Tools.md`.

---

## Implementation Options

**Option A: Sequential (Recommended)**
- Week 1: Prompt #37 (layouts, per-viewport transforms, hotkeys)
- Week 2: Prompt #38 (recording) then #39 (measurements)
- Pros: Simplest integration/testing path; each prompt stabilizes before adding the next.

**Option B: Parallel (Advanced)**
- Developer 1: Prompt #37
- Developer 2: Prompts #38 and #39
- Pros: Faster calendar time; requires early API alignment on AppState additions and rendering hooks.

**Option C: Essential Only**
- Implement Prompt #37 only (8h)
- Pros: Biggest UX impact with minimal scope; recording/measurements can follow later.

---

## How to Use These Prompts

1) Read the prompt document fully (objective, structures, sample code).  
2) Align AppState changes across prompts before coding.  
3) Implement and unit-check each sub-section in the prompt.  
4) Run the matching checklist (16/12/15 cases) and fix regressions immediately.  
5) Log user-facing changes for the release notes and quick-start docs.  

---

## Files You Will Touch

- `src/app/main_imgui.cpp` (AppState, event handling, rendering, UI panels)
- `src/app` or `src/render` helpers for viewport layout, capture, export
- `configs/` or `dist/` assets if bundling `jo_gif.h` or ffmpeg docs (optional)
- New lightweight helpers for CSV export (prompt #39) if not already present

Keep edits local to existing patterns established in Phase 2.75A-C.

---

## Testing Checklists (43 cases total)

| Prompt | Focus | Cases |
|--------|-------|-------|
| #37 | Layout switching, per-viewport zoom/pan, visualization modes, active highlighting, performance | 16 |
| #38 | Start/stop recording, framerate/resolution controls, GIF/MP4/PNG exports, progress UI, error/fallback handling | 12 |
| #39 | Area/perimeter accuracy, annotation workflow, history panel, CSV export, toggles and deletion | 15 |

Recommended order for manual testing:
1. Prompt #37 cases (single -> quad, hover activation, sync toggles, perf smoke test).  
2. Prompt #38 cases (record at 30 fps, MP4 with/without FFmpeg, GIF size under 50 MB, PNG sequence paths).  
3. Prompt #39 cases (shoelace correctness on rectangle/triangle, annotation color/size, CSV import in Excel/LibreOffice).  

---

## Success Criteria (Per Prompt)

**Prompt #37: Multi-Viewport Split**
- Alt+1/2/3/4 switches layouts instantly.  
- Each viewport keeps independent zoom/pan unless sync is enabled.  
- Visualization mode per viewport (Field/Material/Overlay/Magnitude).  
- Active viewport clearly outlined; hover selects active.  
- Minimum size fallback to single view avoids unusable panes.  
- No frame-rate regression beyond acceptable limits.

**Prompt #38: Animation Recording**
- Start/Stop controls record the active viewport(s) at target fps.  
- GIF export succeeds for a 10s clip at 30 fps with expected size.  
- MP4 export works when FFmpeg is available; clean fallback when absent.  
- PNG sequence exports numbered frames to disk.  
- Progress bar and status text reflect capture and encode phases.  
- Auto-play opens the resulting file when enabled.

**Prompt #39: Advanced Measurement Tools**
- Area tool closes polygons and reports correct area/perimeter.  
- Perimeter displayed in meters using grid spacing.  
- Text annotations render at correct grid coordinates with chosen color/size.  
- History panel lists distances, areas, annotations with delete/visibility toggles.  
- CSV export opens cleanly in spreadsheet tools with correct headers/values.  
- Performance remains smooth with 50+ stored measurements.

---

## Known Issues and Solutions

| Issue | Likelihood | Mitigation |
|-------|------------|------------|
| Four viewports drop fps on older GPUs | Medium | Render inactive panes at reduced resolution; throttle redraws when transforms unchanged |
| FFmpeg missing on user machines | High | Detect presence; fall back to GIF/PNG export; document install path; optionally bundle binary on Windows |
| Large GIF output (>50 MB) | Medium | Default to 640x480 @ 20-30 fps; warn when estimated size is high; allow smaller palettes |
| Annotation overlap at high zoom | Low | Provide hide/show toggles and small padding; keep font size adjustable |
| CSV locale differences (decimal separator) | Low | Use dot decimal and simple ASCII headers; user can import with delimiter settings |

---

## Progress Tracking Template

- [ ] Prompt #37 complete (16/16 tests)  
- [ ] Prompt #38 complete (12/12 tests)  
- [ ] Prompt #39 complete (15/15 tests)  
- [ ] Regression pass on existing Phase 2.75A-C behaviors (zoom, pan, dockspace, ruler)  

---

## Ready-to-Implement Checklist

- Prompts #37-39 understood and AppState diffs reconciled.  
- Build remains unchanged (no new libraries required except optional FFmpeg/jo_gif header).  
- UI strings and tooltips consistent with existing tone.  
- Release notes placeholders added for Phase 2.75D.  

*Author: Planning Assistant*  
