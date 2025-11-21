# Phase 2.75D: 5-Minute Smoke Test

**Date:** 2025-11-21
**Build:** emwave_imgui.exe (Nov 21 01:30, 932 KB)
**Config:** configs/waveguide.json
**Objective:** Verify Phase 2.75D features are functional before committing code

---

## How to Run

```bash
cd c:/projects/emwave-c/build-imgui/Release
./emwave_imgui.exe --config ../../configs/waveguide.json
```

**Expected startup output:**
```
Grid: 512x256 (0.50m x 0.25m) CFL=0.90
Sweep: 3 pts, 6.000e+09 Hz to 8.500e+09 Hz, 400 steps/pt
Run mode: fixed (1000 steps)
Boundary: CPML
Materials: 3 rectangles, Sources: 1, Ports: 0
```

✅ **VERIFIED:** App launches without crash

---

## Test 1: Multi-Viewport Layouts (Prompt #37) — 2 minutes

### 1.1 Layout Switching (30 seconds)
- [ ] Press **Alt+1** → Single full viewport appears
- [ ] Press **Alt+2** → Window splits horizontally (viewport A top, B bottom)
- [ ] Press **Alt+3** → Window splits vertically (viewport A left, B right)
- [ ] Press **Alt+4** → Quad view (A/B top, C/D bottom)
- [ ] Check log: Each keypress should log "Layout: Single/Horizontal/Vertical/Quad"

**Expected Result:** Layouts switch instantly, no flickering or crashes

### 1.2 Active Viewport Highlighting (30 seconds)
- [ ] In Quad view (Alt+4), hover mouse over each viewport
- [ ] Verify active viewport shows **blue outline** (IM_COL32(80,160,255,255))
- [ ] Inactive viewports show **gray outline** (IM_COL32(80,80,80,180))
- [ ] Top-left corner of each viewport shows label: "A: Field", "B: Field", etc.

**Expected Result:** Active viewport is clearly indicated

### 1.3 Independent Zoom/Pan (1 minute)
- [ ] In Quad view, click on viewport A
- [ ] Scroll mouse wheel → zoom changes in A only
- [ ] Middle-mouse drag → pan changes in A only
- [ ] Click viewport B, zoom → verify A remains unchanged
- [ ] Open View menu → Check "Sync Zoom" ✅
- [ ] Zoom in any viewport → all viewports zoom together

**Expected Result:** Independent transforms work; sync mode works

---

## Test 2: Animation Recording (Prompt #38) — 1.5 minutes

### 2.1 Recording Panel Opens (15 seconds)
- [ ] Click **Tools → Animation Recording...** (or equivalent menu)
- [ ] Panel opens with controls:
  - Duration slider (1-60s)
  - Framerate slider (10-60 fps)
  - Resolution scale slider (0.5×-2.0×)
  - Format combo (GIF/MP4/PNG Sequence)
  - Auto-play checkbox
  - Start/Stop buttons
  - Progress bar + status text

**Expected Result:** Panel appears without errors

### 2.2 Quick Recording Test (1 minute)
- [ ] Set Duration = 2 seconds
- [ ] Set Framerate = 10 fps
- [ ] Set Format = GIF
- [ ] Click **Start Recording**
- [ ] Status shows "Recording: 0/20 frames" (2s × 10fps = 20 frames)
- [ ] Progress bar fills to 100%
- [ ] Status changes to "Processing..."
- [ ] After 1-2 seconds, status shows "Complete: recordings/emwave_YYYYMMDD_HHMMSS.gif"
- [ ] Check that file exists in recordings/ folder

**Expected Result:** GIF created successfully, no crashes

### 2.3 Stop Button (15 seconds)
- [ ] Start another recording (5s @ 30fps)
- [ ] After 1 second, click **Stop Recording**
- [ ] Recording halts early
- [ ] Partial frames saved or gracefully discarded

**Expected Result:** Stop works without crash

---

## Test 3: Advanced Measurements (Prompt #39) — 1.5 minutes

### 3.1 Area Tool (45 seconds)
- [ ] Press **Shift+A** key
- [ ] Log shows "Area mode: ON"
- [ ] Left-click 3 times in viewport to create triangle vertices
- [ ] Right-click to close polygon
- [ ] Polygon fills with green semi-transparent color
- [ ] Polygon outline appears in bright green
- [ ] Label at centroid shows area (m²) and perimeter (m)

**Expected Result:** Area calculation displays correctly

### 3.2 Text Annotations (30 seconds)
- [ ] Press **Shift+T** key
- [ ] Popup "AddAnnotation" opens with:
  - Text input field
  - Color picker (default yellow)
  - Size slider (8-28 px, default 14)
  - Add / Cancel buttons
- [ ] Type "Test Note" in text field
- [ ] Click **Add**
- [ ] Yellow text "Test Note" appears at cursor grid position with black background

**Expected Result:** Annotation created and visible

### 3.3 Measurement History Panel (15 seconds)
- [ ] Open **Tools → Measurement History** (or equivalent)
- [ ] Panel lists:
  - Distances section (if ruler tool was used)
  - Areas section (shows triangle from step 3.1)
  - Annotations section (shows "Test Note")
- [ ] Each entry has delete button
- [ ] Checkboxes: Show Distances, Show Areas, Show Annotations
- [ ] Buttons: Export to CSV, Clear All

**Expected Result:** History panel displays measurements correctly

---

## Critical Blockers (Stop & Report If Found)

1. **Crash on layout switch** → STOP, report stack trace
2. **Recording fails silently** → Check for error message, verify FFmpeg presence
3. **Area tool crashes on polygon close** → STOP, report error
4. **Annotations don't appear** → Check if text is visible/rendered

---

## Pass Criteria

- ✅ 0/15 critical tests pass → **All systems functional**
- ⚠️ 12-14/15 pass → **Minor issues, proceed with caution**
- ❌ <12/15 pass → **Major regressions, do NOT commit**

---

## Post-Test Actions

### If All Tests Pass:
1. Close the application
2. Review this checklist for any warnings
3. Proceed to commit Phase 2.75D implementation

### If Tests Fail:
1. Note which tests failed
2. Check logs for error messages
3. Report issues before committing

---

## Time Budget

- Test 1 (Multi-Viewport): 2 min
- Test 2 (Recording): 1.5 min
- Test 3 (Measurements): 1.5 min
- **Total:** 5 minutes

---

## Notes

- Application launched successfully on 2025-11-21 01:30
- Config loaded: waveguide.json (512×256 grid)
- Build size: 932 KB
- No compilation errors reported
- Phase 2.75D implementation: 1,134 lines added to main_imgui.cpp

*This smoke test validates critical functionality only. Detailed testing will occur during real usage (bug fixing stage).*
