# Phase 2.75D: 5-Minute Smoke Test (ImGui)

Date: 2025-11-21  
Build: `build-imgui/Release/emwave_imgui.exe`  
Config for quick checks: `configs/waveguide.json`

Note: After pressing Clear, simulation sometimes stops rendering until restart or Apply Run Settings + restart (pending fix).

---

## How to Run
```powershell
cd C:\projects\emwave-c\build-imgui\Release
.\emwave_imgui.exe --config ..\..\configs\waveguide.json
```
Expected startup (example):
```
Grid: 512x256 (...) CFL=0.90
Sweep: ...
```

---

## Test 1 - Viewports & Transforms (2 min)
- [x] Alt+1/2/4 switch layouts (Single / Horizontal / Quad). Quad default: A=Ez, B=Hx, C=Hy, D=|S|.
- [x] Hover each pane: active shows blue outline; labels reflect channel (e.g., "A: Ez").
- [x] Scroll zoom per-pane; Alt+Left-drag pans. Verify only the active pane moves.
- [x] Enable Sync Zoom and Sync Pan (toolbar "Cfg"); zoom/pan one pane -> all panes mirror.
- [x] Toggle Grid/Sources/Vectors in "Cfg" for one pane; only that pane updates. Vectors appear on H/S channels.

Status: PASS (layouts, highlight, per-pane zoom/pan, sync zoom/pan, per-pane overlays OK).

---

## Test 2 - Recording (1.5 min)
- [x] Open Tools -> Animation Recording.
- [x] Set Duration=2s, FPS=10, Format=GIF (or MP4). Start -> status counts frames -> Processing -> Saved: `recordings\...`.
- [x] GIF/MP4 encode inside app. If FFmpeg is found (PATH/`FFMPEG_PATH`/`./ffmpeg/bin/ffmpeg.exe`), MP4+GIF use it; otherwise PNG sequence is produced (status shows failure code). Last FFmpeg command is logged to `recordings/ffmpeg_last.txt`.
- [x] Start a 5s record, click Stop after ~1s -> stops cleanly (no crash).

Status: PASS (frames visible; GIF/MP4/PNG produced in-app; PNG fallback when FFmpeg missing/fails).

---

## Test 3 - Measurements (1.5 min)
- [x] Area tool: Shift+A, click 3+ vertices, Right-click to close -> filled polygon with area/perimeter label.
- [x] Annotation: Shift+T -> popup -> add text -> shows at cursor. Note: Always visible; prefer hide unless hovered on area/context.
- [x] Measurement History: lists distances/areas/annotations; Clear Areas works; Export CSV writes `measurements.csv`. Note: CSV should be numbered (avoid overwrite).

Status: PASS (area, annotations, history, clear, CSV). Numbered CSV export implemented; annotation hover-only toggle added.

---

## Hotkeys Sanity (30 sec)
- [x] Ruler toggle R (measure two clicks, log entry). Request: brief popup showing the measured distance.
- [x] Layout hotkeys Alt+1/2/4 (Alt+3 optional; not critical).
- [x] Help F1 toggles overlay (optional). Background overlays are suppressed while F1 is shown.

Status: Ruler OK (popup added); Alt+1/2/4 OK; Alt+3 not used/critical. F1 overlay background cleaned.

---

## Pass/Fail Notes
- Pass if: no crashes; per-pane zoom/pan/channel work; recording produces visible output; measurements/annotations render; history/CSV OK.
- Current blockers: Clear may stall sim until restart.

---

## Quick Log Checklist
- App launch
- Layout switch
- Per-pane zoom/pan
- Sync zoom/pan
- Vectors
- Recording (visible output)
- Area tool
- Annotation
- History/CSV
