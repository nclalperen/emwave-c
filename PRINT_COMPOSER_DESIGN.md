# Print Composer / Animation Export Draft

## Goals
- Provide a Matlab-style composer window to assemble exports (animations or static captures) without recording the live UI chrome.
- Let users choose sources (field viewports, regions, scopes, FFT, legends, measurements/annotations) and lay them out across one or more pages.
- Offer per-page resolution/format/output, per-item styling, and clean off-screen rendering piped directly to FFmpeg (with PNG-sequence fallback).

## User Flow (High-Level)
1) Open **Print Composer** window.
2) Choose a template (e.g., Field+Legend, Field+Scope, Field+FFT+Legend) or start **Blank**.
3) Add sources from the palette (field viewport, region capture, scope, FFT, legend, measurements/annotations, Smith/S-params). Each add spawns a placeholder on the active page.
4) Arrange placeholders on the page canvas (drag/resize, snap/align). Select an item to tweak binding (which viewport/channel/region), overlays, style, and render scale.
5) Set page properties (resolution, background, format override). Global defaults remain in the top bar.
6) Preview the composed output (play/pause/seek). Optionally hide composer controls during preview/export.
7) Export: stream frames off-screen to FFmpeg (or PNG sequence fallback). Per-page outputs; combined video optional/low priority.

## Window Layout
- **Top bar**: global defaults (format GIF/MP4/PNG seq; FPS; duration; output dir/name), play/pause/seek, memory/time estimate, “Hide UI during export” toggle, countdown toggle.
- **Page area**: main canvas (zoom/pan), optional thumbnail strip (collapsible) for quick page switching; page tabs with right-click for per-page settings.
- **Left palette**: source list with add buttons: Field viewport, Region capture (drag from live viewport), Scope (time), FFT, Legend, Measurements/Annotations, Smith chart, S-params plot, Cursor/Ruler overlay.
- **Right properties** (contextual): for selected item—binding (viewport/channel/region), overlays on/off (legend, measurements, annotations, cursor/ruler), style (bg/border/padding/opacity), render scale, fit-to-page/width.

## Templates
- Presets: Field+Legend, Field+Scope, Field+FFT+Legend, Blank.
- User can still rearrange/resize after loading a template.

## Page & Item Settings
- Page (per tab right-click): resolution (px), background color, format override (mp4/gif/png seq), margins/grid snap toggle, duplicate/delete.
- Item: data binding (which viewport/channel or which region capture), overlay toggles, style, render scale (e.g., 1.25x for sharper downscale), fit-to-page/width buttons, align tools (left/right/top/bottom/center/distribute).
- Region capture: “Add region from viewport…” lets user drag a rectangle on a live viewport; becomes a reusable source.

## Preview
- Shows the composed page only (no live UI chrome). Page selector, play/pause/seek, optional “preview at output resolution.”
- Countdown before export (e.g., 3s) if enabled.

## Export Pipeline (Technical)
- Off-screen render targets per item at requested resolution; composite into the page texture per frame.
- Dynamic layers re-render per frame (field, scope, FFT); static layers cached (legend, unchanging annotations).
- Stream frames directly to FFmpeg via spawn/argv (no shell; avoid %04d issues). PNG sequence fallback if FFmpeg rc != 0.
- Log last FFmpeg command/rc to `recordings/ffmpeg_last.txt`. Hide composer controls during final export if toggle is on.
- Per-page outputs; combined multi-page video optional and off by default.

## Static Snapshots (outside composer)
- When sim paused (spacebar), left-click on viewport/scope/FFT/legend opens a context menu: Export PNG/SVG/Copy. This provides quick static captures without entering the composer.

## Safety & UX
- Snap/grid + align tools; optional guides. Fit-to-page/width per item.
- Estimated memory/time shown; warn if caching too large; auto-stream to FFmpeg to avoid RAM bloat.
- File naming: label + page suffix (e.g., `_page1`, `_page2`); per-page override via tab context menu.

## Layout JSON (headless import)
- `--ch-layout=path.json` loads a saved layout before applying CLI overrides (format/fps/frames/output/resolution/template).
- Expected shape (permissive parser): top-level `"pages"` array. Page keys: `name`, `res_w`, `res_h`, `output_name`, `output_format` (0=BMP,1=PNG seq,2=MP4,3=GIF), `fps`, `frames`, `video_kbps`, `transparent_bg` (0/1), and `bg` as `[r,g,b,a]` or `bg_r/bg_g/bg_b/bg_a`.
- Each page may contain `"items"`: keys include `id` (optional), `type` (string: `field|legend|scope|fft|region|meas|smith` or numeric enum), `viewport_idx` (or `viewport`/`vp`), `pos` `[x,y]`, `size` `[w,h]` (or `w`,`h`), and `region_norm` `[x0,y0,x1,y1]` for region captures. Unknown types default to `field`. Missing ids auto-increment.
- Example:
  ```json
  {
    "pages": [
      {
        "name": "Page 1",
        "res_w": 1280,
        "res_h": 720,
        "output_name": "composer_page_1",
        "output_format": 1,
        "fps": 30,
        "frames": 60,
        "video_kbps": 4000,
        "transparent_bg": 0,
        "bg": [0.05, 0.05, 0.06, 1.0],
        "items": [
          { "type": "field", "viewport": 0, "pos": [80,80], "size": [720,480] },
          { "type": "legend", "pos": [840,80], "size": [280,180] }
        ]
      }
    ]
  }
  ```

## Open Questions / Decisions
- Thumbnail strip: horizontal (collapsible) vs vertical—default horizontal; switchable in settings.
- Do we support PDF/static page export now or later? (nice-to-have for legends/plots)
- Default items to add when opening (e.g., active field + legend + scope) or start blank? 
- Combined multi-page video: keep as low-priority toggle or omit initially.
