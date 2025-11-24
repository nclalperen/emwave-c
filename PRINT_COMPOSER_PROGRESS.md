# Print Composer Progress

## Done in this pass
- Field items: off-screen render (zoom/pan/viz/grid/sources/vectors/material overlay/outlines) into item-sized textures; falls back to placeholders on failure.
- Scope items: off-screen render with auto y-scale; FFT items: off-screen render of scope FFT (dB, 80 dB span).
- Legend items: off-screen render of the help legend; Measurements items: summary counts.
- Snap-to-grid with adjustable step; canvas grid overlay; multi-select with Ctrl+click; Select All/Clear; align (left/right/top/bottom/center X/Y) and distribute (X/Y) with snap.
- Page switching: explicit page list plus “+ Add Page” / Duplicate / Delete; templates (Blank, Field+Legend, Field+Scope, Field+FFT+Legend) with one-click apply.
- Region source: “Add Region” renders a viewport crop defined by normalized bounds; Region pick UX: “Pick Region”/“Re-pick” to drag a box on a live viewport (optional snap to cells), auto-adding/updating a Region item and sizing to aspect.
- Export all pages: one-click exports every page; per-page output name/format picker (BMP/PNG/MP4/GIF). PNG is real (embedded writer). MP4/GIF attempt ffmpeg rawvideo streaming, otherwise fall back to PNG with a log note.
- Page export writes to `recordings/<output_name>.<ext>`; animation export added (steps the simulation per frame; MP4/GIF via ffmpeg; PNG sequence fallback).
- Viewport overlays (sources/material blocks/region picker/paint ring) now render inside the viewport overlay window with proper clipping, preventing bleed into other UI windows.
- Animation export now clears/recreates the `_frames` folder before writing and uses the chosen container filename (mp4/gif) so ffmpeg always receives a valid output path with extension.

## Notes
- Field rendering uses the current simulation state (`sim->step_Ez_absmax` for scaling) and the viewport’s visualization mode. Scope rendering auto-computes y-scale from the current scope buffer. Render state is restored after export.

## Still pending vs. the design
- Preview player/controls, combined multi-page output, thumbnail strip, and a cleaner pipeline for multi-page animations.
