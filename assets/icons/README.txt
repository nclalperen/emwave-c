Rasterize the SVG toolbar icons before use.

If ImageMagick is installed:
  magick ../../emwave_icons/icon_*.svg -resize 40x40 icon_*.png

Place the generated PNGs in this folder. Mapping:
- icon_source_add.png
- icon_block_add.png
- icon_region_pick.png
- icon_snap_magnet.png
- icon_paint_brush.png
- icon_view_reset.png
- icon_layout_preset.png
- icon_composer_cam.png
- icon_help.png
