import os
import zipfile

# Configuration for the icons (Monoline style)
ICON_SIZE = 64
STROKE_WIDTH = 4
STROKE_COLOR = "white"
FILL_COLOR = "none"


# Helper to wrap paths in standard SVG header
def create_svg(content: str) -> str:
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {ICON_SIZE} {ICON_SIZE}" fill="{FILL_COLOR}" stroke="{STROKE_COLOR}" stroke-width="{STROKE_WIDTH}" stroke-linecap="round" stroke-linejoin="round">
{content}
</svg>'''


# Icon Definitions (SVG Paths)
icons = {
    "icon_source_add.svg": """
        <circle cx="16" cy="48" r="5" fill="white" stroke="none"/>
        <path d="M26 38 L54 10 M38 48 A 24 24 0 0 1 56 26 M46 52 A 36 36 0 0 1 60 30" />
        <path d="M10 10 L20 10 M15 5 L15 15" stroke-width="3"/> """,

    "icon_block_add.svg": """
        <path d="M32 6 L54 18 L54 46 L32 58 L10 46 L10 18 Z" />
        <path d="M32 6 L32 30 M54 18 L32 30 L10 18" />
        <path d="M8 8 L18 8 M13 3 L13 13" stroke-width="3"/> """,

    "icon_region_pick.svg": """
        <rect x="16" y="16" width="40" height="40" rx="2" stroke-dasharray="6 6" opacity="0.8"/>
        <path d="M4 4 L20 20 M4 20 L20 4" stroke-width="3" /> <circle cx="12" cy="12" r="9" />
    """,

    "icon_snap_magnet.svg": """
        <path d="M10 20 A 10 10 0 0 1 30 20 L30 36 A 6 6 0 0 1 18 36 L18 54 A 6 6 0 0 1 6 54 L6 36 A 6 6 0 0 1 10 20" transform="rotate(-45 32 32) translate(5,5)"/>
        <path d="M50 14 L58 6 M44 8 L46 4 M56 16 L60 14" /> """,

    "icon_paint_brush.svg": """
        <path d="M48 4 L60 16 L40 36 L28 24 Z" />
        <path d="M28 24 L10 54 L4 58 L8 52 L26 34" />
        <line x1="38" y1="18" x2="50" y2="30" /> """,

    "icon_view_reset.svg": """
        <path d="M6 28 L32 6 L58 28" />
        <path d="M12 34 L12 56 L24 56 L24 40 L40 40 L40 56 L52 56 L52 34" />
    """,

    "icon_layout_preset.svg": """
        <rect x="4" y="8" width="56" height="48" rx="2" />
        <line x1="16" y1="8" x2="16" y2="56" />
        <line x1="44" y1="8" x2="44" y2="56" />
        <line x1="16" y1="44" x2="44" y2="44" />
    """,

    "icon_composer_cam.svg": """
        <circle cx="32" cy="32" r="26" />
        <circle cx="32" cy="32" r="10" />
        <path d="M32 6 L32 12 M32 52 L32 58 M6 32 L12 32 M52 32 L58 32" />
    """,

    "icon_help.svg": """
        <circle cx="32" cy="32" r="28" />
        <path d="M24 24 A 10 10 0 0 1 42 20 C 44 26 32 34 32 40" />
        <circle cx="32" cy="50" r="0.5" fill="white" stroke-width="5"/>
    """
}

# Execution
output_dir = "emwave_icons"
zip_name = "emwave_icons.zip"

os.makedirs(output_dir, exist_ok=True)

print(f"Generating {len(icons)} icons in '{output_dir}'...")

with zipfile.ZipFile(zip_name, 'w') as zipf:
    for filename, svg_content in icons.items():
        full_path = os.path.join(output_dir, filename)
        with open(full_path, "w", encoding="utf-8") as f:
            f.write(create_svg(svg_content))
        zipf.write(full_path, arcname=filename)
        print(f"  - Created {filename}")

print(f"\nDone! Icons zipped into: {zip_name}")
print("Tip: rasterize with ImageMagick, e.g.:")
print("  magick emwave_icons/icon_*.svg -resize 40x40 assets/icons/")
