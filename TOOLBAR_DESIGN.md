# Toolbar Actions Design

This document defines the intended behavior for the left toolbar icons in the ImGui UI. The key principle is:

> Toolbar icons perform actions or enter modes; docked panels are for detailed configuration and management.

Panels remain docked and accessible via the View menu and tab bars. Toolbar tools should never feel like duplicate "open panel" buttons.

---

## 1. Source Tool (icon: source_add)

**Goal:** Place a new source by clicking in the viewport, with a quick configuration dialog that supports both presets and full expressions.

### Flow

1. User clicks the **Source** icon on the toolbar.
2. App enters **place-source mode**:
   - Cursor shows a circular marker over the viewport (snapped to cells).
   - HUD text: `Click to place source, right-click or Esc to cancel`.
3. User left-clicks on the field:
   - We capture grid coordinates `(i, j)` for the active viewport.
   - Simulation is paused while configuring.
   - Open a modal **New Source** dialog.
4. User configures the source and presses **Add** (or **Cancel**).
5. New source is created and appears selected in the **Sources** panel for further editing.

### New Source dialog

Top strip (always visible):

- **Position:** `Cell: (i, j)` (read-only, from click).
- **Field:** radio `Ez / Hx / Hy`.
- **Tabs:** `Basic` | `Expression`.
- Buttons: `Add`, `Cancel` (optionally `Add & Next` to keep placing).

#### Basic tab

First impression for most users; drives a generated expression behind the scenes.

- **Type:** combo `CW / Gaussian / Ricker / Pulse / Custom`.
- **Amplitude:** numeric.
- **Frequency:** default = current simulation frequency.
- **Phase:** numeric (degrees or radians).
- Optional preview line: `Ez(t) ≈ A sin(2π f t + φ)`.

Behavior:

- Changes map into an internal expression template (e.g., `A * sin(2*pi*f*t + phi)`).
- When switching to **Expression** tab, these values are used to prefill the expression and the template selector.

#### Expression tab

Full power; supports arbitrary user expressions.

- **Template** combo, e.g.:
  - `A * sin(omega * t)`
  - `Gaussian pulse`
  - `Ricker wavelet`
  - `CW burst`
  - `Custom (blank)`
- **Expression** multi-line text box.

Behavior:

- Selecting a template pre-fills the expression using current Basic values (A, f, phase).
- User may freely edit the expression (e.g. `10*sin(10*omega*pi)`).
- If the expression diverges from what Basic would generate, mark Basic `Type` as `Custom` and stop auto-overwriting.
- Optionally provide a **Reset to Basic** button to re-sync expression with Basic parameters.

### After Add

- Create a new source using:
  - Position (i, j)
  - Field (Ez/Hx/Hy)
  - Expression string from the Expression tab
  - Any existing source config fields (e.g., enabled flag, port mapping) set to defaults.
- Select the new source in the **Sources** panel and make that panel visible.
- User can still edit all advanced fields in the Sources panel (the toolbar tool is just the fast entry point).

---

## 2. Block / Material Tool (icon: block_add)

**Goal:** Let the user draw a new material block directly on the canvas, then choose material and basic parameters quickly, while keeping full block editing in **Materials / Blocks**.

### Flow

1. User clicks **Block** icon.
2. App enters **draw-block mode**:
   - Cursor indicates block drawing (e.g., crosshair or rectangle outline).
   - HUD text: `Click and drag to draw block, right-click or Esc to cancel`.
3. User drag-selects a rectangle in the viewport.
4. On mouse release, open a **New Block** dialog with the region coordinates prefilled.
5. User configures material and presses **Add**. New block appears in the Materials / Blocks panel selected.

### New Block dialog

Top section:

- **Region:** show normalized or cell coordinates `(x0, y0) → (x1, y1)` (read-only, from drag).

Tabs: `Basic` | `Advanced`.

#### Basic tab

- **Material preset:** combo listing materials from the library (metal/dielectric entries).
- **Type:** `Dielectric / PEC / PMC` (overrides for special cases).
- **Epsilon (epsr):** numeric, prefilled from material preset.
- **Sigma:** numeric, prefilled from material preset.

Behavior:

- Picking a preset copies its properties (epsr, sigma, color, etc.) into the block config.

#### Advanced tab

Optional, for future extension.

- Custom per-cell expressions or parameter fields (e.g., spatially varying epsr).
- Time-dependent behavior tied to a function (analogous to source expression).

### After Add

- Create block in the simulation domain, add it to Materials / Blocks.
- Focus the new block entry in the Materials / Blocks panel for further tweaks.

---

## 3. Region Tool (icon: region_pick)

**Goal:** Provide a general-purpose region selection tool used by multiple features (measurements, composer region capture, advanced grid/domain operations).

### Flow

1. User clicks **Region** icon.
2. App enters **region-select mode**:
   - HUD: `Click and drag to select region, right-click or Esc to cancel`.
   - A rectangle is previewed while dragging.
3. On release:
   - Region `(x0, y0) → (x1, y1)` is stored in an `AppState` field `active_region` in normalized coordinates.
   - Show a small popover **Region Actions**:
     - `Measure area`
     - `Use as composer region`
     - `Set as refinement region` (future grid features)
     - `Cancel`
4. The chosen action delegates to the relevant subsystem:
   - **Measure area:** create an area measurement using this polygon / rectangle.
   - **Composer region:** mark this region as the default region for region-based field exports.
   - Others as we add features.

### Panel interactions

- Grid & Domain / Scenes / Composer windows read the `active_region` to offer “use current region” actions.
- Region tool itself does not open panels; it only updates state and offers a concise action menu.

---

## 4. Snap Tool (icon: snap_magnet)

**Goal:** Toggle and configure snapping for cursor-based tools (source placement, block drawing, region selection, measurements).

### Behavior

- Simple click: toggle global `snap_enabled` on/off.
- Right-click or long-click: open a **Snap Settings** popover:
  - `Snap to grid cells` (default)
  - `Snap step`: numeric (cells or meters)
  - Future: `Snap to probes`, `Snap to vertices`, etc.

Effect:

- When `snap_enabled` is true, all tools that place things (sources, blocks, regions, measurements) quantize coordinates using the configured step.

---

## 5. Paint Tool (icon: paint_brush)

**Goal:** Make material painting an explicit tool with a quick material picker, while leaving the Materials / Blocks panel as the full editor.

### Flow

1. Click **Paint** icon:
   - Toggle `paint_mode` on/off.
   - When turning ON:
     - If no material is selected, open a **Quick Material Picker** popover:
       - Recently used materials (colors + names).
       - A `Browse...` entry that opens the full Material Browser / Materials panel.
2. In paint mode:
   - Cursor shows the paint marker already in use in the app.
   - Clicks paint cells with the chosen material.
3. Click Paint icon again or press its hotkey to exit paint mode.

Panel relation:

- Materials / Blocks / Material Browser remain the place to define and tweak materials.

---

## 6. Reset View Tool (icon: view_reset)

**Goal:** Quickly reset the camera for the active viewport (and optionally other viewports) without touching panels.

### Behavior

- Single left click:
  - Reset the **active viewport** pan and zoom to default (`kDefaultViewportZoom`, pan = 0).
- Right-click or Shift-click:
  - Reset **all viewports** to defaults.

Panel relation:

- No panel opens; this is a pure viewport action.

---

## 7. Layout Tool (icon: layout_preset)

**Goal:** Quickly switch between layout presets and restore a clean dock arrangement.

### Behavior

- Click: cycle through presets in a fixed order, e.g.:
  - Beginner → Power User → Canvas First → Analysis → Beginner...
  - Each click reapplies `apply_layout_preset` (resetting docks).
- Right-click: open a small popover listing presets:
  - `Beginner`
  - `Power User`
  - `Canvas First`
  - `Analysis`
  - `Reset current` (reapply current preset)

Panels:

- This tool never opens/closes individual panels; it only reorganizes them.

---

## 8. Composer Tool (icon: composer_cam)

**Goal:** Treat the Print Composer as a dedicated export workspace; this tool provides fast entry based on the active viewport.

### Behavior

- Click:
  - Open the **Print Composer** window if not already open, and bring it to front.
  - If the composer has no pages yet, create a default page using a template bound to the **active viewport** (e.g., Field + Legend, 1280x720).
  - Optionally set the active page output name to something like `composer_<scene_name>`.
- Right-click:
  - Popover with quick export options:
    - `Open composer`
    - `Export default PNG (field only)`
    - `Export animation (use last composer settings)`

Panels:

- All detailed layout is still done inside the Print Composer window; toolbar just gives a fast entry/shortcut.

---

## 9. Help Tool (icon: help)

**Goal:** Provide immediate orientation (shortcuts, modes, status) without hunting through menus.

### Behavior

- Click:
  - Toggle a **Help Overlay** over the viewport:
    - Lists key shortcuts (run/pause, paint, region, source tool, composer, etc.).
    - Describes current auto-rescale mode and active tools.
  - The overlay is non-interactive and can be dismissed by clicking the Help icon again or pressing Esc.
- Right-click:
  - Bring the **Log** panel to front, for more detailed information.

Panels:

- Log remains the detailed event history; overlay is the quick reference.

---

## Implementation Notes / Plan (high-level)

When implementing:

1. Convert toolbar icon behavior from panel toggles to tool actions and modes as defined above.
2. Ensure each tool updates `AppState` fields (e.g., `place_source_mode`, `draw_block_mode`, `active_region`, `snap_enabled`, `paint_mode`) so they integrate cleanly with event handling and rendering.
3. Keep panel visibility controlled via View menu and tab bars; toolbar should never be the *only* way to reach a panel.
4. For modal dialogs (New Source, New Block), pause simulation while open and resume as appropriate.
5. Add minimal HUD text for active tools so users always know what the next click will do.
6. Wire layout tool to the existing `apply_layout_preset` function, and ensure it’s safe to call repeatedly.
7. Implement the Help overlay as a simple ImGui window rendered on top of the viewport area.

This design should make the toolbar feel like a professional, action-oriented control surface, while the docked panels remain the primary place for configuration and analysis.

Current implementation status (toolbar):

- The Source tool is implemented: clicking the Source icon enters a place-source mode, shows a marker over the active viewport, and opens a **New Source** modal on click. The modal provides Basic and Expression tabs and creates a new runtime source plus config entry.
- The initial version supports **Add** and **Cancel** (no `Add & Next` yet); advanced tweaking still happens in the Sources panel.
