# UI Modernization Brainstorm - Phase 2.75
**Date:** 2025-11-20
**Status:** Planning Phase
**Goal:** Transform emwave-c into a modern, professional simulation tool with Blender-like UX

---

## Current Architecture Analysis

### What We Have ✅
1. **ImGui + SDL2 Hybrid:**
   - SDL2 renders FDTD field to window texture
   - ImGui overlays UI panels on top
   - Docking **ALREADY ENABLED** (`ImGuiConfigFlags_DockingEnable`)

2. **Fixed Layout:**
   - Window: 1920×1080 minimum, **non-resizable**
   - Scale: Fixed at 2 (1 cell = 2 pixels)
   - Panels: Calculated fixed widths (left/right/bottom columns)

3. **Viewport Tracking:**
   - `app.viewport_pos` / `app.viewport_size` calculated
   - Mouse → grid cell conversion works
   - **No zoom or pan controls**

### Current Constraints ❌
1. **Window not resizable** - SDL_WINDOW_SHOWN only (no SDL_WINDOW_RESIZABLE)
2. **Fixed scale = 2** - No zoom in/out
3. **No pan** - Field always centered
4. **Fixed panel layout** - Calculated positions, not user-customizable
5. **SDL renders behind ImGui** - Can't use ImGui::Image for field display

---

## User Requirements

### Explicit Requests:
1. ✅ **Resizable window**
2. ✅ **Complete resizable, modern UI**
3. ✅ **Different screen with magnification (like Blender)**
4. ✅ **Open to more suggestions**

---

## Proposed Improvements - Categories

---

## 1. WINDOW RESIZABILITY (Essential - Quick Win)

### A) Make Window Resizable
**Effort:** ~1 hour
**Impact:** HIGH

**Changes:**
```c
// src/ui/ui_render.c line 422-423
ctx->window = SDL_CreateWindow(title,
                               SDL_WINDOWPOS_CENTERED,
                               SDL_WINDOWPOS_CENTERED,
                               width, height,
                               SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);  // ADD THIS FLAG
```

**Benefits:**
- Users can resize to fit their monitor
- Essential for multi-monitor setups
- Expected behavior for modern apps

**Additional Work:**
- Handle `SDL_WINDOWEVENT_RESIZED` event
- Reflow ImGui layout on resize
- Maintain aspect ratio (optional)

---

## 2. ZOOM & PAN (Critical - Blender-like)

### A) Mouse Wheel Zoom
**Effort:** ~4 hours
**Impact:** HIGH - Game changer for usability

**Implementation:**
```cpp
// In event loop (main_imgui.cpp)
if (event.type == SDL_MOUSEWHEEL) {
    if (app.viewport_valid && mouse_in_viewport) {
        float zoom_delta = event.wheel.y * 0.1f;  // 10% per wheel tick

        // Zoom toward mouse cursor
        float mouse_norm_x = (mx - app.viewport_pos.x) / app.viewport_size.x;
        float mouse_norm_y = (my - app.viewport_pos.y) / app.viewport_size.y;

        scale *= (1.0f + zoom_delta);
        scale = clamp(scale, 0.5f, 32.0f);  // 0.5× to 32× zoom

        // Adjust pan to keep mouse position steady
        // (zoom-to-cursor like Blender)
    }
}
```

**Features:**
- Zoom in: Mouse wheel up (or Ctrl + +)
- Zoom out: Mouse wheel down (or Ctrl + -)
- Zoom range: 0.5× (half size) to 32× (32× magnification)
- **Zoom to cursor** - Blender behavior (zooms toward mouse position)

**Benefits:**
- Inspect fine details (high zoom)
- See entire domain (low zoom)
- Standard CAD/Blender UX

---

### B) Middle Mouse Drag Pan
**Effort:** ~3 hours
**Impact:** MEDIUM-HIGH

**Implementation:**
```cpp
// Pan state in AppState
float viewport_offset_x = 0.0f;
float viewport_offset_y = 0.0f;
bool panning = false;
ImVec2 pan_start_mouse;
ImVec2 pan_start_offset;

// In event loop
if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_MIDDLE) {
    app.panning = true;
    app.pan_start_mouse = ImVec2(mx, my);
    app.pan_start_offset = ImVec2(app.viewport_offset_x, app.viewport_offset_y);
}

if (app.panning && event.type == SDL_MOUSEMOTION) {
    ImVec2 delta = ImVec2(mx - app.pan_start_mouse.x, my - app.pan_start_mouse.y);
    app.viewport_offset_x = app.pan_start_offset.x + delta.x;
    app.viewport_offset_y = app.pan_start_offset.y + delta.y;
}

if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_MIDDLE) {
    app.panning = false;
}
```

**Controls:**
- **Middle mouse drag** - Pan viewport (Blender standard)
- **Shift + Right mouse drag** - Alternative pan (for laptop users)
- **Arrow keys** - Keyboard pan (10 pixels per press)
- **Home key** - Reset zoom/pan to fit

**Benefits:**
- Navigate large domains
- Focus on regions of interest
- Standard Blender/CAD workflow

---

### C) Zoom Controls UI
**Effort:** ~2 hours
**Impact:** MEDIUM

**UI Elements:**
```
Viewport toolbar (bottom-right overlay):
  [−] [Reset] [+]   Zoom: 2.0×   Pan: (0, 0)
```

**Features:**
- Zoom slider: 0.5× to 32×
- Reset button: Fit to viewport
- Zoom indicator: "2.0×" text display
- Pan coordinates: "(50, 30)" pixel offset

---

## 3. MODERN DOCKING UI (High Impact - Uses Existing Feature!)

### A) Full Docking Layout
**Effort:** ~6 hours
**Impact:** VERY HIGH - Professional transformation

**Why This Matters:**
- **Docking is ALREADY ENABLED** (line 2128: `ImGuiConfigFlags_DockingEnable`)
- Currently NOT USED - panels are in fixed columns
- Unlocking this = Blender-like customization

**Implementation:**
```cpp
// Remove fixed layout calculation (lines 2856-2988)
// Replace with DockSpace

ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(),
                              ImGuiDockNodeFlags_PassthruCentralNode);

// Make all panels dockable windows (not child windows)
if (app.show_sources_panel) {
    ImGui::Begin("Sources", &app.show_sources_panel);
    draw_sources_panel_content(...);
    ImGui::End();
}
```

**Benefits:**
- **Drag panels anywhere** - like Blender
- **Split viewport** - horizontal/vertical
- **Tabbed panels** - group related controls
- **Save layouts** - per-project or per-task
- **Multi-monitor** - undock panels to second screen

**Default Layouts:**
1. **Beginner Layout:**
   ```
   ┌──────────────┬─────────────────────┐
   │              │  Simulation         │
   │              │  Controls           │
   │   VIEWPORT   ├─────────────────────┤
   │              │  Sources            │
   │              │  Blocks             │
   └──────────────┴─────────────────────┘
   ```

2. **Power User Layout:**
   ```
   ┌────────┬──────────┬──────────────┐
   │Sources │          │ Material     │
   │Blocks  │          │ Legend       │
   ├────────┤ VIEWPORT ├──────────────┤
   │Scope   │          │ Run Settings │
   │FFT     │          │ Grid Setup   │
   └────────┴──────────┴──────────────┘
   ```

3. **Analysis Layout:**
   ```
   ┌──────────────┬─────────────────────┐
   │   VIEWPORT   │   S-Parameter Plot  │
   ├──────────────┼─────────────────────┤
   │  Scope/FFT   │   Smith Chart       │
   └──────────────┴─────────────────────┘
   ```

---

### B) ImGui Viewports (Multi-Window)
**Effort:** ~2 hours
**Impact:** MEDIUM (for advanced users)

**Enable:**
```cpp
io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
```

**Benefits:**
- Undock panels to separate OS windows
- Multi-monitor workflows
- Second screen for analysis plots

**Caveats:**
- SDL renderer only draws to main window
- Field viewport must stay in main window
- Panels can be separate windows

---

## 4. VIEWPORT RENDERING OVERHAUL (Advanced)

### Option A) Keep SDL Renderer (Current - Simpler)
**Effort:** ~8 hours (zoom/pan only)
**Pros:**
- Hardware accelerated
- Direct pixel manipulation
- Existing code works

**Cons:**
- Can't dock viewport easily
- Rendering happens "behind" ImGui
- Complex zoom/pan calculations

**Implementation:**
- Add zoom/pan transform to SDL rendering
- Adjust `render->scale` dynamically
- Offset rendering by `viewport_offset`

---

### Option B) Switch to ImGui::Image (Modern - Complex)
**Effort:** ~20 hours (full rewrite)
**Pros:**
- Viewport becomes ImGui window
- Can dock/undock field display
- Cleaner architecture

**Cons:**
- Need OpenGL/DirectX texture
- ImGui backend must support textures
- SDL → ImGui texture bridge required

**Implementation:**
1. Create SDL texture for field rendering
2. Convert SDL_Texture → ImTextureID (OpenGL/DX11)
3. Display with `ImGui::Image(texture, size)`
4. Viewport becomes dockable ImGui window

**Recommendation:** Start with Option A, migrate to B later if needed

---

## 5. MODERN UI POLISH (Quality of Life)

### A) Floating Toolbar (Viewport Overlay)
**Effort:** ~3 hours
**Impact:** MEDIUM

**Example:**
```
Viewport (with transparent overlay):
┌─────────────────────────────────────┐
│ ☰ View  🔍 Zoom  📐 Grid  🎨 Theme │ <- Toolbar
│                                     │
│       [Electromagnetic Field]       │
│                                     │
│                 [−] 2.0× [+] ◎      │ <- Zoom controls
└─────────────────────────────────────┘
```

**Features:**
- Transparent background
- Quick access: Grid toggle, colormap, theme
- Zoom controls: -/+/Reset
- View mode: Field/Material/Overlay
- Minimap (optional): Small overview of entire domain

---

### B) Status Bar (Bottom)
**Effort:** ~2 hours
**Impact:** LOW-MEDIUM

```
┌─────────────────────────────────────────────────────────┐
│ Step: 1234 | Time: 4.12 ns | FPS: 60 | Zoom: 2.0×     │
└─────────────────────────────────────────────────────────┘
```

**Info:**
- Simulation step/time
- Render FPS
- Current zoom level
- Mouse grid coordinates (when hovering viewport)
- Selected material/source

---

### C) Dark/Light Theme Toggle
**Effort:** ~1 hour (already exists!)
**Impact:** LOW (already functional)

**Current:** Keyboard 'B' cycles themes
**Add:** Button in toolbar or View menu

---

### D) ImGui Style Customization
**Effort:** ~4 hours
**Impact:** MEDIUM (visual polish)

**Themes:**
1. **Dark Professional** (current - good)
2. **Light Mode** (for presentations)
3. **Blender-inspired** (match Blender colors)
4. **High Contrast** (accessibility)

**Custom styles:**
- Rounded corners (modern look)
- Larger fonts (4K displays)
- Color accents (highlight active panels)

---

## 6. VIEWPORT ENHANCEMENTS

### A) Grid Overlay Improvements
**Effort:** ~3 hours
**Impact:** MEDIUM

**Features:**
- **Adaptive grid** - Denser at high zoom, sparser at low zoom
- **Grid labels** - Show coordinates (0.00m, 0.01m, etc.)
- **Rulers** - Top/left edges with measurements
- **Origin marker** - Cross at (0,0)

---

### B) Measurement Tools
**Effort:** ~6 hours
**Impact:** MEDIUM (professional feature)

**Tools:**
1. **Ruler Tool (R key):**
   - Click two points
   - Shows distance in meters/cells
   - Shows angle
   - Displays on viewport

2. **Area Tool:**
   - Click to define polygon
   - Shows area in m²
   - Useful for material coverage

3. **Coordinate Display:**
   - Mouse cursor shows (x, y) in meters
   - Shows field value at cursor
   - Bottom-left viewport corner

---

### C) Minimap (Optional)
**Effort:** ~8 hours
**Impact:** LOW-MEDIUM

**Example:**
```
Bottom-right corner:
┌──────┐
│▓▓░░░░│ <- Minimap (64×64)
│▓▓░░░░│
│░░░░░░│
│ [■]  │ <- View rectangle
└──────┘
```

**Features:**
- Shows entire domain
- Red rectangle = current view
- Click to jump to location
- Useful for large simulations (>1000×1000)

---

## 7. KEYBOARD SHORTCUTS & USABILITY

### A) Zoom/Pan Shortcuts
**Effort:** ~2 hours
**Impact:** HIGH

**Shortcuts:**
| Key | Action |
|-----|--------|
| **Mouse Wheel** | Zoom in/out (to cursor) |
| **Ctrl + +** | Zoom in (centered) |
| **Ctrl + -** | Zoom out (centered) |
| **Ctrl + 0** | Reset zoom to fit |
| **Middle Mouse Drag** | Pan viewport |
| **Shift + Right Drag** | Pan viewport (alt) |
| **Arrow Keys** | Pan 10 pixels |
| **Shift + Arrow** | Pan 50 pixels |
| **Home** | Reset zoom + pan |
| **F** | Frame selected (zoom to selection) |

---

### B) Quick Actions Menu
**Effort:** ~3 hours
**Impact:** MEDIUM

**Example:** Right-click in viewport
```
┌─────────────────────────┐
│ Add Source Here         │
│ Add Material Block Here │
│ Measure Distance...     │
│─────────────────────────│
│ Zoom to Fit             │
│ Reset View              │
│─────────────────────────│
│ Copy Coordinates        │
│ Properties...           │
└─────────────────────────┘
```

---

## 8. ADVANCED FEATURES (Future)

### A) Multi-Viewport (Blender-style)
**Effort:** ~16 hours
**Impact:** HIGH (for advanced users)

**Layout:**
```
┌──────────────┬──────────────┐
│  Viewport 1  │  Viewport 2  │
│  (Field Ez)  │  (Field Hx)  │
├──────────────┼──────────────┤
│  Viewport 3  │  Viewport 4  │
│  (Material)  │  (FFT)       │
└──────────────┴──────────────┘
```

**Uses:**
- Compare Ez vs. Hx simultaneously
- Field vs. Material overlay comparison
- Scope + FFT + Field at once

---

### B) 3D Visualization (Future)
**Effort:** ~40 hours
**Impact:** VERY HIGH (publication quality)

**Features:**
- OpenGL 3D field rendering
- Rotate/zoom in 3D space
- Volumetric rendering
- Isosurface visualization
- Export to video

**Note:** Major undertaking, consider external tool integration (Paraview)

---

### C) Animation Recording
**Effort:** ~8 hours
**Impact:** MEDIUM

**Features:**
- Record simulation to video (MP4)
- Frame-by-frame export (PNG sequence)
- Animated GIF export
- Configurable framerate/quality

---

## IMPLEMENTATION PRIORITY

### Phase 2.75A: Essential Modernization (~16h)
**Goal:** Resizable, zoom/pan, basic Blender-like UX

1. ✅ **Window Resizability** (1h)
2. ✅ **Mouse Wheel Zoom** (4h)
3. ✅ **Middle Mouse Pan** (3h)
4. ✅ **Viewport Toolbar** (3h)
5. ✅ **Zoom Controls UI** (2h)
6. ✅ **Status Bar** (2h)
7. ✅ **Keyboard Shortcuts** (2h)

**Outcome:** Professional viewport interaction like Blender

---

### Phase 2.75B: Docking UI (~12h)
**Goal:** Fully customizable layout

1. ✅ **Full DockSpace** (6h)
2. ✅ **Panel → Window conversion** (3h)
3. ✅ **Default layouts** (2h)
4. ✅ **Layout save/load** (2h)

**Outcome:** User-customizable interface

---

### Phase 2.75C: Polish (~16h)
**Goal:** Production-ready UX

1. ✅ **Grid overlay improvements** (3h)
2. ✅ **Measurement tools** (6h)
3. ✅ **Quick actions menu** (3h)
4. ✅ **ImGui style customization** (4h)

**Outcome:** Professional, publication-quality tool

---

### Phase 2.75D: Advanced (Optional - ~24h+)
**Goal:** Power-user features

1. 🔲 **Multi-viewport** (16h)
2. 🔲 **Minimap** (8h)
3. 🔲 **Animation recording** (8h)
4. 🔲 **ImGui::Image migration** (20h)

**Outcome:** Research-grade simulation platform

---

## TECHNICAL DECISIONS

### Decision 1: SDL Renderer vs. ImGui::Image
**Recommendation:** Keep SDL renderer for Phase 2.75, migrate later

**Rationale:**
- SDL renderer works well
- Zoom/pan achievable with transform math
- ImGui::Image requires backend rewrite
- Can migrate incrementally later

---

### Decision 2: Docking Layout Strategy
**Recommendation:** Full DockSpace with default layouts

**Rationale:**
- Docking already enabled but unused
- Unlocks massive UX improvement
- Standard in modern apps (Blender, Visual Studio)
- Users expect customization

---

### Decision 3: Zoom Range
**Recommendation:** 0.5× (half) to 32× (32× magnification)

**Rationale:**
- 0.5× - See large domains (1000×1000 cells)
- 1× to 4× - Normal working range
- 8× to 32× - Inspect fine details (single-cell resolution)
- Matches Blender/CAD zoom ranges

---

### Decision 4: Pan Control
**Recommendation:** Middle mouse drag (primary), Shift+Right (fallback)

**Rationale:**
- Middle mouse = Blender standard
- Shift+Right = Laptop users (no middle button)
- Arrow keys = Precision adjustment
- Multiple input methods = accessibility

---

## COMPATIBILITY CONSIDERATIONS

### Cross-Platform:
- ✅ Windows: All features supported
- ✅ Linux: All features supported
- ✅ macOS: Requires testing (SDL2 behavior)

### Performance:
- Zoom/pan: Negligible overhead (<1ms)
- Docking: Built into ImGui (optimized)
- SDL scaling: Hardware accelerated

### Backward Compatibility:
- Save layouts to `imgui.ini` (per-user)
- Config JSON unchanged
- Old SDL2 frontend still works

---

## USER EXPERIENCE IMPROVEMENTS

### Before (Current):
```
❌ Fixed 1920×1080 window
❌ Fixed 2× scale (no zoom)
❌ Can't inspect fine details
❌ Can't see large domains
❌ Fixed panel layout
❌ No customization
```

### After Phase 2.75A (Essential):
```
✅ Resizable window (fit any monitor)
✅ Zoom 0.5× to 32× (mouse wheel)
✅ Pan viewport (middle mouse)
✅ Zoom to cursor (Blender-like)
✅ Keyboard shortcuts
✅ Viewport toolbar
✅ Status bar
```

### After Phase 2.75B (Docking):
```
✅ Drag panels anywhere
✅ Split viewport
✅ Tabbed panels
✅ Save custom layouts
✅ Multi-monitor support
✅ Blender-like flexibility
```

### After Phase 2.75C (Polish):
```
✅ Grid with labels/rulers
✅ Measurement tools (ruler, area)
✅ Right-click quick actions
✅ Custom themes
✅ Professional appearance
```

---

## RECOMMENDATION SUMMARY

### Immediate (Must-Have):
1. **Window resizability** (1h) - Trivial, huge impact
2. **Mouse wheel zoom** (4h) - Expected feature
3. **Middle mouse pan** (3h) - Essential for large domains

### High Priority (Should-Have):
4. **Full DockSpace** (6h) - Unlocks docking (already enabled!)
5. **Viewport toolbar** (3h) - Modern UX
6. **Keyboard shortcuts** (2h) - Power user efficiency

### Medium Priority (Nice-to-Have):
7. **Grid improvements** (3h) - Professional appearance
8. **Measurement tools** (6h) - Research feature
9. **Default layouts** (2h) - User convenience

### Future (Optional):
10. **Multi-viewport** (16h) - Advanced users
11. **Animation export** (8h) - Publications
12. **ImGui::Image migration** (20h) - Architectural improvement

---

## ESTIMATED TOTAL EFFORT

| Phase | Features | Hours |
|-------|----------|-------|
| **2.75A: Essential** | Resizable, zoom, pan, shortcuts | 16h |
| **2.75B: Docking** | Full customizable UI | 12h |
| **2.75C: Polish** | Grid, tools, themes | 16h |
| **2.75D: Advanced** | Multi-viewport, recording | 24h+ |
| **TOTAL (A+B+C)** | **Production-ready modern UI** | **44h** |

---

## NEXT STEPS

### For User:
1. **Review this document** - Which features resonate?
2. **Prioritize** - Essential/Nice-to-have/Future
3. **Choose approach:**
   - **a)** Implement Phase 2.75A only (16h - quick win)
   - **b)** Implement Phase 2.75A + 2.75B (28h - full modernization)
   - **c)** All phases 2.75A+B+C (44h - production-ready)
   - **d)** Custom selection (pick individual features)

### For Implementation:
- Create detailed prompts for each phase
- Break into 2-4 hour tasks
- Test on Windows/Linux/macOS
- Document keyboard shortcuts
- Update README with new features

---

**Status:** ✅ BRAINSTORMING COMPLETE - AWAITING USER DECISION

*This is your planning assistant's comprehensive analysis. Choose your path forward!*
