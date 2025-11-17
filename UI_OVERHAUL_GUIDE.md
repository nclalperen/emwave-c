# 🎨 UI/UX Overhaul Guide for emwave-c

## 🎯 Current Situation

**Problem**: UI code is mixed with simulation in a 2000+ line monolithic file
**Solution**: Clean modular architecture with complete UI isolation

---

## ✨ What This Architecture Enables

### 1. **Drop-In UI Replacements**

The new architecture lets you **completely replace the UI** without touching simulation code:

```
Simulation Core (fdtd_core.c)
       ↓
   [Interface]
       ↓
   Choose ONE:
   ├── SDL2 (current) ────→ ui_render.c + ui_controls.c
   ├── Dear ImGui     ────→ ui_imgui.cpp
   ├── Qt/QML        ────→ ui_qt.cpp
   ├── Web (WebGL)   ────→ ui_web.js + canvas
   └── Headless      ────→ No UI at all (batch processing)
```

---

## 🎨 UI Replacement Options

### Option 1: **Modern SDL2 UI** (Easiest)
Keep SDL2 but modernize the interface

**Pros:**
- Minimal code changes
- Keep existing renderer
- Fast iteration

**Improvements:**
- Better layout system
- Modern color schemes
- Dockable panels
- Touch-optimized controls
- Material library browser

**Effort:** Low (1-2 weeks)

---

### Option 2: **Dear ImGui** (Recommended for Desktop)
Modern immediate-mode GUI - Perfect for technical/scientific apps

**Pros:**
- Dockable windows out of the box
- Rich widget library
- Built-in plotting
- Very fast development
- Still uses SDL2 for rendering backend

**What You Get:**
```cpp
// Example ImGui code
void render_ui(SimulationState* sim) {
    ImGui::Begin("Simulation Control");

    ImGui::SliderFloat("Frequency (GHz)", &freq_ghz, 0.001, 5.0);

    if (ImGui::Button("Run")) paused = false;
    ImGui::SameLine();
    if (ImGui::Button("Pause")) paused = true;

    ImGui::PlotLines("Ez Field", scope_data, scope_size);

    ImGui::End();
}
```

**Effort:** Medium (2-4 weeks)

**Resources:**
- https://github.com/ocornut/imgui
- https://github.com/ocornut/imgui/wiki/Useful-Extensions

---

### Option 3: **Qt/QML** (Professional Desktop App)
Full-featured GUI framework

**Pros:**
- Native look & feel
- Professional appearance
- Rich component library
- Built-in OpenGL integration
- Excellent documentation

**What You Get:**
- Professional menu systems
- Standard dialogs
- File management
- Printing support
- Help system

**Effort:** High (4-8 weeks)

**Resources:**
- https://www.qt.io/
- Qt Creator IDE

---

### Option 4: **Web Interface** (Cloud-Ready)
Run in browser with WebGL visualization

**Pros:**
- No installation needed
- Cross-platform (mobile too!)
- Modern web UI libraries
- Cloud deployment ready

**Tech Stack:**
```
Backend: C/C++ compiled to WebAssembly
Frontend: React/Vue + Three.js/WebGL
Communication: WebSockets or REST API
```

**What You Get:**
- Access from any device
- Share simulations via URL
- Collaborative features possible
- Mobile-friendly

**Effort:** High (6-12 weeks)

---

## 🛠️ Implementation Roadmap

### Phase 1: Complete Modularization (Current)
- [x] Design module interfaces
- [x] Create header files
- [x] Extract code into modules
- [x] Test and verify

**Status:** Modular build and UI split complete

---

## 🧠 Prompt Library for Complete UI/UX Redesigns

Use the following prompt packs verbatim (or as starting points) when ideating in Midjourney, Stable Diffusion, Figma AI, Galileo, v0, or any other design co-pilot. Each one encodes the instrumentation, physics context, and ergonomic requirements of emwave-c so the generated concepts remain grounded in the simulator's capabilities.

### Prompt 1 – Vision & Design Language
```
Design a premium desktop scientific visualization dashboard for “emwave-c”, a real-time 2D electromagnetic FDTD simulator. Focus on a dark, high-contrast laboratory aesthetic with subtle neon gradients inspired by wave interference. Pair bold typography (Space Grotesk / IBM Plex) with ultra-crisp vector icons, depth via frosted glass cards, and glowing accent lines that suggest field propagation. Include a persistent frequency readout (GHz) and grid size indicator. Convey precision, calm control, and zero visual clutter while keeping the viewport hero-sized.
```

### Prompt 2 – Main Workspace Layout
```
Create the primary simulation screen for emwave-c. Center a 2D field viewport rendered as a heatmap with smooth interpolation, adaptive color bar, and ghosted grid. Surround it with dockable panels: left column for paint tools (PEC, PMC, dielectric), brushes, and layer visibility; right column for live material legend, boundary toggles (Mur vs CPML), and probe summaries. Add resizable bottom tray housing transport controls (Play/Pause, Step, Reset), timestep scrubber, and FPS indicator. The UI must support mouse + tablet input with large touch targets and haptic-style feedback cues.
```

### Prompt 3 – Instrumentation & Scopes
```
Visualize a high-end instrumentation stack for emwave-c featuring dual probes, waveform scopes, and FFT charts. Use synchronized cards that show Probe A/B amplitude vs time, peak/hold annotations, logarithmic frequency axes, and export buttons. Integrate tooltips for Auto-P99 and smoothing options. Allow the cards to undock and float over the viewport while maintaining minimal latency aesthetics (micro animations, neon scanlines). Emphasize data legibility on dark backgrounds.
```

### Prompt 4 – Control Console & Scene Browser
```
Design a right-side control console that surfaces configuration JSONs from /configs, scene presets (waveguide, cpw_filter, cavities), and runtime overrides (nx, ny, tmax). Include expandable accordions for source settings (CW, Gaussian, Ricker), frequency sliders with numeric entry, and boundary/material toggles. Provide a mini preview thumbnail per scene, drag-and-drop ordering, and a “batch/headless” switch. The console should feel like a mission panel with glowing outlines, magnetic snapping, and accessible keyboard shortcuts displayed inline.
```

### Prompt 5 – Onboarding, Status & Help Overlay
```
Create an in-app onboarding overlay for emwave-c that highlights controls (Space pause/resume, M boundary toggle, numeric paint modes). Use contextual callouts anchored to UI elements, a condensed cheat sheet, and animated arrows showing gesture paths. Include system status badges (solver stable, boundary reflections < -60 dB, export complete) in the lower-right corner. Blend the overlay seamlessly with the dark neon theme using translucent panels and soft drop shadows.
```

Mix and match these prompts to cover alternate styles (e.g., lighter palettes, web-first layouts) while keeping the simulator’s scientific instrumentation front and center.

---

### Phase 2A: Modernize SDL2 (Quick Wins)
**Timeline:** 1-2 weeks

**Tasks:**
1. **Layout System** (3 days)
   - Define panel regions
   - Make resizable
   - Add docking

2. **Color Scheme** (2 days)
   - Dark mode / Light mode
   - Color palette picker
   - Syntax highlighting for values

3. **Better Controls** (4 days)
   - Grouped sliders
   - Number input boxes
   - Preset buttons
   - Material library

4. **Improved Visualization** (3 days)
   - False-color schemes (jet, viridis, plasma)
   - Vector field arrows
   - Multiple color maps
   - Transparency controls

**Result:** Much better UX with minimal effort

---

### Phase 2B: Dear ImGui (Recommended Next)
**Timeline:** 2-4 weeks

**Week 1: Setup**
- Integrate Dear ImGui
- Setup SDL2 backend
- Basic window system

**Week 2: Core UI**
- Main menu bar
- Control panels
- Parameter sliders
- Run/pause controls

**Week 3: Visualization**
- Field heatmap in ImGui window
- Plotting with ImPlot
- Oscilloscope display
- FFT visualization

**Week 4: Polish**
- Docking layouts
- Preset management
- Help system
- Keybindings editor

**Result:** Professional-grade UI

---

### Phase 2C: Qt/QML or Web (Long-term)
**Timeline:** 4-12 weeks

Only pursue if:
- Need professional distribution
- Want native integration
- Building product, not research tool
- Have UI/UX designer on team

---

## 📐 Recommended UI Layout (Any Framework)

```
┌─────────────────────────────────────────────────────┐
│  Menu Bar: File | Edit | View | Simulation | Help   │
├──────────────┬──────────────────────┬───────────────┤
│              │                      │               │
│  Toolbox     │   Main Viewport      │  Properties   │
│              │                      │               │
│  ┌─────────┐│   [Field Viz]        │ Frequency     │
│  │Materials││                      │ [==========]  │
│  │ - PEC   ││                      │               │
│  │ - PMC   ││                      │ Sources       │
│  │ - Diel  ││                      │ ☑ Source 1    │
│  └─────────┘│                      │ ☐ Source 2    │
│             │                      │               │
│  ┌─────────┐│                      │ Materials     │
│  │Sources  ││                      │ Paint: [PEC▼] │
│  │ - CW    ││                      │               │
│  │ - Pulse ││                      │ Boundaries    │
│  │ - Ricker││                      │ Type: [CPML▼] │
│  └─────────┘│                      │               │
│             │                      │ [Run]  [Pause]│
├──────────────┴──────────────────────┴───────────────┤
│  Timeline / Oscilloscope / FFT                      │
│  [~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]        │
└─────────────────────────────────────────────────────┘
```

---

## 🎨 Modern Design Principles

### Colors
- **Dark Mode Default** - Easier on eyes for long sessions
- **Accent Colors** - Blue for active, Orange for warnings
- **Semantic Colors** - Red=error, Green=success, Yellow=caution

### Layout
- **Dockable Panels** - User can arrange workspace
- **Collapsible Sections** - Hide what's not needed
- **Consistent Spacing** - 4px/8px/16px grid

### Interactions
- **Hover Feedback** - Show tooltips
- **Keyboard Shortcuts** - Power user friendly
- **Undo/Redo** - Non-destructive editing
- **Real-time Preview** - See changes immediately

---

## 🚀 Quick Start: ImGui Integration

Here's how to add Dear ImGui to your modular codebase:

### 1. Add ImGui to project
```cmake
# CMakeLists.txt
add_subdirectory(external/imgui)
target_link_libraries(emwave PRIVATE imgui)
```

### 2. Create ui_imgui.cpp
```cpp
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include "ui_controls.h"
#include "fdtd_core.h"

void render_imgui_frame(SimulationState* sim, UIState* ui) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    // Your UI code here
    ImGui::Begin("Controls");
    ImGui::SliderFloat("Frequency", &sim->freq, FREQ_MIN, FREQ_MAX);
    if (ImGui::Button("Run")) ui->paused = false;
    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
```

### 3. Update main loop
```cpp
// Just replace render_frame() call
render_imgui_frame(sim, ui);
```

**That's it!** Simulation code unchanged.

---

## 📦 Preset System (Any UI)

One feature you'll want regardless of UI choice:

```cpp
// Preset structure
typedef struct {
    char name[64];
    double freq;
    int source_count;
    Source sources[MAX_SRC];
    // ... other params
} SimulationPreset;

// Save/Load
void preset_save(const char* filename, SimulationState* sim);
void preset_load(const char* filename, SimulationState* sim);

// Built-in presets
const SimulationPreset PRESETS[] = {
    {"Dipole Antenna", 1e9, ...},
    {"Waveguide Resonance", 10e9, ...},
    {"Metamaterial Demo", 5e9, ...},
};
```

---

## 🎓 Learning Resources

### Dear ImGui
- Official: https://github.com/ocornut/imgui
- Gallery: https://github.com/ocornut/imgui/issues/5243
- Video: "Dear ImGui Tutorial" on YouTube

### Qt/QML
- Official: https://doc.qt.io/
- Tutorials: https://www.qt.io/learn

### Web/WebGL
- Three.js: https://threejs.org/
- WebAssembly: https://webassembly.org/

---

## 💡 My Recommendation

**Start with Phase 2A (Modernize SDL2)**
- Quick wins
- Low risk
- Immediate improvement
- Learn what you want in UI

**Then Phase 2B (Dear ImGui)**
- Best effort/reward ratio
- Modern, professional look
- Fast development
- Still works on all your target platforms

**Skip Qt/Web unless:**
- Building a product (not research)
- Have UI/UX designer
- Need mobile support
- Have 2+ months for UI work

---

## ✅ Action Items for You

1. **Finish modularization** (this week)
   - Extract code into modules
   - Get it compiling
   - Verify it works

2. **Try modernized SDL2** (next week)
   - Better layout
   - Dark mode
   - Grouped controls
   - See if it's enough

3. **Prototype Dear ImGui** (following week)
   - Download Dear ImGui
   - Build basic window
   - Compare with SDL2 version

4. **Decide on direction** (end of month)
   - SDL2 modernization enough?
   - or Dear ImGui worth the effort?
   - or stick with current?

---

**The beauty of this architecture: You can try all options without breaking anything!**

Your simulation code is completely isolated. The UI is just a "frontend" now. 🎨
