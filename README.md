# emwave-c

**A real-time 2D FDTD electromagnetic wave simulator, written from scratch in C.**

[![CI](https://github.com/nclalperen/emwave-c/actions/workflows/ci.yml/badge.svg)](https://github.com/nclalperen/emwave-c/actions/workflows/ci.yml)
![Language](https://img.shields.io/badge/language-C11-blue)
![Platforms](https://img.shields.io/badge/platforms-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)

emwave-c solves Maxwell's equations on a 2D Yee grid using the finite-difference
time-domain (FDTD) method and renders the field live at interactive frame rates.
You can paint materials onto the grid, drop sources, place probes and ports, and
watch waves propagate, reflect, and interfere in real time — or run the same
scenes headlessly from the command line to script S-parameter sweeps.

![emwave-c demo](recordings/manual_test.gif)

## Features

**Solver**
- 2D TMz FDTD core on a Yee grid, plain C11 with no external math libraries
- CPML and Mur absorbing boundaries (< −60 dB reflection with the CPML presets)
- OpenMP-parallelized field updates when available; single-threaded fallback
- Materials: dielectrics (εr, σ), PEC, and PMC regions, plus an 11-preset
  material library (copper, FR4, …)
- Sources: CW, Gaussian pulse, Ricker wavelet, and user-defined math
  expressions evaluated by a built-in expression parser
- Measurement ports and frequency sweeps that produce S21 curves as CSV

**Interactive front-ends**
- `emwave` — SDL2 UI with paint tools, dual probes, live scope and FFT panes,
  color-mapped field view, and CSV export
- `emwave_imgui` — Dear ImGui + ImPlot front-end with dockable panels, a
  material-block editor, measurement tools, layout presets, and animation
  recording
- `emwave_cli` — headless runner for batch jobs, CI, and parameter sweeps;
  no SDL dependency at all

**Scene system**
- Scenes are plain JSON (`configs/`), covering grid size, materials, sources,
  ports, boundaries, and sweep settings — documented in
  [`configs/SCHEMA.md`](configs/SCHEMA.md)
- Bundled examples: a dielectric-loaded waveguide and a coplanar-waveguide
  filter set up for S-parameter sweeps ([`docs/SCENES.md`](docs/SCENES.md))
- CLI flags override any JSON value, so one scene file serves many runs

## Quick start

Dear ImGui is vendored as a git submodule, so clone with submodules (or
fetch them into an existing clone):

```bash
git clone --recurse-submodules https://github.com/nclalperen/emwave-c.git
# existing clone: git submodule update --init
```

### Linux / macOS

```bash
# Ubuntu: sudo apt-get install -y cmake ninja-build libsdl2-dev libsdl2-ttf-dev
# macOS:  brew install cmake ninja sdl2 sdl2_ttf
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/emwave
```

### Windows (Visual Studio + vcpkg)

```powershell
.\scripts\setup.ps1     # one-time dependency helper (optional if vcpkg is installed)
.\build-msvc.ps1        # add -Clean for a fresh configure
.\build\Release\emwave.exe
```

For the Dear ImGui front-end, use `.\build-imgui.ps1` and run
`.\build-imgui\Release\emwave_imgui.exe`. MSYS2/MinGW instructions and
troubleshooting live in [`docs/QUICK_START_WINDOWS.md`](docs/QUICK_START_WINDOWS.md)
and [`docs/BUILDING_WINDOWS.md`](docs/BUILDING_WINDOWS.md).

### Headless runs and S-parameter sweeps

```bash
# Fixed-step run of the waveguide scene
./build/emwave_cli --config configs/waveguide.json --run-steps=2000 --profile

# Frequency sweep of the CPW filter → prints S21 per frequency (CSV)
./build/emwave_cli --config configs/cpw_filter.json --run-mode=sweep
```

Useful flags: `--boundary=cpml|mur`, `--probe-log=PATH`, `--profile`.
CLI flags always win over values in the JSON scene.

## Controls (SDL front-end)

| Key | Action |
| --- | --- |
| `Space` / `R` / `C` | Pause · reset · clear fields |
| Mouse | Place and drag sources; paint materials in paint modes |
| `F` | Cycle source type (CW → Gaussian → Ricker) |
| `M` | Toggle Mur / CPML boundaries |
| `S` | Toggle S-parameter ports |
| `B`, `V`, `K` | Theme, accent palette, and field colormap |
| `F1` | Onboarding overlay with the full cheat sheet |

## Architecture

```
src/core/    FDTD updates, CPML boundaries, materials, sources,
             expression parser, ports/analysis, JSON config loader
src/ui/      shared rendering/layout/controls for the SDL front-end
src/app/     entry points: main_new.c (SDL), main_imgui.cpp (ImGui), main_cli.c
include/     public headers shared between modules
configs/     JSON scenes + schema
tests/unit/  Check-based unit tests for the solver and config layers
third_party/ vendored ImGui, ImPlot, jsmn, fonts
```

The solver is a standalone static library (`emwave_core`) with zero UI
dependencies — the SDL and ImGui apps, the CLI, and the unit tests all link
against the same core. Everything builds through a single CMake project;
Windows dependencies are pinned via the vcpkg manifest (`vcpkg.json`).

## Testing

```bash
cd build && ctest --output-on-failure
```

Unit tests cover the FDTD core allocation paths, config loading/validation,
runtime config, analysis instrumentation, and the material library. CI builds
the project and runs the full test suite plus headless CLI smoke runs (fixed
step and sweep mode) on every push and pull request.

## Documentation

- [`docs/README.md`](docs/README.md) — detailed build notes and usage guide
- [`configs/SCHEMA.md`](configs/SCHEMA.md) — JSON scene schema reference
- [`docs/SCENES.md`](docs/SCENES.md) — bundled example scenes
- [`docs/`](docs/) — design docs and work-item history for the UI overhaul

## Roadmap

- Multi-threaded/tiled field updates beyond OpenMP loops
- More bundled scenes (patch antenna, photonic crystal)
- 2.5D extensions and richer dispersive material models

---

Issues and PRs with new scenes, instrumentation ideas, or solver improvements
are welcome.
