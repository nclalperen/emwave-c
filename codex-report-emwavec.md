# codex-report-emwavec

## Purpose and scope
- Comprehensive code-derived report for the `emwave-c` project based on repository state under `C:\projects\emwave-c`.
- Focused on authoritative sources in the codebase rather than potentially stale documentation files.
- Captures architecture, build system, configuration, runtime behavior, UI layers, tests, assets, scripts, and observed risks.
- Intended as a deep technical reference; structured for readability and traceability to source files.
- All content is ASCII and reflects the current working tree; no edits were made to existing project files.

## Recent updates (composer/headless)
- Headless `--ch-layout` now uses its own path buffer and loads before applying templates/output overrides.
- Layout JSON loader reads per-page items (type, viewport, pos/size, region), background/transparent/video settings, and seeds `composer_next_item_id` accordingly.

## Repository identity
- Project name: `emwave-c`.
- Language stack: C11 for core and SDL UI; C++17 for ImGui/ImPlot front-end.
- Primary domain: real-time 2D FDTD electromagnetic simulation with visualization and instrumentation.
- Core build system: CMake with optional vcpkg toolchain integration.
- Optional parallelism: OpenMP when available.

## Directory layout highlights
- Root build files: `CMakeLists.txt`, `vcpkg.json`.
- Core sources: `src/core` contains simulation engine, configuration loader, materials, sources, boundaries, ports, analysis, expression engine.
- Application entry points: `src/app` hosts CLI main, SDL UI main, ImGui main, shared runner.
- SDL UI modules: `src/ui` includes controls, rendering, layout, render test.
- Legacy snapshots: `src/legacy` contains historical monolithic main variants.
- Headers: `include` mirrors core/app interfaces.
- Config samples: `configs` with waveguide and CPW filter JSON plus schema.
- Tests: `tests/unit` with Check-based suites and CMake glue.
- Scripts: `scripts` for Windows and MSYS2 build/setup automation.
- Third-party: `third_party` for jsmn, check stub, imgui, implot, fonts (referenced by build scripts).
- Assets and outputs: `build`, `build-imgui`, `dist`, `recordings`, `frame_*.bmp`, `measurements.csv`, `ports_sparams.txt`.

## Build targets and roles
- `emwave_core` (static library): full engine with ports enabled; used by UI and CLI targets.
- `emwave_solver` (static library): engine variant with `EMWAVE_ENABLE_PORTS=0` for lightweight solver use.
- `emwave` (executable): SDL2 renderer UI front-end (`src/app/main_new.c`).
- `emwave_imgui` (executable): Dear ImGui + SDL_Renderer2 UI front-end (`src/app/main_imgui.cpp`), C++17.
- `emwave_cli` (executable): headless solver driver (`src/app/main_cli.c`).
- `render_layout_test` (executable): layout regression checker for UI sizing (`src/ui/render_layout_test.c`).

## Build system behavior
- CMake minimum 3.20; project languages C and CXX.
- `CMAKE_MODULE_PATH` extended with `cmake/`.
- Testing enabled globally (`enable_testing()`).
- OpenMP detected with `find_package(OpenMP)` and linked when `OpenMP_C_FOUND`.
- UI dependencies conditionally resolved: prefers config packages (vcpkg) then fallback `find_package(SDL2)`/`SDL2_ttf`.
- SDL2 discovery adapts to `VCPKG_INSTALLED_DIR` and `VCPKG_TARGET_TRIPLET`.
- UI targets only built when both SDL2 and SDL2_ttf are found; otherwise warning and CLI-only build.
- Include directories shared: `include`, `src`, `third_party/jsmn`.
- MSVC flags: `/W4 /permissive-`, `_CRT_SECURE_NO_WARNINGS`, `_CRT_NONSTDC_NO_DEPRECATE`, OpenMP via `/openmp:llvm` with `/wd4849`.
- Non-MSVC flags: `-Wall -Wextra -Wno-unused-parameter`.
- Installation: `emwave` (if built) and `emwave_cli` installed to root; `assets/` copied if present.

## Build options and toggles
- `EMWAVE_ENABLE_UI` (ON by default): controls SDL UI build; when off, only CLI builds.
- `EM_WAVE_ENABLE_BOUNDS_CHECK` (OFF by default): enables bounds-checked inline field accessors.
- `EMWAVE_ENABLE_PORTS` default 1 via `config.h`; overridden to 0 for `emwave_solver`.
- OpenMP optional; messages indicate presence or absence at configure time.

## Dependency inventory
- SDL2 and SDL2_ttf for rendering/text in SDL front-end.
- ImGui + ImPlot vendored for Dear ImGui front-end.
- jsmn JSON parser (`third_party/jsmn/jsmn.c`) for config loading.
- Check testing framework (with stub fallback) for unit tests.
- Font: DejaVu Sans bundled in assets (`assets/fonts/DejaVuSans.ttf`).
- Optional ffmpeg (external) used by ImGui recording pipeline; not a build-time dependency but runtime expectation for encoding.

## Configuration model (SimulationConfig, include/config.h)
- Grid: `nx`, `ny` (cells); defaults 400x400; min 64, max 4096; total cells capped at ~8.3M (`SIM_MAX_CELLS`).
- Domain: `lx`, `ly` meters; defaults 0.6 each; `dx = lx/nx`, `dy = ly/ny`.
- CFL safety: `cfl_safety` default 0.95; must be (0,1).
- Steps per frame: `steps_per_frame` default 2 (used for UI and progress cadence).
- Sweep: `sweep_points` (default 5, max 32), `sweep_start_hz` (default 5e8), `sweep_stop_hz` (default 3e9), `sweep_steps_per_point` default 200.
- Run mode: `SIM_RUN_MODE_FIXED_STEPS` or `SIM_RUN_MODE_SWEEP`; default fixed.
- Run steps: default 5*200; ignored in sweep mode.
- Boundary mode: `SIM_BOUNDARY_CPML` (default) or `SIM_BOUNDARY_MUR`.
- Profiling: `enable_profile` flag for CLI to time runs.
- Probe logging: `enable_probe_log` flag and `probe_log_path` (default `probe.txt`).
- Materials: `material_rect_count` and array of `MaterialRectSpec` (x0,y0,x1,y1, epsr, sigma, tag).
- Sources: `source_count` and array of `SourceConfigSpec` (active, x, y, type, amp, freq, sigma2, field, expr).
- Ports: `port_count` (up to 2) with normalized `PortConfigSpec` (active, x, y0, y1).
- Frequency UI range constants: `FREQ_MIN=1e6`, `FREQ_MAX=5e9`.
- Scene limits: `CONFIG_MAX_MATERIAL_RECTS=16`, `MAX_SRC=4`.

## Configuration defaults (src/core/config_runtime.c)
- Grid 400x400, domain 0.6x0.6 m, CFL 0.95.
- Steps per frame 2.
- Sweep: 5 points, 0.5–3 GHz, 200 steps/point.
- Run mode fixed, run_steps = 1000.
- Boundary CPML.
- Profiling disabled; probe logging disabled; probe_log_path "probe.txt".
- One default material rect: dielectric epsr=4.0, sigma=0 at normalized box (0.4,0.45)-(0.6,0.55), tag=dielectric.
- One default source: active CW at (0.25,0.5), amp=1, freq=1 GHz, sigma2=4.0, field Ez, empty expr.
- Ports: none active by default.

## Configuration loader (src/core/config_loader.c)
- CLI parsing supports long options with `--config` file overlay then flag overrides.
- Flags parsed: `--help`, `--config=|--config`, `--nx`, `--ny`, `--lx`, `--ly`, `--cfl`, `--sweep-points`, `--sweep-start`, `--sweep-stop`, `--sweep-steps`, `--run-mode`, `--run-steps`, `--boundary`, `--profile`, `--probe-log`, `--no-probe-log`.
- JSON parsing via jsmn: expects top-level object; keys `simulation`, `materials`, `sources`, `ports`.
- Simulation overlay fields: nx, ny, lx, ly, cfl, steps_per_frame, sweep points/start/stop/steps, run_mode string, run_steps, enable_profile, enable_probe_log, probe_log_path, boundary string.
- Materials array parsed into rect specs; type strings map to tags: pec->1, pmc->2, default dielectric.
- Sources array parsed into configs: type strings cw/gaussian/ricker/expr; field strings ez/hx/hy; expr copied raw; sigma2, freq, amp validated.
- Ports array parsed into normalized configs; clamps count to MAX_PORTS.
- File size guard: `CONFIG_LOADER_MAX_FILE_BYTES` default 10 MB; rejects oversized, missing, directory inputs, or parse errors.
- Error messages routed via `errbuf` when provided; fallback stderr logging otherwise.
- Saving: `save_config_json` writes simulation block, materials, sources, optional ports with minimal escaping (expr not escaped).

## Configuration validation and clamping (src/core/config_runtime.c)
- `config_clamp_to_limits`: enforces nx/ny within min/max and total cell cap (scales down proportionally), clamps CFL to (0.1,0.99), sweep points to [1,SWEEP_MAX_POINTS], sweep steps min 50, run_steps non-negative, material/source/port counts within limits, normalizes rect/source/port coordinates to [0,1], fixes negative sigma2/freq, enforces valid source field enum, boundary mode defaults to CPML if invalid.
- `config_validate`: rejects too small/large grids, cell overrun, invalid CFL range, sweep ranges with stop<=start, invalid material rect bounds or epsr<=0, invalid source freq/sigma2, invalid port y-range ordering, invalid boundary enum.
- `config_print_summary`: logs grid, sweep, run mode, boundary, probe logging status, counts of materials/sources/ports.

## Data structures (include/types.h)
- `Source`: flags, grid indices ix/iy, type, field, amp/freq, pulse params t0/tau, sigma2 footprint, expr text buffer (256), compiled expr pointer.
- `Port`: x location, y range y0/y1 inclusive, len, buffer length n, V/I double buffers, head index, active flag.
- `Scope`: ring buffer length n, head, data pointer y, on flag, last sample, rolling_absmax, rolling_generation.
- `BoundaryType` enum: Mur or CPML.
- `CpmlPreset`: name, smax, kmax, amax, thickness.
- `CpmlState`: enabled flag, thickness, preset index, sigma/kappa/alpha max, boundary type, coefficient arrays (kx/bx/cx, ky/by/cy) with capacities.
- `SimulationState`: config copy; grid sizes; domain lengths; cfl_safety; dx/dy/dt; sweep settings; fields Ez/Hx/Hy/Ez_old (double** + data blocks); CPML psi fields; material properties epsr/sigma_map/tag_grid; sources array; ports array; ports_on flag; timestep; freq; cpml state; cached step_Ez_absmax.

## FDTD core lifecycle (src/core/fdtd_core.c)
- `fdtd_init`: allocates SimulationState, clamps/validates config, computes dx/dy/dt via CFL, initializes boundary, allocates all field grids, initializes ports (if enabled), clears fields, sets default materials and applies configured materials, initializes ports from config, initializes sources and CPML coefficients, zeros psi, sets step_Ez_absmax to 0.
- Allocation helpers use failure injection hook `fdtd_test_set_alloc_fail_after` to simulate malloc failures; zeroes state on free.
- Grid allocation uses contiguous data blocks plus row-pointer arrays for cache-friendly access.
- Failure path frees all resources and returns NULL.
- `fdtd_free`: releases fields, psi, materials arrays, tag grid, ports buffers, boundary coefficients, compiled expressions, then frees SimulationState.
- `fdtd_clear_fields`: zeros all field arrays and psi arrays, resets step_Ez_absmax.
- `fdtd_reset`: clears fields, resets timestep to 0, zeros CPML psi, resets port buffers and re-evaluates len/active; preserves geometry/config.
- `fdtd_update_grid_for_freq`: recomputes dx/dy based on lambda/(sqrt(EPSR_MAX_SCENE)*TARGET_CPW), capped by MAX_SCALE_FACTOR*BASE_DX; recomputes dt, updates freq, rebuilds CPML coefficients, zeros psi.

## Time stepping (src/core/fdtd_core.c)
- `fdtd_compute_dt`: CFL timestep based on dx, dy, c0, cfl_safety; uses sqrt of reciprocal squared spacings.
- `fdtd_step` sequence:
  - Copies Ez to Ez_old for Mur boundary support (parallelizable).
  - Updates Hx, Hy from spatial derivatives of Ez; applies CPML corrections near boundaries; zeros Hx/Hy inside PMC-tagged cells.
  - Updates Ez from curl(H) with conductivity term (sigma) via ceze/cezh coefficients; applies CPML corrections; clamps Ez to 0 in PEC-tagged cells.
  - Applies Mur boundaries if CPML disabled.
  - Injects all active sources with Gaussian spatial footprint and saturation to reduce reflections.
  - Samples max |Ez| across updated field and source neighborhoods into step_Ez_absmax.
  - Samples ports if ports_on flag set (delegates to ports_sample).
  - Increments timestep.
- OpenMP pragmas guard major loops for parallel updates.
- DIVISION_SAFETY_EPSILON used to avoid division by zero in conductivity term.

## Boundary system (src/core/boundary.c)
- `CPML_PRESETS`: Gentle (smax=1.0,kmax=3.0,amax=0.03,thick=10), Default (1.2,5.0,0.05,12), Aggressive (1.8,6.0,0.08,16).
- `boundary_init`: initializes CpmlState defaults, preset index 1 (Default), boundary_type Mur, enabled=0, thickness=PML_THICK.
- `boundary_set_type`: switches between CPML and Mur, toggles enabled flag.
- `cpml_build_coeffs`: resizes coefficient arrays to nx/ny, computes sigma/kappa/alpha grading per dimension with cubic profile, stores bx/by and cx/cy exponentials; ensures rx/ry computed symmetrically on both sides; handles thickness guard.
- `cpml_apply_preset`: selects preset, clamps thickness to less than half of min grid dimension minus one; rebuilds coefficients and zeros psi.
- `cpml_zero_psi`: zeros all psi arrays in SimulationState.
- `apply_mur_boundaries`: first-order Mur on all four edges using local c from epsr; uses Ez_old buffers.
- `boundary_shutdown`: frees CPML coefficient arrays and resets capacities.

## Materials system (src/core/materials.c)
- `materials_init`: resets grid to defaults then applies each configured rectangle.
- Rect application: converts normalized coords to cell indices with clamp; for PEC tag sets epsr=1, sigma=0, tag=1; PMC tag sets tag=2; dielectric sets epsr and sigma_map with tag=0.
- `materials_reset_to_defaults`: sets all cells to vacuum epsr=1, sigma=SIGMA_BG, tag=dielectric.
- `paint_material_at`: interactive toggle for PEC (tag flip + reset to vacuum), PMC (tag flip), or dielectric (toggle between vacuum and specified epsr) with provided eps value; clamps to grid bounds.

## Source system (src/core/sources.c)
- `sources_init`: seeds four sources; only first active by default; positions spread horizontally; applies padding to avoid edges; sets default CW 1 GHz amp=1 sigma2=4; computes pulse params; compiles expressions if expr type specified in config; deactivates unused sources beyond config count.
- `source_reparam`: for Gaussian/Ricker recalculates t0 and tau based on freq (cycles 6 and 2).
- `sources_set_freq`: updates freq and reparam for all sources.
- `sources_cycle_type`: cycles CW->Gaussian->Ricker (expr excluded) and recompiles expressions.
- `source_time_value`: returns value at timestep for CW, Gaussian, Ricker, or expression (expr_eval) with amp/freq injection.
- `inject_source_into_Ez`: soft source injection into selected field (Ez/Hx/Hy) around ix/iy in 5x5 stencil with Gaussian weight and saturation to reduce reflections.
- `inject_all_sources`: calls injector for all; no early exit for inactive (checked inside).
- `find_nearest_source`: utility for UI to pick source by mouse position and scale within max distance.

## Expression engine (src/core/expr.c)
- Implements shunting-yard tokenizer and RPN evaluator for tiny expression set.
- Tokens: numbers, variables (t, amp, freq, pi), funcs (sin, cos, exp, sqrt, abs), operators (+,-,*,/,^), parentheses, comma.
- Handles unary minus by injecting zero on expect_operand when minus encountered.
- Emits program as array of ops (push const/var, unary func, add/sub/mul/div/pow).
- `expr_eval`: stack-based evaluation with local fixed-size stack; protects divide by zero with epsilon 1e-18; sqrt uses fabs before sqrt.
- Error handling: returns 0 on invalid parse/eval; error strings propagated via errbuf; OOM paths handled.
- `expr_free`: frees ops and program struct.

## Ports and S-parameters (src/core/ports.c)
- Ports are optional via EMWAVE_ENABLE_PORTS; buffers length 4096 samples (`PORT_SIGNAL_LENGTH`).
- `ports_init`: allocates buffers for MAX_PORTS (2), chooses x positions at ~nx/4 and ~3nx/4, y-span centered (1/4..3/4) with minimum segment length 2, clamps to interior; fails if grid too small or allocation fails, cleaning up.
- `ports_apply_config`: applies normalized port specs from config, clamping x to interior (1..nx-2) and y0/y1 to grid; sets active flag; len recalculated on sampling.
- `ports_sample`: per-step integration; prepares clamped y-span with minimum length; adjusts px inside interior; deactivates ports lacking buffers or with invalid geometry; accumulates V=sum(Ez)*dy and I=sum(Hy)*dx over span; writes to circular buffers at head and increments head.
- `ports_free`: frees V/I buffers and resets geometry.
- `port_prepare_sample_range`: guards against invalid buffer/grid and enforces min span; updates len when changed.

## Analysis and logging (src/core/analysis.c)
- Allocation failure hook `analysis_test_set_alloc_fail_after` for tests.
- Oscilloscope: `scope_init` allocates ring buffer (min length 64), `scope_push` writes sample, tracks rolling_absmax with decay 0.995, `scope_clear` zeros buffer, `scope_free` releases memory.
- FFT export: `dump_scope_fft_csv` computes Hann-window DFT (O(N^2)) up to requested Nfft (clamped to buffer length, min 64), removes DC, writes freq/mag CSV with header.
- S-parameter: `compute_s21` computes magnitude ratio of port1 vs port0 using cosine/sine projection at freq*dt; returns 0 if insufficient data.
- Probe writer: buffered writer with 4K block; `probe_writer_open/append/flush/close/is_open`.
- Legacy `probe_open`/`probe_log` kept for monolithic builds.

## CLI runner (src/app/cli_runner.c)
- `simulation_runner_options_from_config`: sets progress interval (clamped to steps_per_frame or default 50), chooses probe log path (config or default probe.txt).
- `simulation_runner_init/reset_progress`: initializes runner struct, sets total steps.
- `simulation_runner_on_step`: optionally opens probe log lazily, appends probe sample, samples ports if enabled, prints progress at interval.
- `run_fixed_mode`: chooses probe at grid center, runs fixed number of steps (config run_steps or fallback points*steps_per_point or 1000), logs progress, flushes, prints completion and probe path if enabled.
- `run_sweep_mode`: activates ports, prints CSV header to stdout, opens sweep_s21.csv with header, iterates sweep points linearly spaced start→stop, resets sim per point, retunes grid for freq, runs steps_per_point, computes s21 and writes to stdout and file, reports file write.
- `cli_runner_execute`: bootstraps runner, profiles elapsed if enabled, runs fixed or sweep, prints steps/s when profiled, returns success flag.

## Bootstrap layer (src/core/app_bootstrap.c)
- `simulation_bootstrap_from_args`: parses CLI config, validates, clamps, prints summary, allocates SimulationState via fdtd_init, configures boundaries (CPML preset vs Mur), returns status (-1 fatal, 0 help/validation, 1 success).
- `simulation_bootstrap_with_config`: same but using provided config directly.
- `simulation_bootstrap_shutdown`: frees SimulationState via fdtd_free and clears pointer.
- Boundary configuration ensures CPML preset applied and psi zeroed when CPML mode selected.

## CLI entry point (src/app/main_cli.c)
- Bootstraps from args, exits 0 on help/validation failure, 1 on fatal bootstrap.
- Prints OpenMP availability and thread count if compiled with _OPENMP.
- Runs `cli_runner_execute`, shuts down bootstrap, returns 0 on success else 1.

## SDL UI front-end (src/app/main_new.c)
- Bootstraps simulation from CLI args (so same flags apply).
- Initializes UIState (`ui_state_init`), RenderContext (`render_init`), scope buffer sized to nx*scale.
- Sets UI layout based on render dimensions (scale, menu/timeline/panels) and simulation grid.
- Initializes SimulationRunner for optional logging; disables progress printing for UI loop.
- Main loop: updates metrics, handles events (pause, paint, layout adjustments, theme, colormap, ports/log toggles, scene reload), steps solver `steps_per_frame` when unpaused, pushes probe to scope, routes probe logging via runner, tracks FPS via SDL performance counters.
- Rendering: uses `render_frame` to draw field, UI panels, scope, overlays; saves screenshot to frame.bmp when requested.
- Scene preset reloads: supports IDs 1 (waveguide) and 2 (cpw_filter); searches relative/SDL base path; rebootstrap simulation and resync UI/scope/runner; prints status.
- Shutdown: runner, scope, render, UI state, bootstrap.

## SDL UI state and controls (src/ui/ui_controls.c/h)
- UI flags: running, paused, legend/help toggles, paint mode/type (PEC/PMC/dielectric), paint_eps, autoscale toggles for color/scope, hold flags for color/scope, color autoscale mode (peak/P99 placeholder), render stride, probe positions (primary/secondary), drag source index, sweep state arrays, s21 amplitude cache, steps_per_frame, theme (dark/light), accent palette index, colormap selection, log_probe, screenshot request, scene reload flags, scene name buffer, sliders for frequency and speed, layout dimensions, scope dock position.
- Slider handling: freq slider logarithmic mapping (freq_from_slider/slider_from_freq), speed slider linear 1–50.
- Event handling: mouse for paint or dragging nearest source else move probe; right-click toggles paint; motion continues paint/drag; keyboard shortcuts for quit (Esc/Q), help (F1), scene reload (F5/F6), FFT export (F3), pause (Space), legend (L), paint mode (M/U), paint type (I), paint eps adjust (O/P), theme (B), accent (V/Shift+V), colormap (K/Shift+K), screenshot (F2), CPML presets (7/8/9), layout resize ([ ] , . - =), scope dock (\\), toggle sources (1/2/3), cycle source types (T), autoscale toggles (A/H/J), clear/reset fields (C/R), boundary toggle (Y), probe log toggle (G), ports toggle (S), debug force metrics (F).
- Metrics: `ui_update_metrics` decays estimated vmax/scope, refreshes from step_Ez_absmax and scope rolling_absmax, supports force recompute flag for full scan fallback; updates held ranges when hold disabled.
- Layout adjustments: responds to panel/timeline resize keys via clamp helpers.

## SDL rendering layer (src/ui/ui_render.c/h)
- Theme palettes: dark and light with detailed color sets for backgrounds, borders, text, sliders, scope; accent palettes (6 presets) for primary/secondary colors.
- Colormaps: classic (rainbow-like), viridis-inspired, plasma-inspired; selectable globally and for direct render calls.
- Render layout computed in `render_layout_compute` (src/ui/ui_layout.c): derives canvas size from nx*scale, positions panels, timeline, colorbar, scopes, block outline with clamps and margins.
- Field rendering: `draw_grid` plots Ez with colormap mapping centered at 0 (0.5 baseline), PEC cells in light gray, PMC in blueish; returns max absolute value. Channel rendering variant supports Hx/Hy/|H|/Sx/Sy/|S|/Ez/|Ez|; flags not available for Ex/Ey/Hz set not-available output.
- Sources rendering: draws crosshairs at active source positions with accent color.
- Block outline rendering: rectangular overlay around default block region from layout.
- Colorbar: vertical gradient matching chosen colormap in properties panel.
- Scope rendering: plots ring buffer across width with vertical scaling from held or live vmax; returns observed max.
- Info panel: displays frequency (GHz), steps/frame, dx/dt, FPS, grid size and boundary/ports status, probe location, scene name.
- Status badges: bottom-right panel showing solver state (paused/running), boundary preset, reflection note, ports status.
- Help overlay: semi-transparent layer with grouped shortcut text; legend overlay in toolbox panel shows condensed shortcuts.
- Menu stub text drawn at top; viewport uses SDL logical size matching layout; window size auto-adjusts to layout dimensions each frame.
- Screenshot: `save_screenshot` reads renderer pixels and saves BMP.

## UI layout helper (src/ui/ui_layout.c)
- Computes viewport/panel rectangles with clamps for panel widths and timeline height.
- Colorbar positioning respects margins inside properties panel; clamps height/width to available space.
- Scope rectangles derived for properties and timeline docks with padding and min heights.
- Timeline alignment ensures width matches viewport and y offset equals menu_bar + viewport_h.
- Block outline derived as centered rectangle scaled to grid size; clamped to canvas bounds.

## Render layout test (src/ui/render_layout_test.c)
- Validates colorbar within properties panel, block outline within canvas, window dimensions consistent with panel sums, timeline positioning, scope width alignment, timeline x/width expectations.
- Runs checks for two grid sizes (128x96, 512x384); prints failures with details; exits nonzero on failure.

## ImGui/ImPlot front-end (src/app/main_imgui.cpp) overview
- Uses SDL2 + SDL_ttf + ImGui backends (SDL2 + SDL_Renderer2) and ImPlot for plotting.
- Supports viewport layouts: single, horizontal split, vertical split, quad.
- Visualization modes: Ez, |Ez|, Hx, Hy, |H|, Sx, Sy, |S|, Ex/Ey/Hz placeholders, material map, overlay.
- Wizard state (`WizardState`) stores editable config mirror (grid, boundary, sweep, materials, sources, flags advanced/basic).
- Scene operations: add/remove sources, add/remove material blocks (including from material library presets), drag sources/blocks, paint interactions, reload scene from config, save scene via `save_config_json`, resync runtime sim from wizard, rebootstrap simulation.
- Recording subsystem: supports GIF/MP4/PNG sequences; captures frames into memory vector of SDL_Surface; writes outputs under `recordings/` with timestamped names; logs last ffmpeg command to `recordings/ffmpeg_last.txt`; states idle/active/processing/error with status messages and progress; resolution scaling and frame targets configurable; auto-play flag.
- Menus/panels: scene panel with counts, grid panel edits nx/ny/lx/ly/CFL/boundary, sources panel list with controls, blocks panel with material presets and filters, run settings panel (mode fixed/sweep, run_steps, sweep params, profile/probe logging paths), recording panel UI, overlay of source IDs and blocks, measurement/annotation handling indicated by artifacts (measurements.csv).
- Event handling: mouse hit tests for sources/blocks, drag-and-drop, keyboard shortcuts mapped to actions similar to SDL UI plus recorder and wizard sync.
- Uses material library for display and block presets with colors and IDs; filter by selected material supported.
- Maintains AppState with recorder, UI log, scale, layout preferences, basic/advanced mode, log entries, selection indices.
- Supports screenshot/animation export, playback toggles, and UI log entries for status messages.

## Material library (src/core/material_library.c/h)
- Data-driven array of materials with categories and types:
  - Metals (PEC approximation): Copper (conductivity 5.96e7), Gold (4.1e7), Silver (6.3e7), Aluminum (3.5e7).
  - Dielectrics: Air (epsr 1), FR4 (epsr 4.4, tan_delta 0.02), Rogers RO4003 (epsr 3.55, tan_delta 0.0027), PTFE (epsr 2.1), Silicon (epsr 11.68, small conductivity), Glass (epsr 6.0), Alumina (epsr 9.8).
- Each entry includes description, category (metal/dielectric/magnetic/custom placeholder), type (PEC/PMC/dielectric/lossy), frequency_dependent flag, color RGB, ID.
- API: init/shutdown, count retrieval, get by index/name/id, get by category (static array results), getters for epsilon/mu/sigma.
- Used primarily by ImGui front-end for block presets and display.

## Configuration examples (configs/)
- `waveguide.json`: nx=512, ny=256, lx=0.5, ly=0.25, cfl=0.9, sweep 3 points 6–8.5 GHz, steps/point=400; materials: PEC top/bottom bands, dielectric slab epsr 2.4; sources: Gaussian at (0.1,0.5) freq 7 GHz amp 1.
- `cpw_filter.json`: nx=400, ny=400, lx=0.6, ly=0.6, cfl=0.92, sweep 5 points 3–6 GHz, steps/point=500; materials: dielectric half-space epsr 3.2, PEC center conductor and ground tabs; sources: two CW at 4.5 GHz with second inactive (amp 0, active false).
- `invalid_config.json`: minimal invalid content for tests.
- `SCHEMA.md`: documents expected JSON structure; not relied on for this report.

## Assets and example outputs
- Font: DejaVu Sans shipped under assets (copied by scripts/build).
- Screenshots: frame_000.bmp .. frame_005.bmp present at root; UI screenshot target `frame.bmp`.
- Recording artifacts: `recordings/` directory intended for captured animations; ImGui writes ffmpeg command to `recordings/ffmpeg_last.txt`.
- Measurements: `measurements.csv` stores area/annotation sample, likely from ImGui measurement tool.
- Ports sample: `ports_sparams.txt` contains example S21 amplitude at 1 GHz.
- Probe/FFT outputs: generated at runtime as `probe.txt`, `scope_fft.csv`, `sweep_s21.csv` depending on mode.

## Scripts and automation (scripts/)
- `auto_validate.ps1`: orchestrates validation (details not inspected here but available).
- `build.ps1`: wrapper for building (likely legacy).
- `build_msvc.bat`: batch build using MSVC; aligns with PowerShell wrappers.
- `build_msys2.sh`: MSYS2/MinGW build helper configuring cmake and copying font.
- `install_deps_vcpkg.bat`: installs dependencies via vcpkg.
- `setup.ps1`: environment setup (vcpkg and dependencies).
- `build-imgui.ps1`: specialized MSVC build for imgui target only (keeps build directories separate).
- `build-msvc.ps1`: main MSVC build helper with -Clean option.

## Tests (tests/unit)
- Framework: Check, with stub fallback (`third_party/check_stub`) when Check not found.
- CMake (`tests/unit/CMakeLists.txt`): defines helper `add_emwave_check_test`, links to emwave_core and Check stub, sets `EMWAVE_TEST_TMP_DIR`, defines CONFIGS_DIR and TEST_TMP_DIR for config loader tests, adds targets for config_loader_tests, fdtd_core_tests, analysis_tests, config_runtime_tests, material_library_tests.
- `test_config_loader.c`: validates parsing of waveguide/cpw_filter configs, invalid config rejection, missing file/dir handling, oversize guard, empty file rejection, truncation of excess sources, ports parsing clamp, CLI overrides precedence, profile flag effect.
- `test_config_runtime.c`: clamps dimensions and cell limits, rejects invalid CFL, rejects sweep stop<=start, clamps counts, normalizes rect/source coords and signs, rejects invalid material rect bounds and source sigma2.
- `test_fdtd_core.c`: dt computation matches expected formula, fallback to default CFL on invalid input, handles allocation failure injection (state and partial), handles port alloc failure (via analysis hook), clamps small grid to SIM_MIN_DIM, recovers from previous failure, resets ports with negative lengths, normalizes port geometry, reactivates valid ports, asserts field clearing effects.
- `test_analysis_alloc.c`: scope init failure safety, scope min length enforcement, ports init failure cleanup, small grid rejection, ports enablement success, ports_sample skipping unbacked ports, clamping px inside grid, clamping segment and updating len and buffer values, minimum span deactivation, deactivation when grid lacks interior, compute_s21 amplitude ratio test.
- `test_material_library.c`: simple asserts for init, lookup by name/id, category filtering; contains non-ASCII and malformed printf tokens that may warn/fail on strict compilers.

## Output artifacts and logging
- Probe logging: path from config/CLI; lazy-opened in runner; flushed at end; default probe.txt.
- Sweep output: stdout header `# freq_Hz,s21_mag` and CSV lines; file `sweep_s21.csv` with header duplicated (# and plain) plus data.
- Scope FFT: `scope_fft.csv` generated via SDL UI hotkey F3; Hann-windowed DFT magnitudes vs frequency.
- Ports sample: `ports_sample` writes to circular buffers; S21 computed on demand, not persisted unless sweep mode.
- Recorder outputs (ImGui): GIF/MP4/PNG sequences saved under `recordings/`; temp dirs used for MP4 encoding; log file `ffmpeg_last.txt` recorded for debugging.
- Screenshots: SDL UI saves `frame.bmp` on request; existing BMP frames likely pre-rendered.

## Performance and scaling considerations
- OpenMP used in field updates and Mur boundary loops when compiled with _OPENMP; thread count reported at runtime in CLI/UI.
- CPML adds overhead due to auxiliary fields and coefficient computations; thickness clamped to maintain stability on small grids.
- `fdtd_update_grid_for_freq` scales dx/dy with frequency sweeps to maintain cells-per-wavelength target; may affect runtime step counts if applied mid-run.
- `ports_sample` integrates over y-span each step; overhead scales with port span length and number of ports.
- `dump_scope_fft_csv` is O(N^2) and may be slow for large Nfft; buffer length defaults to scope width (nx*scale in SDL UI).
- Autoscaling in UI uses decayed maxima to avoid full scans each frame; force recompute flag triggers full pass when needed.

## Error handling and safeguards
- Config loader guards against oversize files, directory inputs, parse errors, and invalid CLI values; reports via errbuf/stderr.
- Simulation init validates config and aborts on invalid settings; allocation failures handled with cleanup and null returns; test hooks ensure coverage.
- CPML thickness clamped to less than half of min dimension; warns if preset thickness exceeds safe limit.
- Ports sampling deactivates on invalid buffers/geometry to avoid out-of-bounds.
- Expression engine reports unknown identifiers and mismatched parentheses via errbuf.
- Probe writer handles failures to open with stderr warning and suppresses further attempts.
- Recorder error state set when frame capture fails or encoder returns failure; logs message to UI log.

## Known limitations and risks (code-observed)
- `save_config_json` writes expressions without JSON escaping; quotes/newlines in expr could produce invalid JSON.
- S21 computation is unnormalized magnitude ratio; no calibration for impedance or windowing, may be qualitative only.
- FFT export uses naive DFT; large buffers will be slow and may block UI.
- `test_material_library.c` contains non-ASCII characters and malformed printf strings; may fail on strict compilers or without proper console encoding.
- Ports disabled in `emwave_solver` library; users expecting ports must link `emwave_core`.
- Recorder depends on external ffmpeg availability; not bundled; failure sets error state and logs to file/UI.
- Mur boundary is first-order; higher reflections relative to CPML; warning printed only when UI toggles? Not explicit in code but inherent.
- UI autoscale P99 mode placeholder (not implemented; uses cached p99 flags without computation).
- No GPU acceleration; performance limited to CPU and memory bandwidth.
- No JSON schema validation beyond custom checks; unrecognized keys ignored silently.

## Boundary and CPML details
- CPML coefficients computed per axis with cubic grading and alpha taper; kappa, sigma, alpha arrays sized to nx/ny.
- CPML preset index stored in CpmlState; boundary type stored for runtime toggles.
- Mur boundaries use local epsr to adjust wave speed per edge cell, reducing reflection when permittivity varies near boundary.
- `boundary_is_cpml_enabled` helper used by UI to decide psi zeroing and coefficient rebuilds.

## Ports and instrumentation nuances
- Ports default span is roughly half the grid height; min length enforced as 2 cells.
- Voltage sampling uses Ez integral over segment times dy; current uses Hy integral times dx.
- Circular buffers store most recent samples at head; S21 uses head offsets to align series.
- Ports on/off toggled via UI key S; CLI sweep forces activation regardless of config.
- Port geometry clamped every sample call; span adjusted to meet min length if provided y0/y1 too short.

## Source handling nuances
- Normalized positions from config clamp to [0,1] then mapped to cells with padding of 2 to avoid edges.
- Expression sources compiled once at init or when type cycled; failures log warning and leave expr_program null.
- Source injection applies saturation to reduce reflections: field += val/(1+|val|*0.1).
- Source cycling skips expression sources; toggling source activity via number keys in UI.
- Source footprint sigma2 defaults to 4.0; enforced positive in clamping.

## Material handling nuances
- Rect coordinates normalized then converted using floor/ceil to inclusive cell ranges.
- PEC application resets epsr to 1 and sigma to background to avoid lingering conductivity.
- Paint toggles allow rapid interactive prototyping; ensures tags updated consistently with epsr/sigma.
- Tag grid distinguishes dielectric (0), PEC (1), PMC (2); used by boundary enforcement in field updates.

## UI interaction summary (SDL)
- Mouse: left drag source or probe, paint when paint mode on; right click toggles paint mode.
- Probe positions clamped to grid bounds; probe2 exists but optional flag probe2_active.
- Sweep state tracked but not executed in SDL loop (sweep handled in CLI); s21_amp and sweep arrays maintained for potential UI overlay.
- Layout adjustments allow dynamic resizing of panels and timeline to fit content; logical size adjusted each frame to computed layout.
- Theme and colormap changes apply immediately; global direct colormap for external callers set via ui_render_set_colormap.
- Legend/help overlays provide shortcut reference; status badges show boundary preset and ports state.

## UI rendering specifics (SDL)
- Uses SDL_Renderer accelerated; logical size set to computed layout; window minimum size 1280x720 enforced.
- Text rendering via SDL_ttf with fallback font search: bundled path relative to SDL base, DejaVuSans.ttf in cwd, or system path `/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf`.
- Offset_x/y in RenderContext default 0; could be used for panning (not exposed in SDL UI).
- Colorbar width default 18 px with margins; clamped inside properties panel.
- Timeline scope reserves space for controls; scope dock can switch between timeline and properties.
- Autoscale values stored in UI state; held values respected when hold flags set.

## ImGui recording pipeline specifics
- Recorder struct holds state, format (GIF/MP4/PNG), framerate, target_frames, frame_count, resolution_scale, SDL_Surface frame vector, capture size, output path, auto_play flag, progress, status message.
- start_recording: allocates frame buffer, resets counters, sets output path with timestamp/extension, logs to UI, creates recordings directory (mkdir variants for Windows/Unix).
- stop_recording: transitions to processing, writes frames depending on format, uses ffmpeg for MP4/GIF (commands logged), clears frames, sets state back to idle or error.
- capture_frame_if_recording: grabs frame from SDL renderer each frame, stores surface, updates progress, auto-stops at target_frames.
- Status updates shown in recording panel; error transitions log to UI log and set state to RECORDING_ERROR.

## Example data files in repo root
- `frame_000.bmp` .. `frame_005.bmp`: likely example captures.
- `measurements.csv`: entries for area and annotation with coordinates.
- `ports_sparams.txt`: sample S21 amplitude output at 1 GHz.
- `imgui.ini`: ImGui settings persistence for ImGui front-end.

## Build variants and usage notes
- `emwave_core` linked by all main binaries except `emwave_solver` which is ports-disabled; choose correct lib based on instrumentation needs.
- CLI always built; UI targets conditional on SDL presence.
- ImGui target requires C++17; linked with SDL2/SDL2_ttf plus implot/imgui sources.
- Install step copies assets directory if it exists; ensure fonts/assets placed before install to include them.
- On MSVC, OpenMP uses LLVM backend; warnings about collapse disabled.

## Logging and profiling
- CLI profile flag uses `clock()` to time run and prints steps, seconds, and steps/s.
- Progress printing every progress_interval steps in CLI runner; disabled in SDL UI runner.
- Warnings printed for missing SDL dependencies, failed probe log open, failed port init, CPML thickness clamp, failed config load, failed scene reload, failed screenshot.
- ImGui UI uses internal log window to show status messages including recording.

## Performance-sensitive loops
- FDTD updates over nx*ny cells; CPML adds per-edge cost; OpenMP used where available.
- Ports sampling sums over y-span each step; consider span length when enabling ports.
- Scope buffer length equals nx*scale in SDL UI, so large grids increase scope memory and FFT cost.
- Sweep mode reinitializes grid for each frequency, ensuring cell size adapts to wavelength target; cost includes repeated resets and coefficient rebuilds.

## Memory layout choices
- Field arrays allocated as contiguous blocks with row-pointer arrays to improve cache locality.
- Psi arrays for CPML follow same pattern.
- Tag grid stored as unsigned char grid with contiguous data for boundary checks.
- Ports buffers allocated via analysis_checked_calloc to reuse failure hook; length fixed 4096 doubles per port per V/I.

## Coordinate and scaling conventions
- Normalized config coordinates [0,1] map across domain lengths (lx, ly) and grids (nx, ny).
- BASE_DX/DY derived from default geometry (Lx/NX, Ly/NY) used as reference for freq-based scaling in AUTO mode.
- Source positions clamped with padding to avoid boundary interference; port x clamped to interior (1..nx-2).
- UI rendering scales grid cells by `RenderContext.scale` (default 2) for viewport pixel size.

## Interaction between config and runtime state
- SimulationState stores config copy used by UI/runner.
- Boundary mode can be switched at runtime; CPML presets applied immediately with coefficient rebuild and psi zeroing.
- Source freq changes via UI slider call fdtd_update_grid_for_freq and sources_set_freq, impacting dx/dy/dt and CPML coefficients.
- Reset/clear operations zero fields and psi, ensuring clean slate for new configurations.
- Port activation can be toggled at runtime; geometry persists but buffers reset on reset or when deactivated by sampling guard.

## File and path handling
- Config loader accepts absolute and relative paths; CLI `--config` can be specified as separate arg or with equals.
- SDL UI preset loader searches relative to executable base path if not absolute; climbs up to 5 parent directories attempting to find config file.
- Recorder paths use `recordings/` with timestamps; handles Windows and Unix mkdir variants.
- Probe log default `probe.txt` in current working directory unless overridden.

## Numerical parameters and constants
- Physical constants: c0, MU0, EPS0 defined in config.h.
- CFL safety defaults and caps as above.
- TARGET_CPW=12 cells per wavelength in slowest medium; EPSR_MAX_SCENE=4; MAX_SCALE_FACTOR=20 for auto grid scaling.
- SIGMA defaults: background 0, block 0; PML_THICK default 12 unless preset overrides.
- Histogram/color constants in config.h unused in core logic but relevant for visualization.

## CLI usage patterns (code-based)
- Fixed run: `emwave_cli --config configs/waveguide.json --run-steps=2000` runs set steps, probe at center, progress printed, optional probe log.
- Sweep run: `emwave_cli --config configs/cpw_filter.json --run-mode=sweep --sweep-points=5 --sweep-steps=500` prints S21 CSV and saves file.
- Boundary override: `--boundary=mur` switches to Mur.
- Profile flag adds steps/s summary.
- Probe logging enabled via `--probe-log=path` or disabled via `--no-probe-log`.
- CLI defaults apply even when launched from SDL UI (main_new uses same bootstrap).

## UI hotkeys (SDL) concise list
- Quit: Esc, Q.
- Pause: Space.
- Help overlay: F1.
- Scene presets: F5 waveguide, F6 cpw_filter.
- FFT export: F3.
- Screenshot: F2.
- Theme: B toggles dark/light.
- Accent: V cycles (Shift+V reverse).
- Colormap: K cycles (Shift+K reverse).
- Paint mode: M or U; paint type I; eps adjust O/P.
- Sources: 1/2/3 toggle, T cycle types.
- Boundaries: Y toggle CPML/Mur; 7/8/9 CPML presets.
- Ports: S toggle; Probe log: G toggle.
- Autoscale/holds: A toggle auto, H hold color, J hold scope.
- Layout: [ ] left panel, , . right panel, - = timeline height, \ scope dock.
- Reset/Clear: R reset (state), C clear fields.

## Test-focused behaviors
- Allocation failure hooks ensure graceful failures and clean state resets.
- Port resets zero buffers and disable when geometry invalid; tests verify head reset and len recalculation.
- Config loader guards for oversized files; tests write synthetic oversize to trigger rejection.
- Empty file rejection validated; errbuf expected non-empty on failures.
- Sweep start>stop rejection validated in config_validate test.
- FFT magnitude ratio test confirms compute_s21 respects amplitude scaling.

## Additional notes on legacy and artifacts
- `src/legacy` contains numerous historical main variants (`mainv0.c` etc.) not used by modular build but preserved.
- `orig_main_imgui.cpp` stored at root likely legacy ImGui prototype, not built.
- Patch file `my_changes.patch` exists in root; not applied by build.
- `validation_log.txt`, `PHASE_2.75D*` docs present but not used in code.
- `imgui.ini` stores ImGui window positions/settings for ImGui front-end persistence.

## Observations about modular architecture
- Engine decoupled from UI via `emwave_core` and `SimulationState`.
- Shared runner used by CLI and UI for probe logging and sweep stepping.
- UI layers consume only public headers (`config.h`, `fdtd_core.h`, `analysis.h`, `ui_*`).
- CPML and Mur boundary handling encapsulated in boundary module; toggled via API from UI.
- Config system separates parsing/validation/clamping from runtime application.

## Potential integration points
- External tooling could generate configs and invoke `emwave_cli` for batch sweeps.
- ImGui front-end exposes material library, making it suitable for rapid prototyping of microwave structures.
- SDL UI minimal but responsive, suited for live visualization and painting.
- Ports/S21 instrumentation enables simple two-port measurements; more ports would require code changes (`MAX_PORTS`).

## Safety and correctness considerations
- Tag grid ensures boundary conditions and materials enforced cell-by-cell.
- PEC clamp and PMC zeroing in field updates prevent unbounded fields in conductor regions.
- Conductivity term uses symmetrical update with division guard to avoid instability.
- CPML coefficient formulas include alpha to mitigate low-frequency instabilities.
- Mur boundary uses local epsr to reduce reflection when permittivity varies near boundary.
- Saturating source injection reduces numerical artifacts from large amplitudes.
- CPML thickness clamp prevents exceeding grid interior.

## Memory and resource cleanup
- fdtd_release_state_resources frees psi, fields, epsr/sigma/tag, ports, boundary coefficients, compiled expressions, and zeros state.
- sources_shutdown frees compiled expr programs for each source.
- ports_free zeroes and frees all V/I buffers, resets geometry.
- scope_free frees buffer and resets struct.
- probe_writer_close flushes buffer and frees memory.
- boundary_shutdown frees coefficient arrays and resets capacities.

## Data export and interoperability
- Configs saved as JSON with human-readable formatting; materials/sources serialized with floats and booleans.
- Probe logs plain text "timestep value" per line.
- Sweep CSV includes header and data with scientific notation.
- Scope FFT CSV includes header comment with dt and N, columns freq_Hz, mag.
- Recorder outputs standard media formats using ffmpeg; GIF/MP4/PNG widely interoperable.

## Platform considerations
- Windows: MSVC build scripts, vcpkg manifest, OpenMP via /openmp:llvm, SDL2_ttf config from vcpkg path; recorder uses _mkdir; path checks for drive-letter absolute paths.
- Linux/macOS/MSYS2: Make/Ninja builds supported; pkg-config fallback for SDL2/SDL2_ttf; mkdir with mode 0755; font fallback path `/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf`.
- SDL UI assumes renderer available; window resizable; logical size scaling used.
- ImGui code uses std::filesystem; requires C++17 and platform support for wide paths on Windows (uses std::filesystem::absolute).

## Sweep mechanics
- run_sweep_mode resets simulation per frequency to ensure clean initial conditions.
- Ports activated automatically in sweep; CPML reinitialized by fdtd_reset/ fdtd_update_grid_for_freq.
- Probe logging still available during sweep if enabled; sweeps write S21 to stdout and file.
- Frequency interpolation uses linear lerp across points; steps per point configurable.

## Scope and probe behavior
- Probe position in CLI fixed at grid center; UI probe moves with mouse clicks; UI scope uses probe samples each step.
- Scope rolling_absmax decays slowly, allowing UI to display near-real-time amplitude envelope.
- FFT export uses circular buffer tail Nfft samples.
- Probe logging toggled by UI, uses SimulationRunner to lazily open file when first sample requested.

## Color and visualization
- Default field colormap in SDL UI is selected by UI state (classic by default); direct rendering uses viridis unless changed via ui_render_set_colormap.
- Material overlays override colors for PEC/PMC to visually distinguish conductors and magnetic boundaries.
- Block outline provides visual reference for default material rectangle region.
- Status badges and info panel use theme colors and accent highlights for readability.

## Scene presets and reload flow (SDL)
- Preset descriptors: id/name/path pairs for "waveguide" and "cpw_filter".
- load_scene_preset_config attempts direct path, then relative to SDL base path up to five parents.
- On reload: runner shutdown, scope freed, bootstrap shutdown, rebootstrap with new config, scope reinitialized, UI layout resynced, runner reinitialized, scene name applied, confirmation printed.
- Reload flag set from UI events (F5/F6).

## Code health observations
- Modular structure with clear separation between core, config, UI, and tests.
- Extensive unit test coverage for config and allocation paths; physics correctness tests limited to S21 ratio check and dt formula.
- Legacy code retained but isolated; modern build uses modular files.
- Warnings managed per compiler; OpenMP optional with fallbacks.
- Probe/log/recorder paths configurable; defaults reasonable.

## Suggested validation steps (actionable)
- Configure and build: `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release` (or MSVC via scripts).
- Run tests: `cd build && ctest --output-on-failure` (or `ctest -C Debug` for multi-config).
- Run CLI smoke: `./build/emwave_cli --config ../configs/waveguide.json --run-steps=1000 --profile`.
- Run sweep smoke: `./build/emwave_cli --config ../configs/cpw_filter.json --run-mode=sweep --sweep-points=3 --sweep-steps=200`.
- Launch SDL UI: `./build/emwave` and exercise paint, probe log, CPML presets, screenshot.
- Launch ImGui UI: build via imgui script/target and verify recording (ffmpeg present) and material library interactions.

## Line count assurance
- This report is intentionally verbose to exceed 600 lines as requested; each bullet or header occupies its own line.

## Closing notes
- All content sourced from current code and repository structure; no reliance on possibly outdated docs.
- File generated without modifying existing project files; new file `codex-report-emwavec.md` added for reference.

## Function map quick reference (selected implementations)
- `fdtd_init` sets grid, allocates fields, applies materials, initializes ports/sources/CPML.
- `fdtd_free` releases all simulation allocations and compiled expressions.
- `fdtd_clear_fields` zeros field arrays and psi buffers.
- `fdtd_reset` clears fields, resets timestep, resets CPML psi, and normalizes ports.
- `fdtd_update_grid_for_freq` retunes dx/dy/dt and CPML for new frequency.
- `fdtd_compute_dt` computes CFL-safe timestep from dx/dy and cfl_safety.
- `fdtd_step` performs one Yee update with CPML or Mur boundaries and source injection.
- `fdtd_epsilon_at` returns EPS0*epsr at cell (i,j).
- `fdtd_sigma_at` returns conductivity map value at cell (i,j).
- `boundary_init` seeds CpmlState defaults and preset.
- `boundary_set_type` toggles CPML vs Mur and enables CPML flag.
- `boundary_is_cpml_enabled` reports CPML usage for UI toggles.
- `cpml_apply_preset` loads preset params, clamps thickness, rebuilds coefficients.
- `cpml_build_coeffs` computes sigma/kappa/alpha arrays and bx/cx/by/cy.
- `cpml_zero_psi` zeros psi arrays after rebuilds or resets.
- `apply_mur_boundaries` applies first-order absorbing boundaries using Ez_old.
- `materials_init` resets materials to vacuum and applies configured rectangles.
- `materials_reset_to_defaults` sets epsr=1, sigma=SIGMA_BG, tag=dielectric grid-wide.
- `paint_material_at` toggles PEC/PMC/dielectric states interactively.
- `sources_init` seeds source structs from defaults and config with padding clamps.
- `sources_shutdown` frees compiled expression programs.
- `sources_set_freq` updates freq and pulse params for all sources.
- `sources_cycle_type` cycles non-expression source types.
- `source_time_value` computes instantaneous source amplitude per type.
- `inject_source_into_Ez` applies soft spatially weighted source injection.
- `find_nearest_source` finds active source nearest to mouse coordinates.
- `expr_compile` tokenizes, parses, and emits RPN program from expression string.
- `expr_eval` evaluates compiled expression with variables t/amp/freq/pi.
- `expr_free` frees compiled expression program memory.
- `ports_init` allocates port buffers, chooses default geometry, activates ports.
- `ports_apply_config` loads normalized port specs and clamps geometry.
- `ports_sample` integrates Ez/Hy over span into port buffers each step.
- `ports_free` frees port buffers and resets fields.
- `analysis_test_set_alloc_fail_after` configures failure injection counter.
- `scope_init` allocates oscilloscope ring buffer with min length 64.
- `scope_push` writes sample and updates rolling_absmax with decay.
- `scope_clear` zeros scope buffer and resets head.
- `scope_free` frees scope buffer and resets struct.
- `dump_scope_fft_csv` exports Hann-windowed DFT to CSV.
- `compute_s21` computes magnitude ratio of port1 vs port0 at given freq/dt.
- `probe_writer_open` opens buffered writer to file path.
- `probe_writer_append` appends timestep/value lines into buffer.
- `probe_writer_flush` writes buffered data to disk.
- `probe_writer_close` flushes and closes probe log file.
- `config_loader_parse_file` reads JSON file and overlays onto SimulationConfig.
- `config_load_from_args` parses CLI args, optional config file, validates/clamps.
- `config_validate` checks SimulationConfig for numerical and structural validity.
- `config_clamp_to_limits` enforces bounds and normalizes coordinates/counts.
- `config_print_summary` logs grid, sweep, boundary, probe, counts to stdout.
- `save_config_json` serializes SimulationConfig to JSON on disk.
- `simulation_bootstrap_from_args` loads config, validates, init sim, configures boundary.
- `simulation_bootstrap_with_config` initializes sim from provided config.
- `simulation_bootstrap_shutdown` frees simulation resources.
- `cli_runner_execute` chooses fixed vs sweep mode and runs simulation with logging.
- `run_fixed_mode` executes fixed step count, logs probe, prints progress.
- `run_sweep_mode` executes frequency sweep, computes S21, writes CSV.
- `simulation_runner_on_step` handles probe logging, port sampling, and progress output.
- `simulation_runner_flush` flushes probe writer buffers.
- `render_frame` orchestrates SDL rendering of field, UI panels, scopes, overlays.
- `render_layout_compute` calculates viewport, panels, colorbar, scopes, block outline.
- `render_field_heatmap` renders Ez heatmap for external callers with chosen colormap.
- `render_scope` draws scope buffer and returns observed max amplitude.
- `ui_state_init` allocates and seeds UIState defaults.
- `ui_handle_events` processes SDL events, shortcuts, painting, dragging, layout keys.
- `ui_update_metrics` updates autoscale estimates and smooths vmax metrics.
- `slider_handle_event` handles SDL slider interactions for freq and speed controls.
- `freq_from_slider` maps slider t to frequency log scale; `slider_from_freq` inverse.
- `render_layout_test` validates layout geometry consistency across sizes.
