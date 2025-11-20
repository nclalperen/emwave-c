# emwave-c

Real-time 2D finite-difference time-domain (FDTD) electromagnetic simulator written in portable C with an SDL2 UI. Builds cleanly on Windows (MSVC, MSYS2), Linux, and macOS via CMake.

## Highlights
- 512x512 grids sustain 60-120 FPS on modern CPUs and remain interactive down to <16 ms UI latency.
- CPML and Mur boundary presets deliver <-60 dB reflections at the target frequency, plus knobs for manual tuning.
- Built-in instrumentation: dual probes, live scope/FFT panes, Auto-P99/Peak/Hold display modes, paint/edit tools, and CSV exports.
- Reproducible builds: MSVC flow uses the bundled vcpkg manifest (SDL2, SDL2_ttf, OpenMP); Linux/macOS use upstream packages with the same CMake configuration.
- Portable outputs: a self-contained Windows `emwave.exe` (~1-2 MB plus SDL DLLs) and a cross-platform `build/emwave`.

## Quick start

### Windows (Visual Studio + vcpkg)
```powershell
cd C:\projects\emwave-c
.\scripts\setup.ps1    # one-time dependency helper, optional if vcpkg already installed
.\build-msvc.ps1       # add -Clean for a fresh configure
.\build\Release\emwave.exe
```
The PowerShell wrapper locates Visual Studio, loads the MSVC environment, configures CMake with the local `vcpkg.json`, and produces `build\Release\emwave.exe` along with the required SDL DLLs and the bundled font.

Need the experimental Dear ImGui front-end? Run the dedicated helper to target `emwave_imgui` and keep its build outputs isolated:

```powershell
cd C:\projects\emwave-c
.\build-imgui.ps1           # supports -Clean / -Config Release|Debug
.\build-imgui\Release\emwave_imgui.exe
```
The script uses the same MSVC + vcpkg detection but only builds the `emwave_imgui` target so you can switch between the default SDL UI and the ImGui prototype without reconfiguring the primary build directory.

### Windows (MSYS2 / MinGW-w64)
```bash
pacman -Syu --needed base-devel mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-cmake mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf
cd /c/projects/emwave-c
./scripts/build_msys2.sh        # or run the cmake commands below manually
./build/emwave.exe
```
The helper script configures with `-G "MinGW Makefiles"` in Release mode and copies the font asset automatically.

### Linux / macOS
```bash
sudo apt-get install -y cmake ninja-build libsdl2-dev libsdl2-ttf-dev   # Ubuntu example
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/emwave
```
macOS users can substitute Homebrew (`brew install cmake ninja sdl2 sdl2_ttf`) with the same CMake invocation.

For deeper troubleshooting or alternative env setups, see `QUICK_START_WINDOWS.md` and `BUILDING_WINDOWS.md`.

## Detailed build notes

### Visual Studio (manual CMake flow)
```powershell
cd C:\projects\emwave-c
mkdir build
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DCMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake
cmake --build build --config Release
```
OpenMP is detected automatically; if you see `OpenMP not found` ensure the Desktop C++ workload is installed.

### MSYS2 manual commands
```bash
cd /c/projects/emwave-c
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/emwave.exe
```
Run the *MSYS2 MinGW 64-bit* shell (not the MSYS shell) so GCC, SDL2, and OpenMP headers resolve correctly.

### Linux/macOS options
- **Ninja (recommended):**
  ```bash
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
  cmake --build build
  ```
- **Makefiles:**
  ```bash
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build
  ```
Set `-DUSE_OPENMP=OFF` if building on a toolchain without OpenMP.

## Running tests
After configuring with CMake, run the unit test suite from the build directory:
```bash
cd build
ctest --output-on-failure          # single-config generators (Ninja/Make)
ctest -C Debug --output-on-failure # multi-config (Visual Studio/MSVC)
```
The tests cover configuration loading/validation, FDTD core allocation paths, and analysis/ports instrumentation. CI runs the same `ctest` invocation plus a headless smoke run of `emwave_cli` on every push and pull request.

## Headless CLI usage
The `emwave_cli` binary runs the solver without any SDL/UI dependencies so you can script sweeps and batch runs:
```bash
./build/emwave_cli --config configs/waveguide.json --run-steps=2000
./build/emwave_cli --config configs/cpw_filter.json --run-mode=sweep
```
Key flags:
- `--config=PATH` loads a JSON scene (see `configs/SCHEMA.md`).
- `--run-mode=fixed|sweep` selects a fixed-step run or an S-parameter sweep (default: fixed).
- `--run-steps=N` sets total steps when `run-mode=fixed` (otherwise falls back to sweep-derived count).
- `--boundary=cpml|mur` chooses the absorbing boundary model (default: CPML).
- `--profile` prints a short profiling summary (steps, wall time, steps/s).
- `--probe-log=PATH` enables probe logging to a text file; `--no-probe-log` disables it.
CLI flags always override values coming from the JSON config.

## Running the simulator
After building, launch from the build directory or project root:
```powershell
.\build\Release\emwave.exe     # Windows
```
```bash
./build/emwave                 # Linux/macOS/MSYS2
```

### Material blocks (ImGui front-end)
- Launch `build-imgui/Release/emwave_imgui.exe`, open the **Materials / Blocks** panel, and use `+ Add Block` for a default dielectric square or `+ Add Block from Material` to start from a library preset (Copper, FR4, etc.).
- Select a block to edit its material/type and the normalized coordinates `x0/y0/x1/y1`; the panel shows the corresponding meters using the current domain size, and enforces a minimum span of 0.01 in each direction.
- Materials drop-down comes from the Phase 2 library (11 presets). Switching materials swaps in the correct `epsr`, `sigma`, and PEC/PMC tags; dielectric fields expose editable `epsr`/`sigma`, while PEC/PMC hide them.
- Filtering options: pick a material in the browser, enable **Filter by selected material**, and optionally toggle **Auto-filter on selection** to focus the list on matching blocks. Blocks and paint can coexist; paint mode stays independent.
- Example microstrip stack (normalized coordinates): Copper ground (0.00,0.00)-(0.10,0.01), FR4 substrate (0.00,0.01)-(0.10,0.03), Copper trace (0.045,0.03)-(0.055,0.031). Save with `Ctrl+S` and reload to restore blocks.

### Configuration-driven runs
- JSON scenes live under `configs/`. The schema is documented in `configs/SCHEMA.md`.
- Override runtime parameters from the CLI; flags always win over JSON values.
  ```bash
  ./build/emwave --config configs/waveguide.json
  ./build/emwave --config configs/cpw_filter.json --nx=768 --ny=512 --tmax=2000
  ```
- Output artifacts:
  - `probe.txt` - probe samples for each timestep.
  - `probe_fft.csv` / `scope_fft.csv` - FFT/export triggered from the UI.
  - `sweep_s21.csv` - generated when running the S-parameter sweep tool.

### Controls at a glance
- `Space` pause/resume, `R` reset, `C` clear fields.
- `M` toggle Mur vs CPML boundaries.
- Mouse click/drag places or moves sources; numeric keys enter paint modes for PEC/PMC/dielectrics.
- `F` cycles source types (CW, Gaussian, Ricker). Arrow keys still nudge frequency when paused; the bottom timeline tray hosts sliders for frequency and steps/frame for precise control.
- `L` toggles the on-screen legend; `H` freezes the color scale; `J` freezes the scope scale; `S` toggles S-parameter ports.
- `[` `]` shrink/grow the left toolbox column, `,` `.` resize the right-side properties column.
- `-`/`=` change the height of the bottom timeline/scope tray, and `\` toggles whether the scope docks in that tray or lives in the properties column.
- `B` flips between dark and light UI themes; `V` (Shift+V reverses) cycles the neon accent palettes, and `K` (Shift+K reverses) cycles the field colormap (classic, viridis-style, plasma-style).
- `F1` shows the in-app onboarding overlay with a condensed keyboard cheat sheet and status badges.

## Assets and fonts
- The UI ships with DejaVu Sans (`third_party/fonts/DejaVuSans.ttf`). The build copies it into `assets/fonts/DejaVuSans.ttf` so runtime never relies on system fonts.
- Windows builds drop SDL2/SDL2_ttf DLLs plus the font next to the executable; copying `build/Release` (or the `dist` folder) to another machine keeps the app portable.

## Repository layout
- `src/` - simulation core (FDTD updates, CPML, material system), instrumentation, and SDL2 UI.
- `include/` - headers shared between modules.
- `configs/` - JSON scene templates plus the schema.
- `scripts/` - helper PowerShell/Bash scripts for setup and builds.
- `tests/` - regression checks for solver math and config parsing.
- `third_party/` - vendored assets such as fonts and small libraries.
- `build/` - generated binaries and intermediate files (ignored by Git).

## Additional documentation
- `QUICK_START_WINDOWS.md` - step-by-step PowerShell/MSYS2 instructions.
- `BUILDING_WINDOWS.md` - exhaustive Windows toolchain coverage.
- `MODULAR_ARCHITECTURE_SUMMARY.md` and `MODULAR_BUILD_STATUS.md` - deeper context on the code layout and refactors.
- `UI_OVERHAUL_GUIDE.md` - notes on potential UI replacements and future improvements.
- `SCENES.md` - short descriptions of the bundled example configs.
- `COMPILATION_STATUS.md` & `BUILD_SUCCESS.md` - recent validation reports.

emwave-c is actively developed; feel free to file issues or PRs with improvements, scenes, and instrumentation ideas.
