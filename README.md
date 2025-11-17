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

## Running the simulator
After building, launch from the build directory or project root:
```powershell
.\build\Release\emwave.exe     # Windows
```
```bash
./build/emwave                 # Linux/macOS/MSYS2
```

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
- `L` toggles the on-screen legend; `H` freezes the color scale; `S` freezes the scope scale.
- `[` `]` shrink/grow the left toolbox column, `,` `.` resize the right-side properties column.
- `-`/`=` change the height of the bottom timeline/scope tray, and `\` toggles whether the scope docks in that tray or lives in the properties column.
- `B` flips between dark and light UI themes; `V` (Shift+V reverses) cycles the neon accent palettes if you want a different mood for the highlighted data.

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
- `COMPILATION_STATUS.md` & `BUILD_SUCCESS.md` - recent validation reports.

emwave-c is actively developed; feel free to file issues or PRs with improvements, scenes, and instrumentation ideas.
