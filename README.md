# emwave-c

Real-time 2D FDTD EM simulator (C + SDL2). Cross-platform via CMake.

## Key achievements
- Real-time 512x512 @ 60-120 FPS; 1024x768 @ 30-60 FPS; <16 ms UI latency.
- CPML with presets + manual tuning; reflections < -60 dB at target f.
- Instrumentation: dual probes, scope + FFT, Auto-P99/Peak/Hold, paint/edit tools.
- Portable Windows .exe (~1-2 MB), reproducible builds (Linux/Windows).

## Building

### Windows (MSVC + vcpkg)

```
cd C:\projects\emwave-c
.\build-msvc.ps1 -Clean   # optional
```

The script:
- Detects your Visual Studio install and sets up the MSVC environment
- Configures with CMake using the bundled vcpkg manifest (SDL2/SDL2_ttf/OpenMP)
- Builds `build\Release\emwave.exe`

Run the freshly built binary from `build\Release`.

### Linux / macOS

Install dependencies (Ubuntu example):

```
sudo apt-get update
sudo apt-get install -y cmake ninja-build libsdl2-dev libsdl2-ttf-dev
```

Configure and build:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable is `build/emwave`. On macOS you can use Homebrew packages with the same CMake invocation.

### Bundled font assets

The SDL GUI now ships with the [DejaVu Sans](https://dejavu-fonts.github.io/) typeface to avoid relying on host-specific font locations. The font binary lives in `third_party/fonts/DejaVuSans.ttf` together with its Bitstream Vera-derived license and is copied into `assets/fonts/DejaVuSans.ttf` during the build/install steps. At runtime the renderer always loads from that assets directory, so the UI is consistent across Windows, macOS, and Linux.

### Sample configurations

Runtime parameters, material layouts, and source placements can be scripted via
JSON. The schema is documented in [`configs/SCHEMA.md`](configs/SCHEMA.md).

Two ready-to-run examples are included:

```bash
./build/emwave --config configs/waveguide.json
./build/emwave --config configs/cpw_filter.json
```

Command-line overrides always win, so `--nx=768 --config=...` scales an existing
scene without editing the JSON file.

### CI status

GitHub Actions runs the same CMake/Ninja flow on `ubuntu-latest`, ensuring SDL2/SDL2_ttf + OpenMP builds succeed on each push/PR.
