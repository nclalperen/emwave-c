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

### CI status

GitHub Actions runs the same CMake/Ninja flow on `ubuntu-latest`, ensuring SDL2/SDL2_ttf + OpenMP builds succeed on each push/PR.
