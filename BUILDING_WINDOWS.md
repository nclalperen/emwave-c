# Building emwave-c on Windows

This guide provides step-by-step instructions for building the emwave-c FDTD electromagnetic simulator on Windows.

## Prerequisites

### Option 1: MSYS2/MinGW (Recommended)

1. **Install MSYS2** from https://www.msys2.org/
   - Download and run the installer
   - Follow the installation wizard (default location: `C:\msys64`)

2. **Update MSYS2 packages**
   ```bash
   pacman -Syu
   # Close terminal when prompted, then reopen and run:
   pacman -Su
   ```

3. **Install build tools and dependencies**
   ```bash
   pacman -S --needed base-devel mingw-w64-x86_64-toolchain
   pacman -S mingw-w64-x86_64-cmake
   pacman -S mingw-w64-x86_64-SDL2
   pacman -S mingw-w64-x86_64-SDL2_ttf
   ```

4. **Add MSYS2 to PATH**
   - Add `C:\msys64\mingw64\bin` to your system PATH
   - Or use MSYS2 MinGW 64-bit shell for building

### Option 2: Visual Studio (MSVC)

1. **Install Visual Studio 2019 or later**
   - Download from https://visualstudio.microsoft.com/
   - Select "Desktop development with C++" workload
   - Include CMake tools for Windows

2. **Install vcpkg for dependencies**
   ```cmd
   git clone https://github.com/Microsoft/vcpkg.git
   cd vcpkg
   .\bootstrap-vcpkg.bat
   .\vcpkg integrate install
   ```

3. **Install SDL2 libraries**
   ```cmd
   .\vcpkg install sdl2:x64-windows
   .\vcpkg install sdl2-ttf:x64-windows
   ```

## Building with MSYS2/MinGW

1. **Open MSYS2 MinGW 64-bit terminal**

2. **Navigate to project directory**
   ```bash
   cd /c/projects/emwave-c
   ```

3. **Create build directory**
   ```bash
   mkdir -p build
   cd build
   ```

4. **Configure with CMake**
   ```bash
   cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
   ```

5. **Build the project**
   ```bash
   cmake --build . --config Release
   ```

6. **Run the executable**
   ```bash
   ./emwave.exe
   ```

## Building with Visual Studio (MSVC)

### Method 1: Using CMake GUI

1. **Open CMake GUI**
   - Set source directory to `C:\projects\emwave-c`
   - Set build directory to `C:\projects\emwave-c\build`

2. **Configure**
   - Click "Configure"
   - Select "Visual Studio 16 2019" (or your version)
   - Platform: x64
   - Specify vcpkg toolchain file if prompted:
     ```
     C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
     ```

3. **Generate**
   - Click "Generate" to create Visual Studio solution

4. **Build**
   - Open `emwave.sln` in Visual Studio
   - Build → Build Solution (F7)
   - Or use command line:
     ```cmd
     cmake --build build --config Release
     ```

### Method 2: Using Visual Studio Developer Command Prompt

1. **Open "x64 Native Tools Command Prompt for VS"**

2. **Navigate to project**
   ```cmd
   cd C:\projects\emwave-c
   ```

3. **Configure and build**
   ```cmd
   mkdir build
   cd build
   cmake .. -G "Visual Studio 16 2019" -A x64 ^
     -DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
   cmake --build . --config Release
   ```

## OpenMP Support

### MSYS2/MinGW
OpenMP is included by default with GCC. CMake will automatically detect it.

### MSVC
OpenMP is included with Visual Studio. CMake will automatically enable it with the `/openmp` flag.

To verify OpenMP is enabled, check the CMake output:
```
-- OpenMP found - parallel processing enabled
```

If you see:
```
-- OpenMP not found - sequential processing only
```

Then ensure you have the OpenMP libraries installed for your toolchain.

## Troubleshooting

### SDL2 not found (MSYS2)

**Error:** `Could NOT find SDL2 (missing: SDL2_LIBRARY SDL2_INCLUDE_DIR)`

**Solution:**
```bash
pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf
```

### SDL2 not found (MSVC)

**Error:** `Could NOT find SDL2`

**Solution:** Ensure vcpkg toolchain file is specified:
```cmd
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Missing DejaVuSans.ttf

**Error:** `TTF_OpenFont failed`

**Solution:** The build automatically copies `arial.ttf` from Windows fonts directory. If this fails:
1. Check that `C:\Windows\Fonts\arial.ttf` exists
2. Or manually copy any TTF font to the executable directory and rename to `DejaVuSans.ttf`

### OpenMP not found

**MSYS2:** Install OpenMP support:
```bash
pacman -S mingw-w64-x86_64-openmp
```

**MSVC:** Ensure you selected "Desktop development with C++" workload in Visual Studio installer.

### Permission denied on arial.ttf copy

**Symptom:** Build succeeds but font copy fails

**Solution:** Run build as administrator, or manually copy the font after build:
```cmd
copy C:\Windows\Fonts\arial.ttf build\Release\DejaVuSans.ttf
```

## Performance Optimization

### Release Build (Recommended)
Always build in Release mode for optimal performance:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Release mode enables:
- Compiler optimizations (-O3 or /O2)
- OpenMP parallelization
- Inline function expansion
- Fast math operations

Debug builds will run 5-10x slower and should only be used for development.

### OpenMP Thread Count
By default, OpenMP uses all available CPU cores. To limit threads:

**Windows (cmd):**
```cmd
set OMP_NUM_THREADS=4
emwave.exe
```

**MSYS2 (bash):**
```bash
export OMP_NUM_THREADS=4
./emwave.exe
```

## File Structure After Build

```
emwave-c/
├── build/
│   ├── emwave.exe          # Executable (MinGW)
│   ├── DejaVuSans.ttf      # Auto-copied font
│   └── Release/            # MSVC build output
│       ├── emwave.exe
│       └── DejaVuSans.ttf
├── include/                # Header files
├── src/                    # Source files
└── CMakeLists.txt
```

## Next Steps

After successful build:

1. **Run the simulation**
   ```bash
   cd build
   ./emwave.exe
   ```

2. **Test basic functionality**
   - Space: Pause/Resume
   - Click: Add sources
   - R: Reset simulation
   - M: Toggle Mur/CPML boundaries

3. **Check performance**
   - Look for "OpenMP enabled with N threads" message
   - Monitor FPS in the info panel
   - Expected: 30-60 FPS on modern CPUs

4. **Review probe output**
   - Check `probe.txt` for field measurements
   - Use this for validation against original version

## Modular Architecture

The current build uses the new modular architecture:

**Core simulation modules** (implemented):
- `src/fdtd_core.c` - FDTD engine
- `src/boundary.c` - CPML/Mur boundaries
- `src/sources.c` - Wave sources
- `src/materials.c` - Material properties
- `src/analysis.c` - Measurement tools

**UI modules** (deferred):
- Still using monolithic `src/main.c` with embedded UI
- Future: Will be extracted to `ui_render.c` and `ui_controls.c`

All Phase 1 fixes (CPML indexing, buffer safety, NULL checks) and Phase 2 optimizations (OpenMP, caching) are included.
