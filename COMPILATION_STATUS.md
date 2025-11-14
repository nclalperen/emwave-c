# Compilation Status - MSVC on Windows

> **Status update:** The modular SDL front-end (`src/main_new.c`) is now the only
> UI entry point that participates in the build. The legacy monolithic
> implementation has been archived at
> [src/legacy/main_monolithic.c](src/legacy/main_monolithic.c) for historical
> reference and is no longer compiled. The notes below are preserved for context
> on the retired file.

## Current State

### ✅ Completed Fixes

1. **vcpkg manifest mode configured** - SDL2 and SDL2_ttf successfully installed via vcpkg.json
2. **Modular source files fixed for MSVC**:
   - [src/fdtd_core.c](src/fdtd_core.c) - OpenMP loops fixed (pre-declared variables)
   - [src/boundary.c](src/boundary.c) - OpenMP loops fixed
   - [include/config.h](include/config.h) - Changed `BASE_DX`/`BASE_DY` to macros (MSVC doesn't allow non-constant static initializers)

3. **CMakeLists.txt MSVC optimizations**:
   - Added `_CRT_SECURE_NO_WARNINGS` to disable fopen deprecation warnings
   - Added `/wd4849` to suppress OpenMP collapse clause warnings
   - Added `/openmp:experimental` for better OpenMP support in MSVC

### ⚠️ Remaining Issues in main.c

The monolithic [src/main.c](src/main.c) still has compilation errors because it contains:

1. **OpenMP loop variable declaration issues** - Multiple loops use C99-style `for (int i=...)` inside `#pragma omp` directives
2. **OpenMP collapse(2) and reduction directives** - MSVC OpenMP 2.0 doesn't fully support these
3. **OpenMP max reduction** - Line 1928 uses `reduction(max:vmax)` which requires `/openmp:llvm` flag

## Recommended Solutions

### Option 1: Quick Fix - Disable OpenMP in main.c for MSVC

Add conditional compilation around OpenMP pragmas in main.c:

```c
// Replace:
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
for (int i=0;i<NX;++i) for (int j=0;j<NY;++j) Ez_old[i][j]=Ez[i][j];

// With:
#if defined(_OPENMP) && !defined(_MSC_VER)
#pragma omp parallel for collapse(2)
#endif
for (int i=0;i<NX;++i) for (int j=0;j<NY;++j) Ez_old[i][j]=Ez[i][j];
```

This disables OpenMP only in main.c for MSVC, while keeping it enabled in the modular files (which are already fixed).

### Option 2: Fix All OpenMP Loops in main.c (More Work)

Pre-declare loop variables and remove collapse(2) directives throughout main.c:

**Lines requiring fixes:**
- Line 1716-1718: Ez_old save loop
- Line 1722-1724: H field update loop
- Line 1763-1765: E field update loop
- Line 1796-1798: Mur boundary loop (bottom/top)
- Line 1810-1812: Mur boundary loop (left/right)
- Line 1928-1930: vmax reduction loop (needs `/openmp:llvm` for reduction)
- Line 1946-1948: Histogram loop

**Example fix pattern:**
```c
// Before:
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
for (int i=0;i<NX;++i) for (int j=0;j<NY;++j){
    // ...
}

// After:
int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
for (i=0; i<NX; ++i) {
    for (j=0; j<NY; ++j){
        // ...
    }
}
```

### Option 3: Use Modular main_new.c (Best Long-Term)

Extract UI code from main.c and use the clean [src/main_new.c](src/main_new.c) template that calls the already-fixed modular functions. This avoids all OpenMP issues in main.c entirely.

## MSVC-Specific OpenMP Limitations

**MSVC OpenMP 2.0 does NOT support:**
- `collapse(n)` clause on parallel for loops (warning C4849)
- `reduction(max:var)` without `/openmp:llvm` (error C7660)
- C99-style loop variable declarations in OpenMP loops (error C3015)

**Workarounds:**
- Remove `collapse(2)` and use `private(j)` instead
- Add `/openmp:llvm` flag for reduction support
- Pre-declare all loop variables before `#pragma omp`

## Files Status

### ✅ Ready to Compile (MSVC-compatible)
- [src/fdtd_core.c](src/fdtd_core.c)
- [src/boundary.c](src/boundary.c)
- [src/sources.c](src/sources.c)
- [src/materials.c](src/materials.c)
- [src/analysis.c](src/analysis.c)
- [include/*.h](include/) (all headers)

### ❌ Needs Fixes (main.c)
- [src/main.c](src/main.c) - Lines 1716, 1722, 1763, 1796, 1810, 1928, 1946

## Testing After Fix

Once main.c is fixed, the build should complete successfully:

```powershell
cd C:\projects\emwave-c
.\build-msvc.ps1 -Clean
```

Expected output:
```
✅ SDL2 installed via vcpkg
✅ OpenMP found - parallel processing enabled
✅ All source files compile
✅ Executable: build\Release\emwave.exe
```

## Current Build Command

```powershell
cd C:\projects\emwave-c
.\build-msvc.ps1 -Clean
```

This will:
1. Clean previous build
2. Run vcpkg to install SDL2/SDL2_ttf (already done)
3. Configure CMake with MSVC
4. Attempt to build (currently fails on main.c OpenMP issues)

## Quick Command to Apply Option 1 Fix

To disable OpenMP in main.c for MSVC only, search and replace in main.c:

```powershell
# In PowerShell:
(Get-Content src\main.c) -replace '#ifdef _OPENMP', '#if defined(_OPENMP) && !defined(_MSC_VER)' | Set-Content src\main.c.new
Move-Item -Force src\main.c.new src\main.c
```

Then rebuild:
```powershell
.\build-msvc.ps1 -Clean
```

## Summary

**Current blocker:** src/main.c has MSVC-incompatible OpenMP syntax

**Quickest fix:** Disable OpenMP in main.c for MSVC (`#if defined(_OPENMP) && !defined(_MSC_VER)`)

**Best fix:** Complete UI extraction to use modular architecture (already 80% done)

All modular simulation files are already MSVC-compatible and will compile without errors.
