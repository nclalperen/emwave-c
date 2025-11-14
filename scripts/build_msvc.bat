@echo off
REM ============================================================================
REM emwave-c: Build Script for MSVC on Windows
REM ============================================================================

echo =====================================
echo emwave-c Build Script (MSVC)
echo =====================================
echo.

REM Check if running in Visual Studio environment
where cl.exe >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: MSVC compiler not found in PATH
    echo Please run this script from:
    echo   - "x64 Native Tools Command Prompt for VS"
    echo   - Or "Developer Command Prompt for VS"
    echo.
    echo Alternatively, you can manually run vcvarsall.bat:
    echo   "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
    exit /b 1
)

REM Get project root (script is in scripts\)
cd /d "%~dp0\.."
set PROJECT_ROOT=%cd%
echo Project root: %PROJECT_ROOT%
echo.

REM Clean previous build if requested
if "%1"=="clean" (
    echo Cleaning previous build...
    if exist build rmdir /s /q build
    echo Clean complete
    echo.
)

REM Create build directory
if not exist build mkdir build
cd build

REM Check for vcpkg toolchain
if defined VCPKG_ROOT (
    set TOOLCHAIN_FILE=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake
    echo Using vcpkg toolchain: %TOOLCHAIN_FILE%
    echo.
) else (
    echo WARNING: VCPKG_ROOT not set
    echo If CMake cannot find SDL2, please:
    echo   1. Install vcpkg
    echo   2. Set VCPKG_ROOT environment variable
    echo   3. Install SDL2: vcpkg install sdl2:x64-windows sdl2-ttf:x64-windows
    echo.
    set TOOLCHAIN_FILE=
)

REM Configure
echo Configuring with CMake...
if defined TOOLCHAIN_FILE (
    cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="%TOOLCHAIN_FILE%"
) else (
    cmake .. -G "Visual Studio 17 2022" -A x64
)

if %errorlevel% neq 0 (
    echo.
    echo ERROR: CMake configuration failed
    echo.
    echo Make sure you have:
    echo   1. SDL2 and SDL2_ttf installed via vcpkg
    echo   2. VCPKG_ROOT environment variable set
    echo   3. vcpkg integrate install has been run
    exit /b 1
)

echo.
echo Building...
cmake --build . --config Release

if %errorlevel% neq 0 (
    echo.
    echo ERROR: Build failed
    exit /b 1
)

echo.
echo =====================================
echo Build complete!
echo =====================================
echo.
echo Executable: build\Release\emwave.exe
echo.
echo To run:
echo   cd build\Release
echo   emwave.exe
echo.
