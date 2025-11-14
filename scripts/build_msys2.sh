#!/bin/bash
# ============================================================================
# emwave-c: Build Script for MSYS2/MinGW on Windows
# ============================================================================

set -e  # Exit on any error

echo "====================================="
echo "emwave-c Build Script (MSYS2/MinGW)"
echo "====================================="
echo ""

# Check for required tools
echo "Checking build environment..."
if ! command -v cmake &> /dev/null; then
    echo "ERROR: CMake not found. Install with:"
    echo "  pacman -S mingw-w64-x86_64-cmake"
    exit 1
fi

if ! command -v gcc &> /dev/null; then
    echo "ERROR: GCC not found. Install with:"
    echo "  pacman -S mingw-w64-x86_64-toolchain"
    exit 1
fi

# Check for SDL2
if ! pkg-config --exists sdl2; then
    echo "ERROR: SDL2 not found. Install with:"
    echo "  pacman -S mingw-w64-x86_64-SDL2"
    echo "  pacman -S mingw-w64-x86_64-SDL2_ttf"
    exit 1
fi

echo "Build environment OK"
echo ""

# Get project root (script is in scripts/)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "Project root: $PROJECT_ROOT"
echo ""

# Clean previous build if requested
if [ "$1" == "clean" ]; then
    echo "Cleaning previous build..."
    rm -rf build
    echo "Clean complete"
    echo ""
fi

# Create build directory
mkdir -p build
cd build

# Configure
echo "Configuring with CMake..."
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release

echo ""
echo "Building..."
cmake --build . --config Release

echo ""
echo "====================================="
echo "Build complete!"
echo "====================================="
echo ""
echo "Executable: build/emwave.exe"
echo ""
echo "To run:"
echo "  cd build"
echo "  ./emwave.exe"
echo ""
echo "OpenMP status:"
if [ -f CMakeCache.txt ]; then
    grep "OpenMP" CMakeCache.txt | head -3 || echo "  (OpenMP status not found in cache)"
fi
