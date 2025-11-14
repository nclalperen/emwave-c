# ============================================================================
# emwave-c: PowerShell Build Script with Automatic MSVC Detection
# ============================================================================

param(
    [switch]$Clean,
    [ValidateSet('Debug', 'Release')]
    [string]$Config = 'Release'
)

Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "emwave-c Build Script (PowerShell)" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

$projectRoot = $PSScriptRoot
Set-Location $projectRoot

Write-Host "Project root: $projectRoot" -ForegroundColor Gray
Write-Host "Build config: $Config" -ForegroundColor Gray
Write-Host ""

# Find Visual Studio installation
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"

if (Test-Path $vswhere) {
    Write-Host "Detecting Visual Studio installation..." -ForegroundColor Cyan

    $vsPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath

    if ($vsPath) {
        Write-Host "Found Visual Studio at: $vsPath" -ForegroundColor Green

        # Set up MSVC environment
        $vcvarsall = Join-Path $vsPath "VC\Auxiliary\Build\vcvarsall.bat"

        if (Test-Path $vcvarsall) {
            Write-Host "Setting up MSVC environment..." -ForegroundColor Cyan

            # Import Visual Studio environment variables
            $tempFile = [IO.Path]::GetTempFileName()
            cmd /c "`"$vcvarsall`" x64 && set" > $tempFile

            Get-Content $tempFile | ForEach-Object {
                if ($_ -match '^([^=]+)=(.*)$') {
                    [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2])
                }
            }
            Remove-Item $tempFile

            Write-Host "MSVC environment configured" -ForegroundColor Green
            Write-Host ""
        }
    } else {
        Write-Host "WARNING: Visual Studio not found via vswhere" -ForegroundColor Yellow
    }
} else {
    Write-Host "WARNING: vswhere.exe not found" -ForegroundColor Yellow
}

# Check for CMake
$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if (-not $cmake) {
    Write-Host "ERROR: CMake not found in PATH" -ForegroundColor Red
    Write-Host ""
    Write-Host "Please install CMake from: https://cmake.org/download/" -ForegroundColor Yellow
    Write-Host "Or install Visual Studio with 'Desktop development with C++' workload" -ForegroundColor Yellow
    exit 1
}

Write-Host "Using CMake: $($cmake.Source)" -ForegroundColor Gray
Write-Host ""

# Clean if requested
if ($Clean) {
    Write-Host "Cleaning previous build..." -ForegroundColor Cyan
    if (Test-Path "build") {
        Remove-Item -Recurse -Force "build"
        Write-Host "Clean complete" -ForegroundColor Green
    }
    Write-Host ""
}

# Create build directory
if (-not (Test-Path "build")) {
    New-Item -ItemType Directory -Path "build" | Out-Null
}

Set-Location "build"

# Configure CMake
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Configuring with CMake..." -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

# Check for vcpkg
$toolchainFile = $null
if ($env:VCPKG_ROOT) {
    $toolchainFile = Join-Path $env:VCPKG_ROOT "scripts\buildsystems\vcpkg.cmake"
    if (Test-Path $toolchainFile) {
        Write-Host "Using vcpkg toolchain: $toolchainFile" -ForegroundColor Gray
    } else {
        Write-Host "WARNING: VCPKG_ROOT set but toolchain file not found" -ForegroundColor Yellow
        $toolchainFile = $null
    }
} else {
    Write-Host "INFO: VCPKG_ROOT not set, attempting to find SDL2 via system paths" -ForegroundColor Gray
}

Write-Host ""

# Run CMake configure
$cmakeArgs = @("..", "-DCMAKE_BUILD_TYPE=$Config")

if ($toolchainFile) {
    $cmakeArgs += "-DCMAKE_TOOLCHAIN_FILE=$toolchainFile"
}

# Try to detect generator
$cl = Get-Command cl -ErrorAction SilentlyContinue
if ($cl) {
    Write-Host "Detected MSVC compiler, using Visual Studio generator" -ForegroundColor Gray
    # Let CMake auto-detect Visual Studio version
} else {
    Write-Host "MSVC not in PATH, using default generator" -ForegroundColor Gray
}

& cmake $cmakeArgs

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host "CMake Configuration Failed" -ForegroundColor Red
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host ""

    if (-not $env:VCPKG_ROOT) {
        Write-Host "SDL2 not found. Please install dependencies:" -ForegroundColor Yellow
        Write-Host ""
        Write-Host "Option 1: Run the setup script:" -ForegroundColor White
        Write-Host "  .\scripts\setup.ps1" -ForegroundColor Cyan
        Write-Host ""
        Write-Host "Option 2: Install manually with vcpkg:" -ForegroundColor White
        Write-Host "  cd C:\" -ForegroundColor Cyan
        Write-Host "  git clone https://github.com/microsoft/vcpkg.git" -ForegroundColor Cyan
        Write-Host "  cd vcpkg" -ForegroundColor Cyan
        Write-Host "  .\bootstrap-vcpkg.bat" -ForegroundColor Cyan
        Write-Host "  .\vcpkg install sdl2:x64-windows sdl2-ttf:x64-windows" -ForegroundColor Cyan
        Write-Host "  setx VCPKG_ROOT C:\vcpkg" -ForegroundColor Cyan
        Write-Host ""
    }

    exit 1
}

# Build
Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Building..." -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

cmake --build . --config $Config

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host "Build Failed" -ForegroundColor Red
    Write-Host "=====================================" -ForegroundColor Red
    exit 1
}

# Success
Write-Host ""
Write-Host "=====================================" -ForegroundColor Green
Write-Host "Build Successful!" -ForegroundColor Green
Write-Host "=====================================" -ForegroundColor Green
Write-Host ""

# Find executable
$exePaths = @(
    ".\emwave.exe",
    ".\$Config\emwave.exe",
    "..\emwave.exe"
)

$exePath = $null
foreach ($path in $exePaths) {
    if (Test-Path $path) {
        $exePath = (Resolve-Path $path).Path
        break
    }
}

if ($exePath) {
    Write-Host "Executable: $exePath" -ForegroundColor White
    Write-Host ""
    Write-Host "To run:" -ForegroundColor White
    Write-Host "  cd $(Split-Path $exePath -Parent)" -ForegroundColor Cyan
    Write-Host "  .\$(Split-Path $exePath -Leaf)" -ForegroundColor Cyan
} else {
    Write-Host "Executable location:" -ForegroundColor White
    Write-Host "  build\$Config\emwave.exe" -ForegroundColor Gray
}

Write-Host ""
Write-Host "OpenMP Status:" -ForegroundColor White
if (Test-Path "CMakeCache.txt") {
    $openmpInfo = Select-String -Path "CMakeCache.txt" -Pattern "OpenMP" | Select-Object -First 3
    if ($openmpInfo) {
        $openmpInfo | ForEach-Object { Write-Host "  $($_.Line)" -ForegroundColor Gray }
    } else {
        Write-Host "  (OpenMP status not found)" -ForegroundColor Gray
    }
}

Write-Host ""
