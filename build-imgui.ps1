# ============================================================================
# emwave-c: PowerShell Build Script for the Dear ImGui Front-End (MSVC)
# ============================================================================

param(
    [switch]$Clean,
    [ValidateSet('Debug', 'Release')]
    [string]$Config = 'Release'
)

Write-Host "=====================================" -ForegroundColor Magenta
Write-Host "emwave-c ImGui Build Script (PowerShell)" -ForegroundColor Magenta
Write-Host "=====================================" -ForegroundColor Magenta
Write-Host ""

$projectRoot = $PSScriptRoot
$buildDirName = "build-imgui"
$buildDir = Join-Path $projectRoot $buildDirName
Set-Location $projectRoot

Write-Host "Project root : $projectRoot" -ForegroundColor Gray
Write-Host "Build folder : $buildDirName" -ForegroundColor Gray
Write-Host "Build config : $Config" -ForegroundColor Gray
Write-Host ""

$envConfigured = $false

# Find Visual Studio installation so MSVC tools are available
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"

if (Test-Path $vswhere) {
    Write-Host "Detecting Visual Studio installation..." -ForegroundColor Cyan

    $vsPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath

    if ($vsPath) {
        Write-Host "Found Visual Studio at: $vsPath" -ForegroundColor Green

        # Prefer VsDevShell (PowerShell) to avoid long cmd lines
        $vsDevShell = Join-Path $vsPath "Common7\Tools\Launch-VsDevShell.ps1"
        if (Test-Path $vsDevShell) {
            Write-Host "Setting up MSVC environment (VsDevShell)..." -ForegroundColor Cyan
            try {
                & $vsDevShell -Arch x64 -HostArch x64 -SkipAutomaticLocation -NoLogo | Out-Null
                $envConfigured = $true
                Write-Host "MSVC environment configured" -ForegroundColor Green
                Write-Host ""
            } catch {
                Write-Host "VsDevShell failed, falling back to vcvarsall.bat" -ForegroundColor Yellow
            }
        }

        $vcvarsall = Join-Path $vsPath "VC\Auxiliary\Build\vcvarsall.bat"

        if (-not $envConfigured -and (Test-Path $vcvarsall)) {
            Write-Host "Setting up MSVC environment..." -ForegroundColor Cyan

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
            $envConfigured = $true
        }
    } else {
        Write-Host "WARNING: Visual Studio not found via vswhere" -ForegroundColor Yellow
    }
} else {
    Write-Host "WARNING: vswhere.exe not found" -ForegroundColor Yellow
}

# Verify CMake is present
$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if (-not $cmake) {
    Write-Host "ERROR: CMake not found in PATH" -ForegroundColor Red
    Write-Host "Install from https://cmake.org/download/ or add it via Visual Studio." -ForegroundColor Yellow
    exit 1
}

Write-Host "Using CMake: $($cmake.Source)" -ForegroundColor Gray
Write-Host ""

# Clean build directory if requested
if ($Clean -and (Test-Path $buildDir)) {
    Write-Host "Cleaning previous imgui build directory..." -ForegroundColor Cyan
    Remove-Item -Recurse -Force $buildDir
    Write-Host "Clean complete" -ForegroundColor Green
    Write-Host ""
}

if (-not (Test-Path $buildDir)) {
    New-Item -ItemType Directory -Path $buildDir | Out-Null
}

Set-Location $buildDir

# Locate vcpkg toolchain if present for SDL2 dependencies
$toolchainFile = $null
if ($env:VCPKG_ROOT) {
    $toolchainFile = Join-Path $env:VCPKG_ROOT "scripts\buildsystems\vcpkg.cmake"
    if (Test-Path $toolchainFile) {
        Write-Host "Using vcpkg toolchain: $toolchainFile" -ForegroundColor Gray
    } else {
        Write-Host "WARNING: VCPKG_ROOT is set but toolchain file not found" -ForegroundColor Yellow
        $toolchainFile = $null
    }
} else {
    Write-Host "INFO: VCPKG_ROOT not set; CMake will search for SDL2/SDL2_ttf via system paths" -ForegroundColor Yellow
    Write-Host "TIP: Set VCPKG_ROOT to a user-local vcpkg (not under Program Files) to avoid lock contention." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Configuring ImGui target with CMake..." -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

$cmakeArgs = @("..", "-DCMAKE_BUILD_TYPE=$Config", "-DEMWAVE_ENABLE_UI=ON")
if ($toolchainFile) {
    $cmakeArgs += "-DCMAKE_TOOLCHAIN_FILE=$toolchainFile"
}

$cl = Get-Command cl -ErrorAction SilentlyContinue
if ($cl) {
    Write-Host "Detected MSVC compiler, allowing CMake to select the default Visual Studio generator" -ForegroundColor Gray
} else {
    Write-Host "MSVC not detected in PATH; ensure vcvarsall.bat has been executed" -ForegroundColor Yellow
}

& cmake $cmakeArgs
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host "CMake configuration failed for emwave_imgui" -ForegroundColor Red
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host ""
    Write-Host "If SDL2/SDL2_ttf are missing run .\scripts\setup.ps1 or install them via vcpkg/MSYS2." -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Building emwave_imgui (target only)..." -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

cmake --build . --config $Config --target emwave_imgui
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host "emwave_imgui build failed" -ForegroundColor Red
    Write-Host "=====================================" -ForegroundColor Red
    Write-Host ""
    Write-Host "Ensure SDL2 + SDL2_ttf are discoverable and that the SDL UI was enabled during configure." -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Green
Write-Host "ImGui build successful!" -ForegroundColor Green
Write-Host "=====================================" -ForegroundColor Green
Write-Host ""

# Locate executable for convenience
$exePaths = @(
    ".\$Config\emwave_imgui.exe",
    ".\emwave_imgui.exe",
    "..\emwave_imgui.exe"
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
    Write-Host "Run instructions:" -ForegroundColor White
    Write-Host ("  cd {0}" -f (Split-Path $exePath -Parent)) -ForegroundColor Cyan
    Write-Host ("  .\{0}" -f (Split-Path $exePath -Leaf)) -ForegroundColor Cyan
} else {
    Write-Host "Executable expected at: $buildDirName\$Config\emwave_imgui.exe" -ForegroundColor Gray
}

Write-Host ""
