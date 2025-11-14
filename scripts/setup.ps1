# ============================================================================
# emwave-c: Quick Setup Script for Windows (PowerShell)
# ============================================================================

Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "emwave-c Quick Setup for Windows" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

# Check for vcpkg
$vcpkgPaths = @(
    "C:\vcpkg\vcpkg.exe",
    "C:\tools\vcpkg\vcpkg.exe",
    "$env:USERPROFILE\vcpkg\vcpkg.exe",
    "$env:LOCALAPPDATA\vcpkg\vcpkg.exe"
)

$vcpkgExe = $null
foreach ($path in $vcpkgPaths) {
    if (Test-Path $path) {
        $vcpkgExe = $path
        $vcpkgRoot = Split-Path $vcpkgExe -Parent
        Write-Host "Found vcpkg at: $vcpkgRoot" -ForegroundColor Green
        break
    }
}

# Check PATH for vcpkg
if (-not $vcpkgExe) {
    $vcpkgExe = Get-Command vcpkg -ErrorAction SilentlyContinue
    if ($vcpkgExe) {
        $vcpkgRoot = Split-Path $vcpkgExe.Source -Parent
        Write-Host "Found vcpkg in PATH: $vcpkgRoot" -ForegroundColor Green
    }
}

if (-not $vcpkgExe) {
    Write-Host ""
    Write-Host "vcpkg not found!" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Option 1: Install vcpkg (Recommended)" -ForegroundColor White
    Write-Host "  1. Open PowerShell as Administrator" -ForegroundColor Gray
    Write-Host "  2. Run:" -ForegroundColor Gray
    Write-Host "     cd C:\" -ForegroundColor Cyan
    Write-Host "     git clone https://github.com/microsoft/vcpkg.git" -ForegroundColor Cyan
    Write-Host "     cd vcpkg" -ForegroundColor Cyan
    Write-Host "     .\bootstrap-vcpkg.bat" -ForegroundColor Cyan
    Write-Host "     .\vcpkg integrate install" -ForegroundColor Cyan
    Write-Host "     .\vcpkg install sdl2:x64-windows sdl2-ttf:x64-windows" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Option 2: Use MSYS2 instead" -ForegroundColor White
    Write-Host "  1. Install MSYS2 from: https://www.msys2.org/" -ForegroundColor Gray
    Write-Host "  2. Open 'MSYS2 MinGW 64-bit' terminal" -ForegroundColor Gray
    Write-Host "  3. Run:" -ForegroundColor Gray
    Write-Host "     pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf" -ForegroundColor Cyan
    Write-Host "     pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-toolchain" -ForegroundColor Cyan
    Write-Host ""

    $choice = Read-Host "Would you like to install vcpkg now? (y/n)"
    if ($choice -eq 'y' -or $choice -eq 'Y') {
        Write-Host ""
        Write-Host "Installing vcpkg to C:\vcpkg..." -ForegroundColor Cyan

        # Check for git
        $git = Get-Command git -ErrorAction SilentlyContinue
        if (-not $git) {
            Write-Host "ERROR: Git is required. Install from: https://git-scm.com/download/win" -ForegroundColor Red
            exit 1
        }

        # Clone vcpkg
        if (-not (Test-Path "C:\vcpkg")) {
            git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
            if ($LASTEXITCODE -ne 0) {
                Write-Host "ERROR: Failed to clone vcpkg" -ForegroundColor Red
                exit 1
            }
        }

        # Bootstrap
        Set-Location C:\vcpkg
        .\bootstrap-vcpkg.bat
        if ($LASTEXITCODE -ne 0) {
            Write-Host "ERROR: Failed to bootstrap vcpkg" -ForegroundColor Red
            exit 1
        }

        # Integrate
        .\vcpkg integrate install

        $vcpkgExe = "C:\vcpkg\vcpkg.exe"
        $vcpkgRoot = "C:\vcpkg"

        Write-Host "vcpkg installed successfully!" -ForegroundColor Green
    } else {
        Write-Host "Please install vcpkg or MSYS2 manually, then run this script again." -ForegroundColor Yellow
        exit 0
    }
}

# Install SDL2 dependencies
Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Installing SDL2 Dependencies" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

Set-Location $vcpkgRoot

Write-Host "Installing SDL2..." -ForegroundColor Cyan
& $vcpkgExe install sdl2:x64-windows
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to install SDL2" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "Installing SDL2_ttf..." -ForegroundColor Cyan
& $vcpkgExe install sdl2-ttf:x64-windows
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to install SDL2_ttf" -ForegroundColor Red
    exit 1
}

# Set environment variable
Write-Host ""
Write-Host "Setting VCPKG_ROOT environment variable..." -ForegroundColor Cyan
[System.Environment]::SetEnvironmentVariable('VCPKG_ROOT', $vcpkgRoot, [System.EnvironmentVariableTarget]::User)
$env:VCPKG_ROOT = $vcpkgRoot

Write-Host ""
Write-Host "=====================================" -ForegroundColor Green
Write-Host "Setup Complete!" -ForegroundColor Green
Write-Host "=====================================" -ForegroundColor Green
Write-Host ""
Write-Host "Dependencies installed:" -ForegroundColor White
Write-Host "  - SDL2 (x64-windows)" -ForegroundColor Gray
Write-Host "  - SDL2_ttf (x64-windows)" -ForegroundColor Gray
Write-Host ""
Write-Host "Environment:" -ForegroundColor White
Write-Host "  VCPKG_ROOT = $vcpkgRoot" -ForegroundColor Gray
Write-Host ""
Write-Host "Next steps:" -ForegroundColor White
Write-Host "  1. Restart PowerShell (to load VCPKG_ROOT)" -ForegroundColor Gray
Write-Host "  2. cd C:\projects\emwave-c" -ForegroundColor Cyan
Write-Host "  3. .\scripts\build_msvc.bat" -ForegroundColor Cyan
Write-Host ""
Write-Host "Or build manually:" -ForegroundColor White
Write-Host "  mkdir build" -ForegroundColor Cyan
Write-Host "  cd build" -ForegroundColor Cyan
Write-Host "  cmake .. -DCMAKE_TOOLCHAIN_FILE=$vcpkgRoot\scripts\buildsystems\vcpkg.cmake" -ForegroundColor Cyan
Write-Host "  cmake --build . --config Release" -ForegroundColor Cyan
Write-Host ""
