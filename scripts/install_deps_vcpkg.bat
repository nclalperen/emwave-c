@echo off
REM ============================================================================
REM emwave-c: Install Dependencies using vcpkg
REM ============================================================================

echo =====================================
echo emwave-c Dependency Installation
echo Using vcpkg for Windows
echo =====================================
echo.

REM Check if vcpkg is already installed
where vcpkg >nul 2>&1
if %errorlevel% equ 0 (
    echo Found vcpkg in PATH
    goto :install_deps
)

REM Check common vcpkg locations
set VCPKG_PATHS=C:\vcpkg C:\tools\vcpkg %USERPROFILE%\vcpkg %LOCALAPPDATA%\vcpkg

for %%p in (%VCPKG_PATHS%) do (
    if exist "%%p\vcpkg.exe" (
        echo Found vcpkg at %%p
        set VCPKG_ROOT=%%p
        set PATH=%%p;%PATH%
        goto :install_deps
    )
)

echo.
echo vcpkg not found. Installing vcpkg...
echo.

REM Install vcpkg to C:\vcpkg
set INSTALL_DIR=C:\vcpkg
echo Installing to %INSTALL_DIR%...
echo.

REM Check if git is available
where git >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Git is required to install vcpkg
    echo Please install Git from: https://git-scm.com/download/win
    echo Then run this script again.
    pause
    exit /b 1
)

REM Clone vcpkg
if not exist "%INSTALL_DIR%" (
    git clone https://github.com/microsoft/vcpkg.git "%INSTALL_DIR%"
    if %errorlevel% neq 0 (
        echo ERROR: Failed to clone vcpkg repository
        pause
        exit /b 1
    )
)

REM Bootstrap vcpkg
cd /d "%INSTALL_DIR%"
call bootstrap-vcpkg.bat
if %errorlevel% neq 0 (
    echo ERROR: Failed to bootstrap vcpkg
    pause
    exit /b 1
)

REM Integrate with Visual Studio
vcpkg integrate install
if %errorlevel% neq 0 (
    echo WARNING: vcpkg integration failed, continuing anyway...
)

set VCPKG_ROOT=%INSTALL_DIR%
set PATH=%INSTALL_DIR%;%PATH%

echo.
echo vcpkg installation complete!
echo.

:install_deps
echo =====================================
echo Installing SDL2 Dependencies
echo =====================================
echo.

echo Installing SDL2 for x64-windows...
vcpkg install sdl2:x64-windows

if %errorlevel% neq 0 (
    echo ERROR: Failed to install SDL2
    pause
    exit /b 1
)

echo.
echo Installing SDL2_ttf for x64-windows...
vcpkg install sdl2-ttf:x64-windows

if %errorlevel% neq 0 (
    echo ERROR: Failed to install SDL2_ttf
    pause
    exit /b 1
)

echo.
echo =====================================
echo Installation Complete!
echo =====================================
echo.
echo vcpkg root: %VCPKG_ROOT%
echo.
echo Next steps:
echo   1. Set VCPKG_ROOT environment variable (if not already set):
echo      setx VCPKG_ROOT "%VCPKG_ROOT%"
echo.
echo   2. Build the project:
echo      cd C:\projects\emwave-c
echo      scripts\build_msvc.bat
echo.

pause
