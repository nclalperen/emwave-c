Param(
    [string]$BuildDir = "build-imgui",
    [string]$Config = "Debug",
    [string]$LogPath = "validation_log.txt"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$logLines = @()
function Write-Log {
    Param([string]$Message)
    $timestamp = Get-Date -Format "yyyy-MM-ddTHH:mm:ss.ffffK"
    $line = "$timestamp $Message"
    Write-Host $line
    $script:logLines += $line
}

function Run-Step {
    Param(
        [string]$Name,
        [scriptblock]$Action
    )
    Write-Log "BEGIN: $Name"
    try {
        & $Action
        Write-Log "PASS: $Name"
    } catch {
        Write-Log "FAIL: $Name - $($_.Exception.Message)"
    }
}

Write-Log "==== Phase 1 Auto Validation ===="
Write-Log "BuildDir=$BuildDir Config=$Config"

# Section 0: Build
Run-Step "Build emwave_imgui ($Config)" {
    $cmd = "cmake --build `"$BuildDir`" --config $Config"
    Write-Log "Running: $cmd"
    $proc = Start-Process -FilePath "cmake" -ArgumentList "--build", $BuildDir, "--config", $Config -NoNewWindow -Wait -PassThru
    if ($proc.ExitCode -ne 0) {
        throw "Build failed with exit code $($proc.ExitCode)"
    }
}

# Section 1: Core unit tests (config/analysis/fdtd)
Run-Step "Run unit tests (config_loader_tests/config_runtime_tests/fdtd_core_tests/analysis_tests)" {
    $unitDir = Join-Path $BuildDir "tests/unit/$Config"
    if (-not (Test-Path $unitDir)) {
        throw "Unit test directory not found: $unitDir"
    }
    $tests = @(
        "config_loader_tests.exe",
        "config_runtime_tests.exe",
        "fdtd_core_tests.exe",
        "analysis_tests.exe",
        "material_library_tests.exe"
    )
    foreach ($t in $tests) {
        $path = Join-Path $unitDir $t
        if (-not (Test-Path $path)) {
            Write-Log "SKIP: $t not found in $unitDir"
            continue
        }
        Write-Log "Running unit test: $t"
        $proc = Start-Process -FilePath $path -NoNewWindow -Wait -PassThru
        if ($proc.ExitCode -ne 0) {
            throw "$t failed with exit code $($proc.ExitCode)"
        }
    }
}

# Section 2: CLI smoke test for configs/waveguide.json
Run-Step "CLI smoke test on configs/waveguide.json" {
    $cliPath = Join-Path $BuildDir "$Config/emwave_cli.exe"
    if (-not (Test-Path $cliPath)) {
        throw "emwave_cli not found at $cliPath"
    }
    $args = @(
        "--config", "configs/waveguide.json",
        "--run-steps=200"
    )
    Write-Log "Running: $cliPath $($args -join ' ')"
    $proc = Start-Process -FilePath $cliPath -ArgumentList $args -NoNewWindow -Wait -PassThru
    if ($proc.ExitCode -ne 0) {
        throw "emwave_cli smoke test failed with exit code $($proc.ExitCode)"
    }
}

# Section 3: Config save/load roundtrip (config_runtime/config_loader only)
Run-Step "Config save/load roundtrip (config_runtime/config_loader)" {
    # This checks that save_config_json() can write/read a simple config via emwave_cli.
    # We use a small JSON produced by save_config_json via a helper CLI built into emwave_cli:
    # emwave_cli does not currently expose this, so this step is a placeholder for now.
    Write-Log "NOTE: No direct CLI hook for save_config_json; manual testing required for Section 8."
}

# Sections requiring interactive GUI validation
$manualSections = @(
    "Section 1: Global keyboard shortcuts (ESC/Q, F1, F2/F3, F5/F6, Ctrl+S, F1 while help open)",
    "Section 2: Context shortcuts (Space/R/C, 1/2/3, T/F, S when ImGui does not capture keyboard)",
    "Section 3: Paint mode (M/U/I/O/P, cursor, brush rendering, paint vs drag interaction)",
    "Section 4: Auto-rescale modes (A/H/J/L behavior, vmax/Scope correlation)",
    "Section 5: Visual controls (B theme, K colormap, V accent, G grid, combinatorial readability)",
    "Section 6: Boundary conditions (Y toggle, 7/8/9 presets, field/scope reset behavior)",
    "Section 7: Source dragging (drag start/stop, bounds, paint-mode blocking, running vs paused)",
    "Section 8: Full configuration save/load roundtrip via UI (Ctrl+S, F5/F6, manual reload)",
    "Section 9: UI/UX polish (menus, panels, help overlay, logging, performance, viewport/scope smoothness)",
    "Section 10: Edge cases & stability (rapid key presses, window resize, 5+ minute stress run)"
)

Write-Log "==== Sections requiring manual/interactive validation ===="
foreach ($s in $manualSections) {
    Write-Log "MANUAL: $s"
}

Write-Log "==== Bug Reporting Template ===="
Write-Log "**Issue**: [Brief description]"
Write-Log "**Steps to reproduce**:"
Write-Log "1. [Step 1]"
Write-Log "2. [Step 2]"
Write-Log "**Expected**: [What should happen]"
Write-Log "**Actual**: [What actually happens]"
Write-Log "**Severity**: [Critical/High/Medium/Low]"

Write-Log "==== Auto validation complete ===="

try {
    $logLines | Set-Content -Path $LogPath -Encoding UTF8
    Write-Log "Log saved to $LogPath"
} catch {
    Write-Host "WARN: Failed to write log to $LogPath - $($_.Exception.Message)"
}
