# Smoke test for Print Composer headless export.
# Generates a single PNG frame from the sample layout via --ch-layout.

param(
    [string]$OutputName = "smoke_layout",
    [int]$Fps = 30,
    [int]$Frames = 1
)

$root = Split-Path -Parent $PSScriptRoot
$exe = Join-Path $root "build-imgui\\Release\\emwave_imgui.exe"
$layout = Join-Path $root "recordings\\sample_layout.json"
$workdir = Split-Path $exe
$recDir = Join-Path $workdir "recordings"
$outFile = Join-Path $recDir "$OutputName.png"

if (-not (Test-Path $exe)) {
    Write-Error "Missing emwave_imgui.exe at $exe. Build ImGui first (build-imgui.ps1)."
    exit 1
}
if (-not (Test-Path $layout)) {
    Write-Error "Missing layout file at $layout. Run from repo root; sample_layout.json was added under recordings/."
    exit 1
}

if (-not (Test-Path $recDir)) { New-Item -ItemType Directory -Path $recDir -Force | Out-Null }
if (Test-Path $outFile) { Remove-Item $outFile -Force }

Push-Location $workdir
& $exe --composer-headless --ch-layout="$layout" --ch-output="$OutputName" --ch-format=png --ch-frames=$Frames --ch-fps=$Fps
$rc = $LASTEXITCODE
Pop-Location

if ($rc -ne 0) {
    Write-Error "Composer headless export failed (rc=$rc)."
    exit $rc
}

function Test-ImageVariance {
    param([string]$Path)
    Add-Type -AssemblyName System.Drawing
    $bmp = New-Object System.Drawing.Bitmap($Path)
    try {
        $w = $bmp.Width
        $h = $bmp.Height
        $sx = [math]::Max(1, [math]::Floor($w / 64))
        $sy = [math]::Max(1, [math]::Floor($h / 64))
        $p0 = $bmp.GetPixel(0,0)
        for ($y = 0; $y -lt $h; $y += $sy) {
            for ($x = 0; $x -lt $w; $x += $sx) {
                $p = $bmp.GetPixel($x,$y)
                if ($p.ToArgb() -ne $p0.ToArgb()) { return $true }
            }
        }
        return $false
    } finally {
        $bmp.Dispose()
    }
}

$targetPath = $outFile
if ($Frames -gt 1) {
    $frameDir = Join-Path $recDir ($OutputName + "_frames")
    $firstFrame = Join-Path $frameDir "frame_0000.png"
    if (-not (Test-Path $firstFrame)) {
        Write-Error "Expected first frame not found: $firstFrame"
        exit 1
    }
    $targetPath = $firstFrame
}

if (-not (Test-Path $targetPath)) {
    Write-Error "Expected output PNG not found: $targetPath"
    exit 1
}
$info = Get-Item $targetPath
if ($info.Length -le 0) {
    Write-Error "Output PNG is empty: $targetPath"
    exit 1
}
$var = Test-ImageVariance -Path $targetPath
if (-not $var) {
    Write-Error "Output PNG appears uniform (no color variance): $targetPath"
    exit 1
}

Write-Host "Smoke test ok. Output: $targetPath ($([math]::Round($info.Length/1KB,2)) KB)"
