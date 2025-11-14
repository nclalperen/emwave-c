Param(
    [string[]]$CMakeArgs = @("-DCMAKE_BUILD_TYPE=Release")
)

$msysPaths = @("C:\msys64\ucrt64\bin", "C:\msys64\usr\bin")
$env:Path = ($msysPaths -join ";") + ";" + $env:Path

Write-Host "Using PATH override: $($msysPaths -join '; ')"

$configureArgs = @("-S", ".", "-B", "build", "-G", "Ninja") + $CMakeArgs
Write-Host "Running cmake $(($configureArgs -join ' '))"

cmake @configureArgs
if ($LASTEXITCODE -ne 0) {
    throw "Configuration failed"
}

Write-Host "Building..."
cmake --build build -j
if ($LASTEXITCODE -ne 0) {
    throw "Build failed"
}
