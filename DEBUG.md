# Debugging and diagnostics

This document summarizes the knobs that are useful when debugging emwave‑c or
profiling headless runs.

## Build‑time options

### Bounds checking on field access

The core field accessors in `include/fdtd_core.h` can perform cheap bounds
checks in debug builds:

- Define `EMWAVE_BOUNDS_CHECK` at compile time to enable:
  - `fdtd_get_Ez/Hx/Hy` return `0.0` if `i`/`j` are out of range or the field
    pointers are null.
  - In release builds this macro is normally **off** to avoid overhead.

You can enable the checks from CMake:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DEM_WAVE_ENABLE_BOUNDS_CHECK=ON
cmake --build build
```

All targets that use the FDTD core (`emwave_core`, `emwave`, `emwave_cli`, and
the unit tests) will see `EMWAVE_BOUNDS_CHECK=1` and use the guarded helpers.

## Runtime flags (CLI)

The headless CLI entry, `emwave_cli`, accepts several flags that help with
diagnostics:

- `--profile`
  - Enables a short summary at the end of a run:
    - `Profile: <steps> steps in <secs> s (<steps/s> steps/s)`
  - Works for both fixed‑step and sweep runs.

- `--probe-log=PATH` / `--no-probe-log`
  - Enables or disables logging of the probe time series to `PATH`
    (defaults to `probe.txt`).

- `--boundary=cpml|mur`
  - Switches between CPML and Mur boundary conditions for headless tests.

For more CLI details, see the “Headless CLI usage” section in `README.md`.

## Tests and CI

After configuring with CMake, run:

```bash
cd build
ctest --output-on-failure          # Ninja/Make
ctest -C Debug --output-on-failure # Visual Studio/MSVC
```

The test suite currently covers:

- `config_loader_tests` – JSON/CLI overlay, size guards, and error paths.
- `fdtd_core_tests` – timestep calculation, allocation failure paths, port
  reset behavior.
- `analysis_tests` – safe scope/ports allocation and S‑parameter sampling.
- `config_runtime_tests` – clamping/validation of `SimulationConfig`.

CI runs the same `ctest` invocations plus two smoke tests:

- `emwave_cli_smoke` – fixed‑step headless run on `waveguide.json`.
- `emwave_cli_sweep_smoke` – sweep run on `cpw_filter.json`, validated via the
  S21 text output.

