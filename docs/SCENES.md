# Example scenes

emwave-c ships with a couple of JSON scenes under `configs/` that are useful as starting points or regression references.

- **Waveguide (`configs/waveguide.json`)**
  - Rectangular dielectric-loaded waveguide with PEC top/bottom walls.
  - Single Gaussian source near the left, sweeping 6–8.5 GHz.
  - Good for visualizing mode structure and basic reflections.

- **CPW filter (`configs/cpw_filter.json`)**
  - Coplanar-waveguide filter over a dielectric substrate.
  - Two CW sources placed at each end; only the left source is active by default.
  - Simulation block defines conductor rails and ground shields; sweep settings are 3–6 GHz with 5 points.
  - Suitable for running S-parameter sweeps (S21) via `emwave_cli --run-mode=sweep`.

You can copy these configs and tweak `nx`/`ny`, material rectangles, sources, and sweep parameters to build new scenes. The full JSON schema is documented in `configs/SCHEMA.md`.

