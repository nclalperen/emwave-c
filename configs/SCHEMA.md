# emwave JSON configuration schema

Every config file consumed by `--config=` is a JSON document with three optional
sections:

```json
{
  "simulation": { ... },
  "materials": [ ... ],
  "sources": [ ... ]
}
```

## `simulation`

All fields are optional; omitted values fall back to the built-in defaults.

| Key | Type | Description |
| --- | --- | --- |
| `nx`, `ny` | integer | Grid cells along X/Y. Clamped to `[SIM_MIN_DIM, SIM_MAX_DIM]`. |
| `lx`, `ly` | number | Physical domain length in meters. |
| `cfl` | number | CFL safety factor between `0` and `1`. |
| `steps_per_frame` | integer | Update cadence for the UI. |
| `sweep_points` | integer | Number of sweep frequencies to evaluate. |
| `sweep_start_hz`, `sweep_stop_hz` | number | Frequency range in Hz. |
| `sweep_steps_per_point` | integer | Steps to simulate per sweep frequency. |
| `run_mode` | string | Optional. `"fixed"` (default) or `"sweep"`. |
| `run_steps` | integer | When `run_mode` is `"fixed"`, total steps to run. |
| `enable_probe_log` | boolean | Enable writing probe samples to a CSV-like text file. |
| `probe_log_path` | string | Path for probe logging when enabled. Defaults to `"probe.txt"`. |
| `boundary` | string | Optional. `"cpml"` (default) or `"mur"`. |

## `materials`

Array of rectangular regions painted into the Yee grid. Coordinates are unitless
fractions spanning the simulation domain (0.0 → left/bottom boundary, 1.0 → right/top boundary).

Each entry supports the following fields:

| Key | Type | Description |
| --- | --- | --- |
| `type` | string | Optional. `"dielectric"` (default), `"pec"`, or `"pmc"`. |
| `x0`, `y0`, `x1`, `y1` | number | Fractions defining the lower-left and upper-right corners. `x1`/`y1` must be greater than `x0`/`y0`. |
| `epsr` | number | Relative permittivity (dielectric only). Defaults to 1.0. |
| `sigma` | number | Conductivity in S/m (dielectric only). Defaults to 0.0. |

Dielectric rectangles paint `epsr`/`sigma` values. PEC and PMC rectangles only
adjust the boundary tag grid.

## `sources`

Array describing the active sources at start-up. Positions are normalized
fractions (0.0 → left/bottom, 1.0 → right/top) and are clamped away from the
PML to avoid numerical issues.

| Key | Type | Description |
| --- | --- | --- |
| `type` | string | `"cw"`, `"gaussian"`, or `"ricker"`. |
| `x`, `y` | number | Normalized location of the source. |
| `amp` | number | Source amplitude. Defaults to 1.0. |
| `freq` | number | Frequency in Hz. Defaults to 1e9. |
| `sigma2` | number | Spatial variance (controls Gaussian footprint). Defaults to 4.0. |
| `active` | boolean | Whether the source is enabled. Defaults to `true`. |

At most `CONFIG_MAX_MATERIAL_RECTS` material rectangles and `MAX_SRC` sources are
loaded. Excess entries are ignored.
