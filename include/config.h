// ============================================================================
// emwave-c: Configuration Constants
// ============================================================================

#ifndef EMWAVE_CONFIG_H
#define EMWAVE_CONFIG_H

#include <math.h>

/* Optional instrumentation (ports/S-parameters) */
#ifndef EMWAVE_ENABLE_PORTS
#define EMWAVE_ENABLE_PORTS 1
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simulation parameters */
#define STEPS_PER_FRAME 2

/* Runtime grid defaults */
#define NX_DEFAULT 400
#define NY_DEFAULT 400
#define LX_DEFAULT 0.6
#define LY_DEFAULT 0.6

/* Hard limits to guard against GPU/CPU exhaustion */
#define SIM_MIN_DIM 64
#define SIM_MAX_DIM 4096
#define SIM_MAX_CELLS (4096 * 2048)  /* ~8.3 million cells */

/* Preserve legacy macros for older monolithic builds */
#ifndef NX
#define NX NX_DEFAULT
#endif
#ifndef NY
#define NY NY_DEFAULT
#endif

/* Physical constants */
#define c0   299792458.0              /* m/s */
#define MU0  1.256637061435917295e-6  /* H/m (4π×10^-7) */
#define EPS0 8.8541878128e-12         /* F/m */

/* Physical domain (FIXED mode geometry) */
#define Lx LX_DEFAULT        /* meters */
#define Ly LY_DEFAULT        /* meters */
#define BASE_DX (Lx / NX)
#define BASE_DY (Ly / NY)

/* AUTO mode parameters */
#define TARGET_CPW       12.0   /* desired cells/λ in slowest medium */
#define EPSR_MAX_SCENE   4.0    /* slowest medium (inside block) */
#define MAX_SCALE_FACTOR 20.0   /* Δ_auto <= MAX_SCALE_FACTOR * BASE_DX */

/* Conductivity defaults */
#define SIGMA_BG    0.0
#define SIGMA_BLOCK 0.0

/* Rendering and visualization constants */
#define HIST_STRIDE 1                  /* Histogram sampling stride */
#define COLOR_HIST_BINS 256            /* Number of bins for color histogram */
#define UI_DELAY_MS 8                  /* Frame delay in milliseconds */

/* Numerical safety constants */
#define CFL_SAFETY_FACTOR 0.95         /* CFL timestep safety factor */
#define DIVISION_SAFETY_EPSILON 1e-20  /* Small number to prevent division by zero */

/* CPML parameters */
#ifndef PML_THICK
#define PML_THICK 12
#endif

/* Frequency sweep parameters */
#define SWEEP_MAX_POINTS 32

/* Scene layout limits */
#define CONFIG_MAX_MATERIAL_RECTS 16
#define MAX_SRC 4

typedef enum {
    SRC_CW = 0,
    SRC_GAUSS_PULSE = 1,
    SRC_RICKER = 2
} SourceType;

typedef struct {
    double x0;
    double y0;
    double x1;
    double y1;
    double epsr;
    double sigma;
    unsigned char tag; /* 0=dielectric,1=PEC,2=PMC */
} MaterialRectSpec;

typedef struct {
    int active;
    double x;
    double y;
    SourceType type;
    double amp;
    double freq;
    double sigma2;
} SourceConfigSpec;

/* Frequency range for UI */
#define FREQ_MIN 1e6
#define FREQ_MAX 5e9

typedef enum {
    SIM_RUN_MODE_FIXED_STEPS = 0,
    SIM_RUN_MODE_SWEEP = 1
} SimulationRunMode;

typedef enum {
    SIM_BOUNDARY_CPML = 0,
    SIM_BOUNDARY_MUR = 1
} SimulationBoundaryMode;

#define SIM_PROBE_LOG_PATH_MAX 260

typedef struct {
    int nx;
    int ny;
    double lx;
    double ly;
    double cfl_safety;
    int steps_per_frame;
    int sweep_points;
    double sweep_start_hz;
    double sweep_stop_hz;
    int sweep_steps_per_point;
    SimulationRunMode run_mode;
    int run_steps;
    SimulationBoundaryMode boundary_mode;
    int enable_probe_log;
    char probe_log_path[SIM_PROBE_LOG_PATH_MAX];
    int material_rect_count;
    MaterialRectSpec material_rects[CONFIG_MAX_MATERIAL_RECTS];
    int source_count;
    SourceConfigSpec source_configs[MAX_SRC];
} SimulationConfig;

extern const SimulationConfig SIM_CONFIG_DEFAULTS;

#endif /* EMWAVE_CONFIG_H */
