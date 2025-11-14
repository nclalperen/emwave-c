// ============================================================================
// emwave-c: Configuration Constants
// ============================================================================

#ifndef EMWAVE_CONFIG_H
#define EMWAVE_CONFIG_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simulation parameters */
#define STEPS_PER_FRAME 2

/* Grid dimensions (must come before BASE_DX/BASE_DY) */
#define NX 400
#define NY 400

/* Physical constants */
static const double c0   = 299792458.0;              /* m/s */
static const double MU0  = 1.256637061435917295e-6;  /* H/m (4π×1e-7) */
static const double EPS0 = 8.8541878128e-12;         /* F/m */

/* Physical domain (FIXED mode geometry) */
#define Lx 0.6        /* meters */
#define Ly 0.6        /* meters */
#define BASE_DX (Lx / NX)
#define BASE_DY (Ly / NY)

/* AUTO mode parameters */
static const double TARGET_CPW       = 12.0;   /* desired cells/λ in slowest medium */
static const double EPSR_MAX_SCENE   = 4.0;    /* slowest medium (inside block) */
static const double MAX_SCALE_FACTOR = 20.0;   /* Δ_auto ≤ MAX_SCALE_FACTOR * BASE_DX */

/* Conductivity defaults */
static const double SIGMA_BG    = 0.0;
static const double SIGMA_BLOCK = 0.0;

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

/* Frequency range for UI */
#define FREQ_MIN 1e6
#define FREQ_MAX 5e9

#endif /* EMWAVE_CONFIG_H */
