// ============================================================================
// emwave-c: Type Definitions and Structures
// ============================================================================

#ifndef EMWAVE_TYPES_H
#define EMWAVE_TYPES_H

#include <stddef.h>
#include "config.h"

/* Source types */
typedef enum {
    SRC_CW = 0,
    SRC_GAUSS_PULSE = 1,
    SRC_RICKER = 2
} SourceType;

/* Source structure */
#define MAX_SRC 4
typedef struct {
    int active;
    int ix, iy;
    SourceType type;
    double amp, freq;
    double t0, tau;     /* pulse parameters */
    double sigma2;      /* spatial footprint variance in cells^2 */
} Source;

/* Port structure for S-parameter measurements */
#define MAX_PORTS 2
#define PORT_SIGNAL_LENGTH 4096
typedef struct {
    int x;              /* x index (column) where port is sampled */
    int y0, y1;         /* inclusive y-range of the port segment */
    int len;            /* number of cells along port */
    int n;              /* length of the sample buffer */
    double *V;          /* accumulated voltage samples (integral of Ez) */
    double *I;          /* accumulated current samples (integral of H) */
    int head;           /* write index into V/I buffers (circular) */
    int active;         /* whether this port is active */
} Port;

/* Oscilloscope structure */
typedef struct {
    int n, head;
    double *y;
    int on;
    double last;
} Scope;

/* Boundary condition types */
typedef enum {
    BOUNDARY_MUR = 0,
    BOUNDARY_CPML = 1
} BoundaryType;

/* CPML preset structure */
typedef struct {
    const char* name;
    double smax, kmax, amax;
    int thick;
} CpmlPreset;

/* CPML runtime state */
typedef struct {
    int enabled;
    int thickness;
    int preset_idx;
    double sigma_max;
    double kappa_max;
    double alpha_max;
    BoundaryType boundary_type;

    double* kx;
    double* bx;
    double* cx;
    double* ky;
    double* by;
    double* cy;

    int kx_capacity;
    int ky_capacity;
} CpmlState;

/* Simulation state structure */
typedef struct {
    SimulationConfig config;

    /* Runtime grid */
    int nx, ny;
    double lx, ly;
    double cfl_safety;
    double dx, dy, dt;

    /* Sweep configuration */
    int sweep_points;
    double sweep_start_hz;
    double sweep_stop_hz;
    int sweep_steps_per_point;

    /* Field arrays */
    double **Ez;
    double **Hx;
    double **Hy;
    double **Ez_old;

    double *Ez_data;
    double *Hx_data;
    double *Hy_data;
    double *Ez_old_data;

    /* CPML auxiliary fields */
    double **psi_Ezx;
    double **psi_Ezy;
    double **psi_Hyx;
    double **psi_Hxy;

    double *psi_Ezx_data;
    double *psi_Ezy_data;
    double *psi_Hyx_data;
    double *psi_Hxy_data;

    /* Material properties */
    double **epsr;
    double **sigma_map;
    unsigned char **tag_grid;  /* 0=dielectric, 1=PEC, 2=PMC */

    double *epsr_data;
    double *sigma_map_data;
    unsigned char *tag_grid_data;

    /* Sources */
    Source sources[MAX_SRC];

    /* Ports */
    Port ports[MAX_PORTS];
    int ports_on;

    /* Time step counter */
    int timestep;

    /* Frequency */
    double freq;

    /* Boundary / CPML state */
    CpmlState cpml;

} SimulationState;

#endif /* EMWAVE_TYPES_H */
