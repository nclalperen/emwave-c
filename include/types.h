// ============================================================================
// emwave-c: Type Definitions and Structures
// ============================================================================

#ifndef EMWAVE_TYPES_H
#define EMWAVE_TYPES_H

#include <stddef.h>
#include "config.h"  /* For NX, NY */

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

/* Slider UI element */
typedef struct {
    int x, y, w, h;
    double minv, maxv;  /* logical range */
    double value;       /* current logical value */
    int dragging;
} Slider;

/* Autoscale modes */
typedef enum {
    AS_PEAK = 0,
    AS_P99 = 1
} AutoScaleMode;

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

/* Simulation state structure */
typedef struct {
    /* Field arrays */
    double (*Ez)[NY];
    double (*Hx)[NY];
    double (*Hy)[NY];
    double (*Ez_old)[NY];

    /* CPML auxiliary fields */
    double (*psi_Ezx)[NY];
    double (*psi_Ezy)[NY];
    double (*psi_Hyx)[NY];
    double (*psi_Hxy)[NY];

    /* Material properties */
    double epsr[NX][NY];
    double sigma_map[NX][NY];
    unsigned char tag_grid[NX][NY];  /* 0=dielectric, 1=PEC, 2=PMC */

    /* Grid parameters */
    double dx, dy, dt;

    /* Sources */
    Source sources[MAX_SRC];

    /* Ports */
    Port ports[MAX_PORTS];
    int ports_on;

    /* Time step counter */
    int timestep;

    /* Frequency */
    double freq;

} SimulationState;

#endif /* EMWAVE_TYPES_H */
