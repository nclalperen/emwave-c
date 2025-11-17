// =============================================================================
// emwave-c: Port sampling and S-parameter instrumentation helpers
// =============================================================================

#include "ports.h"
#include "config.h"
#include "util.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Allocation helper shared with analysis tests */
extern void* analysis_checked_calloc(size_t count, size_t size);

static const int kPortMinSegmentLen = 2;

static void port_reset(Port* port) {
    if (!port) {
        return;
    }
    port->x = 0;
    port->y0 = 0;
    port->y1 = 0;
    port->len = 0;
    port->n = 0;
    port->head = 0;
    port->active = 0;
    port->V = NULL;
    port->I = NULL;
}

static void port_deactivate(Port* port) {
    if (!port) {
        return;
    }
    if (port->n > 0) {
        size_t bytes = (size_t)port->n * sizeof(double);
        if (port->V) {
            memset(port->V, 0, bytes);
        }
        if (port->I) {
            memset(port->I, 0, bytes);
        }
    }
    port->head = 0;
    port->len = 0;
    port->active = 0;
}

static int choose_port_segment(int ny, int* out_y0, int* out_y1) {
    if (!out_y0 || !out_y1) {
        return 0;
    }
    if (ny < 2) {
        return 0;
    }

    int safe_lo = (ny > 2) ? 1 : 0;
    int safe_hi = (ny > 2) ? ny - 2 : ny - 1;
    if (safe_hi < safe_lo) {
        safe_hi = safe_lo;
    }

    int y0 = util_clamp_int(ny / 4, safe_lo, safe_hi);
    int y1 = util_clamp_int((3 * ny) / 4, y0 + 1, ny - 1);
    if (y1 < y0) {
        y1 = y0;
    }

    int span = y1 - y0 + 1;
    if (span < kPortMinSegmentLen) {
        int deficit = kPortMinSegmentLen - span;
        int grow_lo = (deficit + 1) / 2;
        int grow_hi = deficit / 2;
        y0 = util_clamp_int(y0 - grow_lo, safe_lo, y0);
        y1 = util_clamp_int(y1 + grow_hi, y0, ny - 1);
        span = y1 - y0 + 1;
    }

    if (span < kPortMinSegmentLen) {
        y0 = safe_lo;
        y1 = util_clamp_int(y0 + kPortMinSegmentLen - 1, y0, ny - 1);
        span = y1 - y0 + 1;
    }

    if (span < kPortMinSegmentLen) {
        return 0;
    }

    *out_y0 = y0;
    *out_y1 = y1;
    return 1;
}

/* Initialize ports */
int ports_init(Port* ports, int nx, int ny) {
    int safe_x_lo = 1;
    int safe_x_hi = (nx > 1) ? nx - 2 : 0;

    if (!ports) {
        return 0;
    }

    ports_free(ports);

    if (nx < 2 || ny < 2) {
        fprintf(stderr, "Warning: grid too small for port initialization (%d x %d)\n", nx, ny);
        return 0;
    }

    for (int p = 0; p < MAX_PORTS; p++) {
        int target_x = (p == 0) ? nx / 4 : (3 * nx) / 4;
        int px = util_clamp_int(target_x, safe_x_lo, safe_x_hi);

        int y0 = 0;
        int y1 = 0;
        if (!choose_port_segment(ny, &y0, &y1)) {
            fprintf(stderr, "Warning: Failed to choose port segment for ny=%d\n", ny);
            ports_free(ports);
            return 0;
        }

        int segment_len = y1 - y0 + 1;
        double* vbuf = (double*)analysis_checked_calloc(PORT_SIGNAL_LENGTH, sizeof(double));
        double* ibuf = (double*)analysis_checked_calloc(PORT_SIGNAL_LENGTH, sizeof(double));
        if (!vbuf || !ibuf) {
            free(vbuf);
            free(ibuf);
            fprintf(stderr, "Warning: Failed to allocate port %d buffers\n", p);
            ports_free(ports);
            return 0;
        }
        port_reset(&ports[p]);
        ports[p].x = px;
        ports[p].y0 = y0;
        ports[p].y1 = y1;
        ports[p].len = segment_len;
        ports[p].n = PORT_SIGNAL_LENGTH;
        ports[p].V = vbuf;
        ports[p].I = ibuf;
        ports[p].head = 0;
        ports[p].active = 1;
    }
    return 1;
}

/* Free ports */
void ports_free(Port* ports) {
    if (!ports) {
        return;
    }
    for (int p = 0; p < MAX_PORTS; p++) {
        if (ports[p].V) {
            free(ports[p].V);
        }
        if (ports[p].I) {
            free(ports[p].I);
        }
        port_reset(&ports[p]);
    }
}

static int port_prepare_sample_range(const SimulationState* state, Port* port,
                                     int* out_y0, int* out_y1) {
    if (!state || !port || !out_y0 || !out_y1) {
        return 0;
    }
    if (!port->V || !port->I || port->n <= 0) {
        return 0;
    }
    if (!state->Ez || !state->Hy) {
        return 0;
    }
    if (state->nx < 2 || state->ny < kPortMinSegmentLen) {
        return 0;
    }

    int px = port->x;
    int min_px = 1;
    int max_px = (state->nx > 2) ? state->nx - 2 : 0;
    if (max_px < min_px) {
        return 0;
    }

    int clamped_px = util_clamp_int(px, min_px, max_px);
    if (clamped_px != px) {
        port->x = clamped_px;
        px = clamped_px;
    }

    int y0 = util_clamp_int(port->y0, 0, state->ny - 1);
    int y1 = util_clamp_int(port->y1, y0, state->ny - 1);
    if (y1 < y0) {
        int tmp = y0;
        y0 = y1;
        y1 = tmp;
    }

    int span = y1 - y0 + 1;
    if (span < kPortMinSegmentLen) {
        int deficit = kPortMinSegmentLen - span;
        int grow_top = (deficit + 1) / 2;
        int grow_bottom = deficit - grow_top;
        y0 = util_clamp_int(y0 - grow_bottom, 0, state->ny - 1);
        y1 = util_clamp_int(y1 + grow_top, y0, state->ny - 1);
        span = y1 - y0 + 1;
    }

    if (span < kPortMinSegmentLen) {
        return 0;
    }

    if (port->len != span) {
        port->len = span;
    }

    *out_y0 = y0;
    *out_y1 = y1;
    return 1;
}

/* Sample port voltages and currents */
void ports_sample(SimulationState* state, double dx, double dy) {
    if (!state) {
        return;
    }

    for (int p = 0; p < MAX_PORTS; p++) {
        Port* port = &state->ports[p];
        if (!port->active) continue;

        int y0 = 0;
        int y1 = 0;
        if (!port_prepare_sample_range(state, port, &y0, &y1)) {
            if (port->V && port->I) {
                port_deactivate(port);
            }
            continue;
        }

        int px = port->x;

        double Vsum = 0.0;
        double Isum = 0.0;

        for (int yy = y0; yy <= y1; yy++) {
            Vsum += state->Ez[px][yy];
            Isum += state->Hy[px][yy];
        }

        Vsum *= dy;
        Isum *= dx;

        port->V[port->head] = Vsum;
        port->I[port->head] = Isum;
        port->head = (port->head + 1) % port->n;
    }
}

