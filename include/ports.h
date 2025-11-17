// =============================================================================
// emwave-c: Port sampling and S-parameter instrumentation
// =============================================================================

#ifndef EMWAVE_PORTS_H
#define EMWAVE_PORTS_H

#include "types.h"

/* Port buffer lifecycle */
int ports_init(Port* ports, int nx, int ny);
void ports_free(Port* ports);

/* Sample port voltages and currents into circular buffers */
void ports_sample(SimulationState* state, double dx, double dy);

#endif /* EMWAVE_PORTS_H */

