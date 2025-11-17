// ============================================================================
// emwave-c: Wave Source Management
// ============================================================================

#ifndef EMWAVE_SOURCES_H
#define EMWAVE_SOURCES_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Source initialization and management */
void sources_init(Source* sources, int nx, int ny, const SimulationConfig* cfg);
void sources_shutdown(Source* sources);
void source_reparam(Source* s);
void sources_set_freq(Source* sources, double f);
void sources_cycle_type(Source* sources);

/* Source injection */
double source_time_value(const Source* s, int t, double dt);
void inject_source_into_Ez(Source* s, SimulationState* state, double dt);
void inject_all_sources(SimulationState* state);

/* Source interaction (for UI) */
int find_nearest_source(const Source* sources, int mx, int my, int scale, float max_distance);

#ifdef __cplusplus
}
#endif

#endif /* EMWAVE_SOURCES_H */
