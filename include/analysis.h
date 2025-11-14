// ============================================================================
// emwave-c: Analysis and Measurement Tools
// ============================================================================

#ifndef EMWAVE_ANALYSIS_H
#define EMWAVE_ANALYSIS_H

#include "types.h"
#include <stdio.h>

/* Oscilloscope functions */
void scope_init(Scope* scope, int width);
void scope_free(Scope* scope);
void scope_push(Scope* scope, double v);
void scope_clear(Scope* scope);

/* FFT export */
int dump_scope_fft_csv(const Scope* scope, const char* path, double dt, int Nfft_requested);

/* Port management */
void ports_init(Port* ports, int nx, int ny);
void ports_free(Port* ports);
void ports_sample(SimulationState* state, double dx, double dy);

/* S-parameter computation */
double compute_s21(const Port* ports, double freq, double dt);

/* Probe logging */
FILE* probe_open(const char* filename);
void probe_log(FILE* f, int timestep, double value);

#endif /* EMWAVE_ANALYSIS_H */
