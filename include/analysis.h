// ============================================================================
// emwave-c: Analysis and Measurement Tools
// ============================================================================

#ifndef EMWAVE_ANALYSIS_H
#define EMWAVE_ANALYSIS_H

#include "types.h"
#include <stdio.h>

/* Oscilloscope functions */
int scope_init(Scope* scope, int width);
void scope_free(Scope* scope);
void scope_push(Scope* scope, double v);
void scope_clear(Scope* scope);

/* FFT export */
int dump_scope_fft_csv(const Scope* scope, const char* path, double dt, int Nfft_requested);

/* Port management */
int ports_init(Port* ports, int nx, int ny);
void ports_free(Port* ports);
void ports_sample(SimulationState* state, double dx, double dy);

/* S-parameter computation */
double compute_s21(const Port* ports, double freq, double dt);

/* Probe logging */
typedef struct {
    FILE* file;
    char* buffer;
    size_t capacity;
    size_t size;
} ProbeLogWriter;

int probe_writer_open(ProbeLogWriter* writer, const char* path);
void probe_writer_flush(ProbeLogWriter* writer);
void probe_writer_close(ProbeLogWriter* writer);
void probe_writer_append(ProbeLogWriter* writer, int timestep, double value);
int probe_writer_is_open(const ProbeLogWriter* writer);

FILE* probe_open(const char* filename);
void probe_log(FILE* f, int timestep, double value);

/* Testing hooks */
void analysis_test_set_alloc_fail_after(int count);

#endif /* EMWAVE_ANALYSIS_H */
