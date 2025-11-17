// ============================================================================
// emwave-c: Analysis and Measurement Tools
// ============================================================================

#include "analysis.h"
#include "config.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

static int analysis_alloc_fail_after = -1;
static const int kScopeMinSamples = 64;

static void scope_reset(Scope* scope) {
    if (!scope) {
        return;
    }
    scope->y = NULL;
    scope->n = 0;
    scope->head = 0;
    scope->on = 0;
    scope->last = 0.0;
    scope->rolling_absmax = 0.0;
    scope->rolling_generation = 0;
}

void analysis_test_set_alloc_fail_after(int count) {
    analysis_alloc_fail_after = count;
}

static int analysis_should_fail_allocation(void) {
    if (analysis_alloc_fail_after < 0) {
        return 0;
    }
    if (analysis_alloc_fail_after == 0) {
        return 1;
    }
    analysis_alloc_fail_after--;
    return 0;
}

static void* analysis_checked_malloc(size_t size) {
    if (analysis_should_fail_allocation()) {
        return NULL;
    }
    return malloc(size);
}

void* analysis_checked_calloc(size_t count, size_t size) {
    if (analysis_should_fail_allocation()) {
        return NULL;
    }
    return calloc(count, size);
}

/* Initialize oscilloscope */
int scope_init(Scope* scope, int width) {
    if (!scope) {
        return 0;
    }

    scope_free(scope);
    int buffer_len = width;
    if (buffer_len < kScopeMinSamples) {
        buffer_len = kScopeMinSamples;
    }
    double* buf = (double*)analysis_checked_calloc((size_t)buffer_len, sizeof(double));
    if (!buf) {
        fprintf(stderr, "Warning: Failed to allocate scope buffer\n");
        scope_reset(scope);
        return 0;
    }

    scope->y = buf;
    scope->n = buffer_len;
    scope->head = 0;
    scope->on = 1;
    scope->last = 0.0;
    scope->rolling_absmax = 0.0;
    scope->rolling_generation = 0;
    return 1;
}

/* Free oscilloscope */
void scope_free(Scope* scope) {
    if (!scope) {
        return;
    }
    free(scope->y);
    scope_reset(scope);
}

/* Push value to oscilloscope */
void scope_push(Scope* scope, double v) {
    if (!scope->y || !scope->on) return;

    scope->y[scope->head] = v;
    scope->head = (scope->head + 1) % scope->n;
    scope->last = v;

    double abs_v = fabs(v);
    double decayed = scope->rolling_absmax * 0.995;
    scope->rolling_absmax = (abs_v > decayed) ? abs_v : decayed;
    scope->rolling_generation++;
}

/* Clear oscilloscope */
void scope_clear(Scope* scope) {
    if (!scope->y) return;
    memset(scope->y, 0, sizeof(double) * scope->n);
    scope->head = 0;
    scope->last = 0.0;
    scope->rolling_absmax = 0.0;
    scope->rolling_generation++;
}

/* Export FFT of scope data to CSV */
int dump_scope_fft_csv(const Scope* scope, const char* path, double dt, int Nfft_requested) {
    if (!scope->y || scope->n <= 8 || dt <= 0.0) return 0;

    int N = scope->n;
    int Nfft = Nfft_requested;
    if (Nfft > N) Nfft = N;
    if (Nfft < 64) Nfft = (N < 64 ? N : 64);

    double *x = (double*)analysis_checked_malloc(sizeof(double) * (size_t)Nfft);
    if (!x) {
        fprintf(stderr, "Warning: FFT export failed - memory allocation error\n");
        return 0;
    }

    /* Copy data from circular buffer */
    int start = (scope->head - Nfft + scope->n) % scope->n;
    for (int k = 0; k < Nfft; k++) {
        int idx = (start + k) % scope->n;
        x[k] = scope->y[idx];
    }

    /* Remove DC and apply Hann window */
    double mean = 0.0;
    for (int k = 0; k < Nfft; k++) mean += x[k];
    mean /= (double)Nfft;

    for (int k = 0; k < Nfft; k++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)k / (double)(Nfft - 1)));
        x[k] = (x[k] - mean) * w;
    }

    /* Compute DFT and write to file */
    FILE* f = fopen(path, "w");
    if (!f) {
        free(x);
        return 0;
    }

    fprintf(f, "# FFT of scope (Hann); dt=%.9e; N=%d\n", dt, Nfft);
    fprintf(f, "freq_Hz,mag\n");

    for (int k = 0; k <= Nfft/2; k++) {
        double re = 0.0, im = 0.0;
        double ang_step = -2.0 * M_PI * (double)k / (double)Nfft;

        for (int n = 0; n < Nfft; n++) {
            double ang = ang_step * (double)n;
            re += x[n] * cos(ang);
            im += x[n] * sin(ang);
        }

        double mag = sqrt(re*re + im*im);
        double freq = (double)k / (dt * (double)Nfft);
        fprintf(f, "%.12e,%.12e\n", freq, mag);
    }

    fclose(f);
    free(x);
    return 1;
}

/* Compute S21 parameter */
double compute_s21(const Port* ports, double freq, double dt) {
    int n = ports[0].n;
    if (n <= 0 || !ports[0].V || !ports[1].V) return 0.0;

    double re1 = 0.0, im1 = 0.0;
    double re2 = 0.0, im2 = 0.0;

    int h1 = ports[0].head;
    int h2 = ports[1].head;

    for (int i = 0; i < n; i++) {
        double x1 = ports[0].V[(h1 + i) % n];
        double x2 = ports[1].V[(h2 + i) % n];
        double arg = 2.0 * M_PI * freq * dt * (double)i;
        double ca = cos(arg);
        double sa = sin(arg);

        re1 += x1 * ca;
        im1 -= x1 * sa;
        re2 += x2 * ca;
        im2 -= x2 * sa;
    }

    double mag1 = sqrt(re1 * re1 + im1 * im1);
    double mag2 = sqrt(re2 * re2 + im2 * im2);

    if (mag1 < 1e-30) return 0.0;
    return mag2 / mag1;
}

/* Probe logging helpers */
#define PROBE_WRITER_BLOCK_SIZE 4096

int probe_writer_open(ProbeLogWriter* writer, const char* path) {
    if (!writer || !path || !*path) return 0;
    memset(writer, 0, sizeof(*writer));
    writer->file = fopen(path, "w");
    if (!writer->file) {
        return 0;
    }
    writer->buffer = (char*)malloc(PROBE_WRITER_BLOCK_SIZE);
    if (!writer->buffer) {
        fclose(writer->file);
        writer->file = NULL;
        return 0;
    }
    writer->capacity = PROBE_WRITER_BLOCK_SIZE;
    writer->size = 0;
    return 1;
}

int probe_writer_is_open(const ProbeLogWriter* writer) {
    return writer && writer->file;
}

void probe_writer_flush(ProbeLogWriter* writer) {
    if (!writer || !writer->file || writer->size == 0) {
        return;
    }
    fwrite(writer->buffer, 1, writer->size, writer->file);
    writer->size = 0;
}

void probe_writer_close(ProbeLogWriter* writer) {
    if (!writer) return;
    probe_writer_flush(writer);
    if (writer->file) {
        fclose(writer->file);
        writer->file = NULL;
    }
    if (writer->buffer) {
        free(writer->buffer);
        writer->buffer = NULL;
    }
    writer->capacity = 0;
    writer->size = 0;
}

void probe_writer_append(ProbeLogWriter* writer, int timestep, double value) {
    if (!writer || !writer->file || !writer->buffer) return;

    char line[64];
    int len = snprintf(line, sizeof(line), "%d %.9e\n", timestep, value);
    if (len <= 0) {
        return;
    }
    size_t needed = (size_t)len;
    if (needed >= sizeof(line)) {
        needed = sizeof(line);
    }
    if (needed > writer->capacity) {
        fwrite(line, 1, needed, writer->file);
        return;
    }
    if (writer->size + needed > writer->capacity) {
        probe_writer_flush(writer);
    }
    memcpy(writer->buffer + writer->size, line, needed);
    writer->size += needed;
}

/* Legacy helpers retained for monolithic builds */
FILE* probe_open(const char* filename) {
    return fopen(filename, "w");
}

void probe_log(FILE* f, int timestep, double value) {
    if (f) {
        fprintf(f, "%d %.9e\n", timestep, value);
        if ((timestep & 255) == 0) fflush(f);
    }
}
