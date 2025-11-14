// ============================================================================
// emwave-c: Analysis and Measurement Tools
// ============================================================================

#include "analysis.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Initialize oscilloscope */
void scope_init(Scope* scope, int width) {
    if (scope->y) free(scope->y);

    scope->n = (width > 64 ? width : 64);
    scope->y = (double*)calloc(scope->n, sizeof(double));

    if (!scope->y) {
        fprintf(stderr, "Warning: Failed to allocate scope buffer\n");
        scope->n = 0;
        scope->on = 0;
        return;
    }

    scope->head = 0;
    scope->on = 1;
    scope->last = 0.0;
    scope->rolling_absmax = 0.0;
    scope->rolling_generation = 0;
}

/* Free oscilloscope */
void scope_free(Scope* scope) {
    if (scope->y) {
        free(scope->y);
        scope->y = NULL;
    }
    scope->n = 0;
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

    double *x = (double*)malloc(sizeof(double) * Nfft);
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

static inline int clampi_local(int v, int lo, int hi) {
    if (lo > hi) {
        int tmp = lo;
        lo = hi;
        hi = tmp;
    }
    return v < lo ? lo : (v > hi ? hi : v);
}

/* Initialize ports */
void ports_init(Port* ports, int nx, int ny) {
    int safe_x_lo = 1;
    int safe_x_hi = (nx > 1) ? nx - 2 : 0;
    int safe_y_lo = 1;
    int safe_y_hi = (ny > 1) ? ny - 2 : 0;

    for (int p = 0; p < MAX_PORTS; p++) {
        ports[p].active = 0;
        int target_x = (p == 0) ? nx / 4 : (3 * nx) / 4;
        ports[p].x = clampi_local(target_x, safe_x_lo, safe_x_hi);

        int y0 = clampi_local(ny / 4, safe_y_lo, safe_y_hi);
        int y1 = clampi_local((3 * ny) / 4, y0 + 1, safe_y_hi + 1);
        if (y1 < y0) y1 = y0;

        ports[p].y0 = y0;
        ports[p].y1 = y1;
        ports[p].len = ports[p].y1 - ports[p].y0 + 1;
        ports[p].n = PORT_SIGNAL_LENGTH;
        ports[p].V = (double*)calloc(PORT_SIGNAL_LENGTH, sizeof(double));
        ports[p].I = (double*)calloc(PORT_SIGNAL_LENGTH, sizeof(double));
        ports[p].head = 0;

        if (!ports[p].V || !ports[p].I) {
            fprintf(stderr, "Warning: Failed to allocate port %d buffers\n", p);
            ports[p].active = 0;
        }
    }
}

/* Free ports */
void ports_free(Port* ports) {
    for (int p = 0; p < MAX_PORTS; p++) {
        if (ports[p].V) {
            free(ports[p].V);
            ports[p].V = NULL;
        }
        if (ports[p].I) {
            free(ports[p].I);
            ports[p].I = NULL;
        }
    }
}

/* Sample port voltages and currents */
void ports_sample(SimulationState* state, double dx, double dy) {
    for (int p = 0; p < MAX_PORTS; p++) {
        if (!state->ports[p].active) continue;

        int px = state->ports[p].x;
        if (px <= 0 || px >= state->nx) continue;

        int y0 = clampi_local(state->ports[p].y0, 0, state->ny - 1);
        int y1 = clampi_local(state->ports[p].y1, y0, state->ny - 1);

        double Vsum = 0.0;
        double Isum = 0.0;

        for (int yy = y0; yy <= y1; yy++) {
            Vsum += state->Ez[px][yy];
            Isum += state->Hy[px][yy];
        }

        Vsum *= dy;
        Isum *= dx;

        state->ports[p].V[state->ports[p].head] = Vsum;
        state->ports[p].I[state->ports[p].head] = Isum;
        state->ports[p].head = (state->ports[p].head + 1) % state->ports[p].n;
    }
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
