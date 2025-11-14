// ============================================================================
// emwave-c: Analysis and Measurement Tools
// ============================================================================

#include "analysis.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
}

/* Clear oscilloscope */
void scope_clear(Scope* scope) {
    if (!scope->y) return;
    memset(scope->y, 0, sizeof(double) * scope->n);
    scope->head = 0;
    scope->last = 0.0;
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

/* Initialize ports */
void ports_init(Port* ports) {
    for (int p = 0; p < MAX_PORTS; p++) {
        ports[p].active = 0;
        ports[p].x = (p == 0) ? NX/4 : 3*NX/4;  /* Default positions */
        ports[p].y0 = NY/4;
        ports[p].y1 = 3*NY/4;
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

        double Vsum = 0.0;
        double Isum = 0.0;

        /* Voltage = integral of Ez along port */
        for (int yy = state->ports[p].y0; yy <= state->ports[p].y1; yy++) {
            Vsum += state->Ez[state->ports[p].x][yy];
        }
        Vsum *= dy;

        /* Current = integral of Hy along port */
        for (int yy = state->ports[p].y0; yy <= state->ports[p].y1; yy++) {
            Isum += state->Hy[state->ports[p].x][yy];
        }
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

/* Open probe log file */
FILE* probe_open(const char* filename) {
    return fopen(filename, "w");
}

/* Log probe value */
void probe_log(FILE* f, int timestep, double value) {
    if (f) {
        fprintf(f, "%d %.9e\n", timestep, value);
        if ((timestep & 255) == 0) fflush(f);
    }
}
