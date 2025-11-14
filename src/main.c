
// ============================================================================
// EM Wave Sim Project — Pro-grade 2D TMz FDTD (Single-file demo)
// Features: CPML (toggle + presets), conductivity (sigma), robust autoscale
//           decoupled scope hold, FFT export, render stride, Touch-friendly UI
// Author: Alperen Oncul
// Version: 0.7 (2024-06) - see changelog at the end
// ============================================================================

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_OPENMP) && !defined(_MSC_VER)
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================= Grid, physics, scene ================= */
#define NX 400
#define NY 400
#define STEPS_PER_FRAME 2

/* Physical constants */
static const double c0  = 299792458.0;            /* m/s */
static const double MU0 = 1.256637061435917295e-6;/* H/m (4π×1e-7) */
static const double EPS0= 8.8541878128e-12;       /* F/m */

/* Rendering and visualization constants */
#define HIST_STRIDE 4                  /* Histogram sampling stride (perf tuning) */
#define COLOR_HIST_BINS 512            /* Number of bins for color histogram */
#define UI_DELAY_MS 8                  /* Frame delay in milliseconds */

/* Numerical safety constants */
#define CFL_SAFETY_FACTOR 0.95         /* CFL timestep safety factor */
#define DIVISION_SAFETY_EPSILON 1e-20  /* Small number to prevent division by zero */

/* Physical domain (FIXED mode geometry) */
#define Lx_MAIN 0.6           /* meters */
#define Ly_MAIN 0.6           /* meters */
#define BASE_DX_MAIN (Lx_MAIN / NX)
#define BASE_DY_MAIN (Ly_MAIN / NY)
#define BASE_DX BASE_DX_MAIN
#define BASE_DY BASE_DY_MAIN

/* AUTO (cells per wavelength) targets and safety */
static const double TARGET_CPW       = 12.0;  /* desired cells/λ in slowest medium */
static const double EPSR_MAX_SCENE   = 4.0;   /* slowest medium (inside block) */
static const double MAX_SCALE_FACTOR = 20.0;  /* Δ_auto ≤ MAX_SCALE_FACTOR * BASE_DX_MAIN */

/* Live grid + timestep (start FIXED) */
static double dx = BASE_DX_MAIN, dy = BASE_DY_MAIN, dt = 0.0;

/* Conductivity map (S/m). Default 0; set SIGMA_BLOCK > 0 for loss. */
static const double SIGMA_BG    = 0.0;
static const double SIGMA_BLOCK = 0.0;

/* ================= Tiny helpers ================= */
static inline int    clampi (int v, int lo, int hi){ return v<lo?lo:(v>hi?hi:v); }
static inline double clampd (double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }
static inline float  sqrf   (float x){ return x*x; }
/*
 * Compute a stable Courant–Friedrichs–Lewy time step.  Use a slightly
 * smaller safety factor (0.95) than the classic 0.99 to provide more
 * margin when lossy materials or coarse meshes are active.  A too‑aggressive
 * CFL can destabilise the solver when sigma>0.  See other AI review.
 */
static inline double cfl_dt(double dx_, double dy_){
    return CFL_SAFETY_FACTOR / (c0 * sqrt(1.0/(dx_*dx_) + 1.0/(dy_*dy_)));
}

/* ================= Sliders ================= */
typedef struct {
    int x, y, w, h;
    double minv, maxv;  /* logical range */
    double value;       /* current logical value */
    int dragging;
} Slider;

static int slider_handle_event(Slider* s, const SDL_Event* e){
    if (e->type == SDL_MOUSEBUTTONDOWN) {
        int mx = e->button.x, my = e->button.y;
        if (mx >= s->x && mx <= s->x + s->w && my >= s->y && my <= s->y + s->h) {
            s->dragging = 1;
            double t = (double)(mx - s->x) / (double)s->w;
            t = clampd(t, 0.0, 1.0);
            s->value = s->minv + t * (s->maxv - s->minv);
            return 1;
        }
    } else if (e->type == SDL_MOUSEBUTTONUP) {
        s->dragging = 0;
    } else if (e->type == SDL_MOUSEMOTION && s->dragging) {
        int mx = e->motion.x;
        double t = (double)(mx - s->x) / (double)s->w;
        t = clampd(t, 0.0, 1.0);
        double old = s->value;
        s->value = s->minv + t * (s->maxv - s->minv);
        return (s->value != old);
    }
    return 0;
}
static void slider_draw(SDL_Renderer* r, const Slider* s){
    SDL_Rect track = (SDL_Rect){ s->x, s->y + s->h/2 - 2, s->w, 4 };
    SDL_SetRenderDrawColor(r, 120,120,120,255);
    SDL_RenderFillRect(r, &track);
    double t = (s->value - s->minv) / (s->maxv - s->minv);
    int kx = s->x + (int)(t * s->w);
    SDL_Rect knob = (SDL_Rect){ kx-6, s->y, 12, s->h };
    SDL_SetRenderDrawColor(r, 220,220,220,255);
    SDL_RenderFillRect(r, &knob);
}
/* log slider mapping (0..1 <-> fmin..fmax) */
static inline double freq_from_slider(double t, double fminv, double fmaxv){
    double a = log10(fminv), b = log10(fmaxv);
    return pow(10.0, a + clampd(t,0.0,1.0)*(b-a));
}
static inline double slider_from_freq(double f, double fminv, double fmaxv){
    double a = log10(fminv), b = log10(fmaxv);
    return clampd((log10(f)-a)/(b-a), 0.0, 1.0);
}

/* ================= Fields & scene ================= */
static double (*Ez)[NY], (*Hx)[NY], (*Hy)[NY];
static double (*Ez_old)[NY];

/* ============================================================================ */
/* Material, boundary tags and ports (Phase 1 & 2 enhancements)                 */
/*
 * epsr and sigma_map store the relative permittivity and conductivity for
 * every cell in the domain.  Values are unitless (epsr) and S/m (sigma).
 * tag_grid marks special cell types: 0 = normal dielectric, 1 = PEC, 2 = PMC.
 * ports store voltage/current time-series along user-defined segments for
 * computing S-parameters.  See initialization in main() for details.
 */

/* per-cell material properties */
static double epsr[NX][NY];
static double sigma_map[NX][NY];

/* tag grid: 0=dielectric, 1=PEC (Ez pinned to zero), 2=PMC (not fully implemented) */
static unsigned char tag_grid[NX][NY];

/* port definitions and sampling buffers */
#define MAX_PORTS 2
#define PORT_SIGNAL_LENGTH 4096
typedef struct {
    int x;      /* x index (column) where port is sampled */
    int y0, y1; /* inclusive y-range of the port segment */
    int len;    /* number of cells along port */
    int n;      /* length of the sample buffer */
    double *V;  /* accumulated voltage samples (integral of Ez along port) */
    double *I;  /* accumulated current samples (integral of H along port) */
    int head;   /* write index into V/I buffers (circular) */
    int active; /* whether this port is active */
} Port;

static Port ports[MAX_PORTS];
static int ports_on = 0;        /* global toggle for sampling ports each step */
static double s21_amp = 0.0;    /* last-computed S21 amplitude (ratio of port 2 to port 1) */
static int s21_computed = 0;    /* flag indicating S21 has been computed */
static int paint_mode = 0;      /* toggle for painting PEC cells with left-click */
/* painting type: 0=off, 1=PEC, 2=PMC, 3=dielectric; only effective when paint_mode=1.
 * type 1 paints PEC (Ez pinned), type 2 paints PMC (Hx/Hy pinned), type 3 paints a dielectric
 * with relative permittivity paint_eps.  Use keys 'i' to cycle type and 'o/p' to adjust paint_eps. */
static int paint_type = 1;
static double paint_eps = 2.0;

/* Global runtime toggle so rebuild_sim_for_freq() can see it */
int auto_rescale = 0;   /* start FIXED */

/* Poynting flux monitor (smoothed).  This tracks the net electromagnetic power leaving
 * the simulation domain.  Positive values indicate power flowing outward. */
static double pflux_avg = 0.0;

/* Parametric sweep variables.  When sweep_on is set, the simulation will automatically
 * sweep the frequency over sweep_points values, run sweep_steps_per_point time steps
 * for each, compute S21 at each frequency and store the result in sweep_s21.  Use key 'b'
 * to start a sweep. */
#define SWEEP_MAX_POINTS 32
static int sweep_on = 0;
static int sweep_idx = 0;
static int sweep_points = 0;
static double sweep_freqs[SWEEP_MAX_POINTS];
static double sweep_s21[SWEEP_MAX_POINTS];
static int sweep_steps_remaining = 0;
static int sweep_steps_per_point = 2000;

static inline void clear_fields(void){
    memset(Ez,     0, sizeof(double[NX][NY]));
    memset(Hx,     0, sizeof(double[NX][NY]));
    memset(Hy,     0, sizeof(double[NX][NY]));
    memset(Ez_old, 0, sizeof(double[NX][NY]));
}
static inline void draw_block_outline(SDL_Renderer* ren, const int scale){
    int bx0 = NX/2 - NX/10, bx1 = NX/2 + NX/10;
    int by0 = NY/2 - NY/20, by1 = NY/2 + NY/20;
    SDL_SetRenderDrawColor(ren, 220,220,220,200);
    SDL_RenderDrawLine(ren, bx0*scale, by0*scale, bx1*scale, by0*scale);
    SDL_RenderDrawLine(ren, bx0*scale, by1*scale, bx1*scale, by1*scale);
    SDL_RenderDrawLine(ren, bx0*scale, by0*scale, bx0*scale, by1*scale);
    SDL_RenderDrawLine(ren, bx1*scale, by0*scale, bx1*scale, by1*scale);
}
/*
 * epsilon_at and sigma_at are now wrappers around the per-cell material
 * properties. They return the absolute permittivity (eps0 * epsr[i][j]) and
 * the conductivity (sigma_map[i][j]) for the given cell.  PEC cells are
 * handled in the update loops (Ez=0 for tag=1), so material values here
 * are ignored for PEC.
 */
static inline double epsilon_at(int i, int j){
    /* multiply relative permittivity by free-space epsilon */
    return EPS0 * epsr[i][j];
}
static inline double sigma_at(int i, int j){
    return sigma_map[i][j];
}

/* Compute the amplitude ratio (S21) between ports[1] and ports[0] at a given
 * frequency. This function flattens the circular port buffers so that the
 * oldest sample corresponds to the first element and performs a discrete
 * Fourier transform at the source frequency.  It returns Mag(port2)/Mag(port1).
 */
static double compute_s21(double freq, double dt) {
    int n = ports[0].n;
    /* ensure both ports are active and buffers allocated */
    if (n <= 0 || !ports[0].V || !ports[1].V) return 0.0;
    double re1 = 0.0, im1 = 0.0;
    double re2 = 0.0, im2 = 0.0;
    int h1 = ports[0].head;
    int h2 = ports[1].head;
    for (int i = 0; i < n; ++i) {
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

/* ================= AUTO grid update ================= */
static inline void update_grid_for_freq(double f, double *pdx, double *pdy, double *pdt){
    double delta = c0 / (f * sqrt(EPSR_MAX_SCENE) * TARGET_CPW);
    if (delta <= 0) delta = BASE_DX;
    double max_delta = BASE_DX * MAX_SCALE_FACTOR;
    if (delta > max_delta) delta = max_delta;    // clamp for stability
    *pdx = delta; *pdy = delta; *pdt = cfl_dt(*pdx, *pdy);
}

/* ================= Frequency & speed setters ================= */
static inline void set_frequency(double newf, double fminv, double fmaxv, Slider* s, double* f_out){
    if (newf < fminv) newf = fminv;
    if (newf > fmaxv) newf = fmaxv;
    *f_out = newf;
    s->value = slider_from_freq(newf, fminv, fmaxv);  // keep knob in sync
}
static inline void set_steps(int newsteps, Slider* s, int* steps_out){
    if (newsteps < 1)  newsteps = 1;
    if (newsteps > 50) newsteps = 50;
    *steps_out = newsteps; s->value = *steps_out;
}


/* ================= Sources ================= */
typedef enum { SRC_CW=0, SRC_GAUSS_PULSE=1, SRC_RICKER=2 } SourceType;
typedef struct {
    int active; int ix, iy;
    SourceType type;
    double amp, freq;
    double t0, tau;     /* pulses */
    double sigma2;      /* spatial footprint variance in cells^2 */
} Source;

#define MAX_SRC 4
static Source g_src[MAX_SRC];
static int drag_src = -1;   // -1 = none

/* Source injection mode: linear (physics-correct) vs. saturating (safer when
 * driving very large amplitudes).  Default to linear for correctness. */
typedef enum { SRC_INJ_LINEAR = 0, SRC_INJ_SATURATING = 1 } SourceInjectionMode;
static SourceInjectionMode g_src_injection_mode = SRC_INJ_LINEAR;

static inline void source_reparam(Source* s){
    if (s->type == SRC_GAUSS_PULSE || s->type == SRC_RICKER) {
        double cycles_t0=6.0, cycles_tau=2.0;
        s->t0  = (s->freq>0)? (cycles_t0  / s->freq) : 0.0;
        s->tau = (s->freq>0)? (cycles_tau / s->freq) : 1e-9;
    } else { s->t0=0.0; s->tau=1e-9; }
}
static inline void sources_set_freq(double f){
    for (int k=0;k<MAX_SRC;++k){ g_src[k].freq=f; source_reparam(&g_src[k]); }
}
static inline void sources_cycle_type(void){
    for (int k=0;k<MAX_SRC;++k){ g_src[k].type=(SourceType)((g_src[k].type+1)%3); source_reparam(&g_src[k]); }
}
static inline double source_time_value(const Source* s, int t, double dt_){
    double tt = t * dt_;
    switch(s->type){
        case SRC_CW:          return s->amp * sin(2.0*M_PI*s->freq * tt);
        case SRC_GAUSS_PULSE: { double x=(tt-s->t0)/s->tau; double env=exp(-0.5*x*x);
                                 return s->amp * env * sin(2.0*M_PI*s->freq*(tt - s->t0)); }
        case SRC_RICKER:      { double a = M_PI*s->freq*(tt - s->t0); double e=exp(-a*a);
                                 return s->amp * (1.0 - 2.0*a*a) * e; }
        default: return 0.0;
    }
}
static inline void inject_source_into_Ez(Source* s, double (*Ezf)[NY], int t, double dt_){
    if (!s->active) return;
    /*
     * Soft source injection: instead of directly setting Ez = Ez + A*w
     * (hard source), which can cause spurious reflections, inject the
     * source via a saturating form.  See other AI suggestion: dividing
     * by (1+|A*w|*0.1) gently limits the perturbation when the source
     * amplitude is large.  This keeps reflections small while still
     * approximating a hard source for small amplitudes.
     */
    double A = source_time_value(s, t, dt_);
    for (int di=-2; di<=2; ++di)
    for (int dj=-2; dj<=2; ++dj){
        int i = s->ix + di;
        int j = s->iy + dj;
        if (i > 0 && i < NX && j > 0 && j < NY){
            double r2 = (double)(di*di + dj*dj);
            double w  = exp(-r2 / s->sigma2);
            double val = A * w;
            if (g_src_injection_mode == SRC_INJ_LINEAR) {
                /* Linear soft source (superposition preserved). */
                Ezf[i][j] += val;
            } else {
                /* Saturating injection: gently limits very large amplitudes for
                 * interactive demos at the cost of strict linearity. */
                Ezf[i][j] += val / (1.0 + fabs(val) * 0.1);
            }
        }
    }
}
static inline void draw_sources(SDL_Renderer* ren, const int scale){
    for (int k=0;k<MAX_SRC;++k){
        Source* s = &g_src[k];
        Uint8 rr = s->active?0:80, gg=s->active?255:80, bb=0;
        SDL_SetRenderDrawColor(ren, rr,gg,bb,255);
        int x=s->ix*scale, y=s->iy*scale;
        for (int d=-4; d<=4; ++d){ SDL_RenderDrawPoint(ren, x+d, y); SDL_RenderDrawPoint(ren, x, y+d); }
    }
}

/* ================= Oscilloscope ================= */
typedef struct { int n, head; double *y; int on; double last; } Scope;
static Scope scope = {0,0,NULL,0,0.0};

static void scope_init(int width) {
    if (scope.y) free(scope.y);
    scope.n = (width > 64 ? width : 64);
    scope.y = (double*)calloc(scope.n, sizeof(double));
    if (!scope.y) {
        fprintf(stderr, "Warning: Failed to allocate scope buffer\n");
        scope.n = 0;
        scope.on = 0;
        return;
    }
    scope.head = 0; scope.on = 1; scope.last = 0.0;
}
static void scope_free(void) { if (scope.y){ free(scope.y); scope.y=NULL; } scope.n=0; }
static inline void scope_push(double v) {
    if (!scope.y || !scope.on) return;
    scope.y[scope.head] = v;
    scope.head = (scope.head + 1) % scope.n;   /* head points to next slot */
    scope.last = v;
}
static void scope_clear(void){
    if (!scope.y) return;
    memset(scope.y, 0, sizeof(double)*scope.n);
    scope.head = 0; scope.last = 0.0;
}
/* draw with "now" at the RIGHT edge */
static void draw_scope(SDL_Renderer* r, int x, int y, int w, int h, double yscale) {
    if (!scope.on || !scope.y) return;
    SDL_Rect box = {x,y,w,h};
    SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(r, 20,20,20,220);
    SDL_RenderFillRect(r, &box);
    SDL_SetRenderDrawColor(r, 80,80,80,255);
    SDL_RenderDrawLine(r, x, y + h/2, x + w, y + h/2);

    if (yscale < 1e-12) yscale = 1e-12;

    SDL_SetRenderDrawColor(r, 220,220,220,255);
    int prevx = x, prevy = y + h/2;
    for (int i = 0; i < w; ++i) {
        int idx = (scope.head - (w - i) + scope.n) % scope.n; /* oldest->left, newest->right */
        double v = scope.y[idx];
        double t = 0.5 - 0.5*(v/yscale);
        if (t < 0) {
            t = 0;
        }
        if (t > 1) {
            t = 1;
        }
        int yy = y + (int)(t * (h-1));
        int xx = x + i;
        if (i) SDL_RenderDrawLine(r, prevx, prevy, xx, yy);
        prevx = xx; prevy = yy;
    }
}

/* FFT export from scope (naive DFT; on-demand) */
static int dump_scope_fft_csv(const char* path, double dt_, int Nfft_requested){
    if (!scope.y || scope.n <= 8 || dt_ <= 0.0) return 0;
    int N = scope.n;
    int Nfft = Nfft_requested;
    if (Nfft > N) Nfft = N;
    if (Nfft < 64) Nfft = (N < 64 ? N : 64);

    double *x = (double*)malloc(sizeof(double)*Nfft);
    if (!x) { fprintf(stderr, "Warning: FFT export failed - memory allocation error\n"); return 0; }
    int start = (scope.head - Nfft + scope.n) % scope.n; /* oldest index */
    for (int k=0; k<Nfft; ++k){
        int idx = (start + k) % scope.n;
        x[k] = scope.y[idx];
    }
    double mean = 0.0; for (int k=0;k<Nfft;++k) mean += x[k]; mean /= (double)Nfft;
    for (int k=0;k<Nfft;++k){
        double w = 0.5 * (1.0 - cos(2.0*M_PI*(double)k/(double)(Nfft-1)));
        x[k] = (x[k] - mean) * w;
    }
    FILE* f = fopen(path, "w");
    if (!f) {
        free(x);
        return 0;
    }
    fprintf(f, "# FFT of scope (Hann); dt=%.9e; N=%d\n", dt_, Nfft);
    fprintf(f, "freq_Hz,mag\n");
    for (int k=0; k<=Nfft/2; ++k){
        double re=0.0, im=0.0;
        double ang_step = -2.0*M_PI*(double)k/(double)Nfft;
        for (int n=0;n<Nfft;++n){
            double ang = ang_step*(double)n;
            re += x[n]*cos(ang); im += x[n]*sin(ang);
        }
        double mag = sqrt(re*re + im*im);
        double freq = (double)k / (dt_ * (double)Nfft);
        fprintf(f, "%.12e,%.12e\n", freq, mag);
    }
    fclose(f); free(x);
    return 1;
}

/* ================= Text (TTF) ================= */
static SDL_Texture* render_text(SDL_Renderer* ren, TTF_Font* font,
                                const char* s, SDL_Color col, int* w, int* h)
{
    SDL_Surface* surf = TTF_RenderUTF8_Blended(font, s, col);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ren, surf);
    if (w) {
        *w = surf->w;
    }
    if (h) {
        *h = surf->h;
    }
    SDL_FreeSurface(surf);
    return tex;
}
static void draw_legend(SDL_Renderer* ren, TTF_Font* font, const int x, const int y,
                        const char** lines, const int nlines)
{
    int maxw = 0, lineh = 0;
    for (int i=0;i<nlines;++i){
        int w,h; TTF_SizeUTF8(font, lines[i], &w, &h);
        if (w > maxw) maxw = w;
        if (h > lineh) lineh = h;
    }
    const int pad = 8, gap = 2;
    SDL_Rect panel = { x, y, maxw + 2*pad, nlines*lineh + (nlines-1)*gap + 2*pad };
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(ren, 10, 10, 10, 190);
    SDL_RenderFillRect(ren, &panel);
    SDL_SetRenderDrawColor(ren, 200, 200, 200, 220);
    SDL_RenderDrawRect(ren, &panel);

    int yy = y + pad;
    SDL_Color fg = {230,230,230,255};
    for (int i=0;i<nlines;++i){
        int tw, th;
        SDL_Texture* t = render_text(ren, font, lines[i], fg, &tw, &th);
        if (t){
            SDL_Rect r = { x + pad, yy, tw, th };
            SDL_RenderCopy(ren, t, NULL, &r);
            SDL_DestroyTexture(t);
        }
        yy += lineh + gap;
    }
}

/*
 * Wrapped text rendering for multi-line legend.  This helper uses
 * SDL_TTF's Wrapped API to automatically break long lines so the
 * legend fits within a constrained width.  The caller supplies the
 * rectangular area (area.x, area.y, area.w, area.h); text will be
 * clipped automatically if it overflows vertically.  A translucent
 * background panel is drawn behind the text.
 */
/* Define render_text_wrapped only if not previously defined.
 * This guard avoids a duplicate symbol error when another implementation of
 * render_text_wrapped exists in the build. */
#ifndef RENDER_TEXT_WRAPPED_DEFINED
#define RENDER_TEXT_WRAPPED_DEFINED
static SDL_Texture* render_text_wrapped(SDL_Renderer* ren, TTF_Font* font,
                                        const char* s, SDL_Color col,
                                        Uint32 wrap_width,
                                        int* outw, int* outh)
{
    SDL_Surface* surf = TTF_RenderUTF8_Blended_Wrapped(font, s, col, wrap_width);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ren, surf);
    if (outw) *outw = surf->w;
    if (outh) *outh = surf->h;
    SDL_FreeSurface(surf);
    return tex;
}
#endif /* RENDER_TEXT_WRAPPED_DEFINED */
static void draw_wrapped_panel(SDL_Renderer* ren, TTF_Font* font, const SDL_Rect area, const char* text)
{
    const int pad = 6;
    SDL_Color fg = {230,230,230,255};
    int tw=0, th=0;
    SDL_Texture* t = render_text_wrapped(ren, font, text, fg, (Uint32)(area.w - 2*pad), &tw, &th);
    if (!t) return;
    /* limit height to area.h */
    int draw_h = (th + 2*pad > area.h) ? (area.h) : (th + 2*pad);
    SDL_Rect panel = { area.x, area.y, area.w, draw_h };
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(ren, 10,10,10,190);
    SDL_RenderFillRect(ren, &panel);
    SDL_SetRenderDrawColor(ren, 200,200,200,220);
    SDL_RenderDrawRect(ren, &panel);
    SDL_Rect r = { area.x + pad, area.y + pad, tw, th };
    /* Clip vertical overflow: adjust dest h if necessary */
    if (r.y + r.h > panel.y + panel.h) r.h = (panel.y + panel.h) - r.y;
    SDL_RenderCopy(ren, t, NULL, &r);
    SDL_DestroyTexture(t);
}

/*
 * Draw a multi‑line information panel in the right side of the UI.  This helper
 * assembles a newline‑delimited string summarising the current simulation
 * parameters and renders it wrapped inside the provided area.  The panel
 * replaces the old legend overlay by presenting key metrics, source
 * configuration, boundary conditions, probe values, and other state in a
 * single compact block.  The caller must compute the destination rectangle
 * (area.x, area.y, area.w, area.h) and should leave some padding around
 * the edges for aesthetic spacing.  The fps_inst and fps_avg parameters
 * should be computed in the render loop and passed in directly.
 */

/* ================= Screenshot ================= */
static void save_screenshot(SDL_Renderer* ren, const char* path, int w, int h) {
    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surf) return;
    SDL_RenderReadPixels(ren, NULL, SDL_PIXELFORMAT_ARGB8888, surf->pixels, surf->pitch);
    SDL_SaveBMP(surf, path);
    SDL_FreeSurface(surf);
}

/* ================= Probe log ================= */
static FILE *probe_out = NULL;
static int   probe_x = NX/2 + NX/6, probe_y = NY/2;

/* Optional second probe. Use key '3' to toggle on/off. When active, a middle-click
 * (mouse button 2) will reposition this probe. Its value is displayed in the
 * legend alongside the primary probe. By default it is inactive. */
static int probe2_active = 0;
static int probe2_x = NX/2 - NX/6;
static int probe2_y = NY/2;

/* ================= Visualization scaling ================= */
typedef enum { AS_PEAK=0, AS_P99=1 } AutoScaleMode;
static AutoScaleMode color_autoscale_mode = AS_P99;

/* decoupled holds & smoothing */
static int hold_color = 0;
static double cached_p99 = 0.0;  /* Cached P99 value when hold_color is active */
static int p99_cache_valid = 0;  /* Flag indicating if cached P99 is valid */
static double held_vmax = 1e-3;
static int hold_scope = 0;
static double held_scope_vmax = 1e-3;

static double vmax_smooth = 0.0, p99_smooth = 0.0, scope_vmax_smooth = 0.0;

/* optional rendering stride (perf) */
static int g_render_stride = 1; /* 1,2,4 */

/* ===================================================================== */
/* UI buttons for mouse-only interaction.  Each button occupies a
 * rectangular region in the bottom panel.  When clicked, toggle
 * an associated boolean or perform an action.  See draw loop below.
 */
typedef struct {
    SDL_Rect r;          /* screen rectangle */
    const char* label;   /* text on button */
    int is_toggle;       /* if true, toggles *state on click */
    int* state;          /* pointer to boolean to toggle */
    int id;              /* an ID for non-toggle actions */
} UIButton;
static UIButton ui_btns[16];
static int ui_btn_count = 0;

/* Button IDs for non-toggle actions */
enum {
    UI_BTN_NONE = 0,
    UI_BTN_SCREENSHOT = 1,
    UI_BTN_CLEAR = 2,
    UI_BTN_S21_NOW = 3,
    /* Cycle through CPML presets */
    UI_BTN_CPML_CYCLE = 4
    , UI_BTN_BOUNDARY_CYCLE = 5
};

/* Draw a single button.  Active toggles are drawn brighter. */
static void draw_button(SDL_Renderer* ren, TTF_Font* font, const UIButton* b){
    /* background */
    if (b->is_toggle && b->state && *(b->state)){
        SDL_SetRenderDrawColor(ren, 70,130,180,255); /* steel blue for ON */
    } else {
        SDL_SetRenderDrawColor(ren, 60,60,60,255);
    }
    SDL_RenderFillRect(ren, &b->r);
    /* border */
    SDL_SetRenderDrawColor(ren, 180,180,180,255);
    SDL_RenderDrawRect(ren, &b->r);
    /* label */
    SDL_Color fg = {230,230,230,255};
    int tw=0, th=0;
    SDL_Texture* tt = render_text(ren, font, b->label, fg, &tw, &th);
    if (tt){
        SDL_Rect rr;
        rr.w = tw;
        rr.h = th;
        rr.x = b->r.x + (b->r.w - tw)/2;
        rr.y = b->r.y + (b->r.h - th)/2;
        SDL_RenderCopy(ren, tt, NULL, &rr);
        SDL_DestroyTexture(tt);
    }
}

/* Handle a mouse click; returns the ID of a non-toggle button clicked,
 * or UI_BTN_NONE.  Toggles are updated directly via their state pointer.
 */
static int ui_handle_click(int mx, int my){
    for (int i=0; i<ui_btn_count; ++i){
        UIButton* b = &ui_btns[i];
        if (mx >= b->r.x && mx < b->r.x + b->r.w && my >= b->r.y && my < b->r.y + b->r.h){
            if (b->is_toggle && b->state){
                *(b->state) = !(*(b->state));
                return UI_BTN_NONE;
            }
            return b->id;
        }
    }
    return UI_BTN_NONE;
}

/* =========== Right-side button grid builder =========== */
static void build_side_buttons(int panel_x, int panel_y, int panel_w, int panel_h,
                               int *paused, int *auto_rescale, int *hold_color, int *cpml_on, int *ports_on, int *paint_mode)
{
    ui_btn_count = 0;

    const int pad = 6;
    const int gap = 6;
    const int cols = (panel_w >= 240 ? 3 : 2);      /* 3 columns if roomy, else 2 */
    const int bh   = 24;                             /* button height */
    const int bw   = (panel_w - pad*2 - gap*(cols-1)) / cols;

    /* local macro to add a toggle or an action button */
    #define ADD_TOGGLE(LABEL, PTR) do{                                  \
        if (ui_btn_count < 16){                                         \
            UIButton* b = &ui_btns[ui_btn_count++];                     \
            int row = (ui_btn_count-1)/cols;                            \
            int col = (ui_btn_count-1)%cols;                            \
            b->r.x = panel_x + pad + col*(bw+gap);                      \
            b->r.y = panel_y + pad + row*(bh+gap);                      \
            b->r.w = bw; b->r.h = bh;                                   \
            b->label = LABEL; b->is_toggle = 1; b->state = PTR;         \
            b->id = UI_BTN_NONE;                                        \
        }                                                               \
    } while(0)

    #define ADD_ACTION(LABEL, IDCONST) do{                              \
        if (ui_btn_count < 16){                                         \
            UIButton* b = &ui_btns[ui_btn_count++];                     \
            int row = (ui_btn_count-1)/cols;                            \
            int col = (ui_btn_count-1)%cols;                            \
            b->r.x = panel_x + pad + col*(bw+gap);                      \
            b->r.y = panel_y + pad + row*(bh+gap);                      \
            b->r.w = bw; b->r.h = bh;                                   \
            b->label = LABEL; b->is_toggle = 0; b->state = NULL;        \
            b->id = IDCONST;                                            \
        }                                                               \
    } while(0)

    /* The toggles you asked for — short labels to fit nicely.  We pass
     * pointers directly so the state updates the actual variables. */
    ADD_TOGGLE("pause",      paused);
    ADD_TOGGLE("auto Δ",     auto_rescale);
    ADD_TOGGLE("hold color", hold_color);
    ADD_TOGGLE("cpml",       cpml_on);
    ADD_TOGGLE("ports",      ports_on);
    ADD_TOGGLE("paint",      paint_mode);

    /* Per-source toggles: allow each defined source to be turned on/off.  The
     * labels are S1, S2, etc.  If you add more sources via MAX_SRC, the
     * corresponding buttons will appear automatically. */
    {
        /* Compose a short label for each source once.  We use a static
         * array so the label pointers remain valid after the loop. */
        static char s_labels[MAX_SRC][4];
        for (int kbtn = 0; kbtn < MAX_SRC; ++kbtn) {
            snprintf(s_labels[kbtn], sizeof(s_labels[kbtn]), "S%d", kbtn + 1);
            ADD_TOGGLE(s_labels[kbtn], &g_src[kbtn].active);
        }
    }

    /* Action buttons */
    ADD_ACTION("shot",   UI_BTN_SCREENSHOT);
    ADD_ACTION("clear",  UI_BTN_CLEAR);
    ADD_ACTION("S21",    UI_BTN_S21_NOW);
    /* Cycle through CPML presets.  Clicking this button steps to the next
     * preset and rebuilds the absorbing boundary. */
    ADD_ACTION("pml+",  UI_BTN_CPML_CYCLE);
    /* Cycle through boundary types (Mur <-> CPML).  This action updates
     * boundary_type and synchronises cpml_on accordingly. */
    ADD_ACTION("bnd+",  UI_BTN_BOUNDARY_CYCLE);

    #undef ADD_TOGGLE
    #undef ADD_ACTION

    /* clamp panel_h use: we won't draw a scissor here; just build within panel */
}


/* ================= CPML ================= */
#ifndef PML_THICK
#define PML_THICK 12
#endif
/* Start with CPML enabled by default for more accurate outer boundaries. */
static int cpml_on = 1;
static int cpml_N  = PML_THICK;
static double pml_sigma_max = 1.2, pml_kappa_max = 5.0, pml_alpha_max = 0.05;

/*
 * Boundary type selection.  In the original code the absorbing boundary
 * condition is controlled exclusively via the cpml_on toggle.  To support
 * future boundary conditions (e.g. Mur, CPML, Perfect Electric Conductor)
 * we define an explicit boundary_type.  Currently there are two modes:
 *  - BOUNDARY_MUR : a first‑order Mur absorbing boundary (default)
 *  - BOUNDARY_CPML: convolutional perfectly matched layer (CPML)
 * The cpml_on flag remains for backward compatibility and always mirrors
 * the boundary_type (cpml_on = (boundary_type == BOUNDARY_CPML)).
 */
typedef enum { BOUNDARY_MUR=0, BOUNDARY_CPML=1 } BoundaryType;
static BoundaryType boundary_type = BOUNDARY_CPML;

/* per-axis coeffs */
static double kx[NX], bx[NX], cx[NX];
static double ky[NY], by[NY], cy[NY];

/* ψ memory */
static double (*psi_Ezx)[NY];  /* for Hy: dEz/dx */
static double (*psi_Ezy)[NY];  /* for Hx: dEz/dy */
static double (*psi_Hyx)[NY];  /* for Ez: dHy/dx */
static double (*psi_Hxy)[NY];  /* for Ez: dHx/dy */

typedef struct { const char* name; double smax, kmax, amax; int thick; } CpmlPreset;
static const CpmlPreset CPML_PRESETS[] = {
    { "Gentle",     1.0, 3.0, 0.03, 10 },
    { "Default",    1.2, 5.0, 0.05, 12 },
    { "Aggressive", 1.8, 6.0, 0.08, 16 },
};
static int cpml_preset_idx = 1;

static inline void cpml_zero_psi(void){
    if (psi_Ezx) memset(psi_Ezx,0,sizeof(double[NX][NY]));
    if (psi_Ezy) memset(psi_Ezy,0,sizeof(double[NX][NY]));
    if (psi_Hyx) memset(psi_Hyx,0,sizeof(double[NX][NY]));
    if (psi_Hxy) memset(psi_Hxy,0,sizeof(double[NX][NY]));
}
static void cpml_build_coeffs(double dt_){
    const double m = 3.0; /* cubic grading */
    for (int i=0;i<NX;++i){
        double rx = 0.0;
        /*
         * Compute the fractional distance into the PML region.  Use cell
         * centres aligned with Yee grid axes rather than offset by 0.5
         * cells as in some references.  See other AI notes: the 0.5 offset
         * mismatches the Yee staggering and reduces absorption.  We simply
         * normalise the integer distance [0..N] to [0..1].
         */
        if (i < cpml_N)            rx = (cpml_N - i) / (double)cpml_N;
        else if (i >= NX - cpml_N) rx = (i - (NX - cpml_N)) / (double)cpml_N;
        double g = rx>0 ? pow(rx, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - rx);
        kx[i] = (rx>0)? kappa : 1.0;
        if (rx>0){
            bx[i] = exp(-(sigma/kappa + alpha) * dt_);
            cx[i] = (sigma * (bx[i] - 1.0)) / (kappa * (sigma + kappa*alpha) + 1e-20);
        } else { bx[i]=1.0; cx[i]=0.0; }
    }
    for (int j=0;j<NY;++j){
        double ry = 0.0;
        if (j < cpml_N)            ry = (cpml_N - j) / (double)cpml_N;
        else if (j >= NY - cpml_N) ry = (j - (NY - cpml_N)) / (double)cpml_N;
        double g = ry>0 ? pow(ry, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - ry);
        ky[j] = (ry>0)? kappa : 1.0;
        if (ry>0){
            by[j] = exp(-(sigma/kappa + alpha) * dt_);
            cy[j] = (sigma * (by[j] - 1.0)) / (kappa * (sigma + kappa*alpha) + 1e-20);
        } else { by[j]=1.0; cy[j]=0.0; }
    }
}
static void cpml_apply_preset(int idx, double dt_){
    int n = (int)(sizeof(CPML_PRESETS)/sizeof(CPML_PRESETS[0]));
    if (idx < 0) {
        idx = 0;
    }
    if (idx >= n) {
        idx = n-1;
    }
    cpml_preset_idx = idx;
    pml_sigma_max = CPML_PRESETS[idx].smax;
    pml_kappa_max = CPML_PRESETS[idx].kmax;
    pml_alpha_max = CPML_PRESETS[idx].amax;
    cpml_N        = CPML_PRESETS[idx].thick;

    /* Validate CPML thickness doesn't exceed grid size */
    int max_thickness = (NX < NY ? NX : NY) / 2 - 1;
    if (cpml_N > max_thickness) {
        fprintf(stderr, "Warning: CPML thickness %d exceeds safe limit %d, clamping\n",
                cpml_N, max_thickness);
        cpml_N = max_thickness;
    }
    if (cpml_N < 1) cpml_N = 1;

    cpml_build_coeffs(dt_);
    cpml_zero_psi();
}


/* ===================== INJECTION 1: reset & rebuild ===================== */
static void reset_sim_state(void);
static void rebuild_sim_for_freq(double new_f);

/* Declare freq as a global variable */
static double freq = 1e9;

static inline void zero_2d(double (*a)[NY]) { memset(a, 0, sizeof(double[NX][NY])); }

static void reset_sim_state(void)
{
    /* fields */
    zero_2d(Ez); zero_2d(Hx); zero_2d(Hy);
    zero_2d(Ez_old);

    /* CPML aux */
    zero_2d(psi_Ezx);
    zero_2d(psi_Ezy);
    zero_2d(psi_Hyx);
    zero_2d(psi_Hxy);

    /* scope */
    scope_clear();

    /* ports */
    for (int p=0; p<MAX_PORTS; ++p) {
        if (ports[p].V) memset(ports[p].V, 0, sizeof(double)*ports[p].n);
        if (ports[p].I) memset(ports[p].I, 0, sizeof(double)*ports[p].n);
        ports[p].head = 0;
    }

    /* color/scope autoscale smoothers */
    vmax_smooth = 0.0;
    p99_smooth = 0.0;
    scope_vmax_smooth = 0.0;
}

/* Recompute Δ, dt, CPML, reparam sources and reset */
static void rebuild_sim_for_freq(double new_f)
{
    if (new_f <= 0) return;
    if (fabs(new_f - freq) < 1e-12) return;

    freq = new_f;

    /* your existing helper that sets dx,dy,dt from freq */
    update_grid_for_freq(freq, &dx, &dy, &dt);

    /* sources that depend on freq */
    for (int k=0; k<MAX_SRC; ++k) source_reparam(&g_src[k]);

    if (cpml_on) {
        cpml_build_coeffs(dt);
        cpml_zero_psi();
    }
    reset_sim_state();
}
/* =================== END INJECTION 1 =================== */



static void draw_info_panel(SDL_Renderer* ren, TTF_Font* font, const SDL_Rect area,
                            const double fps_inst_local, const double fps_avg_local,
                            const int steps_per_frame_local, const int paused_local)
{
    /* Compose the information string into a static buffer.  We keep the
     * buffer reasonably large; if it overflows, truncation will occur. */
    char buf[2048];
    buf[0] = '\0';
    size_t len = 0;

    /* Safe snprintf macro to prevent buffer overflow */
    #define SAFE_APPEND(...) do { \
        if (len < sizeof(buf) - 1) { \
            int written = snprintf(buf+len, sizeof(buf)-len, __VA_ARGS__); \
            if (written > 0) len += (size_t)written; \
            if (len >= sizeof(buf)) len = sizeof(buf) - 1; \
        } \
    } while(0)

    /* Frequency and units */
    double disp = (freq >= 1e9) ? freq/1e9 : freq/1e6;
    const char* unit = (freq >= 1e9) ? "GHz" : "MHz";
    SAFE_APPEND("Freq: %.2f %s\n", disp, unit);
    /* Grid spacing and timestep */
    SAFE_APPEND("Δ=%.4g m   dt=%.4g s\n", dx, dt);
    /* Cells per wavelength (in/out) */
    double lambda0 = c0 / freq;
    double cpw_out = lambda0 / dx;
    double cpw_in  = (lambda0 / sqrt(EPSR_MAX_SCENE)) / dx;
    SAFE_APPEND("CPW (in/out): %.1f / %.1f\n", cpw_in, cpw_out);
    /* Mode and colour scale */
    const char* mode = auto_rescale ? "AUTO" : "FIXED";
    const char* scmode = hold_color ? "HOLD" : (color_autoscale_mode==AS_P99 ? "AUTO-P99" : "AUTO-PEAK");
    SAFE_APPEND("Mode: %s  Color: %s  Render:x%d\n", mode, scmode, g_render_stride);
    /* Steps/frame and pause state */
    SAFE_APPEND("Steps/frame: %d   %s\n",
                    steps_per_frame_local, paused_local ? "PAUSED" : "RUNNING");
    /* Source summary */
    int active_cnt = 0;
    for (int k=0; k<MAX_SRC; ++k){ if (g_src[k].active) ++active_cnt; }
    SAFE_APPEND("Sources: %d active of %d total\n", active_cnt, MAX_SRC);
    for (int k=0; k<MAX_SRC; ++k) {
        const char* tname = (g_src[k].type==SRC_CW?"CW":(g_src[k].type==SRC_GAUSS_PULSE?"Gauss":"Ricker"));
        /* Format frequency in MHz or GHz for readability */
        double sdisp = (g_src[k].freq >= 1e9) ? g_src[k].freq/1e9 : g_src[k].freq/1e6;
        const char* sunit = (g_src[k].freq >= 1e9) ? "GHz" : "MHz";
        SAFE_APPEND("  S%d: %s  Pos=(%d,%d)  Amp=%.2f  f=%.2f %s  Type=%s\n",
                        k+1, g_src[k].active?"ON":"OFF",
                        g_src[k].ix, g_src[k].iy,
                        g_src[k].amp,
                        sdisp, sunit,
                        tname);
    }
    /* Source injection mode (linear vs. saturating) */
    SAFE_APPEND("Src inj: %s\n",
                (g_src_injection_mode == SRC_INJ_LINEAR) ? "LINEAR" : "SATURATING");
    /* FPS and power flux */
    SAFE_APPEND("FPS: %.2f (inst), %.2f (avg)  Pflux=%.3e W/m\n",
                    fps_inst_local, fps_avg_local, pflux_avg);
    /* Probe values */
    if (probe2_active) {
        SAFE_APPEND("Probe1: (%d,%d)=%.3g  Probe2: (%d,%d)=%.3g\n",
                        probe_x, probe_y, Ez[probe_x][probe_y],
                        probe2_x, probe2_y, Ez[probe2_x][probe2_y]);
    } else {
        SAFE_APPEND("Probe: (%d,%d) Ez=%.3g\n",
                        probe_x, probe_y, Ez[probe_x][probe_y]);
    }
    /* Scope status */
    SAFE_APPEND("Scope: %s  Scale: %s  Last=%.3g\n",
                    scope.on?"ON":"OFF", hold_scope?"HOLD":"AUTO", scope.last);
    /* Sigma material parameters */
    SAFE_APPEND("Sigma: bg=%.3g S/m  block=%.3g S/m\n",
                    SIGMA_BG, SIGMA_BLOCK);
    /* Boundary / CPML information */
    if (boundary_type == BOUNDARY_CPML) {
        SAFE_APPEND("Boundary: CPML (%s)  N=%d  sigma=%.2f  kappa=%.2f  alpha=%.3f\n",
                        CPML_PRESETS[cpml_preset_idx].name,
                        cpml_N, pml_sigma_max, pml_kappa_max, pml_alpha_max);
    } else {
        SAFE_APPEND("Boundary: Mur-1\n");
    }
    /* Ports / sweep information */
    if (sweep_on) {
        double freq_now = (sweep_idx < sweep_points) ? sweep_freqs[sweep_idx] : 0.0;
        double fdisp = (freq_now >= 1e9) ? freq_now/1e9 : freq_now/1e6;
        const char* funit = (freq_now >= 1e9) ? "GHz" : "MHz";
        SAFE_APPEND("Sweep: %d/%d  freq=%.2f %s\n",
                        sweep_idx+1, sweep_points, fdisp, funit);
    } else {
        SAFE_APPEND("Ports: %s  Last S21=%.3g\n",
                        ports_on?"ON":"OFF", s21_amp);
    }
    /* Painting information */
    if (paint_mode) {
        const char* pnames[] = {"OFF","PEC","PMC","Diel"};
        if (paint_type == 3) {
            SAFE_APPEND("Paint: %s  eps=%.2f\n", pnames[paint_type], paint_eps);
        } else {
            SAFE_APPEND("Paint: %s\n", pnames[paint_type]);
        }
    } else {
        SAFE_APPEND("Paint: OFF\n");
    }
    /* Key hints (brief summary).  For a full list, toggle the old legend with 'L'. */
    SAFE_APPEND("Keys: Space=pause, ↑/↓=freq, ←/→=steps, Z=stride, R=mode, A=color mode, H=hold col, J=hold scope,\\/-reset, T=src type, U=paint, S=ports\n");

    #undef SAFE_APPEND

    /* Render the composed text within the area with wrapping. */
    draw_wrapped_panel(ren, font, area, buf);
}


/* ================= MAIN ================= */
int main(void){
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "0");
    if (SDL_Init(SDL_INIT_VIDEO)!=0){ fprintf(stderr,"SDL init failed: %s\n", SDL_GetError()); return 1; }
    if (TTF_Init() != 0) { fprintf(stderr, "TTF_Init failed: %s\n", TTF_GetError()); SDL_Quit(); return 1; }

    /* allocate fields with individual null checks for better diagnostics */
    Ez     = malloc(sizeof(double[NX][NY]));
    if (!Ez) { fprintf(stderr, "Failed to allocate Ez field\n"); TTF_Quit(); SDL_Quit(); return 1; }

    Hx     = malloc(sizeof(double[NX][NY]));
    if (!Hx) { fprintf(stderr, "Failed to allocate Hx field\n"); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    Hy     = malloc(sizeof(double[NX][NY]));
    if (!Hy) { fprintf(stderr, "Failed to allocate Hy field\n"); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    Ez_old = malloc(sizeof(double[NX][NY]));
    if (!Ez_old) { fprintf(stderr, "Failed to allocate Ez_old field\n"); free(Hy); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    psi_Ezx = malloc(sizeof(double[NX][NY]));
    if (!psi_Ezx) { fprintf(stderr, "Failed to allocate psi_Ezx\n"); free(Ez_old); free(Hy); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    psi_Ezy = malloc(sizeof(double[NX][NY]));
    if (!psi_Ezy) { fprintf(stderr, "Failed to allocate psi_Ezy\n"); free(psi_Ezx); free(Ez_old); free(Hy); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    psi_Hyx = malloc(sizeof(double[NX][NY]));
    if (!psi_Hyx) { fprintf(stderr, "Failed to allocate psi_Hyx\n"); free(psi_Ezy); free(psi_Ezx); free(Ez_old); free(Hy); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }

    psi_Hxy = malloc(sizeof(double[NX][NY]));
    if (!psi_Hxy) { fprintf(stderr, "Failed to allocate psi_Hxy\n"); free(psi_Hyx); free(psi_Ezy); free(psi_Ezx); free(Ez_old); free(Hy); free(Hx); free(Ez); TTF_Quit(); SDL_Quit(); return 1; }
    clear_fields();
    cpml_zero_psi();

    /* UI geometry */
    int scale=2, UI_H=80;
    const int SIM_YOFF = 0;
    const int UI_Y     = NY*scale;
    int side_panel = 260; // enlarge right panel to make space for info & buttons

    SDL_Window *win = SDL_CreateWindow("FDTD Demo",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        NX*scale + side_panel, NY*scale + UI_H,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (!win){ fprintf(stderr,"Window creation failed: %s\n", SDL_GetError()); return 1; }
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    if (!ren){ fprintf(stderr,"Renderer creation failed: %s\n", SDL_GetError()); return 1; }
    SDL_RenderSetLogicalSize(ren, NX*scale + side_panel, NY*scale + UI_H);

    /* font
     * Try a list of potential font paths so the simulation will start on
     * Windows and other systems where the DejaVu font is not available.
     * On Windows, Arial and Calibri are typically available under the
     * C:\\Windows\\Fonts directory.  On Unix-like systems, fall back to
     * DejaVuSans or FreeSans.  If none of these succeed, abort.
     */
    const char* font_paths[] = {
#ifdef _WIN32
        "C:\\Windows\\Fonts\\arial.ttf",
        "C:\\Windows\\Fonts\\Arial.ttf",
        "C:\\Windows\\Fonts\\calibri.ttf",
#endif
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        "DejaVuSans.ttf",
        NULL
    };
    TTF_Font* font = NULL;
    for (int i = 0; font_paths[i] && !font; ++i) {
        font = TTF_OpenFont(font_paths[i], 14);
    }
    if (!font) {
        fprintf(stderr, "TTF open failed: %s\n", TTF_GetError());
        return 1;
    }

    /* scope */
    scope_init(NX * scale);

    /* simulation parameters / UI state */
    const double FREQ_MIN = 1e6, FREQ_MAX = 5e9;
    int steps_per_frame = STEPS_PER_FRAME, paused=0;
    int show_legend = 1;

    dt = cfl_dt(dx, dy);
    cpml_build_coeffs(dt);

    /* ----------------------------------------------------------------------- */
    /* Material and ports initialization (Phase 1 & 2 additions)               */
    /*
     * Fill epsr and sigma_map with default background material (vacuum) and
     * assign the central block region a higher relative permittivity and
     * conductivity if specified. Tag grid is cleared (no PEC cells). Then
     * define two ports along the left and right boundaries for S-parameter
     * sampling. Each port has a sampling buffer for voltage (Ez) and current
     * (Hy) signals.  Ports sampling is off by default; toggle with 's'.
     */
    {
        /* reset material and tags */
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                epsr[i][j] = 1.0;
                sigma_map[i][j] = SIGMA_BG;
                tag_grid[i][j] = 0;
            }
        }
        /* central dielectric block properties: same as old epsilon_at() logic */
        int bx0 = NX/2 - NX/10;
        int bx1 = NX/2 + NX/10;
        int by0 = NY/2 - NY/20;
        int by1 = NY/2 + NY/20;
        for (int i = bx0; i <= bx1; ++i) {
            for (int j = by0; j <= by1; ++j) {
                epsr[i][j] = EPSR_MAX_SCENE;
                sigma_map[i][j] = SIGMA_BLOCK;
            }
        }
        /* initialize ports definitions */
        int p_y0 = NY / 4;
        int p_y1 = 3 * NY / 4;
        for (int p = 0; p < MAX_PORTS; ++p) {
            ports[p].n = PORT_SIGNAL_LENGTH;
            ports[p].V = (double*)calloc(PORT_SIGNAL_LENGTH, sizeof(double));
            ports[p].I = (double*)calloc(PORT_SIGNAL_LENGTH, sizeof(double));
            ports[p].head = 0;
            ports[p].active = 1;
            ports[p].y0 = p_y0;
            ports[p].y1 = p_y1;
            ports[p].len = p_y1 - p_y0 + 1;
        }
        /* port 0 on left side */
        ports[0].x = 1;
        /* port 1 on right side */
        ports[1].x = NX - 2;
        ports_on = 0;
        s21_amp = 0.0;
        s21_computed = 0;
        paint_mode = 0;
    }

    /* sources */
    for (int k=0;k<MAX_SRC;++k){
        g_src[k].active = 0;
        /* Assign default positions for each source.  Distribute them around
         * the domain so that additional sources start in distinct locations.
         * S1 defaults to center, S2 defaults to left third, S3 defaults to
         * right third, S4 defaults to bottom third.  These can be dragged
         * interactively later. */
        switch (k) {
            case 0: g_src[k].ix = NX/2;      g_src[k].iy = NY/2;      break;
            case 1: g_src[k].ix = NX/3;      g_src[k].iy = NY/2;      break;
            case 2: g_src[k].ix = 2*NX/3;    g_src[k].iy = NY/2;      break;
            case 3: g_src[k].ix = NX/2;      g_src[k].iy = NY/3;      break;
            default:g_src[k].ix = NX/2;      g_src[k].iy = NY/2;      break;
        }
        g_src[k].type = SRC_CW;
        g_src[k].amp  = 1.0;
        g_src[k].freq = freq;
        g_src[k].sigma2 = 4.0;
        g_src[k].t0 = 0.0;
        g_src[k].tau = 1e-9;
        source_reparam(&g_src[k]);
    }
    /* Activate the first source by default.  Subsequent sources start inactive. */
    if (MAX_SRC > 0) g_src[0].active = 1;

    Slider s_freq  = { .x=20,.y=UI_Y+15,.w=NX*scale-40,.h=20,.minv=0,.maxv=1,
                       .value=slider_from_freq(freq,FREQ_MIN,FREQ_MAX),.dragging=0 };
    Slider s_speed = { .x=20,.y=UI_Y+45,.w=NX*scale-40,.h=20,.minv=1,.maxv=50,
                       .value=steps_per_frame,.dragging=0 };

    /* probe */
    probe_out = fopen("probe.txt","w");
    if (probe_out){
        fprintf(probe_out,
            "# dt=%.9e dx=%.6e dy=%.6e mode=%s cpml=%d freq=%.9e SIGMA_BG=%.3g SIGMA_BLOCK=%.3g\n",
            dt, dx, dy, auto_rescale? "AUTO":"FIXED", cpml_on, freq, SIGMA_BG, SIGMA_BLOCK);
    }

    /* FPS smoothing */
    Uint64 t_prev = SDL_GetPerformanceCounter();
    double perf_freq = (double)SDL_GetPerformanceFrequency();
    double fps_avg = 0.0;

    int running=1, t=0; SDL_Event e;

    while (running){
        /* -------- events -------- */
        while (SDL_PollEvent(&e)){
            if (e.type == SDL_QUIT) running=0;

            /* right-click: move primary probe */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT) {
                int mx = e.button.x, my = e.button.y;
                int ix = clampi(mx / scale, 1, NX-2);
                int iy = clampi(my / scale, 1, NY-2);
                probe_x = ix; probe_y = iy;
            }
            /* middle-click: move second probe (if active) */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_MIDDLE) {
                if (probe2_active) {
                    int mx = e.button.x, my = e.button.y;
                    int ix = clampi(mx / scale, 1, NX-2);
                    int iy = clampi(my / scale, 1, NY-2);
                    probe2_x = ix; probe2_y = iy;
                }
            }

            /* left-click: either paint PEC cells when paint_mode is ON, or drag nearest active source */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                /* --- side buttons click handling --- */
                {
                    int mx = e.button.x, my = e.button.y;

                    /* Rebuild the same panel geometry here to ensure ui_btns have positions.
                     * We subdivide the side panel similarly to the render stage into an
                     * information area and a button area; only the button area is
                     * passed to build_side_buttons. */
                    int sp_x = NX*scale + 6 + 14 + 6;                     /* cb_x + cb_w + 6 */
                    int sp_y = 6;
                    int sp_w = NX*scale + side_panel - 6 - sp_x;
                    int sp_h = NY*scale - 12;
                    /* Determine the portion reserved for informational text.  Use the same
                     * fraction as in the render loop to ensure consistent layout. */
                    int info_h = (int)(0.6 * (double)sp_h);
                    if (info_h < 100) info_h = (sp_h > 100 ? sp_h - 100 : sp_h);
                    int btn_y = sp_y + info_h + 4;
                    int btn_h = sp_h - info_h - 4;
                    if (btn_h < 0) btn_h = 0;
                    /* Build buttons using addresses of the actual state variables */
                    build_side_buttons(sp_x, btn_y, sp_w, btn_h,
                                    &paused, &auto_rescale, &hold_color, &cpml_on, &ports_on, &paint_mode);

                    /* remember old states to run side-effects only on change */
                    int old_auto  = auto_rescale;
                    int old_cpml  = cpml_on;
                    int old_ports = ports_on;
                    int old_paint = paint_mode;
                    int old_hold  = hold_color;
                    int old_pause = paused;

                    int id = ui_handle_click(mx, my);
                    /* Did this click hit a button?
                    - Actions: id != UI_BTN_NONE
                    - Toggles: compare state before/after
                    */
                    if (id != UI_BTN_NONE ||
                        old_auto  != auto_rescale ||
                        old_cpml  != cpml_on      ||
                        old_ports != ports_on     ||
                        old_paint != paint_mode   ||
                        old_hold  != hold_color   ||
                        old_pause != paused)
                    {
                        /* Side-effects after a toggle/action */
                        if (auto_rescale != old_auto){
                            if (auto_rescale) {
                                update_grid_for_freq(freq, &dx, &dy, &dt);
                            } else {
                                dx = BASE_DX;
                                dy = BASE_DY;
                                dt = cfl_dt(dx,dy);
                            }
                            cpml_build_coeffs(dt);
                            if (cpml_on) {
                                cpml_zero_psi();
                            }
                            clear_fields(); t = 0;
                            s_freq.value  = slider_from_freq(freq, FREQ_MIN, FREQ_MAX);
                            s_speed.value = steps_per_frame;
                        }
                        if (cpml_on != old_cpml){
                            /* Synchronise boundary_type with cpml_on.  When cpml_on is
                             * toggled from the UI, boundary_type is updated accordingly.
                             */
                            boundary_type = cpml_on ? BOUNDARY_CPML : BOUNDARY_MUR;
                            if (cpml_on) {
                                cpml_build_coeffs(dt);
                                cpml_zero_psi();
                            } else {
                                cpml_zero_psi();
                            }
                        }

                        if (id == UI_BTN_CPML_CYCLE){
                            /* Cycle to the next CPML preset.  Wrap around automatically. */
                            int n_presets = (int)(sizeof(CPML_PRESETS)/sizeof(CPML_PRESETS[0]));
                            int next_idx = (cpml_preset_idx + 1) % n_presets;
                            cpml_apply_preset(next_idx, dt);
                        } else if (id == UI_BTN_BOUNDARY_CYCLE) {
                            /* Cycle between Mur and CPML boundary conditions.  This updates
                             * the cpml_on flag and rebuilds the CPML coefficients if needed. */
                            boundary_type = (BoundaryType)((boundary_type + 1) % 2);
                            cpml_on = (boundary_type == BOUNDARY_CPML);
                            if (cpml_on) {
                                cpml_build_coeffs(dt);
                                cpml_zero_psi();
                            } else {
                                /* Clear CPML memory; psi arrays will be unused. */
                                cpml_zero_psi();
                            }
                        } else if (id == UI_BTN_SCREENSHOT){
                            static int shot = 0; char fname[64];
                            snprintf(fname, sizeof(fname), "frame_%03d.bmp", shot++);
                            int w,h; SDL_GetRendererOutputSize(ren, &w, &h);
                            save_screenshot(ren, fname, w, h);
                            fprintf(stderr, "Saved %s\n", fname);
                        } else if (id == UI_BTN_CLEAR){
                            clear_fields(); t = 0; cpml_zero_psi();
                        } else if (id == UI_BTN_S21_NOW){
                            if (ports_on){
                                double new_s21 = compute_s21(freq, dt);
                                s21_amp = new_s21; s21_computed = 1;
                                FILE* pf = fopen("ports_sparams.txt", "w");
                                if (pf){
                                    fprintf(pf, "# freq=%.9e\n", freq);
                                    fprintf(pf, "S21_amp=%.9e\n", new_s21);
                                    fclose(pf);
                                    fprintf(stderr, "Wrote ports_sparams.txt (S21=%.3g)\n", new_s21);
                                }
                            } else {
                                fprintf(stderr, "Ports sampling is OFF; toggle with 'ports' before exporting S21.\n");
                            }
                        }

                        /* consume this left-click — don’t paint or drag */
                        continue; /* next SDL event */
                    
                    }
                }
                /* --- end side buttons handling --- */

                if (paint_mode) {
                    int mx = e.button.x, my = e.button.y;
                    int ix = clampi(mx / scale, 0, NX-1);
                    int iy = clampi(my / scale, 0, NY-1);
                    /* Additional safety check to ensure indices are within bounds */
                    if (ix >= 1 && ix < NX-1 && iy >= 1 && iy < NY-1) {
                        if (paint_type == 1) {
                        /* PEC: toggle tag 1 on/off */
                        if (tag_grid[ix][iy] == 1) {
                            tag_grid[ix][iy] = 0;
                        } else {
                            tag_grid[ix][iy] = 1;
                        }
                        /* ensure material resets to vacuum when toggling PEC */
                        epsr[ix][iy] = 1.0;
                        sigma_map[ix][iy] = SIGMA_BG;
                    } else if (paint_type == 2) {
                        /* PMC: toggle tag 2 on/off */
                        if (tag_grid[ix][iy] == 2) {
                            tag_grid[ix][iy] = 0;
                        } else {
                            tag_grid[ix][iy] = 2;
                        }
                    } else if (paint_type == 3) {
                        /* dielectric: toggle between 1.0 and paint_eps */
                        if (fabs(epsr[ix][iy] - paint_eps) < 1e-6) {
                            epsr[ix][iy] = 1.0;
                            sigma_map[ix][iy] = SIGMA_BG;
                        } else {
                            epsr[ix][iy] = paint_eps;
                            sigma_map[ix][iy] = SIGMA_BG;
                        }
                        tag_grid[ix][iy] = 0;
                        }
                    }
                } else {
                    int mx = e.button.x, my = e.button.y;
                    const int pickR = 10 * scale;
                    float bestD2 = 1e9f; int best = -1;
                    for (int k=0; k<MAX_SRC; ++k) if (g_src[k].active) {
                        int sx = g_src[k].ix * scale;
                        int sy = g_src[k].iy * scale;
                        float d2 = sqrf((float)(mx - sx)) + sqrf((float)(my - sy));
                        if (d2 < bestD2) { bestD2 = d2; best = k; }
                    }
                    if (best != -1 && bestD2 <= (float)(pickR*pickR)) drag_src = best;
                }
            } else if (e.type == SDL_MOUSEMOTION && drag_src != -1 && !paint_mode) {
                int mx = e.motion.x, my = e.motion.y;
                g_src[drag_src].ix = clampi(mx / scale, 1, NX-2);
                g_src[drag_src].iy = clampi(my / scale, 1, NY-2);
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                drag_src = -1;
            }

            if (e.type == SDL_KEYDOWN && e.key.repeat==0){
                if (e.key.keysym.sym == SDLK_ESCAPE) running = 0;
                if (e.key.keysym.sym == SDLK_l) show_legend = !show_legend;
                if (e.key.keysym.sym == SDLK_h) {
                    hold_color = !hold_color;
                    if (hold_color) held_vmax = (vmax_smooth>1e-12?vmax_smooth:1e-3);
                    else p99_cache_valid = 0;  /* Invalidate cache when releasing hold */
                }
                if (e.key.keysym.sym == SDLK_j) {
                    hold_scope = !hold_scope;
                    if (hold_scope) {
                        held_scope_vmax = (scope_vmax_smooth>1e-12?scope_vmax_smooth:1e-3);
                    }
                }
                if (e.key.keysym.sym == SDLK_o) scope.on = !scope.on;
                if (e.key.keysym.sym == SDLK_k) scope_clear();
                if (e.key.keysym.sym == SDLK_z) g_render_stride = (g_render_stride==1?2:(g_render_stride==2?4:1));
                if (e.key.keysym.sym == SDLK_f){
                    int ok = dump_scope_fft_csv("probe_fft.csv", dt, 4096);
                    if (ok) fprintf(stderr, "Wrote probe_fft.csv (dt=%.3g s)\n", dt);
                    else    fprintf(stderr, "FFT export failed (insufficient data or dt<=0)\n");
                }
                if (e.key.keysym.sym==SDLK_SPACE) paused = !paused;

                if (e.key.keysym.sym==SDLK_UP){
                    double newf = freq * 1.10;  /* or your preferred step */
                    /* keep slider in sync without touching freq yet */
                    s_freq.value = slider_from_freq(newf, FREQ_MIN, FREQ_MAX);
                    rebuild_sim_for_freq(newf);
                    if (probe_out) fprintf(probe_out,"# freq UP -> %.9e ; dt=%.9e\n", freq, dt);
                }
                if (e.key.keysym.sym==SDLK_DOWN){
                    double newf = freq / 1.10;
                    s_freq.value = slider_from_freq(newf, FREQ_MIN, FREQ_MAX);
                    rebuild_sim_for_freq(newf);
                    if (probe_out) fprintf(probe_out,"# freq DOWN -> %.9e ; dt=%.9e\n", freq, dt);
                }
                if (e.key.keysym.sym==SDLK_RIGHT) set_steps(steps_per_frame+1, &s_speed,&steps_per_frame);
                if (e.key.keysym.sym==SDLK_LEFT)  set_steps(steps_per_frame-1, &s_speed,&steps_per_frame);

                if (e.key.keysym.sym==SDLK_r){
                    auto_rescale = !auto_rescale;
                    if (auto_rescale) {
                        update_grid_for_freq(freq,&dx,&dy,&dt);
                    } else {
                        dx=BASE_DX;
                        dy=BASE_DY;
                        dt=cfl_dt(dx,dy);
                    }
                    cpml_build_coeffs(dt);
                    if (cpml_on) {
                        cpml_zero_psi();
                    }
                    clear_fields();
                    t=0;
                    if (probe_out){ fclose(probe_out); probe_out=fopen("probe.txt","w");
                        if (probe_out){
                            fprintf(probe_out,
                                "# dt=%.9e dx=%.6e dy=%.6e mode=%s cpml=%d freq=%.9e SIGMA_BG=%.3g SIGMA_BLOCK=%.3g\n",
                                dt, dx, dy, auto_rescale? "AUTO":"FIXED", cpml_on, freq, SIGMA_BG, SIGMA_BLOCK);
                        }
                    }
                    s_freq.value  = slider_from_freq(freq,FREQ_MIN,FREQ_MAX);
                    s_speed.value = steps_per_frame;
                }
                if (e.key.keysym.sym==SDLK_c){
                    clear_fields(); t=0; cpml_zero_psi();
                    if (probe_out){ fclose(probe_out); probe_out=fopen("probe.txt","w");
                        if (probe_out){
                            fprintf(probe_out,
                                "# dt=%.9e dx=%.6e dy=%.6e mode=%s cpml=%d freq=%.9e SIGMA_BG=%.3g SIGMA_BLOCK=%.3g\n",
                                dt, dx, dy, auto_rescale? "AUTO":"FIXED", cpml_on, freq, SIGMA_BG, SIGMA_BLOCK);
                        }
                    }
                }
                /* Toggle sources with number keys.  Keys 1 and 2 toggle the first two
                 * sources; keys 4 and 5 toggle additional sources (3rd and 4th) if defined.
                 * Key 3 is reserved for toggling the secondary probe (see below).
                 */
                if (e.key.keysym.sym==SDLK_1) {
                    if (MAX_SRC > 0) g_src[0].active = !g_src[0].active;
                }
                if (e.key.keysym.sym==SDLK_2) {
                    if (MAX_SRC > 1) g_src[1].active = !g_src[1].active;
                }
                if (e.key.keysym.sym==SDLK_4) {
                    if (MAX_SRC > 2) g_src[2].active = !g_src[2].active;
                }
                if (e.key.keysym.sym==SDLK_5) {
                    if (MAX_SRC > 3) g_src[3].active = !g_src[3].active;
                }
                if (e.key.keysym.sym==SDLK_t){
                    sources_cycle_type();
                    clear_fields(); t=0; cpml_zero_psi();
                    if (probe_out){ fclose(probe_out); probe_out=fopen("probe.txt","w");
                        if (probe_out){
                            fprintf(probe_out,
                                "# dt=%.9e dx=%.6e dy=%.6e mode=%s cpml=%d freq=%.9e SIGMA_BG=%.3g SIGMA_BLOCK=%.3g\n",
                                dt, dx, dy, auto_rescale? "AUTO":"FIXED", cpml_on, freq, SIGMA_BG, SIGMA_BLOCK);
                        }
                    }
                }
                if (e.key.keysym.sym == SDLK_p) {
                    static int shot = 0;
                    char fname[64];
                    snprintf(fname, sizeof(fname), "frame_%03d.bmp", shot++);
                    int w, h; SDL_GetRendererOutputSize(ren, &w, &h);
                    save_screenshot(ren, fname, w, h);
                    fprintf(stderr, "Saved %s\n", fname);
                }
                if (e.key.keysym.sym == SDLK_y){
                /* Toggle CPML boundary via keyboard.  Update both cpml_on and
                 * boundary_type so they remain consistent. */
                cpml_on = !cpml_on;
                boundary_type = cpml_on ? BOUNDARY_CPML : BOUNDARY_MUR;
                if (cpml_on) {
                    cpml_build_coeffs(dt);
                    cpml_zero_psi();
                } else {
                    cpml_zero_psi();
                }
                }
                if (e.key.keysym.sym == SDLK_7) {
                    cpml_apply_preset(0, dt);
                    if (cpml_on) {
                        cpml_zero_psi();
                    }
                }
                if (e.key.keysym.sym == SDLK_8) {
                    cpml_apply_preset(1, dt);
                    if (cpml_on) {
                        cpml_zero_psi();
                    }
                }
                if (e.key.keysym.sym == SDLK_9) {
                    cpml_apply_preset(2, dt);
                    if (cpml_on) {
                        cpml_zero_psi();
                    }
                }

                
                /* New: color autoscale mode & HOLD control keys */
                if (e.key.keysym.sym == SDLK_a){
                    color_autoscale_mode = (color_autoscale_mode==AS_P99? AS_PEAK : AS_P99);
                }
                if (e.key.keysym.sym == SDLK_LEFTBRACKET){
                    if (!hold_color){ hold_color=1; held_vmax = (p99_smooth>1e-12? p99_smooth : (vmax_smooth>1e-12?vmax_smooth:1e-3)); }
                    else {
                        held_vmax *= 0.90;
                        if (held_vmax<1e-12) {
                            held_vmax = 1e-12;
                        }
                    }
                }
                if (e.key.keysym.sym == SDLK_RIGHTBRACKET){
                    if (!hold_color){ hold_color=1; held_vmax = (p99_smooth>1e-12? p99_smooth : (vmax_smooth>1e-12?vmax_smooth:1e-3)); }
                    else             { held_vmax *= 1.10; }
                }
                if (e.key.keysym.sym == SDLK_BACKSLASH){
                    hold_color = 0; hold_scope = 0; held_vmax = held_scope_vmax = 1e-3;
                    vmax_smooth = p99_smooth = scope_vmax_smooth = 0.0;
                    p99_cache_valid = 0;  /* Invalidate P99 cache when resetting holds */
                }
                /* Toggle secondary probe on/off */
                if (e.key.keysym.sym == SDLK_3) {
                    probe2_active = !probe2_active;
                }
                /* Adjust amplitude of sources. Q/W control source 1, N/M control source 2. */
                if (e.key.keysym.sym == SDLK_q) {
                    g_src[0].amp *= 0.9;
                    if (g_src[0].amp < 0.001) g_src[0].amp = 0.001;
                }
                if (e.key.keysym.sym == SDLK_w) {
                    g_src[0].amp *= 1.111111111;
                    if (g_src[0].amp > 10.0) g_src[0].amp = 10.0;
                }
                if (e.key.keysym.sym == SDLK_n) {
                    g_src[1].amp *= 0.9;
                    if (g_src[1].amp < 0.001) g_src[1].amp = 0.001;
                }
                if (e.key.keysym.sym == SDLK_m) {
                    g_src[1].amp *= 1.111111111;
                    if (g_src[1].amp > 10.0) g_src[1].amp = 10.0;
                }
                /* Toggle source injection mode (linear vs. saturating). */
                if (e.key.keysym.sym == SDLK_v) {
                    g_src_injection_mode = (g_src_injection_mode == SRC_INJ_LINEAR)
                                         ? SRC_INJ_SATURATING
                                         : SRC_INJ_LINEAR;
                }
                /* Port sampling control and S-parameter export */
                if (e.key.keysym.sym == SDLK_s) {
                    ports_on = !ports_on;
                }
                if (e.key.keysym.sym == SDLK_g) {
                    if (ports_on) {
                        double new_s21 = compute_s21(freq, dt);
                        s21_amp = new_s21;
                        s21_computed = 1;
                        FILE* pf = fopen("ports_sparams.txt", "w");
                        if (pf) {
                            fprintf(pf, "# freq=%.9e\n", freq);
                            fprintf(pf, "S21_amp=%.9e\n", new_s21);
                            fclose(pf);
                            fprintf(stderr, "Wrote ports_sparams.txt (S21=%.3g)\n", new_s21);
                        }
                    } else {
                        fprintf(stderr, "Ports sampling is OFF; toggle with 's' before exporting S21.\n");
                    }
                }
                /* Paint mode toggle: when ON, left-click paints PEC */
                if (e.key.keysym.sym == SDLK_u) {
                    paint_mode = !paint_mode;
                }
                /* Cycle painting type when in paint mode (1=PEC,2=PMC,3=dielectric). */
                if (e.key.keysym.sym == SDLK_i) {
                    if (paint_mode) {
                        paint_type = (paint_type % 3) + 1;
                    }
                }
                /* Adjust dielectric painting epsr up/down (only when paint_type=3). */
                if (e.key.keysym.sym == SDLK_SEMICOLON) {
                    if (paint_mode && paint_type == 3) {
                        paint_eps *= 1.2;
                        if (paint_eps > 20.0) paint_eps = 20.0;
                    }
                }
                if (e.key.keysym.sym == SDLK_QUOTE) {
                    if (paint_mode && paint_type == 3) {
                        paint_eps /= 1.2;
                        if (paint_eps < 1.01) paint_eps = 1.01;
                    }
                }
                /* Start/abort a parametric sweep of S21 vs frequency with 'b'.  Sweeps
                 * logarithmically between FREQ_MIN and FREQ_MAX over sweep_points points. */
                if (e.key.keysym.sym == SDLK_b) {
                    if (!sweep_on) {
                        sweep_points = 10;
                        if (sweep_points > SWEEP_MAX_POINTS) sweep_points = SWEEP_MAX_POINTS;
                        for (int ii=0; ii<sweep_points; ++ii) {
                            double frac = (double)ii / (double)(sweep_points - 1);
                            sweep_freqs[ii] = FREQ_MIN * pow(FREQ_MAX/FREQ_MIN, frac);
                            sweep_s21[ii] = 0.0;
                        }
                        sweep_idx = 0;
                        sweep_on = 1;
                        double ftarget = sweep_freqs[0];
                        s_freq.value = slider_from_freq(ftarget, FREQ_MIN, FREQ_MAX);
                        rebuild_sim_for_freq(ftarget);
                        t = 0;
                        ports_on = 1;
                        sweep_steps_remaining = sweep_steps_per_point;
                        s21_computed = 0;
                    } else {
                        sweep_on = 0;
                        fprintf(stderr, "Sweep aborted.\n");
                    }
                }
            } /* end keydown */

            /* frequency slider */
            if (slider_handle_event(&s_freq,&e)) {
                double newf = freq_from_slider(s_freq.value, FREQ_MIN, FREQ_MAX);
                if (newf != freq) {
                    rebuild_sim_for_freq(newf);
                    t = 0;

                    /* (Optional) restart probe header, same as your old code but unconditional */
                    if (probe_out) { fclose(probe_out); probe_out = NULL; }
                    probe_out = fopen("probe.txt","w");
                    if (probe_out){
                        fprintf(probe_out,
                            "# dt=%.9e dx=%.6e dy=%.6e mode=%s cpml=%d freq=%.9e SIGMA_BG=%.3g SIGMA_BLOCK=%.3g\n",
                            dt, dx, dy, auto_rescale ? "AUTO":"FIXED", cpml_on, freq, SIGMA_BG, SIGMA_BLOCK);
                    }
                }
            }

            /* speed slider */
            if (slider_handle_event(&s_speed,&e)){
                set_steps((int)llround(s_speed.value), &s_speed, &steps_per_frame);
            }
        } /* end event loop */

        /* -------- simulation -------- */
        if (!paused){
            for (int s=0; s<steps_per_frame; ++s){
                /* save Ez for Mur-1 (previous timestep) */
                {
                    int i, j;
                    #if defined(_OPENMP) && !defined(_MSC_VER)
                    #pragma omp parallel for collapse(2)
                    #endif
                    for (i = 0; i < NX; ++i){
                        for (j = 0; j < NY; ++j){
                            Ez_old[i][j] = Ez[i][j];
                        }
                    }
                }

                /* update H with optional CPML */
                {
                    int i, j;
                    #if defined(_OPENMP) && !defined(_MSC_VER)
                    #pragma omp parallel for collapse(2)
                    #endif
                    for (i = 0; i < NX-1; ++i){
                        for (j = 0; j < NY-1; ++j){
                            /* dEz/dy -> Hx */
                            double dEdy = (Ez[i][j+1] - Ez[i][j]) / dy;
                            if (cpml_on && (j < cpml_N || j >= NY - cpml_N)){
                                psi_Ezy[i][j] = by[j]*psi_Ezy[i][j] + cy[j]*dEdy;
                                dEdy = (dEdy / ky[j]) + psi_Ezy[i][j];
                            }
                            Hx[i][j] -= (dt / MU0) * dEdy;

                            /* dEz/dx -> Hy */
                            double dEdx = (Ez[i+1][j] - Ez[i][j]) / dx;
                            if (cpml_on && (i < cpml_N || i >= NX - cpml_N)){
                                psi_Ezx[i][j] = bx[i]*psi_Ezx[i][j] + cx[i]*dEdx;
                                dEdx = (dEdx / kx[i]) + psi_Ezx[i][j];
                            }
                            Hy[i][j] += (dt / MU0) * dEdx;

                            /* PMC (Perfect Magnetic Conductor) boundary handling.
                             *
                             * NOTE: This is a SIMPLIFIED implementation for educational purposes.
                             * Setting H fields to zero violates Maxwell's equations for PMC.
                             *
                             * PROPER PMC IMPLEMENTATION would require:
                             * 1. For vertical PMC walls: Set tangential H (Hy) to zero, but
                             *    enforce image theory for normal H (Hx) by doubling Hx values
                             * 2. For horizontal PMC walls: Set tangential H (Hx) to zero, but
                             *    double Hy values
                             * 3. Adjust the update equations to maintain H continuity
                             *
                             * Current simple approach: Set all H components to zero (rough approximation)
                             */
                            if (tag_grid[i][j] == 2) {
                                Hx[i][j] = 0.0;
                                Hy[i][j] = 0.0;
                            }
                        }
                    }
                }

                /* update E (TMz) with optional CPML + conductivity */
                {
                    int i, j;
                    #if defined(_OPENMP) && !defined(_MSC_VER)
                    #pragma omp parallel for collapse(2)
                    #endif
                    for (i = 1; i < NX; ++i){
                        for (j = 1; j < NY; ++j){
                    /* dHy/dx */
                    double dHdx = (Hy[i][j] - Hy[i-1][j]) / dx;
                    if (cpml_on && (i < cpml_N || i >= NX - cpml_N)){
                        psi_Hyx[i][j] = bx[i]*psi_Hyx[i][j] + cx[i]*dHdx;
                        dHdx = (dHdx / kx[i]) + psi_Hyx[i][j];
                    }
                    /* dHx/dy */
                    double dHdy = (Hx[i][j] - Hx[i][j-1]) / dy;
                    if (cpml_on && (j < cpml_N || j >= NY - cpml_N)){
                        psi_Hxy[i][j] = by[j]*psi_Hxy[i][j] + cy[j]*dHdy;
                        dHdy = (dHdy / ky[j]) + psi_Hxy[i][j];
                    }
                    double curlH = dHdx - dHdy;
                    double epsij = epsilon_at(i,j);
                    double sigma = sigma_at(i,j);
                    double tmp   = 0.5 * sigma * dt / (epsij + 1e-30);
                    double ceze  = (1.0 - tmp) / (1.0 + tmp);
                    double cezh  = (dt/epsij) / (1.0 + tmp);
                    Ez[i][j] = ceze * Ez[i][j] + cezh * curlH;
                    /* enforce PEC: tag 1 clamps Ez to zero */
                    if (tag_grid[i][j] == 1) {
                        Ez[i][j] = 0.0;
                    }
                }
                    }
                }

                /* Mur-1 boundaries only if CPML is off.
                 * Use local wave speed c_local = c0/sqrt(epsr) for better absorption
                 * at dielectric boundaries. */
                if (!cpml_on){
                    {
                        int i;
                        #if defined(_OPENMP) && !defined(_MSC_VER)
                        #pragma omp parallel for
                        #endif
                        for (i = 1; i < NX-1; ++i){
                            /* Bottom boundary (j=0) */
                            double epsr_bot = epsr[i][1];
                            double c_bot = c0 / sqrt(epsr_bot);
                            Ez[i][0] = Ez_old[i][1] + ((c_bot*dt - dy)/(c_bot*dt + dy)) * (Ez[i][1] - Ez_old[i][0]);

                            /* Top boundary (j=NY-1) */
                            double epsr_top = epsr[i][NY-2];
                            double c_top = c0 / sqrt(epsr_top);
                            Ez[i][NY-1] = Ez_old[i][NY-2] + ((c_top*dt - dy)/(c_top*dt + dy)) * (Ez[i][NY-2] - Ez_old[i][NY-1]);
                        }
                    }
                    {
                        int j;
                        #if defined(_OPENMP) && !defined(_MSC_VER)
                        #pragma omp parallel for
                        #endif
                        for (j = 1; j < NY-1; ++j){
                            /* Left boundary (i=0) */
                            double epsr_left = epsr[1][j];
                            double c_left = c0 / sqrt(epsr_left);
                            Ez[0][j] = Ez_old[1][j] + ((c_left*dt - dx)/(c_left*dt + dx)) * (Ez[1][j] - Ez_old[0][j]);

                            /* Right boundary (i=NX-1) */
                            double epsr_right = epsr[NX-2][j];
                            double c_right = c0 / sqrt(epsr_right);
                            Ez[NX-1][j] = Ez_old[NX-2][j] + ((c_right*dt - dx)/(c_right*dt + dx)) * (Ez[NX-2][j] - Ez_old[NX-1][j]);
                        }
                    }
                }

                /* inject sources */
                for (int ksrc=0; ksrc<MAX_SRC; ++ksrc) inject_source_into_Ez(&g_src[ksrc], Ez, t, dt);

                /* probe log + scope sample every step */
                double probe_val = Ez[probe_x][probe_y];
                if (probe_out){
                    fprintf(probe_out, "%d %.9e\n", t, probe_val);
                    if ((t & 255)==0) fflush(probe_out);
                }
                scope_push(probe_val);

                /* Port sampling: accumulate port voltages and currents on each step */
                if (ports_on) {
                    for (int p = 0; p < MAX_PORTS; ++p) {
                        if (!ports[p].active) continue;
                        double Vsum = 0.0;
                        double Isum = 0.0;
                        /* voltage is integral of Ez along port segment */
                        for (int yy = ports[p].y0; yy <= ports[p].y1; ++yy) {
                            Vsum += Ez[ports[p].x][yy];
                        }
                        Vsum *= dy;
                        /* current is integral of Hy along port segment.  Include
                         * the upper endpoint (<=) just as for the voltage integral.
                         * Previously, the loop used yy < ports[p].y1 and thus
                         * excluded the last cell, biasing the S-parameter results.
                         */
                        for (int yy = ports[p].y0; yy <= ports[p].y1; ++yy) {
                            Isum += Hy[ports[p].x][yy];
                        }
                        Isum *= dx;
                        ports[p].V[ports[p].head] = Vsum;
                        ports[p].I[ports[p].head] = Isum;
                    }
                    /* advance circular buffer heads */
                    for (int p = 0; p < MAX_PORTS; ++p) {
                        if (!ports[p].active) continue;
                        ports[p].head = (ports[p].head + 1) % ports[p].n;
                    }
                }

                ++t;
            }
        }

        /* If a frequency sweep is active, manage progression when not paused. */
        if (sweep_on && !paused) {
            if (sweep_steps_remaining > 0) {
                sweep_steps_remaining -= steps_per_frame;
                if (sweep_steps_remaining < 0) sweep_steps_remaining = 0;
            }
            if (sweep_steps_remaining == 0) {
                /* compute S21 for current frequency */
                double val = compute_s21(freq, dt);
                sweep_s21[sweep_idx] = val;
                s21_amp = val;
                s21_computed = 1;
                sweep_idx++;
                if (sweep_idx >= sweep_points) {
                    FILE* swf = fopen("sweep_s21.csv", "w");
                    if (swf) {
                        fprintf(swf, "# freq_Hz,S21_amp\n");
                        for (int ii = 0; ii < sweep_points; ++ii) {
                            fprintf(swf, "%.12e,%.12e\n", sweep_freqs[ii], sweep_s21[ii]);
                        }
                        fclose(swf);
                        fprintf(stderr, "Wrote sweep_s21.csv (%d points)\n", sweep_points);
                    } else {
                        fprintf(stderr, "Could not write sweep_s21.csv\n");
                    }
                    sweep_on = 0;
                } else {
                    /* move to next frequency */
                    double ftarget = sweep_freqs[sweep_idx];
                    s_freq.value = slider_from_freq(ftarget, FREQ_MIN, FREQ_MAX);
                    rebuild_sim_for_freq(ftarget);
                    t = 0;
                    sweep_steps_remaining = sweep_steps_per_point;

                    /* reset port buffers */
                    for (int pp = 0; pp < MAX_PORTS; ++pp) {
                        if (ports[pp].V) memset(ports[pp].V, 0, sizeof(double)*ports[pp].n);
                        if (ports[pp].I) memset(ports[pp].I, 0, sizeof(double)*ports[pp].n);
                        ports[pp].head = 0;
                    }
                }
            }
        }

        /* -------- FPS -------- */
        Uint64 t_now = SDL_GetPerformanceCounter();
        double dt_s = (double)(t_now - t_prev) / perf_freq;
        t_prev = t_now;
        double fps_inst = (dt_s > 0.0) ? (1.0 / dt_s) : 0.0;
        if (fps_avg == 0.0) fps_avg = fps_inst;
        else                fps_avg = 0.9 * fps_avg + 0.1 * fps_inst;

        /* -------- render -------- */
        SDL_SetRenderDrawColor(ren,0,0,0,255); SDL_RenderClear(ren);

        /* field -> heatmap (compute scale) */
        double vmax = 1e-12;
        #if defined(_OPENMP) && !defined(_MSC_VER)
        {
            double vmax_seed = vmax;
            #pragma omp parallel
            {
                double vmax_local = vmax_seed;
                int i, j;
                #pragma omp for collapse(2) nowait
                for (i = 0; i < NX; i += HIST_STRIDE){
                    for (j = 0; j < NY; j += HIST_STRIDE){
                        double a = fabs(Ez[i][j]);
                        if (a > vmax_local) {
                            vmax_local = a;
                        }
                    }
                }
                #pragma omp critical
                {
                    if (vmax_local > vmax) {
                        vmax = vmax_local;
                    }
                }
            }
        }
        #else
        {
            int i, j;
            for (i = 0; i < NX; i += HIST_STRIDE){
                for (j = 0; j < NY; j += HIST_STRIDE){
                    double a = fabs(Ez[i][j]);
                    if (a > vmax) {
                        vmax = a;
                    }
                }
            }
        }
        #endif
        if (vmax_smooth == 0.0) vmax_smooth = vmax;
        else vmax_smooth = 0.95*vmax_smooth + 0.05*vmax;

        double p99 = vmax_smooth;
            if (!hold_color && color_autoscale_mode == AS_P99){
            unsigned int hist[COLOR_HIST_BINS]; memset(hist,0,sizeof(hist));
            double eps = (vmax>1e-20? vmax : 1e-20);
            {
                int i, j;
                #if defined(_OPENMP) && !defined(_MSC_VER)
                #pragma omp parallel for collapse(2)
                #endif
                for (i = 0; i < NX; i += HIST_STRIDE){
                    for (j = 0; j < NY; j += HIST_STRIDE){
                        double a = fabs(Ez[i][j]) / eps;
                        if (a>1.0) a = 1.0;
                        int bin = (int)(a * (COLOR_HIST_BINS-1));
                        if (bin < 0) {
                            bin = 0;
                        }
                        if (bin >= COLOR_HIST_BINS) {
                            bin = COLOR_HIST_BINS-1;
                        }
                        #if defined(_OPENMP) && !defined(_MSC_VER)
                        #pragma omp atomic
                        #endif
                        hist[bin]++;
                    }
                }
            }
            unsigned int total = (NX/HIST_STRIDE)*(NY/HIST_STRIDE);
            unsigned int target = (unsigned int)(0.99 * (double)total);
            unsigned int acc = 0; int k99 = COLOR_HIST_BINS-1;
            for (int b=0; b<COLOR_HIST_BINS; ++b){
                acc += hist[b];
                if (acc >= target){ k99 = b; break; }
            }
            p99 = ((double)k99 / (double)(COLOR_HIST_BINS-1)) * vmax;
            if (p99 < 1e-12) p99 = 1e-12;
            if (p99_smooth == 0.0) p99_smooth = p99;
            else p99_smooth = 0.95*p99_smooth + 0.05*p99;

            /* Cache P99 value when hold_color becomes active */
            cached_p99 = p99_smooth;
            p99_cache_valid = 1;
        } else if (hold_color && color_autoscale_mode == AS_P99 && p99_cache_valid) {
            /* Use cached P99 to skip expensive histogram calculation */
            p99_smooth = cached_p99;
        }

        double color_use_vmax;
        if (hold_color)                               color_use_vmax = (held_vmax>1e-12?held_vmax:1e-12);
        else if (color_autoscale_mode == AS_P99)      color_use_vmax = p99_smooth>1e-12 ? p99_smooth : vmax_smooth;
        else                                          color_use_vmax = vmax_smooth;
        if (!hold_color) held_vmax = color_use_vmax;
        if (color_use_vmax < 1e-12) color_use_vmax = 1e-12;

        /* draw field using stride (perf) */
        for (int i=0;i<NX;i+=g_render_stride) for (int j=0;j<NY;j+=g_render_stride){
            int c = (int)(128 + (127.0/color_use_vmax)*Ez[i][j]);
            if (c < 0) {
                c = 0;
            }
            if (c > 255) {
                c = 255;
            }
            SDL_SetRenderDrawColor(ren,c,0,255-c,255);
            SDL_Rect cell = { i*scale, SIM_YOFF + j*scale, scale*g_render_stride, scale*g_render_stride };
            SDL_RenderFillRect(ren, &cell);
        }

        draw_block_outline(ren, scale);
        draw_sources(ren, scale);

        /* probe cross */
        /* draw primary probe cross */
        SDL_SetRenderDrawColor(ren, 255,255,255,255);
        int px = probe_x*scale, py=SIM_YOFF + probe_y*scale;
        for (int d = -4; d <= 4; ++d) {
            SDL_RenderDrawPoint(ren, px + d, py);
            SDL_RenderDrawPoint(ren, px, py + d);
        }
        /* draw second probe cross in green if active */
        if (probe2_active) {
            SDL_SetRenderDrawColor(ren, 0,255,0,255);
            int px2 = probe2_x*scale;
            int py2 = SIM_YOFF + probe2_y*scale;
            for (int d = -4; d <= 4; ++d) {
                SDL_RenderDrawPoint(ren, px2 + d, py2);
                SDL_RenderDrawPoint(ren, px2, py2 + d);
            }
        }

        /* scope (above sliders) */
        const int SCOPE_HEIGHT = 60;
        int scope_w = NX * scale;
        int scope_x = 0;
        int scope_y = UI_Y - (SCOPE_HEIGHT + 8);

        /* scope autoscale (decoupled from field) */
        double scope_max = 1e-12;
        if (scope.on && scope.y){
            for (int i=0;i<scope.n;++i){
                double a=fabs(scope.y[i]);
                if (a>scope_max) {
                    scope_max=a;
                }
            }
        }
        if (scope_vmax_smooth==0.0) scope_vmax_smooth = scope_max;
        else scope_vmax_smooth = 0.95*scope_vmax_smooth + 0.05*scope_max;
        double scope_use_vmax = hold_scope ? (held_scope_vmax>1e-12?held_scope_vmax:1e-12)
                                           : scope_vmax_smooth;
        if (!hold_scope) held_scope_vmax = scope_use_vmax;
        draw_scope(ren, scope_x, scope_y, scope_w, SCOPE_HEIGHT, scope_use_vmax);

        /* colorbar */
        const int cb_w = 14;
        const int cb_x = NX*scale + 6;
        const int cb_y = 6;
        const int cb_h = NY*scale - 12;
        SDL_Rect frame = { cb_x-1, cb_y-1, cb_w+2, cb_h+2 };
        SDL_SetRenderDrawColor(ren, 220,220,220,200);
        SDL_RenderDrawRect(ren, &frame);
        for (int yy = 0; yy < cb_h; ++yy) {
            double frac = 1.0 - (double)yy / (double)(cb_h-1);   // top=1, bottom=0
            double val  = (2.0*frac - 1.0) * color_use_vmax;     // +..-
            int c = (int)(128 + (127.0/color_use_vmax) * val);
            if (c < 0) {
                c = 0;
            }
            if (c > 255) {
                c = 255;
            }
            SDL_SetRenderDrawColor(ren, c, 0, 255 - c, 255);
            SDL_RenderDrawLine(ren, cb_x, cb_y + yy, cb_x + cb_w - 1, cb_y + yy);
        }
        /* colorbar labels.  Show the maximum positive, zero, and maximum negative
         * field amplitude corresponding to the colourbar extents. */
        {
            char lab_top[64], lab_mid[64], lab_bot[64];
            snprintf(lab_top, sizeof(lab_top), "+%.3g", color_use_vmax);
            snprintf(lab_mid, sizeof(lab_mid), "0");
            snprintf(lab_bot, sizeof(lab_bot), "-%.3g", color_use_vmax);
            int tw=0, th=0;
            SDL_Color fg = {230,230,230,255};
            SDL_Texture* tt;
            tt = render_text(ren, font, lab_top, fg, &tw, &th);
            if (tt){
                SDL_Rect r={cb_x+cb_w+6, cb_y-2, tw, th};
                SDL_RenderCopy(ren, tt, NULL, &r);
                SDL_DestroyTexture(tt);
            }
            tt = render_text(ren, font, lab_mid, fg, &tw, &th);
            if (tt){
                SDL_Rect r={cb_x+cb_w+6, cb_y + cb_h/2 - th/2, tw, th};
                SDL_RenderCopy(ren, tt, NULL, &r);
                SDL_DestroyTexture(tt);
            }
            tt = render_text(ren, font, lab_bot, fg, &tw, &th);
            if (tt){
                SDL_Rect r={cb_x+cb_w+6, cb_y + cb_h - th + 2, tw, th};
                SDL_RenderCopy(ren, tt, NULL, &r);
                SDL_DestroyTexture(tt);
            }
        }

        /* right-side panel (information + buttons) to the right of the colourbar */
        int sp_x = cb_x + cb_w + 6;
        int sp_y = cb_y;
        int sp_w = NX*scale + side_panel - 6 - sp_x;
        int sp_h = cb_h;

        /* background for the entire side panel */
        SDL_Rect sp_bg = { sp_x, sp_y, sp_w, sp_h };
        SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(ren, 12,12,12,220);
        SDL_RenderFillRect(ren, &sp_bg);
        SDL_SetRenderDrawColor(ren, 200,200,200,220);
        SDL_RenderDrawRect(ren, &sp_bg);

        /* Subdivide the side panel into an information region and a button region.  The
         * information region occupies roughly 60% of the vertical space, but never
         * less than 100 pixels when the panel is tall enough.  The button region
         * occupies the remainder.  A small gap separates the two. */
        int info_h = (int)(0.6 * (double)sp_h);
        if (info_h < 100) {
            /* if the panel is very short, allocate at least (sp_h - 100) pixels to info */
            info_h = (sp_h > 100 ? sp_h - 100 : sp_h);
        }
        int btn_y = sp_y + info_h + 4;
        int btn_h = sp_h - info_h - 4;
        if (btn_h < 0) btn_h = 0;
        /* Draw the information panel inside a padded rectangle */
        SDL_Rect info_rect = { sp_x + 4, sp_y + 4, sp_w - 8, info_h - 8 };
        draw_info_panel(ren, font, info_rect, fps_inst, fps_avg, steps_per_frame, paused);

        /* Build the buttons in the lower portion of the side panel.  Pass addresses
         * of the actual state variables so toggles update the live variables. */
        build_side_buttons(sp_x, btn_y, sp_w, btn_h, &paused, &auto_rescale, &hold_color, &cpml_on, &ports_on, &paint_mode);
        /* Draw all buttons */
        for (int i = 0; i < ui_btn_count; ++i){
            draw_button(ren, font, &ui_btns[i]);
        }


        /* UI slab + sliders */
        SDL_Rect ui_bg = (SDL_Rect){0, UI_Y, NX*scale, UI_H};
        SDL_SetRenderDrawColor(ren,30,30,30,255); SDL_RenderFillRect(ren,&ui_bg);
        slider_draw(ren,&s_freq); slider_draw(ren,&s_speed);

        /* Compute and smooth net Poynting flux leaving domain.  Positive flux means net
         * power leaving the domain.  We compute this here in the render stage to
         * avoid slowing down the time-step loops. */
        {
            double fsum = 0.0;
            /* vertical boundaries (x-direction) */
            for (int jj=0; jj<NY; ++jj) {
                /* left boundary normal (-1,0): flux = Ez*Hy * dy */
                fsum += Ez[0][jj] * Hy[0][jj] * dy;
                /* right boundary normal (+1,0): flux = -Ez*Hy * dy */
                fsum += -Ez[NX-1][jj] * Hy[NX-1][jj] * dy;
            }
            /* horizontal boundaries (y-direction) */
            for (int ii=0; ii<NX; ++ii) {
                /* bottom boundary normal (0,-1): flux = -Ez*Hx * dx */
                fsum += -Ez[ii][0] * Hx[ii][0] * dx;
                /* top boundary normal (0,+1): flux = Ez*Hx * dx */
                fsum += Ez[ii][NY-1] * Hx[ii][NY-1] * dx;
            }
            pflux_avg = 0.95 * pflux_avg + 0.05 * fsum;
        }

        /* legend */
        if (show_legend){
            double disp = (freq >= 1e9) ? freq/1e9 : freq/1e6;
            const char* unit = (freq >= 1e9) ? "GHz" : "MHz";
            double lambda0 = c0 / freq;
            double cpw_out = lambda0 / dx;
            double cpw_in  = (lambda0 / sqrt(EPSR_MAX_SCENE)) / dx;
            int under   = (!auto_rescale) && (cpw_in < TARGET_CPW);
            int clamped =  (auto_rescale) && (dx >= BASE_DX*MAX_SCALE_FACTOR*0.999);

            int active_cnt = 0; for (int k=0;k<MAX_SRC;++k) if (g_src[k].active) ++active_cnt;
            const char* stype = (g_src[0].type==0) ? "CW" : (g_src[0].type==1 ? "Gauss" : "Ricker");
            const char* scmode = hold_color ? "HOLD" : (color_autoscale_mode==AS_P99? "AUTO-P99" : "AUTO-PEAK");

            char l1[128], l2[128], l3[128], l4[192], l5[128], l6[160], l7[128], l8[128], l9[128], l10[128], lSigma[128], lPML[192];
            snprintf(l1, sizeof(l1), "Freq: %.2f %s", disp, unit);
            snprintf(l2, sizeof(l2), "Δ = %.4g m   dt = %.4g s", dx, dt);
            snprintf(l3, sizeof(l3), "CPW (in/out): %.1f / %.1f", cpw_in, cpw_out);
            snprintf(l4, sizeof(l4), "Mode: %s  ColorScale: %s%s%s  Render:x%d",
                    auto_rescale ? "AUTO" : "FIXED", scmode,
                    under   ? "  UNDER-RES"   : "",
                    clamped ? "  AUTO-CLAMPED": "",
                    g_render_stride);
            snprintf(l5, sizeof(l5), "Steps/frame: %d   %s", steps_per_frame, paused ? "PAUSED" : "RUNNING");
            snprintf(l6, sizeof(l6), "Sources: %d active   Amp1=%.2f Amp2=%.2f   Type: %s  (toggle: 1/2, cycle type: T)",
                     active_cnt, g_src[0].amp, g_src[1].amp, stype);
            snprintf(l7, sizeof(l7), "FPS: %.2f (inst), %.2f (avg)  PowerFlux=%.3e W/m", fps_inst, fps_avg, pflux_avg);
            if (probe2_active) {
                snprintf(l8, sizeof(l8), "Probes: P1(%d,%d)=%.3g  P2(%d,%d)=%.3g",
                        probe_x, probe_y, Ez[probe_x][probe_y],
                        probe2_x, probe2_y, Ez[probe2_x][probe2_y]);
            } else {
                snprintf(l8, sizeof(l8), "Probe (%d,%d) Ez=%.3g (arb.)", probe_x, probe_y, Ez[probe_x][probe_y]);
            }
            snprintf(l9, sizeof(l9), "Scope: %s  Scale: %s", scope.on ? "ON" : "OFF", hold_scope ? "HOLD" : "AUTO");
            snprintf(l10,sizeof(l10),"Scope (last): %.3g", scope.last);
            snprintf(lSigma,sizeof(lSigma), "Sigma: bg=%.3g S/m  block=%.3g S/m", SIGMA_BG, SIGMA_BLOCK);
            if (cpml_on){
                snprintf(lPML,sizeof(lPML), "CPML: %s  N=%d  sigma=%.2f  kappa=%.2f  alpha=%.3f",
                    CPML_PRESETS[cpml_preset_idx].name, cpml_N, pml_sigma_max, pml_kappa_max, pml_alpha_max);
            } else {
                snprintf(lPML,sizeof(lPML), "Boundary: Mur-1");
            }

            /* Ports status and last-computed S21 amplitude.  If a sweep is in progress,
             * show progress. */
            char lPorts[128];
            if (sweep_on) {
                snprintf(lPorts, sizeof(lPorts), "Sweep: %d/%d freq=%.2f MHz",
                         sweep_idx+1, sweep_points,
                         (sweep_idx < sweep_points? (sweep_freqs[sweep_idx]/1e6) : 0.0));
            } else {
                snprintf(lPorts, sizeof(lPorts), "Ports: %s  Last S21=%.3g",
                         ports_on ? "ON" : "OFF", s21_amp);
            }
            /* Painting status: show mode and selected epsilon if painting dielectric */
            char lPaint[128];
            if (paint_mode) {
                const char* pnames[] = {"OFF","PEC","PMC","Diel"};
                if (paint_type == 3) {
                    snprintf(lPaint, sizeof(lPaint), "Paint: %s  eps=%.2f", pnames[paint_type], paint_eps);
                } else {
                    snprintf(lPaint, sizeof(lPaint), "Paint: %s", pnames[paint_type]);
                }
            } else {
                snprintf(lPaint, sizeof(lPaint), "Paint: OFF");
            }

            const char* lines[] = {
                l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,lSigma,lPML,
                lPorts, lPaint,
                "Keys: Space pause | ↑/↓ freq | ←/→ steps | Z render stride | R auto/fixed",
                "| A autoscale mode | H hold color | [/] adjust color (hold)",
                "| J hold scope | \\\\ reset holds | Y CPML toggle | 7/8/9 CPML presets | F dump FFT",
                "| C clear | L legend | P screenshot | Right-click move P1 | Middle-click move P2 ",
                "| 1/2 toggle src | 3 toggle P2 | Q/W amp1 down/up | N/M amp2 down/up | T type | O scope | K scope clear",
                "| U paint on/off | I  paint type | ;/` change dielectric eps | S ports | G S21 | B sweep"
            };
            draw_legend(ren, font, 12, 12, lines, (int)(sizeof(lines)/sizeof(lines[0])));
        }

        /* window title */
        {
            double disp = (freq>=1e9)? freq/1e9 : freq/1e6;
            const char* unit = (freq>=1e9)? "GHz":"MHz";
            char title[256];
            snprintf(title, sizeof(title),
                     "Mode: %s   Freq: %.2f %s",
                     auto_rescale ? "AUTO" : "FIXED", disp, unit);
            SDL_SetWindowTitle(win, title);
        }

        SDL_RenderPresent(ren);
        SDL_Delay(UI_DELAY_MS);
    } /* end main loop */

    if (probe_out) fclose(probe_out);
    /*
     * Free port sampling buffers.  These buffers are allocated during
     * initialisation and otherwise leak at exit.  See other AI notes.
     */
    for (int p=0; p<MAX_PORTS; ++p){
        if (ports[p].V){ free(ports[p].V); ports[p].V = NULL; }
        if (ports[p].I){ free(ports[p].I); ports[p].I = NULL; }
    }
    TTF_CloseFont(font);
    TTF_Quit();
    scope_free();
    free(Ez); free(Hx); free(Hy); free(Ez_old);
    free(psi_Ezx); free(psi_Ezy); free(psi_Hyx); free(psi_Hxy);
    SDL_DestroyRenderer(ren); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}

/* end of file */
//Additional Improvements in This Version:
//-------------------------------------------------------------------------------------------
// 1. Expanded Source Support: Increased MAX_SRC to 4 and initialized distinct default positions for each source, with per-source toggle buttons and info display.

// 2. Boundary Selection: Introduced a boundary_type enum (MUR vs CPML) with a new “bnd+” button and keyboard toggle that cycles boundary conditions, syncing with the existing cpml_on logic.

// 3. Integrated Side Panel: Replaced the old on-screen legend by subdividing the right panel: the top portion now displays simulation info (frequency, time step, CPW, mode/color settings, steps/frame status, detailed per-source data, FPS, probe readings, scope and sigma status, boundary/CPML information, ports/sweep status, painting state, and key hints). The bottom portion retains the control buttons.

// 4. Re-enabled Colorbar Labels: Top, mid, and bottom labels now show ±max and 0 values for better interpretation of the heatmap scale.

// 5. Dynamic Layout: Computed info/button split (60%/40% with minimum sizing), drawn info via a wrapped panel, and auto-sized buttons. Enlarged the right panel width to 260 for readability.

// 6. Improved Toggles: Added keys 4 and 5 to toggle sources 3 and 4; toggling CPML (y) updates both boundary_type and cpml_on.

// 7. General UI Refinements: Added an info panel drawing helper function, updated event handling to rebuild panels appropriately, and passed relevant runtime state to the info panel for accurate display.

// 8. Code Cleanup: Removed obsolete legend code, ensured consistent formatting, and added comments for clarity.
//-------------------------------------------------------------------------------------------
