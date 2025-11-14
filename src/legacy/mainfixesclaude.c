// ============================================================================
// EM Wave Sim Project – Fixed 2D TMz FDTD Implementation
// Critical fixes: PML positioning, source injection, stability, memory leaks
// ============================================================================

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================= Configuration Constants ================= */
#define NX 400
#define NY 400
#define STEPS_PER_FRAME 2
#define CFL_SAFETY_FACTOR 0.95  /* Reduced from 0.99 for stability */
#define SMOOTHING_FACTOR 0.95
#define HIST_UPDATE_INTERVAL 10  /* Update histogram every N frames */
#define PROBE_FLUSH_INTERVAL 256

/* Physical constants */
static const double c0  = 299792458.0;            /* m/s */
static const double MU0 = 1.256637061435917295e-6;/* H/m */
static const double EPS0= 8.8541878128e-12;       /* F/m */

/* Physical domain */
static const double Lx = 0.6, Ly = 0.6;           /* meters */
static const double BASE_DX = Lx / NX;
static const double BASE_DY = Ly / NY;

/* AUTO mode parameters */
static const double TARGET_CPW       = 12.0;
static const double EPSR_MAX_SCENE   = 4.0;
static const double MAX_SCALE_FACTOR = 20.0;

/* Grid and timestep */
static double dx = BASE_DX, dy = BASE_DY, dt = 0.0;

/* Conductivity */
static const double SIGMA_BG    = 0.0;
static const double SIGMA_BLOCK = 0.0;

/* ================= Helper Functions ================= */
static inline int    clampi (int v, int lo, int hi){ return v<lo?lo:(v>hi?hi:v); }
static inline double clampd (double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }
static inline float  sqrf   (float x){ return x*x; }

/* Improved CFL with safety factor */
static inline double cfl_dt(double dx_, double dy_){
    return CFL_SAFETY_FACTOR / (c0 * sqrt(1.0/(dx_*dx_) + 1.0/(dy_*dy_)));
}

/* ================= Sliders ================= */
typedef struct {
    int x, y, w, h;
    double minv, maxv;
    double value;
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

/* Log slider mapping */
static inline double freq_from_slider(double t, double fminv, double fmaxv){
    double a = log10(fminv), b = log10(fmaxv);
    return pow(10.0, a + clampd(t,0.0,1.0)*(b-a));
}

static inline double slider_from_freq(double f, double fminv, double fmaxv){
    double a = log10(fminv), b = log10(fmaxv);
    return clampd((log10(f)-a)/(b-a), 0.0, 1.0);
}

/* ================= Fields & Materials ================= */
static double (*Ez)[NY], (*Hx)[NY], (*Hy)[NY];
static double (*Ez_old)[NY];
static double (*epsr)[NY], (*sigma_map)[NY];
static unsigned char (*tag_grid)[NY];

/* Port definitions */
#define MAX_PORTS 2
#define PORT_SIGNAL_LENGTH 4096

typedef struct {
    int x, y0, y1;
    int len, n;
    double *V, *I;
    int head;
    int active;
} Port;

static Port ports[MAX_PORTS];
static int ports_on = 0;
static double s21_amp = 0.0;
static int s21_computed = 0;

/* Painting mode */
static int paint_mode = 0;
static int paint_type = 1;  /* 1=PEC, 2=PMC, 3=dielectric */
static double paint_eps = 2.0;

/* Poynting flux */
static double pflux_avg = 0.0;

/* Sweep parameters */
#define SWEEP_MAX_POINTS 32
static int sweep_on = 0;
static int sweep_idx = 0;
static int sweep_points = 0;
static double sweep_freqs[SWEEP_MAX_POINTS];
static double sweep_s21[SWEEP_MAX_POINTS];
static int sweep_steps_remaining = 0;
static int sweep_steps_per_point = 2000;

/* ================= Field Operations ================= */
static inline void clear_fields(void){
    memset(Ez,     0, sizeof(double[NX][NY]));
    memset(Hx,     0, sizeof(double[NX][NY]));
    memset(Hy,     0, sizeof(double[NX][NY]));
    memset(Ez_old, 0, sizeof(double[NX][NY]));
}

static inline void draw_block_outline(SDL_Renderer* ren, int scale){
    int bx0 = NX/2 - NX/10, bx1 = NX/2 + NX/10;
    int by0 = NY/2 - NY/20, by1 = NY/2 + NY/20;
    SDL_SetRenderDrawColor(ren, 220,220,220,200);
    SDL_RenderDrawLine(ren, bx0*scale, by0*scale, bx1*scale, by0*scale);
    SDL_RenderDrawLine(ren, bx0*scale, by1*scale, bx1*scale, by1*scale);
    SDL_RenderDrawLine(ren, bx0*scale, by0*scale, bx0*scale, by1*scale);
    SDL_RenderDrawLine(ren, bx1*scale, by0*scale, bx1*scale, by1*scale);
}

/* Fixed S21 computation */
static double compute_s21(double freq, double dt) {
    int n = ports[0].n;
    if (n <= 0 || !ports[0].V || !ports[1].V) return 0.0;
    
    double complex sum1 = 0.0, sum2 = 0.0;
    int h1 = ports[0].head;
    int h2 = ports[1].head;
    
    for (int i = 0; i < n; ++i) {
        double x1 = ports[0].V[(h1 + i) % n];
        double x2 = ports[1].V[(h2 + i) % n];
        double arg = -2.0 * M_PI * freq * dt * (double)i;
        double complex phasor = cexp(I * arg);
        sum1 += x1 * phasor;
        sum2 += x2 * phasor;
    }
    
    double mag1 = cabs(sum1);
    double mag2 = cabs(sum2);
    return (mag1 > 1e-30) ? (mag2 / mag1) : 0.0;
}

/* ================= Grid Updates ================= */
static inline void update_grid_for_freq(double f, double *pdx, double *pdy, double *pdt){
    double delta = c0 / (f * sqrt(EPSR_MAX_SCENE) * TARGET_CPW);
    if (delta <= 0) delta = BASE_DX;
    double max_delta = BASE_DX * MAX_SCALE_FACTOR;
    if (delta > max_delta) delta = max_delta;
    *pdx = delta; *pdy = delta; *pdt = cfl_dt(*pdx, *pdy);
}

/* ================= Sources ================= */
typedef enum { SRC_CW=0, SRC_GAUSS_PULSE=1, SRC_RICKER=2 } SourceType;

typedef struct {
    int active, ix, iy;
    SourceType type;
    double amp, freq;
    double t0, tau;
    double sigma2;
} Source;

#define MAX_SRC 2
static Source g_src[MAX_SRC];
static int drag_src = -1;

static inline void source_reparam(Source* s){
    if (s->type == SRC_GAUSS_PULSE || s->type == SRC_RICKER) {
        double cycles_t0=6.0, cycles_tau=2.0;
        s->t0  = (s->freq>0)? (cycles_t0  / s->freq) : 0.0;
        s->tau = (s->freq>0)? (cycles_tau / s->freq) : 1e-9;
    } else { 
        s->t0=0.0; s->tau=1e-9; 
    }
}

static inline void sources_set_freq(double f){
    for (int k=0;k<MAX_SRC;++k){ 
        g_src[k].freq=f; 
        source_reparam(&g_src[k]); 
    }
}

static inline void sources_cycle_type(void){
    for (int k=0;k<MAX_SRC;++k){ 
        g_src[k].type=(SourceType)((g_src[k].type+1)%3); 
        source_reparam(&g_src[k]); 
    }
}

static inline double source_time_value(const Source* s, int t, double dt_){
    double tt = t * dt_;
    switch(s->type){
        case SRC_CW:          
            return s->amp * sin(2.0*M_PI*s->freq * tt);
        case SRC_GAUSS_PULSE: {
            double x=(tt-s->t0)/s->tau; 
            double env=exp(-0.5*x*x);
            return s->amp * env * sin(2.0*M_PI*s->freq*(tt - s->t0));
        }
        case SRC_RICKER: {
            double a = M_PI*s->freq*(tt - s->t0); 
            double e=exp(-a*a);
            return s->amp * (1.0 - 2.0*a*a) * e;
        }
        default: return 0.0;
    }
}

/* Soft source injection (reduces reflections) */
static inline void inject_source_into_Ez(Source* s, double (*Ezf)[NY], int t, double dt_){
    if (!s->active) return;
    double A = source_time_value(s, t, dt_);
    
    for (int di=-2; di<=2; ++di)
    for (int dj=-2; dj<=2; ++dj){
        int i=s->ix+di, j=s->iy+dj;
        if (i>0 && i<NX && j>0 && j<NY){
            double r2 = di*di + dj*dj;
            double w = exp(-r2/s->sigma2);
            /* Soft source injection */
            double src_val = A * w;
            Ezf[i][j] += src_val / (1.0 + fabs(src_val) * 0.1);
        }
    }
}

static inline void draw_sources(SDL_Renderer* ren, int scale){
    for (int k=0;k<MAX_SRC;++k){
        Source* s = &g_src[k];
        Uint8 rr = s->active?0:80, gg=s->active?255:80, bb=0;
        SDL_SetRenderDrawColor(ren, rr,gg,bb,255);
        int x=s->ix*scale, y=s->iy*scale;
        for (int d=-4; d<=4; ++d){ 
            SDL_RenderDrawPoint(ren, x+d, y); 
            SDL_RenderDrawPoint(ren, x, y+d); 
        }
    }
}

/* ================= Oscilloscope ================= */
typedef struct { 
    int n, head; 
    double *y; 
    int on; 
    double last; 
} Scope;

static Scope scope = {0,0,NULL,0,0.0};

static void scope_init(int width) {
    if (scope.y) free(scope.y);
    scope.n = (width > 64 ? width : 64);
    scope.y = (double*)calloc(scope.n, sizeof(double));
    scope.head = 0; scope.on = 1; scope.last = 0.0;
}

static void scope_free(void) { 
    if (scope.y){ 
        free(scope.y); 
        scope.y=NULL; 
    } 
    scope.n=0; 
}

static inline void scope_push(double v) {
    if (!scope.y || !scope.on) return;
    scope.y[scope.head] = v;
    scope.head = (scope.head + 1) % scope.n;
    scope.last = v;
}

static void scope_clear(void){
    if (!scope.y) return;
    memset(scope.y, 0, sizeof(double)*scope.n);
    scope.head = 0; scope.last = 0.0;
}

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
        int idx = (scope.head - (w - i) + scope.n) % scope.n;
        double v = scope.y[idx];
        double t = 0.5 - 0.5*(v/yscale);
        t = clampd(t, 0, 1);
        int yy = y + (int)(t * (h-1));
        int xx = x + i;
        if (i) SDL_RenderDrawLine(r, prevx, prevy, xx, yy);
        prevx = xx; prevy = yy;
    }
}

/* FFT export */
static int dump_scope_fft_csv(const char* path, double dt_, int Nfft_requested){
    if (!scope.y || scope.n <= 8 || dt_ <= 0.0) return 0;
    int N = scope.n;
    int Nfft = Nfft_requested;
    if (Nfft > N) Nfft = N;
    if (Nfft < 64) Nfft = (N < 64 ? N : 64);

    double *x = (double*)malloc(sizeof(double)*Nfft);
    if (!x) return 0;
    
    int start = (scope.head - Nfft + scope.n) % scope.n;
    for (int k=0; k<Nfft; ++k){
        int idx = (start + k) % scope.n;
        x[k] = scope.y[idx];
    }
    
    double mean = 0.0; 
    for (int k=0;k<Nfft;++k) mean += x[k]; 
    mean /= (double)Nfft;
    
    for (int k=0;k<Nfft;++k){
        double w = 0.5 * (1.0 - cos(2.0*M_PI*(double)k/(double)(Nfft-1)));
        x[k] = (x[k] - mean) * w;
    }
    
    FILE* f = fopen(path, "w"); 
    if (!f){ free(x); return 0; }
    fprintf(f, "# FFT of scope (Hann); dt=%.9e; N=%d\n", dt_, Nfft);
    fprintf(f, "freq_Hz,mag\n");
    
    for (int k=0; k<=Nfft/2; ++k){
        double re=0.0, im=0.0;
        double ang_step = -2.0*M_PI*(double)k/(double)Nfft;
        for (int n=0;n<Nfft;++n){
            double ang = ang_step*(double)n;
            re += x[n]*cos(ang); 
            im += x[n]*sin(ang);
        }
        double mag = sqrt(re*re + im*im);
        double freq = (double)k / (dt_ * (double)Nfft);
        fprintf(f, "%.12e,%.12e\n", freq, mag);
    }
    fclose(f); 
    free(x);
    return 1;
}

/* ================= Text Rendering ================= */
static SDL_Texture* render_text(SDL_Renderer* ren, TTF_Font* font,
                                const char* s, SDL_Color col, int* w, int* h)
{
    if (!font) return NULL;  /* Fixed: check font validity */
    SDL_Surface* surf = TTF_RenderUTF8_Blended(font, s, col);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ren, surf);
    if (w) *w = surf->w; 
    if (h) *h = surf->h;
    SDL_FreeSurface(surf);
    return tex;
}

static void draw_legend(SDL_Renderer* ren, TTF_Font* font, int x, int y,
                        const char** lines, int nlines)
{
    if (!font) return;  /* Fixed: check font validity */
    
    int maxw = 0, lineh = 0;
    for (int i=0;i<nlines;++i){
        int w,h; 
        TTF_SizeUTF8(font, lines[i], &w, &h);
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

/* ================= Screenshot ================= */
static void save_screenshot(SDL_Renderer* ren, const char* path, int w, int h) {
    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surf) return;
    SDL_RenderReadPixels(ren, NULL, SDL_PIXELFORMAT_ARGB8888, surf->pixels, surf->pitch);
    SDL_SaveBMP(surf, path);
    SDL_FreeSurface(surf);
}

/* ================= Probe ================= */
static FILE *probe_out = NULL;
static int   probe_x = NX/2 + NX/6, probe_y = NY/2;
static int probe2_active = 0;
static int probe2_x = NX/2 - NX/6;
static int probe2_y = NY/2;

/* ================= Visualization ================= */
typedef enum { AS_PEAK=0, AS_P99=1 } AutoScaleMode;
#define COLOR_HIST_BINS 512
#define HIST_STRIDE 4

static AutoScaleMode color_autoscale_mode = AS_P99;
static int hold_color = 0;
static double held_vmax = 1e-3;
static int hold_scope = 0;
static double held_scope_vmax = 1e-3;
static double vmax_smooth = 0.0, p99_smooth = 0.0, scope_vmax_smooth = 0.0;
static int g_render_stride = 1;
static int hist_frame_counter = 0;  /* Fixed: only update histogram periodically */

/* ================= CPML ================= */
#define PML_THICK 12
static int cpml_on = 0;
static int cpml_N  = PML_THICK;
static double pml_sigma_max = 1.2, pml_kappa_max = 5.0, pml_alpha_max = 0.05;

/* PML coefficients */
static double kx[NX], bx[NX], cx[NX];
static double ky[NY], by[NY], cy[NY];

/* PML auxiliary fields */
static double (*psi_Ezx)[NY];
static double (*psi_Ezy)[NY];
static double (*psi_Hyx)[NY];
static double (*psi_Hxy)[NY];

typedef struct { 
    const char* name; 
    double smax, kmax, amax; 
    int thick; 
} CpmlPreset;

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

/* Fixed PML coefficients for proper Yee grid positioning */
static void cpml_build_coeffs(double dt_){
    const double m = 3.0;
    
    /* X-direction coefficients (for Ex, Hy fields) */
    for (int i=0;i<NX;++i){
        double rx = 0.0;
        if (i < cpml_N)            
            rx = (cpml_N - i) / (double)cpml_N;  /* Fixed: removed 0.5 offset */
        else if (i >= NX - cpml_N) 
            rx = (i - (NX - cpml_N - 1)) / (double)cpml_N;
            
        double g = rx>0 ? pow(rx, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - rx);
        
        kx[i] = kappa;
        if (rx>0){
            double denom = kappa + (sigma + kappa*alpha)*dt_/(2.0*EPS0);
            bx[i] = exp(-(sigma/kappa + alpha) * dt_);
            cx[i] = (sigma * (bx[i] - 1.0)) / (denom + 1e-30);
        } else { 
            bx[i]=1.0; cx[i]=0.0; 
        }
    }
    
    /* Y-direction coefficients (for Ey, Hx fields) */
    for (int j=0;j<NY;++j){
        double ry = 0.0;
        if (j < cpml_N)            
            ry = (cpml_N - j) / (double)cpml_N;  /* Fixed: removed 0.5 offset */
        else if (j >= NY - cpml_N) 
            ry = (j - (NY - cpml_N - 1)) / (double)cpml_N;
            
        double g = ry>0 ? pow(ry, m) : 0.0;
        double sigma = pml_sigma_max * g;
        double kappa = 1.0 + (pml_kappa_max - 1.0) * g;
        double alpha = pml_alpha_max * (1.0 - ry);
        
        ky[j] = kappa;
        if (ry>0){
            double denom = kappa + (sigma + kappa*alpha)*dt_/(2.0*EPS0);
            by[j] = exp(-(sigma/kappa + alpha) * dt_);
            cy[j] = (sigma * (by[j] - 1.0)) / (denom + 1e-30);
        } else { 
            by[j]=1.0; cy[j]=0.0; 
        }
    }
}

static void cpml_apply_preset(int idx, double dt_){
    int n = (int)(sizeof(CPML_PRESETS)/sizeof(CPML_PRESETS[0]));
    if (idx < 0) idx = 0; 
    if (idx >= n) idx = n-1;
    cpml_preset_idx = idx;
    pml_sigma_max = CPML_PRESETS[idx].smax;
    pml_kappa_max = CPML_PRESETS[idx].kmax;
    pml_alpha_max = CPML_PRESETS[idx].amax;
    cpml_N        = CPML_PRESETS[idx].thick;
    cpml_build_coeffs(dt_);
    cpml_zero_psi();
}

/* ================= Resource Cleanup ================= */
static void cleanup_resources(void) {
    /* Fixed: proper cleanup of all resources */
    for (int p = 0; p < MAX_PORTS; ++p) {
        if (ports[p].V) { free(ports[p].V); ports[p].V = NULL; }
        if (ports[p].I) { free(ports[p].I); ports[p].I = NULL; }
    }
    
    if (probe_out) { fclose(probe_out); probe_out = NULL; }
    
    scope_free();
    
    if (Ez) { free(Ez); Ez = NULL; }
    if (Hx) { free(Hx); Hx = NULL; }
    if (Hy) { free(Hy); Hy = NULL; }
    if (Ez_old) { free(Ez_old); Ez_old = NULL; }
    if (epsr) { free(epsr); epsr = NULL; }
    if (sigma_map) { free(sigma_map); sigma_map = NULL; }
    if (tag_grid) { free(tag_grid); tag_grid = NULL; }
    
    if (psi_Ezx) { free(psi_Ezx); psi_Ezx = NULL; }
    if (psi_Ezy) { free(psi_Ezy); psi_