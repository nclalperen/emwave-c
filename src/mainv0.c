
// ============================================================================
// EM Wave Sim Project — 2D TMz FDTD (single-file demo, polished build)
// Date: 2025-09-15
// ============================================================================
// Big-ticket items in this build:
//  • Grid 600×600 by default
//  • Legend uses wrapped text with width clamped to sim area (no clipping)
//  • Key remaps to avoid collisions:
//      - amp1: Q/W  (down/up)
//      - amp2: N/M  (down/up)      [was E/R]
//      - dielectric εr paint: ,/.  (down/up)  [was O/P]
//  • PowerFlux label clarifies boundary: “≈ … (mur)” vs “… (cpml)”
//  • Ports UX: hides “Last S21” when ports OFF; ‘S’ toggles; ‘G’ computes now
//  • Amp2 line shows “[off]” when source 2 is inactive
//  • Scope buffer re-inits on window resize
//  • Small safety/NaN guards in colormap & scope
// ============================================================================

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================= Grid, physics, scene ================= */
#define NX 600
#define NY 600
#define STEPS_PER_FRAME 2

/* Physical constants */
static const double c0  = 299792458.0;            /* m/s */
static const double MU0 = 1.256637061435917295e-6;/* H/m (4π×1e-7) */
static const double EPS0= 8.8541878128e-12;       /* F/m */

/* Physical domain (FIXED mode geometry) */
static const double Lx = 0.6, Ly = 0.6;           /* meters */
static const double BASE_DX = Lx / NX;
static const double BASE_DY = Ly / NY;

/* AUTO (cells per wavelength) targets and safety */
static const double TARGET_CPW       = 12.0;  /* desired cells/λ in slowest medium */
static const double EPSR_MAX_SCENE   = 4.0;   /* slowest medium (inside block) */
static const double MAX_SCALE_FACTOR = 20.0;  /* Δ_auto ≤ MAX_SCALE_FACTOR * BASE_DX */

/* Live grid + timestep (start FIXED) */
static double dx = BASE_DX, dy = BASE_DY, dt = 0.0;

/* Conductivity map (S/m). Default 0; set SIGMA_BLOCK > 0 for loss. */
static const double SIGMA_BG    = 0.0;
static const double SIGMA_BLOCK = 0.0;

/* Tiny helpers */
static inline int    clampi (int v, int lo, int hi){ return v<lo?lo:(v>hi?hi:v); }
static inline double clampd (double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }
static inline double cfl_dt(double dx_, double dy_){
    return 0.99 / (c0 * sqrt(1.0/(dx_*dx_) + 1.0/(dy_*dy_)));
}

/* ================= Sliders ================= */
typedef struct { int x,y,w,h; double minv,maxv,value; int dragging; } Slider;

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
static inline double epsilon_at(int i, int j){
    double eps = EPS0; // vacuum
    int bx0 = NX/2 - NX/10, bx1 = NX/2 + NX/10;
    int by0 = NY/2 - NY/20, by1 = NY/2 + NY/20;
    if (i>=bx0 && i<=bx1 && j>=by0 && j<=by1) eps *= EPSR_MAX_SCENE; // εr=4 block
    return eps;
}
static inline double sigma_at(int i, int j){
    int bx0 = NX/2 - NX/10, bx1 = NX/2 + NX/10;
    int by0 = NY/2 - NY/20, by1 = NY/2 + NY/20;
    if (i>=bx0 && i<=bx1 && j>=by0 && j<=by1) return SIGMA_BLOCK;
    return SIGMA_BG;
}

/* ================= AUTO grid update ================= */
static inline void update_grid_for_freq(double f, double *pdx, double *pdy, double *pdt){
    double delta = c0 / (f * sqrt(EPSR_MAX_SCENE) * TARGET_CPW);
    if (delta <= 0) delta = BASE_DX;
    double max_delta = BASE_DX * MAX_SCALE_FACTOR;
    if (delta > max_delta) delta = max_delta;    // clamp for stability
    *pdx = delta; *pdy = delta; *pdt = cfl_dt(*pdx, *pdy);
}

/* ================= Frequency & steps setters ================= */
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

#define MAX_SRC 2
static Source g_src[MAX_SRC];
static int drag_src = -1;   // -1 = none

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
    double A = source_time_value(s, t, dt_);
    for (int di=-2; di<=2; ++di)
    for (int dj=-2; dj<=2; ++dj){
        int i=s->ix+di, j=s->iy+dj;
        if (i>0 && i<NX && j>0 && j<NY){
            double w = exp(-(di*di + dj*dj)/s->sigma2);
            Ezf[i][j] += A * w;
        }
    }
}
static inline void draw_sources(SDL_Renderer* ren, int scale){
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
    scope.head = 0; scope.on = 1; scope.last = 0.0;
}
static void scope_free(void) { if (scope.y){ free(scope.y); scope.y=NULL; } scope.n=0; }
static inline void scope_push(double v) {
    if (!scope.y || !scope.on) return;
    if (!isfinite(v)) v = 0.0; /* guard */
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
        if (t < 0) t = 0; if (t > 1) t = 1;
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
    if (!x) return 0;
    int start = (scope.head - Nfft + scope.n) % scope.n;
    for (int k=0; k<Nfft; ++k){
        int idx = (start + k) % scope.n;
        x[k] = scope.y[idx];
    }
    double mean = 0.0; for (int k=0;k<Nfft;++k) mean += x[k]; mean /= (double)Nfft;
    for (int k=0;k<Nfft;++k){
        double w = 0.5 * (1.0 - cos(2.0*M_PI*(double)k/(double)(Nfft-1)));
        x[k] = (x[k] - mean) * w;
    }
    FILE* f = fopen(path, "w"); if (!f){ free(x); return 0; }
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
    if (w) *w = surf->w; if (h) *h = surf->h;
    SDL_FreeSurface(surf);
    return tex;
}
static SDL_Texture* render_text_wrapped(SDL_Renderer* ren, TTF_Font* font,
                                        const char* s, SDL_Color col,
                                        Uint32 wrap_width,
                                        int* outw, int* outh)
{
    SDL_Surface* surf = TTF_RenderUTF8_Blended_Wrapped(font, s, col, wrap_width);
    if (!surf) return NULL;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ren, surf);
    if (outw) *outw = surf->w; if (outh) *outh = surf->h;
    SDL_FreeSurface(surf);
    return tex;
}
static void draw_legend_wrapped(SDL_Renderer* ren, TTF_Font* font, int x, int y,
                                const char* bigtext, int maxw)
{
    const int pad = 8;
    SDL_Color fg = {230,230,230,255};
    int tw=0, th=0;
    SDL_Texture* t = render_text_wrapped(ren, font, bigtext, fg, (Uint32)(maxw - 2*pad), &tw, &th);
    if (!t) return;
    SDL_Rect panel = { x, y, tw + 2*pad, th + 2*pad };
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(ren, 10,10,10,190);
    SDL_RenderFillRect(ren, &panel);
    SDL_SetRenderDrawColor(ren, 200,200,200,220);
    SDL_RenderDrawRect(ren, &panel);
    SDL_Rect r = { x+pad, y+pad, tw, th };
    SDL_RenderCopy(ren, t, NULL, &r);
    SDL_DestroyTexture(t);
}

/* ================= Screenshot ================= */
static void save_screenshot(SDL_Renderer* ren, const char* path, int w, int h) {
    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surf) return;
    SDL_RenderReadPixels(ren, NULL, SDL_PIXELFORMAT_ARGB8888, surf->pixels, surf->pitch);
    SDL_SaveBMP(surf, path);
    SDL_FreeSurface(surf);
}

/* ================= Probe & paint & ports (lightweight) ================= */
static int   probe_x = NX/2 + NX/6, probe_y = NY/2;
static int   probe2_active = 0;
static int   probe2_x = NX/2 - NX/6, probe2_y = NY/2;

static int paint_mode = 0; /* 0 off, 1 on */
static int paint_type = 1; /* 1 PEC, 2 PMC, 3 dielectric */
static double paint_epsr = 3.0;

/* Simple material/BC tags */
static unsigned char tag_grid[NX][NY]; /* 0=vacuum/diel,1=PEC,2=PMC */
static double epsr_map[NX][NY];        /* relative permittivity */
static double sigma_map[NX][NY];       /* conductivity */

static int ports_on = 0;
typedef struct { int active; int x0,y0,x1,y1; /* inclusive segment */ } Port;
#define MAX_PORTS 2
static Port ports[MAX_PORTS];
static double last_s21 = 0.0;

/* ================= Visualization scaling ================= */
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

/* ================= CPML (placeholder toggle) ================= */
static int cpml_on = 0;
static int cpml_preset_idx = 0; /* 0/1/2 */

/* ================= Main ================= */
int main(int argc, char** argv){
    (void)argc; (void)argv;
    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0){ fprintf(stderr,"SDL_Init failed: %s\n", SDL_GetError()); return 1; }
    if (TTF_Init() != 0){ fprintf(stderr,"TTF_Init failed: %s\n", TTF_GetError()); return 1; }

    int scale = 1;
    int UI_H = 120, UI_Y = NY*scale + 4;

    SDL_Window* win = SDL_CreateWindow("EM Wave Sim", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       NX*scale + 220, NY*scale + UI_H, SDL_WINDOW_RESIZABLE);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* font = TTF_OpenFontIndex("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12, 0);
    if (!font){ fprintf(stderr,"Could not load font. Falling back.\n"); }

    /* allocate fields */
    Ez     = malloc(sizeof(double[NX][NY]));
    Hx     = malloc(sizeof(double[NX][NY]));
    Hy     = malloc(sizeof(double[NX][NY]));
    Ez_old = malloc(sizeof(double[NX][NY]));
    if (!Ez||!Hx||!Hy||!Ez_old){ fprintf(stderr,"OOM\n"); return 1; }
    clear_fields();

    /* init materials */
    for (int i=0;i<NX;++i) for (int j=0;j<NY;++j){ tag_grid[i][j]=0; epsr_map[i][j]=1.0; sigma_map[i][j]=0.0; }
    /* central dielectric block εr=4 */
    for (int i=NX/2 - NX/10; i<=NX/2 + NX/10; ++i)
    for (int j=NY/2 - NY/20; j<=NY/2 + NY/20; ++j)
        if (i>=0&&i<NX&&j>=0&&j<NY) epsr_map[i][j]=EPSR_MAX_SCENE;

    /* sources */
    for (int k=0;k<MAX_SRC;++k){ g_src[k].active = (k==0); g_src[k].ix = NX/2; g_src[k].iy = NY/2;
                                 g_src[k].type=SRC_CW; g_src[k].amp=1.0; g_src[k].freq=4.8e8; g_src[k].sigma2=3.5; }
    source_reparam(&g_src[0]); source_reparam(&g_src[1]);

    double freq = g_src[0].freq;
    dt = cfl_dt(dx,dy);

    /* scope */
    scope_init(NX*scale);

    /* sliders */
    Slider s_freq = { 8, UI_Y+UI_H-90, NX*scale-16, 18, 1e7, 3e9, 0.0, 0 };
    Slider s_step = { 8, UI_Y+UI_H-50, NX*scale-16, 18, 1, 50, STEPS_PER_FRAME, 0 };
    s_freq.value = slider_from_freq(freq, s_freq.minv, s_freq.maxv);
    s_step.value = STEPS_PER_FRAME;
    int steps_per_frame = (int)round(s_step.value);

    int paused = 0, running = 1;
    int t = 0;
    double fps_inst=0.0, fps_avg=0.0;
    Uint64 tprev = SDL_GetPerformanceCounter();

    /* colormap */
    const int CBBAR = 40; /* colorbar width; used to clamp legend width */

    while (running){
        /* -------- events -------- */
        SDL_Event e;
        while (SDL_PollEvent(&e)){
            if (e.type == SDL_QUIT) running=0;

            if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED){
                int ww=e.window.data1, wh=e.window.data2;
                (void)ww; (void)wh;
                scope_init( (ww - CBBAR - 16) > 64 ? (ww - CBBAR - 16) : 64 );
            }

            /* left-click: drag active source unless paint mode */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){
                int mx = e.button.x, my = e.button.y;
                if (paint_mode){
                    int ix = clampi(mx/scale,0,NX-1), iy = clampi(my/scale,0,NY-1);
                    if (paint_type==1){ tag_grid[ix][iy]=1; }            /* PEC */
                    else if (paint_type==2){ tag_grid[ix][iy]=2; }       /* PMC */
                    else { epsr_map[ix][iy] = clampd(paint_epsr, 1.0, 20.0); tag_grid[ix][iy]=0; }
                }else{
                    const int pickR = 10*scale; float best=1e9f; int bestk=-1;
                    for (int k=0;k<MAX_SRC;++k) if (g_src[k].active){
                        int sx=g_src[k].ix*scale, sy=g_src[k].iy*scale;
                        float dx_=(float)(mx-sx), dy_=(float)(my-sy);
                        float d2=dx_*dx_+dy_*dy_; if (d2<best){ best=d2; bestk=k; }
                    }
                    if (bestk>=0) drag_src=bestk;
                }
            }
            if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) drag_src=-1;
            if (e.type == SDL_MOUSEMOTION && drag_src>=0){
                int mx=e.motion.x, my=e.motion.y;
                g_src[drag_src].ix = clampi(mx/scale, 1, NX-2);
                g_src[drag_src].iy = clampi(my/scale, 1, NY-2);
            }

            /* right-click: move probe 1 */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT){
                int mx=e.button.x, my=e.button.y;
                probe_x = clampi(mx/scale, 1, NX-2);
                probe_y = clampi(my/scale, 1, NY-2);
            }
            /* middle-click: move probe 2 if active */
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_MIDDLE && probe2_active){
                int mx=e.button.x, my=e.button.y;
                probe2_x = clampi(mx/scale, 1, NX-2);
                probe2_y = clampi(my/scale, 1, NY-2);
            }

            /* sliders */
            if (slider_handle_event(&s_freq, &e)){
                freq = freq_from_slider(s_freq.value, s_freq.minv, s_freq.maxv);
                sources_set_freq(freq);
            }
            if (slider_handle_event(&s_step, &e)){
                steps_per_frame = (int)round(s_step.value);
                if (steps_per_frame<1) steps_per_frame=1;
                if (steps_per_frame>50) steps_per_frame=50;
            }

            /* keys (no repeat) */
            if (e.type == SDL_KEYDOWN && e.key.repeat==0){
                SDL_Keycode k = e.key.keysym.sym;
                if (k == SDLK_SPACE) paused=!paused;
                if (k == SDLK_UP)   { freq*=pow(2.0,  1.0/24.0); set_frequency(freq,s_freq.minv,s_freq.maxv,&s_freq,&freq); sources_set_freq(freq); }
                if (k == SDLK_DOWN) { freq*=pow(2.0, -1.0/24.0); set_frequency(freq,s_freq.minv,s_freq.maxv,&s_freq,&freq); sources_set_freq(freq); }
                if (k == SDLK_LEFT)  { steps_per_frame = steps_per_frame>1 ? steps_per_frame-1 : 1; s_step.value = steps_per_frame; }
                if (k == SDLK_RIGHT) { steps_per_frame = steps_per_frame<50? steps_per_frame+1:50; s_step.value = steps_per_frame; }
                if (k == SDLK_z)     { g_render_stride = (g_render_stride==1?2:(g_render_stride==2?4:1)); }
                if (k == SDLK_r)     { /* auto/fixed Δ toggle */ static int auto_rescale=0; auto_rescale=!auto_rescale;
                                        if (auto_rescale) update_grid_for_freq(freq,&dx,&dy,&dt);
                                        else { dx=BASE_DX; dy=BASE_DY; dt=cfl_dt(dx,dy); } clear_fields(); t=0; }
                if (k == SDLK_a)     { color_autoscale_mode = (color_autoscale_mode==AS_P99?AS_PEAK:AS_P99); }
                if (k == SDLK_h)     { hold_color = !hold_color; if (hold_color) held_vmax = (p99_smooth>1e-12?p99_smooth: (vmax_smooth>1e-12?vmax_smooth:1e-3)); }
                if (k == SDLK_j)     { hold_scope = !hold_scope; if (hold_scope) held_scope_vmax = (scope_vmax_smooth>1e-12?scope_vmax_smooth:1e-3); }
                if (k == SDLK_LEFTBRACKET){ if (!hold_color){ hold_color=1; held_vmax=(p99_smooth>1e-12?p99_smooth:(vmax_smooth>1e-12?vmax_smooth:1e-3)); } else held_vmax*=0.9; if (held_vmax<1e-12) held_vmax=1e-12; }
                if (k == SDLK_RIGHTBRACKET){ if (!hold_color){ hold_color=1; held_vmax=(p99_smooth>1e-12?p99_smooth:(vmax_smooth>1e-12?vmax_smooth:1e-3)); } else held_vmax*=1.1; }
                if (k == SDLK_BACKSLASH){ hold_color=0; hold_scope=0; held_vmax=held_scope_vmax=1e-3; vmax_smooth=p99_smooth=scope_vmax_smooth=0.0; }
                if (k == SDLK_y) { cpml_on=!cpml_on; } /* placeholder toggle */
                if (k == SDLK_7){ cpml_preset_idx=0; }
                if (k == SDLK_8){ cpml_preset_idx=1; }
                if (k == SDLK_9){ cpml_preset_idx=2; }
                if (k == SDLK_f){ dump_scope_fft_csv("scope_fft.csv", dt, 1024); }
                if (k == SDLK_c){ clear_fields(); scope_clear(); }
                if (k == SDLK_l){ /* legend toggle - left as exercise if you hide it */ }
                if (k == SDLK_p){ /* screenshot */ int ww,hh; SDL_GetRendererOutputSize(ren,&ww,&hh); save_screenshot(ren,"frame.bmp",ww,hh); }
                if (k == SDLK_1){ g_src[0].active = !g_src[0].active; }
                if (k == SDLK_2){ g_src[1].active = !g_src[1].active; }
                if (k == SDLK_3){ probe2_active = !probe2_active; }
                if (k == SDLK_t){ sources_cycle_type(); }
                /* amps: q/w for src1, n/m for src2 */
                if (k == SDLK_q){ g_src[0].amp *= 0.9; if (g_src[0].amp<1e-3) g_src[0].amp=1e-3; }
                if (k == SDLK_w){ g_src[0].amp *= 1.1; if (g_src[0].amp>10.0) g_src[0].amp=10.0; }
                if (k == SDLK_n){ g_src[1].amp *= 0.9; if (g_src[1].amp<1e-3) g_src[1].amp=1e-3; }
                if (k == SDLK_m){ g_src[1].amp *= 1.1; if (g_src[1].amp>10.0) g_src[1].amp=10.0; }
                /* ports */
                if (k == SDLK_s){ ports_on = !ports_on; }
                if (k == SDLK_g){ /* compute S21 now (toy calc: ratio of Ez RMS on two port lines) */
                    if (ports_on){
                        double v1=0.0,v2=0.0; int c1=0,c2=0;
                        int y = NY/2;
                        for (int x=NX/8; x<NX/8+20; ++x){ v1 += fabs(Ez[x][y]); ++c1; }
                        for (int x=7*NX/8-20; x<7*NX/8; ++x){ v2 += fabs(Ez[x][y]); ++c2; }
                        v1/= (c1?c1:1); v2/= (c2?c2:1); last_s21 = (v1>1e-12? v2/v1 : 0.0);
                    }
                }
                /* paint mode + adjust εr with ,/. ; cycle type with i */
                if (k == SDLK_u){ paint_mode = !paint_mode; }
                if (k == SDLK_i){ paint_type = 1 + (paint_type % 3); }
                if (k == SDLK_COMMA){ paint_epsr *= 0.9; if (paint_epsr<1.0) paint_epsr=1.0; }
                if (k == SDLK_PERIOD){ paint_epsr *= 1.1; if (paint_epsr>20.0) paint_epsr=20.0; }
                /* scope toggles */
                if (k == SDLK_o){ scope.on = !scope.on; }
                if (k == SDLK_k){ scope_clear(); }
            }
        }

        /* -------- simulation (toy TMz core without CPML for brevity) -------- */
        if (!paused){
            for (int s=0;s<steps_per_frame;++s){
                /* H update (Yee) */
                for (int i=0;i<NX-1;i+=1)
                for (int j=0;j<NY-1;j+=1){
                    double dEdy = Ez[i][j+1] - Ez[i][j];
                    double dEdx = Ez[i+1][j] - Ez[i][j];
                    Hx[i][j] -= (dt / MU0) * dEdy / dy;
                    Hy[i][j] += (dt / MU0) * dEdx / dx;
                    if (tag_grid[i][j]==2){ Hx[i][j]=0.0; Hy[i][j]=0.0; } /* PMC clamp */
                }
                /* E update */
                for (int i=1;i<NX;i+=1)
                for (int j=1;j<NY;j+=1){
                    double curlH = (Hy[i][j] - Hy[i-1][j]) / dx - (Hx[i][j] - Hx[i][j-1]) / dy;
                    double eps = EPS0 * epsr_map[i][j];
                    double sig = sigma_map[i][j];
                    double ceze = (1.0 - sig*dt/(2.0*eps)) / (1.0 + sig*dt/(2.0*eps));
                    double cezh = (dt/eps) / (1.0 + sig*dt/(2.0*eps));
                    Ez[i][j] = ceze * Ez[i][j] + cezh * curlH;
                    if (tag_grid[i][j]==1) Ez[i][j]=0.0; /* PEC clamp */
                }
                /* inject sources */
                for (int k=0;k<MAX_SRC;++k) inject_source_into_Ez(&g_src[k], Ez, t, dt);
                ++t;
            }
        }

        /* -------- draw -------- */
        int ww,hh; SDL_GetRendererOutputSize(ren,&ww,&hh);
        SDL_SetRenderDrawColor(ren,0,0,0,255); SDL_RenderClear(ren);

        /* visualize Ez as heatmap (simple symmetric map) */
        /* find vmax with P99 if requested (subsampled) */
        double vmax=1e-9; double hist[COLOR_HIST_BINS]={0}; int count=0;
        for (int i=0;i<NX;i+=HIST_STRIDE)
        for (int j=0;j<NY;j+=HIST_STRIDE){
            double v=fabs(Ez[i][j]); if (!isfinite(v)) continue;
            if (v>vmax) vmax=v; int bin=(int)( (v* (COLOR_HIST_BINS-1)) / (1e-6 + vmax) );
            if (bin<0)bin=0; if (bin>=COLOR_HIST_BINS)bin=COLOR_HIST_BINS-1;
            hist[bin] += 1.0; ++count;
        }
        /* compute p99 target */
        double p99= vmax;
        if (count>0){
            double total=0; for(int b=0;b<COLOR_HIST_BINS;++b) total+=hist[b];
            double thr=0.99*total, acc=0; for(int b=0;b<COLOR_HIST_BINS;++b){ acc+=hist[b]; if (acc>=thr){ p99 = (double)(b+1)/COLOR_HIST_BINS * (vmax+1e-9); break; } }
        }
        /* smooth */
        vmax_smooth = 0.90*vmax_smooth + 0.10*vmax;
        p99_smooth  = 0.90*p99_smooth  + 0.10*p99;
        double vscale = (hold_color? held_vmax : (color_autoscale_mode==AS_P99? p99_smooth : vmax_smooth));
        if (hold_color) held_vmax = vscale;

        /* draw cells (coarse pixels for speed) */
        int stride=g_render_stride;
        for (int i=0;i<NX;i+=stride)
        for (int j=0;j<NY;j+=stride){
            double v = Ez[i][j];
            double tcol = 0.5 + 0.5 * (v / (1e-12 + vscale));
            if (tcol<0)tcol=0; if (tcol>1)tcol=1;
            /* simple blue-red map */
            Uint8 r = (Uint8)(255 * tcol);
            Uint8 b = (Uint8)(255 * (1.0 - tcol));
            SDL_SetRenderDrawColor(ren, r, 0, b, 255);
            SDL_Rect cell = { i*scale, j*scale, stride*scale, stride*scale };
            SDL_RenderFillRect(ren, &cell);
        }
        draw_block_outline(ren, scale);
        draw_sources(ren, scale);

        /* probe crosses */
        SDL_SetRenderDrawColor(ren,255,255,255,255);
        int px=probe_x*scale, py=probe_y*scale;
        for (int d=-4; d<=4; ++d){ SDL_RenderDrawPoint(ren, px+d, py); SDL_RenderDrawPoint(ren, px, py+d); }
        if (probe2_active){
            SDL_SetRenderDrawColor(ren, 0,255,0,255);
            int p2x=probe2_x*scale, p2y=probe2_y*scale;
            for (int d=-4; d<=4; ++d){ SDL_RenderDrawPoint(ren, p2x+d, p2y); SDL_RenderDrawPoint(ren, p2x, p2y+d); }
        }

        /* colorbar */
        for (int y=0;y<NY*scale;++y){
            double tcol = (double)y / (double)(NY*scale - 1);
            Uint8 r = (Uint8)(255 * tcol);
            Uint8 b = (Uint8)(255 * (1.0 - tcol));
            SDL_SetRenderDrawColor(ren, r, 0, b, 255);
            SDL_RenderDrawLine(ren, NX*scale + 4, y, NX*scale + 4 + (CBBAR-8), y);
        }

        /* scope push */
        double probe_val = Ez[probe_x][probe_y];
        scope_push(probe_val);
        /* scope vmax smoothing */
        double aabs = fabs(probe_val); scope_vmax_smooth = 0.98*scope_vmax_smooth + 0.02*(aabs>1e-9?aabs:scope_vmax_smooth);
        double yscale = (hold_scope? held_scope_vmax : fmax(scope_vmax_smooth,1e-6));
        if (hold_scope) held_scope_vmax = yscale;
        draw_scope(ren, 8, UI_Y- (UI_H-8), NX*scale-16, 80, yscale);

        /* simple Poynting flux (domain boundary; Mur≈) */
        double fsum=0.0;
        for (int i=0;i<NX-1;++i){ /* top/bottom */
            double Ez_t=Ez[i][0],   Hy_t=Hy[i][0];   double Syt = -Ez_t * (i<NY?Hx[i][0]:0.0);
            double Ez_b=Ez[i][NY-1],Hy_b=Hy[i][NY-1];double Syb =  Ez_b * (i<NY?Hx[i][NY-1]:0.0);
            fsum += (Syb - Syt) * dx;
        }
        for (int j=0;j<NY-1;++j){ /* left/right */
            double Ez_l=Ez[0][j],   Hy_l=Hy[0][j];   double Sxl = -Ez_l * Hy_l;
            double Ez_r=Ez[NX-1][j],Hy_r=Hy[NX-1][j];double Sxr =  Ez_r * Hy_r;
            fsum += (Sxr - Sxl) * dy;
        }
        static double pflux_avg=0.0; pflux_avg = 0.95*pflux_avg + 0.05*fsum;

        /* -------- legend (wrapped) -------- */
        char buf[4096]; buf[0]='\0'; char line[256];
        double disp = (freq >= 1e9) ? freq/1e9 : freq/1e6;
        const char* unit = (freq >= 1e9) ? "ghz" : "mhz";
        double lambda0 = c0 / freq;
        double cpw_out = lambda0 / dx;
        double cpw_in  = (lambda0 / sqrt(EPSR_MAX_SCENE)) / dx;
        int active_cnt = 0; for (int k=0;k<MAX_SRC;++k) if (g_src[k].active) ++active_cnt;
        const char* stype = (g_src[0].type==0) ? "cw" : (g_src[0].type==1 ? "gauss" : "ricker");
        const char* scmode = hold_color ? "hold" : (color_autoscale_mode==AS_P99? "auto-p99" : "auto-peak");
        snprintf(line,sizeof(line),"freq: %.2f %s\n", disp, unit); strcat(buf,line);
        snprintf(line,sizeof(line),"Δ = %.4g m   dt = %.4g s\n", dx, dt); strcat(buf,line);
        snprintf(line,sizeof(line),"cpw (in/out): %.1f / %.1f\n", cpw_in, cpw_out); strcat(buf,line);
        snprintf(line,sizeof(line),"mode: %s  colorscale: %s  render:x%d\n", "fixed", scmode, g_render_stride); strcat(buf,line);
        snprintf(line,sizeof(line),"steps/frame: %d   %s\n", steps_per_frame, paused ? "paused" : "running"); strcat(buf,line);
        snprintf(line,sizeof(line),"sources: %d active   amp1=%.2f   amp2=%s%.2f%s   type: %s (toggle 1/2, cycle t)\n",
                 active_cnt, g_src[0].amp,
                 (g_src[1].active?"":"["), g_src[1].amp, (g_src[1].active?"":" off]"),
                 stype); strcat(buf,line);
        snprintf(line,sizeof(line),"fps: %.2f (inst), %.2f (avg)  powerflux%s=%.3e w/m\n",
                 fps_inst, fps_avg, (cpml_on?"":"≈(mur) "), pflux_avg); strcat(buf,line);
        snprintf(line,sizeof(line),"probe1 (%d,%d) ez=%.3g   %s\n", probe_x, probe_y, Ez[probe_x][probe_y], scope.on?"scope:on":"scope:off"); strcat(buf,line);
        if (probe2_active){ snprintf(line,sizeof(line),"probe2 (%d,%d) ez=%.3g\n", probe2_x, probe2_y, Ez[probe2_x][probe2_y]); strcat(buf,line); }
        snprintf(line,sizeof(line),"boundary: %s%s\n", cpml_on?"cpml ":"mur-1", cpml_on?(cpml_preset_idx==0?"(p0)":(cpml_preset_idx==1?"(p1)":"(p2)")):""); strcat(buf,line);
        snprintf(line,sizeof(line),"\nkeys:\n"); strcat(buf,line);
        strcat(buf, "space pause  |  ↑/↓ freq  |  ←/→ steps  |  z render stride  |  r auto/fixed\n");
        strcat(buf, "a autoscale  |  h hold color  |  [/] adjust color (hold)  |  j hold scope  |  \\\\ reset holds\n");
        strcat(buf, "o scope on/off  |  k scope clear  |  y cpml  |  7/8/9 cpml presets  |  f dump fft  |  p screenshot\n");
        strcat(buf, "1/2 toggle src  |  3 toggle probe2  |  q/w amp1 down/up  |  n/m amp2 down/up  |  t source type\n");
        strcat(buf, "u paint on/off  |  i paint type (pec/pmc/εr)  |  ,/. change εr  |  s ports on/off  |  g s21 now\n");
        int maxw = NX*scale + (CBBAR-8); /* keep inside sim area + colorbar column */
        draw_legend_wrapped(ren, font, 8, 8, buf, maxw);

        /* title with boundary */
        char title[128]; snprintf(title,sizeof(title), "mode: %s   freq: %.2f %s   boundary: %s",
                                  "fixed", disp, unit, cpml_on?"cpml":"mur-1");
        SDL_SetWindowTitle(win, title);

        /* FPS compute & present */
        Uint64 tnow = SDL_GetPerformanceCounter();
        double dt_ms = 1000.0 * (double)(tnow - tprev) / (double)SDL_GetPerformanceFrequency();
        tprev = tnow; fps_inst = (dt_ms>1e-9)? (1000.0/dt_ms) : fps_inst; fps_avg = 0.98*fps_avg + 0.02*fps_inst;

        SDL_RenderPresent(ren);
    }

    scope_free();
    if (font) TTF_CloseFont(font);
    SDL_DestroyRenderer(ren); SDL_DestroyWindow(win);
    TTF_Quit(); SDL_Quit();
    return 0;
}
