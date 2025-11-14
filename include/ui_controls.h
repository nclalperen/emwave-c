// ============================================================================
// emwave-c: UI Controls and Event Handling
// This module handles ALL user input
// Perfect isolation for future UI/UX overhaul
// ============================================================================

#ifndef EMWAVE_UI_CONTROLS_H
#define EMWAVE_UI_CONTROLS_H

#include "types.h"
#include <SDL2/SDL.h>

/* UI state structure */
typedef struct {
    int running;
    int paused;
    int show_legend;
    int paint_mode;
    int paint_type;
    double paint_eps;
    int auto_rescale;

    /* Color and scope holds */
    int hold_color;
    int hold_scope;
    double held_vmax;
    double held_scope_vmax;

    /* Autoscale modes */
    AutoScaleMode color_autoscale_mode;

    /* Smoothed values */
    double vmax_smooth;
    double p99_smooth;
    double scope_vmax_smooth;

    /* Cached values */
    double cached_p99;
    int p99_cache_valid;

    /* Render stride */
    int render_stride;

    /* Probe positions */
    int probe_x, probe_y;
    int probe2_x, probe2_y;
    int probe2_active;

    /* Dragging state */
    int drag_src;  /* -1 = none, 0-3 = source index */

    /* Sweep state */
    int sweep_on;
    int sweep_idx;
    int sweep_points;
    double sweep_freqs[SWEEP_MAX_POINTS];
    double sweep_s21[SWEEP_MAX_POINTS];
    int sweep_steps_remaining;
    int sweep_steps_per_point;

    /* S-parameter */
    double s21_amp;
    int s21_computed;

    /* Sliders */
    Slider freq_slider;
    Slider speed_slider;

} UIState;

/* Initialization */
UIState* ui_state_init(void);
void ui_state_free(UIState* state);

/* Event handling */
int ui_handle_events(UIState* ui, SimulationState* sim, Scope* scope);

/* Slider interaction */
int slider_handle_event(Slider* s, const SDL_Event* e);

/* Slider value conversion helpers */
double freq_from_slider(double t, double fmin, double fmax);
double slider_from_freq(double f, double fmin, double fmax);

#endif /* EMWAVE_UI_CONTROLS_H */
