// ============================================================================
// emwave-c: SDL Input Handling and UI State Management
// ============================================================================

#include "ui_controls.h"

#include "config.h"
#include "fdtd_core.h"
#include "materials.h"
#include "sources.h"
#include "analysis.h"
#include "boundary.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static double clampd(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static int clampi_local(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

UIState* ui_state_init(void) {
    UIState* ui = (UIState*)calloc(1, sizeof(UIState));
    if (!ui) {
        return NULL;
    }

    ui->running = 1;
    ui->paused = 0;
    ui->show_legend = 1;
    ui->paint_mode = 0;
    ui->paint_type = 1;
    ui->paint_eps = 2.0;
    ui->auto_rescale = 1;
    ui->hold_color = 0;
    ui->hold_scope = 0;
    ui->color_autoscale_mode = AS_P99;
    ui->render_stride = 1;
    ui->probe_x = NX / 2;
    ui->probe_y = NY / 2;
    ui->probe2_active = 0;
    ui->drag_src = -1;
    ui->sweep_on = 0;
    ui->sweep_idx = 0;
    ui->sweep_points = 0;
    ui->sweep_steps_per_point = 2000;
    ui->steps_per_frame = STEPS_PER_FRAME;
    ui->scale = 2;
    ui->ui_height = 90;
    ui->side_panel_width = 240;
    ui->log_probe = 0;

    ui->freq_slider.minv = 0.0;
    ui->freq_slider.maxv = 1.0;
    ui->freq_slider.value = 0.5;
    ui->freq_slider.dragging = 0;

    ui->speed_slider.minv = 1.0;
    ui->speed_slider.maxv = 50.0;
    ui->speed_slider.value = (double)ui->steps_per_frame;
    ui->speed_slider.dragging = 0;

    return ui;
}

void ui_state_free(UIState* state) {
    if (state) {
        free(state);
    }
}

void ui_state_set_layout(UIState* ui, int scale, int ui_height, int side_panel_width) {
    if (!ui) return;
    ui->scale = scale;
    ui->ui_height = ui_height;
    ui->side_panel_width = side_panel_width;

    int canvas_height = NY * scale;
    int slider_width = NX * scale - 40;
    if (slider_width < 100) slider_width = 100;

    ui->freq_slider.x = 20;
    ui->freq_slider.y = canvas_height + 15;
    ui->freq_slider.w = slider_width;
    ui->freq_slider.h = 22;

    ui->speed_slider.x = 20;
    ui->speed_slider.y = canvas_height + 50;
    ui->speed_slider.w = slider_width;
    ui->speed_slider.h = 22;
}

void ui_state_sync_with_sim(UIState* ui, const SimulationState* sim) {
    if (!ui || !sim) return;
    ui->freq_slider.value = slider_from_freq(sim->freq, FREQ_MIN, FREQ_MAX);
    ui->speed_slider.value = (double)ui->steps_per_frame;
}

void ui_update_metrics(UIState* ui, const SimulationState* sim, const Scope* scope) {
    if (!ui || !sim) return;

    double vmax = 0.0;
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double val = fabs(sim->Ez[i][j]);
            if (val > vmax) vmax = val;
        }
    }

    if (ui->vmax_smooth == 0.0) ui->vmax_smooth = vmax;
    ui->vmax_smooth = 0.9 * ui->vmax_smooth + 0.1 * vmax;

    if (!ui->hold_color) {
        ui->held_vmax = ui->vmax_smooth;
    }

    if (scope && scope->y && scope->n > 0) {
        double smax = 0.0;
        for (int k = 0; k < scope->n; ++k) {
            double val = fabs(scope->y[k]);
            if (val > smax) smax = val;
        }
        if (ui->scope_vmax_smooth == 0.0) ui->scope_vmax_smooth = smax;
        ui->scope_vmax_smooth = 0.9 * ui->scope_vmax_smooth + 0.1 * smax;
        if (!ui->hold_scope) {
            ui->held_scope_vmax = ui->scope_vmax_smooth;
        }
    }
}

static void apply_paint(UIState* ui, SimulationState* sim, int mx, int my, int scale) {
    if (!ui->paint_mode) return;
    if (mx < 0 || my < 0) return;
    if (mx >= NX * scale || my >= NY * scale) return;
    int gx = mx / scale;
    int gy = my / scale;
    paint_material_at(sim, gx, gy, ui->paint_type, ui->paint_eps);
}

static void update_frequency_from_slider(UIState* ui, SimulationState* sim) {
    double freq = freq_from_slider(ui->freq_slider.value, FREQ_MIN, FREQ_MAX);
    fdtd_update_grid_for_freq(sim, freq);
    sources_set_freq(sim->sources, freq);
    ui->held_vmax = 0.0;
}

static void update_steps_from_slider(UIState* ui) {
    ui->steps_per_frame = (int)round(clampd(ui->speed_slider.value, 1.0, 50.0));
    ui->speed_slider.value = (double)ui->steps_per_frame;
}

int ui_handle_events(UIState* ui, SimulationState* sim, Scope* scope,
                     int scale, int ui_height, int side_panel_width) {
    if (!ui || !sim) return 0;

    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) {
            ui->running = 0;
            return 0;
        }

        if (slider_handle_event(&ui->freq_slider, &e)) {
            update_frequency_from_slider(ui, sim);
            continue;
        }
        if (slider_handle_event(&ui->speed_slider, &e)) {
            update_steps_from_slider(ui);
            continue;
        }

        if (e.type == SDL_MOUSEBUTTONDOWN) {
            if (e.button.button == SDL_BUTTON_LEFT) {
                if (ui->paint_mode) {
                    apply_paint(ui, sim, e.button.x, e.button.y, scale);
                } else {
                    if (e.button.x < NX * scale && e.button.y < NY * scale) {
                        ui->probe_x = clampi_local(e.button.x / scale, 0, NX - 1);
                        ui->probe_y = clampi_local(e.button.y / scale, 0, NY - 1);
                    }
                }
            } else if (e.button.button == SDL_BUTTON_RIGHT) {
                ui->paint_mode = !ui->paint_mode;
            }
        } else if (e.type == SDL_MOUSEMOTION) {
            if (ui->paint_mode && (e.motion.state & SDL_BUTTON_LMASK)) {
                apply_paint(ui, sim, e.motion.x, e.motion.y, scale);
            }
        }

        if (e.type == SDL_KEYDOWN) {
            SDL_Keycode sym = e.key.keysym.sym;
            switch (sym) {
                case SDLK_ESCAPE:
                case SDLK_q:
                    ui->running = 0;
                    return 0;
                case SDLK_SPACE:
                    ui->paused = !ui->paused;
                    break;
                case SDLK_l:
                    ui->show_legend = !ui->show_legend;
                    break;
                case SDLK_m:
                    ui->paint_mode = !ui->paint_mode;
                    break;
                case SDLK_1:
                    ui->paint_type = 1;
                    break;
                case SDLK_2:
                    ui->paint_type = 2;
                    break;
                case SDLK_3:
                    ui->paint_type = 3;
                    break;
                case SDLK_o:
                    ui->paint_eps = clampd(ui->paint_eps * 0.9, 1.0, 20.0);
                    break;
                case SDLK_p:
                    ui->paint_eps = clampd(ui->paint_eps * 1.1, 1.0, 20.0);
                    break;
                case SDLK_a:
                    ui->auto_rescale = !ui->auto_rescale;
                    if (!ui->auto_rescale) {
                        ui->held_vmax = ui->vmax_smooth;
                    }
                    break;
                case SDLK_h:
                    ui->hold_color = !ui->hold_color;
                    if (ui->hold_color) {
                        ui->held_vmax = ui->vmax_smooth;
                    }
                    break;
                case SDLK_j:
                    ui->hold_scope = !ui->hold_scope;
                    if (ui->hold_scope) {
                        ui->held_scope_vmax = ui->scope_vmax_smooth;
                    }
                    break;
                case SDLK_c:
                    sources_cycle_type(sim->sources);
                    break;
                case SDLK_r:
                    fdtd_reset(sim);
                    scope_clear(scope);
                    cpml_zero_psi(sim);
                    break;
                case SDLK_y:
                    cpml_on = !cpml_on;
                    boundary_type = cpml_on ? BOUNDARY_CPML : BOUNDARY_MUR;
                    if (cpml_on) {
                        cpml_build_coeffs(sim->dt);
                        cpml_zero_psi(sim);
                    }
                    break;
                case SDLK_g:
                    ui->log_probe = !ui->log_probe;
                    break;
                case SDLK_s:
                    sim->ports_on = !sim->ports_on;
                    break;
                default:
                    break;
            }
        }
    }

    return 1;
}

int slider_handle_event(Slider* s, const SDL_Event* e) {
    if (!s) return 0;
    if (e->type == SDL_MOUSEBUTTONDOWN) {
        int mx = e->button.x;
        int my = e->button.y;
        if (mx >= s->x && mx <= s->x + s->w && my >= s->y && my <= s->y + s->h) {
            s->dragging = 1;
            double t = (double)(mx - s->x) / (double)s->w;
            s->value = s->minv + clampd(t, 0.0, 1.0) * (s->maxv - s->minv);
            return 1;
        }
    } else if (e->type == SDL_MOUSEBUTTONUP) {
        s->dragging = 0;
    } else if (e->type == SDL_MOUSEMOTION && s->dragging) {
        int mx = e->motion.x;
        double t = (double)(mx - s->x) / (double)s->w;
        double new_val = s->minv + clampd(t, 0.0, 1.0) * (s->maxv - s->minv);
        if (fabs(new_val - s->value) > 1e-6) {
            s->value = new_val;
            return 1;
        }
    }
    return 0;
}

double freq_from_slider(double t, double fmin, double fmax) {
    double a = log10(fmin);
    double b = log10(fmax);
    return pow(10.0, a + clampd(t, 0.0, 1.0) * (b - a));
}

double slider_from_freq(double f, double fmin, double fmax) {
    double a = log10(fmin);
    double b = log10(fmax);
    double v = (log10(f) - a) / (b - a);
    return clampd(v, 0.0, 1.0);
}
