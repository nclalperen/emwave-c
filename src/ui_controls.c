// ============================================================================
// emwave-c: SDL Input Handling and UI State Management
// ============================================================================

#include "ui_controls.h"

#include "config.h"
#include "ui_layout.h"
#include "fdtd_core.h"
#include "materials.h"
#include "sources.h"
#include "analysis.h"
#include "boundary.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static int clampi_local(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static double compute_full_field_max(const SimulationState* sim) {
    if (!sim || sim->nx <= 0 || sim->ny <= 0) return 0.0;
    double vmax = 0.0;
    for (int i = 0; i < sim->nx; ++i) {
        for (int j = 0; j < sim->ny; ++j) {
            double val = fabs(sim->Ez[i][j]);
            if (val > vmax) vmax = val;
        }
    }
    return vmax;
}

static double compute_scope_abs_max(const Scope* scope) {
    if (!scope || !scope->y || scope->n <= 0) return 0.0;
    double vmax = 0.0;
    for (int k = 0; k < scope->n; ++k) {
        double val = fabs(scope->y[k]);
        if (val > vmax) vmax = val;
    }
    return vmax;
}

static void ui_request_metrics_refresh(UIState* ui) {
    if (ui) {
        ui->force_metrics_recompute = 1;
    }
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
    ui->probe_x = 0;
    ui->probe_y = 0;
    ui->probe2_active = 0;
    ui->drag_src = -1;
    ui->sweep_on = 0;
    ui->sweep_idx = 0;
    ui->sweep_points = 0;
    ui->sweep_steps_per_point = 2000;
    ui->steps_per_frame = STEPS_PER_FRAME;
    ui->scale = RENDER_DEFAULT_SCALE;
    ui->ui_height = RENDER_DEFAULT_UI_HEIGHT;
    ui->side_panel_width = RENDER_DEFAULT_SIDE_PANEL;
    ui->log_probe = 0;

    ui->freq_slider.minv = 0.0;
    ui->freq_slider.maxv = 1.0;
    ui->freq_slider.value = 0.5;
    ui->freq_slider.dragging = 0;

    ui->speed_slider.minv = 1.0;
    ui->speed_slider.maxv = 50.0;
    ui->speed_slider.value = (double)ui->steps_per_frame;
    ui->speed_slider.dragging = 0;

    ui->render_field_vmax = 0.0;
    ui->render_scope_vmax = 0.0;
    ui->render_metrics_valid = 0;
    ui->estimated_field_vmax = 0.0;
    ui->estimated_scope_vmax = 0.0;
    ui->field_sample_timestep = -1;
    ui->scope_sample_generation = 0;
    ui->force_metrics_recompute = 0;
    ui->debug_force_metrics = 0;

    return ui;
}

void ui_state_free(UIState* state) {
    if (state) {
        free(state);
    }
}

void ui_state_set_layout(UIState* ui, int scale, int ui_height, int side_panel_width,
                         int nx, int ny) {
    if (!ui) return;
    ui->scale = scale;
    ui->ui_height = ui_height;
    ui->side_panel_width = side_panel_width;

    int canvas_height = ny * scale;
    int slider_width = nx * scale - 40;
    if (slider_width < 100) slider_width = 100;

    ui->freq_slider.x = 20;
    ui->freq_slider.y = canvas_height + 15;
    ui->freq_slider.w = slider_width;
    ui->freq_slider.h = 22;

    ui->speed_slider.x = 20;
    ui->speed_slider.y = canvas_height + 50;
    ui->speed_slider.w = slider_width;
    ui->speed_slider.h = 22;

    ui_request_metrics_refresh(ui);
}

void ui_state_sync_with_sim(UIState* ui, const SimulationState* sim) {
    if (!ui || !sim) return;
    ui->freq_slider.value = slider_from_freq(sim->freq, FREQ_MIN, FREQ_MAX);
    ui->speed_slider.value = (double)ui->steps_per_frame;
    ui->probe_x = sim->nx / 2;
    ui->probe_y = sim->ny / 2;
    ui->field_sample_timestep = sim->timestep;
}

void ui_update_metrics(UIState* ui, const SimulationState* sim, const Scope* scope) {
    if (!ui || !sim) return;

    const double decay = 0.995;

    if (ui->field_sample_timestep != sim->timestep) {
        double decayed = ui->estimated_field_vmax * decay;
        double sample = fabs(sim->step_Ez_absmax);
        if (sample > decayed) decayed = sample;
        ui->estimated_field_vmax = decayed;
        ui->field_sample_timestep = sim->timestep;
    } else {
        ui->estimated_field_vmax *= decay;
    }

    if (scope) {
        if (ui->scope_sample_generation != scope->rolling_generation) {
            double decayed = ui->estimated_scope_vmax * decay;
            double sample = fabs(scope->rolling_absmax);
            if (sample > decayed) decayed = sample;
            ui->estimated_scope_vmax = decayed;
            ui->scope_sample_generation = scope->rolling_generation;
        } else {
            ui->estimated_scope_vmax *= decay;
        }
    } else {
        ui->estimated_scope_vmax *= decay;
    }

    int force_full = ui->force_metrics_recompute || ui->debug_force_metrics;

    double vmax = 0.0;
    if (ui->render_metrics_valid && !force_full) {
        vmax = ui->render_field_vmax;
    } else if (!force_full && ui->estimated_field_vmax > 0.0) {
        vmax = ui->estimated_field_vmax;
    } else {
        vmax = compute_full_field_max(sim);
        ui->estimated_field_vmax = vmax;
    }

    if (ui->vmax_smooth == 0.0) ui->vmax_smooth = vmax;
    ui->vmax_smooth = 0.9 * ui->vmax_smooth + 0.1 * vmax;

    if (!ui->hold_color) {
        ui->held_vmax = ui->vmax_smooth;
    }

    double smax = 0.0;
    int have_scope_data = 0;
    if (ui->render_metrics_valid && !force_full) {
        smax = ui->render_scope_vmax;
        have_scope_data = 1;
    } else if (!force_full && ui->estimated_scope_vmax > 0.0) {
        smax = ui->estimated_scope_vmax;
        have_scope_data = 1;
    } else if (scope && scope->y && scope->n > 0) {
        smax = compute_scope_abs_max(scope);
        have_scope_data = 1;
        ui->estimated_scope_vmax = smax;
    }

    if (have_scope_data) {
        if (ui->scope_vmax_smooth == 0.0) ui->scope_vmax_smooth = smax;
        ui->scope_vmax_smooth = 0.9 * ui->scope_vmax_smooth + 0.1 * smax;
        if (!ui->hold_scope) {
            ui->held_scope_vmax = ui->scope_vmax_smooth;
        }
    }

    ui->render_metrics_valid = 0;
    if (ui->force_metrics_recompute && !ui->debug_force_metrics) {
        ui->force_metrics_recompute = 0;
    }
}

static void apply_paint(UIState* ui, SimulationState* sim, int mx, int my, int scale) {
    if (!ui->paint_mode) return;
    if (mx < 0 || my < 0) return;
    if (mx >= sim->nx * scale || my >= sim->ny * scale) return;
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
                    int max_x = sim->nx * scale;
                    int max_y = sim->ny * scale;
                    int mx = e.button.x;
                    int my = e.button.y;
                    if (mx < max_x && my < max_y) {
                        /* Prefer dragging nearest active source; fall back to moving probe */
                        int idx = find_nearest_source(sim->sources, mx, my, scale, 10.0f * scale);
                        if (idx >= 0) {
                            ui->drag_src = idx;
                        } else {
                            ui->probe_x = clampi_local(mx / scale, 0, sim->nx - 1);
                            ui->probe_y = clampi_local(my / scale, 0, sim->ny - 1);
                        }
                    }
                }
            } else if (e.button.button == SDL_BUTTON_RIGHT) {
                ui->paint_mode = !ui->paint_mode;
            }
        } else if (e.type == SDL_MOUSEMOTION) {
            if (ui->paint_mode && (e.motion.state & SDL_BUTTON_LMASK)) {
                apply_paint(ui, sim, e.motion.x, e.motion.y, scale);
            } else if (!ui->paint_mode && ui->drag_src >= 0 && (e.motion.state & SDL_BUTTON_LMASK)) {
                int mx = e.motion.x;
                int my = e.motion.y;
                int ix = clampi_local(mx / scale, 1, sim->nx - 2);
                int iy = clampi_local(my / scale, 1, sim->ny - 2);
                sim->sources[ui->drag_src].ix = ix;
                sim->sources[ui->drag_src].iy = iy;
            }
        } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
            ui->drag_src = -1;
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
                case SDLK_u:
                    ui->paint_mode = !ui->paint_mode;
                    break;
                case SDLK_i:
                    ui->paint_type = (ui->paint_type % 3) + 1;
                    break;
                case SDLK_o:
                    ui->paint_eps = clampd(ui->paint_eps * 0.9, 1.0, 20.0);
                    break;
                case SDLK_p:
                    ui->paint_eps = clampd(ui->paint_eps * 1.1, 1.0, 20.0);
                    break;

                /* Source controls: toggle individual sources and cycle type */
                case SDLK_1:
                    sim->sources[0].active = !sim->sources[0].active;
                    break;
                case SDLK_2:
                    if (MAX_SRC > 1) sim->sources[1].active = !sim->sources[1].active;
                    break;
                case SDLK_3:
                    if (MAX_SRC > 2) sim->sources[2].active = !sim->sources[2].active;
                    break;
                case SDLK_t:
                    sources_cycle_type(sim->sources);
                    fdtd_clear_fields(sim);
                    scope_clear(scope);
                    cpml_zero_psi(sim);
                    ui_request_metrics_refresh(ui);
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
                    fdtd_clear_fields(sim);
                    scope_clear(scope);
                    cpml_zero_psi(sim);
                    ui_request_metrics_refresh(ui);
                    break;
                case SDLK_r:
                    fdtd_reset(sim);
                    scope_clear(scope);
                    cpml_zero_psi(sim);
                    ui_request_metrics_refresh(ui);
                    break;
                case SDLK_y: {
                    BoundaryType type = boundary_get_type(sim);
                    type = (type == BOUNDARY_CPML) ? BOUNDARY_MUR : BOUNDARY_CPML;
                    boundary_set_type(sim, type);
                    if (boundary_is_cpml_enabled(sim)) {
                        cpml_build_coeffs(sim);
                        cpml_zero_psi(sim);
                    }
                    break;
                }
                case SDLK_g:
                    ui->log_probe = !ui->log_probe;
                    break;
                case SDLK_s:
                    sim->ports_on = !sim->ports_on;
                    break;
                case SDLK_f:
                    ui->debug_force_metrics = !ui->debug_force_metrics;
                    if (!ui->debug_force_metrics) {
                        ui_request_metrics_refresh(ui);
                    }
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
