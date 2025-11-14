// ============================================================================
// emwave-c: Headless CLI sweep runner
// ============================================================================

#include "config.h"
#include "types.h"
#include "fdtd_core.h"
#include "boundary.h"
#include "materials.h"
#include "analysis.h"
#include "sources.h"

#include <stdio.h>

static void setup_ports(SimulationState* sim) {
    sim->ports_on = 1;
    for (int p = 0; p < MAX_PORTS; ++p) {
        sim->ports[p].active = 1;
        sim->ports[p].head = 0;
    }
}

static void init_world(SimulationState* sim) {
    materials_init(sim);
    ports_init(sim->ports, sim->nx, sim->ny);
    setup_ports(sim);
    boundary_set_type(sim, BOUNDARY_CPML);
    cpml_apply_preset(sim, cpml_get_preset_index(sim));
    cpml_zero_psi(sim);
}

int main(void) {
    SimulationState* sim = fdtd_init();
    if (!sim) {
        fprintf(stderr, "Failed to initialize simulation\n");
        return 1;
    }

    init_world(sim);

    const int sweep_points = 8;
    const double f_start = 5e8;
    const double f_stop = 3e9;

    printf("# freq_Hz,s21_mag\n");
    for (int idx = 0; idx < sweep_points; ++idx) {
        double t = (double)idx / (double)(sweep_points - 1);
        double freq = f_start + t * (f_stop - f_start);

        fdtd_reset(sim);
        fdtd_update_grid_for_freq(sim, freq);
        cpml_build_coeffs(sim);
        cpml_zero_psi(sim);

        int total_steps = 4000;
        for (int step = 0; step < total_steps; ++step) {
            fdtd_step(sim);
            if (sim->ports_on) {
                ports_sample(sim, sim->dx, sim->dy);
            }
        }

        double s21 = compute_s21(sim->ports, freq, sim->dt);
        printf("%.9e,%.6e\n", freq, s21);
    }

    ports_free(sim->ports);
    fdtd_free(sim);
    return 0;
}
