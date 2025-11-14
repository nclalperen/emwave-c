// =============================================================================
// emwave-c: Simulation configuration loader
// =============================================================================

#ifndef EMWAVE_CONFIG_LOADER_H
#define EMWAVE_CONFIG_LOADER_H

#include <stddef.h>
#include "config.h"

#ifndef CONFIG_LOADER_MAX_FILE_BYTES
#define CONFIG_LOADER_MAX_FILE_BYTES (10 * 1024 * 1024)
#endif

int config_loader_parse_file(const char* path, SimulationConfig* out_config,
                             char* errbuf, size_t errbuf_len);
int config_load_from_args(int argc, char** argv, SimulationConfig* out_config);
void config_print_summary(const SimulationConfig* cfg);
int config_validate(const SimulationConfig* cfg, char* errbuf, size_t errbuf_len);
void config_clamp_to_limits(SimulationConfig* cfg);

#endif /* EMWAVE_CONFIG_LOADER_H */
