#ifndef MONODOMAIN_SOLVER_H
#define MONODOMAIN_SOLVER_H

#include "../../include/config_parser.h"
#include "../../include/core_definitions.h"
#include "../logger/logger.h"
#include "../../include/auxfuncs.h"
#include "../numerical_methods/numerical_methods.h"

#ifdef __cplusplus
extern "C" {
#endif

int runMonodomainSimulationSerial(const SimulationConfig *config);
int runMonodomainSimulationOpenMP(const SimulationConfig *config);
int runMonodomainSimulationCUDA(const SimulationConfig *config);

#ifdef __cplusplus
}
#endif

#endif // MONODOMAIN_SOLVER_H