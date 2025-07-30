#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "numerical_methods_helpers.h"
#include "../cell_models/cell_models_headers.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef void (*numerical_method_t)(const SimulationConfig *, Measurement *, const real *, const CellModelSolver *, real *, real *, const real *, const real *, const real *);

// Serial methods
void runSSIADI(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runOSADI(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runFE(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runDO(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runHV(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runCG(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);

// OpenMP methods
void runSSIADI_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runOSADI_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runFE_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runDO_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runHV_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runCG_OMP(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);

// CUDA methods
void runSSIADI_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runOSADI_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runFE_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runDO_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runDG_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runDY_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runHV_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runCS_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);
void runCG_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver_CUDA, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy);

// Function to get the appropriate numerical method
numerical_method_t get_numerical_method(const ExecutionMode *exec_mode, const NumericalMethod *method);

#ifdef __cplusplus
}
#endif

#endif // NUMERICAL_METHODS_H