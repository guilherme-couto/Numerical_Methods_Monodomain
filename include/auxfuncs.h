#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "core_definitions.h"
#include "config_parser.h"
#include "../src/logger/logger.h"

#ifdef __cplusplus
extern "C" {
#endif

void populateDiagonalThomasAlgorithm(real *la, real *lb, real *lc, int N, real phi);
void prefactorizeThomas(const real *la, const real *lb, const real *lc, real *c_prime, real *denominator, const int sysSize);
int createDirectories(char *dir_path, bool remove_old_files);
void initializeTimeArray(real *timeArray, int M, real dt);
void initializeMeasurement(Measurement *measurement);
int populateStimuli(SimulationConfig *config);
void saveCopyOfSimulationConfig(const char *ini_file_path, const char *output_dir);
void tridiagonalSystemSolver_x(const int Nx, real *rhs, real *result, real *c_prime, real *d_prime, const real phi_x, const ElementProperties *elements, const int sys_idx);
void tridiagonalSystemSolver_y(const int Ny, real *rhs, real *result, real *c_prime, real *d_prime, const real phi_y, const ElementProperties *elements, const int sys_idx, const int Nx);
int saveSimulationInfos(const SimulationConfig *config, const Measurement *measurement);
real calculateNorm2Error(real *Vm, real *exact, int Nx, int Ny, real totalTime, real delta_x, real delta_y, real Lx, real Ly);
void initializeElementsProperties(const SimulationConfig *config, ElementProperties *elements_properties);
void initializeVariableWithExactSolution(real *Var, int Nx, int Ny, real delta_x, real delta_y, real Lx, real Ly);
void initializeVariableWithValue(real *Var, int Nx, int Ny, real value);
void initializeVariableFromFile(real *Var, int Nx, int Ny, char *filename, real delta_x, real delta_y, char *varName, real reference_dx, real reference_dy, real Lx, real Ly);
void shiftVariableToLeft(real *Var, int Nx, int Ny, real length, real delta_x, real delta_y, real initValue, char *varName);

#ifdef __cplusplus
}
#endif

#ifdef CABLEEQ

void initialize1DVariableWithValue(real *Var, int N, real value);
void initialize1DVariableFromFile(real *Var, int N, char *filename, real delta_x, char *varName, real reference_dx, real Lx);
void shift1DVariableToLeft(real *Var, int N, real length, real delta_x, real initValue, char *varName);

#endif // CABLEEQ

#endif // AUXFUNCS_H