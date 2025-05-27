#ifndef CELL_MODELS_H
#define CELL_MODELS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"

#ifdef __cplusplus
extern "C"
{
#endif

    // Types definition for the cell model solvers
    typedef void (*initialize_t)(real *, real *, const int, const int);
    typedef void (*get_actual_sV_t)(real *, const real *, const int);
    typedef real (*compute_dVmdt_t)(const real, const real *);
    typedef void (*update_sVtilde_t)(real *, const real, const real *, const real);
    typedef void (*update_sV_t)(real *, const real *, const real, const real *, const real, const int);

    // Structure to hold the cell model solver information
    typedef struct
    {
        int n_state_vars;
        real chiCm;
        real activation_thershold;
        initialize_t initialize;
        get_actual_sV_t get_actual_sV;
        compute_dVmdt_t compute_dVmdt;
        update_sVtilde_t update_sVtilde;
        update_sV_t update_sV;
    } CellModelSolver;

    // Instantiate the cell model solvers
    extern const CellModelSolver AFHN_SOLVER;
    extern const CellModelSolver TT2_SOLVER;
    extern const CellModelSolver MV_SOLVER;

    // Get the cell model solver based on the cell model type
    const CellModelSolver *get_solver_struct(const CellModel *cell_model);

    #if defined(__CUDACC__)
    #define MAX_NSV 16
    __device__ void select_get_actual_sV(const CellModel cell_model, real *actualsV, const real *sV, const int idx);
    __device__ real select_compute_dVmdt(const CellModel cell_model, const real Vm, const real *sV);
    __device__ void select_update_sVtilde(const CellModel cell_model, real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t);
    __device__ void select_update_sV(const CellModel cell_model, real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx);
    #endif // __CUDACC__

#ifdef __cplusplus
}
#endif

#endif // CELL_MODELS_H