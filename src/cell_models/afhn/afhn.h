#ifndef AFHN_SOLVER_H
#define AFHN_SOLVER_H

#include "../../../include/core_definitions.h"
#include "../cell_models.h"

// Constants for the AFHN model
#define AFHN_chi 1.0e3f // cm^-1
#define AFHN_Cm 1.0e-3f // mF * cm^-2

#define AFHN_NSV 1
#define AFHN_ACTIVATION_THRESHOLD 10.0f

#ifdef __cplusplus
extern "C"
{
#endif

    extern const CellModelSolver AFHN_SOLVER;

    static void initialize_AFHN(real *Vm, real *sV, const int Nx, const int Ny)
    {
        // Initial conditions
        const real Vm_init = 0.0f;
        const real W_init = 0.0f;

        // Initialize Vm and W with the initial condition
        for (int idx = 0, idx_sv = 0; idx < Nx * Ny; idx++, idx_sv = idx)
        {
            Vm[idx] = Vm_init;
            sV[idx_sv] = W_init;
        }
    }

    // Model parameters - Based on Gerardo_Giorda 2007
    // #define sigma 1.2e-3f // omega^-1 * cm^-1
    // #define chi 1.0e3f    // cm^-1
    // #define Cm 1.0e-3f    // mF * cm^-2

#define G 1.5f      // omega^-1 * cm^-2
#define eta1 4.4f   // omega^-1 * cm^-1
#define eta2 0.012f // dimensionless
#define eta3 1.0f   // dimensionless
#define vth 13.0f   // mV
#define vp 100.0f   // mV

#ifdef __cplusplus
}
#endif

#endif // AFHN_SOLVER_H