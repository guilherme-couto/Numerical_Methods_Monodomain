#ifndef MV_SOLVER_H
#define MV_SOLVER_H

#include "../../../include/core_definitions.h"
#include "../cell_models.h"

// Constants for the MV model
#define MV_chi 1.0f
#define MV_Cm 1.0f

#define MV_NSV 3                     // Number of state variables in the MV model
#define MV_ACTIVATION_THRESHOLD 0.8f // Activation threshold for the MV model

// Options: ENDO, M, EPI, PB, TNNP -> default is ENDO
#if !defined(MCELL) && !defined(EPI) && !defined(ENDO) && !defined(PB) && !defined(TNNP)
#define ENDO
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    extern const CellModelSolver MV_SOLVER;

    static void initialize_MV(real *Vm, real *sV, const int Nx, const int Ny)
    {
        // Initial conditions
        const real u_init = 0.0f;
        const real v_init = 1.0f;
        const real w_init = 1.0f;
        const real s_init = 0.0f;

        // Initialize Vm and W with the initial condition
        for (int idx = 0, idx_sv = 0; idx < Nx * Ny; idx++, idx_sv = idx * MV_NSV)
        {
            Vm[idx] = u_init;
            sV[idx_sv] = v_init;
            sV[idx_sv + 1] = w_init;
            sV[idx_sv + 2] = s_init;
        }
    }

    // Rescale Vm -> from Minimal Ventricular paper
    static real rescaleVm(const real Vm)
    {
        return 85.7f * Vm - 84.0f;
    }

    #if defined(__CUDACC__)
    __device__ void d_get_actual_sV_MV(real *actualsV, const real *sV, const int idx);
    __device__ real d_compute_dVmdt_MV(const real Vm, const real *sV);
    __device__ void d_update_sVtilde_MV(real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t);
    __device__ void d_update_sV_MV(real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx);
    #endif 

    // Parameters - Based on Minimal Ventricular model
    // Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub

    // #define Dtilde = 0.65f * 1e-3; // cm^2/s

#ifdef EPI

#define u_o 0.0f
#define u_u 1.55f
#define theta_v 0.3f
#define theta_w 0.13f
#define theta_vminus 0.006f
#define theta_o 0.006f
#define tau_v1minus 60.0f
#define tau_v2minus 1150.0f
#define tau_vplus 1.4506f
#define tau_w1minus 60.0f
#define tau_w2minus 15.0f
#define k_wminus 65.0f
#define u_wminus 0.03f
#define tau_wplus 200.0f
#define tau_fi 0.11f
#define tau_o1 400.0f
#define tau_o2 6.0f
#define tau_so1 30.0181f
#define tau_so2 0.9957f
#define k_so 2.0458f
#define u_so 0.65f
#define tau_s1 2.7342f
#define tau_s2 16.0f
#define k_s 2.0994f
#define u_s 0.9087f
#define tau_si 1.8875f
#define tau_winf 0.07f
#define w_infstar 0.94f

#endif // EPI

#ifdef ENDO

#define u_o 0.0f
#define u_u 1.56f
#define theta_v 0.3f
#define theta_w 0.13f
#define theta_vminus 0.2f
#define theta_o 0.006f
#define tau_v1minus 75.0f
#define tau_v2minus 10.0f
#define tau_vplus 1.4506f
#define tau_w1minus 6.0f
#define tau_w2minus 140.0f
#define k_wminus 200.0f
#define u_wminus 0.016f
#define tau_wplus 280.0f
#define tau_fi 0.1f
#define tau_o1 470.0f
#define tau_o2 6.0f
#define tau_so1 40.0f
#define tau_so2 1.2f
#define k_so 2.0f
#define u_so 0.65f
#define tau_s1 2.7342f
#define tau_s2 2.0f
#define k_s 2.0994f
#define u_s 0.9087f
#define tau_si 2.9013f
#define tau_winf 0.0273f
#define w_infstar 0.78f

#endif // ENDO

#ifdef MCELL

#define u_o 0.0f
#define u_u 1.61f
#define theta_v 0.3f
#define theta_w 0.13f
#define theta_vminus 0.1f
#define theta_o 0.005f
#define tau_v1minus 80.0f
#define tau_v2minus 1.4506f
#define tau_vplus 1.4506f
#define tau_w1minus 70.0f
#define tau_w2minus 8.0f
#define k_wminus 200.0f
#define u_wminus 0.016f
#define tau_wplus 280.0f
#define tau_fi 0.078f
#define tau_o1 410.0f
#define tau_o2 7.0f
#define tau_so1 91.0f
#define tau_so2 0.8f
#define k_so 2.1f
#define u_so 0.6f
#define tau_s1 2.7342f
#define tau_s2 4.0f
#define k_s 2.0994f
#define u_s 0.9087f
#define tau_si 3.3849f
#define tau_winf 0.01f
#define w_infstar 0.5f

#endif // MCELL

#ifdef PB

#define u_o 0.0f
#define u_u 1.45f
#define theta_v 0.35f
#define theta_w 0.13f
#define theta_vminus 0.175f
#define theta_o 0.006f
#define tau_v1minus 10.0f
#define tau_v2minus 1150.0f
#define tau_vplus 1.4506f
#define tau_w1minus 140.0f
#define tau_w2minus 6.25f
#define k_wminus 65.0f
#define u_wminus 0.015f
#define tau_wplus 326.0f
#define tau_fi 0.105f
#define tau_o1 400.0f
#define tau_o2 6.0f
#define tau_so1 30.0181f
#define tau_so2 0.9957f
#define k_so 2.0458f
#define u_so 0.65f
#define tau_s1 2.7342f
#define tau_s2 16.0f
#define k_s 2.0994f
#define u_s 0.9087f
#define tau_si 1.8875f
#define tau_winf 0.175f
#define w_infstar 0.9f

#endif // PB

#ifdef TNNP

#define u_o 0.0f
#define u_u 1.58f
#define theta_v 0.3f
#define theta_w 0.015f
#define theta_vminus 0.015f
#define theta_o 0.006f
#define tau_v1minus 60.0f
#define tau_v2minus 1150.0f
#define tau_vplus 1.4506f
#define tau_w1minus 70.0f
#define tau_w2minus 20.0f
#define k_wminus 65.0f
#define u_wminus 0.03f
#define tau_wplus 280.0f
#define tau_fi 0.11f
#define tau_o1 6.0f
#define tau_o2 6.0f
#define tau_so1 43.0f
#define tau_so2 0.2f
#define k_so 2.0f
#define u_so 0.65f
#define tau_s1 2.7342f
#define tau_s2 3.0f
#define k_s 2.0994f
#define u_s 0.9087f
#define tau_si 2.8723f
#define tau_winf 0.07f
#define w_infstar 0.94f

#endif // TNNP

#ifdef __cplusplus
}
#endif

#endif // MV_SOLVER_H
