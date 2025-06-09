#ifndef AFHN_CUH
#define AFHN_CUH

#if defined __CUDACC__

#ifdef __cplusplus
extern "C"
{
#endif

#include "afhn.h"

// AFHN solver functions
__forceinline__ __device__ void d_get_actual_sV_AFHN(real *actualsV, const real *sV, const int idx)
{
    const int idx_sv = idx * AFHN_NSV;
    actualsV[0] = sV[idx_sv];
}

__forceinline__ __device__ real d_compute_dVmdt_AFHN(const real Vm, const real *sV)
{
    const real W = sV[0];
    return (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
}

__forceinline__ __device__ void d_update_sVtilde_AFHN(real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t)
{
    const real W = rhs_sV[0];
    const real dWdt = (eta2 * ((Vm / vp) - (eta3 * W)));

    sVtilde[0] = W + delta_t * dWdt;
}

__forceinline__ __device__ void d_update_sV_AFHN(real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    const real W = rhs_sV[0];
    const real dSdt_W = dSdt_sV[0];
    const real dWdt = (eta2 * ((dSdt_Vm / vp) - (eta3 * dSdt_W)));

    const int idx_sv = idx * AFHN_NSV;
    sV[idx_sv] = W + delta_t * dWdt;
}

#ifdef __cplusplus
}
#endif

#endif // __CUDACC__
#endif // AFHN_CUH

