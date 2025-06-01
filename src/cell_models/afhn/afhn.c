#include "afhn.h"

// AFHN solver functions
static void get_actual_sV_AFHN(real *actualsV, const real *sV, const int idx)
{
    const int idx_sv = idx * AFHN_NSV;
    actualsV[0] = sV[idx_sv];
}

static const real compute_dVmdt_AFHN(const real Vm, const real *sV)
{
    const real W = sV[0];
    return (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
}

static void update_sVtilde_AFHN(real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t)
{
    const real W = rhs_sV[0];
    const real dWdt = (eta2 * ((Vm / vp) - (eta3 * W)));

    sVtilde[0] = W + delta_t * dWdt;
}

static void update_sV_AFHN(real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    const real W = rhs_sV[0];
    const real dSdt_W = dSdt_sV[0];
    const real dWdt = (eta2 * ((dSdt_Vm / vp) - (eta3 * dSdt_W)));

    const int idx_sv = idx * AFHN_NSV;
    sV[idx_sv] = W + delta_t * dWdt;
}

// Instantiate the AFHN model
const CellModelSolver AFHN_SOLVER = {
    .n_state_vars = AFHN_NSV,
    .denom_chiCm = 1.0f / (AFHN_chi * AFHN_Cm),
    .activation_thershold = AFHN_ACTIVATION_THRESHOLD,
    .initialize = initialize_AFHN,
    .get_actual_sV = get_actual_sV_AFHN,
    .compute_dVmdt = compute_dVmdt_AFHN,
    .update_sVtilde = update_sVtilde_AFHN,
    .update_sV = update_sV_AFHN,
};

// real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma)
// {
//     // Calculate coefficients for the ADI method
//     const real diff_coeff = sigma / (AFHN_Cm * AFHN_chi);
//     real exactVm = (exp(-t)) * cos(_PI * x / Lx) * cos(_PI * y / Ly);
//     real reaction = dVmdt(exactVm, W);
//     return (exactVm * (-(AFHN_chi * AFHN_Cm) + 2.0f * (sigma / (AFHN_chi * AFHN_Cm)) * _PI * _PI / (Lx * Ly))) + (AFHN_chi * reaction);
// }