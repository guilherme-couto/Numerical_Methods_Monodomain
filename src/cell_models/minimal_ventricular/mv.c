#include "mv.h"

// MV solver functions
static void get_actual_sV_MV(real *actualsV, const real *sV, const int idx)
{
    const int idx_sv = idx * MV_NSV;
    for (int i = 0; i < MV_NSV; ++i)
    {
        actualsV[i] = sV[idx_sv + i];
    }
}

static real compute_dVmdt_MV(const real Vm, const real *sV)
{
    // Extract state variables
    const real v = sV[0];
    const real w = sV[1];
    const real s = sV[2];

    // Auxiliary variables
    const real Vm_theta_w = Vm - theta_w;
    const real Vm_theta_o = Vm - theta_o;
    const real Vm_theta_v = Vm - theta_v;

    const bool Htheta_w = Vm_theta_w > 0.0f;
    const bool Htheta_o = Vm_theta_o > 0.0f;
    const bool Htheta_v = Vm_theta_v > 0.0f;

    const real tau_o = Htheta_o ? tau_o2 : tau_o1;
    const real tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vm - u_so))) * 0.5f;

    // Currents
    const real J_fi = -v * Htheta_v * (Vm_theta_v) * (u_u - Vm) / tau_fi;
    const real J_so = ((Vm - u_o) * (1.0f - Htheta_w) / tau_o) + (Htheta_w / tau_so);
    const real J_si = -Htheta_w * w * s / tau_si;

    // Compute the derivative of Vm
    return J_fi + J_so + J_si;
}

static void update_sVtilde_MV(real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t)
{
    // Extract state variables
    const real v = rhs_sV[0];
    const real w = rhs_sV[1];
    const real s = rhs_sV[2];

    // Auxiliary variables
    const real Vm_theta_w = Vm - theta_w;
    const real Vm_theta_o = Vm - theta_o;
    const real Vm_theta_v = Vm - theta_v;
    const real Vm_theta_vminus = Vm - theta_vminus;

    const bool Htheta_w = Vm_theta_w > 0.0f;
    const bool Htheta_o = Vm_theta_o > 0.0f;
    const bool Htheta_v = Vm_theta_v > 0.0f;
    const bool Htheta_vminus = Vm_theta_vminus > 0.0f;

    const real tau_vminus = Htheta_vminus ? tau_v2minus : tau_v1minus;

    // Auxiliary variables for Rush-Larsen or Forward Euler
    const real v_inf = !Htheta_vminus;
    const real commom_denominator_v = 1.0f / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
    const real tau_v_RL = tau_vplus * tau_vminus * commom_denominator_v;
    const real v_inf_RL = tau_vplus * v_inf * !Htheta_v * commom_denominator_v;

    const real w_inf = !Htheta_o * (1.0f - (Vm / tau_winf)) + Htheta_o * w_infstar;
    const real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (Vm - u_wminus))) * 0.5f;
    const real commom_denominator_w = 1.0f / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
    const real tau_w_RL = tau_wplus * tau_wminus * commom_denominator_w;
    const real w_inf_RL = tau_wplus * w_inf * !Htheta_w * commom_denominator_w;

    const real tau_s = !Htheta_w * tau_s1 + Htheta_w * tau_s2;
    const real s_inf_RL = (1.0f + tanh(k_s * (Vm - u_s))) * 0.5f;

    // Calculate approximations with Rush-Larsen or Forward Euler using half time step
    sVtilde[0] = (tau_v_RL > 1.0e-10f)
                     ? v_inf_RL - (v_inf_RL - v) * exp(-delta_t / tau_v_RL)
                     : v + delta_t * (!Htheta_v * (v_inf - v) / tau_vminus - Htheta_v * v / tau_vplus);

    sVtilde[1] = (tau_w_RL > 1.0e-10f)
                     ? w_inf_RL - (w_inf_RL - w) * exp(-delta_t / tau_w_RL)
                     : w + delta_t * (!Htheta_w * (w_inf - w) / tau_wminus - Htheta_w * w / tau_wplus);

    sVtilde[2] = (tau_s > 1.0e-10f)
                     ? s_inf_RL - (s_inf_RL - s) * exp(-delta_t / tau_s)
                     : s + delta_t * (s_inf_RL - s) / tau_s;
}

static void update_sV_MV(real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    // Create a local copy of the state variables to avoid aliasing
    const real local_rhs[MV_NSV] = {rhs_sV[0], rhs_sV[1], rhs_sV[2]};
    const real local_dSdt[MV_NSV] = {dSdt_sV[0], dSdt_sV[1], dSdt_sV[2]};

    // Extract state variables
    const real dSdt_v = local_dSdt[0];
    const real dSdt_w = local_dSdt[1];
    const real dSdt_s = local_dSdt[2];
    const real v = local_rhs[0];
    const real w = local_rhs[1];
    const real s = local_rhs[2];

    // Auxiliary variables
    const real Vm_theta_w = dSdt_Vm - theta_w;
    const real Vm_theta_o = dSdt_Vm - theta_o;
    const real Vm_theta_v = dSdt_Vm - theta_v;
    const real Vm_theta_vminus = dSdt_Vm - theta_vminus;

    const bool Htheta_w = Vm_theta_w > 0.0f;
    const bool Htheta_o = Vm_theta_o > 0.0f;
    const bool Htheta_v = Vm_theta_v > 0.0f;
    const bool Htheta_vminus = Vm_theta_vminus > 0.0f;

    const real tau_vminus = Htheta_vminus ? tau_v2minus : tau_v1minus;

    // Auxiliary variables for Rush-Larsen or Forward Euler
    const real v_inf = !Htheta_vminus;
    const real commom_denominator_v = 1.0f / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
    const real tau_v_RL = tau_vplus * tau_vminus * commom_denominator_v;
    const real v_inf_RL = tau_vplus * v_inf * !Htheta_v * commom_denominator_v;

    const real w_inf = !Htheta_o * (1.0f - (dSdt_Vm / tau_winf)) + Htheta_o * w_infstar;
    const real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (dSdt_Vm - u_wminus))) * 0.5f;
    const real commom_denominator_w = 1.0f / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
    const real tau_w_RL = (tau_wplus * tau_wminus) * commom_denominator_w;
    const real w_inf_RL = (tau_wplus * w_inf * !Htheta_w) * commom_denominator_w;

    const real tau_s = !Htheta_w * tau_s1 + Htheta_w * tau_s2;
    const real s_inf_RL = (1.0f + tanh(k_s * (dSdt_Vm - u_s))) * 0.5f;

    // Calculate approximations with Rush-Larsen or Forward Euler using half time step
    const int idx_sv = idx * MV_NSV;
    sV[idx_sv] = (tau_v_RL > 1.0e-10f)
                     ? v_inf_RL - (v_inf_RL - v) * exp(-delta_t / tau_v_RL)
                     : v + delta_t * (!Htheta_v * (v_inf - dSdt_v) / tau_vminus - Htheta_v * dSdt_v / tau_vplus);

    sV[idx_sv + 1] = (tau_w_RL > 1.0e-10f)
                         ? w_inf_RL - (w_inf_RL - w) * exp(-delta_t / tau_w_RL)
                         : w + delta_t * (!Htheta_w * (w_inf - dSdt_w) / tau_wminus - Htheta_w * dSdt_w / tau_wplus);

    sV[idx_sv + 2] = (tau_s > 1.0e-10f)
                         ? s_inf_RL - (s_inf_RL - s) * exp(-delta_t / tau_s)
                         : s + delta_t * (s_inf_RL - dSdt_s) / tau_s;
}

// Instantiate the MV model
const CellModelSolver MV_SOLVER = {
    .n_state_vars = MV_NSV,
    .chiCm = MV_chi * MV_Cm,
    .activation_thershold = MV_ACTIVATION_THRESHOLD,
    .initialize = initialize_MV,
    .get_actual_sV = get_actual_sV_MV,
    .compute_dVmdt = compute_dVmdt_MV,
    .update_sVtilde = update_sVtilde_MV,
    .update_sV = update_sV_MV,
};