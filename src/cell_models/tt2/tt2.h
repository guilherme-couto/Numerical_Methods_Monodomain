#ifndef TT2_SOLVER_H
#define TT2_SOLVER_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"
#include "../../logger/logger.h"

// Options: ENDO, MCELL, EPI -> default is ENDO
#if !defined(MCELL) && !defined(EPI) && !defined(ENDO)
#define ENDO
#endif // MCELL and EPI are mutually exclusive

#ifdef __cplusplus
extern "C"
{
#endif

    void solveMonodomainTT2(const SimulationConfig *config, Measurement *measurement, const real *time_array);

    // Model parameters - Based on Ten Tusscher 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
    // from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
    // and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D

    // #define chi 1400.0f  // Surface-to-volume ratio (cm^-1)
    // #define Cm 0.185f    // Membrane capacitance (uF/cm^2)
    // #define sigma 1.171f // Diffusion coefficient (cm^2/s)

#define R 8314.472f      // Universal gas constant (J/(kmol*K))
#define T 310.0f         // Temperature (K)
#define F 96485.3415f    // Faraday's constant (C/mol)
#define RTONF 26.713761f // R*T/F
#define FONRT 0.037434f  // F/(R*T)

// Intracellular volumes
#define V_C 0.016404f    // Cellular volume -> (???) [16404 um^3]
#define V_SR 0.001094f   // Sarcoplasmic reticulum volume -> (???) [1094 um^3]
#define V_SS 0.00005468f // Subsarcolemmal space volume -> (???) [54.68 um^3]

// External concentrations
#define K_o 5.4f    // Extracellular potassium (K+) concentration -> mM
#define Na_o 140.0f // Extracellular sodium (Na+) concentration -> mM
#define Ca_o 2.0f   // Extracellular calcium (Ca++) concentration -> mM

// Parameters for currents
#define G_Na 14.838f // Maximal I_Na (sodium current) conductance -> nS/pF
#define G_K1 5.405f  // Maximal I_K1 (late rectifier potassium current) conductance -> nS/pF

#if defined(EPI) || defined(MCELL)

#define G_to 0.294f // Maximal I_to (transient outward potassium current) conductance -> nS/pF (epi and M cells)

#endif // EPI || MCELL

#ifdef ENDO

#define G_to 0.073f // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)

#endif // ENDO

#define G_Kr 0.153f // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF

#if defined(EPI) || defined(ENDO)

#define G_Ks 0.392f // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)

#endif // EPI || ENDO

#ifdef MCELL

#define G_Ks 0.098f // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells)

#endif // MCELL

#define p_KNa 0.03f        // Relative I_Ks permeability to Na+ over K+ -> dimensionless
#define G_CaL 3.98e-5f     // Maximal I_CaL (L-type calcium current) conductance -> cm/ms/uF
#define k_NaCa 1000.0f     // Maximal I_NaCa (Na+/Ca++ exchanger current) -> pA/pF
#define gamma_I_NaCa 0.35f // Voltage dependence parameter of I_NaCa -> dimensionless
#define K_mCa 1.38f        // Half-saturation constant of I_NaCa for intracellular Ca++ -> mM
#define K_mNa_i 87.5f      // Half-saturation constant of I_NaCa for intracellular Na+ -> mM
#define k_sat 0.1f         // Saturation factor for I_NaCa -> dimensionless
#define alpha 2.5f         // Factor enhancing outward nature of I_NaCa -> dimensionless
#define P_NaK 2.724f       // Maximal I_NaK (Na+/K+ pump current) -> pA/pF
#define K_mK 1.0f          // Half-saturation constant of I_NaK for Ko -> mM
#define K_mNa 40.0f        // Half-saturation constant of I_NaK for intracellular Na+ -> mM
#define G_pK 0.0146f       // Maximal I_pK (plateau potassium current) conductance -> nS/pF
#define G_pCa 0.1238f      // Maximal I_pCa (plateau calcium current) conductance -> nS/pF
#define K_pCa 0.0005f      // Half-saturation constant of I_pCa for intracellular Ca++ -> mM
#define G_bNa 0.00029f     // Maximal I_bNa (sodium background current) conductance -> nS/pF
#define G_bCa 0.000592f    // Maximal I_bCa (calcium background current) conductance -> nS/pF

// Intracellular calcium flux dynamics
#define V_maxup 0.006375f // Maximal I_up -> mM/ms
#define K_up 0.00025f     // Half-saturation constant of I_up -> mM
#define V_rel 0.102f      // Maximal I_rel conductance -> mM/ms
#define k1_prime 0.15f    // R to O and RI to I I_rel transition rate -> mM^-2*ms^-1
#define k2_prime 0.045f   // O to I  and R to RI I_rel transition rate -> mM^-1*ms^-1
#define k3 0.06f          // O to R and I to RI I_rel transition rate -> ms^-1
#define k4 0.005f         // I to O and RI to I I_rel transition rate -> ms^-1
#define EC 1.5f           // Half-saturation constant of k_Ca_SR -> mM
#define max_SR 2.5f       // Maximum value of k_Ca_SR -> dimensionless
#define min_SR 1.0f       // Minimum value of k_Ca_SR -> dimensionless
#define V_leak 0.00036f   // Maximal I_leak conductance -> mM/ms
#define V_xfer 0.0038f    // Maximal I_xfer conductance -> mM/ms

// Calcium buffering dynamics
#define bufC 0.2f        // Total cytoplasmic buffer concentration -> mM
#define K_bufC 0.001f    // Half-saturation constant of cytoplasmic buffers -> mM
#define bufSR 10.0f      // Total sarcoplasmic reticulum buffer concentration -> mM
#define K_bufSR 0.3f     // Half-saturation constant of sarcoplasmic reticulum buffers -> mM
#define bufSS 0.4f       // Total subspace buffer concentration -> mM
#define K_bufSS 0.00025f // Half-saturation constant of subspace buffer -> mM

    // Initial conditions

#ifdef EPI

    // Initial conditions for EPI cells
    // from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/
    // const real Vm_init -85.23f      // Initial membrane potential -> mV
    // const real X_r1_init 0.00621f   // Initial rapid time-dependent potassium current Xr1 gate -> dimensionless
    // const real X_r2_init 0.4712f    // Initial rapid time-dependent potassium current Xr2 gate -> dimensionless
    // const real X_s_init 0.0095f     // Initial slow time-dependent potassium current Xs gate -> dimensionless
    // const real m_init 0.00172f      // Initial fast sodium current m gate -> dimensionless
    // const real h_init 0.7444f       // Initial fast sodium current h gate -> dimensionless
    // const real j_init 0.7045f       // Initial fast sodium current j gate -> dimensionless
    // const real d_init 3.373e-5f     // Initial L-type calcium current d gate -> dimensionless
    // const real f_init 0.7888f       // Initial L-type calcium current f gate -> dimensionless
    // const real f2_init 0.9755f      // Initial L-type calcium current f2 gate -> dimensionless
    // const real fCaSS_init 0.9953f   // Initial L-type calcium current fCaSS gate -> dimensionless
    // const real s_init 0.999998f     // Initial transient outward current s gate -> dimensionless
    // const real r_init 2.42e-8f      // Initial transient outward current r gate -> dimensionless
    // const real Ca_i_init 0.000126f  // Initial intracellular Ca++ concentration -> mM
    // const real Ca_SR_init 3.64f     // Initial sarcoplasmic reticulum Ca++ concentration -> mM
    // const real Ca_SS_init 0.00036f  // Initial subspace Ca++ concentration -> mM
    // const real R_prime_init 0.9073f // Initial ryanodine receptor -> dimensionless
    // const real Na_i_init 8.604f     // Initial intracellular Na+ concentration -> mM
    // const real K_i_init 136.89f     // Initial intracellular K+ concentration -> mM

#endif // EPI

#if defined(ENDO) || defined(MCELL)

    // Initial conditions for ENDO and MCELL cells
    // from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/
    // const real Vm_init -86.2f      // Initial membrane potential -> mV
    // const real X_r1_init 0.0f      // Initial rapid time-dependent potassium current Xr1 gate -> dimensionless
    // const real X_r2_init 1.0f      // Initial rapid time-dependent potassium current Xr2 gate -> dimensionless
    // const real X_s_init 0.0f       // Initial slow time-dependent potassium current Xs gate -> dimensionless
    // const real m_init 0.0f         // Initial fast sodium current m gate -> dimensionless
    // const real h_init 0.75f        // Initial fast sodium current h gate -> dimensionless
    // const real j_init 0.75f        // Initial fast sodium current j gate -> dimensionless
    // const real d_init 0.0f         // Initial L-type calcium current d gate -> dimensionless
    // const real f_init 1.0f         // Initial L-type calcium current f gate -> dimensionless
    // const real f2_init 1.0f        // Initial L-type calcium current f2 gate -> dimensionless
    // const real fCaSS_init 1.0f     // Initial L-type calcium current fCaSS gate -> dimensionless
    // const real s_init 1.0f         // Initial transient outward current s gate -> dimensionless
    // const real r_init 0.0f         // Initial transient outward current r gate -> dimensionless
    // const real Ca_i_init 0.00007f  // Initial intracellular Ca++ concentration -> mM
    // const real Ca_SR_init 1.3f     // Initial sarcoplasmic reticulum Ca++ concentration -> mM
    // const real Ca_SS_init 0.00007f // Initial subspace Ca++ concentration -> mM
    // const real R_prime_init 1.0f   // Initial ryanodine receptor -> dimensionless
    // const real Na_i_init 7.67f     // Initial intracellular Na+ concentration -> mM
    // const real K_i_init 138.3f     // Initial intracellular K+ concentration -> mM

#endif // ENDO || MCELL

#ifdef __cplusplus
}
#endif

#endif // TT2_SOLVER_H