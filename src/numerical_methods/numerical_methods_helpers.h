#ifndef NUMERICAL_METHODS_HELPERS_H
#define NUMERICAL_METHODS_HELPERS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../logger/logger.h"
#include "../cell_models/cell_models_headers.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define HARMONIC_MEAN(a, b) (2.0f * (a) * (b) / ((a) + (b) + 1e-12f))
#define ARITHMETIC_MEAN(a, b) (((a) + (b)) * 0.5f)

// Function to get the stimulus value at a given time and position
MODIFIERS real get_stimulus_value(const real actualTime, const int i, const int j, const Stimulus *stimuli, const int numberOfStimuli)
{
    #pragma unroll
    for (int k = 0; k < numberOfStimuli; k++)
    {
        const Stimulus *s = &stimuli[k];

        if (actualTime < s->start_time || actualTime > s->start_time + s->duration)
            continue;

        if (j < s->x_discretized.min || j > s->x_discretized.max)
            continue;

        if (i < s->y_discretized.min || i > s->y_discretized.max)
            continue;

        return s->amplitude;
    }
    return 0.0f;
}

// Function to get the active stimuli
MODIFIERS int update_and_get_num_active_stimuli(const real actualTime, const Stimulus *stimuli, const int numberOfStimuli, Stimulus *active_stimuli)
{
    int num_active_stimuli = 0;
    for (int k = 0; k < numberOfStimuli; k++)
    {
        const Stimulus *s = &stimuli[k];
        if (actualTime >= s->start_time && actualTime <= s->start_time + s->duration)
            active_stimuli[num_active_stimuli++] = *s;
    }
    return num_active_stimuli;
}

// Function to handle the saving of frames
static inline void handle_frame_saving(const char *pathToSaveData, const char *file_extension,
                                       const save_function_t save_function, const int timeStepCounter,
                                       const real *Vm, const int Nx, const int Ny,
                                       const real delta_x, const real delta_y, const real actualTime)
{
    static char file_path[MAX_STRING_SIZE];
    snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
    save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
}

// Function to handle the velocity measurement
static inline void handle_velocity_measurement(const real Vm_x0, const real Vm_x1, real *t0, real *t1, const real thereshold,
                                               bool *aux_stim_velocity_flag, bool *stim_velocity_measured,
                                               const real actualTime, const real x0, const real x1, real *stim_velocity)
{
    if (!*aux_stim_velocity_flag)
    {
        if (Vm_x0 > thereshold)
        {
            *t0 = actualTime;
            *aux_stim_velocity_flag = true;
        }
    }
    else
    {
        if (Vm_x1 > thereshold)
        {
            *t1 = actualTime;
            *stim_velocity = ((x1 - x0) / (*t1 - *t0)) * 1000.0f; // cm/ms
            *stim_velocity_measured = true;
            INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, *stim_velocity);
        }
    }
}

MODIFIERS real compute_diffusion_term_x_no_rotation(const real *Vm, const real *Dxx,
                                                    const int i, const int j, const int Nx, const int Ny,
                                                    const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left  = (j > 0);
    const bool right = (j < Nx - 1);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];

    // Diffusion coefficients at center
    const real Dxx_c = Dxx[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_l = left  ? HARMONIC_MEAN(Dxx_c, Dxx[idx_l])  : Dxx_c;
    const real Dxx_r = right ? HARMONIC_MEAN(Dxx_c, Dxx[idx_r])  : Dxx_c;

    const real denom_delta_x = 1.0f / delta_x;

    // Fluxes
    const real J_l = left  * ((Dxx_l * denom_delta_x) * (Vm_c - Vm_l));
    const real J_r = right * ((Dxx_r * denom_delta_x) * (Vm_r - Vm_c));

    // Divergence of flux -> nabla dot J
    return (J_r - J_l) * denom_delta_x;
}

MODIFIERS real compute_diffusion_term_x_with_rotation(const real *Vm, const real *Dxx, const real *Dxy,
                                                      const int i, const int j, const int Nx, const int Ny,
                                                      const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Diagonal neighbors with mirroring
    const int idx_ne = (top && right) ? idx + Nx + 1 : 
                      (top && !right) ? idx + Nx :
                      (!top && right) ? idx + 1 :
                      idx - Nx - 1;
    
    const int idx_nw = (top && left) ? idx + Nx - 1 :
                      (top && !left) ? idx + Nx :
                      (!top && left) ? idx - 1 :
                      idx - Nx + 1;
    
    const int idx_se = (bottom && right) ? idx - Nx + 1 :
                      (bottom && !right) ? idx - Nx :
                      (!bottom && right) ? idx + 1 :
                      idx + Nx - 1;
    
    const int idx_sw = (bottom && left) ? idx - Nx - 1 :
                      (bottom && !left) ? idx - Nx :
                      (!bottom && left) ? idx - 1 :
                      idx + Nx + 1;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];
    const real Vm_ne = Vm[idx_ne];
    const real Vm_nw = Vm[idx_nw];
    const real Vm_se = Vm[idx_se];
    const real Vm_sw = Vm[idx_sw];

    // Diffusion coefficients at center
    const real Dxx_c = Dxx[idx];
    const real Dxy_c = Dxy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right ? HARMONIC_MEAN(Dxx_c, Dxx[idx_r]) : Dxx_c;
    const real Dxx_l = left ? HARMONIC_MEAN(Dxx_c, Dxx[idx_l]) : Dxx_c;

    const real Dxy_r = right ? HARMONIC_MEAN(Dxy_c, Dxy[idx_r]) : Dxy_c;
    const real Dxy_l = left ? HARMONIC_MEAN(Dxy_c, Dxy[idx_l]) : Dxy_c;

    const real VmtVmb = Vm_t - Vm_b;
    const real denom_delta_x  = 1.0f / delta_x;
    const real denom_4delta_y = 0.25f * delta_y;

    // Fluxes
    const real J_r = right * ((Dxx_r * denom_delta_x) * (Vm_r - Vm_c) + (Dxy_r * denom_4delta_y) * ((Vm_ne - Vm_se) + VmtVmb));
    const real J_l = left  * ((Dxx_l * denom_delta_x) * (Vm_c - Vm_l) + (Dxy_l * denom_4delta_y) * (VmtVmb + (Vm_nw - Vm_sw)));

    // Divergence of flux -> nabla dot J
    return (J_r - J_l) * denom_delta_x;
}

MODIFIERS real compute_diffusion_term_y_no_rotation(const real *Vm, const real *Dyy,
                                                    const int i, const int j, const int Nx, const int Ny,
                                                    const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];

    // Diffusion coefficients at center
    const real Dyy_c = Dyy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dyy_t = top    ? HARMONIC_MEAN(Dyy_c, Dyy[idx_t]) : Dyy_c;
    const real Dyy_b = bottom ? HARMONIC_MEAN(Dyy_c, Dyy[idx_b]) : Dyy_c;

    const real denom_delta_y = 1.0f / delta_y;

    // Fluxes
    const real J_t = top    * ((Dyy_t * denom_delta_y) * (Vm_t - Vm_c));
    const real J_b = bottom * ((Dyy_b * denom_delta_y) * (Vm_c - Vm_b));

    // Divergence of flux -> nabla dot J
    return (J_t - J_b) * denom_delta_y;
}

MODIFIERS real compute_diffusion_term_y_with_rotation(const real *Vm, const real *Dyy, const real *Dxy,
                                                      const int i, const int j, const int Nx, const int Ny,
                                                      const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Diagonal neighbors with mirroring
    const int idx_ne = (top && right) ? idx + Nx + 1 : 
                      (top && !right) ? idx + Nx :
                      (!top && right) ? idx + 1 :
                      idx - Nx - 1;
    
    const int idx_nw = (top && left) ? idx + Nx - 1 :
                      (top && !left) ? idx + Nx :
                      (!top && left) ? idx - 1 :
                      idx - Nx + 1;
    
    const int idx_se = (bottom && right) ? idx - Nx + 1 :
                      (bottom && !right) ? idx - Nx :
                      (!bottom && right) ? idx + 1 :
                      idx + Nx - 1;
    
    const int idx_sw = (bottom && left) ? idx - Nx - 1 :
                      (bottom && !left) ? idx - Nx :
                      (!bottom && left) ? idx - 1 :
                      idx + Nx + 1;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];
    const real Vm_ne = Vm[idx_ne];
    const real Vm_nw = Vm[idx_nw];
    const real Vm_se = Vm[idx_se];
    const real Vm_sw = Vm[idx_sw];

    // Diffusion coefficients at center
    const real Dyy_c = Dyy[idx];
    const real Dxy_c = Dxy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dyy_t = top    ? HARMONIC_MEAN(Dyy_c, Dyy[idx_t]) : Dyy_c;
    const real Dyy_b = bottom ? HARMONIC_MEAN(Dyy_c, Dyy[idx_b]) : Dyy_c;

    const real Dxy_t = top    ? HARMONIC_MEAN(Dxy_c, Dxy[idx_t]) : Dxy_c;
    const real Dxy_b = bottom ? HARMONIC_MEAN(Dxy_c, Dxy[idx_b]) : Dxy_c;

    const real VmrVml = Vm_r - Vm_l;
    const real denom_delta_y  = 1.0f / delta_y;
    const real denom_4delta_x = 0.25f * delta_x;

    // Fluxes
    const real J_t = top    * ((Dyy_t * denom_delta_y) * (Vm_t - Vm_c) + (Dxy_t * denom_4delta_x) * ((Vm_ne - Vm_nw) + VmrVml));
    const real J_b = bottom * ((Dyy_b * denom_delta_y) * (Vm_c - Vm_b) + (Dxy_b * denom_4delta_x) * (VmrVml + (Vm_se - Vm_sw)));

    // Divergence of flux -> nabla dot J
    return (J_t - J_b) * denom_delta_y;
}

MODIFIERS real compute_diffusion_term_xy(const bool is_aligned, 
                                         const real *Vm, const real *Dxx, const real *Dyy, const real *Dxy,
                                         const int i, const int j, const int Nx, const int Ny,
                                         const real delta_x, const real delta_y)
{
    if (is_aligned) return 0.0f; // No diffusion in the xy direction for aligned fibers

    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left = (j > 0);
    const bool right = (j < Nx - 1);
    const bool top = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Diagonal neighbors with mirroring
    const int idx_ne = (top && right) ? idx + Nx + 1 : 
                      (top && !right) ? idx + Nx :
                      (!top && right) ? idx + 1 :
                      idx - Nx - 1;
    
    const int idx_nw = (top && left) ? idx + Nx - 1 :
                      (top && !left) ? idx + Nx :
                      (!top && left) ? idx - 1 :
                      idx - Nx + 1;
    
    const int idx_se = (bottom && right) ? idx - Nx + 1 :
                      (bottom && !right) ? idx - Nx :
                      (!bottom && right) ? idx + 1 :
                      idx + Nx - 1;
    
    const int idx_sw = (bottom && left) ? idx - Nx - 1 :
                      (bottom && !left) ? idx - Nx :
                      (!bottom && left) ? idx - 1 :
                      idx + Nx + 1;

    // Vm values
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];
    const real Vm_ne = Vm[idx_ne];
    const real Vm_nw = Vm[idx_nw];
    const real Vm_se = Vm[idx_se];
    const real Vm_sw = Vm[idx_sw];

    // Diffusion coefficients at center
    const real Dxy_c = Dxy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxy_r = right  ? HARMONIC_MEAN(Dxy_c, Dxy[idx_r]) : Dxy_c;
    const real Dxy_l = left   ? HARMONIC_MEAN(Dxy_c, Dxy[idx_l]) : Dxy_c;
    const real Dxy_t = top    ? HARMONIC_MEAN(Dxy_c, Dxy[idx_t]) : Dxy_c;
    const real Dxy_b = bottom ? HARMONIC_MEAN(Dxy_c, Dxy[idx_b]) : Dxy_c;

    const real VmtVmb = Vm_t - Vm_b;
    const real VmrVml = Vm_r - Vm_l;
    const real denom_delta_x = 1.0f / delta_x;
    const real denom_delta_y = 1.0f / delta_y;
    const real denom_4delta_y = 0.25f * delta_y;
    const real denom_4delta_x = 0.25f * delta_x;

    // Fluxes
    const real J_r = right  * ((Dxy_r * denom_4delta_y) * ((Vm_ne - Vm_se) + VmtVmb));
    const real J_l = left   * ((Dxy_l * denom_4delta_y) * (VmtVmb + (Vm_nw - Vm_sw)));
    const real J_t = top    * ((Dxy_t * denom_4delta_x) * ((Vm_ne - Vm_nw) + VmrVml));
    const real J_b = bottom * ((Dxy_b * denom_4delta_x) * (VmrVml + (Vm_se - Vm_sw)));

    // Divergence of flux -> nabla dot J
    return ((J_r - J_l) * denom_delta_x) + ((J_t - J_b) * denom_delta_y);
}

// Function to compute the diffusion term in 2D considering no rotation
MODIFIERS real compute_diffusion_term_no_rotation(const real *Vm, const real *Dxx, const real *Dyy,
                                                  const int i, const int j, const int Nx, const int Ny,
                                                  const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];

    // Diffusion coefficients at center
    const real Dxx_c = Dxx[idx];
    const real Dyy_c = Dyy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right  ? HARMONIC_MEAN(Dxx_c, Dxx[idx_r]) : Dxx_c;
    const real Dxx_l = left   ? HARMONIC_MEAN(Dxx_c, Dxx[idx_l]) : Dxx_c;
    const real Dyy_t = top    ? HARMONIC_MEAN(Dyy_c, Dyy[idx_t]) : Dyy_c;
    const real Dyy_b = bottom ? HARMONIC_MEAN(Dyy_c, Dyy[idx_b]) : Dyy_c;

    const real denom_delta_x = 1.0f / delta_x;
    const real denom_delta_y = 1.0f / delta_y;

    // Fluxes
    const real J_r = right  * ((Dxx_r * denom_delta_x) * (Vm_r - Vm_c));
    const real J_l = left   * ((Dxx_l * denom_delta_x) * (Vm_c - Vm_l));
    const real J_t = top    * ((Dyy_t * denom_delta_y) * (Vm_t - Vm_c));
    const real J_b = bottom * ((Dyy_b * denom_delta_y) * (Vm_c - Vm_b));

    // Divergence of flux -> nabla dot J
    return ((J_r - J_l) * denom_delta_x) + ((J_t - J_b) * denom_delta_y);
}

MODIFIERS real compute_diffusion_term_with_rotation(const real *Vm, const real *Dxx, const real *Dyy, const real *Dxy,
                                                    const int i, const int j, const int Nx, const int Ny,
                                                    const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Neighbor indices with mirroring
    const int idx_l = left   ? idx - 1 : idx;
    const int idx_r = right  ? idx + 1 : idx;
    const int idx_t = top    ? idx + Nx : idx;
    const int idx_b = bottom ? idx - Nx : idx;

    // Diagonal neighbors with mirroring
    const int idx_ne = (top && right) ? idx + Nx + 1 : 
                      (top && !right) ? idx + Nx :
                      (!top && right) ? idx + 1 :
                      idx - Nx - 1;
    
    const int idx_nw = (top && left) ? idx + Nx - 1 :
                      (top && !left) ? idx + Nx :
                      (!top && left) ? idx - 1 :
                      idx - Nx + 1;
    
    const int idx_se = (bottom && right) ? idx - Nx + 1 :
                      (bottom && !right) ? idx - Nx :
                      (!bottom && right) ? idx + 1 :
                      idx + Nx - 1;
    
    const int idx_sw = (bottom && left) ? idx - Nx - 1 :
                      (bottom && !left) ? idx - Nx :
                      (!bottom && left) ? idx - 1 :
                      idx + Nx + 1;

    // Diagonal neighbors with mirroring
    // const int idx_ne = (top && right) ? idx + Nx + 1 : 
    //                   (top && !right) ? idx + Nx - 1 :
    //                   (!top && right) ? idx - Nx + 1 :
    //                   idx - Nx - 1;
    
    // const int idx_nw = (top && left) ? idx + Nx - 1 :
    //                   (top && !left) ? idx + Nx + 1 :
    //                   (!top && left) ? idx - Nx - 1 :
    //                   idx - Nx + 1;
    
    // const int idx_se = (bottom && right) ? idx - Nx + 1 :
    //                   (bottom && !right) ? idx - Nx - 1 :
    //                   (!bottom && right) ? idx + Nx + 1 :
    //                   idx + Nx - 1;
    
    // const int idx_sw = (bottom && left) ? idx - Nx - 1 :
    //                   (bottom && !left) ? idx - Nx + 1 :
    //                   (!bottom && left) ? idx + Nx - 1 :
    //                   idx + Nx + 1;

    // Vm values
    const real Vm_c = Vm[idx];
    const real Vm_l = Vm[idx_l];
    const real Vm_r = Vm[idx_r];
    const real Vm_t = Vm[idx_t];
    const real Vm_b = Vm[idx_b];
    const real Vm_ne = Vm[idx_ne];
    const real Vm_nw = Vm[idx_nw];
    const real Vm_se = Vm[idx_se];
    const real Vm_sw = Vm[idx_sw];

    // Diffusion coefficients at center
    const real Dxx_c = Dxx[idx];
    const real Dyy_c = Dyy[idx];
    const real Dxy_c = Dxy[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right ? HARMONIC_MEAN(Dxx_c, Dxx[idx_r]) : Dxx_c;
    const real Dxy_r = right ? HARMONIC_MEAN(Dxy_c, Dxy[idx_r]) : Dxy_c;

    const real Dxx_l = left ? HARMONIC_MEAN(Dxx_c, Dxx[idx_l]) : Dxx_c;
    const real Dxy_l = left ? HARMONIC_MEAN(Dxy_c, Dxy[idx_l]) : Dxy_c;

    const real Dyy_t = top ? HARMONIC_MEAN(Dyy_c, Dyy[idx_t]) : Dyy_c;
    const real Dxy_t = top ? HARMONIC_MEAN(Dxy_c, Dxy[idx_t]) : Dxy_c;

    const real Dyy_b = bottom ? HARMONIC_MEAN(Dyy_c, Dyy[idx_b]) : Dyy_c;
    const real Dxy_b = bottom ? HARMONIC_MEAN(Dxy_c, Dxy[idx_b]) : Dxy_c;

    const real VmtVmb = Vm_t - Vm_b;
    const real VmrVml = Vm_r - Vm_l;
    const real denom_delta_x = 1.0f / delta_x;
    const real denom_delta_y = 1.0f / delta_y;
    const real denom_4delta_y = 0.25f * delta_y;
    const real denom_4delta_x = 0.25f * delta_x;

    // Fluxes
    const real J_r = right  ? ((Dxx_r * denom_delta_x) * (Vm_r - Vm_c) + (Dxy_r * denom_4delta_y) * ((Vm_ne - Vm_se) + VmtVmb)) : 0.0f;
    const real J_l = left   ? ((Dxx_l * denom_delta_x) * (Vm_c - Vm_l) + (Dxy_l * denom_4delta_y) * (VmtVmb + (Vm_nw - Vm_sw))) : 0.0f;
    const real J_t = top    ? ((Dyy_t * denom_delta_y) * (Vm_t - Vm_c) + (Dxy_t * denom_4delta_x) * ((Vm_ne - Vm_nw) + VmrVml)) : 0.0f;
    const real J_b = bottom ? ((Dyy_b * denom_delta_y) * (Vm_c - Vm_b) + (Dxy_b * denom_4delta_x) * (VmrVml + (Vm_se - Vm_sw))) : 0.0f;

    // Divergence of flux -> nabla dot J
    return ((J_r - J_l) * denom_delta_x) + ((J_t - J_b) * denom_delta_y);
}

MODIFIERS real compute_diffusion_term_with_rotation_3D(const real *Vm,
                                             const real *Dxx, const real *Dyy, const real *Dzz,
                                             const real *Dxy, const real *Dxz, const real *Dyz,
                                             const int i, const int j, const int k,
                                             const int Nx, const int Ny, const int Nz,
                                             const real dx, const real dy, const real dz)
{
    const int idx = (k * Ny + i) * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);
    const bool up     = (k < Nz - 1);
    const bool down   = (k > 0);

    // Main neighbors with mirroring
    const int idx_l = left   ? idx - 1         : idx;
    const int idx_r = right  ? idx + 1         : idx;
    const int idx_t = top    ? idx + Nx        : idx;
    const int idx_b = bottom ? idx - Nx        : idx;
    const int idx_u = up     ? idx + Nx*Ny     : idx;
    const int idx_d = down   ? idx - Nx*Ny     : idx;

    // Diagonals in XY plane
    const int idx_tr = (top && right)  ? idx + Nx + 1 :
                       (top && !right) ? idx + Nx     :
                       (!top && right) ? idx + 1      :
                       idx - Nx - 1;

    const int idx_tl = (top && left)   ? idx + Nx - 1 :
                       (top && !left)  ? idx + Nx     :
                       (!top && left)  ? idx - 1      :
                       idx - Nx + 1;

    const int idx_br = (bottom && right)  ? idx - Nx + 1 :
                       (bottom && !right) ? idx - Nx     :
                       (!bottom && right) ? idx + 1      :
                       idx + Nx - 1;

    const int idx_bl = (bottom && left)   ? idx - Nx - 1 :
                       (bottom && !left)  ? idx - Nx     :
                       (!bottom && left)  ? idx - 1      :
                       idx + Nx + 1;

    // Diagonals in XZ plane
    const int idx_ur = (up && right)    ? idx + Nx*Ny + 1 :
                       (up && !right)   ? idx + Nx*Ny     :
                       (!up && right)   ? idx + 1         :
                       idx - Nx*Ny - 1;

    const int idx_ul = (up && left)     ? idx + Nx*Ny - 1 :
                       (up && !left)    ? idx + Nx*Ny     :
                       (!up && left)    ? idx - 1         :
                       idx - Nx*Ny + 1;

    const int idx_dr = (down && right)  ? idx - Nx*Ny + 1 :
                       (down && !right) ? idx - Nx*Ny     :
                       (!down && right) ? idx + 1         :
                       idx + Nx*Ny - 1;

    const int idx_dl = (down && left)   ? idx - Nx*Ny - 1 :
                       (down && !left)  ? idx - Nx*Ny     :
                       (!down && left)  ? idx - 1         :
                       idx + Nx*Ny + 1;

    // Diagonals in YZ plane
    const int idx_ut = (up && top)     ? idx + Nx*Ny + Nx :
                       (up && !top)    ? idx + Nx*Ny      :
                       (!up && top)    ? idx + Nx         :
                       idx - Nx*Ny - Nx;

    const int idx_ub = (up && bottom)  ? idx + Nx*Ny - Nx :
                       (up && !bottom) ? idx + Nx*Ny      :
                       (!up && bottom) ? idx - Nx         :
                       idx - Nx*Ny + Nx;

    const int idx_dt = (down && top)   ? idx - Nx*Ny + Nx :
                       (down && !top)  ? idx - Nx*Ny      :
                       (!down && top)  ? idx + Nx         :
                       idx + Nx*Ny - Nx;

    const int idx_db = (down && bottom) ? idx - Nx*Ny - Nx :
                        (down && !bottom)? idx - Nx*Ny     :
                        (!down && bottom)? idx - Nx        :
                        idx + Nx*Ny + Nx;

    // Vm values
    const real Vm_c  = Vm[idx];
    const real Vm_l  = Vm[idx_l];
    const real Vm_r  = Vm[idx_r];
    const real Vm_b  = Vm[idx_b];
    const real Vm_t  = Vm[idx_t];
    const real Vm_d  = Vm[idx_d];
    const real Vm_u  = Vm[idx_u];

    const real Vm_tr = Vm[idx_tr];
    const real Vm_tl = Vm[idx_tl];
    const real Vm_br = Vm[idx_br];
    const real Vm_bl = Vm[idx_bl];

    const real Vm_ur = Vm[idx_ur];
    const real Vm_ul = Vm[idx_ul];
    const real Vm_dr = Vm[idx_dr];
    const real Vm_dl = Vm[idx_dl];

    const real Vm_ut = Vm[idx_ut];
    const real Vm_ub = Vm[idx_ub];
    const real Vm_dt = Vm[idx_dt];
    const real Vm_db = Vm[idx_db];

    // Diffusion coefficients at center
    const real Dxx_c = Dxx[idx], Dyy_c = Dyy[idx], Dzz_c = Dzz[idx];
    const real Dxy_c = Dxy[idx], Dxz_c = Dxz[idx], Dyz_c = Dyz[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right ? HARMONIC_MEAN(Dxx_c, Dxx[idx_r]) : Dxx_c;
    const real Dxx_l = left ? HARMONIC_MEAN(Dxx_c, Dxx[idx_l]) : Dxx_c;

    const real Dyy_t = top ? HARMONIC_MEAN(Dyy_c, Dyy[idx_t]) : Dyy_c;
    const real Dyy_b = bottom ? HARMONIC_MEAN(Dyy_c, Dyy[idx_b]) : Dyy_c;

    const real Dzz_u = up ? HARMONIC_MEAN(Dzz_c, Dzz[idx_u]) : Dzz_c;
    const real Dzz_d = down ? HARMONIC_MEAN(Dzz_c, Dzz[idx_d]) : Dzz_c;

    const real Dxy_r = right ? HARMONIC_MEAN(Dxy_c, Dxy[idx_r]) : Dxy_c;
    const real Dxy_l = left ? HARMONIC_MEAN(Dxy_c, Dxy[idx_l]) : Dxy_c;
    const real Dxy_t = top ? HARMONIC_MEAN(Dxy_c, Dxy[idx_t]) : Dxy_c;
    const real Dxy_b = bottom ? HARMONIC_MEAN(Dxy_c, Dxy[idx_b]) : Dxy_c;

    const real Dxz_r = right ? HARMONIC_MEAN(Dxz_c, Dxz[idx_r]) : Dxz_c;
    const real Dxz_l = left ? HARMONIC_MEAN(Dxz_c, Dxz[idx_l]) : Dxz_c;
    const real Dxz_u = up ? HARMONIC_MEAN(Dxz_c, Dxz[idx_u]) : Dxz_c;
    const real Dxz_d = down ? HARMONIC_MEAN(Dxz_c, Dxz[idx_d]) : Dxz_c;

    const real Dyz_t = top ? HARMONIC_MEAN(Dyz_c, Dyz[idx_t]) : Dyz_c;
    const real Dyz_b = bottom ? HARMONIC_MEAN(Dyz_c, Dyz[idx_b]) : Dyz_c;
    const real Dyz_u = up ? HARMONIC_MEAN(Dyz_c, Dyz[idx_u]) : Dyz_c;
    const real Dyz_d = down ? HARMONIC_MEAN(Dyz_c, Dyz[idx_d]) : Dyz_c;

    // Precomputed finite difference terms
    const real denom_dx = 1.0f / dx;
    const real denom_dy = 1.0f / dy;
    const real denom_dz = 1.0f / dz;
    const real denom_4dx = 0.25f * dx;
    const real denom_4dy = 0.25f * dy;
    const real denom_4dz = 0.25f * dz;

    const real Vm_t_minus_b = Vm_t - Vm_b;
    const real Vm_r_minus_l = Vm_r - Vm_l;
    const real Vm_u_minus_d = Vm_u - Vm_d;

    // Fluxes
    const real J_r = right ? ((Dxx_r * denom_dx) * (Vm_r - Vm_c)
                              + (Dxy_r * denom_4dy) * ((Vm_tr - Vm_br) + Vm_t_minus_b)
                              + (Dxz_r * denom_4dz) * ((Vm_ur - Vm_dr) + Vm_u_minus_d)) : 0.0f;

    const real J_l = left ? ((Dxx_l * denom_dx) * (Vm_c - Vm_l)
                              + (Dxy_l * denom_4dy) * (Vm_t_minus_b + (Vm_tl - Vm_bl))
                              + (Dxz_l * denom_4dz) * (Vm_u_minus_d + (Vm_ul - Vm_dl))) : 0.0f;

    const real J_t = top ? ((Dyy_t * denom_dy) * (Vm_t - Vm_c)
                              + (Dxy_t * denom_4dx) * ((Vm_tr - Vm_tl) + Vm_r_minus_l)
                              + (Dyz_t * denom_4dz) * ((Vm_ut - Vm_dt) + Vm_u_minus_d)) : 0.0f;

    const real J_b = bottom ? ((Dyy_b * denom_dy) * (Vm_c - Vm_b)
                              + (Dxy_b * denom_4dx) * (Vm_r_minus_l + (Vm_br - Vm_bl))
                              + (Dyz_b * denom_4dz) * (Vm_u_minus_d + (Vm_ub - Vm_db))) : 0.0f;

    const real J_u = up ? ((Dzz_u * denom_dz) * (Vm_u - Vm_c)
                              + (Dxz_u * denom_4dx) * ((Vm_ur - Vm_ul) + Vm_r_minus_l)
                              + (Dyz_u * denom_4dy) * ((Vm_ut - Vm_ub) + Vm_t_minus_b)) : 0.0f;

    const real J_d = down ? ((Dzz_d * denom_dz) * (Vm_c - Vm_d)
                              + (Dxz_d * denom_4dx) * (Vm_r_minus_l + (Vm_dr - Vm_dl))
                              + (Dyz_d * denom_4dy) * (Vm_t_minus_b + (Vm_dt - Vm_db))) : 0.0f;

    // Divergence of flux: ∇·J
    return ((J_r - J_l) * denom_dx) +
           ((J_t - J_b) * denom_dy) +
           ((J_u - J_d) * denom_dz);
}


MODIFIERS real select_compute_diffusion_term_x(const bool is_aligned,
                                               const real *Vm, const real *Dxx, const real *Dxy,
                                               const int i, const int j, const int Nx, const int Ny,
                                               const real delta_x, const real delta_y)
{
    real diffusion_term;

    diffusion_term = is_aligned
                            ? compute_diffusion_term_x_no_rotation(Vm, Dxx, i, j, Nx, Ny, delta_x, delta_y)
                            : compute_diffusion_term_x_with_rotation(Vm, Dxx, Dxy, i, j, Nx, Ny, delta_x, delta_y);

    return diffusion_term;
}

MODIFIERS real select_compute_diffusion_term_y(const bool is_aligned,
                                               const real *Vm, const real *Dyy, const real *Dxy,
                                               const int i, const int j, const int Nx, const int Ny,
                                               const real delta_x, const real delta_y)
{
    real diffusion_term;

    diffusion_term = is_aligned
                            ? compute_diffusion_term_y_no_rotation(Vm, Dyy, i, j, Nx, Ny, delta_x, delta_y)
                            : compute_diffusion_term_y_with_rotation(Vm, Dyy, Dxy, i, j, Nx, Ny, delta_x, delta_y);

    return diffusion_term;
}

MODIFIERS real select_compute_diffusion_term(const bool is_aligned,
                                             const real *Vm, const real *Dxx, const real *Dyy, const real *Dxy,
                                             const int i, const int j, const int Nx, const int Ny,
                                             const real delta_x, const real delta_y)
{
    real diffusion_term;

    diffusion_term = is_aligned
                            ? compute_diffusion_term_no_rotation(Vm, Dxx, Dyy, i, j, Nx, Ny, delta_x, delta_y)
                            : compute_diffusion_term_with_rotation(Vm, Dxx, Dyy, Dxy, i, j, Nx, Ny, delta_x, delta_y);

    return diffusion_term;
}

#if defined(__CUDACC__)

__forceinline__ __device__ void select_get_actual_sV(const CellModel cell_model, real *actualsV, const real *sV, const int idx)
{
    switch (cell_model)
    {
    case CELL_MODEL_AFHN:
        d_get_actual_sV_AFHN(actualsV, sV, idx);
        break;
    case CELL_MODEL_MV:
        d_get_actual_sV_MV(actualsV, sV, idx);
        break;
    default:
        printf("Error: Unsupported cell model\n");
        break;
    }
}

__forceinline__ __device__ real select_compute_dVmdt(const CellModel cell_model, const real Vm, const real *sV)
{
    switch (cell_model)
    {
    case CELL_MODEL_AFHN:
        return d_compute_dVmdt_AFHN(Vm, sV);
    case CELL_MODEL_MV:
        return d_compute_dVmdt_MV(Vm, sV);
    default:
        printf("Error: Unsupported cell model\n");
        return 0.0f;
    }
}

__forceinline__ __device__ void select_update_sVtilde(const CellModel cell_model, real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t)
{
    switch (cell_model)
    {
    case CELL_MODEL_AFHN:
        d_update_sVtilde_AFHN(sVtilde, Vm, rhs_sV, delta_t);
        break;
    case CELL_MODEL_MV:
        d_update_sVtilde_MV(sVtilde, Vm, rhs_sV, delta_t);
        break;
    default:
        printf("Error: Unsupported cell model\n");
        break;
    }
}

__forceinline__ __device__ void select_update_sV(const CellModel cell_model, real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    switch (cell_model)
    {
    case CELL_MODEL_AFHN:
        d_update_sV_AFHN(sV, rhs_sV, dSdt_Vm, dSdt_sV, delta_t, idx);
        break;
    case CELL_MODEL_MV:
        d_update_sV_MV(sV, rhs_sV, dSdt_Vm, dSdt_sV, delta_t, idx);
        break;
    default:
        printf("Error: Unsupported cell model\n");
        break;
    }
}

#endif // __CUDACC__

#ifdef __cplusplus
}
#endif

#endif // NUMERICAL_METHODS_HELPERS_H