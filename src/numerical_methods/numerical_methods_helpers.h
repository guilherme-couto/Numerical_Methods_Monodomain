#ifndef NUMERICAL_METHODS_HELPERS_H
#define NUMERICAL_METHODS_HELPERS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../logger/logger.h"

#ifdef __cplusplus
extern "C" {
#endif

// Function to get the stimulus value at a given time and position
STATIC_MODIFIER real get_stimulus_value(const real actualTime, const int i, const int j, const Stimulus *stimuli, const int numberOfStimuli)
{
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

// Funtion to determine the limits of the grid and the boundaries condition
STATIC_MODIFIER int lim(const int num, const int N)
{
    return num == -1 ? 1 : (num == N ? N - 2 : num);
}

STATIC_MODIFIER real harmonic_mean(const real s1, const real s2)
{
    return (s1 + s2) > 0.0f ? 2.0f * s1 * s2 / (s1 + s2) : 0.0f;
}

STATIC_MODIFIER real compute_diffusion_term_x_axis(const real *Vm, const ElementProperties *elements,
                                                   const int i, const int j, const int Nx, const int Ny,
                                                   const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left  = (j > 0);
    const bool right = (j < Nx - 1);
    const bool top   = (i < Ny - 1);
    const bool bottom = (i > 0);

    const bool ne = top && right;
    const bool nw = top && left;
    const bool se = bottom && right;
    const bool sw = bottom && left;

    // Indexing
    const int idx_l = i * Nx + (j - 1);
    const int idx_r = i * Nx + (j + 1);
    const int idx_t = (i + 1) * Nx + j;
    const int idx_b = (i - 1) * Nx + j;

    const int idx_ne = (i + 1) * Nx + (j + 1);
    const int idx_nw = (i + 1) * Nx + (j - 1);
    const int idx_se = (i - 1) * Nx + (j + 1);
    const int idx_sw = (i - 1) * Nx + (j - 1);

    // Vm values
    const real Vm_c  = Vm[idx];
    const real Vm_l  = left   ? Vm[idx_l]  : 0.0f;
    const real Vm_r  = right  ? Vm[idx_r]  : 0.0f;
    const real Vm_t  = top    ? Vm[idx_t]  : 0.0f;
    const real Vm_b  = bottom ? Vm[idx_b]  : 0.0f;
    
    const real Vm_ne = ne     ? Vm[idx_ne] : 0.0f;
    const real Vm_nw = nw     ? Vm[idx_nw] : 0.0f;
    const real Vm_se = se     ? Vm[idx_se] : 0.0f;
    const real Vm_sw = sw     ? Vm[idx_sw] : 0.0f;

    // Diffusion coefficients at center
    const ElementProperties center = elements[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right  ? harmonic_mean(center.D_xx, elements[idx_r].D_xx) : center.D_xx;
    const real Dxy_r = right  ? harmonic_mean(center.D_xy, elements[idx_r].D_xy) : center.D_xy;

    const real Dxx_l = left   ? harmonic_mean(center.D_xx, elements[idx_l].D_xx) : center.D_xx;
    const real Dxy_l = left   ? harmonic_mean(center.D_xy, elements[idx_l].D_xy) : center.D_xy;

    // Fluxes
    const real J_r = right ? (
        (Dxx_r / delta_x) * (Vm_r - Vm_c)
        +(Dxy_r / (4.0f * delta_y)) * ((Vm_ne - Vm_se) + (Vm_t - Vm_b))
    ) : 0.0f;

    const real J_l = left ? (
        (Dxx_l / delta_x) * (Vm_c - Vm_l)
        +(Dxy_l / (4.0f * delta_y)) * ((Vm_t - Vm_b) + (Vm_nw - Vm_sw))
    ) : 0.0f;

    // Divergence of flux -> nabla dot J
    const real x_term = (J_r - J_l) / delta_x;

    return x_term;
}

STATIC_MODIFIER real compute_diffusion_term_y_axis(const real *Vm, const ElementProperties *elements,
                                                   const int i, const int j, const int Nx, const int Ny,
                                                   const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    const bool ne = top && right;
    const bool nw = top && left;
    const bool se = bottom && right;
    const bool sw = bottom && left;

    // Indexing
    const int idx_l = i * Nx + (j - 1);
    const int idx_r = i * Nx + (j + 1);
    const int idx_t = (i + 1) * Nx + j;
    const int idx_b = (i - 1) * Nx + j;

    const int idx_ne = (i + 1) * Nx + (j + 1);
    const int idx_nw = (i + 1) * Nx + (j - 1);
    const int idx_se = (i - 1) * Nx + (j + 1);
    const int idx_sw = (i - 1) * Nx + (j - 1);

    // Vm values
    const real Vm_c  = Vm[idx];
    const real Vm_l  = left   ? Vm[idx_l]  : 0.0f;
    const real Vm_r  = right  ? Vm[idx_r]  : 0.0f;
    const real Vm_t  = top    ? Vm[idx_t]  : 0.0f;
    const real Vm_b  = bottom ? Vm[idx_b]  : 0.0f;
    
    const real Vm_ne = ne     ? Vm[idx_ne] : 0.0f;
    const real Vm_nw = nw     ? Vm[idx_nw] : 0.0f;
    const real Vm_se = se     ? Vm[idx_se] : 0.0f;
    const real Vm_sw = sw     ? Vm[idx_sw] : 0.0f;

    // Diffusion coefficients at center
    const ElementProperties center = elements[idx];

    // Interpolated coefficients using harmonic mean
    const real Dyy_t = top    ? harmonic_mean(center.D_yy, elements[idx_t].D_yy) : center.D_yy;
    const real Dxy_t = top    ? harmonic_mean(center.D_xy, elements[idx_t].D_xy) : center.D_xy;

    const real Dyy_b = bottom ? harmonic_mean(center.D_yy, elements[idx_b].D_yy) : center.D_yy;
    const real Dxy_b = bottom ? harmonic_mean(center.D_xy, elements[idx_b].D_xy) : center.D_xy;

    // Fluxes
    const real J_t = top ? (
        (Dyy_t / delta_y) * (Vm_t - Vm_c)
        +(Dxy_t / (4.0f * delta_x)) * ((Vm_ne - Vm_nw) + (Vm_r - Vm_l))
    ) : 0.0f;

    const real J_b = bottom ? (
        (Dyy_b / delta_y) * (Vm_c - Vm_b)
        +(Dxy_b / (4.0f * delta_x)) * ((Vm_r - Vm_l) + (Vm_se - Vm_sw))
    ) : 0.0f;

    // Divergence of flux -> nabla dot J
    const real y_term = (J_t - J_b) / delta_y;

    return y_term;
}

// Function to compute the diffusion term in 2D considering no rotation
STATIC_MODIFIER real compute_diffusion_term_no_rotation(const real *Vm, const ElementProperties *elements,
                                                        const int i, const int j, const int Nx, const int Ny,
                                                        const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    // Indexing
    const int idx_l = i * Nx + (j - 1);
    const int idx_r = i * Nx + (j + 1);
    const int idx_t = (i + 1) * Nx + j;
    const int idx_b = (i - 1) * Nx + j;

    // Vm values
    const real Vm_c  = Vm[idx];
    const real Vm_l  = left   ? Vm[idx_l]  : 0.0f;
    const real Vm_r  = right  ? Vm[idx_r]  : 0.0f;
    const real Vm_t  = top    ? Vm[idx_t]  : 0.0f;
    const real Vm_b  = bottom ? Vm[idx_b]  : 0.0f;

    // Diffusion coefficients at center
    const ElementProperties center = elements[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right  ? harmonic_mean(center.D_xx, elements[idx_r].D_xx)  : center.D_xx;
    const real Dxx_l = left   ? harmonic_mean(center.D_xx, elements[idx_l].D_xx)  : center.D_xx;
    const real Dyy_t = top    ? harmonic_mean(center.D_yy, elements[idx_t].D_yy)  : center.D_yy;
    const real Dyy_b = bottom ? harmonic_mean(center.D_yy, elements[idx_b].D_yy)  : center.D_yy;

    // Fluxes
    const real J_r = right  ? ((Dxx_r / delta_x) * (Vm_r - Vm_c))  : 0.0f;
    const real J_l = left   ? ((Dxx_l / delta_x) * (Vm_c - Vm_l))  : 0.0f;
    const real J_t = top    ? ((Dyy_t / delta_y) * (Vm_t - Vm_c))  : 0.0f;
    const real J_b = bottom ? ((Dyy_b / delta_y) * (Vm_c - Vm_b))  : 0.0f;

    // Divergence of flux -> nabla dot J
    const real x_term = (J_r - J_l) / delta_x;
    const real y_term = (J_t - J_b) / delta_y;

    return x_term + y_term;
}

STATIC_MODIFIER real compute_diffusion_term_anisotropic(const real *Vm, const ElementProperties *elements,
                                                        const int i, const int j, const int Nx, const int Ny,
                                                        const real delta_x, const real delta_y)
{
    const int idx = i * Nx + j;

    // Neighbor validation
    const bool left   = (j > 0);
    const bool right  = (j < Nx - 1);
    const bool top    = (i < Ny - 1);
    const bool bottom = (i > 0);

    const bool ne = top && right;
    const bool nw = top && left;
    const bool se = bottom && right;
    const bool sw = bottom && left;

    // Indexing
    const int idx_l = i * Nx + (j - 1);
    const int idx_r = i * Nx + (j + 1);
    const int idx_t = (i + 1) * Nx + j;
    const int idx_b = (i - 1) * Nx + j;

    const int idx_ne = (i + 1) * Nx + (j + 1);
    const int idx_nw = (i + 1) * Nx + (j - 1);
    const int idx_se = (i - 1) * Nx + (j + 1);
    const int idx_sw = (i - 1) * Nx + (j - 1);

    // Vm values
    const real Vm_c  = Vm[idx];
    const real Vm_l  = left   ? Vm[idx_l]  : 0.0f;
    const real Vm_r  = right  ? Vm[idx_r]  : 0.0f;
    const real Vm_t  = top    ? Vm[idx_t]  : 0.0f;
    const real Vm_b  = bottom ? Vm[idx_b]  : 0.0f;
    
    const real Vm_ne = ne     ? Vm[idx_ne] : 0.0f;
    const real Vm_nw = nw     ? Vm[idx_nw] : 0.0f;
    const real Vm_se = se     ? Vm[idx_se] : 0.0f;
    const real Vm_sw = sw     ? Vm[idx_sw] : 0.0f;

    // Diffusion coefficients at center
    const ElementProperties center = elements[idx];

    // Interpolated coefficients using harmonic mean
    const real Dxx_r = right  ? harmonic_mean(center.D_xx, elements[idx_r].D_xx) : center.D_xx;
    const real Dxy_r = right  ? harmonic_mean(center.D_xy, elements[idx_r].D_xy) : center.D_xy;

    const real Dxx_l = left   ? harmonic_mean(center.D_xx, elements[idx_l].D_xx) : center.D_xx;
    const real Dxy_l = left   ? harmonic_mean(center.D_xy, elements[idx_l].D_xy) : center.D_xy;

    const real Dyy_t = top    ? harmonic_mean(center.D_yy, elements[idx_t].D_yy) : center.D_yy;
    const real Dxy_t = top    ? harmonic_mean(center.D_xy, elements[idx_t].D_xy) : center.D_xy;

    const real Dyy_b = bottom ? harmonic_mean(center.D_yy, elements[idx_b].D_yy) : center.D_yy;
    const real Dxy_b = bottom ? harmonic_mean(center.D_xy, elements[idx_b].D_xy) : center.D_xy;

    // Fluxes
    const real J_r = right ? (
        (Dxx_r / delta_x) * (Vm_r - Vm_c)
        +(Dxy_r / (4.0f * delta_y)) * ((Vm_ne - Vm_se) + (Vm_t - Vm_b))
    ) : 0.0f;

    const real J_l = left ? (
        (Dxx_l / delta_x) * (Vm_c - Vm_l)
        +(Dxy_l / (4.0f * delta_y)) * ((Vm_t - Vm_b) + (Vm_nw - Vm_sw))
    ) : 0.0f;

    const real J_t = top ? (
        (Dyy_t / delta_y) * (Vm_t - Vm_c)
        +(Dxy_t / (4.0f * delta_x)) * ((Vm_ne - Vm_nw) + (Vm_r - Vm_l))
    ) : 0.0f;

    const real J_b = bottom ? (
        (Dyy_b / delta_y) * (Vm_c - Vm_b)
        +(Dxy_b / (4.0f * delta_x)) * ((Vm_r - Vm_l) + (Vm_se - Vm_sw))
    ) : 0.0f;

    // Divergence of flux -> nabla dot J
    const real x_term = (J_r - J_l) / delta_x;
    const real y_term = (J_t - J_b) / delta_y;

    return x_term + y_term;
}

// Function to get the active stimuli
static int update_and_get_num_active_stimuli(const real actualTime, const Stimulus *stimuli, const int numberOfStimuli, Stimulus *active_stimuli)
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
static void handle_frame_saving(const char *pathToSaveData, const char *file_extension,
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
static void handle_velocity_measurement(const real Vm_x0, const real Vm_x1, real *t0, real *t1, const real thereshold,
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

#ifdef __cplusplus
}
#endif

#endif // NUMERICAL_METHODS_HELPERS_H