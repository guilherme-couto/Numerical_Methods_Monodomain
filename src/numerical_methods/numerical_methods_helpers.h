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

// Function to compute the diffusion term in 2D considering isotropic diffusion
STATIC_MODIFIER real compute_diffusion_term(const real *Vm, const int i, const int j, const int Nx, const int Ny,
                                            const real diff_coeff_x, const real diff_coeff_y, const real phi_x, const real phi_y)
{
    // Get the neighboring values
    const int idx = i * Nx + j;
    const int idx_left = i * Nx + lim(j - 1, Nx);
    const int idx_right = i * Nx + lim(j + 1, Nx);
    const int idx_top = lim(i + 1, Ny) * Nx + j;
    const int idx_bottom = lim(i - 1, Ny) * Nx + j;

    const real Vm_center = Vm[idx];
    const real Vm_left = Vm[idx_left];
    const real Vm_right = Vm[idx_right];
    const real Vm_top = Vm[idx_top];
    const real Vm_bottom = Vm[idx_bottom];

    // Compute the diffusion term
    const real x_term = Vm_left - 2.0f * Vm_center + Vm_right;
    const real y_term = Vm_top - 2.0f * Vm_center + Vm_bottom;

    return diff_coeff_x * phi_x * x_term +
           diff_coeff_y * phi_y * y_term;
}

// Funtion to determine the limits of the grid and the boundaries condition
STATIC_MODIFIER bool validate_boundary(const int num, const int N)
{
    // TODO: check if this is correct and second order accurate
    return num == -1 ? false : (num == N ? false : true);
}

// Function to compute the diffusion term in 2D considering anisotropic diffusion with cross-terms
STATIC_MODIFIER real compute_diffusion_term_anisotropic(const real *Vm, const int i, const int j, const int Nx, const int Ny,
                                                        const real diff_coeff_x, const real diff_coeff_y, const real diff_coeff_xy,
                                                        const real phi_x, const real phi_y, const real phi_xy)
{
    // Get the neighboring values (TODO: check if this is correct and second order accurate)
   
    const bool is_left_valid = validate_boundary(j - 1, Nx);
    const bool is_right_valid = validate_boundary(j + 1, Nx);
    const bool is_top_valid = validate_boundary(i + 1, Ny);
    const bool is_bottom_valid = validate_boundary(i - 1, Ny);
    const bool is_northeast_valid = validate_boundary(i + 1, Nx) && validate_boundary(j + 1, Ny);
    const bool is_northwest_valid = validate_boundary(i + 1, Nx) && validate_boundary(j - 1, Ny);
    const bool is_southeast_valid = validate_boundary(i - 1, Nx) && validate_boundary(j + 1, Ny);
    const bool is_southwest_valid = validate_boundary(i - 1, Nx) && validate_boundary(j - 1, Ny);

    const int idx = i * Nx + j;
    const int idx_left = i * Nx + (j - 1);
    const int idx_right = i * Nx + (j + 1);
    const int idx_top = (i + 1) * Nx + j;
    const int idx_bottom = (i - 1) * Nx + j;
    const int idx_northeast = (i + 1) * Nx + (j + 1);
    const int idx_northwest = (i + 1) * Nx + (j - 1);
    const int idx_southeast = (i - 1) * Nx + (j + 1);
    const int idx_southwest = (i - 1) * Nx + (j - 1);

    const real Vm_center = Vm[idx];
    const real Vm_left = is_left_valid ? Vm[idx_left] : 0.0f;
    const real Vm_right = is_right_valid ? Vm[idx_right] : 0.0f;
    const real Vm_top = is_top_valid ? Vm[idx_top] : 0.0f;
    const real Vm_bottom = is_bottom_valid ? Vm[idx_bottom] : 0.0f;
    const real Vm_northeast = is_northeast_valid ? Vm[idx_northeast] : 0.0f;
    const real Vm_northwest = is_northwest_valid ? Vm[idx_northwest] : 0.0f;
    const real Vm_southeast = is_southeast_valid ? Vm[idx_southeast] : 0.0f;
    const real Vm_southwest = is_southwest_valid ? Vm[idx_southwest] : 0.0f;

    // Compute the diffusion term
    const real x_term = Vm_left - 2.0f * Vm_center + Vm_right;
    const real y_term = Vm_top - 2.0f * Vm_center + Vm_bottom;
    const real xy_term = Vm_northeast - Vm_northwest - Vm_southeast + Vm_southwest;

    return diff_coeff_x  * phi_x  * x_term +
           diff_coeff_y  * phi_y  * y_term +
           diff_coeff_xy * phi_xy * xy_term;
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