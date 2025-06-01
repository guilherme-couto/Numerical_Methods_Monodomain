// ==============================================
//                Douglas Scheme
// ==============================================

#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"

#ifndef DO_THETA
#define DO_THETA 0.5f
#endif // DO_THETA

void runDO(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const ElementProperties *elements)
{
    // Unpack configuration parameters
    const int M = config->M;
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const int numberOfStimuli = config->stimulus_count;
    const Stimulus *stimuli = config->stimuli;

    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *saveFunctionName = config->save_function_name;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Get the solver functions
    const real denom_chiCm = cell_model_solver->denom_chiCm;
    const real activation_threshold = cell_model_solver->activation_thershold;
    const get_actual_sV_t get_actual_sV = cell_model_solver->get_actual_sV;
    const compute_dVmdt_t compute_dVmdt = cell_model_solver->compute_dVmdt;
    const update_sVtilde_t update_sVtilde = cell_model_solver->update_sVtilde;
    const update_sV_t update_sV = cell_model_solver->update_sV;

    // Measure velocity variables
    real stim_velocity, t0, t1;
    const real x0 = config->Lx / 3.0f;
    const real x1 = 2.0f * x0;
    const int idx_x0 = round(x0 / delta_x) + 1;
    const int idx_x1 = round(x1 / delta_x) + 1;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, num_active_stimuli;
    int idx, idx_left, idx_right, idx_top, idx_bottom;

    // Auxiliary variables for the operations
    Stimulus *active_stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    real diff_term, stim, actualVm;
    real Vmtilde, dVmdt, prevY;
    real *actualsV = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    real *sVtilde = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    real *Y = (real *)malloc(Nx * Ny * sizeof(real));
    real *auxVm = (real *)malloc(Nx * Ny * sizeof(real));

    // Calculate auxiliary coefficients for Thomas algorithm
    const real thomas_coeff_x = (DO_THETA * delta_t * denom_chiCm) / (delta_x * delta_x);
    const real thomas_coeff_y = (DO_THETA * delta_t * denom_chiCm) / (delta_y * delta_y);

    // Auxiliary arrays for Thomas algorithm
    real *c_prime, *d_prime, *ls_rhs, *result;
    if (Nx > Ny)
    {
        c_prime = (real *)malloc(Nx * sizeof(real));
        d_prime = (real *)malloc(Nx * sizeof(real));
        ls_rhs = (real *)malloc(Nx * sizeof(real));
        result = (real *)malloc(Nx * sizeof(real));
    }
    else
    {
        c_prime = (real *)malloc(Ny * sizeof(real));
        d_prime = (real *)malloc(Ny * sizeof(real));
        ls_rhs = (real *)malloc(Ny * sizeof(real));
        result = (real *)malloc(Ny * sizeof(real));
    }

    // Variables for time measurement
    real startTime = 0.0f;
    real startExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;
    real startLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;

    SIMPLEMSG("");
    INFOMSG("Starting simulation with Do (SERIAL)...\n");

    // Main time loop
    startExecutionTime = omp_get_wtime();

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // Update the active stimuli
        num_active_stimuli = update_and_get_num_active_stimuli(actualTime, stimuli, numberOfStimuli, active_stimuli);

        // =================================================
        //  Calculate Y0 and Update ODEs
        // =================================================
        startTime = omp_get_wtime();

        diff_term = 0.0f;
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[idx];
                get_actual_sV(actualsV, sV, idx);

                // Stimulation
                stim = (num_active_stimuli > 0)
                           ? (get_stimulus_value(actualTime, i, j, active_stimuli, num_active_stimuli))
                           : (0.0f);

                // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*dt*(diffusion + R(Vmn, sV))
                diff_term = compute_diffusion_term_anisotropic(Vm, elements, i, j, Nx, Ny, delta_x, delta_y);
                dVmdt = compute_dVmdt(actualVm, actualsV);
                Vmtilde = actualVm + 0.5f * delta_t * ((diff_term * denom_chiCm) + stim - dVmdt);

                // Calculate approximation for state variables using Vmtilde and update them
                update_sVtilde(sVtilde, actualVm, actualsV, 0.5f * delta_t);
                update_sV(sV, actualsV, Vmtilde, sVtilde, delta_t, idx);

                // Update Y with the new approximation -> Y_0 = U_n-1 + dt * F(t_n-1, U_n-1)
                Y[idx] = actualVm + delta_t * ((diff_term * denom_chiCm) + stim - dVmdt);
            }
        }

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // =================================================
        //  First step -> Result goes to d_RHS
        //  diffusion implicit and explicit in x
        // =================================================
        startTime = omp_get_wtime();

        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along x
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                prevY = Y[idx];
                diff_term = compute_diffusion_term_x_axis(Vm, elements, i, j, Nx, Ny, delta_x, delta_y);

                // Y_1 = Y_0 + theta * dt * F_1(t_n-1, U_n-1) -> 0 represents the first step and 1 represents the x-axis
                ls_rhs[j] = prevY - DO_THETA * delta_t * diff_term * denom_chiCm;
            }

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver_x(Nx, ls_rhs, result, c_prime, d_prime, thomas_coeff_x, elements, i);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                auxVm[idx] = result[j];
            }

            elapsedTime1stLS += omp_get_wtime() - startLSTime;
        }

        // =================================================
        //  Second step -> Result goes to Vm
        //  diffusion implicit and explicit in y
        // =================================================
        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along y
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                prevY = auxVm[idx];
                diff_term = compute_diffusion_term_y_axis(Vm, elements, i, j, Nx, Ny, delta_x, delta_y);

                // Y_2 = Y_1 + theta * dt * F_2(t_n-1, U_n-1) -> 1 represents the x-axis and 2 represents the y-axis
                ls_rhs[i] = prevY - DO_THETA * delta_t * diff_term * denom_chiCm;
            }

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver_y(Ny, ls_rhs, result, c_prime, d_prime, thomas_coeff_y, elements, j, Nx);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                Vm[idx] = result[i];
            }

            elapsedTime2ndLS += omp_get_wtime() - startLSTime;
        }

        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        if (saveFrames && (timeStepCounter % frameSaveRate == 0))
        {
            startTime = omp_get_wtime();
            handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, Vm, Nx, Ny, delta_x, delta_y, actualTime);
            elapsedSaveFramesTime += omp_get_wtime() - startTime;
        }

        // Measure velocity if needed
        if (measureVelocity && !stim_velocity_measured)
        {
            startTime = omp_get_wtime();
            handle_velocity_measurement(Vm[idx_x0], Vm[idx_x1], &t0, &t1, activation_threshold, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
            elapsedMeasureVelocityTime += omp_get_wtime() - startTime;
        }

        // Update time step counter
        timeStepCounter++;
    }

    elapsedExecutionTime = omp_get_wtime() - startExecutionTime;

    INFOMSG("Simulation done!\n");
    SIMPLEMSG("");

    // Update measurement structure
    measurement->elapsedExecutionTime = elapsedExecutionTime;
    measurement->elapsedTime1stPart = elapsedTime1stPart;
    measurement->elapsedTime2ndPart = elapsedTime2ndPart;
    measurement->elapsedTime1stLS = elapsedTime1stLS;
    measurement->elapsedTime2ndLS = elapsedTime2ndLS;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Free allocated memory
    free(active_stimuli);
    free(actualsV);
    free(sVtilde);
    free(Y);
    free(auxVm);
    free(c_prime);
    free(d_prime);
    free(ls_rhs);
    free(result);
}