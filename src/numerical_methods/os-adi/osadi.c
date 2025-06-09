// ==============================================
//             Operator-Splitting ADI
// ==============================================

#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"

void runOSADI(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy)
{
    // Unpack configuration parameters
    const int M = config->M;
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const bool is_aligned = config->is_fiber_aligned;
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
    int idx;

    // Auxiliary variables for the operations
    Stimulus *active_stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    real stim, actualVm;
    real *actualsV = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    real *partRHS = (real *)malloc(Nx * Ny * sizeof(real));

    // Calculate auxiliary coefficients for Thomas algorithm
    const real thomas_coeff_x = (delta_t * denom_chiCm) / (delta_x * delta_x);
    const real thomas_coeff_y = (delta_t * denom_chiCm) / (delta_y * delta_y);

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
    INFOMSG("Starting simulation with OSADI (SERIAL)...\n");

    // Main time loop
    startExecutionTime = omp_get_wtime();

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // Update the active stimuli
        num_active_stimuli = update_and_get_num_active_stimuli(actualTime, stimuli, numberOfStimuli, active_stimuli);

        // =================================================
        //  Compute Reaction and Update ODEs
        // =================================================
        startTime = omp_get_wtime();

        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Get the actual Vm and sV
                actualVm = Vm[idx];
                get_actual_sV(actualsV, sV, idx);

                // Stimulation
                stim = (num_active_stimuli > 0)
                           ? (get_stimulus_value(actualTime, i, j, active_stimuli, num_active_stimuli))
                           : (0.0f);

                // Calculate part of the RHS of the following linear systems with Forward Euler
                partRHS[idx] = delta_t * (stim - compute_dVmdt(actualVm, actualsV));

                // Update state variables
                update_sV(sV, actualsV, actualVm, actualsV, delta_t, idx);
            }
        }

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // =================================================
        //  Calculate Vm at n+1/2 -> Result goes to Vm
        // =================================================
        startTime = omp_get_wtime();

        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                ls_rhs[i] = Vm[idx] + 0.5f * partRHS[idx];
            }

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver_y(Ny, ls_rhs, result, c_prime, d_prime, thomas_coeff_y, Dyy, j, Nx);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                Vm[idx] = result[i];
            }

            elapsedTime1stLS += omp_get_wtime() - startLSTime;
        }

        // =================================================
        //  Calculate Vm at n+1 -> Result goes to Vm
        // =================================================
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                ls_rhs[j] = Vm[idx] + 0.5f * partRHS[idx];
            }

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver_x(Nx, ls_rhs, result, c_prime, d_prime, thomas_coeff_x, Dxx, i);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                Vm[idx] = result[j];
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
    free(partRHS);
    free(c_prime);
    free(d_prime);
    free(ls_rhs);
    free(result);
}