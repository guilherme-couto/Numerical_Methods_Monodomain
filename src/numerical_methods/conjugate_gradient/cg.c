// ==============================================
//         Conjugate Gradient
// ==============================================

#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"

MODIFIERS real dot(const real *a, const real *b, const int N)
{
    real sum = 0.0;
    for (int i = 0; i < N; ++i)
        sum += a[i] * b[i];
    return sum;
}

// Function to apply the diffusion operator to prepare the Ax vector
MODIFIERS void apply_diffusion_operator(const real *x, real *Ax, const real *Dxx, const real *Dyy, const real *Dxy,
                                        const int Nx, const int Ny, const real dx, const real dy, const real dt, const real denom_chiCm, const bool is_aligned)
{
    // const real dx2 = dx * dx;
    // const real dy2 = dy * dy;

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            const int idx = i * Nx + j;

            real lap = select_compute_diffusion_term(is_aligned, x, Dxx, Dyy, Dxy, i, j, Nx, Ny, dx, dy);

            Ax[idx] = x[idx] - dt * lap * denom_chiCm;
        }
    }
}

void runCG(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy)
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
    int i, j, idx;

    // Auxiliary variables for the operations
    real diff_term, stim, actualVm;
    const int max_iterations = 1000;
    const real tolerance = 1e-6;
    const int total_points = Nx * Ny;
    real *actualsV = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    real *RHS = (real *)malloc(total_points * sizeof(real));
    real *r = (real *)malloc(total_points * sizeof(real));
    real *p = (real *)malloc(total_points * sizeof(real));
    real *Ap = (real *)malloc(total_points * sizeof(real));

    // Variables for time measurement
    real startTime = 0.0f;
    real startExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;

    SIMPLEMSG("");
    INFOMSG("Starting simulation with CG (SERIAL)...\n");

    // Main time loop
    startExecutionTime = omp_get_wtime();

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // =================================================
        //  Calculate Approxs. and Update ODEs
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
                stim = get_stimulus_value(actualTime, i, j, stimuli, numberOfStimuli);

                // Prepare the RHS of the following linear system
                RHS[idx] = actualVm + delta_t * (stim - compute_dVmdt(actualVm, actualsV));

                // Update state variables
                update_sV(sV, actualsV, actualVm, actualsV, delta_t, idx);
            }
        }

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // =================================================
        //  Solve Linear System with CG
        // =================================================
        startTime = omp_get_wtime();
        
        // Initialize residual r₀ = b - A·x₀
        apply_diffusion_operator(Vm, Ap, Dxx, Dyy, Dxy, Nx, Ny, delta_x, delta_y, delta_t, denom_chiCm, is_aligned);
        for (int i = 0; i < total_points; ++i)
            r[i] = RHS[i] - Ap[i];

        memcpy(p, r, total_points * sizeof(real));
        real rs_old = dot(r, r, total_points);

        for (int iter = 0; iter < max_iterations; ++iter)
        {
            apply_diffusion_operator(p, Ap, Dxx, Dyy, Dxy, Nx, Ny, delta_x, delta_y, delta_t, denom_chiCm, is_aligned);
            real alphak = rs_old / dot(p, Ap, total_points);

            for (int i = 0; i < total_points; ++i)
                Vm[i] += alphak * p[i];
            for (int i = 0; i < total_points; ++i)
                r[i] -= alphak * Ap[i];

            real rs_new = dot(r, r, total_points);
            if (rs_new < tolerance * tolerance) break;
            

            real betak = rs_new / rs_old;
            for (int i = 0; i < total_points; ++i)
                p[i] = r[i] + betak * p[i];
            rs_old = rs_new;
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
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Free allocated memory
    free(actualsV);
    free(RHS);
    free(r);
    free(p);
    free(Ap);
}