#include "monodomain.h"

int runMonodomainSimulationSerial(const SimulationConfig *config)
{
    // Allocate and populate time array
    real *time_array = (real *)malloc(config->M * sizeof(real));
    initializeTimeArray(time_array, config->M, config->dt);

    // Structure for measurement
    Measurement measurement;
    initializeMeasurement(&measurement);

    // Define which cell model solver to use
    const CellModelSolver *cell_model_solver = get_solver_struct(&config->cell_model);
    if (cell_model_solver == NULL)
    {
        ERRORMSG("Invalid cell model selected.");
        free(time_array);
        return -1;
    }

    // Allocate and initialize state variables arrays
    int total_points = config->Nx * config->Ny;
    real *Vm = (real *)malloc(total_points * sizeof(real));
    real *sV = (real *)malloc(total_points * cell_model_solver->n_state_vars * sizeof(real));
    cell_model_solver->initialize(Vm, sV, config->Nx, config->Ny);

    // Allocate and initialize diffusion coefficients arrays
    real *Dxx = (real *)malloc(total_points * sizeof(real));
    real *Dyy = (real *)malloc(total_points * sizeof(real));
    real *Dxy = (real *)malloc(total_points * sizeof(real));
    if (Dxx == NULL || Dyy == NULL || Dxy == NULL)
    {
        ERRORMSG("Failed to allocate memory for diffusion coefficients.");
        free(time_array);
        free(Vm);
        free(sV);
        return -2;
    }
    initializeProperties(config, Dxx, Dyy, Dxy);

    // Run the simulation based on the selected method
    numerical_method_t run_method = get_numerical_method(&config->method);
    if (run_method == NULL)
    {
        ERRORMSG("Invalid numerical method selected.");
        free(time_array);
        free(Vm);
        free(sV);
        free(Dxx);
        free(Dyy);
        free(Dxy);
        return -3;
    }
    
    // Run the selected method
    run_method(config, &measurement, time_array, cell_model_solver, Vm, sV, Dxx, Dyy, Dxy);

    // Save last frame
    if (config->save_last_frame)
    {
        real startTime = omp_get_wtime();
        static char file_path[MAX_STRING_SIZE];
        snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", config->output_dir, config->M, config->file_extension);
        config->save_function(file_path, Vm, config->Nx, config->Ny, config->dx, config->dy);
        measurement.elapsedSaveFramesTime += (omp_get_wtime() - startTime);
        SUCCESSMSG("Last frame (%.2f ms) saved to %s\n", config->M * config->dt, file_path);
    }

    // Save last state
    if (config->save_last_state)
    {
        real startTime = omp_get_wtime();
        // TODO: save as binary
        static char file_path[MAX_STRING_SIZE];
        snprintf(file_path, MAX_STRING_SIZE, "%s/state_%05d.dat", config->output_dir, config->M);
        measurement.elapsedSaveStateTime += (omp_get_wtime() - startTime);
        SUCCESSMSG("Last state (%.2f ms) saved to %s\n", config->M * config->dt, file_path);
    }

    // Free memory
    free(time_array);
    free(Vm);
    free(sV);
    free(Dxx);
    free(Dyy);
    free(Dxy);

    // Save simulation information
    saveSimulationInfos(config, &measurement);

    return 0;
}

// void baseForRunMonodomain(const SimulationConfig *config)
// {
//     // Unpack configuration parameters
//     const ExecutionMode exec_mode = config->exec_mode;
//     const CellModel cell_model = config->cell_model;
//     const NumericalMethod method = config->method;
//     const real delta_t = config->dt;
//     const real delta_x = config->dx;
//     const real delta_y = config->dy;
//     const real sigma = config->sigma;
//     const real totalTime = config->total_time;
//     const real Lx = config->Lx;
//     const real Ly = config->Ly;
//     const int numberOfStimuli = config->stimulus_count;
//     const Stimulus *stimuli = config->stimuli;
//     const int frameSaveRate = config->frame_save_rate;
//     const int numberOfThreads = config->number_of_threads;
//     const char *pathToSaveData = config->output_dir;
//     const char *pathToRestoreStateFiles = config->path_to_restore_state_files;
//     const char *saveFunctionName = config->save_function_name;
//     const save_function_t save_function = config->save_function;
//     const bool shiftState = config->shift_state;
//     const bool saveFrames = config->save_frames;
//     const bool saveLastFrame = config->save_last_frame;
//     const bool saveLastState = config->save_last_state;
//     const bool measureVelocity = config->measure_velocity;

//     // Number of steps
//     int M = round(totalTime / delta_t); // Number of time steps
//     int Nx = round(Lx / delta_x) + 1;   // Spatial steps in x
//     int Ny = round(Ly / delta_y) + 1;   // Spatial steps in y

//     printf("Nx = %d\n", Nx);
//     printf("Ny = %d\n", Ny);
//     printf("Points in the domain = %d\n", Nx * Ny);

//     // Allocate and populate time array
//     real *time_array = (real *)malloc(M * sizeof(real));
//     initializeTimeArray(time_array, M, delta_t);

//     // Allocate arrays for variables (2D matrices will be flattened) and populate them

// #ifdef TT2

//     real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_r1 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_r2 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_s = (real *)malloc(Nx * Ny * sizeof(real));
//     real *m = (real *)malloc(Nx * Ny * sizeof(real));
//     real *h = (real *)malloc(Nx * Ny * sizeof(real));
//     real *_j = (real *)malloc(Nx * Ny * sizeof(real));
//     real *d = (real *)malloc(Nx * Ny * sizeof(real));
//     real *f = (real *)malloc(Nx * Ny * sizeof(real));
//     real *f2 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *fCaSS = (real *)malloc(Nx * Ny * sizeof(real));
//     real *s = (real *)malloc(Nx * Ny * sizeof(real));
//     real *r = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_i = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_SR = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_SS = (real *)malloc(Nx * Ny * sizeof(real));
//     real *R_prime = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Na_i = (real *)malloc(Nx * Ny * sizeof(real));
//     real *K_i = (real *)malloc(Nx * Ny * sizeof(real));

//     initializeVariableWithValue(Vm, Nx, Ny, Vm_init);
//     initializeVariableWithValue(X_r1, Nx, Ny, X_r1_init);
//     initializeVariableWithValue(X_r2, Nx, Ny, X_r2_init);
//     initializeVariableWithValue(X_s, Nx, Ny, X_s_init);
//     initializeVariableWithValue(m, Nx, Ny, m_init);
//     initializeVariableWithValue(h, Nx, Ny, h_init);
//     initializeVariableWithValue(_j, Nx, Ny, j_init);
//     initializeVariableWithValue(d, Nx, Ny, d_init);
//     initializeVariableWithValue(f, Nx, Ny, f_init);
//     initializeVariableWithValue(f2, Nx, Ny, f2_init);
//     initializeVariableWithValue(fCaSS, Nx, Ny, fCaSS_init);
//     initializeVariableWithValue(s, Nx, Ny, s_init);
//     initializeVariableWithValue(r, Nx, Ny, r_init);
//     initializeVariableWithValue(Ca_i, Nx, Ny, Ca_i_init);
//     initializeVariableWithValue(Ca_SR, Nx, Ny, Ca_SR_init);
//     initializeVariableWithValue(Ca_SS, Nx, Ny, Ca_SS_init);
//     initializeVariableWithValue(R_prime, Nx, Ny, R_prime_init);
//     initializeVariableWithValue(Na_i, Nx, Ny, Na_i_init);
//     initializeVariableWithValue(K_i, Nx, Ny, K_i_init);

// #endif // TT2

// #ifdef RESTORE_STATE

//     printf("\n");
//     printf("Restoring state variables...\n");

//     // Initialize variables with a solution
//     real real_ref_dx = 0.0001f;
//     real real_def_dy = 0.0001f;

//     char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE);

// #ifdef TT2

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r1.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(X_r1, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r1", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(X_r2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r2", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_s.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(X_s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_s", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframem.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(m, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "m", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeh.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(h, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "h", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframej.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(_j, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "j", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframed.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(d, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "d", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(f, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(f2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f2", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframefCaSS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(fCaSS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "fCaSS", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframer.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(r, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "r", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(Ca_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_i", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SR.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(Ca_SR, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SR", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(Ca_SS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SS", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeR_prime.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(R_prime, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "R_prime", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeNa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(Na_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Na_i", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeK_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initializeVariableFromFile(K_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "K_i", real_ref_dx, real_def_dy);

// #endif // TT2

//     free(pathToRestoreStateFiles);

// #endif // RESTORE_STATE

// #ifdef SHIFT_STATE

//     printf("\n");
//     printf("Shifting state variables...\n");

//     // Shift variables
//     real lengthToShift = 0.7f; // Length to shift in cm

// #ifdef TT2

//     shiftVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, Vm_init, "Vm");
//     shiftVariableToLeft(X_r1, Nx, Ny, lengthToShift, delta_x, delta_y, X_r1_init, "X_r1");
//     shiftVariableToLeft(X_r2, Nx, Ny, lengthToShift, delta_x, delta_y, X_r2_init, "X_r2");
//     shiftVariableToLeft(X_s, Nx, Ny, lengthToShift, delta_x, delta_y, X_s_init, "X_s");
//     shiftVariableToLeft(m, Nx, Ny, lengthToShift, delta_x, delta_y, m_init, "m");
//     shiftVariableToLeft(h, Nx, Ny, lengthToShift, delta_x, delta_y, h_init, "h");
//     shiftVariableToLeft(_j, Nx, Ny, lengthToShift, delta_x, delta_y, j_init, "j");
//     shiftVariableToLeft(d, Nx, Ny, lengthToShift, delta_x, delta_y, d_init, "d");
//     shiftVariableToLeft(f, Nx, Ny, lengthToShift, delta_x, delta_y, f_init, "f");
//     shiftVariableToLeft(f2, Nx, Ny, lengthToShift, delta_x, delta_y, f2_init, "f2");
//     shiftVariableToLeft(fCaSS, Nx, Ny, lengthToShift, delta_x, delta_y, fCaSS_init, "fCaSS");
//     shiftVariableToLeft(s, Nx, Ny, lengthToShift, delta_x, delta_y, s_init, "s");
//     shiftVariableToLeft(r, Nx, Ny, lengthToShift, delta_x, delta_y, r_init, "r");
//     shiftVariableToLeft(Ca_i, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_i_init, "Ca_i");
//     shiftVariableToLeft(Ca_SR, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SR_init, "Ca_SR");
//     shiftVariableToLeft(Ca_SS, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SS_init, "Ca_SS");
//     shiftVariableToLeft(R_prime, Nx, Ny, lengthToShift, delta_x, delta_y, R_prime_init, "R_prime");
//     shiftVariableToLeft(Na_i, Nx, Ny, lengthToShift, delta_x, delta_y, Na_i_init, "Na_i");
//     shiftVariableToLeft(K_i, Nx, Ny, lengthToShift, delta_x, delta_y, K_i_init, "K_i");

// #endif // TT2

// #endif // SHIFT_STATE

//     // Populate auxiliary arrays for Thomas algorithm
//     real phi_x = delta_t / (delta_x * delta_x);
//     real phi_y = delta_t / (delta_y * delta_y);
//     real tau = 0.5f; // Used for calculating the explicit diffusion term on the right-hand side of the ADI method
//     real diff_coeff;

// #ifdef TT2

//     diff_coeff = sigma / chi;

// #endif // TT2

//     // Allocate auxiliaries 2D arrays for variables
//     real *partRHS = (real *)malloc(Nx * Ny * sizeof(real));

//     // Auxiliary variables for the time loop
//     int timeStepCounter = 0;
//     real actualTime = 0.0f;
//     int index;

//     // Start measuring the execution time
//     printf("\n");
//     printf("Starting simulation...\n");
//     startExecutionTime = omp_get_wtime();

//     // Main time loop
//     while (timeStepCounter < M)
//     {
//         // Get time step
//         actualTime = time_array[timeStepCounter];

//         // Start measuring time of 1st part
//         startTime = omp_get_wtime();

//         // =================================================
//         //  Calculate Approxs. and Update ODEs           
//         // =================================================
//         real diff_term = 0.0f;
//         for (int i = 0; i < Ny; i++)
//         {
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;

// #ifdef TT2
//                 // TODO: Implement TT2 for MONODOMAIN in SERIAL mode
// #endif // TT2

//             }
//         }

//         // End of the 1st part of the time step
//         finishTime = omp_get_wtime();
//         elapsedTime1stPart += finishTime - startTime;

//         // Start measuring time of 2nd part
//         startTime = omp_get_wtime();

// #if defined(SSIADI) || defined(THETASSIADI)

//         // =================================================
//         //  Calculate Vm at n+1/2 -> Result goes to RHS  
//         //  diffusion implicit in y and explicit in x    
//         // =================================================
//         for (int j = 0; j < Nx; j++)
//         {
//             // Calculate the RHS of the linear system with the explicit diffusion term along x
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 real actualVm = Vm[index];
//                 diff_term = diff_coeff * tau * phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]);
//                 LS_b_y[i] = actualVm + diff_term + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 1st LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

//             // Update with the result
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 RHS[index] = result_y[i];
//             }

//             // End measuring time of 1st LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime1stLS += finishLSTime - startLSTime;
//         }

//         // =================================================
//         //  Calculate Vm at n+1 -> Result goes to Vm     
//         //  diffusion implicit in x and explicit in y    
//         // =================================================
//         for (int i = 0; i < Ny; i++)
//         {
//             // Calculate the RHS of the linear system with the explicit diffusion term along y
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 real actualVm = RHS[index];
//                 diff_term = diff_coeff * tau * phi_y * (RHS[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + RHS[(lim(i + 1, Ny) * Nx + j)]);
//                 LS_b_x[j] = actualVm + diff_term + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 2nd LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

//             // Update with the result
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = result_x[j];
//             }

//             // End measuring time of 2nd LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime2ndLS += finishLSTime - startLSTime;
//         }

// #endif // SSIADI || THETASSIADI

// #ifdef OSADI

//         // =================================================
//         //  Calculate Vm at n+1/2 -> Result goes to Vm     
//         // =================================================
//         for (int j = 0; j < Nx; j++)
//         {
//             // Calculate the RHS of the linear system
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 LS_b_y[i] = Vm[index] + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 1st LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

//             // Update with the result
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = result_y[i];
//             }

//             // End measuring time of 1st LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime1stLS += finishLSTime - startLSTime;
//         }

//         // =================================================
//         //  Calculate Vm at n+1 -> Result goes to Vm       
//         // =================================================
//         for (int i = 0; i < Ny; i++)
//         {
//             // Calculate the RHS of the linear system
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 LS_b_x[j] = Vm[index] + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 2nd LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

//             // Update with the result
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = result_x[j];
//             }

//             // End measuring time of 2nd LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime2ndLS += finishLSTime - startLSTime;
//         }

// #endif // OSADI

// #ifdef FE

//         // ===================
//         //  Update Vm       
//         // ===================
//         for (int i = 0; i < Ny; i++)
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = partRHS[index];
//             }

// #endif // FE

//         // End of the 2nd part of the time step
//         finishTime = omp_get_wtime();
//         elapsedTime2ndPart += finishTime - startTime;

//         // Save frame if needed
//         if (saveFrames)
//         {
//             startSaveFramesTime = omp_get_wtime();

//             // If save frames is true and time step is multiple of frame save rate
//             if (timeStepCounter % frameSaveRate == 0)
//             {
//                 // Save frame
//                 snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
//                 save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
//                 SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
//             }

//             finishSaveFramesTime = omp_get_wtime();
//             elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
//         }

//         // Measure velocity if needed
//         if (measureVelocity)
//         {
//             startMeasureVelocityTime = omp_get_wtime();

//             // Calculate stim velocity
//             if (!stim_velocity_measured)
//             {
//                 real point_potential = 0.0f;
//                 if (!aux_stim_velocity_flag)
//                 {
//                     point_potential = Vm[idx0];
//                     if (cell_model == CELL_MODEL_MV)
//                         point_potential = rescaleVm(point_potential);

//                     if (point_potential > 10.0f)
//                     {
//                         t0 = actualTime;
//                         aux_stim_velocity_flag = true;
//                     }
//                 }
//                 else
//                 {
//                     point_potential = Vm[idx1];
//                     if (cell_model == CELL_MODEL_MV)
//                         point_potential = rescaleVm(point_potential);

//                     if (point_potential > 10.0f)
//                     {
//                         t1 = actualTime;
//                         stim_velocity = (x1 - x0) / (t1 - t0); // cm/ms
//                         stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
//                         stim_velocity_measured = true;
//                         INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, stim_velocity);
//                     }
//                 }
//             }

//             finishMeasureVelocityTime = omp_get_wtime();
//             elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;
//         }

//         // Update time step counter
//         timeStepCounter++;
//     }

//     finishExecutionTime = omp_get_wtime();
//     elapsedExecutionTime += finishExecutionTime - startExecutionTime;

//     printf("Simulation done!\n");
//     printf("\n");

//     // Save last frame
//     snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", pathToSaveData, M, file_extension);
//     save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
//     SUCCESSMSG("Last frame (time %.2f ms) saved to %s\n", M * delta_t, file_path);

//     // Write infos to file
//     srand(time(NULL));
//     int random_number = (rand() % 500) + 1;
//     char infosFilePath[MAX_STRING_SIZE];
//     snprintf(infosFilePath, MAX_STRING_SIZE, "%s/infos%d.txt", pathToSaveData, random_number);
//     FILE *fpInfos = fopen(infosFilePath, "w");
//     fprintf(fpInfos, "EXECUTION TYPE = %s\n", exec_mode);
//     fprintf(fpInfos, "PRECISION = %s\n", REAL_TYPE);
//     fprintf(fpInfos, "PROBLEM = %s\n", PROBLEM);
//     fprintf(fpInfos, "CELL MODEL = %s\n", cell_model);
//     fprintf(fpInfos, "METHOD = %s\n", method);
//     fprintf(fpInfos, "\n");

//     fprintf(fpInfos, "DOMAIN LENGTH IN X = %.4g cm\n", Lx);
//     fprintf(fpInfos, "DOMAIN LENGTH IN Y = %.4g cm\n", Ly);
//     fprintf(fpInfos, "TOTAL TIME = %.4g ms\n", totalTime);
//     fprintf(fpInfos, "\n");

//     fprintf(fpInfos, "delta_t = %.5g ms (%d time steps)\n", delta_t, M);
//     fprintf(fpInfos, "delta_x = %.5g cm (%d um) (%d space steps in x)\n", delta_x, CM_TO_UM(delta_x), Nx);
//     fprintf(fpInfos, "delta_y = %.5g cm (%d um) (%d space steps in y)\n", delta_y, CM_TO_UM(delta_y), Ny);
//     fprintf(fpInfos, "TOTAL POINTS IN DOMAIN = %d\n", Nx * Ny);
//     fprintf(fpInfos, "DIFFUSION COEFFICIENT = %.8g\n", diff_coeff);
//     fprintf(fpInfos, "NUMBER OF STIMULI = %d\n", numberOfStimuli);
//     for (int i = 0; i < numberOfStimuli; i++)
//     {
//         fprintf(fpInfos, "STIMULUS %d: START TIME = %.5g ms\n", i + 1, stimuli[i].start_time);
//     }

//     if (measureVelocity)
//     {
//         fprintf(fpInfos, "\n");
//         fprintf(fpInfos, "STIMULUS S1 VELOCITY = %.4g cm/s\n", stim_velocity);
//     }

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "TIME OF THE FIRST PART = %.5g s\n", elapsedTime1stPart);
//     fprintf(fpInfos, "TIME OF THE SECOND PART = %.5g s\n", elapsedTime2ndPart);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     fprintf(fpInfos, "TIME TO SOLVE THE 1st LINEAR SYSTEM = %.5g s\n", elapsedTime1stLS);
//     fprintf(fpInfos, "TIME TO SOLVE THE 2nd LINEAR SYSTEM = %.5g s\n", elapsedTime2ndLS);

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

//     if (measureVelocity)
//         fprintf(fpInfos, "TIME TO MEASURE VELOCITY = %.5g s\n", elapsedMeasureVelocityTime);

//     if (saveFrames)
//         fprintf(fpInfos, "TIME TO SAVE FRAMES = %.5g s\n", elapsedSaveFramesTime);

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "SIMULATION TOTAL EXECUTION TIME = %.5g s\n", elapsedExecutionTime);

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "OUTPUT DIRECTORY = %s\n", pathToSaveData);
//     fclose(fpInfos);

//     INFOMSG("Simulation total execution time = %.5g s\n", elapsedExecutionTime);
//     SUCCESSMSG("Simulation infos saved to %s\n", infosFilePath);

//     if (saveLastState)
//     {
//         // TODO: save as binary
//         snprintf(file_path, MAX_STRING_SIZE, "%s/state_%05d.dat", pathToSaveData, M);
//     }

//     // Free memory
//     free(time_array);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     free(c_prime_x);
//     free(d_prime_x);
//     free(LS_b_x);
//     free(result_x);
//     free(c_prime_y);
//     free(d_prime_y);
//     free(LS_b_y);
//     free(result_y);

//     free(la_x);
//     free(lb_x);
//     free(lc_x);
//     free(la_y);
//     free(lb_y);
//     free(lc_y);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2)

//     free(RHS);

// #endif // SSIADI || THETASSIADI || THETASSIRK2

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

//     free(partRHS);

// #ifdef TT2

//     free(Vm);
//     free(X_r1);
//     free(X_r2);
//     free(X_s);
//     free(m);
//     free(h);
//     free(j);
//     free(d);
//     free(f);
//     free(f2);
//     free(fCaSS);
//     free(s);
//     free(r);
//     free(R_prime);
//     free(Ca_i);
//     free(Ca_SR);
//     free(Ca_SS);
//     free(Na_i);
//     free(K_i);

// #endif // TT2

// #ifdef MV

//     free(Vm);
//     free(v);
//     free(w);
//     free(s);

// #endif // MV

//     return;
// }

// void runSimulationSerial(const SimulationConfig *config)
// {
//     // Unpack configuration parameters
//     const ExecutionMode exec_mode = config->exec_mode;
//     const CellModel cell_model = config->cell_model;
//     const NumericalMethod method = config->method;
//     const real delta_t = config->dt;
//     const real delta_x = config->dx;
//     const real delta_y = config->dy;
//     const real totalTime = config->total_time;
//     const real Lx = config->Lx;
//     const real Ly = config->Ly;
//     const int numberOfStimuli = config->stimulus_count;
//     const Stimulus *stimuli = config->stimuli;
//     const int frameSaveRate = config->frame_save_rate;
//     const int numberOfThreads = config->number_of_threads;
//     const char *pathToSaveData = config->output_dir;
//     const char *pathToRestoreStateFiles = config->path_to_restore_state_files;
//     const char *saveFunctionName = config->save_function_name;
//     const save_function_t save_function = config->save_function;
//     const bool shiftState = config->shift_state;
//     const bool saveFrames = config->save_frames;
//     const bool saveLastFrame = config->save_last_frame;
//     const bool saveLastState = config->save_last_state;
//     const bool measureVelocity = config->measure_velocity;

//     // Number of steps
//     int M = round(totalTime / delta_t); // Number of time steps
//     int Nx = round(Lx / delta_x) + 1;   // Spatial steps in x
//     printf("Nx = %d\n", Nx);

// #ifndef CABLEEQ

//     int Ny = round(Ly / delta_y) + 1; // Spatial steps in y
//     printf("Ny = %d\n", Ny);
//     printf("Points in the domain = %d\n", Nx * Ny);

// #endif // not CABLEEQ

//     // Allocate and populate time array
//     real *time_array = (real *)malloc(M * sizeof(real));
//     initializeTimeArray(time_array, M, delta_t);

// #ifndef CABLEEQ

//     // Allocate 2D arrays for variables
//     real *partRHS = (real *)malloc(Nx * Ny * sizeof(real));

// #if defined(SSIADI) || defined(THETASSIADI)

//     real *RHS = (real *)malloc(Nx * Ny * sizeof(real));

// #endif // SSIADI || THETASSIADI

// #else // if def CABLEEQ

//     // Allocate 1D arrays for variables
//     real *partRHS = (real *)malloc(Nx * sizeof(real));
//     real *AP = (real *)malloc(M * sizeof(real));

// #endif // not CABLEEQ

// #ifdef MONODOMAIN

// #ifdef AFHN

//     real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
//     real *W = (real *)malloc(Nx * Ny * sizeof(real));

// #endif // AFHN

// #ifdef TT2

//     real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_r1 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_r2 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *X_s = (real *)malloc(Nx * Ny * sizeof(real));
//     real *m = (real *)malloc(Nx * Ny * sizeof(real));
//     real *h = (real *)malloc(Nx * Ny * sizeof(real));
//     real *_j = (real *)malloc(Nx * Ny * sizeof(real));
//     real *d = (real *)malloc(Nx * Ny * sizeof(real));
//     real *f = (real *)malloc(Nx * Ny * sizeof(real));
//     real *f2 = (real *)malloc(Nx * Ny * sizeof(real));
//     real *fCaSS = (real *)malloc(Nx * Ny * sizeof(real));
//     real *s = (real *)malloc(Nx * Ny * sizeof(real));
//     real *r = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_i = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_SR = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Ca_SS = (real *)malloc(Nx * Ny * sizeof(real));
//     real *R_prime = (real *)malloc(Nx * Ny * sizeof(real));
//     real *Na_i = (real *)malloc(Nx * Ny * sizeof(real));
//     real *K_i = (real *)malloc(Nx * Ny * sizeof(real));

// #endif // TT2

// #ifdef MV

//     real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
//     real *v = (real *)malloc(Nx * Ny * sizeof(real));
//     real *w = (real *)malloc(Nx * Ny * sizeof(real));
//     real *s = (real *)malloc(Nx * Ny * sizeof(real));

// #endif // MV

// #endif // MONODOMAIN

// #ifdef CABLEEQ

// #ifdef AFHN

//     real *Vm = (real *)malloc(Nx * sizeof(real));
//     real *W = (real *)malloc(Nx * sizeof(real));

// #endif // AFHN

// #ifdef TT2

//     real *Vm = (real *)malloc(Nx * sizeof(real));
//     real *X_r1 = (real *)malloc(Nx * sizeof(real));
//     real *X_r2 = (real *)malloc(Nx * sizeof(real));
//     real *X_s = (real *)malloc(Nx * sizeof(real));
//     real *m = (real *)malloc(Nx * sizeof(real));
//     real *h = (real *)malloc(Nx * sizeof(real));
//     real *_j = (real *)malloc(Nx * sizeof(real));
//     real *d = (real *)malloc(Nx * sizeof(real));
//     real *f = (real *)malloc(Nx * sizeof(real));
//     real *f2 = (real *)malloc(Nx * sizeof(real));
//     real *fCaSS = (real *)malloc(Nx * sizeof(real));
//     real *s = (real *)malloc(Nx * sizeof(real));
//     real *r = (real *)malloc(Nx * sizeof(real));
//     real *Ca_i = (real *)malloc(Nx * sizeof(real));
//     real *Ca_SR = (real *)malloc(Nx * sizeof(real));
//     real *Ca_SS = (real *)malloc(Nx * sizeof(real));
//     real *R_prime = (real *)malloc(Nx * sizeof(real));
//     real *Na_i = (real *)malloc(Nx * sizeof(real));
//     real *K_i = (real *)malloc(Nx * sizeof(real));

// #endif // TT2

// #ifdef MV

//     real *Vm = (real *)malloc(Nx * sizeof(real));
//     real *v = (real *)malloc(Nx * sizeof(real));
//     real *w = (real *)malloc(Nx * sizeof(real));
//     real *s = (real *)malloc(Nx * sizeof(real));

// #endif // MV

// #endif // CABLEEQ

//     // Initialize variables
// #ifdef MONODOMAIN

// #ifdef AFHN

// #ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

//     initialize2DVariableWithExactSolution(Vm, Nx, Ny, delta_x, delta_y);

// #else // if not def CONVERGENCE_ANALYSIS_FORCING_TERM

//     initialize2DVariableWithValue(Vm, Nx, Ny, Vm_init);

// #endif // CONVERGENCE_ANALYSIS_FORCING_TERM

//     initialize2DVariableWithValue(W, Nx, Ny, W_init);

// #endif // AFHN

// #ifdef TT2

//     initialize2DVariableWithValue(Vm, Nx, Ny, Vm_init);
//     initialize2DVariableWithValue(X_r1, Nx, Ny, X_r1_init);
//     initialize2DVariableWithValue(X_r2, Nx, Ny, X_r2_init);
//     initialize2DVariableWithValue(X_s, Nx, Ny, X_s_init);
//     initialize2DVariableWithValue(m, Nx, Ny, m_init);
//     initialize2DVariableWithValue(h, Nx, Ny, h_init);
//     initialize2DVariableWithValue(_j, Nx, Ny, j_init);
//     initialize2DVariableWithValue(d, Nx, Ny, d_init);
//     initialize2DVariableWithValue(f, Nx, Ny, f_init);
//     initialize2DVariableWithValue(f2, Nx, Ny, f2_init);
//     initialize2DVariableWithValue(fCaSS, Nx, Ny, fCaSS_init);
//     initialize2DVariableWithValue(s, Nx, Ny, s_init);
//     initialize2DVariableWithValue(r, Nx, Ny, r_init);
//     initialize2DVariableWithValue(Ca_i, Nx, Ny, Ca_i_init);
//     initialize2DVariableWithValue(Ca_SR, Nx, Ny, Ca_SR_init);
//     initialize2DVariableWithValue(Ca_SS, Nx, Ny, Ca_SS_init);
//     initialize2DVariableWithValue(R_prime, Nx, Ny, R_prime_init);
//     initialize2DVariableWithValue(Na_i, Nx, Ny, Na_i_init);
//     initialize2DVariableWithValue(K_i, Nx, Ny, K_i_init);

// #endif // TT2

// #ifdef MV

//     initialize2DVariableWithValue(Vm, Nx, Ny, u_init);
//     initialize2DVariableWithValue(v, Nx, Ny, v_init);
//     initialize2DVariableWithValue(w, Nx, Ny, w_init);
//     initialize2DVariableWithValue(s, Nx, Ny, s_init);

// #endif // MV

// #endif // MONODOMAIN

// #ifdef CABLEEQ

// #ifdef AFHN

//     initialize1DVariableWithValue(Vm, Nx, Vm_init);
//     initialize1DVariableWithValue(W, Nx, W_init);

// #endif // AFHN

// #ifdef TT2

//     initialize1DVariableWithValue(Vm, Nx, Vm_init);
//     initialize1DVariableWithValue(X_r1, Nx, X_r1_init);
//     initialize1DVariableWithValue(X_r2, Nx, X_r2_init);
//     initialize1DVariableWithValue(X_s, Nx, X_s_init);
//     initialize1DVariableWithValue(m, Nx, m_init);
//     initialize1DVariableWithValue(h, Nx, h_init);
//     initialize1DVariableWithValue(_j, Nx, j_init);
//     initialize1DVariableWithValue(d, Nx, d_init);
//     initialize1DVariableWithValue(f, Nx, f_init);
//     initialize1DVariableWithValue(f2, Nx, f2_init);
//     initialize1DVariableWithValue(fCaSS, Nx, fCaSS_init);
//     initialize1DVariableWithValue(s, Nx, s_init);
//     initialize1DVariableWithValue(r, Nx, r_init);
//     initialize1DVariableWithValue(Ca_i, Nx, Ca_i_init);
//     initialize1DVariableWithValue(Ca_SR, Nx, Ca_SR_init);
//     initialize1DVariableWithValue(Ca_SS, Nx, Ca_SS_init);
//     initialize1DVariableWithValue(R_prime, Nx, R_prime_init);
//     initialize1DVariableWithValue(Na_i, Nx, Na_i_init);
//     initialize1DVariableWithValue(K_i, Nx, K_i_init);

// #endif // TT2

// #ifdef MV

//     initialize1DVariableWithValue(Vm, Nx, u_init);
//     initialize1DVariableWithValue(v, Nx, v_init);
//     initialize1DVariableWithValue(w, Nx, w_init);
//     initialize1DVariableWithValue(s, Nx, s_init);

// #endif // MV

// #endif // CABLEEQ

// #ifdef RESTORE_STATE

//     printf("\n");
//     printf("Restoring state variables...\n");

//     // Initialize variables with a solution
//     real real_ref_dx = 0.002f;

// #ifndef CABLEEQ

//     real real_def_dy = 0.002f;

// #endif // not CABLEEQ

//     char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE);

// #ifdef CABLEEQ

// #ifdef AFHN

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Vm, Nx, pathToRestoreStateFiles, delta_x, "Vm", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeW.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(W, Nx, pathToRestoreStateFiles, delta_x, "W", real_ref_dx);

// #endif // AFHN

// #ifdef TT2

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Vm, Nx, pathToRestoreStateFiles, delta_x, "Vm", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r1.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(X_r1, Nx, pathToRestoreStateFiles, delta_x, "X_r1", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(X_r2, Nx, pathToRestoreStateFiles, delta_x, "X_r2", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_s.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(X_s, Nx, pathToRestoreStateFiles, delta_x, "X_s", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframem.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(m, Nx, pathToRestoreStateFiles, delta_x, "m", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeh.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(h, Nx, pathToRestoreStateFiles, delta_x, "h", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframej.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(_j, Nx, pathToRestoreStateFiles, delta_x, "j", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframed.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(d, Nx, pathToRestoreStateFiles, delta_x, "d", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(f, Nx, pathToRestoreStateFiles, delta_x, "f", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(f2, Nx, pathToRestoreStateFiles, delta_x, "f2", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframefCaSS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(fCaSS, Nx, pathToRestoreStateFiles, delta_x, "fCaSS", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(s, Nx, pathToRestoreStateFiles, delta_x, "s", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframer.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(r, Nx, pathToRestoreStateFiles, delta_x, "r", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Ca_i, Nx, pathToRestoreStateFiles, delta_x, "Ca_i", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SR.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Ca_SR, Nx, pathToRestoreStateFiles, delta_x, "Ca_SR", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Ca_SS, Nx, pathToRestoreStateFiles, delta_x, "Ca_SS", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeR_prime.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(R_prime, Nx, pathToRestoreStateFiles, delta_x, "R_prime", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeNa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Na_i, Nx, pathToRestoreStateFiles, delta_x, "Na_i", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeK_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(K_i, Nx, pathToRestoreStateFiles, delta_x, "K_i", real_ref_dx);

// #endif // TT2

// #ifdef MV

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(Vm, Nx, pathToRestoreStateFiles, delta_x, "Vm", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframev.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(v, Nx, pathToRestoreStateFiles, delta_x, "v", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframew.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(w, Nx, pathToRestoreStateFiles, delta_x, "w", real_ref_dx);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize1DVariableFromFile(s, Nx, pathToRestoreStateFiles, delta_x, "s", real_ref_dx);

// #endif // MV

// #endif // CABLEEQ

// #ifdef MONODOMAIN

// #ifdef AFHN

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeW.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(W, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "W", real_ref_dx, real_def_dy);

// #endif // AFHN

// #ifdef TT2

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r1.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(X_r1, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r1", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(X_r2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r2", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_s.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(X_s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_s", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframem.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(m, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "m", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeh.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(h, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "h", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframej.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(_j, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "j", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframed.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(d, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "d", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(f, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef2.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(f2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f2", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframefCaSS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(fCaSS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "fCaSS", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframer.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(r, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "r", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Ca_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_i", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SR.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Ca_SR, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SR", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SS.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Ca_SS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SS", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeR_prime.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(R_prime, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "R_prime", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeNa_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Na_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Na_i", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeK_i.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(K_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "K_i", real_ref_dx, real_def_dy);

// #endif // TT2

// #ifdef MV

//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframev.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(v, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "v", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframew.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(w, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "w", real_ref_dx, real_def_dy);
//     snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
//     initialize2DVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);

// #endif // MV

// #endif // MONODOMAIN

//     free(pathToRestoreStateFiles);

// #endif // RESTORE_STATE

// #ifdef SHIFT_STATE

//     printf("\n");
//     printf("Shifting state variables...\n");

//     // Shift variables
//     real lengthToShift = 0.5f;

// #ifdef CABLEEQ

//     shift1DVariableToLeft(Vm, Nx, lengthToShift, delta_x, Vm_init, "Vm");

// #ifdef AFHN

//     shift1DVariableToLeft(W, Nx, lengthToShift, delta_x, W_init, "W");

// #endif // AFHN

// #ifdef TT2

//     shift1DVariableToLeft(X_r1, Nx, lengthToShift, delta_x, X_r1_init, "X_r1");
//     shift1DVariableToLeft(X_r2, Nx, lengthToShift, delta_x, X_r2_init, "X_r2");
//     shift1DVariableToLeft(X_s, Nx, lengthToShift, delta_x, X_s_init, "X_s");
//     shift1DVariableToLeft(m, Nx, lengthToShift, delta_x, m_init, "m");
//     shift1DVariableToLeft(h, Nx, lengthToShift, delta_x, h_init, "h");
//     shift1DVariableToLeft(j, Nx, lengthToShift, delta_x, j_init, "j");
//     shift1DVariableToLeft(d, Nx, lengthToShift, delta_x, d_init, "d");
//     shift1DVariableToLeft(f, Nx, lengthToShift, delta_x, f_init, "f");
//     shift1DVariableToLeft(f2, Nx, lengthToShift, delta_x, f2_init, "f2");
//     shift1DVariableToLeft(fCaSS, Nx, lengthToShift, delta_x, fCaSS_init, "fCaSS");
//     shift1DVariableToLeft(s, Nx, lengthToShift, delta_x, s_init, "s");
//     shift1DVariableToLeft(r, Nx, lengthToShift, delta_x, r_init, "r");
//     shift1DVariableToLeft(Ca_i, Nx, lengthToShift, delta_x, Ca_i_init, "Ca_i");
//     shift1DVariableToLeft(Ca_SR, Nx, lengthToShift, delta_x, Ca_SR_init, "Ca_SR");
//     shift1DVariableToLeft(Ca_SS, Nx, lengthToShift, delta_x, Ca_SS_init, "Ca_SS");
//     shift1DVariableToLeft(R_prime, Nx, lengthToShift, delta_x, R_prime_init, "R_prime");
//     shift1DVariableToLeft(Na_i, Nx, lengthToShift, delta_x, Na_i_init, "Na_i");
//     shift1DVariableToLeft(K_i, Nx, lengthToShift, delta_x, K_i_init, "K_i");

// #endif // TT2

// #ifdef MV

//     shift1DVariableToLeft(Vm, Nx, lengthToShift, delta_x, u_init, "Vm");
//     shift1DVariableToLeft(v, Nx, lengthToShift, delta_x, v_init, "v");
//     shift1DVariableToLeft(w, Nx, lengthToShift, delta_x, w_init, "w");
//     shift1DVariableToLeft(s, Nx, lengthToShift, delta_x, s_init, "s");

// #endif // MV

// #endif // CABLEEQ

// #ifdef MONODOMAIN

// #ifdef AFHN

//     shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, Vm_init, "Vm");
//     shift2DVariableToLeft(W, Nx, Ny, lengthToShift, delta_x, delta_y, W_init, "W");

// #endif // AFHN

// #ifdef TT2

//     shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, Vm_init, "Vm");
//     shift2DVariableToLeft(X_r1, Nx, Ny, lengthToShift, delta_x, delta_y, X_r1_init, "X_r1");
//     shift2DVariableToLeft(X_r2, Nx, Ny, lengthToShift, delta_x, delta_y, X_r2_init, "X_r2");
//     shift2DVariableToLeft(X_s, Nx, Ny, lengthToShift, delta_x, delta_y, X_s_init, "X_s");
//     shift2DVariableToLeft(m, Nx, Ny, lengthToShift, delta_x, delta_y, m_init, "m");
//     shift2DVariableToLeft(h, Nx, Ny, lengthToShift, delta_x, delta_y, h_init, "h");
//     shift2DVariableToLeft(_j, Nx, Ny, lengthToShift, delta_x, delta_y, j_init, "j");
//     shift2DVariableToLeft(d, Nx, Ny, lengthToShift, delta_x, delta_y, d_init, "d");
//     shift2DVariableToLeft(f, Nx, Ny, lengthToShift, delta_x, delta_y, f_init, "f");
//     shift2DVariableToLeft(f2, Nx, Ny, lengthToShift, delta_x, delta_y, f2_init, "f2");
//     shift2DVariableToLeft(fCaSS, Nx, Ny, lengthToShift, delta_x, delta_y, fCaSS_init, "fCaSS");
//     shift2DVariableToLeft(s, Nx, Ny, lengthToShift, delta_x, delta_y, s_init, "s");
//     shift2DVariableToLeft(r, Nx, Ny, lengthToShift, delta_x, delta_y, r_init, "r");
//     shift2DVariableToLeft(Ca_i, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_i_init, "Ca_i");
//     shift2DVariableToLeft(Ca_SR, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SR_init, "Ca_SR");
//     shift2DVariableToLeft(Ca_SS, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SS_init, "Ca_SS");
//     shift2DVariableToLeft(R_prime, Nx, Ny, lengthToShift, delta_x, delta_y, R_prime_init, "R_prime");
//     shift2DVariableToLeft(Na_i, Nx, Ny, lengthToShift, delta_x, delta_y, Na_i_init, "Na_i");
//     shift2DVariableToLeft(K_i, Nx, Ny, lengthToShift, delta_x, delta_y, K_i_init, "K_i");

// #endif // TT2

// #ifdef MV

//     shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, u_init, "Vm");
//     shift2DVariableToLeft(v, Nx, Ny, lengthToShift, delta_x, delta_y, v_init, "v");
//     shift2DVariableToLeft(w, Nx, Ny, lengthToShift, delta_x, delta_y, w_init, "w");
//     shift2DVariableToLeft(s, Nx, Ny, lengthToShift, delta_x, delta_y, s_init, "s");

// #endif // MV

// #endif // MONODOMAIN

// #endif // SHIFT_STATE

//     // Populate auxiliary arrays for Thomas algorithm
//     real phi_x = delta_t / (delta_x * delta_x);
//     real phi_y = delta_t / (delta_y * delta_y);
//     real tau = 0.5f; // Used for calculating the explicit diffusion term on the right-hand side of the ADI method
//     real diff_coeff;

// #ifdef AFHN

//     diff_coeff = sigma / (Cm * chi);

// #endif // AFHN

// #ifdef TT2

//     diff_coeff = sigma / chi;

// #endif // TT2

// #ifdef MV

//     diff_coeff = Dtilde;

// #endif // MV

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     // Auxiliary arrays for Thomas algorithm
//     real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
//     real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
//     real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal

//     real *c_prime_x = (real *)malloc(Nx * sizeof(real));
//     real *d_prime_x = (real *)malloc(Nx * sizeof(real));
//     real *LS_b_x = (real *)malloc(Nx * sizeof(real));
//     real *result_x = (real *)malloc(Nx * sizeof(real));

// #ifndef CABLEEQ

//     real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
//     real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
//     real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal

//     real *c_prime_y = (real *)malloc(Ny * sizeof(real));
//     real *d_prime_y = (real *)malloc(Ny * sizeof(real));
//     real *LS_b_y = (real *)malloc(Ny * sizeof(real));
//     real *result_y = (real *)malloc(Ny * sizeof(real));

// #endif // not CABLEEQ

// #if defined(SSIADI)

//     populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, 0.5f * phi_x * diff_coeff);

// #ifndef CABLEEQ

//     populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, 0.5f * phi_y * diff_coeff);

// #endif // not CABLEEQ

// #endif // SSIADI

// #if defined(THETASSIADI) || defined(THETASSIRK2)

//     populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, THETA * phi_x * diff_coeff);

// #ifndef CABLEEQ

//     populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, THETA * phi_y * diff_coeff);

// #endif // not CABLEEQ

//     tau = 1.0f - THETA;

// #endif // THETASSIADI || THETARK2

// #ifdef OSADI

//     populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, phi_x * diff_coeff);

// #ifndef CABLEEQ

//     populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, phi_y * diff_coeff);

// #endif // not CABLEEQ

// #endif // OSADI

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

// #ifndef CONVERGENCE_ANALYSIS_FORCING_TERM

// #if defined(MONODOMAIN) || defined(CABLEEQ)

// #endif // MONODOMAIN || CABLEEQ

// #endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

// #ifdef CABLEEQ

//     // Choose cell at 0.5 cm to measure Action Potential
//     int APCellIndex = round(0.5f / delta_x);

// #endif // CABLEEQ

//     // Create directories
//     char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE);
//     createDirectories(delta_t, delta_x, delta_y, pathToSaveData);

// #ifdef SAVE_FRAMES

//     // Save frames variables
//     char framesPath[MAX_STRING_SIZE];
//     FILE *fpFrames;
//     snprintf(framesPath, MAX_STRING_SIZE, "%s/frames.txt", pathToSaveData);
//     fpFrames = fopen(framesPath, "w");
//     real startSaveFramesTime, finishSaveFramesTime, elapsedSaveFramesTime = 0.0f;

// #endif // SAVE_FRAMES

// #ifdef MEASURE_VELOCITY

//     // Measure velocity
//     real stim_velocity = 0.0;
//     real x0 = Lx / 3.0f;
//     real x1 = 2.0f * x0;
//     int idx0 = round(x0 / delta_x) + 1;
//     int idx1 = round(x1 / delta_x) + 1;
//     real t0 = 0.0;
//     real t1 = 0.0;
//     bool aux_stim_velocity_flag = false;
//     bool stim_velocity_measured = false;
//     real startMeasureVelocityTime, finishMeasureVelocityTime, elapsedMeasureVelocityTime = 0.0f;

// #endif // MEASURE_VELOCITY

//     // Variables for measuring the execution time
//     real startExecutionTime, finishExecutionTime, elapsedExecutionTime = 0.0f;
//     real startTime, finishTime, elapsedTime1stPart, elapsedTime2ndPart = 0.0f;

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     real startLSTime, finishLSTime, elapsedTime1stLS, elapsedTime2ndLS = 0.0f;

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

//     // Auxiliary variables for the time loop
//     int timeStepCounter = 0;
//     real actualTime = 0.0f;
//     int index;

//     // Start measuring the execution time
//     printf("\n");
//     printf("Starting simulation...\n");
//     startExecutionTime = omp_get_wtime();

// #if defined(LINMONO) || defined(MONODOMAIN)

//     while (timeStepCounter < M)
//     {
//         // Get time step
//         actualTime = time_array[timeStepCounter];

//         // Start measuring time of 1st part
//         startTime = omp_get_wtime();

//         // =================================================
//         //  Calculate Approxs. and Update ODEs           
//         // =================================================
//         real diff_term = 0.0f;
//         for (int i = 0; i < Ny; i++)
//         {
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;

// #ifdef LINMONO

//                 diff_term = diff_coeff * 0.5f * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * Vm[index] + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * Vm[index] + Vm[(lim(i + 1, Ny) * Nx + j)]));
//                 real x = j * delta_x;
//                 real y = i * delta_y;
//                 real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t)) / (chi * Cm);
//                 real reac_term = G * Vm[index] / Cm;
//                 real Vmtilde = Vm[index] + diff_term + (0.5f * delta_t * (for_term - reac_term));

//                 // Preparing part of the RHS of the following linear systems
//                 real reac_tilde_term = G * Vmtilde / Cm;
//                 partRHS[index] = delta_t * (for_term - reac_tilde_term);

// #endif // LINMONO

// #ifdef MONODOMAIN

// #ifdef AFHN

//                 // Calculate the explicit part of the RHS, including the diffusion term in both directions
//                 real actualVm = Vm[index];
//                 real actualW = W[index];
//                 real RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

// #ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

//                 // Calculate forcing term
//                 real x = j * delta_x;
//                 real y = i * delta_y;
//                 real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t), actualW) / (chi * Cm);

// #if defined(SSIADI) || defined(THETASSIADI)

//                 diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
//                 real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (for_term - RHS_Vm_term));

//                 // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
//                 real Wtilde = actualW + (0.5f * delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))));
//                 real RHS_Vmtilde_term = (G * Vmtilde * (1.0f - (Vmtilde / vth)) * (1.0f - (Vmtilde / vp))) + (eta1 * Vmtilde * Wtilde);
//                 partRHS[index] = delta_t * (for_term - RHS_Vmtilde_term);

//                 // Update state variables
//                 W[index] = actualW + delta_t * (eta2 * ((Vmtilde / vp) - (eta3 * Wtilde))); // with RK2 -> sV+1 = sV + dt*R(Vm*, W*)

// #endif // SSIADI || THETASSIADI

// #if defined(OSADI)

//                 // Calculate part of the RHS of the following linear systems with Forward Euler
//                 partRHS[index] = delta_t * (for_term - RHS_Vm_term);

//                 // Update state variables
//                 W[index] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))); // with Forward Euler -> sV+1 = sV + dt*R(Vmn, sV)

// #endif // OSADI

// #else // if not def CONVERGENCE_ANALYSIS_FORCING_TERM

//                 // Stimulation
//                 real stim = 0.0f;

// #pragma unroll
//                 for (int si = 0; si < numberOfStimuli; si++)
//                 {
//                     if (actualTime >= stimuli[si].start_time && actualTime <= stimuli[si].start_time + stimuli[si].duration && j >= stimuli[si].x_discretized.min && j <= stimuli[si].x_discretized.max && i >= stimuli[si].y_discretized.min && i <= stimuli[si].y_discretized.max)
//                     {
//                         stim = stimuli[si].amplitude;
//                         break;
//                     }
//                 }

// #if defined(SSIADI) || defined(THETASSIADI)

//                 // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, sV)
//                 diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
//                 real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

//                 // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
//                 real Wtilde = actualW + (0.5f * delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))));
//                 real RHS_Vmtilde_term = (G * Vmtilde * (1.0f - (Vmtilde / vth)) * (1.0f - (Vmtilde / vp))) + (eta1 * Vmtilde * Wtilde);
//                 partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

//                 // Update state variables
//                 W[index] = actualW + delta_t * (eta2 * ((Vmtilde / vp) - (eta3 * Wtilde))); // with RK2 -> sV+1 = sV + dt*R(Vm*, W*)

// #endif // SSIADI || THETASSIADI

// #if defined(OSADI)

//                 // Calculate part of the RHS of the following linear systems with Forward Euler
//                 partRHS[index] = delta_t * (stim - RHS_Vm_term);

//                 // Update state variables
//                 W[index] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))); // with Forward Euler -> sV+1 = sV + dt*R(Vmn, sV)

// #endif // OSADI

// #ifdef FE

//                 // Update variables explicitly
//                 diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
//                 partRHS[index] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

//                 W[index] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW)));

// #endif // FE

// #endif // CONVERGENCE_ANALYSIS_FORCING_TERM

// #endif // AFHN

// #ifdef TT2
//                 // TODO: Implement TT2 for MONODOMAIN in SERIAL mode
// #endif // TT2

// #ifdef MV

//                 // Calculate the explicit part of the RHS, including the diffusion term in both directions
//                 real actualVm = Vm[index];
//                 real actualv = v[index];
//                 real actualw = w[index];
//                 real actuals = s[index];

//                 real stim = 0.0f;

// #pragma unroll
//                 for (int si = 0; si < numberOfStimuli; si++)
//                 {
//                     if (actualTime >= stimuli[si].start_time && actualTime <= stimuli[si].start_time + stimuli[si].duration && j >= stimuli[si].x_discretized.min && j <= stimuli[si].x_discretized.max && i >= stimuli[si].y_discretized.min && i <= stimuli[si].y_discretized.max)
//                     {
//                         stim = stimuli[si].amplitude;
//                         break;
//                     }
//                 }

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2)

//                 diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || FE

//                 // Calculate RHS of the equations
//                 // Auxiliary variables
//                 real Htheta_w = (actualVm - theta_w > 0.0f) ? 1.0f : 0.0f;
//                 real Htheta_o = (actualVm - theta_o > 0.0f) ? 1.0f : 0.0f;
//                 real Htheta_v = (actualVm - theta_v > 0.0f) ? 1.0f : 0.0f;
//                 real Htheta_vminus = (actualVm - theta_vminus > 0.0f) ? 1.0f : 0.0f;

//                 real tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
//                 real tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (actualVm - u_so))) * 0.5f;
//                 real tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

//                 // Currents
//                 real J_fi = -actualv * Htheta_v * (actualVm - theta_v) * (u_u - actualVm) / tau_fi;
//                 real J_so = ((actualVm - u_o) * (1.0f - Htheta_w) / tau_o) + (Htheta_w / tau_so);
//                 real J_si = -Htheta_w * actualw * actuals / tau_si;

//                 // RHS of the state variables
//                 real RHS_Vm_term = J_fi + J_so + J_si;

// #if defined(SSIADI) || defined(THETASSIADI)

//                 // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
//                 real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

//                 // Auxiliary variables for Rush-Larsen or Forward Euler
//                 real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
//                 real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
//                 real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

//                 real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
//                 real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
//                 real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
//                 real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

//                 real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
//                 real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

//                 // Calculate approximations with Rush-Larsen or Forward Euler using half time step
//                 real vtilde, wtilde, stilde;
//                 (tau_v_RL > 1.0e-10)
//                     ? (vtilde = v_inf_RL - (v_inf_RL - actualv) * exp(-0.5f * delta_t / tau_v_RL))
//                     : (vtilde = actualv + 0.5f * delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

//                 (tau_w_RL > 1.0e-10)
//                     ? (wtilde = w_inf_RL - (w_inf_RL - actualw) * exp(-0.5f * delta_t / tau_w_RL))
//                     : (wtilde = actualw + 0.5f * delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

//                 (tau_s > 1.0e-10)
//                     ? (stilde = s_inf_RL - (s_inf_RL - actuals) * exp(-0.5f * delta_t / tau_s))
//                     : (stilde = actuals + 0.5f * delta_t * (s_inf_RL - actuals) / tau_s);

//                 // Calculate RHS of the equations with approximations
//                 // Auxiliary variables
//                 Htheta_w = (Vmtilde - theta_w > 0.0f) ? 1.0f : 0.0f;
//                 Htheta_o = (Vmtilde - theta_o > 0.0f) ? 1.0f : 0.0f;
//                 Htheta_v = (Vmtilde - theta_v > 0.0f) ? 1.0f : 0.0f;
//                 Htheta_vminus = (Vmtilde - theta_vminus > 0.0f) ? 1.0f : 0.0f;

//                 tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
//                 tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so))) * 0.5f;
//                 tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

//                 // Currents
//                 real J_fi_tilde = -actualv * Htheta_v * (Vmtilde - theta_v) * (u_u - Vmtilde) / tau_fi;
//                 real J_so_tilde = ((Vmtilde - u_o) * (1.0f - Htheta_w) / ((1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2)) + (Htheta_w / (tau_so1 + (((tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so)))) * 0.5f)));
//                 real J_si_tilde = -Htheta_w * wtilde * stilde / tau_si;

//                 // Update partRHS
//                 real RHS_Vmtilde_term = J_fi_tilde + J_so_tilde + J_si_tilde;
//                 partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

//                 // Update auxiliary variables for Rush-Larsen or Forward Euler with approximations
//                 v_inf = ((Vmtilde < theta_vminus) ? 1.0f : 0.0f);
//                 tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
//                 v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

//                 w_inf = ((1.0f - Htheta_o) * (1.0f - (Vmtilde / tau_winf)) + Htheta_o * w_infstar);
//                 tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (Vmtilde - u_wminus))) * 0.5f;
//                 tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
//                 w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

//                 tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
//                 s_inf_RL = (1.0f + tanh(k_s * (Vmtilde - u_s))) * 0.5f;

//                 // Update state variables with Rush-Larsen or RK2 (Heun's method) using approximations
//                 (tau_v_RL > 1.0e-10)
//                     ? (v[index] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
//                     : (v[index] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - vtilde) / tau_vminus - Htheta_v * vtilde / tau_vplus);

//                 (tau_w_RL > 1.0e-10)
//                     ? (w[index] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
//                     : (w[index] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - wtilde) / tau_wminus - Htheta_w * wtilde / tau_wplus);

//                 (tau_s > 1.0e-10)
//                     ? (s[index] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
//                     : (s[index] = actuals + delta_t * (s_inf_RL - stilde) / tau_s);

// #endif // SSIADI || THETASSIADI

// #ifdef OSADI

//                 // Update partRHS
//                 partRHS[index] = delta_t * (stim - RHS_Vm_term);

//                 // Auxiliary variables for Rush-Larsen or Forward Euler
//                 real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
//                 real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
//                 real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

//                 real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
//                 real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
//                 real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
//                 real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

//                 real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
//                 real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

//                 // Update state variables with Rush-Larsen or Forward Euler
//                 (tau_v_RL > 1.0e-10)
//                     ? (v[index] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
//                     : (v[index] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

//                 (tau_w_RL > 1.0e-10)
//                     ? (w[index] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
//                     : (w[index] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

//                 (tau_s > 1.0e-10)
//                     ? (s[index] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
//                     : (s[index] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

// #endif // OSADI

// #ifdef FE

//                 // Update partRHS (auxiliary variable) with Forward Euler
//                 diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
//                 partRHS[index] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

//                 // Auxiliary variables for Rush-Larsen or Forward Euler
//                 real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
//                 real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
//                 real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

//                 real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
//                 real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
//                 real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
//                 real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

//                 real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
//                 real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

//                 // Update state variables with Rush-Larsen or Forward Euler
//                 (tau_v_RL > 1.0e-10)
//                     ? (v[index] = v_inf_RL + (actualv - v_inf_RL) * exp(-delta_t / tau_v_RL))
//                     : (v[index] = actualv + delta_t * ((1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus));

//                 (tau_w_RL > 1.0e-10)
//                     ? (w[index] = w_inf_RL + (actualw - w_inf_RL) * exp(-delta_t / tau_w_RL))
//                     : (w[index] = actualw + delta_t * ((1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus));

//                 (tau_s > 1.0e-10)
//                     ? (s[index] = s_inf_RL + (actuals - s_inf_RL) * exp(-delta_t / tau_s))
//                     : (s[index] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

// #endif // FE

// #endif // MV

// #endif // MONODOMAIN
//             }
//         }

//         // End of the 1st part of the time step
//         finishTime = omp_get_wtime();
//         elapsedTime1stPart += finishTime - startTime;

//         // Start measuring time of 2nd part
//         startTime = omp_get_wtime();

// #if defined(SSIADI) || defined(THETASSIADI)

//         // =================================================
//         //  Calculate Vm at n+1/2 -> Result goes to RHS  
//         //  diffusion implicit in y and explicit in x    
//         // =================================================
//         for (int j = 0; j < Nx; j++)
//         {
//             // Calculate the RHS of the linear system with the explicit diffusion term along x
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 real actualVm = Vm[index];
//                 diff_term = diff_coeff * tau * phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]);
//                 LS_b_y[i] = actualVm + diff_term + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 1st LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

//             // Update with the result
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 RHS[index] = result_y[i];
//             }

//             // End measuring time of 1st LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime1stLS += finishLSTime - startLSTime;
//         }

//         // =================================================
//         //  Calculate Vm at n+1 -> Result goes to Vm     
//         //  diffusion implicit in x and explicit in y    
//         // =================================================
//         for (int i = 0; i < Ny; i++)
//         {
//             // Calculate the RHS of the linear system with the explicit diffusion term along y
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 real actualVm = RHS[index];
//                 diff_term = diff_coeff * tau * phi_y * (RHS[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + RHS[(lim(i + 1, Ny) * Nx + j)]);
//                 LS_b_x[j] = actualVm + diff_term + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 2nd LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

//             // Update with the result
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = result_x[j];
//             }

//             // End measuring time of 2nd LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime2ndLS += finishLSTime - startLSTime;
//         }

// #endif // SSIADI || THETASSIADI

// #ifdef OSADI

//         // =================================================
//         //  Calculate Vm at n+1/2 -> Result goes to Vm     
//         // =================================================
//         for (int j = 0; j < Nx; j++)
//         {
//             // Calculate the RHS of the linear system
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 LS_b_y[i] = Vm[index] + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 1st LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

//             // Update with the result
//             for (int i = 0; i < Ny; i++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = result_y[i];
//             }

//             // End measuring time of 1st LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime1stLS += finishLSTime - startLSTime;
//         }

//         // =================================================
//         //  Calculate Vm at n+1 -> Result goes to Vm       
//         // =================================================
//         for (int i = 0; i < Ny; i++)
//         {
//             // Calculate the RHS of the linear system
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 LS_b_x[j] = Vm[index] + 0.5f * partRHS[index];
//             }

//             // Start measuring time of 2nd LS
//             startLSTime = omp_get_wtime();

//             // Solve the linear system
//             tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

//             // Update with the result
//             for (int j = 0; j < Nx; j++)
//             {
//                 Vm[index] = result_x[j];
//             }

//             // End measuring time of 2nd LS
//             finishLSTime = omp_get_wtime();
//             elapsedTime2ndLS += finishLSTime - startLSTime;
//         }

// #endif // OSADI

// #ifdef FE

//         // ===================
//         //  Update Vm       
//         // ===================
//         for (int i = 0; i < Ny; i++)
//             for (int j = 0; j < Nx; j++)
//             {
//                 index = i * Nx + j;
//                 Vm[index] = partRHS[index];
//             }

// #endif // FE

//         // End of the 2nd part of the time step
//         finishTime = omp_get_wtime();
//         elapsedTime2ndPart += finishTime - startTime;

// #ifdef SAVE_FRAMES

//         startSaveFramesTime = omp_get_wtime();

//         // If save frames is true and time step is multiple of frame save rate
//         if (timeStepCounter % frameSaveRate == 0)
//         {
//             // Save frame
//             fprintf(fpFrames, "%lf\n", actualTime);
//             saveFrame(fpFrames, Vm, Nx, Ny);
//             SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
//         }

//         finishSaveFramesTime = omp_get_wtime();
//         elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;

// #endif // SAVE_FRAMES

// #ifndef CONVERGENCE_ANALYSIS_FORCING_TERM

// #ifdef MEASURE_VELOCITY

//         startMeasureVelocityTime = omp_get_wtime();

//         // Calculate stim velocity
//         if (!stim_velocity_measured)
//         {
//             real point_potential = 0.0f;
//             if (!aux_stim_velocity_flag)
//             {
//                 point_potential = Vm[idx0];
// #ifdef MV
//                 point_potential = rescaleVm(point_potential);
// #endif // MV
//                 if (point_potential > 10.0f)
//                 {
//                     t0 = actualTime;
//                     aux_stim_velocity_flag = true;
//                 }
//             }
//             else
//             {
//                 point_potential = Vm[idx1];
// #ifdef MV
//                 point_potential = rescaleVm(point_potential);
// #endif // MV
//                 if (point_potential > 10.0f)
//                 {
//                     t1 = actualTime;
//                     stim_velocity = (x1 - x0) / (t1 - t0); // cm/ms
//                     stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
//                     stim_velocity_measured = true;
//                     INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, stim_velocity);
//                 }
//             }
//         }

//         finishMeasureVelocityTime = omp_get_wtime();
//         elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;

// #endif // MEASURE_VELOCITY

// #endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

//         // Update time step counter
//         timeStepCounter++;
//     }

// #endif // LINMONO || MONODOMAIN

// #ifdef CABLEEQ

//     if (strcmp(method, "theta-RK2") == 0)
//     {
//         while (timeStepCounter < M)
//         {
//             // Get time step
//             actualTime = time_array[timeStepCounter];

//             // Get info for Action Potential
//             AP[timeStepCounter] = Vm[APCellIndex];

//             // =================================================
//             //  Calcula Approx.                              
//             // =================================================
//             real diff_term = 0.0f;
//             for (int i = 0; i < Nx; i++)
//             {

// #ifdef AFHN

//                 real actualVm = Vm[i];
//                 real actualW = W[i];

//                 diff_term = diff_coeff * phi_x * (Vm[lim(i - 1, Nx)] - 2.0f * actualVm + Vm[lim(i + 1, Nx)]);
//                 real RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

//                 // Stimulation
//                 real stim = 0.0f;
//                 for (int si = 0; si < numberOfStimuli; si++)
//                 {
//                     if (actualTime >= stimuli[si].start_time && actualTime <= stimuli[si].start_time + stimuli[si].duration && i >= stimuli[si].x_discretized.min && i <= stimuli[si].x_discretized.max)
//                     {
//                         stim = stimuli[si].amplitude;
//                         break;
//                     }
//                 }

//                 // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
//                 real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

//                 // Calculate W approximation
//                 real RHS_W_term = (eta2 * ((actualVm / vp) - (eta3 * actualW)));
//                 real Wtilde = actualW + (0.5f * delta_t * RHS_W_term);

//                 // Preparing part of the RHS of the following linear systems
//                 real RHS_Vmtilde_term = (G * Vmtilde * (1.0f - (Vmtilde / vth)) * (1.0f - (Vmtilde / vp))) + (eta1 * Vmtilde * Wtilde);
//                 partRHS[i] = delta_t * (stim - RHS_Vmtilde_term);

//                 // Update sV+1 with RK2 -> sV+1 = sV + dt*R(Vm*, W*)
//                 W[i] = actualW + delta_t * (eta2 * ((Vmtilde / vp) - (eta3 * Wtilde)));

// #endif // AFHN

// #ifdef TT2

//                 real actualVm = Vm[i];

//                 diff_term = diff_coeff * phi_x * (Vm[lim(i - 1, Nx)] - 2.0f * actualVm + Vm[lim(i + 1, Nx)]);

//                 // State variables
//                 real actualX_r1 = X_r1[i];
//                 real actualX_r2 = X_r2[i];
//                 real actualX_s = X_s[i];
//                 real actualm = m[i];
//                 real actualh = h[i];
//                 real actualj = j[i];
//                 real actuald = d[i];
//                 real actualf = f[i];
//                 real actualf2 = f2[i];
//                 real actualfCaSS = fCaSS[i];
//                 real actuals = s[i];
//                 real actualr = r[i];
//                 real actualCa_i = Ca_i[i];
//                 real actualCa_SR = Ca_SR[i];
//                 real actualCa_SS = Ca_SS[i];
//                 real actualR_prime = R_prime[i];
//                 real actualNa_i = Na_i[i];
//                 real actualK_i = K_i[i];

//                 // Auxiliary variables
//                 real VmENa = actualVm - (RTONF * log(Na_o / actualNa_i));
//                 real E_K = (RTONF * log(K_o / actualK_i));
//                 real VmEK = actualVm - E_K;
//                 real alpha_K1 = 0.1f / (1.0f + exp(0.06f * (actualVm - E_K - 200.0f)));
//                 real beta_K1 = (3.0f * exp(0.0002f * (actualVm - E_K + 100.0f)) + exp(0.1f * (actualVm - E_K - 10.0f))) / (1.0f + exp(-0.5f * (actualVm - E_K)));
//                 real E_Ks = (RTONF * log((K_o + p_KNa * Na_o) / (actualK_i + p_KNa * actualNa_i)));
//                 real E_Ca = 0.5f * RTONF * log(Ca_o / actualCa_i);

//                 // Currents
//                 real INa = G_Na * (actualm * actualm * actualm) * actualh * actualj * VmENa;
//                 real IbNa = G_bNa * VmENa;
//                 real IK1 = G_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * VmEK;
//                 real Ito = G_to * actualr * actuals * VmEK;
//                 real IKr = G_Kr * sqrt(K_o / 5.4f) * actualX_r1 * actualX_r2 * VmEK;
//                 real IKs = G_Ks * actualX_s * actualX_s * (actualVm - E_Ks);
//                 real ICaL; // !!!
//                 (actualVm < 15.0f - 1.0e-5f)
//                     ? (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 4.0f * (actualVm - 15.0f) * (F * F) * (0.25f * actualCa_SS * exp(2.0f * (actualVm - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualVm - 15.0f) * FONRT) - 1.0f)))
//                     : (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 2.0f * F * (0.25f * actualCa_SS - Ca_o));
//                 real INaK = ((((p_KNa * K_o) / (K_o + K_mK)) * actualNa_i) / (actualNa_i + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualVm * FONRT))) + (0.0353f * exp(((-actualVm) * FONRT))));
//                 real INaCa; // !!!
//                 INaCa = (k_NaCa * ((exp((gamma_I_NaCa * actualVm * FONRT)) * (actualNa_i * actualNa_i * actualNa_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualVm * FONRT)) * (Na_o * Na_o * Na_o) * actualCa_i * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualVm * FONRT)))));
//                 real IpCa = (G_pCa * actualCa_i) / (K_pCa + actualCa_i);
//                 real IpK = (G_pK * VmEK) / (1.0f + exp((25.0f - actualVm) / 5.98f));
//                 real IbCa = G_bCa * (actualVm - E_Ca);

//                 // RHS of the main equation
//                 real RHS_Vm_term = INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;

//                 real stim = 0.0f;
//                 for (int si = 0; si < numberOfStimuli; si++)
//                 {
//                     if (actualTime >= stimuli[si].start_time && actualTime <= stimuli[si].start_time + stimuli[si].duration && i >= stimuli[si].x_discretized.min && i <= stimuli[si].x_discretized.max)
//                     {
//                         stim = stimuli[si].amplitude;
//                         break;
//                     }
//                 }

//                 // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
//                 real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

//                 // Preparing part of the RHS of the following linear systems
//                 // Calculate approximation for state variables
//                 // Rush-Larsen method - auxiliary variables
//                 real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualVm) / 7.0f));
//                 real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualVm) / 10.0f));
//                 real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualVm) / 11.5f));

//                 real X_r2_inf = 1.0f / (1.0f + exp((actualVm + 88.0f) / 24.0f));
//                 real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualVm) / 20.0f));
//                 real beta_X_r2 = 1.12f / (1.0f + exp((actualVm - 60.0f) / 20.0f));

//                 real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualVm) / 14.0f));
//                 real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualVm) / 6.0f));
//                 real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualVm) / 15.0f));

//                 real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualVm) / 9.03f)) * (1.0f + exp((-56.86f - actualVm) / 9.03f)));
//                 real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualVm) / 5.0f));
//                 real beta_m = 0.1f / (1.0f + exp((actualVm + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVm - 50.0f) / 200.0f)));

//                 real h_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
//                 real alpha_h;
//                 (actualVm < -40.0f)
//                     ? (alpha_h = 0.057f * exp(-(actualVm + 80.0f) / 6.8f))
//                     : (alpha_h = 0.0f);
//                 real beta_h;
//                 (actualVm < -40.0f)
//                     ? (beta_h = 2.7f * exp(0.079f * actualVm) + 3.1f * 1.0e5f * exp(0.3485f * actualVm))
//                     : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualVm + 10.66f) / -11.1f))));

//                 real j_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
//                 real alpha_j;
//                 (actualVm < -40.0f)
//                     ? (alpha_j = ((-25428.0f * exp(0.2444f * actualVm) - (6.948e-6f * exp((-0.04391f) * actualVm))) * (actualVm + 37.78f)) / (1.0f + exp(0.311f * (actualVm + 79.23f))))
//                     : (alpha_j = 0.0f);
//                 real beta_j;
//                 (actualVm < -40.0f)
//                     ? (beta_j = (0.02424f * exp(-0.01052f * actualVm)) / (1.0f + exp(-0.1378f * (actualVm + 40.14f))))
//                     : (beta_j = (0.6f * exp(0.057f * actualVm)) / (1.0f + exp(-0.1f * (actualVm + 32.0f))));

//                 real inf = 1.0f / (1.0f + exp((-8.0f - actualVm) / 7.5f));
//                 real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualVm) / 13.0f)) + 0.25f;
//                 real beta_d = 1.4f / (1.0f + exp((actualVm + 5.0f) / 5.0f));
//                 real gamma_d = 1.0f / (1.0f + exp((50.0f - actualVm) / 20.0f));

//                 real f_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 7.0f));
//                 real alpha_f = 1102.5f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 225.0f);
//                 real beta_f = 200.0f / (1.0f + exp((13.0f - actualVm) / 10.0f));
//                 real gamma_f = 180.0f / (1.0f + exp((actualVm + 30.0f) / 10.0f)) + 20.0f;

//                 real f2_inf = 0.67f / (1.0f + exp((actualVm + 35.0f) / 7.0f)) + 0.33f;
//                 real alpha_f2; // !!!
//                 alpha_f2 = 562.0f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 240.0f);
//                 real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualVm) / 10.0f));
//                 real gamma_f2; // !!!
//                 gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualVm) / 10.0f));

//                 real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
//                 real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

// #if defined(EPI) || defined(MCELL)

//                 real s_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 5.0f));
//                 real tau_s = 85.0f * exp(-(actualVm + 45.0f) * (actualVm + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVm - 20.0f) / 5.0f)) + 3.0f;

// #endif // EPI || MCELL
// #ifdef ENDO

//                 real s_inf = 1.0f / (1.0f + exp((actualVm + 28.0f) / 5.0f));
//                 real tau_s = 1000.0f * exp(-(actualVm + 67.0f) * (actualVm + 67.0f) / 1000.0f) + 8.0f;

// #endif // ENDO

//                 real r_inf = 1.0f / (1.0f + exp((20.0f - actualVm) / 6.0f));
//                 real tau_r = 9.5f * exp(-(actualVm + 40.0f) * (actualVm + 40.0f) / 1800.0f) + 0.8f;

//                 // Explicit method - auxiliary variables
//                 real Ileak = V_leak * (actualCa_SR - actualCa_i);
//                 real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
//                 real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
//                 real k1 = k1_prime / k_CaSR;
//                 real k2 = k2_prime * k_CaSR;
//                 real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
//                 real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
//                 real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
//                 real Ca_i_bufC; // !!!
//                 Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
//                 real Ca_SR_bufSR; // !!!
//                 Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
//                 real Ca_SS_bufSS; // !!!
//                 Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

//                 // Explicit method - RHS of the state variables
//                 real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
//                 real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
//                 real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
//                 real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
//                 real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
//                 real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

//                 // Rush-Larsen method - update approximations
//                 real X_r1tilde = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-(0.5f * delta_t) / (alpha_X_r1 * beta_X_r1));
//                 real X_r2tilde = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-(0.5f * delta_t) / (alpha_X_r2 * beta_X_r2));
//                 real X_stilde = X_s_inf - (X_s_inf - actualX_s) * exp(-(0.5f * delta_t) / (alpha_X_s * beta_X_s + 80.0f));
//                 real mtilde = m_inf - (m_inf - actualm) * exp(-(0.5f * delta_t) / (alpha_m * beta_m));
//                 real htilde = h_inf - (h_inf - actualh) * exp(-(0.5f * delta_t) * (alpha_h + beta_h));
//                 real jtilde = j_inf - (j_inf - actualj) * exp(-(0.5f * delta_t) * (alpha_j + beta_j));
//                 real dtilde = inf - (inf - actuald) * exp(-(0.5f * delta_t) / (alpha_d * beta_d + gamma_d));
//                 real ftilde = f_inf - (f_inf - actualf) * exp(-(0.5f * delta_t) / (alpha_f + beta_f + gamma_f));
//                 real f2tilde = f2_inf - (f2_inf - actualf2) * exp(-(0.5f * delta_t) / (alpha_f2 + beta_f2 + gamma_f2));
//                 real fCaSStilde = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-(0.5f * delta_t) / tau_fCaSS);
//                 real stilde = s_inf - (s_inf - actuals) * exp(-(0.5f * delta_t) / tau_s);
//                 real rtilde = r_inf - (r_inf - actualr) * exp(-(0.5f * delta_t) / tau_r);

//                 // Explicit method - update approximations
//                 real R_primetilde = actualR_prime + (0.5f * delta_t * RHS_R_prime_term);
//                 real Ca_itilde = actualCa_i + (0.5f * delta_t * RHS_Ca_i_term);
//                 real Ca_SRtilde = actualCa_SR + (0.5f * delta_t * RHS_Ca_SR_term);
//                 real Ca_SStilde = actualCa_SS + (0.5f * delta_t * RHS_Ca_SS_term);
//                 real Na_itilde = actualNa_i + (0.5f * delta_t * RHS_Na_i_term);
//                 real K_itilde = actualK_i + (0.5f * delta_t * RHS_K_i_term);

//                 // Auxiliary variables with Vmtilde
//                 real VmENatilde = Vmtilde - (RTONF * log(Na_o / Na_itilde));
//                 real E_Ktilde = (RTONF * log(K_o / K_itilde));
//                 real VmEKtilde = Vmtilde - E_Ktilde;
//                 real alpha_K1tilde = 0.1f / (1.0f + exp(0.06f * (Vmtilde - E_Ktilde - 200.0f)));
//                 real beta_K1tilde = (3.0f * exp(0.0002f * (Vmtilde - E_Ktilde + 100.0f)) + exp(0.1f * (Vmtilde - E_Ktilde - 10.0f))) / (1.0f + exp(-0.5f * (Vmtilde - E_Ktilde)));
//                 real E_Kstilde = (RTONF * log((K_o + p_KNa * Na_o) / (K_itilde + p_KNa * Na_itilde)));
//                 real E_Catilde = 0.5f * RTONF * log(Ca_o / Ca_itilde);

//                 // Currents with Vmtilde
//                 real INatilde = G_Na * (mtilde * mtilde * mtilde) * htilde * jtilde * VmENatilde;
//                 real IbNatilde = G_bNa * VmENatilde;
//                 real IK1tilde = G_K1 * (alpha_K1tilde / (alpha_K1tilde + beta_K1tilde)) * VmEKtilde;
//                 real Itotilde = G_to * rtilde * stilde * VmEKtilde;
//                 real IKrtilde = G_Kr * sqrt(K_o / 5.4f) * X_r1tilde * X_r2tilde * VmEKtilde;
//                 real IKstilde = G_Ks * X_stilde * X_stilde * (Vmtilde - E_Kstilde);
//                 real ICaLtilde; // !!!
//                 (Vmtilde < 15.0f - 1.0e-5f)
//                     ? (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 4.0f * (Vmtilde - 15.0f) * (F * F) * (0.25f * Ca_SStilde * exp(2.0f * (Vmtilde - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (Vmtilde - 15.0f) * FONRT) - 1.0f)))
//                     : (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 2.0f * F * (0.25f * Ca_SStilde - Ca_o));
//                 real INaKtilde = ((((p_KNa * K_o) / (K_o + K_mK)) * Na_itilde) / (Na_itilde + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * Vmtilde * FONRT))) + (0.0353f * exp(((-Vmtilde) * FONRT))));
//                 real INaCatilde; // !!!
//                 INaCatilde = (k_NaCa * ((exp((gamma_I_NaCa * Vmtilde * FONRT)) * (Na_itilde * Na_itilde * Na_itilde) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * Vmtilde * FONRT)) * (Na_o * Na_o * Na_o) * Ca_itilde * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*Vmtilde * FONRT)))));
//                 real IpCatilde = (G_pCa * Ca_itilde) / (K_pCa + Ca_itilde);
//                 real IpKtilde = (G_pK * VmEKtilde) / (1.0f + exp((25.0f - Vmtilde) / 5.98f));
//                 real IbCatilde = G_bCa * (Vmtilde - E_Catilde);

//                 // part of RHS of the main equation with Vmtilde
//                 real RHS_Vmtilde_term = INatilde + IbNatilde + IK1tilde + Itotilde + IKrtilde + IKstilde + ICaLtilde + INaKtilde + INaCatilde + IpCatilde + IpKtilde + IbCatilde;
//                 partRHS[i] = delta_t * (stim - RHS_Vmtilde_term);

//                 // Update state variables
//                 // RHS of the state variables with tilde approximations
//                 // Rush-Larsen method - auxiliary variables
//                 real X_r1_inftilde = 1.0f / (1.0f + exp((-26.0f - Vmtilde) / 7.0f));
//                 real alpha_X_r1tilde = 450.0f / (1.0f + exp((-45.0f - Vmtilde) / 10.0f));
//                 real beta_X_r1tilde = 6.0f / (1.0f + exp((30.0f + Vmtilde) / 11.5f));

//                 real X_r2_inftilde = 1.0f / (1.0f + exp((Vmtilde + 88.0f) / 24.0f));
//                 real alpha_X_r2tilde = 3.0f / (1.0f + exp((-60.0f - Vmtilde) / 20.0f));
//                 real beta_X_r2tilde = 1.12f / (1.0f + exp((Vmtilde - 60.0f) / 20.0f));

//                 real X_s_inftilde = 1.0f / (1.0f + exp((-5.0f - Vmtilde) / 14.0f));
//                 real alpha_X_stilde = 1400.0f / sqrt(1.0f + exp((5.0f - Vmtilde) / 6.0f));
//                 real beta_X_stilde = 1.0f / (1.0f + exp((-35.0f + Vmtilde) / 15.0f));

//                 real m_inftilde = 1.0f / ((1.0f + exp((-56.86 - Vmtilde) / 9.03f)) * (1.0f + exp((-56.86f - Vmtilde) / 9.03f)));
//                 real alpha_mtilde = 1.0f / (1.0f + exp((-60.0f - Vmtilde) / 5.0f));
//                 real beta_mtilde = 0.1f / (1.0f + exp((Vmtilde + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((Vmtilde - 50.0f) / 200.0f)));

//                 real h_inftilde = 1.0f / ((1.0f + exp((Vmtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vmtilde + 71.55f) / 7.43f)));
//                 real alpha_htilde;
//                 (Vmtilde < -40.0f)
//                     ? (alpha_htilde = 0.057f * exp(-(Vmtilde + 80.0f) / 6.8f))
//                     : (alpha_htilde = 0.0f);
//                 real beta_htilde;
//                 (Vmtilde < -40.0f)
//                     ? (beta_htilde = 2.7f * exp(0.079f * Vmtilde) + 3.1f * 1.0e5f * exp(0.3485f * Vmtilde))
//                     : (beta_htilde = 0.77f / (0.13f * (1.0f + exp((Vmtilde + 10.66f) / -11.1f))));

//                 real j_inftilde = 1.0f / ((1.0f + exp((Vmtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vmtilde + 71.55f) / 7.43f)));
//                 real alpha_jtilde;
//                 (Vmtilde < -40.0f)
//                     ? (alpha_jtilde = ((-25428.0f * exp(0.2444f * Vmtilde) - (6.948e-6f * exp((-0.04391f) * Vmtilde))) * (Vmtilde + 37.78f)) / (1.0f + exp(0.311f * (Vmtilde + 79.23f))))
//                     : (alpha_jtilde = 0.0f);
//                 real beta_jtilde;
//                 (Vmtilde < -40.0f)
//                     ? (beta_jtilde = (0.02424f * exp(-0.01052f * Vmtilde)) / (1.0f + exp(-0.1378f * (Vmtilde + 40.14f))))
//                     : (beta_jtilde = (0.6f * exp(0.057f * Vmtilde)) / (1.0f + exp(-0.1f * (Vmtilde + 32.0f))));

//                 real d_inftilde = 1.0f / (1.0f + exp((-8.0f - Vmtilde) / 7.5f));
//                 real alpha_dtilde = 1.4f / (1.0f + exp((-35.0f - Vmtilde) / 13.0f)) + 0.25f;
//                 real beta_dtilde = 1.4f / (1.0f + exp((Vmtilde + 5.0f) / 5.0f));
//                 real gamma_dtilde = 1.0f / (1.0f + exp((50.0f - Vmtilde) / 20.0f));

//                 real f_inftilde = 1.0f / (1.0f + exp((Vmtilde + 20.0f) / 7.0f));
//                 real alpha_ftilde = 1102.5f * exp(-(Vmtilde + 27.0f) * (Vmtilde + 27.0f) / 225.0f);
//                 real beta_ftilde = 200.0f / (1.0f + exp((13.0f - Vmtilde) / 10.0f));
//                 real gamma_ftilde = 180.0f / (1.0f + exp((Vmtilde + 30.0f) / 10.0f)) + 20.0f;

//                 real f2_inftilde = 0.67f / (1.0f + exp((Vmtilde + 35.0f) / 7.0f)) + 0.33f;
//                 real alpha_f2tilde; // !!!
//                 alpha_f2tilde = 562.0f * exp(-(Vmtilde + 27.0f) * (Vmtilde + 27.0f) / 240.0f);
//                 real beta_f2tilde = 31.0f / (1.0f + exp((25.0f - Vmtilde) / 10.0f));
//                 real gamma_f2tilde; // !!!
//                 gamma_f2tilde = 80.0f / (1.0f + exp((30.0f + Vmtilde) / 10.0f));

//                 real fCaSS_inftilde = 0.6f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 0.4f;
//                 real tau_fCaSStilde = 80.0f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 2.0f;

// #if defined(EPI) || defined(MCELL)

//                 real s_inftilde = 1.0f / (1.0f + exp((Vmtilde + 20.0f) / 5.0f));
//                 real tau_stilde = 85.0f * exp(-(Vmtilde + 45.0f) * (Vmtilde + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((Vmtilde - 20.0f) / 5.0f)) + 3.0f;

// #endif // EPI || MCELL
// #ifdef ENDO

//                 real s_inftilde = 1.0f / (1.0f + exp((Vmtilde + 28.0f) / 5.0f));
//                 real tau_stilde = 1000.0f * exp(-(Vmtilde + 67.0f) * (Vmtilde + 67.0f) / 1000.0f) + 8.0f;

// #endif // ENDO

//                 real r_inftilde = 1.0f / (1.0f + exp((20.0f - Vmtilde) / 6.0f));
//                 real tau_rtilde = 9.5f * exp(-(Vmtilde + 40.0f) * (Vmtilde + 40.0f) / 1800.0f) + 0.8f;

//                 // Explicit method - auxiliary variables
//                 real Ileaktilde = V_leak * (Ca_SRtilde - Ca_itilde);
//                 real Iuptilde = V_maxup / (1.0f + (K_up * K_up) / (Ca_itilde * Ca_itilde));
//                 real k_CaSRtilde = max_SR - (max_SR - min_SR) / (1.0f + (EC / Ca_SRtilde) * (EC / Ca_SRtilde));
//                 real k1tilde = k1_prime / k_CaSRtilde;
//                 real k2tilde = k2_prime * k_CaSRtilde;
//                 real Otilde = k1tilde * Ca_SStilde * Ca_SStilde * R_primetilde / (k3 + k1tilde * Ca_SStilde * Ca_SStilde);
//                 real Ireltilde = V_rel * Otilde * (Ca_SRtilde - Ca_SStilde);
//                 real Ixfertilde = V_xfer * (Ca_SStilde - Ca_itilde);
//                 real Ca_i_bufCtilde; // !!!
//                 Ca_i_bufCtilde = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + Ca_itilde) * (K_bufC + Ca_itilde)));
//                 real Ca_SR_bufSRtilde; // !!!
//                 Ca_SR_bufSRtilde = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + Ca_SRtilde) * (K_bufSR + Ca_SRtilde)));
//                 real Ca_SS_bufSStilde; // !!!
//                 Ca_SS_bufSStilde = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + Ca_SStilde) * (K_bufSS + Ca_SStilde)));

//                 // Explicit method - RHS of the state variables
//                 real RHS_R_primetilde_term = (k4 * (1.0f - R_primetilde)) - k2tilde * Ca_SStilde * R_primetilde;
//                 real RHS_Ca_itilde_term = Ca_i_bufCtilde * ((((Ileaktilde - Iuptilde) * V_SR / V_C) + Ixfertilde) - (((IbCatilde + IpCatilde) - 2.0f * INaCatilde) * Cm / (2.0f * V_C * F)));
//                 real RHS_Ca_SRtilde_term = Ca_SR_bufSRtilde * (Iuptilde - Ileaktilde - Ireltilde);
//                 real RHS_Ca_SStilde_term = Ca_SS_bufSStilde * (((-ICaLtilde * Cm / (2.0 * V_SS * F)) + (Ireltilde * V_SR / V_SS)) - (Ixfertilde * V_C / V_SS));
//                 real RHS_Na_itilde_term = -((INatilde + IbNatilde + 3.0f * INaKtilde + 3.0f * INaCatilde) * Cm / (V_C * F));
//                 real RHS_K_itilde_term = -((IK1tilde + Itotilde + IKrtilde + IKstilde - 2.0f * INaKtilde + IpKtilde + stim) * Cm / (V_C * F));

//                 // Rush-Larsen method - update approximations
//                 X_r1[i] = X_r1_inftilde - (X_r1_inftilde - actualX_r1) * exp(-delta_t / (alpha_X_r1tilde * beta_X_r1tilde));
//                 X_r2[i] = X_r2_inftilde - (X_r2_inftilde - actualX_r2) * exp(-delta_t / (alpha_X_r2tilde * beta_X_r2tilde));
//                 X_s[i] = X_s_inftilde - (X_s_inftilde - actualX_s) * exp(-delta_t / (alpha_X_stilde * beta_X_stilde + 80.0f));
//                 m[i] = m_inftilde - (m_inftilde - actualm) * exp(-delta_t / (alpha_mtilde * beta_mtilde));
//                 h[i] = h_inftilde - (h_inftilde - actualh) * exp(-delta_t * (alpha_htilde + beta_htilde));
//                 j[i] = j_inftilde - (j_inftilde - actualj) * exp(-delta_t * (alpha_jtilde + beta_jtilde));
//                 d[i] = d_inftilde - (d_inftilde - actuald) * exp(-delta_t / (alpha_dtilde * beta_dtilde + gamma_dtilde));
//                 f[i] = f_inftilde - (f_inftilde - actualf) * exp(-delta_t / (alpha_ftilde + beta_ftilde + gamma_ftilde));
//                 f2[i] = f2_inftilde - (f2_inftilde - actualf2) * exp(-delta_t / (alpha_f2tilde + beta_f2tilde + gamma_f2tilde));
//                 fCaSS[i] = fCaSS_inftilde - (fCaSS_inftilde - actualfCaSS) * exp(-delta_t / tau_fCaSStilde);
//                 s[i] = s_inftilde - (s_inftilde - actuals) * exp(-delta_t / tau_stilde);
//                 r[i] = r_inftilde - (r_inftilde - actualr) * exp(-delta_t / tau_rtilde);

//                 // Explicit method - update approximations
//                 R_prime[i] = actualR_prime + (delta_t * RHS_R_primetilde_term);
//                 Ca_i[i] = actualCa_i + (delta_t * RHS_Ca_itilde_term);
//                 Ca_SR[i] = actualCa_SR + (delta_t * RHS_Ca_SRtilde_term);
//                 Ca_SS[i] = actualCa_SS + (delta_t * RHS_Ca_SStilde_term);
//                 Na_i[i] = actualNa_i + (delta_t * RHS_Na_itilde_term);
//                 K_i[i] = actualK_i + (delta_t * RHS_K_itilde_term);

// #endif // TT2
//             }

//             // =================================================
//             //  Calculate Vm at n + 1 -> Result goes to Vm   
//             // =================================================
//             // Calculate the RHS
//             for (int i = 0; i < Nx; i++)
//             {
//                 real actualVm = Vm[i];
//                 diff_term = diff_coeff * phi_x * (Vm[lim(i - 1, Nx)] - 2.0f * actualVm + Vm[lim(i + 1, Nx)]);

//                 LS_b_x[i] = actualVm + (1.0f - THETA) * diff_term + partRHS[i];
//             }

//             // Solve the tridiagonal system
//             tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

//             // Update Vm
//             for (int i = 0; i < Nx; i++)
//             {
//                 Vm[i] = result_x[i];
//             }

// #ifdef SAVE_FRAMES

//             startSaveFramesTime = omp_get_wtime();

//             // If save frames is true and time step is multiple of frame save rate
//             if (timeStepCounter % frameSaveRate == 0)
//             {
//                 // Save frame
//                 fprintf(fpFrames, "%lf\n", actualTime);
//                 saveFrame(fpFrames, Vm, Nx);
//                 SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
//             }

//             endSaveFramesTime = omp_get_wtime();
//             saveFramesElapsedTime += endSaveFramesTime - startSaveFramesTime;

// #endif // SAVE_FRAMES

// #ifdef MEASURE_VELOCITY

//             startMeasureVelocityTime = omp_get_wtime();

//             // Calculate stim velocity
//             if (!stim_velocity_measured)
//             {
//                 if (!aux_stim_velocity_flag)
//                 {
//                     if (Vm[idx0] > 10.0f)
//                     {
//                         t0 = actualTime;
//                         aux_stim_velocity_flag = true;
//                     }
//                 }
//                 else
//                 {
//                     if (Vm[idx1] > 10.0f)
//                     {
//                         t1 = actualTime;
//                         stim_velocity = (x1 - x0) / (t1 - t0); // cm/ms
//                         stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
//                         stim_velocity_measured = true;
//                         INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, stim_velocity);
//                     }
//                 }
//             }

//             finishMeasureVelocityTime = omp_get_wtime();
//             elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;

// #endif // MEASURE_VELOCITY

//             // Update time step counter
//             timeStepCounter++;
//         }
//     }

// #endif // CABLEEQ

//     finishExecutionTime = omp_get_wtime();
//     elapsedExecutionTime += finishExecutionTime - startExecutionTime;

// #ifdef SAVE_FRAMES

//     fprintf(fpFrames, "%lf\n", actualTime);

// #ifndef CABLEEQ

//     saveFrame(fpFrames, Vm, Nx, Ny);

// #else // if def CABLEEQ

//     saveFrame(fpFrames, Vm, Nx);

// #endif // CABLEEQ

//     SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
//     fclose(fpFrames);

// #endif // SAVE_FRAMES

//     printf("Simulation done!\n");
//     printf("\n");

// #if defined(CONVERGENCE_ANALYSIS_FORCING_TERM) && !defined(CABLEEQ)

//     real *exact = (real *)malloc(Nx * Ny * sizeof(real));
//     for (int i = 0; i < Ny; i++)
//         exact[i] = (real *)malloc(Nx * sizeof(real));

//     real norm2error = calculateNorm2Error(Vm, exact, Nx, Ny, totalTime, delta_x, delta_y);

//     char exactFilePath[MAX_STRING_SIZE];
//     snprintf(exactFilePath, MAX_STRING_SIZE, "%s/exact.txt", pathToSaveData);
//     FILE *fpExact = fopen(exactFilePath, "w");

//     char errorsFilePath[MAX_STRING_SIZE];
//     snprintf(errorsFilePath, MAX_STRING_SIZE, "%s/errors.txt", pathToSaveData);
//     FILE *fpErrors = fopen(errorsFilePath, "w");

//     for (int i = 0; i < Ny; i++)
//     {
//         for (int j = 0; j < Nx; j++)
//         {
//             index = i * Nx + j;
//             fprintf(fpExact, "%e ", exact[index]);
//             fprintf(fpErrors, "%e ", abs(Vm[index] - exact[index]));
//         }
//         fprintf(fpExact, "\n");
//         fprintf(fpErrors, "\n");
//     }

//     for (int i = 0; i < Ny; i++)
//         free(exact[i]);
//     free(exact);

//     SUCCESSMSG("Exact solution saved to %s\n", exactFilePath);
//     SUCCESSMSG("Errors saved to %s\n", errorsFilePath);
//     fclose(fpExact);
//     fclose(fpErrors);

// #endif // CONVERGENCE_ANALYSIS_FORCING_TERM && !CABLEEQ

//     // Write infos to file
//     srand(time(NULL));
//     int random_number = (rand() % 500) + 1;
//     char infosFilePath[MAX_STRING_SIZE];
//     snprintf(infosFilePath, MAX_STRING_SIZE, "%s/infos%d.txt", pathToSaveData, random_number);
//     FILE *fpInfos = fopen(infosFilePath, "w");
//     fprintf(fpInfos, "EXECUTION TYPE = %s\n", exec_mode);
//     fprintf(fpInfos, "PRECISION = %s\n", REAL_TYPE);
//     fprintf(fpInfos, "PROBLEM = %s\n", PROBLEM);
//     fprintf(fpInfos, "CELL MODEL = %s\n", cell_model);
//     fprintf(fpInfos, "METHOD = %s\n", method);

// #ifdef THETA

//     fprintf(fpInfos, "theta = %.2f\n", THETA);

// #endif // THETA

//     fprintf(fpInfos, "\n");

// #ifdef CABLEEQ

//     fprintf(fpInfos, "CABLE LENGTH = %.4g cm\n", Lx);

// #else // if not def CABLEEQ

//     fprintf(fpInfos, "DOMAIN LENGTH IN X = %.4g cm\n", Lx);
//     fprintf(fpInfos, "DOMAIN LENGTH IN Y = %.4g cm\n", Ly);

// #endif // CABLEEQ

//     fprintf(fpInfos, "TOTAL TIME = %.4g ms\n", totalTime);
//     fprintf(fpInfos, "\n");

//     fprintf(fpInfos, "delta_t = %.5g ms (%d time steps)\n", delta_t, M);

// #ifdef CABLEEQ

//     fprintf(fpInfos, "delta_x = %.5g cm (%d space steps)\n", delta_x, Nx);

// #else // if not def CABLEEQ

//     fprintf(fpInfos, "delta_x = %.5g cm (%d um) (%d space steps in x)\n", delta_x, CM_TO_UM(delta_x), Nx);
//     fprintf(fpInfos, "delta_y = %.5g cm (%d um) (%d space steps in y)\n", delta_y, CM_TO_UM(delta_y), Ny);
//     fprintf(fpInfos, "TOTAL POINTS IN DOMAIN = %d\n", Nx * Ny);

// #endif // CABLEEQ

// #if defined(MONODOMAIN) || defined(CABLEEQ)

//     fprintf(fpInfos, "DIFFUSION COEFFICIENT = %.8g\n", diff_coeff);
//     fprintf(fpInfos, "NUMBER OF STIMULI = %d\n", numberOfStimuli);
//     for (int i = 0; i < numberOfStimuli; i++)
//     {
//         fprintf(fpInfos, "STIMULUS %d: START TIME = %.5g ms\n", i + 1, stimuli[i].start_time);
//     }

// #endif // MONODOMAIN || CABLEEQ

// #ifdef MEASURE_VELOCITY

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "STIMULUS S1 VELOCITY = %.4g cm/s\n", stim_velocity);

// #endif // MEASURE_VELOCITY

// #ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "NORM-2 ERROR = %lf\n", norm2error);

// #endif // CONVERGENCE_ANALYSIS_FORCING_TERM

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "TIME OF THE FIRST PART = %.5g s\n", elapsedTime1stPart);
//     fprintf(fpInfos, "TIME OF THE SECOND PART = %.5g s\n", elapsedTime2ndPart);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     fprintf(fpInfos, "TIME TO SOLVE THE 1st LINEAR SYSTEM = %.5g s\n", elapsedTime1stLS);
//     fprintf(fpInfos, "TIME TO SOLVE THE 2nd LINEAR SYSTEM = %.5g s\n", elapsedTime2ndLS);

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

// #ifdef MEASURE_VELOCITY

//     fprintf(fpInfos, "TIME TO MEASURE VELOCITY = %.5g s\n", elapsedMeasureVelocityTime);

// #endif // MEASURE_VELOCITY

// #ifdef SAVE_FRAMES

//     fprintf(fpInfos, "TIME TO SAVE FRAMES = %.5g s\n", elapsedSaveFramesTime);

// #endif // SAVE_FRAMES

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "SIMULATION TOTAL EXECUTION TIME = %.5g s\n", elapsedExecutionTime);

//     fprintf(fpInfos, "\n");
//     fprintf(fpInfos, "OUTPUT DIRECTORY = %s\n", pathToSaveData);
//     fclose(fpInfos);

//     INFOMSG("Simulation total execution time = %.5g s\n", elapsedExecutionTime);
//     SUCCESSMSG("Simulation infos saved to %s\n", infosFilePath);

// #ifdef SAVE_LAST_FRAME

//     // Save last frame
//     char lastFrameFilePath[MAX_STRING_SIZE];
//     snprintf(lastFrameFilePath, MAX_STRING_SIZE, "%s/lastframe.txt", pathToSaveData);
//     FILE *fpLast = fopen(lastFrameFilePath, "w");

// #ifndef CABLEEQ

//     saveFrame(fpLast, Vm, Nx, Ny);

// #else // if def CABLEEQ

//     saveFrame(fpLast, Vm, Nx);

// #endif // CABLEEQ

//     SUCCESSMSG("Last frame saved to %s\n", lastFrameFilePath);
//     fclose(fpLast);

// #endif // SAVE_LAST_FRAME

// #ifdef SAVE_LAST_STATE

// #ifdef AFHN

//     char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathW[MAX_STRING_SIZE];
//     snprintf(lastFrameFilePathVm, MAX_STRING_SIZE, "%s/lastframeVm.txt", pathToSaveData);
//     snprintf(lastFrameFilePathW, MAX_STRING_SIZE, "%s/lastframeW.txt", pathToSaveData);

//     FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
//     FILE *fpLastW = fopen(lastFrameFilePathW, "w");

// #ifndef CABLEEQ

//     for (int i = 0; i < Ny; i++)
//     {
//         for (int j = 0; j < Nx; j++)
//         {
//             index = i * Nx + j;
//             fprintf(fpLastVm, "%e ", Vm[index]);
//             fprintf(fpLastW, "%e ", W[index]);
//         }
//         fprintf(fpLastVm, "\n");
//         fprintf(fpLastW, "\n");
//     }

// #else // if def CABLEEQ

//     for (int i = 0; i < Nx; i++)
//     {
//         fprintf(fpLastVm, "%e ", Vm[i]);
//         fprintf(fpLastW, "%e ", W[i]);
//     }

// #endif // CABLEEQ

//     SUCCESSMSG("Last Vm frame saved to %s\n", lastFrameFilePathVm);
//     SUCCESSMSG("Last W frame saved to %s\n", lastFrameFilePathW);
//     fclose(fpLastVm);
//     fclose(fpLastW);

// #endif // AFHN

// #ifdef TT2

//     char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathX_r1[MAX_STRING_SIZE], lastFrameFilePathX_r2[MAX_STRING_SIZE], lastFrameFilePathX_s[MAX_STRING_SIZE],
//         lastFrameFilePathm[MAX_STRING_SIZE], lastFrameFilePathh[MAX_STRING_SIZE], lastFrameFilePathj[MAX_STRING_SIZE], lastFrameFilePathd[MAX_STRING_SIZE], lastFrameFilePathf[MAX_STRING_SIZE],
//         lastFrameFilePathf2[MAX_STRING_SIZE], lastFrameFilePathfCaSS[MAX_STRING_SIZE], lastFrameFilePaths[MAX_STRING_SIZE], lastFrameFilePathr[MAX_STRING_SIZE],
//         lastFrameFilePathR_prime[MAX_STRING_SIZE], lastFrameFilePathCa_i[MAX_STRING_SIZE], lastFrameFilePathCa_SR[MAX_STRING_SIZE], lastFrameFilePathCa_SS[MAX_STRING_SIZE],
//         lastFrameFilePathNa_i[MAX_STRING_SIZE], lastFrameFilePathK_i[MAX_STRING_SIZE];
//     snprintf(lastFrameFilePathVm, MAX_STRING_SIZE, "%s/lastframeVm.txt", pathToSaveData);
//     snprintf(lastFrameFilePathX_r1, MAX_STRING_SIZE, "%s/lastframeX_r1.txt", pathToSaveData);
//     snprintf(lastFrameFilePathX_r2, MAX_STRING_SIZE, "%s/lastframeX_r2.txt", pathToSaveData);
//     snprintf(lastFrameFilePathX_s, MAX_STRING_SIZE, "%s/lastframeX_s.txt", pathToSaveData);
//     snprintf(lastFrameFilePathm, MAX_STRING_SIZE, "%s/lastframem.txt", pathToSaveData);
//     snprintf(lastFrameFilePathh, MAX_STRING_SIZE, "%s/lastframeh.txt", pathToSaveData);
//     snprintf(lastFrameFilePathj, MAX_STRING_SIZE, "%s/lastframej.txt", pathToSaveData);
//     snprintf(lastFrameFilePathd, MAX_STRING_SIZE, "%s/lastframed.txt", pathToSaveData);
//     snprintf(lastFrameFilePathf, MAX_STRING_SIZE, "%s/lastframef.txt", pathToSaveData);
//     snprintf(lastFrameFilePathf2, MAX_STRING_SIZE, "%s/lastframef2.txt", pathToSaveData);
//     snprintf(lastFrameFilePathfCaSS, MAX_STRING_SIZE, "%s/lastframefCaSS.txt", pathToSaveData);
//     snprintf(lastFrameFilePaths, MAX_STRING_SIZE, "%s/lastframes.txt", pathToSaveData);
//     snprintf(lastFrameFilePathr, MAX_STRING_SIZE, "%s/lastframer.txt", pathToSaveData);
//     snprintf(lastFrameFilePathR_prime, MAX_STRING_SIZE, "%s/lastframeR_prime.txt", pathToSaveData);
//     snprintf(lastFrameFilePathCa_i, MAX_STRING_SIZE, "%s/lastframeCa_i.txt", pathToSaveData);
//     snprintf(lastFrameFilePathCa_SR, MAX_STRING_SIZE, "%s/lastframeCa_SR.txt", pathToSaveData);
//     snprintf(lastFrameFilePathCa_SS, MAX_STRING_SIZE, "%s/lastframeCa_SS.txt", pathToSaveData);
//     snprintf(lastFrameFilePathNa_i, MAX_STRING_SIZE, "%s/lastframeNa_i.txt", pathToSaveData);
//     snprintf(lastFrameFilePathK_i, MAX_STRING_SIZE, "%s/lastframeK_i.txt", pathToSaveData);

//     FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
//     FILE *fpLastX_r1 = fopen(lastFrameFilePathX_r1, "w");
//     FILE *fpLastX_r2 = fopen(lastFrameFilePathX_r2, "w");
//     FILE *fpLastX_s = fopen(lastFrameFilePathX_s, "w");
//     FILE *fpLastm = fopen(lastFrameFilePathm, "w");
//     FILE *fpLasth = fopen(lastFrameFilePathh, "w");
//     FILE *fpLastj = fopen(lastFrameFilePathj, "w");
//     FILE *fpLastd = fopen(lastFrameFilePathd, "w");
//     FILE *fpLastf = fopen(lastFrameFilePathf, "w");
//     FILE *fpLastf2 = fopen(lastFrameFilePathf2, "w");
//     FILE *fpLastfCaSS = fopen(lastFrameFilePathfCaSS, "w");
//     FILE *fpLasts = fopen(lastFrameFilePaths, "w");
//     FILE *fpLastr = fopen(lastFrameFilePathr, "w");
//     FILE *fpLastR_prime = fopen(lastFrameFilePathR_prime, "w");
//     FILE *fpLastCa_i = fopen(lastFrameFilePathCa_i, "w");
//     FILE *fpLastCa_SR = fopen(lastFrameFilePathCa_SR, "w");
//     FILE *fpLastCa_SS = fopen(lastFrameFilePathCa_SS, "w");
//     FILE *fpLastNa_i = fopen(lastFrameFilePathNa_i, "w");
//     FILE *fpLastK_i = fopen(lastFrameFilePathK_i, "w");

// #ifndef CABLEEQ

//     for (int i = 0; i < Ny; i++)
//     {
//         for (int j = 0; j < Nx; j++)
//         {
//             index = i * Nx + j;
//             fprintf(fpLastVm, "%e ", Vm[index]);
//             fprintf(fpLastX_r1, "%e ", X_r1[index]);
//             fprintf(fpLastX_r2, "%e ", X_r2[index]);
//             fprintf(fpLastX_s, "%e ", X_s[index]);
//             fprintf(fpLastm, "%e ", m[index]);
//             fprintf(fpLasth, "%e ", h[index]);
//             fprintf(fpLastj, "%e ", _j[index]);
//             fprintf(fpLastd, "%e ", d[index]);
//             fprintf(fpLastf, "%e ", f[index]);
//             fprintf(fpLastf2, "%e ", f2[index]);
//             fprintf(fpLastfCaSS, "%e ", fCaSS[index]);
//             fprintf(fpLasts, "%e ", s[index]);
//             fprintf(fpLastr, "%e ", r[index]);
//             fprintf(fpLastR_prime, "%e ", R_prime[index]);
//             fprintf(fpLastCa_i, "%e ", Ca_i[index]);
//             fprintf(fpLastCa_SR, "%e ", Ca_SR[index]);
//             fprintf(fpLastCa_SS, "%e ", Ca_SS[index]);
//             fprintf(fpLastNa_i, "%e ", Na_i[index]);
//             fprintf(fpLastK_i, "%e ", K_i[index]);
//         }
//         fprintf(fpLastVm, "\n");
//         fprintf(fpLastX_r1, "\n");
//         fprintf(fpLastX_r2, "\n");
//         fprintf(fpLastX_s, "\n");
//         fprintf(fpLastm, "\n");
//         fprintf(fpLasth, "\n");
//         fprintf(fpLastj, "\n");
//         fprintf(fpLastd, "\n");
//         fprintf(fpLastf, "\n");
//         fprintf(fpLastf2, "\n");
//         fprintf(fpLastfCaSS, "\n");
//         fprintf(fpLasts, "\n");
//         fprintf(fpLastr, "\n");
//         fprintf(fpLastR_prime, "\n");
//         fprintf(fpLastCa_i, "\n");
//         fprintf(fpLastCa_SR, "\n");
//         fprintf(fpLastCa_SS, "\n");
//         fprintf(fpLastNa_i, "\n");
//         fprintf(fpLastK_i, "\n");
//     }

// #else // if def CABLEEQ

//     for (int i = 0; i < Nx; i++)
//     {
//         fprintf(fpLastVm, "%e ", Vm[i]);
//         fprintf(fpLastX_r1, "%e ", X_r1[i]);
//         fprintf(fpLastX_r2, "%e ", X_r2[i]);
//         fprintf(fpLastX_s, "%e ", X_s[i]);
//         fprintf(fpLastm, "%e ", m[i]);
//         fprintf(fpLasth, "%e ", h[i]);
//         fprintf(fpLastj, "%e ", j[i]);
//         fprintf(fpLastd, "%e ", d[i]);
//         fprintf(fpLastf, "%e ", f[i]);
//         fprintf(fpLastf2, "%e ", f2[i]);
//         fprintf(fpLastfCaSS, "%e ", fCaSS[i]);
//         fprintf(fpLasts, "%e ", s[i]);
//         fprintf(fpLastr, "%e ", r[i]);
//         fprintf(fpLastR_prime, "%e ", R_prime[i]);
//         fprintf(fpLastCa_i, "%e ", Ca_i[i]);
//         fprintf(fpLastCa_SR, "%e ", Ca_SR[i]);
//         fprintf(fpLastCa_SS, "%e ", Ca_SS[i]);
//         fprintf(fpLastNa_i, "%e ", Na_i[i]);
//         fprintf(fpLastK_i, "%e ", K_i[i]);
//     }

// #endif // CABLEEQ

//     SUCCESSMSG("Last Vm frame saved to %s\n", lastFrameFilePathVm);
//     SUCCESSMSG("Last X_r1 frame saved to %s\n", lastFrameFilePathX_r1);
//     SUCCESSMSG("Last X_r2 frame saved to %s\n", lastFrameFilePathX_r2);
//     SUCCESSMSG("Last X_s frame saved to %s\n", lastFrameFilePathX_s);
//     SUCCESSMSG("Last m frame saved to %s\n", lastFrameFilePathm);
//     SUCCESSMSG("Last h frame saved to %s\n", lastFrameFilePathh);
//     SUCCESSMSG("Last j frame saved to %s\n", lastFrameFilePathj);
//     SUCCESSMSG("Last d frame saved to %s\n", lastFrameFilePathd);
//     SUCCESSMSG("Last f frame saved to %s\n", lastFrameFilePathf);
//     SUCCESSMSG("Last f2 frame saved to %s\n", lastFrameFilePathf2);
//     SUCCESSMSG("Last fCaSS frame saved to %s\n", lastFrameFilePathfCaSS);
//     SUCCESSMSG("Last s frame saved to %s\n", lastFrameFilePaths);
//     SUCCESSMSG("Last r frame saved to %s\n", lastFrameFilePathr);
//     SUCCESSMSG("Last R_prime frame saved to %s\n", lastFrameFilePathR_prime);
//     SUCCESSMSG("Last Ca_i frame saved to %s\n", lastFrameFilePathCa_i);
//     SUCCESSMSG("Last Ca_SR frame saved to %s\n", lastFrameFilePathCa_SR);
//     SUCCESSMSG("Last Ca_SS frame saved to %s\n", lastFrameFilePathCa_SS);
//     SUCCESSMSG("Last Na_i frame saved to %s\n", lastFrameFilePathNa_i);
//     SUCCESSMSG("Last K_i frame saved to %s\n", lastFrameFilePathK_i);

//     fclose(fpLastVm);
//     fclose(fpLastX_r1);
//     fclose(fpLastX_r2);
//     fclose(fpLastX_s);
//     fclose(fpLastm);
//     fclose(fpLasth);
//     fclose(fpLastj);
//     fclose(fpLastd);
//     fclose(fpLastf);
//     fclose(fpLastf2);
//     fclose(fpLastfCaSS);
//     fclose(fpLasts);
//     fclose(fpLastr);
//     fclose(fpLastR_prime);
//     fclose(fpLastCa_i);
//     fclose(fpLastCa_SR);
//     fclose(fpLastCa_SS);
//     fclose(fpLastNa_i);
//     fclose(fpLastK_i);

// #endif // TT2

// #ifdef MV

//     // Save last state
//     char lastStateFilePathVm[MAX_STRING_SIZE], lastStateFilePathv[MAX_STRING_SIZE],
//         lastStateFilePathw[MAX_STRING_SIZE], lastStateFilePaths[MAX_STRING_SIZE];
//     snprintf(lastStateFilePathVm, MAX_STRING_SIZE, "%s/laststateVm.txt", pathToSaveData);
//     snprintf(lastStateFilePathv, MAX_STRING_SIZE, "%s/laststatev.txt", pathToSaveData);
//     snprintf(lastStateFilePathw, MAX_STRING_SIZE, "%s/laststatew.txt", pathToSaveData);
//     snprintf(lastStateFilePaths, MAX_STRING_SIZE, "%s/laststates.txt", pathToSaveData);

//     FILE *fpLastStateVm = fopen(lastStateFilePathVm, "w");
//     FILE *fpLastStatev = fopen(lastStateFilePathv, "w");
//     FILE *fpLastStatew = fopen(lastStateFilePathw, "w");
//     FILE *fpLastStates = fopen(lastStateFilePaths, "w");

// #ifndef CABLEEQ

//     for (int i = 0; i < Ny; i++)
//     {
//         for (int j = 0; j < Nx; j++)
//         {
//             index = i * Nx + j;
//             fprintf(fpLastStateVm, "%e ", Vm[index]);
//             fprintf(fpLastStatev, "%e ", v[index]);
//             fprintf(fpLastStatew, "%e ", w[index]);
//             fprintf(fpLastStates, "%e ", s[index]);
//         }
//         fprintf(fpLastStateVm, "\n");
//         fprintf(fpLastStatev, "\n");
//         fprintf(fpLastStatew, "\n");
//         fprintf(fpLastStates, "\n");
//     }

// #else // if def CABLEEQ

//     for (int i = 0; i < Nx; i++)
//     {
//         fprintf(fpLastStateVm, "%e ", Vm[i]);
//         fprintf(fpLastStatev, "%e ", v[i]);
//         fprintf(fpLastStatew, "%e ", w[i]);
//         fprintf(fpLastStates, "%e ", s[i]);
//     }

// #endif // CABLEEQ

//     SUCCESSMSG("Last Vm state saved to %s\n", lastStateFilePathVm);
//     SUCCESSMSG("Last v state saved to %s\n", lastStateFilePathv);
//     SUCCESSMSG("Last w state saved to %s\n", lastStateFilePathw);
//     SUCCESSMSG("Last s state saved to %s\n", lastStateFilePaths);

//     fclose(fpLastStateVm);
//     fclose(fpLastStatev);
//     fclose(fpLastStatew);
//     fclose(fpLastStates);

// #endif // MV

// #endif // SAVE_LAST_STATE

// #ifdef CABLEEQ

//     // Save Action Potential
//     char APFilePath[MAX_STRING_SIZE];
//     snprintf(APFilePath, MAX_STRING_SIZE, "%s/AP.txt", pathToSaveData);
//     FILE *fpAP = fopen(APFilePath, "w");

//     for (int i = 0; i < M; i++)
//     {
//         fprintf(fpAP, "%e ", AP[i]);
//     }

//     SUCCESSMSG("Action Potential saved to %s\n", APFilePath);
//     fclose(fpAP);

// #endif // CABLEEQ

//     // Free memory
//     free(time_array);
//     free(pathToSaveData);

// #ifndef CABLEEQ

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     free(c_prime_y);
//     free(d_prime_y);
//     free(LS_b_y);
//     free(result_y);

//     free(la_y);
//     free(lb_y);
//     free(lc_y);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2)

//     free(RHS);

// #endif // SSIADI || THETASSIADI || THETASSIRK2

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

// #endif // not CABLEEQ

//     free(Vm);
//     free(partRHS);

// #if defined(SSIADI) || defined(THETASSIADI) || defined(THETASSIRK2) || defined(OSADI)

//     free(c_prime_x);
//     free(d_prime_x);
//     free(LS_b_x);
//     free(result_x);

//     free(la_x);
//     free(lb_x);
//     free(lc_x);

// #endif // SSIADI || THETASSIADI || THETASSIRK2 || OSADI

// #if defined(MONODOMAIN) || defined(CABLEEQ)

// #ifdef CABLEEQ

//     free(AP);

// #endif // CABLEEQ

// #ifdef AFHN

//     free(W);

// #endif // AFHN

// #ifdef TT2

//     free(X_r1);
//     free(X_r2);
//     free(X_s);
//     free(m);
//     free(h);
//     free(j);
//     free(d);
//     free(f);
//     free(f2);
//     free(fCaSS);
//     free(s);
//     free(r);
//     free(R_prime);
//     free(Ca_i);
//     free(Ca_SR);
//     free(Ca_SS);
//     free(Na_i);
//     free(K_i);

// #endif // TT2

// #ifdef MV

//     free(v);
//     free(w);
//     free(s);

// #endif // MV

// #endif // MONODOMAIN || CABLEEQ

//     return;
// }