#include "monodomain.h"

int runMonodomainSimulationCUDA(const SimulationConfig *config)
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

    // If the number os state variables is greater than MAX_NSV, print an error message and exit
    if (cell_model_solver->n_state_vars > MAX_NSV)
    {
        ERRORMSG("Number of state variables (%d) exceeds the maximum allowed (%d).", cell_model_solver->n_state_vars, MAX_NSV);
        free(time_array);
        return -2;
    }

    // If Nx or Ny is greater than MAX_SYS_SIZE, print an error message and exit
    if (config->Nx > MAX_SYS_SIZE || config->Ny > MAX_SYS_SIZE)
    {
        ERRORMSG("Nx (%d) or Ny (%d) is greater than the maximum system allowed for Thomas (%d).", config->Nx, config->Ny, MAX_SYS_SIZE);
        free(time_array);
        return -2;
    }
    
    // Allocate and initialize state variables arrays
    int total_points = config->Nx * config->Ny;
    real *Vm = (real *)malloc(total_points * sizeof(real));
    real *sV = (real *)malloc(total_points * cell_model_solver->n_state_vars * sizeof(real));
    if (config->restore_state)
        restore_simulation_state(config->restore_path, Vm, sV, config->Nx, config->Ny, cell_model_solver->n_state_vars);
    else
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
        return -3;
    }
    initializeProperties(config, Dxx, Dyy, Dxy);
    
    // Run the simulation based on the selected method
    numerical_method_t run_method = get_numerical_method(&config->exec_mode, &config->method);
    if (run_method == NULL)
    {
        ERRORMSG("Invalid CUDA numerical method selected.");
        free(time_array);
        free(Vm);
        free(sV);
        free(Dxx);
        free(Dyy);
        free(Dxy);
        return -4;
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
        measurement.elapsedSaveFramesTime += omp_get_wtime() - startTime;
        SUCCESSMSG("Last frame (%.2f ms) saved to %s\n", config->M * config->dt, file_path);
        snprintf(file_path, MAX_STRING_SIZE, "%s/frames/lastframe.%s", config->output_dir, config->file_extension);
        config->save_function(file_path, Vm, config->Nx, config->Ny, config->dx, config->dy);
    }

    // Save last state
    if (config->save_state)
    {
        real startTime = omp_get_wtime();
        static char file_path[MAX_STRING_SIZE];
        snprintf(file_path, MAX_STRING_SIZE, "%s/simulation.dat", config->output_dir);
        save_simulation_state(file_path, Vm, sV, config->Nx, config->Ny, cell_model_solver->n_state_vars);
        measurement.elapsedSaveStateTime += omp_get_wtime() - startTime;
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