// ==============================================
//           Hundsdorfer-Verwer Scheme
// ==============================================

#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"
#include "../../cell_models/cell_models.h"

#ifndef HV_THETA
#define HV_THETA 0.25f
#endif // HV_THETA

static __global__ void computeY0AndUpdateSV(const int Nx, const int Ny, const real delta_t, const real delta_x, const real delta_y,
                                            const real actualTime, const int num_active_stimuli, const Stimulus *d_active_stimuli,
                                            const real *d_Vm, real *d_sV, real *d_Y0, real *d_F, const ElementProperties *d_elements, const CellModel cell_model, const real denom_chiCm)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Declare auxiliary arrays
        real d_actualsV[MAX_NSV];
        real d_sVtilde[MAX_NSV];

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        const real actualVm_center = d_Vm[idx];
        select_get_actual_sV(cell_model, d_actualsV, d_sV, idx);

        // Stimulation
        const real stim = get_stimulus_value(actualTime, i, j, d_active_stimuli, num_active_stimuli);

        // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*dt*(diffusion + R(Vmn, sVn))
        const real diff_term = compute_diffusion_term_anisotropic(d_Vm, d_elements, i, j, Nx, Ny, delta_x, delta_y);
        const real dVmdt = select_compute_dVmdt(cell_model, actualVm_center, d_actualsV);
        const real Vmtilde = actualVm_center + 0.5f * delta_t * ((diff_term * denom_chiCm) + stim - dVmdt);

        // Calculate approximation for state variables using Vmtilde and update them
        select_update_sVtilde(cell_model, d_sVtilde, actualVm_center, d_actualsV, 0.5f * delta_t);
        select_update_sV(cell_model, d_sV, d_actualsV, Vmtilde, d_sVtilde, delta_t, idx);

        // Update d_Y0 with the new approximation -> Y_0 = U_n-1 + dt * F(t_n-1, U_n-1)
        const real computeF = (diff_term * denom_chiCm) + stim - dVmdt;
        d_F[idx] = computeF; // Store the computed F for later use
        d_Y0[idx] = actualVm_center + delta_t * computeF; // Update Y0 with the new approximation
    }
}

static __global__ void computeY0tilde(const int Nx, const int Ny, const real delta_t, const real delta_x, const real delta_y,
                                      const real actualTime, const int num_active_stimuli, const Stimulus *d_active_stimuli,
                                      const real *d_Y, real *d_sV_tn, real *d_Y0, real *d_F, real *d_Y0tilde, const ElementProperties *d_elements, const CellModel cell_model, const real denom_chiCm)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Declare auxiliary arrays
        real d_actualsV_tn[MAX_NSV];

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        const real actualY0_center = d_Y0[idx];
        const real Y_tn_center = d_Y[idx];
        const real actualF_center = d_F[idx];
        select_get_actual_sV(cell_model, d_actualsV_tn, d_sV_tn, idx);

        // Stimulation
        const real stim = get_stimulus_value(actualTime, i, j, d_active_stimuli, num_active_stimuli);

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        const real diff_term = compute_diffusion_term_anisotropic(d_Y, d_elements, i, j, Nx, Ny, delta_x, delta_y);
        const real dVmdt = select_compute_dVmdt(cell_model, Y_tn_center, d_actualsV_tn);
        const real compute_F_tn = (diff_term * denom_chiCm) + stim - dVmdt;

        // Update d_Y0tilde with the new approximation -> Y_0tilde = Y0 + 0.5 * dt * (F(t_n, Y_k) - F(t_n-1, U_n-1))
        d_Y0tilde[idx] = actualY0_center + 0.5f * delta_t * (compute_F_tn - actualF_center);
    }
}

static __global__ void prepareRHS_x(const int Nx, const int Ny, const real delta_x, const real delta_y,
                                    const real delta_t, const real *d_Vm, const real *d_prevY,
                                    real *d_RHS, const ElementProperties *d_elements, const real denom_chiCm)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Calculate the RHS of the linear system with explicit diffusion term along x
        const real prevY_center = d_prevY[idx];
        const real diff_term = compute_diffusion_term_x_axis(d_Vm, d_elements, i, j, Nx, Ny, delta_x, delta_y);

        // Y_1 = Y_0 + theta * dt * F_1(t_n-1, U_n-1) -> 0 represents the first step and 1 represents the x-axis
        d_RHS[idx] = prevY_center - HV_THETA * delta_t * diff_term * denom_chiCm;
    }
}

static __global__ void prepareRHS_y(const int Nx, const int Ny, const real delta_x, const real delta_y,
                                    const real delta_t, const real *d_Vm, const real *d_prevY,
                                    real *d_RHS, const ElementProperties *d_elements, const real denom_chiCm)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Calculate the RHS of the linear system with explicit diffusion term along y
        const real prevY_center = d_prevY[idx];
        const real diff_term = compute_diffusion_term_y_axis(d_Vm, d_elements, i, j, Nx, Ny, delta_x, delta_y);

        // Y_2 = Y_1 + theta * dt * F_2(t_n-1, U_n-1) -> 1 represents the x-axis and 2 represents the y-axis
        d_RHS[idx] = prevY_center - HV_THETA * delta_t * diff_term * denom_chiCm;
    }
}

static __global__ void parallelThomas_x(const int numSys, const int sysSize, real *d_rhs,
                                        const real phi_x, const ElementProperties *d_elements)
{
    // Obtain the index of the thread - each thread will handle a system
    const int sysIdx = blockIdx.x * blockDim.x + threadIdx.x;
    const int offset = sysIdx * sysSize;

    if (sysIdx < numSys)
    {
        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];

        real denom;
        real coeff = phi_x * d_elements[offset].D_xx;
        real lalc = -coeff;
        real lb = 1.0f + coeff; // Coefficient for the first element

        c_prime[0] = lalc / lb;
        d_prime[0] = d_rhs[offset] / lb;

        for (int i = 1; i < sysSize; i++)
        {
            coeff = phi_x * d_elements[offset + i].D_xx;

            lalc = -coeff;
            lb = (i < sysSize - 1) ? 1.0f + 2.0f * coeff : 1.0f + coeff; // Last element has a different coefficient
            denom = 1.0f / (lb - c_prime[i - 1] * lalc);

            c_prime[i] = lalc * denom;
            d_prime[i] = (d_rhs[offset + i] - d_prime[i - 1] * lalc) * denom;
        }

        d_rhs[offset + sysSize - 1] = d_prime[sysSize - 1];

        for (int i = sysSize - 2; i >= 0; i--)
            d_rhs[offset + i] = d_prime[i] - c_prime[i] * d_rhs[offset + i + 1];
    }
}

static __global__ void parallelThomas_y(const int numSys, const int sysSize, real *d_rhs,
                                        const real phi_y, const ElementProperties *d_elements)
{
    // Obtain the index of the thread - each thread will handle a system
    const int sysIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (sysIdx < numSys)
    {
        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];

        real denom;
        real coeff = phi_y * d_elements[sysIdx].D_yy;
        real lalc = -coeff;
        real lb = 1.0f + coeff; // Coefficient for the first element

        c_prime[0] = lalc / lb;
        d_prime[0] = d_rhs[sysIdx] / lb;

        for (int i = 1; i < sysSize; i++)
        {
            coeff = phi_y * d_elements[sysIdx + numSys * i].D_yy;

            lalc = -coeff;
            lb = (i < sysSize - 1) ? 1.0f + 2.0f * coeff : 1.0f + coeff; // Last element has a different coefficient
            denom = 1.0f / (lb - c_prime[i - 1] * lalc);

            c_prime[i] = lalc * denom;
            d_prime[i] = (d_rhs[sysIdx + numSys * i] - d_prime[i - 1] * lalc) * denom;
        }

        d_rhs[sysIdx + numSys * (sysSize - 1)] = d_prime[sysSize - 1];

        for (int i = sysSize - 2; i >= 0; i--)
            d_rhs[sysIdx + numSys * i] = d_prime[i] - c_prime[i] * d_rhs[sysIdx + numSys * (i + 1)];
    }
}

void runHV_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, const ElementProperties *elements)
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

    const CellModel cell_model = config->cell_model;
    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Get the solver functions
    const real denom_chiCm = cell_model_solver->denom_chiCm;
    const real activation_threshold = cell_model_solver->activation_thershold;

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

    // Create device variables, allocate memory on device, and copy data
    real *d_Vm, *d_sV;
    Stimulus *d_stimuli;
    ElementProperties *d_elements;
    const int total_points = Nx * Ny;

    CUDA_CALL(cudaMalloc(&d_Vm, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_sV, total_points * cell_model_solver->n_state_vars * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMalloc(&d_elements, total_points * sizeof(ElementProperties)));

    CUDA_CALL(cudaMemcpy(d_Vm, Vm, total_points * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_sV, sV, total_points * cell_model_solver->n_state_vars * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_elements, elements, total_points * sizeof(ElementProperties), cudaMemcpyHostToDevice));

    // Auxiliary variables for the operations
    real *d_Y, *d_RHS, *d_Y0, *d_F;
    CUDA_CALL(cudaMalloc(&d_Y, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Y0, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_F, total_points * sizeof(real)));

    // Calculate auxiliary coefficients for Thomas algorithm
    const real thomas_coeff_x = (HV_THETA * delta_t * denom_chiCm) / (delta_x * delta_x);
    const real thomas_coeff_y = (HV_THETA * delta_t * denom_chiCm) / (delta_y * delta_y);

    // CUDA grid and block allocation
    // Device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Number of SMs and minimum number of blocks to maximize the parallelism
    const int numSMs = prop.multiProcessorCount;
    const int minBlocks = 2 * numSMs;

    // Print information
    SIMPLEMSG("");
    INFOMSG("Device name: %s (%d SMs)\n", prop.name, numSMs);

    // Calculate the number of blocks and threads for the full domain kernels
    dim3 fullDomainBlockSize(FULL_DOMAIN_BLOCK_SIZE_X, FULL_DOMAIN_BLOCK_SIZE_Y);
    dim3 fullDomainGridSize((Nx + fullDomainBlockSize.x - 1) / fullDomainBlockSize.x, (Ny + fullDomainBlockSize.y - 1) / fullDomainBlockSize.y);

    // Adjust the number of blocks
    if (fullDomainGridSize.x * fullDomainGridSize.y < minBlocks)
        fullDomainGridSize.x = (minBlocks + fullDomainGridSize.y - 1) / fullDomainGridSize.y;

    // Print information
    SIMPLEMSG("");
    INFOMSG("For full domain kernels:\n");
    INFOMSG("Block size: %d x %d threads (total %d threads per block)\n", fullDomainBlockSize.x, fullDomainBlockSize.y, fullDomainBlockSize.x * fullDomainBlockSize.y);
    INFOMSG("Grid size: %d x %d blocks (total %d blocks, total %d threads)\n", fullDomainGridSize.x, fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y * fullDomainBlockSize.x * fullDomainBlockSize.y);

    // Calculate the number of blocks and threads for individual directions of ADI that will be used in Thomas kernel
    const int gridSize_x = (Nx + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;
    const int gridSize_y = (Ny + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;

    // Print information
    SIMPLEMSG("");
    INFOMSG("For Thomas kernel:\n");
    INFOMSG("Grid size for x: %d blocks (%d threads per block, total %d threads)\n", gridSize_x, THOMAS_KERNEL_BLOCK_SIZE, gridSize_x * THOMAS_KERNEL_BLOCK_SIZE);
    INFOMSG("Grid size for y: %d blocks (%d threads per block, total %d threads)\n", gridSize_y, THOMAS_KERNEL_BLOCK_SIZE, gridSize_y * THOMAS_KERNEL_BLOCK_SIZE);

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
    real elapsedTimeApproximation = 0.0f;
    real elapsedTime3rdLS = 0.0f;
    real elapsedTime4thLS = 0.0f;

    SIMPLEMSG("");
    INFOMSG("Starting simulation with HV (CUDA)...\n");

    // Main time loop
    startExecutionTime = omp_get_wtime();

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // =================================================
        //  Calculate Y0 and Update ODEs
        // =================================================
        startTime = omp_get_wtime();

        // Launch kernel to compute the Y0 approximation (including cross derivative terms and reaction) and update state variables
        computeY0AndUpdateSV<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, delta_x, delta_y, actualTime, numberOfStimuli, d_stimuli,
                                                                          d_Vm, d_sV, d_Y0, d_F, d_elements, cell_model, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // =================================================
        //  First step -> Result goes to d_RHS
        //  diffusion implicit and explicit in x
        // =================================================
        startTime = omp_get_wtime();

        prepareRHS_x<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_x, delta_y, delta_t, d_Vm, d_Y0, d_RHS, d_elements, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_x<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_RHS, thomas_coeff_x, d_elements);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stLS += omp_get_wtime() - startLSTime;

        // =================================================
        //  Second step -> Result goes to d_Y
        //  diffusion implicit and explicit in y
        // =================================================
        prepareRHS_y<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_x, delta_y, delta_t, d_Vm, d_RHS, d_Y, d_elements, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_y<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_Y, thomas_coeff_y, d_elements);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime2ndLS += omp_get_wtime() - startLSTime;

        // =================================================
        //  Calculate Y0tilde -> result goes to d_Vm
        // =================================================
        startLSTime = omp_get_wtime();

        computeY0tilde<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, delta_x, delta_y, actualTime, numberOfStimuli, d_stimuli,
                                                                    d_Y, d_sV, d_Y0, d_F, d_Vm, d_elements, cell_model, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTimeApproximation += omp_get_wtime() - startLSTime;

        // =================================================
        //  Third step -> Result goes to d_RHS
        //  diffusion implicit and explicit in x
        // =================================================
        prepareRHS_x<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_x, delta_y, delta_t, d_Y, d_Vm, d_RHS, d_elements, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_x<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_RHS, thomas_coeff_x, d_elements);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime3rdLS += omp_get_wtime() - startLSTime;

        // =================================================
        //  Fourth step -> Result goes to d_Vm
        //  diffusion implicit and explicit in y
        // =================================================
        prepareRHS_y<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_x, delta_y, delta_t, d_Y, d_RHS, d_Vm, d_elements, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_y<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_Vm, thomas_coeff_y, d_elements);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime4thLS += omp_get_wtime() - startLSTime;

        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        if (saveFrames && (timeStepCounter % frameSaveRate == 0))
        {
            startTime = omp_get_wtime();
            CUDA_CALL(cudaMemcpy(Vm, d_Vm, total_points * sizeof(real), cudaMemcpyDeviceToHost));
            handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, Vm, Nx, Ny, delta_x, delta_y, actualTime);
            elapsedSaveFramesTime += omp_get_wtime() - startTime;
        }

        // Measure velocity if needed
        if (measureVelocity && !stim_velocity_measured)
        {
            startTime = omp_get_wtime();

            // Copy only the Vm values needed for velocity measurement
            real Vmidx_x0, Vmidx_x1;
            CUDA_CALL(cudaMemcpy(&Vmidx_x0, &d_Vm[idx_x0], sizeof(real), cudaMemcpyDeviceToHost));
            CUDA_CALL(cudaMemcpy(&Vmidx_x1, &d_Vm[idx_x1], sizeof(real), cudaMemcpyDeviceToHost));
            handle_velocity_measurement(Vmidx_x0, Vmidx_x1, &t0, &t1, activation_threshold, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);

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
    measurement->elapsedTimeApproximation = elapsedTimeApproximation;
    measurement->elapsedTime3rdLS = elapsedTime3rdLS;
    measurement->elapsedTime4thLS = elapsedTime4thLS;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Copy results back to host
    CUDA_CALL(cudaMemcpy(Vm, d_Vm, total_points * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(sV, d_sV, total_points * cell_model_solver->n_state_vars * sizeof(real), cudaMemcpyDeviceToHost));

    // Free allocated memory
    CUDA_CALL(cudaFree(d_Vm));
    CUDA_CALL(cudaFree(d_sV));
    CUDA_CALL(cudaFree(d_stimuli));
    CUDA_CALL(cudaFree(d_elements));
    CUDA_CALL(cudaFree(d_Y));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_Y0));
    CUDA_CALL(cudaFree(d_F));
}