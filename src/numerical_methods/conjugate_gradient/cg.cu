// ==============================================
//              Conjugate Gradient
// ==============================================

#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"
#include "../../cell_models/cell_models.h"
#include <thrust/inner_product.h>
#include <thrust/device_vector.h>

static __global__ void prepareRHSAndUpdateSV(const int Nx, const int Ny, const real delta_t,
                                             const real actualTime, const int num_active_stimuli, const Stimulus *d_active_stimuli,
                                             const real *d_Vm, real *d_sV, real *d_RHS, const CellModel cell_model, const real denom_chiCm)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Declare auxiliary array
        real d_actualsV[MAX_NSV];

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        const real actualVm = d_Vm[idx];
        select_get_actual_sV(cell_model, d_actualsV, d_sV, idx);

        // Stimulation
        const real stim = get_stimulus_value(actualTime, i, j, d_active_stimuli, num_active_stimuli);

        // Prepare the RHS of the following linear system
        d_RHS[idx] = actualVm + delta_t * (stim - select_compute_dVmdt(cell_model, actualVm, d_actualsV));

        // Update state variables
        select_update_sV(cell_model, d_sV, d_actualsV, actualVm, d_actualsV, delta_t, idx);
    }
}

static __global__ void apply_diff_operator(const real *x, real *Ax, const real *Dxx, const real *Dyy, const real *Dxy,
                                           const int Nx, const int Ny, const real delta_x, const real delta_y, const real delta_t, const real denom_chiCm, const bool is_aligned)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Compute the diffusion term without rotation
        const real lap = select_compute_diffusion_term(is_aligned, x, Dxx, Dyy, Dxy, i, j, Nx, Ny, delta_x, delta_y);

        // Update Ax with the diffusion term
        Ax[idx] = x[idx] - delta_t * lap * denom_chiCm;   
    }
}

static __global__ void compute_residual(const real *b, const real *Ax, real *r, const int Nx, const int Ny)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
        r[idx] = b[idx] - Ax[idx]; 
}

static __global__ void update_x_and_residual(const real *p, const real alphak, real *x, const real *Ap, real *r, const int Nx, const int Ny)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Update xₖ₊₁ = xₖ + αₖ·pₖ
        x[idx] += alphak * p[idx];

        // Update residual rₖ₊₁ = rₖ - αₖ·A·pₖ
        r[idx] -= alphak * Ap[idx];
    }
}

static __global__ void update_direction(const real *r, const real *p, const real betak, real *new_p, const int Nx, const int Ny)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
        new_p[idx] = r[idx] + betak * p[idx];
}

void runCG_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array,
                const CellModelSolver *cell_model_solver, real *Vm, real *sV, const real *Dxx, const real *Dyy, const real *Dxy)
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
    real *d_Dxx, *d_Dyy, *d_Dxy;
    const int total_points = Nx * Ny;

    CUDA_CALL(cudaMalloc(&d_Vm, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_sV, total_points * cell_model_solver->n_state_vars * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMalloc(&d_Dxx, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Dyy, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Dxy, total_points * sizeof(real)));

    CUDA_CALL(cudaMemcpy(d_Vm, Vm, total_points * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_sV, sV, total_points * cell_model_solver->n_state_vars * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Dxx, Dxx, total_points * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Dyy, Dyy, total_points * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Dxy, Dxy, total_points * sizeof(real), cudaMemcpyHostToDevice));

    // Auxiliary variables for the operations
    const int max_iterations = 1000; // Maximum number of iterations for CG
    const real tolerance = 1e-6; // Tolerance for convergence
    real *d_RHS, *d_r, *d_p, *d_Ap;
    CUDA_CALL(cudaMalloc(&d_RHS, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_r, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_p, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ap, total_points * sizeof(real)));

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

    // Variables for time measurement
    real startTime = 0.0f;
    real startExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;

    SIMPLEMSG("");
    INFOMSG("Starting simulation with CG (CUDA)...\n");

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

        // Launch kernel to compute the reaction term and update state variables
        prepareRHSAndUpdateSV<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, actualTime, numberOfStimuli, d_stimuli,
                                                                           d_Vm, d_sV, d_RHS, cell_model, denom_chiCm);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // =================================================
        //  Solve Linear System with CG
        // =================================================
        startTime = omp_get_wtime();

        // Initialize residual r₀ = b - A·x₀
        apply_diff_operator<<<fullDomainGridSize, fullDomainBlockSize>>>(d_Vm, d_Ap, d_Dxx, d_Dyy, d_Dxy, Nx, Ny, delta_x, delta_y, delta_t, denom_chiCm, is_aligned);
        CUDA_CALL(cudaDeviceSynchronize());

        compute_residual<<<fullDomainGridSize, fullDomainBlockSize>>>(d_RHS, d_Ap, d_r, Nx, Ny);
        CUDA_CALL(cudaDeviceSynchronize());

        // Initialize direction p₀ = r₀
        CUDA_CALL(cudaMemcpy(d_p, d_r, total_points * sizeof(real), cudaMemcpyDeviceToDevice));

        // Initialize variables for CG
        real rs_old = thrust::inner_product(
                                            thrust::device_pointer_cast(d_r),
                                            thrust::device_pointer_cast(d_r) + total_points,
                                            thrust::device_pointer_cast(d_r),
                                            0.0);
        
        // Iterate until convergence or maximum iterations
        for (int iter = 0; iter < max_iterations; ++iter)
        {
            // Compute A·pₖ
            apply_diff_operator<<<fullDomainGridSize, fullDomainBlockSize>>>(d_p, d_Ap, d_Dxx, d_Dyy, d_Dxy, Nx, Ny, delta_x, delta_y, delta_t, denom_chiCm, is_aligned);
            CUDA_CALL(cudaDeviceSynchronize());

            // Compute the step size αₖ = rₖᵀ·rₖ / (pₖᵀ·A·pₖ)
            real p_Ap = thrust::inner_product(
                                            thrust::device_pointer_cast(d_p),
                                            thrust::device_pointer_cast(d_p) + total_points,
                                            thrust::device_pointer_cast(d_Ap),
                                            0.0);
            real alphak = rs_old / p_Ap;

            // Update xₖ₊₁ = xₖ + αₖ·pₖ and rₖ₊₁ = rₖ - αₖ·A·pₖ
            update_x_and_residual<<<fullDomainGridSize, fullDomainBlockSize>>>(d_p, alphak, d_Vm, d_Ap, d_r, Nx, Ny);
            CUDA_CALL(cudaDeviceSynchronize());

            real rs_new = thrust::inner_product(
                                            thrust::device_pointer_cast(d_r),
                                            thrust::device_pointer_cast(d_r) + total_points,
                                            thrust::device_pointer_cast(d_r),
                                            0.0);

            // Check for convergence
            if (rs_new < tolerance * tolerance) break;

            // Compute the new direction pₖ₊₁ = rₖ₊₁ + βₖ·pₖ
            real betak = rs_new / rs_old;
            update_direction<<<fullDomainGridSize, fullDomainBlockSize>>>(d_r, d_p, betak, d_p, Nx, Ny);
            CUDA_CALL(cudaDeviceSynchronize());

            // Update the old residual
            rs_old = rs_new;
        }

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
    CUDA_CALL(cudaFree(d_Dxx));
    CUDA_CALL(cudaFree(d_Dyy));
    CUDA_CALL(cudaFree(d_Dxy));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_r));
    CUDA_CALL(cudaFree(d_p));
    CUDA_CALL(cudaFree(d_Ap));
}
