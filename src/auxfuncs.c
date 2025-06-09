#include "../include/auxfuncs.h"

// Populate diagonals for Thomas algorithm
// la -> Subdiagonal
// lb -> Diagonal
// lc -> Superdiagonal
// N -> Number of elements
// coeff -> Coefficient for the diagonals
void populateDiagonalThomasAlgorithm(real *la, real *lb, real *lc, int N, real coeff)
{
    for (int i = 0; i < N; i++)
    {
        la[i] = -coeff;
        lb[i] = 1.0f + 2.0f * coeff;
        lc[i] = -coeff;
    }

    // Set the first and last elements of the diagonals, boundary conditions
    lc[0] = lc[0] + la[0];
    la[N - 1] = la[N - 1] + lc[N - 1];
    la[0] = 0.0f;
    lc[N - 1] = 0.0f;
}

// TODO: Correct this function
void prefactorizeThomas(const real *la, const real *lb, const real *lc, real *c_prime, real *denominator, const int sysSize)
{
    c_prime[0] = lc[0] / lb[0];
    denominator[0] = 1.0f / lb[0];

    for (int i = 1; i < sysSize; i++)
    {
        real denom = 1.0f / (lb[i] - c_prime[i - 1] * la[i]);
        if (i < sysSize - 1)
            c_prime[i] = lc[i] * denom;
        denominator[i] = denom;
    }
}

int createDirectories(char *dir_path, bool remove_old_files)
{
    char path[MAX_STRING_SIZE];
    char frames_path[MAX_STRING_SIZE];

    // Check for path length to avoid buffer overflow
    if (strlen(dir_path) >= MAX_STRING_SIZE)
    {
        ERRORMSG("Path length exceeds maximum size of %d characters.\n", MAX_STRING_SIZE);
        return -1;
    }

    // Copy and normalize path separators to platform-specific style
    size_t len = strlen(dir_path);
    for (size_t i = 0; i < len; i++)
    {
        if (dir_path[i] == '/' || dir_path[i] == '\\')
        {
            path[i] = PATH_SEPARATOR;
        }
        else
        {
            path[i] = dir_path[i];
        }
    }
    path[len] = '\0'; // Ensure null-termination

    // Traverse the path one segment at a time and create directories if needed
    for (char *p = path + 1; *p; p++)
    {
        if (*p == PATH_SEPARATOR)
        {
            char saved = *p;
            *p = '\0'; // Temporarily terminate the string to isolate the subpath

            // Attempt to create the directory
            // If mkdir fail because it already exists (error EEXIST), we ignore the error
            if (mkdir(path, 0777) != 0 && errno != EEXIST)
            {
                perror("mkdir");
                return -1;
            }

            *p = saved; // Restore the original character
        }
    }

    // Create the final (full) directory
    if (mkdir(path, 0777) != 0 && errno != EEXIST)
    {
        perror("mkdir");
        return -1;
    }

    // If remove_old_files is true, remove the old files in the directory, including the frames subdirectory
    if (remove_old_files)
    {
        // Open the directory
        DIR *dir = opendir(path);
        if (dir == NULL)
        {
            perror("opendir");
            return -1;
        }

        // Check if "frames" subdirectory exists and remove its contents
        snprintf(frames_path, sizeof(frames_path), "%s/frames", path);
        DIR *frames_dir = opendir(frames_path);
        if (frames_dir != NULL)
        {
            // Remove files inside "frames" subdirectory
            struct dirent *entry;
            while ((entry = readdir(frames_dir)) != NULL)
            {
                // Skip the current and parent directory entries
                if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
                    continue;

                // Construct the full path to the file
                char full_path[MAX_STRING_SIZE];
                snprintf(full_path, sizeof(full_path), "%s/%s", frames_path, entry->d_name);

                // Remove the file
                if (remove(full_path) != 0)
                {
                    perror("remove");
                    closedir(frames_dir);
                    closedir(dir);
                    return -1;
                }
            }
            closedir(frames_dir);

            // Remove the "frames" directory itself
            if (rmdir(frames_path) != 0 && errno != ENOENT)
            {
                perror("rmdir");
                closedir(dir);
                return -1;
            }
        }

        struct dirent *entry;
        while ((entry = readdir(dir)) != NULL)
        {
            // Skip the current and parent directory entries
            if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
                continue;

            // Construct the full path to the file
            char full_path[MAX_STRING_SIZE];
            snprintf(full_path, sizeof(full_path), "%s/%s", path, entry->d_name);

            // Remove the file
            if (remove(full_path) != 0)
            {
                perror("remove");
                closedir(dir);
                return -1;
            }
        }

        closedir(dir);
    }

    // Check if there is a path separator at the end of the dir_path, if so, remove it
    size_t path_len = strlen(path);
    if (path_len > 0 && (path[path_len - 1] == '/' || path[path_len - 1] == '\\'))
    {
        path[path_len - 1] = '\0';
    }

    // Check if the path is empty after removing the separator
    if (strlen(path) == 0)
    {
        ERRORMSG("Error: The path is empty after removing the separator.\n");
        return -1;
    }

    // Safely copy the path to the original dir_path using snprintf
    snprintf(dir_path, MAX_STRING_SIZE, "%s", path);

    // Create the frames subdirectory
    snprintf(frames_path, sizeof(frames_path), "%s/frames", dir_path);
    if (mkdir(frames_path, 0777) != 0 && errno != EEXIST)
    {
        perror("mkdir");
        return -1;
    }

    return 0;
}

void initializeTimeArray(real *timeArray, int M, real dt)
{
    for (int i = 0; i < M; i++)
        timeArray[i] = i * dt;
}

void initializeMeasurement(Measurement *measurement)
{
    measurement->elapsedExecutionTime = 0.0f;
    measurement->elapsedTime1stPart = 0.0f;
    measurement->elapsedTime2ndPart = 0.0f;
    measurement->elapsedTime1stLS = 0.0f;
    measurement->elapsedTime2ndLS = 0.0f;
    measurement->elapsedSaveFramesTime = 0.0f;
    measurement->elapsedMeasureVelocityTime = 0.0f;
    measurement->elapsedSaveStateTime = 0.0f;
    measurement->stimVelocity = 0.0f;
}

int populateStimuli(SimulationConfig *config)
{
    // Unpack parameters from the config
    real delta_x = config->dx;
    real delta_y = config->dy;
    int numberOfStimuli = config->stimulus_count;

    for (int i = 0; i < numberOfStimuli; i++)
    {
        // Get the stimulus and update the discretized limits
        Stimulus *stim = &config->stimuli[i];
        if (stim == NULL)
        {
            ERRORMSG("Error in stimulus %d: stimulus is NULL\n", i);
            return -1;
        }
        if (stim->x_range.min >= stim->x_range.max)
        {
            ERRORMSG("Error in stimulus %d: x_min >= x_max\n", i);
            return -1;
        }
        if (stim->y_range.min >= stim->y_range.max)
        {
            ERRORMSG("Error in stimulus %d: y_min >= y_max\n", i);
            return -1;
        }

        stim->x_discretized.max = round(stim->x_range.max / delta_x);
        stim->x_discretized.min = round(stim->x_range.min / delta_x);
        stim->y_discretized.max = round(stim->y_range.max / delta_y);
        stim->y_discretized.min = round(stim->y_range.min / delta_y);
    }

    // Check if the stimulus is within the simulation domain
    for (int i = 0; i < numberOfStimuli; i++)
    {
        Stimulus *stim = &config->stimuli[i];
        if (stim == NULL)
        {
            ERRORMSG("Error in stimulus %d: stimulus is NULL\n", i);
            return -1;
        }
        if (stim->x_discretized.min < 0 || stim->x_discretized.max >= config->Nx)
        {
            ERRORMSG("Error in stimulus %d: x_discretized out of bounds\n", i);
            return -1;
        }
        if (stim->y_discretized.min < 0 || stim->y_discretized.max >= config->Ny)
        {
            ERRORMSG("Error in stimulus %d: y_discretized out of bounds\n", i);
            return -1;
        }
    }

    return 0;
}

void saveCopyOfSimulationConfig(const char *ini_file_path, const char *output_dir)
{
    // Copy the configuration file to the output directory
    char command[MAX_STRING_SIZE];
    snprintf(command, sizeof(command), "cp %s %s/simulation_config.ini", ini_file_path, output_dir);
    int result = system(command);
    if (result != 0)
    {
        ERRORMSG("Error copying simulation configuration file\n");
    }
}

// Basic Thomas algorithm for solving tridiagonal systems
void tridiagonalSystemSolver_x(const int Nx, real *rhs, real *result, real *c_prime, real *d_prime, const real phi_x, const real *Dxx, const int sys_idx)
{
    real denom;
    real coeff = phi_x * Dxx[sys_idx];
    real lalc = -coeff;
    real lb = 1.0f + coeff; // Coefficient for the first element

    c_prime[0] = lalc / lb;
    d_prime[0] = rhs[0] / lb;

    for (int i = 1; i < Nx; i++)
    {
        coeff = phi_x * Dxx[sys_idx + i];
        
        lalc = -coeff;
        lb = (i < Nx - 1) ? 1.0f + 2.0f * coeff : 1.0f + coeff; // Last element has a different coefficient
        denom = 1.0f / (lb - c_prime[i - 1] * lalc);
        
        c_prime[i] = lalc * denom;
        d_prime[i] = (rhs[i] - d_prime[i - 1] * lalc) * denom;
    }

    result[Nx - 1] = d_prime[Nx - 1];

    for (int i = Nx - 2; i >= 0; i--)
    {
        result[i] = d_prime[i] - c_prime[i] * result[i + 1];
    }
}

void tridiagonalSystemSolver_y(const int Ny, real *rhs, real *result, real *c_prime, real *d_prime, const real phi_y, const real *Dyy, const int sys_idx, const int Nx)
{
    real denom;
    real coeff = phi_y * Dyy[sys_idx];
    real lalc = -coeff;
    real lb = 1.0f + coeff; // Coefficient for the first element

    c_prime[0] = lalc / lb;
    d_prime[0] = rhs[0] / lb;

    for (int i = 1; i < Ny; i++)
    {
        coeff = phi_y * Dyy[sys_idx + Nx * i];
        
        lalc = -coeff;
        lb = (i < Ny - 1) ? 1.0f + 2.0f * coeff : 1.0f + coeff; // Last element has a different coefficient
        denom = 1.0f / (lb - c_prime[i - 1] * lalc);
        
        c_prime[i] = lalc * denom;
        d_prime[i] = (rhs[i] - d_prime[i - 1] * lalc) * denom;
    }

    result[Ny - 1] = d_prime[Ny - 1];

    for (int i = Ny - 2; i >= 0; i--)
    {
        result[i] = d_prime[i] - c_prime[i] * result[i + 1];
    }
}

int saveSimulationInfos(const SimulationConfig *config, const Measurement *measurement)
{
    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE, "%s/simulation_infos.txt", config->output_dir);

    FILE *fpInfos = fopen(infosFilePath, "w");
    if (fpInfos == NULL)
    {
        ERRORMSG("Error opening file %s\n", infosFilePath);
        return -1;
    }

    fprintf(fpInfos, "EXECUTION TYPE = %s\n", executionModeToString(config->exec_mode));
    fprintf(fpInfos, "PRECISION = %s\n", REAL_TYPE);
    fprintf(fpInfos, "PROBLEM = %s\n", equationTypeToString(config->equation_type));
    fprintf(fpInfos, "CELL MODEL = %s\n", cellModelToString(config->cell_model));
    fprintf(fpInfos, "METHOD = %s\n", numericalMethodToString(config->method));
    fprintf(fpInfos, "\n");

    fprintf(fpInfos, "DOMAIN LENGTH IN X = %.4g cm\n", config->Lx);
    fprintf(fpInfos, "DOMAIN LENGTH IN Y = %.4g cm\n", config->Ly);
    fprintf(fpInfos, "TOTAL TIME = %.4g ms\n", config->total_time);
    fprintf(fpInfos, "\n");

    fprintf(fpInfos, "delta_t = %.5g ms (%d time steps)\n", config->dt, config->M);
    fprintf(fpInfos, "delta_x = %.5g cm (%d um) (%d space steps in x)\n", config->dx, CM_TO_UM(config->dx), config->Nx);
    fprintf(fpInfos, "delta_y = %.5g cm (%d um) (%d space steps in y)\n", config->dy, CM_TO_UM(config->dy), config->Ny);
    fprintf(fpInfos, "TOTAL POINTS IN DOMAIN = %d\n", config->Nx * config->Ny);
    fprintf(fpInfos, "SIGMA LONGITUDINAL = %.8g\n", config->sigma_l);
    fprintf(fpInfos, "SIGMA TRANSVERSAL = %.8g\n", config->sigma_t);
    fprintf(fpInfos, "FIBER ORIENTATION = %.2f degrees\n", config->fiber_orientation);
    fprintf(fpInfos, "NUMBER OF STIMULI = %d\n", config->stimulus_count);
    for (int i = 0; i < config->stimulus_count; i++)
    {
        fprintf(fpInfos, "STIMULUS %d: START TIME = %.5g ms\n", i + 1, config->stimuli[i].start_time);
    }

    if (config->measure_velocity)
    {
        fprintf(fpInfos, "\n");
        fprintf(fpInfos, "STIMULUS S1 VELOCITY = %.4g cm/s\n", measurement->stimVelocity);
    }

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "TIME OF THE FIRST PART = %.5g s\n", measurement->elapsedTime1stPart);
    fprintf(fpInfos, "TIME OF THE SECOND PART = %.5g s\n", measurement->elapsedTime2ndPart);

    if (config->method == METHOD_OSADI || config->method == METHOD_SSIADI || config->method == METHOD_DO || config->method == METHOD_HV)
    {
        fprintf(fpInfos, "TIME TO SOLVE THE 1st LINEAR SYSTEM = %.5g s\n", measurement->elapsedTime1stLS);
        fprintf(fpInfos, "TIME TO SOLVE THE 2nd LINEAR SYSTEM = %.5g s\n", measurement->elapsedTime2ndLS);
    }

    if (config->method == METHOD_HV)
    {
        fprintf(fpInfos, "TIME TO COMPUTE THE APPROXIMATION = %.5g s\n", measurement->elapsedTimeApproximation);
        fprintf(fpInfos, "TIME TO SOLVE THE 3rd LINEAR SYSTEM = %.5g s\n", measurement->elapsedTime3rdLS);
        fprintf(fpInfos, "TIME TO SOLVE THE 4th LINEAR SYSTEM = %.5g s\n", measurement->elapsedTime4thLS);
    }

    if (config->measure_velocity)
        fprintf(fpInfos, "TIME TO MEASURE VELOCITY = %.5g s\n", measurement->elapsedMeasureVelocityTime);

    if (config->save_frames)
        fprintf(fpInfos, "TIME TO SAVE FRAMES = %.5g s\n", measurement->elapsedSaveFramesTime);

    if (config->save_last_state)
        fprintf(fpInfos, "TIME TO SAVE LAST STATE = %.5g s\n", measurement->elapsedSaveStateTime);

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "SIMULATION TOTAL EXECUTION TIME = %.5g s\n", measurement->elapsedExecutionTime);

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "OUTPUT DIRECTORY = %s\n", config->output_dir);
    fclose(fpInfos);

    SIMPLEMSG("");
    INFOMSG("Simulation total execution time = %.5g s\n", measurement->elapsedExecutionTime);
    SUCCESSMSG("Simulation infos saved to %s\n", infosFilePath);

    return 0;
}

// Exact solution when using forcing term
real exactSolution(real t, real x, real y, real Lx, real Ly)
{
    return (exp(-t)) * cos(_PI * x / Lx) * cos(_PI * y / Ly);
}

real calculateNorm2Error(real *Vm, real *exact, int Nx, int Ny, real totalTime, real delta_x, real delta_y, real Lx, real Ly)
{
    real x, y;
    real solution;
    real sum = 0.0;
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            x = j * delta_x;
            y = i * delta_y;
            solution = exactSolution(totalTime, x, y, Lx, Ly);
            exact[index] = solution;
            sum += ((Vm[index] - solution) * (Vm[index] - solution));
        }
    }
    return sqrt(sum / (Nx * Ny));
}

void initializeProperties(const SimulationConfig *config, real *Dxx, real *Dyy, real *Dxy)
{
    // Unpack parameters from the config
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real sigma_l = config->sigma_l;
    const real sigma_t = config->sigma_t;
    const real fiber_orientation_rad = config->fiber_orientation * _PI / 180.0f;

    for (int idx = 0; idx < Nx * Ny; idx++)
    {
        // Calculate the diffusion coefficients based on the fiber orientation
        // D_xx = (sigma_l - sigma_t) * cos^2(theta) + sigma_t
        // D_yy = (sigma_l - sigma_t) * sin^2(theta) + sigma_t
        // D_xy = (sigma_l - sigma_t) * sin(theta) * cos(theta)
        Dxx[idx] = (sigma_l - sigma_t) * cos(fiber_orientation_rad) * cos(fiber_orientation_rad) + sigma_t;
        Dyy[idx] = (sigma_l - sigma_t) * sin(fiber_orientation_rad) * sin(fiber_orientation_rad) + sigma_t;
        Dxy[idx] = (sigma_l - sigma_t) * sin(fiber_orientation_rad) * cos(fiber_orientation_rad);
    }
}

void initializeVariableWithExactSolution(real *Var, int Nx, int Ny, real delta_x, real delta_y, real Lx, real Ly)
{
    real x, y;
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            y = i * delta_y;
            x = j * delta_x;
            Var[index] = exactSolution(0.0f, x, y, Lx, Ly);
        }
    }
}

void initializeVariableWithValue(real *Var, int Nx, int Ny, real value)
{
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            Var[index] = value;
        }
    }
}

void initializeVariableFromFile(real *Var, int Nx, int Ny, char *filename, real delta_x, real delta_y, char *varName, real reference_dx, real reference_dy, real Lx, real Ly)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        ERRORMSG("Error opening file %s\n", filename);
        exit(1);
    }

    int baseNx = round(Lx / reference_dx) + 1;
    int baseNy = round(Ly / reference_dy) + 1;
    int rate_x = round(delta_x / reference_dx);
    int rate_y = round(delta_y / reference_dy);

    int sizeFile = 0;
    int sizeVar = 0;
    real value;

    INFOMSG("Reading file %s to initialize variable with rate_x=%d and rate_y=%d\n", filename, rate_x, rate_y);
    for (int i = 0; i < baseNy; i++)
    {
        for (int j = 0; j < baseNx; j++)
        {
            // Read value from file to variable
            // If i and j are multiples of rate, read value to Var
            fscanf(file, FSCANF_REAL, &value);
            if (i % rate_y == 0 && j % rate_x == 0)
            {
                Var[sizeVar] = value;
                if (isnan(value))
                {
                    ERRORMSG("At var index %d, file index %d, value is NaN\n", sizeVar, i * baseNx + j);
                    exit(1);
                }
                sizeVar++;
            }
            sizeFile++;
        }
    }
    fclose(file);
    INFOMSG("Variable %s initialized with %d values from the %d values in file\n", varName, sizeVar, sizeFile);
}

void shiftVariableToLeft(real *Var, int Nx, int Ny, real length, real delta_x, real delta_y, real initValue, char *varName)
{
    real *temp = (real *)malloc(Nx * sizeof(real));
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            temp[j] = Var[index];
        }

        int lengthIndex = round(length / delta_x) + 1;
        for (int j = 0; j < Nx - lengthIndex; j++)
        {
            index = i * Nx + j;
            Var[index] = temp[j + lengthIndex];
        }
        for (int j = Nx - lengthIndex; j < Nx; j++)
        {
            index = i * Nx + j;
            Var[index] = initValue;
        }
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

#ifdef CABLEEQ

void initialize1DVariableWithValue(real *Var, int N, real value)
{
    for (int i = 0; i < N; i++)
        Var[i] = value;
}

void initialize1DVariableFromFile(real *Var, int N, char *filename, real delta_x, char *varName, real reference_dx, real Lx)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        ERRORMSG("Error opening file %s\n", filename);
        exit(1);
    }

    int baseN = round(Lx / reference_dx) + 1;
    int rate = round(delta_x / reference_dx);

    int sizeFile = 0;
    int sizeVar = 0;
    real value;

    INFOMSG("Reading file %s to initialize variable with a rate of %d\n", filename, rate);

    for (int i = 0; i < baseN; i++)
    {
        // Read value from file to variable
        // If i is multiple of rate, read value to Var
        fscanf(file, FSCANF_REAL, &value);
        if (i % rate == 0)
        {
            Var[(int)(i / rate)] = value;
            if (isnan(value))
            {
                ERRORMSG("At var index [%d], file index %d, value is NaN\n", (int)(i / rate), i);
                exit(1);
            }
            sizeVar++;
        }
        sizeFile++;
    }
    fclose(file);

    SUCCESSMSG("Variable %s initialized with %d values from the %d values in file\n", varName, sizeVar, sizeFile);
}

void shift1DVariableToLeft(real *Var, int N, real length, real delta_x, real initValue, char *varName)
{
    real *temp = (real *)malloc(N * sizeof(real));
    for (int i = 0; i < N; i++)
    {
        temp[i] = Var[i];
    }

    int lengthIndex = round(length / delta_x) + 1;
    for (int i = 0; i < N - lengthIndex; i++)
    {
        Var[i] = temp[i + lengthIndex];
    }
    for (int i = N - lengthIndex; i < N; i++)
    {
        Var[i] = initValue;
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

#endif // CABLEEQ