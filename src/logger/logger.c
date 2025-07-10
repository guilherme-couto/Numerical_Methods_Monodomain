#include "logger.h"
#include <time.h>

FILE *log_file;

void init_logger(const char *filename)
{
    log_file = fopen(filename, "w");
    if (!log_file)
    {
        fprintf(stderr, "Error opening log file '%s'\n", filename);
        return;
    }
}

void close_logger(void)
{
    if (log_file)
    {
        fclose(log_file);
        log_file = NULL;
    }
}

void log_message(const char *color, const char *prefix, const char *fmt, ...)
{
    va_list args;

    // Terminal (with colors)
    printf("%s%s%s", color, prefix, RESET);
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    // Log file (without colors)
    if (log_file)
    {
        fprintf(log_file, "%s", prefix);
        va_start(args, fmt);
        vfprintf(log_file, fmt, args);
        va_end(args);
    }
}

static const char *get_current_timestamp()
{
    static char buffer[32];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", t);
    return buffer;
}

void log_simulation_header(const SimulationConfig *config, const char *config_filename)
{
    if (!log_file || !config)
        return;

    fprintf(log_file, "================================================================\n");
    fprintf(log_file, "                         Simulation Log\n");
    fprintf(log_file, "----------------------------------------------------------------\n");
    fprintf(log_file, " Start Time             : %s\n", get_current_timestamp());
    fprintf(log_file, " Configuration File     : %s\n", config_filename);
    fprintf(log_file, " Output Directory       : %s\n", config->output_dir);
    fprintf(log_file, "\n");

    fprintf(log_file, " Execution Mode         : %s\n", executionModeToString(config->exec_mode));
    fprintf(log_file, " Equation Type          : %s\n", equationTypeToString(config->equation_type));
    fprintf(log_file, " Cell Model             : %s\n", cellModelToString(config->cell_model));
    fprintf(log_file, " Numerical Method       : %s\n", numericalMethodToString(config->method));
    fprintf(log_file, " Real Precision         : %s\n", REAL_TYPE);
    fprintf(log_file, "\n");

    fprintf(log_file, " Save Frames            : %s\n", config->save_frames ? "Yes" : "No");
    fprintf(log_file, " Save Last Frame        : %s\n", config->save_last_frame ? "Yes" : "No");
    fprintf(log_file, " Save Last State        : %s\n", config->save_state ? "Yes" : "No");
    fprintf(log_file, " Shift State            : %s\n", config->shift_state ? "Yes" : "No");
    fprintf(log_file, " Measure Velocity       : %s\n", config->measure_velocity ? "Yes" : "No");
    if (config->save_frames)
        fprintf(log_file, " Frame Save Rate        : %d\n", config->frame_save_rate);
    if (config->exec_mode == EXEC_OPENMP)
        fprintf(log_file, " Number of Threads      : %d\n", config->number_of_threads);
    if (config->restore_path[0] != '\0')
        fprintf(log_file, " Restore State Path     : %s\n", config->restore_path);
    fprintf(log_file, " Initial State Mode     : %s\n", config->init_mode);
    fprintf(log_file, "\n");

    fprintf(log_file, " Grid Size (Nx x Ny)    : %d x %d\n", config->Nx, config->Ny);
    fprintf(log_file, " Domain Size (Lx x Ly)  : %.2f x %.2f\n", config->Lx, config->Ly);
    fprintf(log_file, " Spatial Step (dx, dy)  : %.4g, %.4g\n", config->dx, config->dy);
    fprintf(log_file, " Time Step (dt)         : %.6g\n", config->dt);
    fprintf(log_file, " Total Time             : %.2f\n", config->total_time);
    fprintf(log_file, " Number of Steps (M)    : %d\n", config->M);
    fprintf(log_file, " Fiber Orientation      : %.2f degrees\n", config->fiber_orientation);
    fprintf(log_file, " Is Fiber Aligned       : %s\n", config->is_fiber_aligned ? "Yes" : "No");
    fprintf(log_file, " Sigma Longitudinal     : %.8g\n", config->sigma_l);
    fprintf(log_file, " Sigma Transversal      : %.8g\n", config->sigma_t);
    fprintf(log_file, "\n");

    for (int i = 0; i < config->stimulus_count; i++)
    {
        fprintf(log_file, " Stimulus %d:\n", i + 1);
        fprintf(log_file, "   Amplitude  : %.4g\n", config->stimuli[i].amplitude);
        fprintf(log_file, "   Duration   : %.4g\n", config->stimuli[i].duration);
        fprintf(log_file, "   Start Time : %.4g\n", config->stimuli[i].start_time);
        fprintf(log_file, "   X Range    : [%.4g, %.4g]\n", config->stimuli[i].x_range.min, config->stimuli[i].x_range.max);
        fprintf(log_file, "   Y Range    : [%.4g, %.4g]\n", config->stimuli[i].y_range.min, config->stimuli[i].y_range.max);
        if (i < config->stimulus_count - 1)
            fprintf(log_file, "\n");
    }

    fprintf(log_file, "================================================================\n\n");

    fflush(log_file);
}

void log_machine_info()
{
    if (!log_file)
        return;

    fprintf(log_file, "================================================================\n");
    fprintf(log_file, "                         Machine Info\n");
    fprintf(log_file, "----------------------------------------------------------------\n");

#ifdef _WIN32
    fprintf(log_file, " Machine information not supported on Windows systems.\n");
#else
    struct utsname uname_info;
    if (uname(&uname_info) != 0)
    {
        fprintf(log_file, " Failed to retrieve system information.\n");
        return;
    }

    char hostname[256];
    if (gethostname(hostname, sizeof(hostname)) != 0)
        snprintf(hostname, sizeof(hostname), "Unknown");

    long num_cores = sysconf(_SC_NPROCESSORS_ONLN);

    struct sysinfo s_info;
    if (sysinfo(&s_info) != 0)
    {
        fprintf(log_file, " Failed to retrieve memory info.\n");
        return;
    }

    double total_mem_gb = s_info.totalram / (1024.0 * 1024.0 * 1024.0);
    double free_mem_gb = s_info.freeram / (1024.0 * 1024.0 * 1024.0);

    fprintf(log_file, " Hostname              : %s\n", hostname);
    fprintf(log_file, " OS                    : %s\n", uname_info.sysname);
    fprintf(log_file, " Architecture          : %s\n", uname_info.machine);
    fprintf(log_file, " CPU Cores             : %ld\n", num_cores);
    fprintf(log_file, " Total Memory (GB)     : %.2f\n", total_mem_gb);
    fprintf(log_file, " Free Memory (GB)      : %.2f\n", free_mem_gb);
#endif

#ifdef USE_OPENMP
    int max_threads = omp_get_max_threads();
    int num_procs   = omp_get_num_procs();
    int nested      = omp_get_nested();

    fprintf(log_file, "---------------------------- OpenMP -----------------------------\n");
    fprintf(log_file, " Available Threads     : %d\n", num_procs);
    fprintf(log_file, " Max Threads Allowed   : %d\n", max_threads);
    fprintf(log_file, " Nested Parallelism    : %s\n", nested ? "Enabled" : "Disabled");
#endif

    fprintf(log_file, "================================================================\n\n");

    fflush(log_file);
}