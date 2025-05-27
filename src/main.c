#include "../include/config_parser.h"
#include "../include/core_definitions.h"
#include "../include/auxfuncs.h"
#include "../src/solvers/monodomain.h"
#include "../src/logger/logger.h"
#include <signal.h>

SimulationConfig config;
bool config_loaded = false;

void handle_sigint(int sig)
{
    ERRORMSG("Simulation interrupted by user (SIGINT)\n");

    if (config_loaded)
    {
        free_simulation_config(&config);
    }

    close_logger();
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    // Handle SIGINT (Ctrl+C)
    signal(SIGINT, handle_sigint);

    // Verify arguments
    if (argc < 2)
    {
        printf("Usage: %s <configuration_file>", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Then proceed with normal loading
    int config_status = load_simulation_config(argv[1], &config);
    if (config_status != 0)
    {
        printf("Failed to load configuration (error code: %d)\n", config_status);
        free_simulation_config(&config);
        return -1;
    }
    config_loaded = true;

    // Create directories
    if (createDirectories(config.output_dir, config.remove_old_files) != 0)
    {
        printf("Error creating directories\n");
        free_simulation_config(&config);
        return -1;
    }

    // Populate the stimuli array
    if (populateStimuli(&config) != 0)
    {
        printf("Failed to populate stimuli\n");
        free_simulation_config(&config);
        return -1;
    }

    // Initialize the logger
    char log_filename[MAX_STRING_SIZE];
    snprintf(log_filename, sizeof(log_filename), "%s/output.log", config.output_dir);
    init_logger(log_filename);
    atexit(close_logger); // Ensure logger is closed on exit

    // Log the machine info and simulation header
    log_machine_info();
    #ifdef USE_CUDA
    log_device_info(log_file);
    #endif // USE_CUDA
    log_simulation_header(&config, argv[1]);

    // Save a copy of the configuration file
    saveCopyOfSimulationConfig(argv[1], config.output_dir);

    // Check the execution mode and the equation type to determine the simulation type
    int simulation_result = -1;
    if (config.exec_mode == EXEC_SERIAL)
    {
        if (config.equation_type == EQUATION_MONODOMAIN)
        {
            simulation_result = runMonodomainSimulationSerial(&config);
        }
        else
        {
            ERRORMSG("Unsupported equation type");
            free_simulation_config(&config);
            return -1;
        }
    }

#ifdef USE_CUDA

    else if (config.exec_mode == EXEC_CUDA)
    {
        if (config.equation_type == EQUATION_MONODOMAIN)
        {
            simulation_result = runMonodomainSimulationCUDA(&config);
        }
        else
        {
            ERRORMSG("Unsupported equation type");
            free_simulation_config(&config);
            return -1;
        }
    }

#endif // USE_CUDA

    // Check if the simulation was successful
    if (simulation_result != 0)
    {
        ERRORMSG("Simulation failed with error code: %d\n", simulation_result);
        free_simulation_config(&config);
        return -1;
    }

    // Free the configuration
    free_simulation_config(&config);

    return EXIT_SUCCESS;
}
