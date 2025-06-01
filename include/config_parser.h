#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include "core_definitions.h"
#include "save_functions.h"

#ifdef __cplusplus
extern "C"
{
#endif

    // Execution modes
    typedef enum
    {
        EXEC_SERIAL,
        EXEC_OPENMP,
        EXEC_CUDA,
        EXEC_INVALID
    } ExecutionMode;

    // Equation types
    typedef enum
    {
        EQUATION_MONODOMAIN,
        EQUATION_INVALID
    } EquationType;

    // Cell models
    typedef enum
    {
        CELL_MODEL_AFHN,
        CELL_MODEL_TT2,
        CELL_MODEL_MV,
        CELL_MODEL_INVALID
    } CellModel;

    // Numerical methods
    typedef enum
    {
        METHOD_OSADI,
        METHOD_SSIADI,
        METHOD_FE,
        METHOD_DO,
        METHOD_HV,
        METHOD_INVALID
    } NumericalMethod;

    typedef struct
    {
        // Basic parameters
        ExecutionMode exec_mode;
        EquationType equation_type;
        CellModel cell_model;
        NumericalMethod method;

        // Numeric parameters
        real dt;
        real dx;
        real dy;

        real sigma_l;           // Longitudinal diffusion coefficient
        real sigma_t;           // Transversal diffusion coefficient
        real fiber_orientation; // Fiber angle with respect to the x-axis (in degrees)

        // Simulation parameters
        real total_time;
        real Lx;
        real Ly;
        int M;  // Number of time steps
        int Nx; // Number of spatial steps in x
        int Ny; // Number of spatial steps in y

        // Stimuli parameters
        Stimulus *stimuli;
        int stimulus_count;

        int frame_save_rate;
        int number_of_threads;

        char output_dir[MAX_STRING_SIZE];
        bool remove_old_files;
        char path_to_restore_state_files[MAX_STRING_SIZE];

        char save_function_name[MAX_STRING_SIZE];
        save_function_t save_function;
        char file_extension[4];

        // Boolean flags (using stdbool.h)
        bool shift_state;
        bool save_frames;
        bool save_last_frame;
        bool save_last_state;
        bool measure_velocity;

        // String parameters
        char init_mode[20];
    } SimulationConfig;

    int load_simulation_config(const char *filename, SimulationConfig *config);
    bool validate_simulation_config(const SimulationConfig *config);
    void free_simulation_config(SimulationConfig *config);

    const char *executionModeToString(ExecutionMode mode);
    const char *equationTypeToString(EquationType type);
    const char *cellModelToString(CellModel model);
    const char *numericalMethodToString(NumericalMethod method);

#ifdef __cplusplus
}
#endif

#endif // CONFIG_PARSER_H