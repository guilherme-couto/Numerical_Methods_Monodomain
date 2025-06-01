#include "../include/config_parser.h"
#include "../external/inih/ini.h"
#include "../external/tinyexpr/tinyexpr.h"

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

// Convert string to lowercase for case-insensitive comparison
static void strtolower(char *str)
{
    for (; *str; ++str)
        *str = tolower(*str);
}

// Function to evaluate expressions using tinyexpr and not considering inline comments
static double te_interp_no_comments(const char *expression, int *error)
{
    // Remove inline comments beginning with ';' or '#'
    char *comment_start = strpbrk(expression, ";#");
    if (comment_start)
    {
        size_t length = comment_start - expression;
        char *clean_expression = (char *)malloc(length + 1);
        if (!clean_expression)
        {
            fprintf(stderr, "Memory allocation error\n");
            return NAN;
        }
        strncpy(clean_expression, expression, length);
        clean_expression[length] = '\0';

        // Evaluate the expression
        double result = te_interp(clean_expression, error);
        free(clean_expression);
        return result;
    }
    else
    {
        // No comments found, evaluate the expression directly
        return te_interp(expression, error);
    }
}

static int config_parser_handler(void *user, const char *section, const char *name, const char *value)
{
    SimulationConfig *config = (SimulationConfig *)user;
    char lower_value[MAX_STRING_SIZE];

    // Make case-insensitive comparison easier
    strncpy(lower_value, value, sizeof(lower_value));
    strtolower(lower_value);

    // Handle the "simulation" section
    if (MATCH("simulation", "execution_mode"))
    {
        if (strstr(lower_value, "gpu") || strstr(lower_value, "cuda"))
            config->exec_mode = EXEC_CUDA;
        else if (strstr(lower_value, "openmp"))
            config->exec_mode = EXEC_OPENMP;
        else if (strstr(lower_value, "serial"))
            config->exec_mode = EXEC_SERIAL;
        else
            config->exec_mode = EXEC_INVALID;
    }
    else if (MATCH("simulation", "equation_type"))
    {
        if (strstr(lower_value, "monodomain"))
            config->equation_type = EQUATION_MONODOMAIN;
        else
            config->equation_type = EQUATION_INVALID;
    }
    else if (MATCH("simulation", "cell_model"))
    {
        if (strstr(lower_value, "afhn"))
            config->cell_model = CELL_MODEL_AFHN;
        else if (strstr(lower_value, "tt2"))
            config->cell_model = CELL_MODEL_TT2;
        else if (strstr(lower_value, "mv"))
            config->cell_model = CELL_MODEL_MV;
        else
            config->cell_model = CELL_MODEL_INVALID;
    }
    else if (MATCH("simulation", "numerical_method"))
    {
        if (strstr(lower_value, "osadi"))
            config->method = METHOD_OSADI;
        else if (strstr(lower_value, "ssiadi"))
            config->method = METHOD_SSIADI;
        else if (strstr(lower_value, "fe"))
            config->method = METHOD_FE;
        else if (strstr(lower_value, "do"))
            config->method = METHOD_DO;
        else if (strstr(lower_value, "hv"))
            config->method = METHOD_HV;
        else
            config->method = METHOD_INVALID;
    }
    else if (MATCH("simulation", "save_function"))
    {
        strncpy(config->save_function_name, value, sizeof(config->save_function_name));
        config->save_function = get_save_function(config->save_function_name);
        strncpy(config->file_extension, get_file_extension(config->save_function_name), sizeof(config->file_extension));
    }

    else if (MATCH("simulation", "dt"))
        config->dt = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "dx"))
        config->dx = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "dy"))
        config->dy = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "sigma_l"))
        config->sigma_l = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "sigma_t"))
        config->sigma_t = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "fiber_orientation"))
        config->fiber_orientation = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "total_time"))
        config->total_time = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "length_x"))
        config->Lx = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "length_y"))
        config->Ly = te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "save_rate"))
        config->frame_save_rate = (int)te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "number_of_threads"))
        config->number_of_threads = (int)te_interp_no_comments(value, 0);
    else if (MATCH("simulation", "output_dir"))
        strncpy(config->output_dir, value, sizeof(config->output_dir));
    else if (MATCH("simulation", "remove_old_files"))
        config->remove_old_files = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "path_to_restore_state_files"))
    {
        strncpy(config->path_to_restore_state_files, value, sizeof(config->path_to_restore_state_files));
        strncpy(config->init_mode, "Restore State", sizeof(config->init_mode));
    }
    else if (MATCH("simulation", "shift_state"))
        config->shift_state = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_frames"))
        config->save_frames = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_last_frame"))
        config->save_last_frame = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_last_state"))
        config->save_last_state = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "measure_velocity"))
        config->measure_velocity = (strstr(lower_value, "true") != NULL);

    // Handle the "stimuli" section
    else if (strncmp(section, "stim_", 5) == 0)
    {
        int stim_index = atoi(section + 5) - 1;

        // Validate index
        if (stim_index < 0)
        {
            printf("Invalid stimulus index in section: %s -> Stimulus index must be a positive integer\n", section);
            return 0;
        }

        // Resize the stimuli array if necessary
        if (stim_index >= config->stimulus_count)
        {
            Stimulus *new_stim = realloc(config->stimuli, (stim_index + 1) * sizeof(Stimulus));
            if (!new_stim)
            {
                printf("Failed to allocate memory for stimuli\n");
                return 0;
            }
            config->stimuli = new_stim;

            // Initialize new stimulus
            memset(&config->stimuli[stim_index], 0, sizeof(Stimulus));
            config->stimulus_count = stim_index + 1;
        }

        // Set stimulus parameters
        Stimulus *stim = &config->stimuli[stim_index];
        if (strcmp(name, "amplitude") == 0)
            stim->amplitude = te_interp_no_comments(value, 0);
        else if (strcmp(name, "start_time") == 0)
            stim->start_time = te_interp_no_comments(value, 0);
        else if (strcmp(name, "duration") == 0)
            stim->duration = te_interp_no_comments(value, 0);
        else if (strcmp(name, "x_min") == 0)
            stim->x_range.min = te_interp_no_comments(value, 0);
        else if (strcmp(name, "x_max") == 0)
            stim->x_range.max = te_interp_no_comments(value, 0);
        else if (strcmp(name, "y_min") == 0)
            stim->y_range.min = te_interp_no_comments(value, 0);
        else if (strcmp(name, "y_max") == 0)
            stim->y_range.max = te_interp_no_comments(value, 0);
        else
            printf("Unknown configuration parameter: %s.%s -> IGNORING...\n", section, name);
    }
    return 1;
}

int load_simulation_config(const char *filename, SimulationConfig *config)
{
    printf("Loading config from: %s\n", filename);

    memset(config, 0, sizeof(SimulationConfig));

    // Initialize default values
    config->exec_mode = EXEC_INVALID;
    config->equation_type = EQUATION_INVALID;
    config->cell_model = CELL_MODEL_INVALID;
    config->method = METHOD_INVALID;
    config->dt = -1.0f;
    config->dx = -1.0f;
    config->dy = -1.0f;
    config->sigma_l = -1.0f;
    config->sigma_t = -1.0f;
    config->fiber_orientation = 0.0f;
    config->total_time = -1.0;
    config->Lx = -1.0f;
    config->Ly = -1.0f;
    config->frame_save_rate = -1;
    config->number_of_threads = -1;
    config->shift_state = false;
    config->save_frames = true;
    config->save_last_frame = true;
    config->save_last_state = false;
    config->measure_velocity = true;
    config->stimuli = NULL;
    config->stimulus_count = -1;
    config->save_function = NULL;
    config->save_function_name[0] = '\0';
    config->file_extension[0] = '\0';
    config->output_dir[0] = '\0';
    config->remove_old_files = true;
    config->path_to_restore_state_files[0] = '\0';
    strncpy(config->init_mode, "Initial Condition", sizeof(config->init_mode));

    int result = ini_parse(filename, config_parser_handler, config);
    if (result < 0)
    {
        printf("Failed to open config file: %s\n", filename);
        return -1;
    }
    if (result > 0)
    {
        printf("Error parsing config file at line %d\n", result);
        return -2;
    }

    if (validate_simulation_config(config) == false)
    {
        return -3;
    }

    // Number of steps
    int M = round(config->total_time / config->dt); // Number of time steps
    int Nx = round(config->Lx / config->dx) + 1;    // Spatial steps in x
    int Ny = round(config->Ly / config->dy) + 1;    // Spatial steps in y
    config->M = M;
    config->Nx = Nx;
    config->Ny = Ny;

    printf("Configuration loaded successfully\n");

    return 0;
}

bool validate_simulation_config(const SimulationConfig *config)
{
    bool valid = true;

    if (config->exec_mode == EXEC_INVALID)
    {
        printf("Invalid execution mode specified\n");
        valid = false;
    }
    if (config->equation_type == EQUATION_INVALID)
    {
        printf("Invalid equation type specified\n");
        valid = false;
    }
    if (config->cell_model == CELL_MODEL_INVALID)
    {
        printf("Invalid cell model specified\n");
        valid = false;
    }
    if (config->method == METHOD_INVALID)
    {
        printf("Invalid numerical method specified\n");
        valid = false;
    }
    if (config->dt <= 0)
    {
        printf("Time step (dt) must be positive\n");
        valid = false;
    }
    if (config->dx <= 0 || config->dy <= 0)
    {
        printf("Spatial steps (dx, dy) must be positive\n");
        valid = false;
    }
    if (config->sigma_l <= 0 || config->sigma_t <= 0)
    {
        printf("Diffusion coefficients (sigma_l, sigma_t) must be positive\n");
        valid = false;
    }
    if (config->fiber_orientation < 0 || config->fiber_orientation > 180)
    {
        printf("Fiber orientation must be between 0 and 180 degrees\n");
        valid = false;
    }
    if (config->total_time <= 0)
    {
        printf("Total time must be positive\n");
        valid = false;
    }
    if (config->Lx <= 0 || config->Ly <= 0)
    {
        printf("Simulation domain dimensions (Lx, Ly) must be positive\n");
        valid = false;
    }
    if (config->save_function == NULL && (config->save_last_frame || config->save_frames))
    {
        printf("Save function must be specified for saving frames\n");
        valid = false;
    }
    if (config->number_of_threads <= 0 && config->exec_mode == EXEC_OPENMP)
    {
        printf("Number of threads must be positive for OpenMP execution\n");
        valid = false;
    }
    if (config->frame_save_rate <= 0 && config->save_frames)
    {
        printf("Frame save rate must be positive for saving frames\n");
        valid = false;
    }
    if (config->output_dir[0] == '\0')
    {
        printf("Output directory must be specified\n");
        valid = false;
    }
    for (int i = 0; i < config->stimulus_count; i++)
    {
        if (config->stimuli[i].start_time < 0)
        {
            printf("Stimulus begin time must be non-negative\n");
            valid = false;
        }
        if (config->stimuli[i].duration <= 0)
        {
            printf("Stimulus duration must be positive\n");
            valid = false;
        }
        if (config->stimuli[i].x_range.min < 0 || config->stimuli[i].x_range.max < 0 ||
            config->stimuli[i].y_range.min < 0 || config->stimuli[i].y_range.max < 0)
        {
            printf("Stimulus spatial range must be non-negative\n");
            valid = false;
        }
    }

    return valid;
}

void free_simulation_config(SimulationConfig *config)
{
    if (config->stimuli)
    {
        free(config->stimuli);
        config->stimuli = NULL;
    }
    config->stimulus_count = 0;
}

const char *executionModeToString(ExecutionMode mode)
{
    switch (mode)
    {
    case EXEC_SERIAL:
        return "SERIAL";
    case EXEC_OPENMP:
        return "OPENMP";
    case EXEC_CUDA:
        return "CUDA";
    case EXEC_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *equationTypeToString(EquationType type)
{
    switch (type)
    {
    case EQUATION_MONODOMAIN:
        return "MONODOMAIN";
    case EQUATION_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *cellModelToString(CellModel model)
{
    switch (model)
    {
    case CELL_MODEL_AFHN:
        return "AFHN";
    case CELL_MODEL_TT2:
        return "TT2";
    case CELL_MODEL_MV:
        return "MV";
    case CELL_MODEL_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *numericalMethodToString(NumericalMethod method)
{
    switch (method)
    {
    case METHOD_OSADI:
        return "OSADI";
    case METHOD_SSIADI:
        return "SSIADI";
    case METHOD_FE:
        return "FE";
    case METHOD_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}