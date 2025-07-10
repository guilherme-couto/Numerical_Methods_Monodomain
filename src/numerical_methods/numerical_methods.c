#include "numerical_methods.h"

// Map of numerical methods and their corresponding solvers
static const struct
{
    const NumericalMethod method;
    numerical_method_t solver_SERIAL;
    numerical_method_t solver_OMP;
    numerical_method_t solver_CUDA;
} numerical_method_map[] = {
    {METHOD_SSIADI, runSSIADI, runSSIADI_OMP, runSSIADI_CUDA},
    {METHOD_OSADI, runOSADI, NULL, runOSADI_CUDA},
    {METHOD_FE, runFE, NULL, runFE_CUDA},
    {METHOD_DO, runDO, NULL, runDO_CUDA},
    {METHOD_HV, NULL, NULL, runHV_CUDA},
    {METHOD_CG, runCG, NULL, runCG_CUDA},
    {METHOD_INVALID, NULL, NULL, NULL}};

numerical_method_t get_numerical_method(const ExecutionMode *exec_mode, const NumericalMethod *method)
{
    if (method == NULL)
    {
        printf("Error: numerical method is NULL\n");
        return NULL;
    }

    for (int i = 0; numerical_method_map[i].method != METHOD_INVALID; i++)
        if (numerical_method_map[i].method == *method)
        {
            if (*exec_mode == EXEC_SERIAL)
                return numerical_method_map[i].solver_SERIAL;
            else if (*exec_mode == EXEC_OPENMP)
                return numerical_method_map[i].solver_OMP;
            else if (*exec_mode == EXEC_CUDA)
                return numerical_method_map[i].solver_CUDA;
        }
    return NULL;
}