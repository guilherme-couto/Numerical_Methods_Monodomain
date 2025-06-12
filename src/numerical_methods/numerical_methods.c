#include "numerical_methods.h"

// Map of numerical methods and their corresponding solvers
static const struct
{
    const NumericalMethod method;
    numerical_method_t solver;
} numerical_method_map[] = {
    {METHOD_SSIADI, runSSIADI},
    {METHOD_OSADI, runOSADI},
    {METHOD_FE, runFE},
    {METHOD_DO, runDO},
    // {METHOD_HV, runHV}, // HV method is not implemented yet
    {METHOD_CG, runCG},
    {METHOD_INVALID, NULL}};

numerical_method_t get_numerical_method(const NumericalMethod *method)
{
    if (method == NULL)
    {
        printf("Error: numerical method is NULL\n");
        return NULL;
    }

    for (int i = 0; numerical_method_map[i].solver != NULL; i++)
        if (numerical_method_map[i].method == *method)
            return numerical_method_map[i].solver;
    return NULL;
}

// Map of numerical methods and their corresponding solvers for CUDA
static const struct
{
    const NumericalMethod method;
    numerical_method_t solver;
} numerical_method_map_CUDA[] = {
    {METHOD_SSIADI, runSSIADI_CUDA},
    {METHOD_OSADI, runOSADI_CUDA},
    {METHOD_FE, runFE_CUDA},
    {METHOD_DO, runDO_CUDA},
    {METHOD_HV, runHV_CUDA},
    {METHOD_CG, runCG_CUDA},
    {METHOD_INVALID, NULL}};

numerical_method_t get_numerical_method_CUDA(const NumericalMethod *method)
{
    if (method == NULL)
    {
        printf("Error: numerical method is NULL\n");
        return NULL;
    }

    for (int i = 0; numerical_method_map_CUDA[i].solver != NULL; i++)
        if (numerical_method_map_CUDA[i].method == *method)
            return numerical_method_map_CUDA[i].solver;
    return NULL;
}
