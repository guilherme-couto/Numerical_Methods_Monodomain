#ifndef CORE_DEFINITIONS_H
#define CORE_DEFINITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#ifdef USE_CUDA

#include <cuda_runtime.h>

#endif // USE_CUDA

#ifdef __cplusplus
extern "C" {
#endif

// Define path separator based on OS
#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path) // On Windows, _mkdir ignores mode
#define PATH_SEPARATOR '\\'
#else
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/sysinfo.h>
#define PATH_SEPARATOR '/'
#endif

// Define maximum string size
#define MAX_STRING_SIZE 200

#if defined(__CUDACC__)
#define FORCE_INLINE __forceinline__
#define STATIC_MODIFIER static __device__
#else
#define FORCE_INLINE inline
#define STATIC_MODIFIER static
#endif

// Define block size for GPU
#define FULL_DOMAIN_BLOCK_SIZE_X 16
#define FULL_DOMAIN_BLOCK_SIZE_Y 16
#define THOMAS_KERNEL_BLOCK_SIZE 32

#define MAX_SYS_SIZE 1024
#define NUMTHREADS 6

// Convert CM to UM
#define CM_TO_UM(x) ((int)(x * 1.0e4))

// Define real type
typedef double real;
#define REAL_TYPE "double"
#define STRTOREAL strtod
#define FSCANF_REAL "%le"
#define PRINTF_REAL "%lf"
#ifdef USE_FLOAT
typedef float real;
#define REAL_TYPE "float"
#define STRTOREAL strtof
#define FSCANF_REAL "%e"
#define PRINTF_REAL "%f"
#endif

#define _PI 3.14159265358979323846f

// Define measurement structure
typedef struct
{
    real elapsedExecutionTime;
    real elapsedTime1stPart;
    real elapsedTime2ndPart;
    real elapsedTime1stLS;
    real elapsedTime2ndLS;
    real elapsedSaveFramesTime;
    real elapsedMeasureVelocityTime;
    real elapsedSaveStateTime;
    real stimVelocity;
} Measurement;

// Define stimulus structure
typedef struct
{
    real amplitude;
    real start_time;
    real duration;
    struct
    {
        real min;
        real max;
    } x_range, y_range, x_discretized, y_discretized;
} Stimulus;

// Define cell phenotypes
typedef enum
{
    CELL_PHENOTYPE_EPI,
    CELL_PHENOTYPE_ENDO,
    CELL_PHENOTYPE_M,
    CELL_PHENOTYPE_INVALID
} CellPhenotype;

// Define a structure for element properties
typedef struct
{
    CellPhenotype cell_phenotype;
    real D_xx; // Diffusion coefficient in x direction
    real D_yy; // Diffusion coefficient in y direction
    real D_xy; // Diffusion coefficient in xy direction
} ElementProperties;

// Define CUDA error checking
#define CUDA_CALL(call)                                                                                                          \
    do                                                                                                                           \
    {                                                                                                                            \
        cudaError_t error = call;                                                                                                \
        if (error != cudaSuccess)                                                                                                \
        {                                                                                                                        \
            fprintf(stderr, "\n\033[31m[x] \033[0mCUDA error at %s:%d - %s\n\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            return;                                                                                                              \
        }                                                                                                                        \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif // CORE_DEFINITIONS_H