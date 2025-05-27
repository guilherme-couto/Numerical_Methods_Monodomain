#include "logger.h"

__host__ void log_device_info(FILE *log_file)
{
    if (!log_file)
        return;

    fprintf(log_file, "================================================================\n");
    fprintf(log_file, "                         Device Info\n");
    fprintf(log_file, "----------------------------------------------------------------\n");

    int runtime_version = 0;
    cudaRuntimeGetVersion(&runtime_version);

    int device_count = 0;
    cudaGetDeviceCount(&device_count);

    fprintf(log_file, "----------------------------- CUDA ------------------------------\n");
    fprintf(log_file, " CUDA Runtime Version  : %d.%d\n", runtime_version / 1000, (runtime_version % 1000) / 10);
    fprintf(log_file, " CUDA Devices Found    : %d\n", device_count);

    cudaDeviceProp prop;
    for (int i = 0; i < device_count; ++i)
    {           
        cudaGetDeviceProperties(&prop, i);
        fprintf(log_file, "  Device %d: %s\n", i, prop.name);
        fprintf(log_file, "    Compute Capability : %d.%d\n", prop.major, prop.minor);
        fprintf(log_file, "    Global Memory (GB) : %.2f\n", prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
        fprintf(log_file, "    Multiprocessors    : %d\n", prop.multiProcessorCount);
        fprintf(log_file, "    Clock Rate (MHz)   : %.2f\n", prop.clockRate / 1000.0);
    }

    fprintf(log_file, "================================================================\n\n");

    fflush(log_file);
}
