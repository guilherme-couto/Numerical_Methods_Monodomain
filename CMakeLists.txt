cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0104 NEW)

set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_C_COMPILER} CACHE FILEPATH "Host compiler used by NVCC")

project(MonodomainSimulation C CXX CUDA)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# === Set Release build and flags ===
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops -fno-trapping-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops -fno-trapping-math")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 --use_fast_math --fmad=true -Xcompiler \"-O3 -march=native -ffast-math -funroll-loops -fno-trapping-math\" ")

# === Build options ===
option(USE_CUDA "Enable CUDA support" ON)
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_FLOAT "Use float precision instead of double" OFF)

# === Directories ===
set(INC_DIR "${CMAKE_SOURCE_DIR}/include")
set(EXT_DIR "${CMAKE_SOURCE_DIR}/external")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin")
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")

set(CELL_MODELS_DIR "${SRC_DIR}/cell_models")
set(NUMERICAL_METHODS_DIR "${SRC_DIR}/numerical_methods")
set(LOGGER_DIR "${SRC_DIR}/logger")
set(SOLVERS_DIR "${SRC_DIR}/solvers")


# === Include paths ===
include_directories(${INC_DIR} ${EXT_DIR} ${CELL_MODELS_DIR} ${NUMERICAL_METHODS_DIR} ${LOGGER_DIR} ${SOLVERS_DIR})
file(GLOB_RECURSE EXT_HEADERS "${EXT_DIR}/*.h")
file(GLOB_RECURSE CELL_HEADERS "${CELL_MODELS_DIR}/*.h")
file(GLOB_RECURSE NUMERICAL_METHODS_HEADERS "${NUMERICAL_METHODS_DIR}/*.h")
file(GLOB_RECURSE LOGGER_HEADERS "${LOGGER_DIR}/*.h")
file(GLOB_RECURSE SOLVERS_HEADERS "${SOLVERS_DIR}/*.h")

# === Source files ===
file(GLOB SRC_C_FILES "${SRC_DIR}/*.c")
file(GLOB_RECURSE EXT_C_FILES "${EXT_DIR}/*.c")
file(GLOB_RECURSE CELL_C_FILES "${CELL_MODELS_DIR}/*.c")
file(GLOB_RECURSE NUMERICAL_METHODS_C_FILES "${NUMERICAL_METHODS_DIR}/*.c")
file(GLOB_RECURSE LOGGER_C_FILES "${LOGGER_DIR}/*.c")
file(GLOB_RECURSE SOLVERS_C_FILES "${SOLVERS_DIR}/*.c")

set (C_SOURCES
    ${SRC_C_FILES}
    ${EXT_C_FILES}
    ${CELL_C_FILES}
    ${NUMERICAL_METHODS_C_FILES}
    ${LOGGER_C_FILES}
    ${SOLVERS_C_FILES}
)

# === Separate libraries ===
add_library(core_c STATIC ${C_SOURCES})
target_include_directories(core_c PUBLIC ${INC_DIR} ${EXT_DIR})

# === CUDA source files ===
if (USE_CUDA)
    file(GLOB_RECURSE SRC_CU_FILES "${SRC_DIR}/*.cu")
    file(GLOB_RECURSE CELL_CU_FILES "${CELL_MODELS_DIR}/*.cu*")
    file(GLOB_RECURSE NUMERICAL_METHODS_CU_FILES "${NUMERICAL_METHODS_DIR}/*.cu")
    file(GLOB_RECURSE LOGGER_CU_FILES "${LOGGER_DIR}/*.cu")
    file(GLOB_RECURSE SOLVERS_CU_FILES "${SOLVERS_DIR}/*.cu")
    
    set(CUDA_SOURCES ${SRC_CU_FILES} ${CELL_CU_FILES} ${NUMERICAL_METHODS_CU_FILES} ${LOGGER_CU_FILES} ${SOLVERS_CU_FILES})
    
    add_library(core_cuda STATIC ${CUDA_SOURCES})
    target_include_directories(core_cuda PUBLIC ${INC_DIR} ${EXT_DIR})
    target_compile_definitions(core_cuda PRIVATE USE_CUDA)

    set_target_properties(core_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )
endif()

# === Create executable ===
add_executable(monodomain_simulation ${SRC_DIR}/main.c)

# === System libraries ===
target_link_libraries(monodomain_simulation PRIVATE core_c pthread m)

# === Output binary ===
set_target_properties(monodomain_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# === CUDA support ===
if(USE_CUDA)
    enable_language(CUDA)
    message(STATUS "CUDA compiler: ${CMAKE_CUDA_COMPILER}")

    set_target_properties(monodomain_simulation PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    target_compile_definitions(monodomain_simulation PRIVATE USE_CUDA)

    target_link_libraries(monodomain_simulation PRIVATE core_cuda ${CUDA_LIBRARIES} cusparse)
    
    message(STATUS "✅ CUDA support ENABLED")

else()
    message(STATUS "❌ CUDA support DISABLED")
endif()

# === OpenMP support ===
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        target_link_libraries(monodomain_simulation PRIVATE OpenMP::OpenMP_C)
        target_compile_definitions(monodomain_simulation PRIVATE USE_OPENMP)
        target_compile_definitions(core_c PRIVATE USE_OPENMP)
        message(STATUS "✅ OpenMP support ENABLED")
    endif()
else()
    message(STATUS "❌ OpenMP support DISABLED")
endif()

# === Precision ===
if(USE_FLOAT)
    target_compile_definitions(monodomain_simulation PRIVATE USE_FLOAT)
    message(STATUS "🧮 Precision set to: FLOAT")
else()
    target_compile_definitions(monodomain_simulation PRIVATE USE_DOUBLE)
    message(STATUS "🧮 Precision set to: DOUBLE")
endif()

# === CUDA architecture (logging only) ===
if(USE_CUDA AND DEFINED CMAKE_CUDA_ARCHITECTURES)
    message(STATUS "🧠 CUDA architecture: ${CMAKE_CUDA_ARCHITECTURES}")
endif()
