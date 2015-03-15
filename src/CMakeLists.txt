# by default the final executable is next to the source file, to clearly
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wfatal-errors ${OpenMP_CXX_FLAGS}")

list(APPEND CUDA_NVCC_FLAGS -O3 -gencode arch=compute_52,code=sm_52 -lineinfo)

INCLUDE_DIRECTORIES(
    # ${CMAKE_SOURCE_DIR}/lib/local/openblas/include
    # ${CMAKE_SOURCE_DIR}/lib/local/magma/include
    ${CMAKE_SOURCE_DIR}/lib
    ${LIBOPENBLAS_INCLUDE_DIR}
    ${LIBMAGMA_INCLUDE_DIR}
)

# defines all source code files
set(SOURCES
    parameters.cc
    performance.cc
    resources.cc
    main.cu
    simulation.cu
)

# defines all header code files
set(HEADERS
    common.h
    fiberopt.h
    parameters.h
    resources.h
    simulation.h
)

set(KERNELS
    kernels/assemble_system_1D.cu
    kernels/assemble_system_2D.cu
    kernels/assemble_system_3D.cu
    kernels/compute_inner_integral_analytically.cu
    kernels/compute_inner_integral_numerically.cu
    kernels/constants.cu
    kernels/reset_system.cu
    kernels/update_fibers.cu
    kernels/update_fibers_firststep.cu
    kernels/update_velocities_1D.cu
    kernels/update_velocities_2D.cu
    kernels/reset_velocities.cu
)

set(ADDONS
    ../fibers_runner.py
    ../tools/validate.py
    ../tools/validate_mapping.py
    ../tools/gen2.py
)


# add an executable/binary to the project
CUDA_ADD_EXECUTABLE(fibers ${SOURCES})

# it is not necessary to include the headers here for cmake, however this way
# they also appear in the file browser of different IDEs like Xcode/VisualStudio
add_custom_target(fibers_files DEPENDS ${SOURCES} SOURCES ${HEADERS} ${KERNELS} ${ADDONS})

TARGET_LINK_LIBRARIES(fibers ${LIBMAGMA_LIBRARY} ${LIBOPENBLAS_LIBRARY})

# it is important that this comes AFTER target_link_libraries
# no idea why
CUDA_ADD_CUBLAS_TO_TARGET(fibers)
