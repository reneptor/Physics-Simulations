cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(ISING)


set(CMAKE_BUILD_TYPE Debug)

# Set the source directory
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# Add the executable
add_executable(ISING 
    "${SOURCE_DIR}/main.cpp"
    "${SOURCE_DIR}/ising_lattice.cpp"
    "${SOURCE_DIR}/metropolis.cpp"
    "${SOURCE_DIR}/wolff.cpp"
    "${SOURCE_DIR}/swedsen_wang.cpp"
    "${SOURCE_DIR}/store_data.cpp"  
    # "${SOURCE_DIR}/run_multiple.cpp"   Fix this later if necessary.
    "${PROJECT_SOURCE_DIR}/include/ising_model.h"
)

# Include the Eigen headers (not necessary if `find_package()` works)
target_include_directories(ISING PRIVATE "${PROJECT_SOURCE_DIR}/include")

