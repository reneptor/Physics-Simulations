cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(ICING)


set(CMAKE_BUILD_TYPE Debug)

# Set the source directory
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# Add the executable
add_executable(ICING 
    "${SOURCE_DIR}/main.cpp"
    "${SOURCE_DIR}/icing_lattice.cpp"
    "${SOURCE_DIR}/store_data.cpp"
    "${SOURCE_DIR}/run_multiple.cpp"
    "${PROJECT_SOURCE_DIR}/include/icing_model.h"
)

# Include the Eigen headers (not necessary if `find_package()` works)
target_include_directories(ICING PRIVATE "${PROJECT_SOURCE_DIR}/include")


# Link the Eigen3::Eigen target to your executable
# target_link_libraries(COOMATRIX PRIVATE Eigen3::Eigen)
