# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.15)
project(cuik_molmaker VERSION 0.1 LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add option for building tests
option(CUIKMOLMAKER_BUILD_TESTS "Build tests for cuik_molmaker" OFF)

# Set RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

# Add compiler flags
add_compile_options(-O3 -Wall -Wmissing-field-initializers -Wmaybe-uninitialized -Wuninitialized)

# Set Python paths to match RDKit's Python dependency
set(PYTHON_ROOT_DIR $ENV{CONDA_PREFIX})
message(STATUS "PYTHON_ROOT_DIR: ${PYTHON_ROOT_DIR}")

# Set Python paths and hints before FindPython is called
set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
message(STATUS "Using Python executable: ${Python_EXECUTABLE}")

# Make sure FindPython uses our specified executable by setting Python_ROOT
set(Python_ROOT "${CONDA_PREFIX}")

# Explicitly specify Python hints to FindPython
set(CMAKE_FIND_FRAMEWORK NEVER) # Avoid issues with macOS frameworks
set(CMAKE_PREFIX_PATH "${CONDA_PREFIX};${CMAKE_PREFIX_PATH}")

# Find required packages - will use the Python_EXECUTABLE we specified
find_package(Python COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

find_package(Torch REQUIRED) 
message(STATUS "TORCH_INCLUDE_DIRS: ${TORCH_INCLUDE_DIRS}")
find_library(TORCH_PYTHON_LIBRARY torch_python PATH "${TORCH_INSTALL_PREFIX}/lib")

# Enable testing only if CUIKMOLMAKER_BUILD_TESTS is ON
if(CUIKMOLMAKER_BUILD_TESTS)
    message(STATUS "Building with tests enabled")
    enable_testing()
    
    # Include module to fetch Catch2
    include(cmake/FetchCatch2.cmake)
endif()

find_package(pybind11 REQUIRED)

# Find Boost package
find_package(Boost REQUIRED COMPONENTS system serialization iostreams)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")

# Set RDKit paths explicitly
# set(RDKit_DIR "${RDKIT_PYTHON_PATH}/lib/cmake/rdkit")
# message(STATUS "Setting RDKit_DIR to: ${RDKit_DIR}")

find_package(RDKit REQUIRED)
message(STATUS "RDKit_VERSION: ${RDKit_VERSION}")
message(STATUS "RDKit_INCLUDE_DIRS: ${RDKit_INCLUDE_DIRS}")
message(STATUS "RDKit_LIBRARIES: ${RDKit_LIBRARIES}")

set(RDKit_LIBS
    RDKit::GraphMol
    RDKit::MolStandardize
    RDKit::SmilesParse
    )

find_package(Python3 COMPONENTS NumPy REQUIRED)
message(STATUS "NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")

# Create the core library that can be used by pybind11 and tests
add_library(cuik_molmaker_core SHARED
    src/features.cpp
    src/float_features.cpp
    src/one_hot.cpp
)

# Add include directories for core library
target_include_directories(cuik_molmaker_core PUBLIC
    ${TORCH_INCLUDE_DIRS}
    ${RDKIT_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Link against required libraries for core library
target_link_libraries(cuik_molmaker_core PUBLIC
    ${RDKit_LIBS}
    ${TORCH_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${TORCH_PYTHON_LIBRARY}
    ${Boost_LIBRARIES}
)
set_target_properties(cuik_molmaker_core PROPERTIES
    INSTALL_RPATH "$ORIGIN/../../../.."
    BUILD_WITH_INSTALL_RPATH TRUE
)
# Create the Python module using the core library
pybind11_add_module(cuik_molmaker
    src/cuik_molmaker_cpp.cpp
)

# Link the Python module against the core library
target_link_libraries(cuik_molmaker PRIVATE
    cuik_molmaker_core
)

# Set RPATH for the Python module - platform specific
if(APPLE)
    # macOS: Use @loader_path to find dylib in lib/ subdirectory
    set_target_properties(cuik_molmaker PROPERTIES
        INSTALL_RPATH "@loader_path/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
else()
    # Linux: Use $ORIGIN to find .so in lib/ subdirectory
    set_target_properties(cuik_molmaker PROPERTIES
        INSTALL_RPATH "$ORIGIN/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Set installation targets
install(TARGETS cuik_molmaker
    LIBRARY DESTINATION ${Python_SITELIB}
    RUNTIME DESTINATION ${Python_SITELIB}
)

# Build test executables only if CUIKMOLMAKER_BUILD_TESTS is ON
if(CUIKMOLMAKER_BUILD_TESTS)      
    # Catch2 tests
    # Create the main executable for Catch2 tests
    add_executable(catch2_tests 
        tests/cpp/catch2/test_features.cpp
    )
    
    # Link against the core library and Catch2
    target_link_libraries(catch2_tests PRIVATE
        cuik_molmaker_core
        Catch2::Catch2WithMain
    )
    
    # Register Catch2 tests with CTest
    include(CTest)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Catch2/Catch.cmake)
    catch_discover_tests(catch2_tests)
endif()

