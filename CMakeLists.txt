# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.15)
project(cuik_molmaker_cpp VERSION 0.2 LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include required modules
include(CMakeDependentOption)

# Add option for building tests
option(CUIKMOLMAKER_BUILD_TESTS "Build tests for cuik_molmaker" OFF)
option(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT
       "Whether or not to build against the rdkit from pip" OFF)
option(CUIKMOLMAKER_CXX11_ABI
    "Whether or not to build with CXX11 ABI"
    OFF
)
cmake_dependent_option(
  CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR
  "Location of the rdkit library directory for pip installs" ""
  CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT "")
cmake_dependent_option(
  CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR
  "Location of the rdkit include directory for pip installs" ""
  CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT "")
cmake_dependent_option(
  CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR
  "Location of the boost include directory for pip installs" ""
  CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT "")


# Set RPATH settings
# Note: Target-specific RPATH settings are configured below for each library
# to use $ORIGIN relative paths for portability in conda environments
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# Don't set a global CMAKE_INSTALL_RPATH here - let target-specific settings handle it
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

if(CUIKMOLMAKER_CXX11_ABI)
  message(STATUS "Using cxx11 ABI")
  add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=1)
else()
  message(STATUS "Using pre-cxx11 ABI")
  add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)
endif()

# Platform-specific compiler flags
if(MSVC)
    add_compile_options(/O2 /W4)
else()
    add_compile_options(-O3 -Wall -Wmissing-field-initializers -Wmaybe-uninitialized -Wuninitialized)
endif()

# Windows-specific: Tell RDKit we're importing from DLLs
if(WIN32)
    add_compile_definitions(RDKIT_DYN_LINK)
    add_compile_definitions(BOOST_ALL_DYN_LINK)
endif()

# Set Python paths to match RDKit's Python dependency
set(PYTHON_ROOT_DIR $ENV{CONDA_PREFIX})
message(STATUS "PYTHON_ROOT_DIR: ${PYTHON_ROOT_DIR}")

# Set Python paths and hints before FindPython is called
if(WIN32)
    set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/python.exe")
else()
    set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
endif()
message(STATUS "Using Python executable: ${Python_EXECUTABLE}")

# Make sure FindPython uses our specified executable by setting Python_ROOT
set(Python_ROOT "$ENV{CONDA_PREFIX}")

# Explicitly specify Python hints to FindPython
set(CMAKE_FIND_FRAMEWORK NEVER) # Avoid issues with macOS frameworks
set(CMAKE_PREFIX_PATH "${CONDA_PREFIX};${CMAKE_PREFIX_PATH}")

# Find required packages - will use the Python_EXECUTABLE we specified
find_package(Python COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

# Enable testing only if CUIKMOLMAKER_BUILD_TESTS is ON
if(CUIKMOLMAKER_BUILD_TESTS)
    message(STATUS "Building with tests enabled")
    enable_testing()
    
    # Include module to fetch Catch2
    include(cmake/FetchCatch2.cmake)
endif()

# Use Python include directories found by FindPython
set(PYTHON_INCLUDE_DIRS "${Python_INCLUDE_DIRS}")
message(STATUS "Using PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")

find_package(pybind11 2.13 REQUIRED)

if(NOT CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)
    # Find Boost package
    find_package(Boost REQUIRED COMPONENTS serialization iostreams OPTIONAL_COMPONENTS system)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Building against system RDKit from conda")
    find_package(RDKit REQUIRED)
    set(RDKit_LIBS
    RDKit::GraphMol
    RDKit::MolStandardize
    RDKit::SmilesParse
    )

    message(STATUS "RDKit_VERSION: ${RDKit_VERSION}")
    message(STATUS "RDKit_INCLUDE_DIRS: ${RDKit_INCLUDE_DIRS}")
    message(STATUS "RDKit_LIBRARIES: ${RDKit_LIBRARIES}")
else()
    if(NOT CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR)
    message(
    FATAL_ERROR
        "CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR must be set for builds against pip install"
    )
    endif()
    if(NOT CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR)
    message(
    FATAL_ERROR
        "CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR must be set for builds against pip install"
    )
    endif()
    if(NOT CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR)
    message(
    FATAL_ERROR
        "CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR must be set for builds against pip install"
    )
    endif()

    message(STATUS "CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR: ${CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR}")
    message(STATUS "CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR: ${CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR}")
    message(STATUS "CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR: ${CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR}")

  # make a list of all files ine the libdir
  # Normalize path to use forward slashes (CMake prefers this)
  file(TO_CMAKE_PATH "${CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR}" LIBDIR_NORMALIZED)
  
  if(WIN32)
    # On Windows, look for .dll files
    file(GLOB ALL_RDKIT_FILES "${LIBDIR_NORMALIZED}/RDKit*.dll")
    # Select only certain RDKit*.dll files
    set(RDKIT_FILES "")
    foreach(file ${ALL_RDKIT_FILES})
        if(file MATCHES "RDKitGraphMol" OR file MATCHES "RDKitRDGeneral" OR file MATCHES "RDKitSmilesParse" OR file MATCHES "RDKitSubstructMatch")
            list(APPEND RDKIT_FILES "${file}")
        endif()
    endforeach()
  else()
    # On Unix, look for all library files
    file(GLOB RDKIT_FILES "${LIBDIR_NORMALIZED}/*")
  endif()
  
  # for each file, make an imported library with that lib as source
  message(
    STATUS "Searched for RDKit libs in: ${LIBDIR_NORMALIZED}")
  message(STATUS "Found RDKit pip libs: ${RDKIT_FILES}")
  # Populate RDKIT_LIBS with the imported libraries
  foreach(lib ${RDKIT_FILES})
    get_filename_component(libname ${lib} NAME_WE)
    message(STATUS "libname: ${libname}")
    add_library(${libname} SHARED IMPORTED)
    set_target_properties(${libname} PROPERTIES IMPORTED_LOCATION ${lib})
    if(WIN32)
        # In Windows, we need to set the IMPORTED_IMPLIB property to the .lib file
        if(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)
            # Get the path to the matching .lib file by replacing .dll with .lib in the filename
            string(REGEX REPLACE "\\.dll$" ".lib" pip_rdkit_lib_path "${lib}")
            # Normalize the path (for CMake compatibility)
            file(TO_CMAKE_PATH "${pip_rdkit_lib_path}" pip_rdkit_lib_path_normalized)
            message(STATUS "  pip_rdkit_lib: ${pip_rdkit_lib_path_normalized}")
            set_target_properties(${libname} PROPERTIES IMPORTED_IMPLIB "${pip_rdkit_lib_path_normalized}")

        else()
            set_target_properties(${libname} PROPERTIES IMPORTED_IMPLIB "${libname}.lib")
        endif()
    endif()
    # Set include dirs to include both CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR and
    # CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR
    target_include_directories(
      ${libname} SYSTEM INTERFACE ${CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR}
                                  ${CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR})
    list(APPEND RDKit_LIBS ${libname})
  endforeach()
  set(BOOST_INCLUDE_DIRS ${CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR})
  message(STATUS "Using boost libs from pip RDKit ${CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR}")
endif(NOT CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)


# Create the core library that can be used by pybind11 and tests
add_library(cuik_molmaker_core SHARED
    src/features.cpp
    src/float_features.cpp
    src/one_hot.cpp
)

# Add include directories for core library
target_include_directories(cuik_molmaker_core PUBLIC
    ${RDKIT_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

message(STATUS "pybind11_INCLUDE_DIRS: ${pybind11_INCLUDE_DIRS}")

# Link against required libraries for core library
message(STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

if(APPLE)
    message(STATUS "Linking order: MacOS")
    target_link_libraries(cuik_molmaker_core PUBLIC
        ${RDKit_LIBS}
        ${pybind11_LIBRARIES}
        ${PYTHON_LIBRARIES} # PYTHON_LIBRARIES needed for MacOS
        ${Boost_LIBRARIES}
    )
else()
    message(STATUS "Linking order: Linux or Windows")
    target_link_libraries(cuik_molmaker_core PUBLIC
        ${Boost_LIBRARIES}
        ${RDKit_LIBS}
        ${Python_LIBRARIES}
    )
endif()

# On macOS, allow undefined symbols for Python C API functions
# They will be resolved at runtime when loaded by the Python interpreter
if(APPLE)
    target_link_options(cuik_molmaker_core PRIVATE "-undefined" "dynamic_lookup")
endif()

if(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)
    if(UNIX)
        set_target_properties(cuik_molmaker_core PROPERTIES
            INSTALL_RPATH "$ORIGIN/../../rdkit.libs;$ORIGIN/../../../.."
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    elseif(APPLE)
        # TODO: Figure this out for MacOS
        set_target_properties(cuik_molmaker_core PROPERTIES
            INSTALL_RPATH "$ORIGIN/../../../.."
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    else()
        # On Windows, DLLs should be in the same directory or on PATH
        # No RPATH needed, but we can set output directory
        set_target_properties(cuik_molmaker_core PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    endif()
else()
    if(APPLE OR UNIX)
        set_target_properties(cuik_molmaker_core PROPERTIES
            INSTALL_RPATH "$ORIGIN/../../../.."
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    else()
        # TODO: Need to figure this out for Windows
        # On Windows, DLLs should be in the same directory or on PATH
        # No RPATH needed, but we can set output directory
        # set_target_properties(cuik_molmaker_core PROPERTIES
        #     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        # )
    endif(APPLE OR UNIX)
endif(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)

    # Create the Python module using the core library
pybind11_add_module(cuik_molmaker_cpp
    src/cuik_molmaker_cpp.cpp
)

message(STATUS "pybind11_LIBRARIES: ${pybind11_LIBRARIES}")

# Link the Python module against the core library
target_link_libraries(cuik_molmaker_cpp PRIVATE
    cuik_molmaker_core
)

# Set RPATH for the Python module - platform specific
# The compiled module links the RDKit/Boost libraries directly, so its own RPATH
# must point at them. When building against pip RDKit those libraries live in the
# sibling "rdkit.libs" directory (site-packages/rdkit.libs), one level above the
# module (site-packages/cuik_molmaker). Modern toolchains emit DT_RUNPATH, which is
# NOT used to resolve transitive dependencies, so the runpath on the core library
# alone is not sufficient and we must add rdkit.libs here as well.
if(APPLE)
    # macOS: Use @loader_path to find dylib in lib/ subdirectory
    set(rpath_entries "@loader_path/lib")
    if(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)
        list(APPEND rpath_entries "@loader_path/../rdkit.libs")
    endif()

    set_target_properties(cuik_molmaker_cpp PROPERTIES
        INSTALL_RPATH "${rpath_entries}"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
# elseif(WIN32)
    # On Windows, do not set INSTALL_RPATH/BUILD_WITH_INSTALL_RPATH.
    # Instead, ensure that the DLLs are either (1) in the same directory as the .pyd,
    # or (2) on PATH at runtime. This is necessary so cuik_molmaker_cpp can find its dependencies.
    # Optionally, copy dependent DLLs to the output directory if needed.
    # No RPATH equivalent for Windows.
else()
    # Linux: Use $ORIGIN to find .so in lib/ subdirectory
    set(rpath_entries "$ORIGIN/lib")
    if(CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT)
        list(APPEND rpath_entries "$ORIGIN/../rdkit.libs")
    endif()
    list(APPEND rpath_entries "$ORIGIN/../../../")
    set_target_properties(cuik_molmaker_cpp PROPERTIES
        INSTALL_RPATH "${rpath_entries}"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Set installation targets
install(TARGETS cuik_molmaker_cpp
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
        ${Python_LIBRARIES}  # Add this if not already present
        pybind11::embed
    )
    
    # Register Catch2 tests with CTest
    include(CTest)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Catch2/Catch.cmake)
    catch_discover_tests(catch2_tests)
endif()