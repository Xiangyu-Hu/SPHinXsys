#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.16)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

set(SPHINXSYS_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake)
include(Common) # brings macro to current namespace

###### Change settings in this file if you want to build with Simbody and/or oneTBB source code ######
option(SPHINXSYS_2D "Build sphinxsys_2d library" ON)
option(SPHINXSYS_3D "Build sphinxsys_3d library" ON)
option(SPHINXSYS_BUILD_TESTS "Build tests" ON)

set(CMAKE_VERBOSE_MAKEFILE ON) # Enable verbose makefile for whoever use them
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate `compile_commands.json` file which is necessary for code completioner, static analyzer, etc

set(CXX_STANDARD_REQUIRED ON) # C++ standard requirement for targets is NOT optional 
set(CXX_EXTENSIONS OFF) # Disable any non-standard compiler specific C++ extensions.

add_library(sphinxsys_core INTERFACE)
target_compile_features(sphinxsys_core INTERFACE cxx_std_17)

find_package(Simbody CONFIG REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE SimTKmath SimTKcommon SimTKsimbody)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries (sphinxsys_core INTERFACE Eigen3::Eigen)

find_package(TBB CONFIG REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
target_compile_definitions(sphinxsys_core INTERFACE TBB_SUPPRESS_DEPRECATED_MESSAGES)

if(${CMAKE_SYSTEM_NAME} MATCHES Windows)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
#    set(BUILD_SHARED_LIBS OFF)
    target_compile_definitions(sphinxsys_core INTERFACE _USE_MATH_DEFINES)
    target_compile_options(sphinxsys_core INTERFACE /MP)
    target_compile_options(sphinxsys_core INTERFACE /permissive-) # /W4 /WX 
else()
    find_package(Threads REQUIRED)
    target_link_libraries(sphinxsys_core INTERFACE Threads::Threads)
    target_compile_options(sphinxsys_core INTERFACE -Werror)
    # target_compile_options(sphinxsys_core INTERFACE -pedantic -Wall -Wextra) # For when one gets rid of Simbody 
endif()

    target_compile_definitions(sphinxsys_core INTERFACE BOOST_AVAILABLE)
    find_package(Boost REQUIRED COMPONENTS program_options filesystem system)
    target_include_directories(sphinxsys_core INTERFACE ${Boost_INCLUDE_DIRS})
    target_link_directories(sphinxsys_core INTERFACE ${Boost_LIBRARY_DIRS})
    target_link_libraries(sphinxsys_core INTERFACE ${Boost_LIBRARIES})


add_subdirectory(SPHINXsys)
add_subdirectory(modules)

if(SPHINXSYS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/PythonScriptStore/RegressionTest/regression_test_base_tool.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PythonScriptStore/RegressionTest)
        