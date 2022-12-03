#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.16)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

set(SPHINXSYS_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake)
include(Common) # brings macro to current namespace

###### Change settings in this file if you want to build with Simbody and/or oneTBB source code ######
option(SPH_ONLY_STATIC_BUILD "SPH_ONLY_STATIC_BUILD" 0)
option(ONLY_3D "ONLY_3D" 0)
option(BUILD_WITH_SIMBODY "BUILD_WITH_SIMBODY" 0)
option(EMSCRIPTEN "EMSCRIPTEN" 0)
option(BUILD_SPHINXSYS_TESTS "BUILD_SPHINXSYS_TESTS" 1)


set(CMAKE_VERBOSE_MAKEFILE ON) # Enable verbose makefile for whoever use them
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate `compile_commands.json` file which is necessary for code completioner, static analyzer, etc

set(CXX_STANDARD_REQUIRED ON) # C++ standard requirement for targets is NOT optional 
set(CXX_EXTENSIONS OFF) # Disable any non-standard compiler specific C++ extensions.

add_library(sphinxsys_core INTERFACE)
target_compile_features(sphinxsys_core INTERFACE cxx_std_17)

FINDLIBHOME(SIMBODY_HOME)
find_package(SIMBODY REQUIRED)
target_include_directories(sphinxsys_core INTERFACE ${Simbody_INCLUDE_DIR})
target_link_directories(sphinxsys_core INTERFACE ${Simbody_LIB_DIR})
target_link_libraries(sphinxsys_core INTERFACE ${Simbody_LIBRARIES})

find_package(Eigen3 3.4 CONFIG QUIET)
if(TARGET Eigen3::Eigen)
    target_link_libraries (sphinxsys_core INTERFACE Eigen3::Eigen)
else()
    find_path(EIGEN_INCLUDE_DIR Eigen/Eigen HINTS ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/eigen)
    target_include_directories(sphinxsys_core INTERFACE ${EIGEN_INCLUDE_DIR})
endif()

find_package(TBB REQUIRED)
target_include_directories(sphinxsys_core INTERFACE ${TBB_INCLUDE_DIR})
target_link_directories(sphinxsys_core INTERFACE ${TBB_LIB_DIR})
target_link_libraries(sphinxsys_core INTERFACE ${TBB_LIBRARIES})
target_compile_definitions(sphinxsys_core INTERFACE TBB_SUPPRESS_DEPRECATED_MESSAGES)


if(MSVC)
    target_compile_definitions(sphinxsys_core INTERFACE _USE_MATH_DEFINES)
    target_compile_options(sphinxsys_core INTERFACE /MP)
    target_compile_options(sphinxsys_core INTERFACE /permissive-) # /W4 /WX 
endif()
if(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(sphinxsys_core INTERFACE Threads::Threads)
    target_compile_options(sphinxsys_core INTERFACE -Werror)
    # target_compile_options(sphinxsys_core INTERFACE -pedantic -Wall -Wextra) # For when one gets rid of Simbody 
endif()

if(EMSCRIPTEN)
    set(SPH_ONLY_STATIC_BUILD 1)
    set(BUILD_WITH_SIMBODY 1)
    set(ONLY_3D 1)
    set(BUILD_SPHINXSYS_TESTS 1)
    target_compile_definitions(sphinxsys_core INTERFACE -D__SIMBODY_WITHOUT_LAPACK__)
    target_compile_definitions(sphinxsys_core INTERFACE -D__EIGEN__)
    target_compile_definitions(sphinxsys_core INTERFACE -D__EMSCRIPTEN__)

    remove_definitions(-DTBB_2021_2_0)
    set(TBB_ROOT ${SPHINXSYS_PROJECT_DIR}/3rd_party/wasmtbb)
    target_include_directories(sphinxsys_core INTERFACE "${TBB_ROOT}/include")
    set(OBJDIR "${TBB_ROOT}/build/linux_wasm32_emscripten_emscripten_wasm32_release")
    file(GLOB TBB_OBJECTS ${OBJDIR}/*.o)
else()
    target_compile_definitions(sphinxsys_core INTERFACE BOOST_AVAILABLE)
    find_package(Boost REQUIRED COMPONENTS program_options filesystem system)
    target_include_directories(sphinxsys_core INTERFACE ${Boost_INCLUDE_DIRS})
    target_link_directories(sphinxsys_core INTERFACE ${Boost_LIBRARY_DIRS})
    target_link_libraries(sphinxsys_core INTERFACE ${Boost_LIBRARIES})
endif()

if(BUILD_WITH_SIMBODY)
    target_compile_definitions(sphinxsys_core INTERFACE -D__SIMBODY_WITHOUT_LAPACK__)
    target_compile_definitions(sphinxsys_core INTERFACE -D__EIGEN__)
endif()

if(SPH_ONLY_STATIC_BUILD)
    set(BUILD_STATIC_LIBRARIES 1) # Simbody
    set(BUILD_DYNAMIC_LIBRARIES 0) # Simbody
else()
    set(BUILD_SHARED_LIBS 1)
endif()

if(BUILD_WITH_SIMBODY)
    add_subdirectory(${SPHINXSYS_PROJECT_DIR}/3rd_party/simbody)
endif()

add_subdirectory(SPHINXsys)

if(BUILD_SPHINXSYS_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

add_subdirectory(modules)

file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/PythonScriptStore/RegressionTest/regression_test_base_tool.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PythonScriptStore/RegressionTest)
        