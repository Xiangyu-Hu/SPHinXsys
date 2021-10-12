#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.16)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

set(SPHINXSYS_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake)

###### Change settings in this file if you want to build with Simbody and/or oneTBB source code ######
option(STATIC_BUILD "STATIC_BUILD" 0)
option(ONLY_3D "ONLY_3D" 0)
option(BUILD_WITH_SIMBODY "BUILD_WITH_SIMBODY" 0)
option(EMSCRIPTEN "EMSCRIPTEN" 0)
option(BUILD_SPHINXSYS_TESTS "BUILD_SPHINXSYS_TESTS" 1)

if(EMSCRIPTEN)
    set(STATIC_BUILD 1)
    set(BUILD_WITH_SIMBODY 1)
    set(ONLY_3D 1)
    set(BUILD_SPHINXSYS_TESTS 1)
    add_definitions(-D__SIMBODY_WITHOUT_LAPACK__)
    add_definitions(-D__EIGEN__)
    add_definitions(-D__EMSCRIPTEN__)

    remove_definitions(-DTBB_2021_2_0)
    set(TBB_ROOT ${SPHINXSYS_PROJECT_DIR}/3rd_party/wasmtbb)
    include_directories("${TBB_ROOT}/include")
    set(OBJDIR "${TBB_ROOT}/build/linux_wasm32_emscripten_emscripten_wasm32_release")
    file(GLOB TBB_OBJECTS ${OBJDIR}/*.o)
else()
    set(BOOST_AVAILABLE 1)
    add_definitions(-DBOOST_AVAILABLE)
endif()

if(BUILD_WITH_SIMBODY)
    add_definitions(-D__SIMBODY_WITHOUT_LAPACK__)
    add_definitions(-D__EIGEN__)
    set(EIGEN_ROOT ${SPHINXSYS_PROJECT_DIR}/3rd_party/eigen)
    message("${EIGEN_ROOT}")
    include_directories("${EIGEN_ROOT}")
endif()

if(STATIC_BUILD)
    set(BUILD_STATIC_LIBRARIES 1) # Simbody
    set(BUILD_DYNAMIC_LIBRARIES 0) # Simbody
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
            "MinSizeRel" "RelWithDebInfo")
endif()

set(ACTIVATE_CUDA OFF CACHE BOOL "Activate CUDA?")
if(ACTIVATE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
      enable_language(CUDA)
    else()
      message(STATUS "No CUDA support!")
    endif()

    set(BUILD_GPU_ARCH "75" CACHE STRING "Define the CUDA architecture")
    #option(BUILD_GPU_FORCE_INLINES "CUDA 7.5 AND Ubuntu 16.04 needs this enabled" OFF)
    #string(APPEND CMAKE_CUDA_FLAGS " -arch=sm_${BUILD_GPU_ARCH} UAMMD_DEBUG")
    string(APPEND CMAKE_CUDA_FLAGS " -arch=sm_${BUILD_GPU_ARCH}")

endif()

if(MSVC)
    SET(CMAKE_CL_64 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" " /MP")
    add_definitions(-D_USE_MATH_DEFINES)
elseif(EMSCRIPTEN)
    # set(PRE_POST "")
    # #set(EM_FLAGS "${PRE_POST} -fPIC -g -O0 --profiling -DUNIX -D__linux__=1 --bind -std=c++11 -s USE_PTHREADS=1 -s DISABLE_EXCEPTION_CATCHING=0")
    # set(EM_FLAGS "${PRE_POST} -fPIC -g -O3 -DUNIX -D__linux__=1 --bind -std=c++11 -s USE_PTHREADS=1 -s DISABLE_EXCEPTION_CATCHING=0")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EM_FLAGS}")
    # set(EM_FLAGS " -s INITIAL_MEMORY=536870912 -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=0 -s ERROR_ON_UNDEFINED_SYMBOLS=1")
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EM_FLAGS}")
else()
    add_compile_definitions(__linux__=1)
    set(CMAKE_CXX_FLAGS "-w -std=c++11 -DUNIX -pthread")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -ggdb")
    IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ELSE(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
endif()

if(NOT BUILD_WITH_SIMBODY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

add_compile_definitions(TBB_SUPPRESS_DEPRECATED_MESSAGES=1)

message("===========================================")
message("${PROJECT_NAME}, VERSION ${PROJECT_VERSION}")
message("Cmake build type: ${CMAKE_BUILD_TYPE}")
message("CXX flags: ${CMAKE_CXX_FLAGS}")
message("CUDA flags: ${CMAKE_CUDA_FLAGS}")
message("===========================================")

########### GPU - compliler flags for physical properties #############
# 1. Turn ON if the body force is time dependent, e.g. sloshing tank 
option(_TIMEDEPENDENT_BODYFORCE_ "Enable time dependent body force" OFF)
if (${_TIMEDEPENDENT_BODYFORCE_})
    add_definitions(-D_TIMEDEPENDENT_BODYFORCE_)
endif()
#######################################################################

enable_testing()

include(Common)
include(Import3rdpartylibs)

if(BUILD_WITH_SIMBODY)
    set(simbody_root "${SPHINXSYS_PROJECT_DIR}/3rd_party/simbody")
    add_subdirectory(${simbody_root})
endif()

add_subdirectory(SPHINXsys)

if(BUILD_SPHINXSYS_TESTS)
    add_subdirectory(tests)
endif()

if(ACTIVATE_CUDA)
    add_subdirectory(cases_test_gpu)
endif()

add_subdirectory(modules)
