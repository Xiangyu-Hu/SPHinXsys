#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.10)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

set(INCLUDE_DIRS)
# append
set(INCLUDE_DIRS
  ${PLATFORM_INCLUDE_DIRECTORIES}
  ${PROJECT_SOURCE_DIR}/SPHINXsys/src/shared/simbody-include
  ${PROJECT_SOURCE_DIR}/SPHINXsys/src/shared/oneTBB/include)

include_directories(INCLUDE_DIRS)

set(SIMBODY_MAJOR_VERSION 3)
set(SIMBODY_MINOR_VERSION 7)
set(SIMBODY_PATCH_VERSION 0)

# Report the version number to the CMake UI. Don't include the
# build version if it is zero.
set(PATCH_VERSION_STRING)
if(SIMBODY_PATCH_VERSION)
    set(PATCH_VERSION_STRING ".${SIMBODY_PATCH_VERSION}")
endif()

set(SIMBODY_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}${PATCH_VERSION_STRING}"
    CACHE STRING
    "This is the version that will be built (can't be changed in GUI)."
    FORCE)

set(SIMBODY_SONAME_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}"
    CACHE STRING
    "Soname version; appended to names of shared libs
    (can't be changed in GUI)."
    FORCE)

# This is the suffix if we're building and depending on versioned libraries.
set(VN "_${SIMBODY_VERSION}")

set(SimTKCOMMON_MAJOR_VERSION ${SIMBODY_MAJOR_VERSION})
set(SimTKCOMMON_MINOR_VERSION ${SIMBODY_MINOR_VERSION})
set(SimTKCOMMON_PATCH_VERSION ${SIMBODY_PATCH_VERSION})

set(SIMMATH_MAJOR_VERSION ${SIMBODY_MAJOR_VERSION})
set(SIMMATH_MINOR_VERSION ${SIMBODY_MINOR_VERSION})
set(SIMMATH_PATCH_VERSION ${SIMBODY_PATCH_VERSION})

add_definitions(-DSimTK_SIMMATH_MAJOR_VERSION=${SIMMATH_MAJOR_VERSION}
                -DSimTK_SIMMATH_MINOR_VERSION=${SIMMATH_MINOR_VERSION}
        -DSimTK_SIMMATH_PATCH_VERSION=${SIMMATH_PATCH_VERSION})

add_definitions(-DSimTK_SIMBODY_MAJOR_VERSION=${SIMMATH_MAJOR_VERSION}
                -DSimTK_SIMBODY_MINOR_VERSION=${SIMMATH_MINOR_VERSION}
        -DSimTK_SIMBODY_PATCH_VERSION=${SIMMATH_PATCH_VERSION})   
        
add_definitions(-DBUILD_VISUALIZER=off)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

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
else(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -m64 -DUNIX")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -ggdb")
    IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")
    ELSE(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
endif(MSVC)

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
add_subdirectory(SPHINXsys)
add_subdirectory(cases_test)

if(ACTIVATE_CUDA)
    add_subdirectory(cases_test_gpu)
endif()

add_subdirectory(cases_user)

set(BUILD_WITH_BENCHMARK_CASES OFF CACHE BOOL "Include the benchmark cases into the project?")
if(BUILD_WITH_BENCHMARK_CASES)
add_subdirectory(cases_benchmark)
endif()

set(BUILD_WITH_XIANGYU_HU_CASES OFF CACHE BOOL "Include Xiangyu Hu's cases into the project?")
if(BUILD_WITH_XIANGYU_HU_CASES)
add_subdirectory(cases_xyh)
endif()

set(BUILD_WITH_CHI_ZHANG_CASES OFF CACHE BOOL "Include Chi Zhang's cases into the project?")
if(BUILD_WITH_CHI_ZHANG_CASES)
add_subdirectory(cases_cz)
endif()

set(BUILD_WITH_SHUOGUO_ZHANG_CASES OFF CACHE BOOL "Include Shuoguo Zhang's cases into the project?")
if(BUILD_WITH_SHUOGUO_ZHANG_CASES)
add_subdirectory(cases_sgz)
endif()

set(BUILD_WITH_MASSOUD_REZAVAND_CASES OFF CACHE BOOL "Include Massoud Rezavand's cases into the project?")
if(BUILD_WITH_MASSOUD_REZAVAND_CASES)
add_subdirectory(cases_mr)
endif()

set(BUILD_WITH_YONGCHUAN_YU_CASES OFF CACHE BOOL "Include Yongchuan Yu's cases into the project?")
if(BUILD_WITH_YONGCHUAN_YU_CASES)
add_subdirectory(cases_ycy)
endif()

set(BUILD_WITH_YUJIE_ZHU_CASES OFF CACHE BOOL "Include Yujie Zhu's cases into the project?")
if(BUILD_WITH_YUJIE_ZHU_CASES)
add_subdirectory(cases_yjz)
endif()

set(BUILD_WITH_WENBIN_ZHANG_CASES OFF CACHE BOOL "Include Wenbin Zhang's cases into the project?")
if(BUILD_WITH_WENBIN_ZHANG_CASES)
add_subdirectory(cases_wbz)
endif()

set(BUILD_WITH_DONG_WU_CASES OFF CACHE BOOL "Include Dong Wu's cases into the project?")
if(BUILD_WITH_DONG_WU_CASES)
add_subdirectory(cases_dw)
endif()

set(BUILD_WITH_BR_CASES OFF CACHE BOOL "Include BR's cases into the project?")
if(BUILD_WITH_BR_CASES)
add_subdirectory(cases_br)
endif()

set(BUILD_WITH_NA_CASES OFF CACHE BOOL "Include Niklas Aufenanger's cases into the project?")
if(BUILD_WITH_NA_CASES)
add_subdirectory(cases_na)
endif()

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.

if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(PLATFORM_ABI x64)
else()
    set(PLATFORM_ABI x86)
endif()

set(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)
    
# Make the pkgconfig file: select the proper flags depending on the platform
if(WIN32)
    if(${PLATFORM_ABI} MATCHES "x64")
        # win 64 bits
        set(PKGCONFIG_PLATFORM_LIBS "-lliblapack -llibblas")
    else()
        set(PKGCONFIG_PLATFORM_LIBS "-lliblapack -llibblas")
    endif()
elseif(APPLE)
    set(PKGCONFIG_PLATFORM_LIBS     "-llapack -lblas -ldl")
else()
    set(PKGCONFIG_PLATFORM_LIBS     "-llapack -lblas -lpthread -lrt -ldl -lm")
endif()