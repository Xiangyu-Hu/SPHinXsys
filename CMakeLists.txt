#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.10)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

# webassembly build #
set(WASM_BUILD 0)

# oneTBB #
# "off": oneTBB not built by this project
# "static": oneTBB static build
# "dynamic": oneTBB dynamic build
# static oneTBB build does not work currently
set(oneTBB_BUILD "dynamic")
if(oneTBB_BUILD MATCHES "static")
    set(BUILD_SHARED_LIBS OFF)
endif(oneTBB_BUILD MATCHES "static")

## current oneTBB version: v2021.2.0
## set oneTBB testing off in oneTBB CmakeLists

###### BOOST ######
## turn Boost on here
#add_definitions(-DBOOST_AVAILABLE)
###### BOOST ######

###### Simbody ######
set(ONLY_3D 1)
include_directories(${PLATFORM_INCLUDE_DIRECTORIES})
set(SIMBODY_MAJOR_VERSION 3)
set(SIMBODY_MINOR_VERSION 7)
set(SIMBODY_PATCH_VERSION 0)
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
add_definitions(-DSimTK_SimTKCOMMON_MAJOR_VERSION=${SimTKCOMMON_MAJOR_VERSION}
        -DSimTK_SimTKCOMMON_MINOR_VERSION=${SimTKCOMMON_MINOR_VERSION}
        -DSimTK_SimTKCOMMON_PATCH_VERSION=${SimTKCOMMON_PATCH_VERSION})

add_definitions(-DBUILD_VISUALIZER=off)
set(CMAKE_C_FLAGS "-DINTEGER_STAR_8")
###### Simbody ######

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
    add_compile_definitions(__linux__=1)
    #set(CMAKE_CXX_FLAGS "-Wall -m64 -std=c++11 -DUNIX -pthread")
    set(CMAKE_CXX_FLAGS "-Wall -std=c++11 -DUNIX -pthread")
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

if(NOT WASM_BUILD)
    include(Common)
endif()

include(Import3rdpartylibs) # TBB and Boost

add_subdirectory(SPHINXsys)
if(oneTBB_BUILD MATCHES "static" OR oneTBB_BUILD MATCHES "dynamic")
    add_subdirectory(oneTBB)
endif(oneTBB_BUILD MATCHES "static" OR oneTBB_BUILD MATCHES "dynamic")

if(NOT WASM_BUILD)
    add_subdirectory(cases_test)
endif()

if(ACTIVATE_CUDA)
    add_subdirectory(cases_test_gpu)
endif()

if(NOT WASM_BUILD)
    add_subdirectory(cases_user)
endif()

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
