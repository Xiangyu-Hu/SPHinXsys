#---------------------------------------
# SPHINXSYS V.1.0.0
#---------------------------------------
cmake_minimum_required(VERSION 3.10)
project(SPHINXsys VERSION 1.0.0 LANGUAGES CXX)

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
    set(CMAKE_CXX_FLAGS "-Wall -m64 -std=c++11 -DUNIX")
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

##### compliler flags for physical properties #####
option(_RIEMANN_ "Enable Riemann solvers"  ON)

if (${_RIEMANN_})
    add_definitions(-D_RIEMANN_)
endif()
###################################################

enable_testing()

include(Common)
include(Import3rdpartylibs)
add_subdirectory(SPHINXsys)
add_subdirectory(cases_test)

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
