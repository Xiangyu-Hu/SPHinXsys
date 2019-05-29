#-------------------------------------------------------------------------------
# SPHINXSYS
#
# This is the master CMake file that coordinates
# the build of SPHINXSYS. There are four steps:
#
#-------------------------------------------------------------------------------
if(MSVC)
    MESSAGE("Windows and Visual Studio is used!")
else(MSVC)
    if(NOT APPLE)
        MESSAGE(" Linux system is used!")
    else(NOT APPLE)
        MESSAGE("MAC system is used!")
    endif(NOT APPLE)
endif(MSVC)

if(MSVC)
    set(SPHINXSYS_REQUIRED_CMAKE_VERSION 3.1.3)
else(MSVC)
    set(SPHINXSYS_REQUIRED_CMAKE_VERSION 2.8.10)
endif(MSVC)

cmake_minimum_required(VERSION ${SPHINXSYS_REQUIRED_CMAKE_VERSION})


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)


# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
            "MinSizeRel" "RelWithDebInfo")
endif()
#MESSAGE("${CMAKE_BUILD_TYPE}")

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0005 NEW)
endif()
# RPATH
# -----

set(SPHINXSYS_USE_INSTALL_RPATH FALSE)
if(APPLE AND NOT (${CMAKE_VERSION} VERSION_LESS 2.8.12))
    # CMake 2.8.12 introduced the ability to set RPATH for shared libraries on
    # OSX. This helps executables find the libraries they depend on without
    # having to set the DYLD_LIBRARY_PATH environment variable.

    # Consider a library libfoo.dyld and an executable bar.
    # On OSX, libraries have an "install name" that, when linking, is copied
    # into the target (e.g., bar). The "install name" can be the full path to
    # libfoo.dylib, in which case bar will have no trouble finding libfoo.dylib
    # at runtime (since it has the full path to it). This doesn't work if you
    # want to be able to relocate your project. Therefore, it's possible to use
    # special tokens in the "install name" that are evaluated separately. The
    # token "@exectuable_path" is evaluated (at run time) to the full path of
    # the exectuable (e.g., bar) that is trying to load libfoo.dylib. An even
    # more flexible token is "@rpath", which is evaluated to a path (called
    # RPATH) that can be baked into the executable just after compiling or any
    # time before running the executable (using the executable
    # "install_name_tool"). The RPATH stored in executables can also contain
    # "@executable_path", etc.

    # Set the install name of libraries to contain "@rpath".
    # This allows clients of our libraries to point to them however they wish.
    set(CMAKE_MACOSX_RPATH ON)

    # We only need to set RPATH in executables if the libraries are installed
    # into directories that are not already searched by the linker.
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
    # CMake sets isSystemDir to -1 if the lib dir is NOT a system lib dir.
    if("${isSystemDir}" STREQUAL "-1")
        # This variable is used later on to toggle if RPATH should be set for
        # specific targets.
        set(SPHINXSYS_USE_INSTALL_RPATH TRUE)
    endif()
endif()


if(MSVC)
    SET(CMAKE_CL_64 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" " /MP")
else(MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -m64 -std=c++11 -DUNIX")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -ggdb")
    IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")
    ELSE(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
endif(MSVC)
## When building in any of the Release modes, tell gcc/clang to use
## not-quite most agressive optimization.  Here we
## are specifying *all* of the Release flags, overriding CMake's defaults.
## Watch out for optimizer bugs in particular gcc versions!

if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR
   ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

    # C++11
    # -----
    # Using C++11 on OSX requires using libc++ instead of libstd++.
    # libc++ is an implementation of the C++ standard library for OSX.
    if(APPLE)
        if(XCODE)
            set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
            set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
            if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
            endif()
        endif()
    else() # not APPLE
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()

    if(inst_set_to_use)
        string(TOLOWER ${inst_set_to_use} GCC_INST_SET)
        set(GCC_INST_SET "-m${GCC_INST_SET}")
    else()
        set(GCC_INST_SET)
    endif()

    # Get the gcc or clang version number in major.minor.build format
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)

    # Unrolling fixed-count loops was a useful optimization for Simmatrix
    # in earlier gcc versions.
    # Doesn't have a big effect for current compiler crop and may be
    # pushing our luck with optimizer bugs. So let the compilers decide
    # how to handle loops instead.
    ##set(GCC_OPT_ENABLE "-funroll-loops")

    # If you know of optimization bugs that affect Simbody in particular
    # gcc versions, this is the place to turn off those optimizations.
    set(GCC_OPT_DISABLE)

    # C++
    set(BUILD_CXX_FLAGS_DEBUG          "-g ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_RELEASE
      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_RELWITHDEBINFO
    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os ${GCC_INST_SET}")

    # C
    set(BUILD_C_FLAGS_DEBUG            "-g ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_RELEASE
      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_RELWITHDEBINFO
    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os ${GCC_INST_SET}")

    # C++
    set(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    # C
    set(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)

    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

endif()

# Determine which math libraries to use for this platform.
# Intel MKL: mkl_intel_c_dll;mkl_sequential_dll;mkl_core_dll
set(BUILD_USING_OTHER_LAPACK "" CACHE STRING
 "If you have your own Lapack and Blas, put libraries here, separated by semicolons (full paths or paths that are on the (DY)LD_LIBRARY_PATH (UNIX) or PATH (Windows)). See LAPACK_BEING_USED to see what's actually being used.")
if(WIN32)
    set(LAPACK_PLATFORM_DEFAULT lapack;blas)
else()
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND LAPACK_FOUND)
        set(LAPACK_PLATFORM_DEFAULT ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    else()
        message(WARNING "Could not find blas/lapack")
    endif()
endif()
set(LAPACK_BEING_USED ${LAPACK_PLATFORM_DEFAULT} CACHE STRING
    "Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)

if(BUILD_USING_OTHER_LAPACK)
    set(LAPACK_BEING_USED ${BUILD_USING_OTHER_LAPACK} CACHE STRING
"Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)
endif()

if(UNIX)
    if(NOT APPLE)
        set(REALTIME_LIB rt)
        # We do not use pthreads directly anymore, but on Linux, C++11 threads
        # are implemented with pthreads, and we must link to it.
        set(PTHREAD_LIB pthread)
    endif()
    set(MATH_LIBRARY m)
    set(DL_LIBRARY dl)
endif()

set(MATH_LIBS_TO_USE    ${LAPACK_BEING_USED} ${PTHREAD_LIB}
                        ${REALTIME_LIB} ${DL_LIBRARY} ${MATH_LIBRARY})
set(MATH_LIBS_TO_USE_VN ${MATH_LIBS_TO_USE})


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

