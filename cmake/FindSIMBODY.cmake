# FindSimbody.cmake
#
# Try to find the Simbody multibody dynamics package that is part of the
# SimTK biosimulation toolkit. This includes several libraries:
#     libSimTKsimbody
#     libSimTKmath
#     libSimTKcommon
#     libSimTKlapack
#     (and on Windows pthreads too)
# We expect to find all of these in the same place.
#
# This module looks for certain CMake variables on input and behaves 
# accordingly if they are present:
#
#   SimTK_SDK 	
#       If this is set (probably by the nightly build system) then that
#       will be the only place we'll look for Simbody. Otherwise we'll 
#       hunt around.
#   BUILD_USING_NAMESPACE
#       This is used to look for Simbody libraries that are prefixed by
#       a string. An underscore is appended to the "namespace".
#   BUILD_UNVERSIONED_LIBRARIES
#       If this is set we look for the non-versioned Simbody libraries.
#   BUILD_VERSIONED_LIBRARIES
#       If this is set we'll look for Simbody libraries that have the same
#       version number as the version of Molmodel we're building.
#
# Note that both versioned and unversioned libraries may be produced by the
# same build so we have to allow for both the above to be set.
#
# Created January 2011 by Michael Sherman
# Simbios National Center for Physics Based Simulation of Biological Structures
# Stanford University
#
# Once done this will define:
#
#   Simbody_FOUND - Whether search for Simbody libraries and headers succeeded.
#   Simbody_ROOT_DIR - the installation directory; all the pieces must be
#                      found together
#   Simbody_INCLUDE_DIR - location of Simbody.h
#   Simbody_LIB_DIR     - location of libSimbody.{a,so,dylib} or Simbody.lib
#   Simbody_BIN_DIR     - location of VisualizerGUI and .dll's on Windows
#   Simbody_LIBRARIES   - suitable for target_link_libraries(); includes
#                           both optimized and debug libraries if both are
#                           available
#   Simbody_STATIC_LIBRARIES - suitable for target_link_libraries(); includes
#                              both optimized and debug static libraries if
#                              both are available
#
# == Using Simbody ==
#
#     find_package(Simbody)
#     if(Simbody_FOUND)
#         include_directories(${Simbody_INCLUDE_DIR})
#         link_directories(${Simbody_LIB_DIR})
#         target_link_libraries(Foo ${Simbody_LIBRARIES})
#     endif()
#

#set(Simbody_ROOT_DIR "$ENV{SIMBODY_HOME}")
set(Simbody_ROOT_DIR "${SIMBODY_HOME}")
set(Simbody_SEARCH_PATHS "${Simbody_ROOT_DIR}")

# We'll use the main Simbody header as the key to finding the installation
# root directory. We're assuming the header is in 
#    Simbody_INCLUDE_DIR == ${Simbody_ROOT_DIR}/include
# So stripping off the "/include" should give the root directory.

# Force this to be recalculated every time.
#set(Simbody_INCLUDE_DIR "Simbody_INCLUDE_DIR-NOTFOUND" CACHE PATH
#   "The Simbody and SimTK include directory." FORCE)

foreach(pth IN LISTS Simbody_SEARCH_PATHS)
find_path(Simbody_INCLUDE_DIR 
    NAMES "SimTKsimbody.h" "Simbody.h"
    PATHS "${pth}"
    PATH_SUFFIXES "include" "SimTK/include" "simtk/include"
                "Simbody/include" "simbody/include" "include/simbody"
    NO_DEFAULT_PATH
    DOC "Location of top-level installed Simbody header files"
)
endforeach()

# This will only be executed if the first loop fails. We're getting
# desperate!
find_path(Simbody_INCLUDE_DIR 
    NAMES "SimTKsimbody.h" "Simbody.h"
    PATH_SUFFIXES "include" "SimTK/include" "simtk/include"
                "Simbody/include" "simbody/include" "include/simbody"
    DOC "Location of top-level installed Simbody header files"
)

set(Simbody_LIBRARY_LIST SimTKsimbody SimTKmath SimTKcommon)

find_library(Simbody_LIB_DIR_TEMP 
    NAMES ${Simbody_LIBRARY_LIST}
    PATHS ${Simbody_ROOT_DIR}/lib ${Simbody_ROOT_DIR}/lib64 ${Simbody_ROOT_DIR}/lib/x86_64-linux-gnu
    NO_DEFAULT_PATH)

    get_filename_component(Simbody_LIB_DIR ${Simbody_LIB_DIR_TEMP} DIRECTORY)

set(Simbody_LAPACK_LIBRARY_LIST )
set(Simbody_EXTRA_LIBRARY_LIST )

if (WIN32)
    set(Simbody_LAPACK_LIBRARY_LIST liblapack; libblas)
    set(Simbody_EXTRA_LIBRARY_LIST pthreadVC2_x64)
else(WIN32)
    set(Simbody_LAPACK_LIBRARY_LIST lapack; blas)
    if(NOT APPLE)
        set(Simbody_EXTRA_LIBRARY_LIST pthread; rt; dl)
    else(NOT APPLE)
        set(Simbody_EXTRA_LIBRARY_LIST pthread; dl)
    endif(NOT APPLE)
endif(WIN32)

# Find out which of the unversioned libraries are available.
find_library(Simbody_LIBRARY NAMES SimTKsimbody
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody library."
    NO_DEFAULT_PATH)
find_library(Simbody_DEBUG_LIBRARY NAMES SimTKsimbody_d
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody debug library."
    NO_DEFAULT_PATH)

# Set composite Simbody_LIBRARIES variable
set(LIBS)
set(SIMBODYRELEASELIBS)
if(Simbody_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(SIMBODYRELEASELIBS ${SIMBODYRELEASELIBS} "${libname}")
    endforeach()
endif()

set(SIMBODYDEBUGLIBS)
if(Simbody_DEBUG_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(SIMBODYDEBUGLIBS ${SIMBODYDEBUGLIBS} "${libname}_d")
    endforeach()
endif()

IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    if(SIMBODYDEBUGLIBS)
    set(LIBS ${SIMBODYDEBUGLIBS})
    else()
    set(LIBS ${SIMBODYRELEASELIBS})
    endif()
ELSE()
    if(SIMBODYRELEASELIBS)
        set(LIBS ${SIMBODYRELEASELIBS})
    else()
        set(LIBS ${SIMBODYDEBUGLIBS})
    endif()
ENDIF()

if (LIBS)
    foreach(lapack_lib IN LISTS Simbody_LAPACK_LIBRARY_LIST)
        set(LIBS ${LIBS} "${lapack_lib}")
        set(SIMBODYRELEASELIBS ${SIMBODYRELEASELIBS} "${lapack_lib}")
        set(SIMBODYDEBUGLIBS ${SIMBODYDEBUGLIBS} "${lapack_lib}")
    endforeach()
    foreach(extra_lib IN LISTS Simbody_EXTRA_LIBRARY_LIST)
        set(LIBS ${LIBS} "${extra_lib}")
        set(SIMBODYRELEASELIBS ${SIMBODYRELEASELIBS} "${extra_lib}")
        set(SIMBODYDEBUGLIBS ${SIMBODYDEBUGLIBS} "${extra_lib}")
    endforeach()
    set(Simbody_LIBRARIES ${LIBS} CACHE STRING 
        "Simbody dynamic libraries" FORCE)
    set(Simbody_DEBUG_LIBRARIES ${SIMBODYDEBUGLIBS} CACHE STRING 
        "Simbody debug dynamic libraries" FORCE)
    set(Simbody_RELEASE_LIBRARIES ${SIMBODYRELEASELIBS} CACHE STRING 
        "Simbody release dynamic libraries" FORCE)
else()
    set(Simbody_LIBRARIES Simbody_LIBRARIES-NOTFOUND CACHE STRING 
        "Simbody dynamic libraries" FORCE)
endif()

include(FindPackageHandleStandardArgs OPTIONAL)
find_package_handle_standard_args(Simbody DEFAULT_MSG 
    Simbody_INCLUDE_DIR Simbody_LIB_DIR Simbody_LIBRARIES Simbody_DEBUG_LIBRARIES Simbody_RELEASE_LIBRARIES )

unset(Simbody_LIBRARY CACHE)
unset(Simbody_DEBUG_LIBRARY CACHE)
unset(Simbody_LIB_DIR_TEMP CACHE)

# Not all the variables we produced need be returned.

mark_as_advanced(
    Simbody_INCLUDE_DIR
    Simbody_BIN_DIR
    Simbody_LIB_DIR
    Simbody_LIBRARIES 
    Simbody_RELEASE_LIBRARIES
    Simbody_DEBUG_LIBRARIES
)