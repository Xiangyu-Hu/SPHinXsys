# - Find ThreadingBuildingBlocks include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TBB
#    [REQUIRED]             # Fail with error if TBB is not found
#    )                      #
# Once done, this will define
#
#  TBB_FOUND - system has TBB
#  TBB_INCLUDE_DIR - the TBB include directories
#  TBB_LIBRARIES - TBB libraries to be lined, doesn't include malloc or
#                  malloc proxy
#
#  TBB_VERSION_MAJOR - Major Product Version Number
#  TBB_VERSION_MINOR - Minor Product Version Number
#  TBB_INTERFACE_VERSION - Engineering Focused Version Number
#  TBB_COMPATIBLE_INTERFACE_VERSION - The oldest major interface version
#                                     still supported. This uses the engineering
#                                     focused interface version numbers.
#
# This module reads hints about search locations from variables:
#  ENV TBB_ARCH_PLATFORM - for eg. set it to "mic" for Xeon Phi builds
#  ENV TBB_HOME or just TBB_HOME - root directory of tbb installation
#  ENV TBB_BUILD_PREFIX - specifies the build prefix for user built tbb
#                         libraries. Should be specified with ENV TBB_HOME
#                         and optionally...
#  ENV TBB_BUILD_DIR - if build directory is different than ${TBB_HOME}/build
#
#
# Modified by Robert Maynard from the original OGRE source
#
#-------------------------------------------------------------------
# This file is part of the CMake build system for OGRE
#     (Object-oriented Graphics Rendering Engine)
# For the latest info, see http://www.ogre3d.org/
#
# The contents of this file are placed in the public domain. Feel
# free to make use of it in any way you like.
#-------------------------------------------------------------------
#
#=============================================================================
# Copyright 2010-2012 Kitware, Inc.
# Copyright 2012      Rolf Eike Beer <eike@sf-mail.de>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)


#=============================================================================
#  FindTBB helper functions and macros
#

#===============================================
# Do the final processing for the package find.
#===============================================
macro(findpkg_finish PREFIX)
  # skip if already processed during this run
  if (NOT ${PREFIX}_FOUND)
    if (${PREFIX}_INCLUDE_DIR AND ${PREFIX}_LIBRARY)
      set(${PREFIX}_FOUND TRUE)
      set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIR})
      set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARY})
    else ()
      if (${PREFIX}_FIND_REQUIRED AND NOT ${PREFIX}_FIND_QUIETLY)
        message(FATAL_ERROR "Required library ${PREFIX} not found.")
      endif ()
    endif ()

  #mark the following variables as internal variables
  mark_as_advanced(${PREFIX}_INCLUDE_DIR
                    ${PREFIX}_LIBRARY
                    ${PREFIX}_LIBRARY_DEBUG
                    ${PREFIX}_LIBRARY_RELEASE)
  endif ()
endmacro()

#===============================================
# Generate debug names from given release names
#===============================================
macro(get_debug_names PREFIX)
  foreach(i ${${PREFIX}})
    set(${PREFIX}_DEBUG ${${PREFIX}_DEBUG} ${i}d ${i}D ${i}_d ${i}_D ${i}_debug ${i})
  endforeach()
endmacro()

#===============================================
# Couple a set of release AND debug libraries
#===============================================
macro(make_library_set PREFIX)
  if (${PREFIX}_RELEASE AND ${PREFIX}_DEBUG)
    set(${PREFIX} optimized ${${PREFIX}_RELEASE} debug ${${PREFIX}_DEBUG})
  elseif (${PREFIX}_RELEASE)
    set(${PREFIX} ${${PREFIX}_RELEASE})
  elseif (${PREFIX}_DEBUG)
    set(${PREFIX} ${${PREFIX}_DEBUG})
  endif ()
endmacro()


#=============================================================================
#  Now to actually find TBB
#

# Get path, convert backslashes as ${ENV_${var}}
SET(ENV_TBB_HOME, ${TBB_HOME})

# initialize search paths
set(TBB_PREFIX_PATH ${TBB_HOME} ${ENV_TBB_HOME})
set(TBB_INC_SEARCH_PATH "")
set(TBB_LIB_SEARCH_PATH "")

# For Windows, let's assume that the user might be using the precompiled
# TBB packages from the main website. These use a rather awkward directory
# structure (at least for automatically finding the right files) depending
# on platform and compiler, but we'll do our best to accommodate it.
# Not adding the same effort for the precompiled linux builds, though. Those
# have different versions for CC compiler versions and linux kernels which
# will never adequately match the user's setup, so there is no feasible way
# to detect the "best" version to use. The user will have to manually
# select the right files. (Chances are the distributions are shipping their
# custom version of tbb, anyway, so the problem is probably nonexistant.)
if (WIN32 AND MSVC)
  set(COMPILER_PREFIX "vc14")
  if (MSVC_VERSION EQUAL 1400)
    set(COMPILER_PREFIX "vc8")
  elseif(MSVC_VERSION EQUAL 1500)
    set(COMPILER_PREFIX "vc9")
  elseif(MSVC_VERSION EQUAL 1600)
    set(COMPILER_PREFIX "vc10")
  elseif(MSVC_VERSION EQUAL 1700)
    set(COMPILER_PREFIX "vc11")
  elseif(MSVC_VERSION EQUAL 1800)
    set(COMPILER_PREFIX "vc12")
  elseif(MSVC_VERSION EQUAL 1900)
    set(COMPILER_PREFIX "vc14")
  endif ()

  # for each prefix path, add ia32/64\${COMPILER_PREFIX}\lib to the lib search path
  foreach (dir IN LISTS TBB_PREFIX_PATH)
    if (CMAKE_CL_64)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia64/${COMPILER_PREFIX}/lib)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia64/${COMPILER_PREFIX})
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/${COMPILER_PREFIX}/lib)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64/${COMPILER_PREFIX})
    else ()
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/${COMPILER_PREFIX}/lib)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32/${COMPILER_PREFIX})
    endif ()
  endforeach ()
endif ()


# check compiler ABI
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_PREFIX)
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    list(APPEND COMPILER_PREFIX "gcc4.7")
  endif()
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
    list(APPEND COMPILER_PREFIX "gcc4.4")
  endif()
  list(APPEND COMPILER_PREFIX "gcc4.1")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(COMPILER_PREFIX)
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
    list(APPEND COMPILER_PREFIX "gcc4.7")
  endif()
  list(APPEND COMPILER_PREFIX "gcc4.4")
else() # Assume compatibility with 4.4 for other compilers
  list(APPEND COMPILER_PREFIX "gcc4.4")
endif ()

# if platform architecture is explicitly specified
set(TBB_ARCH_PLATFORM $ENV{TBB_ARCH_PLATFORM})
if (TBB_ARCH_PLATFORM)
  foreach (dir IN LISTS TBB_PREFIX_PATH)
    list(APPEND TBB_LIB_SEARCH_PATH ${dir}/${TBB_ARCH_PLATFORM}/lib)
    list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/${TBB_ARCH_PLATFORM})
  endforeach ()
endif ()

foreach (dir IN LISTS TBB_PREFIX_PATH)
  foreach (prefix IN LISTS COMPILER_PREFIX)
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64/${prefix})
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/lib)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/${prefix}/lib)
    else ()
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32/${prefix})
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/lib)
      list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/${prefix}/lib)
    endif ()
  endforeach()
endforeach ()

# add general search paths
foreach (dir IN LISTS TBB_PREFIX_PATH)
  list(APPEND TBB_LIB_SEARCH_PATH ${dir} ${dir}/lib ${dir}/Lib ${dir}/lib/tbb
    ${dir}/Libs)
  list(APPEND TBB_INC_SEARCH_PATH /usr/include/ ${dir}/include ${dir}/Include
    ${dir}/include/tbb)
endforeach ()

set(TBB_LIBRARY_NAMES tbb tbbmalloc tbbmalloc_proxy)
get_debug_names(TBB_LIBRARY_NAMES)

find_path(TBB_INCLUDE_DIR
          NAMES tbb/tbb.h
          PATHS ${TBB_INC_SEARCH_PATH} 
    NO_DEFAULT_PATH)

IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    find_library(TBB_LIB_DIR_TEMP 
    NAMES ${TBB_LIBRARY_NAMES_DEBUG}
    PATHS ${TBB_LIB_SEARCH_PATH}
    NO_DEFAULT_PATH)   
ELSE()
    find_library(TBB_LIB_DIR_TEMP 
    NAMES ${TBB_LIBRARY_NAMES}
    PATHS ${TBB_LIB_SEARCH_PATH}
    NO_DEFAULT_PATH)
ENDIF()

get_filename_component(TBB_LIB_DIR ${TBB_LIB_DIR_TEMP} DIRECTORY)

# Set composite TBB_LIBRARIES variable
set(TBB_LIBS)
set(TBB_RELEASE_LIBS)
foreach(libname IN LISTS TBB_LIBRARY_NAMES)
    set(TBB_RELEASE_LIBS ${TBB_RELEASE_LIBS} "${libname}")
endforeach()

set(TBB_DEBUG_LIBS)
foreach(libname IN LISTS TBB_LIBRARY_NAMES)
    set(TBB_DEBUG_LIBS ${TBB_DEBUG_LIBS} "${libname}_debug")
endforeach()

IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    if(TBB_DEBUG_LIBS)
  set(TBB_LIBS ${TBB_DEBUG_LIBS})
    else()
  set(TBB_LIBS ${TBB_RELEASE_LIBS})
    endif()    
ELSE()
    if(TBB_RELEASE_LIBS)
        set(TBB_LIBS ${TBB_RELEASE_LIBS})
    else()
  set(TBB_LIBS ${TBB_DEBUG_LIBS})
    endif()  
ENDIF()

set(TBB_LIBRARYS ${TBB_LIBS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB DEFAULT_MSG 
  TBB_INCLUDE_DIR TBB_LIB_DIR TBB_LIBRARYS)

#=============================================================================

if(TBB_FOUND)
  set(TBB_LIBRARYS ${TBB_LIBRARYS} )
  set(TBB_INCLUDE_DIR ${TBB_INCLUDE_DIR} )
  set(TBB_LIB_DIR ${TBB_LIB_DIR} )
  set(TBB_DEFINITIONS)
endif()

unset(TBB_LIB_DIR_TEMP CACHE)

mark_as_advanced(
    TBB_INCLUDE_DIR
    TBB_LIB_DIR
    TBB_LIBRARYS
)