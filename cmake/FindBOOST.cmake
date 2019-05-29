# FindBoost.cmake
#

cmake_minimum_required(VERSION 2.8)

set(Boost_ROOT_DIR "${BOOST_HOME}")
set(Boost_SEARCH_PATHS "${Boost_ROOT_DIR}")

foreach(pth IN LISTS Boost_SEARCH_PATHS)
  find_path(Boost_INCLUDE_DIR 
    NAMES "mpi.hpp" 
    PATHS "${pth}"
    PATH_SUFFIXES "include" "include/boost"
    NO_DEFAULT_PATH
    DOC "Location of top-level installed Simbody header files"
  )
endforeach()

# This will only be executed if the first loop fails. We're getting
# desperate!
find_path(Boost_INCLUDE_DIR 
    NAMES "mpi.hpp" 
    PATH_SUFFIXES "include"
    DOC "Location of top-level installed Simbody header files"
)

set(Boost_LIBRARY_LIST boost_filesystem boost_system)

find_library(Boost_LIB_DIR_TEMP 
    NAMES ${Boost_LIBRARY_LIST}
    PATHS ${Boost_ROOT_DIR}/lib ${Boost_ROOT_DIR}/lib64
    NO_DEFAULT_PATH)

get_filename_component(Boost_LIB_DIR ${Boost_LIB_DIR_TEMP} DIRECTORY)

# Find out which of the unversioned libraries are available.

find_library(Boost_LIBRARY NAMES boost_filesystem
    PATHS ${Boost_LIB_DIR}
    DOC "This is the main Boost library."
    NO_DEFAULT_PATH)

# Set composite Simbody_LIBRARIES variable
set(BOOSTLIBS)
if(Boost_LIBRARY)
    foreach(libname IN LISTS Boost_LIBRARY_LIST)
        set(BOOSTLIBS ${BOOSTLIBS} "${libname}")
    endforeach()
endif()

set(Boost_LIBRARIES ${BOOSTLIBS} CACHE STRING 
        "Boost dynamic libraries" FORCE)

include(FindPackageHandleStandardArgs OPTIONAL)
find_package_handle_standard_args(Boost DEFAULT_MSG Boost_INCLUDE_DIR Boost_LIB_DIR Boost_LIBRARIES)

unset(Boost_LIBRARY CACHE)
unset(Simbody_LIB_DIR_TEMP CACHE)

# Not all the variables we produced need be returned.

mark_as_advanced(
    Boost_INCLUDE_DIR
    Boost_LIB_DIR
    Boost_LIBRARIES 
)
