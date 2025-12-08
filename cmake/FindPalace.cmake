# FindPalace.cmake  — Locate and configure the Palace static library for SPHinXsys
#
# Usage in your top-level CMakeLists.txt:
#
#   # Optionally (recommended) set the root of the Palace source/build tree:
#   #   - via CMake cache:
#   #       -DPALACE_ROOT=/path/to/palace
#   #     or
#   #       -DPalace_ROOT=/path/to/palace
#   #   - or via environment variable:
#   #       export PALACE_ROOT=/path/to/palace
#   #
#   find_package(Palace REQUIRED)
#   target_link_libraries(<your_target> PUBLIC Palace::palace)
#
# This module defines:
#   Palace::palace          - INTERFACE target for Palace and all its dependencies
#   Palace_FOUND            - TRUE if Palace was found and configured
#   PALACE_ROOT             - Root directory of the Palace source/build tree
#   PALACE_INCLUDE_DIRS     - Include directories needed to use Palace
#   PALACE_LIB_DIR          - Directory containing libpalace.a and related static libs
#   PALACE_LIBRARIES        - List of libraries linked by Palace::palace

# Guard: if the target already exists, we assume Palace is configured.
if (TARGET Palace::palace)
    set(Palace_FOUND TRUE)
    return()
endif()

# ------------------------------------------------------------------------------
# 1. Determine PALACE_ROOT
# ------------------------------------------------------------------------------

# Preferred CMake convention is <PackageName>_ROOT, i.e. Palace_ROOT.
# We also support PALACE_ROOT for convenience and the PALACE_ROOT environment
# variable. The resolution order is:
#   1) PALACE_ROOT (CMake cache or command line)
#   2) Palace_ROOT (CMake cache or command line)
#   3) $ENV{PALACE_ROOT}
if (NOT PALACE_ROOT)
    if (Palace_ROOT)
        set(PALACE_ROOT "${Palace_ROOT}" CACHE PATH "Palace root directory" FORCE)
    elseif (DEFINED ENV{PALACE_ROOT})
        set(PALACE_ROOT "$ENV{PALACE_ROOT}" CACHE PATH "Palace root directory" FORCE)
    else()
        set(PALACE_ROOT "" CACHE PATH "Palace root directory")
    endif()
endif()

if (NOT PALACE_ROOT)
    message(FATAL_ERROR
        "PALACE_ROOT is not configured, Palace cannot be found.\n"
        "Please add -DPALACE_ROOT=/path/to/palace (or -DPalace_ROOT=/path/to/palace) "
        "to your CMake command, or set the PALACE_ROOT environment variable."
    )
endif()

# Expected superbuild layout:
#   ${PALACE_ROOT}/palace         - Palace source headers (palace/...)
#   ${PALACE_ROOT}/build/include  - Installed headers (mfem.hpp, mfem/mfem.hpp, etc.)
#   ${PALACE_ROOT}/build/lib      - Static libraries (libpalace.a, libmfem.a, ...)
set(PALACE_INCLUDE_DIRS
    ${PALACE_ROOT}
    ${PALACE_ROOT}/palace            # for #include "drivers/...", "fem/...", etc.
    ${PALACE_ROOT}/build/include
)

set(PALACE_LIB_DIR ${PALACE_ROOT}/build/lib)

# Quick sanity check: libpalace.a must exist.
if (NOT EXISTS "${PALACE_LIB_DIR}/libpalace.a")
    message(FATAL_ERROR
        "Cannot find ${PALACE_LIB_DIR}/libpalace.a.\n"
        "Please ensure that PALACE_ROOT points to the root directory "
        "of your Palace source/build tree."
    )
endif()

# ------------------------------------------------------------------------------
# 2. Locate MPI (required by the MPI-enabled MFEM used inside Palace)
# ------------------------------------------------------------------------------

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ------------------------------------------------------------------------------
# 3. Define the INTERFACE target Palace::palace
# ------------------------------------------------------------------------------

# 3.1 Create a pure interface library that carries include dirs and link flags
add_library(palace_interface INTERFACE)

# 3.2 Palace / MFEM header paths
target_include_directories(palace_interface INTERFACE
    ${PALACE_INCLUDE_DIRS}
)

target_compile_definitions(palace_interface INTERFACE
    PALACE_WITH_SLEPC
)

# 3.3 All (static) third-party libraries that a typical Palace application needs.
#     These paths assume the superbuild layout under ${PALACE_LIB_DIR}.
target_link_libraries(palace_interface INTERFACE
    ${PALACE_LIB_DIR}/libpalace.a
    ${PALACE_LIB_DIR}/libslepc.a
    ${PALACE_LIB_DIR}/libpetsc.a

    ${PALACE_LIB_DIR}/libmfem.a
    ${PALACE_LIB_DIR}/libHYPRE.a
    ${PALACE_LIB_DIR}/libsuperlu_dist.a
    ${PALACE_LIB_DIR}/libparmetis.a
    ${PALACE_LIB_DIR}/libmetis.a

    ${PALACE_LIB_DIR}/libsundials_nvecserial.a
    ${PALACE_LIB_DIR}/libsundials_nvecparallel.a
    ${PALACE_LIB_DIR}/libsundials_nvecmpiplusx.a
    ${PALACE_LIB_DIR}/libsundials_cvodes.a
    ${PALACE_LIB_DIR}/libsundials_arkode.a
    ${PALACE_LIB_DIR}/libsundials_kinsol.a
    ${PALACE_LIB_DIR}/libsundials_core.a

    ${PALACE_LIB_DIR}/libgs.a
    ${PALACE_LIB_DIR}/libfmt.a
    ${PALACE_LIB_DIR}/libscn.a


    ${PALACE_LIB_DIR}/libceed.so

    z
    
    # MPI C++ interface target (provides <mpi.h> include dir and MPI link flags)
    MPI::MPI_CXX
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# 3.4 Public alias with the canonical namespace-qualified name
add_library(Palace::palace ALIAS palace_interface)

# ------------------------------------------------------------------------------
# 4. Standard find_package() variables
# ------------------------------------------------------------------------------

set(Palace_FOUND TRUE)
set(PALACE_LIBRARIES palace_interface)

