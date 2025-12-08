# modules/Palace/tests/test_palace_external_call/CMakeLists.txt

# 1. Ensure PALACE_ROOT is defined.
#    This should be provided either via:
#       - cmake -DPALACE_ROOT=/path/to/palace
#       - or the environment variable PALACE_ROOT.
#    The test cannot proceed without this information.
if (NOT PALACE_ROOT)
    message(FATAL_ERROR
        "PALACE_ROOT is not set.\n"
        "Please set environment variable PALACE_ROOT, e.g.:\n"
        "  export PALACE_ROOT=/path/to/palace\n"
        "or pass -DPALACE_ROOT=/path/to/palace when running CMake."
    )
endif()


# 2. If user did not supply PALACE_EXECUTABLE, construct a default one
#    based on the typical build output path inside Palace.
if (NOT PALACE_EXECUTABLE)
    set(PALACE_EXECUTABLE
        "${PALACE_ROOT}/build/bin/palace"
        CACHE FILEPATH "Path to Palace main executable"
    )
endif()

# 3. Define the test executable that directly calls the external palace binary.
add_executable(test_palace_external_call
    test_palace_external_call.cpp
)

# 4. Link Palace interface layer so the executable has access to
#    the Palace static libs and header interfaces.
target_link_libraries(test_palace_external_call
    PRIVATE
        sphinxsys_palace
)

# 5. Embed the path to the Palace executable as a macro into C++,
#    so the code can reference the runtime binary.
target_compile_definitions(test_palace_external_call
    PRIVATE
        PALACE_EXECUTABLE=\"${PALACE_EXECUTABLE}\"
)

# 6. Copy input data directory (e.g., rings.json, mesh/) to the build dir
#    so the executable can run without needing absolute paths.
add_custom_command(TARGET test_palace_external_call POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

