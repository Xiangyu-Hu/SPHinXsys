# modules/Palace/tests/test_3d_palace_rings/CMakeLists.txt

# Build the test executable for the 3D Palace magnetostatic rings example.
add_executable(test_3d_palace_rings
    test_3d_palace_rings.cpp
)

# Link Palace interface module so the case can invoke Palace solver internally.
target_link_libraries(test_3d_palace_rings
    PRIVATE
        sphinxsys_palace
)

# Default path to the rings.json input configuration.
# This file will be copied later into the runtime directory.
set(PALACE_JSON_PATH_DEFAULT "data/rings.json")

# Allow user to override PALACE_JSON_PATH from CMake command line.
if (NOT PALACE_JSON_PATH)
    set(PALACE_JSON_PATH "${PALACE_JSON_PATH_DEFAULT}")
endif()

# Pass the JSON path as a macro into the C++ code.
target_compile_definitions(test_3d_palace_rings
    PRIVATE
        PALACE_JSON_PATH=\"${PALACE_JSON_PATH}\"
)

# Copy data directory (including rings.json + mesh/) so that the test
# can run from the build directory without requiring absolute paths.
add_custom_command(TARGET test_3d_palace_rings POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

