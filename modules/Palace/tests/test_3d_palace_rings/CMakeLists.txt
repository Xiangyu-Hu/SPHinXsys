# modules/Palace/tests/test_3d_palace_rings/CMakeLists.txt

add_executable(test_3d_palace_rings
    test_3d_palace_rings.cpp
)

target_link_libraries(test_3d_palace_rings
    PRIVATE
        sphinxsys_palace
)


set(PALACE_JSON_PATH_DEFAULT "data/rings.json")


if (NOT PALACE_JSON_PATH)
    set(PALACE_JSON_PATH "${PALACE_JSON_PATH_DEFAULT}")
endif()

target_compile_definitions(test_3d_palace_rings
    PRIVATE
        PALACE_JSON_PATH=\"${PALACE_JSON_PATH}\"
)


add_custom_command(TARGET test_3d_palace_rings POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

