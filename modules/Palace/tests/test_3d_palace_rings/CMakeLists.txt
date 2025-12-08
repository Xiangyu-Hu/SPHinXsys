# modules/Palace/tests/test_3d_palace_rings/CMakeLists.txt

add_executable(test_3d_palace_rings
    test_3d_palace_rings.cpp
)

target_link_libraries(test_3d_palace_rings
    PRIVATE
        sphinxsys_palace
)

# 默认 JSON 路径（相对运行目录）
set(PALACE_JSON_PATH_DEFAULT "data/rings.json")

# 允许用户通过 -DPALACE_JSON_PATH=... 覆盖
if (NOT PALACE_JSON_PATH)
    set(PALACE_JSON_PATH "${PALACE_JSON_PATH_DEFAULT}")
endif()

target_compile_definitions(test_3d_palace_rings
    PRIVATE
        PALACE_JSON_PATH=\"${PALACE_JSON_PATH}\"
)

# 运行时需要 data/ 目录（rings.json + mesh）
add_custom_command(TARGET test_3d_palace_rings POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

