# modules/Palace/tests/test_palace_external_call/CMakeLists.txt

# 1. 确认 PALACE_ROOT 已经由 FindPalace.cmake 设置好
#    （来源可以是 -DPALACE_ROOT=... 或环境变量 PALACE_ROOT）
if (NOT PALACE_ROOT)
    message(FATAL_ERROR
        "PALACE_ROOT is not set.\n"
        "Please set environment variable PALACE_ROOT, e.g.:\n"
        "  export PALACE_ROOT=/path/to/palace\n"
        "or pass -DPALACE_ROOT=/path/to/palace when running CMake."
    )
endif()

# 2. 如果用户没有显式指定 PALACE_EXECUTABLE，就从 PALACE_ROOT 推一个默认值
if (NOT PALACE_EXECUTABLE)
    set(PALACE_EXECUTABLE
        "${PALACE_ROOT}/build/bin/palace"
        CACHE FILEPATH "Path to Palace main executable"
    )
endif()

# 3. 定义测试可执行文件
add_executable(test_palace_external_call
    test_palace_external_call.cpp
)

# 4. 链接 sphinxsys_palace（它会把 Palace::palace 的所有库带上）
target_link_libraries(test_palace_external_call
    PRIVATE
        sphinxsys_palace
)

# 5. 把 PALACE_EXECUTABLE 路径以宏的形式传进 C++
target_compile_definitions(test_palace_external_call
    PRIVATE
        PALACE_EXECUTABLE=\"${PALACE_EXECUTABLE}\"
)

# 6. 拷贝 data 目录到构建目录，保证运行时有 rings.json 和 mesh/
add_custom_command(TARGET test_palace_external_call POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

