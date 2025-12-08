# modules/Palace/CMakeLists.txt

# 顶层如果没开 Palace，就直接跳过这个模块
if (NOT SPHINXSYS_WITH_PALACE)
    return()
endif()

# 这里用 INTERFACE 库：不需要源文件，只当“配置转发器”
add_library(sphinxsys_palace INTERFACE
    #palace_interface.cpp
)

# INTERFACE 库，必须用 INTERFACE 关键字
target_link_libraries(sphinxsys_palace
    INTERFACE
    #PUBLIC
        Palace::palace
)

target_include_directories(sphinxsys_palace
    INTERFACE
    #PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# 测试子目录
add_subdirectory(tests)

