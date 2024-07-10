### Python, for ${Python_EXECUTABLE}
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
### Pybind11
find_package(pybind11 CONFIG REQUIRED)


# 根据编译器设置不同的编译标志
if (MSVC)
    # MSVC 编译器设置
    add_compile_options(/W4)
    add_compile_options(/wd4251) # 例如，禁用特定警告
else()
    # GCC 或 Clang 编译器设置
    add_compile_options(-Wall -Wextra -Wno-reorder -Wno-unused-parameter)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake) # main (top) cmake dir

set(CMAKE_VERBOSE_MAKEFILE on)

STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")


# 指定源文件夹和目标文件夹
set(SOURCE_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/gym-flap_fish)
set(DESTINATION_FOLDER ${LIBRARY_OUTPUT_PATH})
# 创建目标文件夹
file(MAKE_DIRECTORY ${DESTINATION_FOLDER})
# 复制源文件夹内容到目标文件夹
file(COPY ${SOURCE_FOLDER} DESTINATION ${DESTINATION_FOLDER})

set(SOURCE_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/gym-flap_fish_test)
set(DESTINATION_FOLDER ${LIBRARY_OUTPUT_PATH})
# 创建目标文件夹
file(MAKE_DIRECTORY ${DESTINATION_FOLDER})
# 复制源文件夹内容到目标文件夹
file(COPY ${SOURCE_FOLDER} DESTINATION ${DESTINATION_FOLDER})

set(SOURCE_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/drl_tool)
set(DESTINATION_FOLDER ${LIBRARY_OUTPUT_PATH})
# 创建目标文件夹
file(MAKE_DIRECTORY ${DESTINATION_FOLDER})
file(COPY ${SOURCE_FOLDER} DESTINATION ${DESTINATION_FOLDER})

aux_source_directory(. DIR_SRCS)
pybind11_add_module(${PROJECT_NAME} ${DIR_SRCS})
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
target_link_libraries(${PROJECT_NAME} PRIVATE sphinxsys_2d)

# test项目:
add_test(NAME ${PROJECT_NAME} COMMAND ${Python3_EXECUTABLE} "${LIBRARY_OUTPUT_PATH}/drl_tool/install.py")
set_tests_properties(${PROJECT_NAME} PROPERTIES WORKING_DIRECTORY "${LIBRARY_OUTPUT_PATH}/drl_tool"
    PASS_REGULAR_EXPRESSION "success")