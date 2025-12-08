# modules/Palace/tests/test_palace_main_orig/CMakeLists.txt

add_executable(test_palace_main_orig
    main.cpp
)


target_link_libraries(test_palace_main_orig
    PRIVATE
        Palace::palace
)


add_custom_command(TARGET test_palace_main_orig POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../test_3d_palace_rings/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
)

