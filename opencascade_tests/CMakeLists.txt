option(SPHINXSYS_BUILD_OPENCASCADE_TESTS  "SPHINXSYS_BUILD_OPENCASCADE_TESTS"  ON)

find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

target_link_libraries(sphinxsys_core INTERFACE GTest::gtest GTest::gtest_main)


if(SPHINXSYS_BUILD_UNIT_TESTS)
	ADD_SUBDIRECTORY(opencascade_tests)
endif()





