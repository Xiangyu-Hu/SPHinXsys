option(BUILD_2D_EXAMPLES "BUILD_2D_EXAMPLES" 1)
option(BUILD_3D_EXAMPLES "BUILD_3D_EXAMPLES" 1)
option(BUILD_UNIT_TESTS "BUILD_UNIT_TESTS" 1)

if(NOT EMSCRIPTEN)
	find_package(GTest REQUIRED)
	include(GoogleTest)

	if(ONLY_3D)
		set(BUILD_2D_EXAMPLES 0)
	endif()

	if(BUILD_2D_EXAMPLES)
		ADD_SUBDIRECTORY(2d_examples)
	endif()

	if(BUILD_UNIT_TESTS)
		ADD_SUBDIRECTORY(unit_tests_src)
	endif()

	if(NOT SPH_ONLY_STATIC_BUILD)
		add_subdirectory(modules)
	endif()
endif()

if(BUILD_3D_EXAMPLES)
	ADD_SUBDIRECTORY(3d_examples)
	ADD_SUBDIRECTORY(webassembly_models)
endif()
