set(GPU_PAPER_BENCHMARK_ROOT "${CMAKE_CURRENT_LIST_DIR}")
set(PAPER_BENCHMARK_COMMON_ROOT
    "${CMAKE_CURRENT_LIST_DIR}/../paper_benchmarks_common")
get_filename_component(GPU_PAPER_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/../.."
                       ABSOLUTE)
set(GPU_PAPER_RESULT_ROOT "${GPU_PAPER_BUILD_ROOT}/benchmark_results")

execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GPU_PAPER_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE GPU_PAPER_GIT_RESULT)
if(NOT GPU_PAPER_GIT_RESULT EQUAL 0)
    set(GPU_PAPER_GIT_COMMIT "unknown")
endif()

if(CMAKE_BUILD_TYPE)
    set(GPU_PAPER_BUILD_ID "${CMAKE_BUILD_TYPE}")
else()
    set(GPU_PAPER_BUILD_ID "multi-config-or-unspecified")
endif()
set(GPU_PAPER_COMPILER
    "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")

function(gpu_paper_add_benchmark target_name dimension)
    set(options)
    set(one_value_args ASSET_DIR)
    cmake_parse_arguments(BENCH "${options}" "${one_value_args}" "" ${ARGN})

    file(GLOB benchmark_sources CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
    add_executable(${target_name} ${benchmark_sources})
    target_link_libraries(${target_name} PRIVATE sphinxsys_${dimension})
    target_compile_features(${target_name} PRIVATE cxx_std_17)
    target_include_directories(${target_name} PRIVATE
        "${PAPER_BENCHMARK_COMMON_ROOT}")
    target_compile_definitions(${target_name} PRIVATE
        "SPHINXSYS_BENCHMARK_GIT_COMMIT=\"${GPU_PAPER_GIT_COMMIT}\""
        "SPHINXSYS_BENCHMARK_BUILD=\"${GPU_PAPER_BUILD_ID}\""
        "SPHINXSYS_BENCHMARK_COMPILER=\"${GPU_PAPER_COMPILER}\""
        "SPHINXSYS_BENCHMARK_RESULT_ROOT=\"${GPU_PAPER_RESULT_ROOT}\"")

    set(benchmark_runtime_dir "${CMAKE_CURRENT_BINARY_DIR}/bin")
    set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${benchmark_runtime_dir}"
        VS_DEBUGGER_WORKING_DIRECTORY "${benchmark_runtime_dir}")

    if(BENCH_ASSET_DIR)
        file(MAKE_DIRECTORY "${benchmark_runtime_dir}/input")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${BENCH_ASSET_DIR}/"
            DESTINATION "${benchmark_runtime_dir}/input")
    endif()
endfunction()

add_subdirectory(verification)
