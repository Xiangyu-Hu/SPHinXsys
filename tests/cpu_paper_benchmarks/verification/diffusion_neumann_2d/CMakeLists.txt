if(SPHINXSYS_USE_SYCL)
    set(cpu_paper_diffusion_assets regression_test_tool_sycl)
else()
    set(cpu_paper_diffusion_assets regression_test_tool_ck)
endif()

cpu_paper_add_benchmark(cpu_paper_diffusion_neumann_2d 2d
    ASSET_DIR ${cpu_paper_diffusion_assets})