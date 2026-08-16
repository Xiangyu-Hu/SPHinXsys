if(SPHINXSYS_USE_SYCL)
    set(gpu_paper_diffusion_assets regression_test_tool_sycl)
else()
    set(gpu_paper_diffusion_assets regression_test_tool_ck)
endif()

gpu_paper_add_benchmark(gpu_paper_diffusion_neumann_2d 2d
    ASSET_DIR ${gpu_paper_diffusion_assets})