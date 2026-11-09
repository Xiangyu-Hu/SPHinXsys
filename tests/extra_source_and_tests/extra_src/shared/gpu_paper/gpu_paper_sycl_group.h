/**
 * @file gpu_paper_sycl_group.h
 * @brief Apply GPU_PAPER_GROUP to SYCL ExecutionInstance without modifying src/.
 * @details Call once at process start (before heavy kernels). getQueue() initially clamps
 *          work-group size to 64; we override afterward from env.
 */
#ifndef SPHINXSYS_GPU_PAPER_SYCL_GROUP_H
#define SPHINXSYS_GPU_PAPER_SYCL_GROUP_H

#include <cstdlib>
#include <iostream>
#include <string>

#if SPHINXSYS_USE_SYCL
#include "implementation_sycl.h"
#endif

namespace sphinxsys_gpu_paper
{

/** Returns the effective work-group size after applying env (0 if CPU / unset path). */
inline size_t applyGpuPaperSyclGroupFromEnv()
{
#if SPHINXSYS_USE_SYCL
    using SPH::execution::execution_instance;
    // Ensure queue exists (constructor path may set default 64).
    auto &q = execution_instance.getQueue();
    (void)q;

    size_t group = execution_instance.getWorkGroupSize();
    if (const char *v = std::getenv("GPU_PAPER_GROUP"))
    {
        const long parsed = std::atol(v);
        if (parsed > 0)
            group = static_cast<size_t>(parsed);
    }

    const auto device = execution_instance.getQueue().get_device();
    const size_t device_max = device.get_info<sycl::info::device::max_work_group_size>();
    // NVIDIA-style limit: registers/block ≤ 65536. Heavy SPH kernels ~64–96 regs/item
    // → safe local size often ≤ 512; paper scan uses ≤256. Cap hard at 512.
    constexpr size_t kPaperSafeMax = 512;
    if (group > device_max)
        group = device_max;
    if (group > kPaperSafeMax)
    {
        std::cout << "[gpu_paper] WARNING: GPU_PAPER_GROUP=" << group
                  << " exceeds paper safe max " << kPaperSafeMax
                  << " (register pressure); clamping.\n";
        group = kPaperSafeMax;
    }
    if (group < 1)
        group = 1;

    execution_instance.setWorkGroupSize(group);
    std::cout << "[gpu_paper] SYCL work_group_size=" << execution_instance.getWorkGroupSize()
              << " (GPU_PAPER_GROUP env, device_max=" << device_max
              << ", paper_safe_max=" << kPaperSafeMax << ")\n";
    return execution_instance.getWorkGroupSize();
#else
    (void)std::getenv("GPU_PAPER_GROUP");
    return 0;
#endif
}

inline std::string gpuPaperGroupMetaString()
{
#if SPHINXSYS_USE_SYCL
    return std::to_string(SPH::execution::execution_instance.getWorkGroupSize());
#else
    if (const char *v = std::getenv("GPU_PAPER_GROUP"))
        return std::string(v);
    return "na_cpu";
#endif
}

} // namespace sphinxsys_gpu_paper

#endif
