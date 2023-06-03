#ifndef SPHINXSYS_EXECUTION_PROXY_HPP
#define SPHINXSYS_EXECUTION_PROXY_HPP

#include "execution_context.hpp"
#include "execution_policy.h"

namespace SPH
{
namespace execution
{
template <typename BaseT, typename KernelT>
class ExecutionProxy
{
  public:
    using Base = BaseT;
    using Kernel = KernelT;

    ExecutionProxy(BaseT *base, KernelT *proxy) : base(base), kernel(proxy) {}

    template <class ExecutionPolicy = ParallelPolicy>
    BaseT *get(const ExecutionPolicy & = par) const
    {
        return base;
    }

    KernelT *get(const ParallelSYCLDevicePolicy &) const
    {
        return kernel;
    }

    template <class ExecutionPolicy>
    void init_memory_access(const Context<ExecutionPolicy> &context)
    {
        if constexpr (std::is_same_v<ExecutionPolicy, ParallelSYCLDevicePolicy>)
            this->init_device_memory_access(context.cgh);
    }

  protected:
    virtual void init_device_memory_access(sycl::handler &cgh) = 0;

    BaseT *base;
    KernelT *kernel;
};

template <typename T>
class NoProxy : public ExecutionProxy<T, T>
{
  public:
    explicit NoProxy(T *base) : ExecutionProxy<T, T>(base, base) {}

  protected:
    void init_device_memory_access(sycl::handler &cgh) override
    {
        static_assert("No device memory access available for this class.");
    }
};
} // namespace execution
} // namespace SPH

#endif // SPHINXSYS_EXECUTION_PROXY_HPP
