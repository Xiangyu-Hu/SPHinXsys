#ifndef SPHINXSYS_EXECUTION_CONTEXT_HPP
#define SPHINXSYS_EXECUTION_CONTEXT_HPP

#include <sycl/sycl.hpp>

namespace SPH
{
namespace execution
{
class ParallelSYCLDevicePolicy;

template <typename ExecutionUnit>
struct Context
{
    template <class... Args>
    explicit Context(Args &...args){};

    template <class... Args>
    explicit Context(Args &&...args){};
};

template <>
struct Context<ParallelSYCLDevicePolicy>
{
    explicit Context(sycl::handler &cgh) : cgh(cgh) {}

    sycl::handler &cgh; /**< Command group handler*/
};
} // namespace execution
} // namespace SPH

#endif // SPHINXSYS_EXECUTION_CONTEXT_HPP
