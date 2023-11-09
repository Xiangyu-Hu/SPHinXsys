#ifndef SPHINXSYS_DEVICEIMPLEMENTATION_H
#define SPHINXSYS_DEVICEIMPLEMENTATION_H

#include "execution_policy.h"
#include <memory>
#include <sycl/sycl.hpp>

namespace SPH::execution
{
template <class Device>
class DeviceImplementation
{
  public:
    DeviceImplementation() = default;

    template <class... Args>
    explicit DeviceImplementation(Args &&...device_args)
        : device_(std::make_shared<Device>(std::forward<Args>(device_args)...)) {}

    using KernelType = Device;

    inline Device *get_ptr() const { return device_.get(); }

    inline sycl::buffer<Device> &get_buffer()
    {
        if (!device_buffer_)
            device_buffer_ = std::make_unique<sycl::buffer<Device>>(device_, 1);
        return *device_buffer_;
    }

  private:
    std::shared_ptr<Device> device_;
    std::shared_ptr<sycl::buffer<Device>> device_buffer_;
};
} // namespace SPH::execution

#endif // SPHINXSYS_DEVICEIMPLEMENTATION_H
