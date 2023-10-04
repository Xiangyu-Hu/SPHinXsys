#ifndef SPHINXSYS_DEVICE_EXECUTABLE_HPP
#define SPHINXSYS_DEVICE_EXECUTABLE_HPP

#include "execution_proxy.hpp"

namespace SPH::execution {
    template<class Base, class Kernel>
    class DeviceExecutable {
    public:
        template<class ...Args>
        DeviceExecutable(Args&& ...args) : device_proxy(args...) {}

        auto &getDeviceProxy() {
            return device_proxy;
        }

    private:
      ExecutionProxy<Base, Kernel> device_proxy;
    };
}

#endif //SPHINXSYS_DEVICE_EXECUTABLE_HPP
