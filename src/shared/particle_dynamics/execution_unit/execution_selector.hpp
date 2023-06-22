#ifndef SPHINXSYS_EXECUTION_SELECTOR_HPP
#define SPHINXSYS_EXECUTION_SELECTOR_HPP

#include <sycl/sycl.hpp>
#include <utility>

#include "execution_proxy.hpp"

namespace SPH::execution {
        template<class T, class ExecutionPolicy>
        class ExecutionSelector {
        public:
            explicit ExecutionSelector(T* localDynamics) : proxy(localDynamics) {}

            auto& getProxy() {
                return proxy;
            }

        private:
            NoProxy<T> proxy;
        };

        template<class T>
        class ExecutionSelector<T, ParallelSYCLDevicePolicy> {
        public:
            explicit ExecutionSelector(T* localDynamics) : local_dynamics(localDynamics),
                buffer_local_dynamics_kernel(sycl::buffer<KernelT>(local_dynamics->getDeviceProxy().getKernel(), 1)) {}

            auto& getProxy() {
                return local_dynamics->getDeviceProxy();
            }

            auto& getBuffer() {
                return buffer_local_dynamics_kernel;
            }

        private:
            using KernelT =
                    std::remove_pointer_t<decltype(std::declval<T>().getDeviceProxy().getKernel())>;
            T* local_dynamics;
            sycl::buffer<KernelT, 1> buffer_local_dynamics_kernel;
        };
    }

#endif //SPHINXSYS_EXECUTION_SELECTOR_HPP
