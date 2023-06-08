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

            void writeBack() {}

        private:
            NoProxy<T> proxy;
        };

        template<class T>
        class ExecutionSelector<T, ParallelSYCLDevicePolicy> {
        public:
            explicit ExecutionSelector(T* localDynamics) : local_dynamics(localDynamics) {}

            auto& getProxy() {
                return local_dynamics->getDeviceProxy();
            }

            void writeBack() {
                local_dynamics->writeBack();
            }

        private:
            T* local_dynamics;
        };
    }

#endif //SPHINXSYS_EXECUTION_SELECTOR_HPP
