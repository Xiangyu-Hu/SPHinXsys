#ifndef SPHINXSYS_EXECUTION_ARGUMENT_HPP
#define SPHINXSYS_EXECUTION_ARGUMENT_HPP

#include <sycl/sycl.hpp>
#include <utility>

#include "execution_proxy.hpp"

namespace SPH {
    namespace execution {
        template<class T, sycl::access_mode access_mode, typename=std::true_type>
        class DeviceVariable {
            using VarType = std::remove_reference_t<T>;
        public:
            DeviceVariable(VarType* var_addr, std::size_t var_size) : var_buffer(var_addr, var_size){}

            auto init_memory_access(sycl::handler& cgh) {
                return var_buffer.template get_access<access_mode>(cgh);
            }

            auto get_host_memory_access() {
                return var_buffer.get_host_access();
            }

        private:
            sycl::buffer<VarType , 1> var_buffer;
        };

        template<class Proxy, typename=std::enable_if_t<std::is_base_of_v<ExecutionProxy<typename Proxy::Base,
                typename Proxy::Kernel>, Proxy>>>
        class DeviceProxyVariable {
        public:
            explicit DeviceProxyVariable(Proxy &proxy) : var_proxy(proxy) {}

            auto init_memory_access(sycl::handler& cgh) {
                var_proxy.get_memory_access(Context<ParallelSYCLDevicePolicy>(cgh), par_sycl);
                return *var_proxy.get(par_sycl);
            }
        private:
            Proxy& var_proxy;
        };

        template<typename BaseT, typename KernelT, typename ...DeviceVariables>
        class DeviceProxy : public ExecutionProxy<BaseT, KernelT> {
        public:
            explicit DeviceProxy(BaseT* base, DeviceVariables& ...variables)
                : ExecutionProxy<BaseT, KernelT>(base, new KernelT()), deviceVariables(variables...) {}

            ~DeviceProxy() {
                delete this->kernel;
            }

        protected:

            void get_memory_access_device(sycl::handler &cgh) override {
                auto build_accessors_argument = [&](DeviceVariables& ...vars) {
                    return std::make_tuple(vars.init_memory_access(cgh)...);
                };
                this->kernel->setAccessors(std::apply(build_accessors_argument, deviceVariables));
            }

        private:
            std::tuple<DeviceVariables&...> deviceVariables;
        };

        template<class T>
        class DeviceDispatcher {
        public:
            explicit DeviceDispatcher(T* localDynamics) : local_dynamics(localDynamics) {}

            auto& getProxy() {
                return local_dynamics->getDeviceProxy();
            }

            void writeBack() {
                local_dynamics->writeBack();
            }

        private:
            T* local_dynamics;
        };

        template<class T>
        class NoDispatcher {
        public:
            explicit NoDispatcher(T* localDynamics) : proxy(localDynamics) {}

            auto& getProxy() {
                return proxy;
            }

            void writeBack() {}

        private:
            NoProxy<T> proxy;
        };
    }
}

#endif //SPHINXSYS_EXECUTION_ARGUMENT_HPP
