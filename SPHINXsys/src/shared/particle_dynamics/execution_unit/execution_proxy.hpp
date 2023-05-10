#ifndef SPHINXSYS_EXECUTION_PROXY_HPP
#define SPHINXSYS_EXECUTION_PROXY_HPP

#include "execution_policy.h"
#include "execution_context.hpp"

namespace SPH {
    namespace execution {
        template<typename BaseT, typename KernelT>
        class ExecutionProxy {
        public:
            using Base = BaseT;
            using Kernel = KernelT;

            ExecutionProxy(BaseT* base, KernelT* proxy) : base(base), kernel(proxy) {}

            template<class ExecutionPolicy = ParallelPolicy>
            BaseT* get(const ExecutionPolicy& = par,
                       std::enable_if_t<std::negation_v<std::is_same<ExecutionPolicy, ParallelSYCLDevicePolicy>>>* = nullptr) {
                return base;
            }

            KernelT* get(const ParallelSYCLDevicePolicy&) {
                return kernel;
            }

            template<class ExecutionPolicy>
            void copy_memory(const ExecutionPolicy&) {
                if constexpr (std::is_same_v<ExecutionPolicy, ParallelSYCLDevicePolicy>)
                    this->copy_memory_to_device();
            }

            template<class ExecutionPolicy>
            void copy_back(const ExecutionPolicy&) {
                if constexpr (std::is_same_v<ExecutionPolicy, ParallelSYCLDevicePolicy>)
                    this->copy_back_from_device();
            }

        protected:
            virtual void copy_memory_to_device() = 0;
            virtual void copy_back_from_device() = 0;

            BaseT* base;
            KernelT* kernel;
        };


        template<typename T>
        class NoProxy : public ExecutionProxy<T, T> {
        public:
            explicit NoProxy(T *base) : ExecutionProxy<T, T>(base, base) {}

            void get_memory_access(sycl::handler& cgh) {
                static_assert("No device memory access available for this class.");
            }

        protected:
            void copy_memory_to_device() override {
                static_assert("No device copy available for this class.");
            }

            void copy_back_from_device() override {
                static_assert("No device copy-back available for this class.");
            }
        };
    }
}

#endif //SPHINXSYS_EXECUTION_PROXY_HPP
