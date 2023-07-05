#ifndef SPHINXSYS_EXECUTION_PROXY_HPP
#define SPHINXSYS_EXECUTION_PROXY_HPP

#include "execution_policy.h"

namespace SPH {
    namespace execution {
        template<typename BaseT, typename KernelT>
        class ExecutionProxy {
        public:
            using Base = BaseT;
            using Kernel = KernelT;

            template<class ...Args>
            ExecutionProxy(BaseT* base, Args&&... kernelArgs) : base(base), kernel(new KernelT(std::forward<Args>(kernelArgs)...)) {}

            ~ExecutionProxy() {
                delete kernel;
            }

            template<class ExecutionPolicy = ParallelPolicy>
            BaseT* get(const ExecutionPolicy& = par) const {
                return base;
            }

            KernelT* get(const ParallelSYCLDevicePolicy&) const {
                return kernel;
            }

            BaseT *getBase() const {
                return base;
            }

            KernelT *getKernel() const {
                return kernel;
            }

        protected:
            BaseT* base;
            KernelT* kernel;
        };


        template<typename T>
        class NoProxy {
        public:
            explicit NoProxy(T *base) : base_(base) {}

            template<class ExecutionPolicy = ParallelPolicy>
            T* get(const ExecutionPolicy& = par) const {
                return base_;
            }
        private:
            T* base_;
        };
    }
}

#endif //SPHINXSYS_EXECUTION_PROXY_HPP
