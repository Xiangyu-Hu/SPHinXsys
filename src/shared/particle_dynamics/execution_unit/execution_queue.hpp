#ifndef SPHINXSYS_EXECUTION_QUEUE_H
#define SPHINXSYS_EXECUTION_QUEUE_H

#include <sycl/sycl.hpp>

namespace SPH::execution {
    class ExecutionQueue {
    public:
        ExecutionQueue(ExecutionQueue const&) = delete;
        void operator=(ExecutionQueue const&) = delete;

        static ExecutionQueue& getInstance() {
            static ExecutionQueue instance;
            return instance;
        }

        sycl::queue &getQueue() {
            if(!sycl_queue)
                sycl_queue = std::make_unique<sycl::queue>(sycl::gpu_selector_v);
            return *sycl_queue;
        }

        auto getWorkGroupSize() const {
            return work_group_size;
        }

        void setWorkGroupSize(size_t workGroupSize) {
            work_group_size = workGroupSize;
        }

        static inline sycl::nd_range<1> getUniformNdRange(size_t global_size, size_t local_size) {
            return {global_size % local_size ? (global_size / local_size + 1) * local_size : global_size , local_size};
        }

        inline sycl::nd_range<1> getUniformNdRange(size_t global_size) {
            // sycl::nd_range is trivially-copyable, no std::move required
            return getUniformNdRange(global_size, work_group_size);
        }

    private:
        ExecutionQueue() : work_group_size(32), sycl_queue() {}

        std::size_t work_group_size;
        std::unique_ptr<sycl::queue> sycl_queue;

    } static &executionQueue = ExecutionQueue::getInstance();
}

#endif //SPHINXSYS_EXECUTION_QUEUE_H
