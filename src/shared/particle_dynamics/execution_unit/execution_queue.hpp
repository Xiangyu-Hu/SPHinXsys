#ifndef SPHINXSYS_EXECUTION_QUEUE_H
#define SPHINXSYS_EXECUTION_QUEUE_H

#include <sycl/sycl.hpp>

namespace SPH::execution {
    class ExecutionQueue {
    public:
        ExecutionQueue() : work_group_size(32), sycl_queue(sycl::gpu_selector_v) {}

        sycl::queue &getQueue() {
            return sycl_queue;
        }

        auto getWorkGroupSize() const {
            return work_group_size;
        }

        void setWorkGroupSize(size_t workGroupSize) {
            work_group_size = workGroupSize;
        }

    private:
        std::size_t work_group_size;
        sycl::queue sycl_queue;
    } static executionQueue;
}

#endif //SPHINXSYS_EXECUTION_QUEUE_H
