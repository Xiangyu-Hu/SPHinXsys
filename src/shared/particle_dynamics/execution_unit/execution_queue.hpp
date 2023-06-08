#ifndef SPHINXSYS_EXECUTION_QUEUE_H
#define SPHINXSYS_EXECUTION_QUEUE_H

#include <sycl/sycl.hpp>

namespace SPH::execution {
    class ExecutionQueue {
    public:
        static sycl::queue &getQueue() {
            static std::unique_ptr<sycl::queue> sycl_queue;
            if (sycl_queue == nullptr) {
                try {
                    static sycl::async_handler error_handler = [](const auto &list_errors) {
                        for (const auto &error: list_errors)
                            std::rethrow_exception(error);
                    };

                    sycl_queue = std::make_unique<sycl::queue>(sycl::gpu_selector_v, error_handler);

                } catch (const sycl::exception &error) {
                    std::cerr << error.what() << std::endl;
                }
            }
            return *sycl_queue;
        }
    };
}

#endif //SPHINXSYS_EXECUTION_QUEUE_H
