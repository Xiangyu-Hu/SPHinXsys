#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_QUEUE_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_QUEUE_HPP

/** Thin SYCL execution helpers for matrix-free A–phi operators (extra_src). */

#if SPHINXSYS_USE_SYCL
#include "implementation_sycl.h"
#include <cstdlib>
#include <sycl/sycl.hpp>
#include <utility>

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

/** Submit one command group on the global execution queue and wait for completion (host-visible state). */
template <typename Fn>
inline void matrixFreeSyclSubmitAndWait(Fn &&fn)
{
    sycl::queue &q = execution::ExecutionInstance::getInstance().getQueue();
    q.submit(std::forward<Fn>(fn));
    q.wait();
}

/** Submit without waiting. Chained submits on the same default in-order queue complete in submission order; call
 *  `matrixFreeSyclQueueWait` (or `matrixFreeSyclSubmitAndWait`) before host accessors that read the results. */
template <typename Fn>
inline void matrixFreeSyclSubmit(Fn &&fn)
{
    sycl::queue &q = execution::ExecutionInstance::getInstance().getQueue();
    q.submit(std::forward<Fn>(fn));
}

inline void matrixFreeSyclQueueWait()
{
    execution::ExecutionInstance::getInstance().getQueue().wait();
}

/** When `EM_APHI_MATRIX_FREE_SYCL_FIELD_VALUES_USM=1`, scalar complex fields may use device USM mirrors (see scalar/Jacobi helpers). */
inline bool matrixFreeSyclFieldValuesUseDeviceUsm()
{
    static const bool use = []() {
        const char *e = std::getenv("EM_APHI_MATRIX_FREE_SYCL_FIELD_VALUES_USM");
        return e != nullptr && e[0] == '1' && e[1] == '\0';
    }();
    return use;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // SPHINXSYS_USE_SYCL

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_QUEUE_HPP
