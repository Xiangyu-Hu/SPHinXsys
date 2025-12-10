#ifndef MESH_ITERATOR_SYCL_HPP
#define MESH_ITERATOR_SYCL_HPP

#include "implementation_sycl.h"
#include "mesh_iterators.h"

namespace SPH
{
//=================================================================================================//
template <typename FunctionOnData>
void package_for(const ParallelDevicePolicy &par_device, UnsignedInt start_index,
                 UnsignedInt num_grid_pkgs, const FunctionOnData &function)
{
    UnsignedInt operations = num_grid_pkgs - start_index;
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue.submit([&](sycl::handler &cgh)
                      { cgh.parallel_for(execution_instance.getUniformNdRange(operations),
                                         [=](sycl::nd_item<1> index)
                                         {
                                             if (index.get_global_id(0) < operations)
                                                 function(index.get_global_id(0) + start_index);
                                         }); })
        .wait_and_throw();
}
//=================================================================================================//
} // namespace SPH
#endif