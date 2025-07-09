#ifndef MESH_ITERATOR_SYCL_HPP
#define MESH_ITERATOR_SYCL_HPP

#include "mesh_iterator.h"

namespace SPH
{
//=================================================================================================//
template <typename FunctionOnData>
void package_parallel_for(const ParallelDevicePolicy &par_device,
                          size_t num_grid_pkgs, const FunctionOnData &function)
{
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue.submit([&](sycl::handler &cgh)
    {
        cgh.parallel_for(execution_instance.getUniformNdRange(num_grid_pkgs),
                         [=](sycl::nd_item<1> index)
        {
            if(index.get_global_id(0) + 2 < num_grid_pkgs)
                function(index.get_global_id(0) + 2); 
        });
    }).wait_and_throw();
}
//=================================================================================================//
} // namespace SPH
#endif