#ifndef UPDATE_CELL_LINKED_LIST_SYCL_HPP
#define UPDATE_CELL_LINKED_LIST_SYCL_HPP

#include "update_cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType>::
    setParticleOffsetListUpperBound(const ParallelDevicePolicy &par_device)
{
    UnsignedInt *host_data_field = setParticleOffsetListUpperBound();
    copyToDevice(host_data_field + number_of_cells_, particle_offset_list_ + number_of_cells_, 1);
}
//=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType>::
    exclusiveScanParticleOffsetList(const ParallelDevicePolicy &par_device)
{
    execution_instance.getQueue()
        .submit(
            [=](sycl::handler &cgh)
            {
                cgh.parallel_for(executionQueue.getUniformNdRange(executionQueue.getWorkGroupSize()),
                                 [=](sycl::nd_item<1> item)
                                 {
                                     if (item.get_group_linear_id() == 0)
                                     {
                                         sycl::joint_exclusive_scan(item.get_group(), particle_offset_list_,
                                                                    particle_offset_list_ + number_of_cells_ + 1,
                                                                    particle_offset_list_, sycl::plus<>{});
                                     }
                                 });
            })
        .wait_and_throw();
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_HPP
