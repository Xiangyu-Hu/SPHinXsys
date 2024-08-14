#ifndef UPDATE_CELL_LINKED_LIST_SYCL_HPP
#define UPDATE_CELL_LINKED_LIST_SYCL_HPP

#include "update_cell_linked_list_sycl.h"

namespace SPH
{
//=================================================================================================//
template <typename MeshType>
void UpdateCellLinkedList<MeshType>::
    clearParticleOffsetList(const ParallelDevicePolicy &par_device)
{
    execution_instance.getQueue()
        .submit(
            [=, number_of_cells_plus_one = number_of_cells_plus_one_,
             particle_id_list = particle_id_list_,
             particle_offset_list = particle_offset_list_,
             current_size_list = current_size_list_](sycl::handler &cgh)
            {
                cgh.parallel_for(execution_instance.getUniformNdRange(number_of_cells_plus_one),
                                 [=](sycl::nd_item<1> item)
                                 {
                                     if (item.get_global_id(0) < number_of_cells_plus_one)
                                     {
                                         UnsignedInt linear_index = item.get_global_id(0);
                                         particle_offset_list[linear_index] = 0;
                                         current_size_list[linear_index] = 0;
                                         particle_id_list[linear_index] = 0;
                                     }
                                 });
            })
        .wait_and_throw();
}
//=================================================================================================//
template <typename MeshType>
void UpdateCellLinkedList<MeshType>::
    incrementCellSize(const ParallelDevicePolicy &par_device)
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    execution_instance.getQueue()
        .submit( // note that particle_id_list_ is used here as temporary storage
            [=, mesh = mesh_, pos = pos_, particle_id_list = particle_offset_list_](sycl::handler &cgh)
            {
                cgh.parallel_for(execution_instance.getUniformNdRange(total_real_particles),
                                 [=](sycl::nd_item<1> item)
                                 {
                                     if (item.get_global_id(0) < total_real_particles)
                                     {
                                         UnsignedInt particle_i = item.get_global_id(0);
                                         UnsignedInt linear_index =
                                             mesh->LinearCellIndexFromPosition(pos[particle_i]);
                                         AtomicUnsignedIntRef<ParallelDevicePolicy>::type
                                             atomic_cell_size(particle_id_list[linear_index]);
                                         ++atomic_cell_size;
                                     }
                                 });
            })
        .wait_and_throw();
}
//=================================================================================================//
template <typename MeshType>
void UpdateCellLinkedList<MeshType>::
    exclusiveScanParticleOffsetList(const ParallelDevicePolicy &par_device)
{
    execution_instance.getQueue()
        .submit( // note that particle_id_list_ is still used here as temporary storage
            [=, particle_id_list = particle_id_list_, particle_offset_list = particle_offset_list_,
             number_of_cells_plus_one = number_of_cells_plus_one_](sycl::handler &cgh)
            {
                cgh.parallel_for(execution_instance.getUniformNdRange(execution_instance.getWorkGroupSize()),
                                 [=](sycl::nd_item<1> item)
                                 {
                                     if (item.get_group_linear_id() == 0)
                                     {
                                         sycl::joint_exclusive_scan(
                                             item.get_group(), particle_id_list,
                                             particle_id_list + number_of_cells_plus_one,
                                             particle_offset_list, sycl::plus<>{});
                                     }
                                 });
            })
        .wait_and_throw();
}
//=================================================================================================//
template <typename MeshType>
void UpdateCellLinkedList<MeshType>::
    updateCellLists(const ParallelDevicePolicy &par_device)
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    execution_instance.getQueue()
        .submit( // note that particle_id_list_ is finally used for its intended purpose
            [=, mesh = mesh_, pos = pos_, current_size_list = current_size_list_,
             particle_id_list = particle_id_list_,
             particle_offset_list = particle_offset_list_](sycl::handler &cgh)
            {
                cgh.parallel_for(execution_instance.getUniformNdRange(total_real_particles),
                                 [=](sycl::nd_item<1> item)
                                 {
                                     if (item.get_global_id(0) < total_real_particles)
                                     {

                                         UnsignedInt particle_i = item.get_global_id(0);
                                         UnsignedInt linear_index =
                                             mesh->LinearCellIndexFromPosition(pos[particle_i]);
                                         AtomicUnsignedIntRef<ParallelDevicePolicy>::type
                                             atomic_current_cell_size(current_size_list[linear_index]);
                                         particle_id_list[particle_offset_list[linear_index] + atomic_current_cell_size++] = particle_i;
                                     }
                                 });
            })
        .wait_and_throw();
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_HPP
