#ifndef UPDATE_CELL_LINKED_LIST_SYCL_HPP
#define UPDATE_CELL_LINKED_LIST_SYCL_HPP

#include "update_cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
template <typename MeshType>
UpdateCellLinkedList<MeshType, ParallelDevicePolicy>::UpdateCellLinkedList(RealBody &real_body)
    : LocalDynamics(real_body), BaseDynamics<void>(real_body),
      mesh_(DynamicsCast<CellLinkedList>(this, &real_body.getCellLinkedList())),
      particles_bound_(particles_->ge)
{
    particles_bound_ = mesh_.NumberOfCells();
    pos_ = particles_->getVariablePositions();
    particle_id_list_ = particles_->getVariableParticleIndex();
    particle_offset_list_ = particles_->getVariableParticleOffset();
    v_total_real_particles_ = particles_->getVariable<indexScalar, Real>("TotalRealParticles");
}
{
    if (existDeviceDataField())
    {
        copyFromDevice(data_field_, device_data_field_, size_);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_HPP
