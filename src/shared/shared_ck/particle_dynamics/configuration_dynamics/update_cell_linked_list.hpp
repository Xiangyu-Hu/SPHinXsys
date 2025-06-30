#ifndef UPDATE_CELL_LINKED_LIST_HPP
#define UPDATE_CELL_LINKED_LIST_HPP

#include "update_cell_linked_list.h"

#include "adaptation.hpp"
#include "base_particles.hpp"
#include "mesh_iterators.hpp"
#include "particle_iterators_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::
    UpdateCellLinkedList(DynamicsIdentifier &identifier)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier), BaseDynamics<void>(),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, identifier.getCellLinkedList())),
      mesh_(cell_linked_list_.getMesh()),
      cell_offset_list_size_(cell_linked_list_.getCellOffsetListSize()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_particle_index_(cell_linked_list_.dvParticleIndex()),
      dv_cell_offset_(cell_linked_list_.dvCellOffset()),
      dv_current_cell_size_(DiscreteVariable<UnsignedInt>("CurrentCellSize", cell_offset_list_size_)),
      identifier_(identifier), kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mesh_(encloser.mesh_), particle_mask_(ex_policy, encloser.identifier_),
      cell_offset_list_size_(encloser.cell_offset_list_size_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      particle_index_(encloser.dv_particle_index_->DelegatedData(ex_policy)),
      cell_offset_(encloser.dv_cell_offset_->DelegatedData(ex_policy)),
      current_cell_size_(encloser.dv_current_cell_size_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
void UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::ComputingKernel::
    clearAllLists(UnsignedInt index_i)
{
    cell_offset_[index_i] = 0;
    current_cell_size_[index_i] = 0;
    particle_index_[index_i] = 0;
}
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
void UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::ComputingKernel::
    incrementCellSize(UnsignedInt index_i)
{
    if (particle_mask_(index_i))
    {
        // Here, particle_index_ takes role of current_cell_size_list_.
        const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
        AtomicRef<UnsignedInt> atomic_cell_size(particle_index_[linear_index]);
        ++atomic_cell_size;
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
void UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::ComputingKernel::
    updateCellList(UnsignedInt index_i)
{
    if (particle_mask_(index_i))
    {
        // Here, particle_index_ takes its original role.
        const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
        AtomicRef<UnsignedInt> atomic_current_cell_size(current_cell_size_[linear_index]);
        particle_index_[cell_offset_[linear_index] + atomic_current_cell_size++] = index_i;
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename DynamicsIdentifier>
void UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();

    particle_for(ExecutionPolicy{},
                 IndexRange(0, this->cell_offset_list_size_),
                 [=](size_t i)
                 { computing_kernel->clearAllLists(i); });

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementCellSize(i); });

    UnsignedInt *particle_index = this->dv_particle_index_->DelegatedData(ExecutionPolicy{});
    UnsignedInt *cell_offset = this->dv_cell_offset_->DelegatedData(ExecutionPolicy{});
    exclusive_scan(ExecutionPolicy{}, particle_index, cell_offset,
                   this->cell_offset_list_size_,
                   typename PlusUnsignedInt<ExecutionPolicy>::type());

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateCellList(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
