#ifndef UPDATE_CELL_LINKED_LIST_HPP
#define UPDATE_CELL_LINKED_LIST_HPP

#include "update_cell_linked_list.h"

#include "adaptation.hpp"
#include "base_particles.hpp"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, typename CellLinkedListType>
UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::UpdateCellLinkedList(RealBody &real_body)
    : LocalDynamics(real_body), BaseDynamics<void>(),
      cell_linked_list_(DynamicCast<CellLinkedListType>(this, real_body.getCellLinkedList())),
      mesh_(cell_linked_list_),
      cell_offset_list_size_(cell_linked_list_.getCellOffsetListSize()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_particle_index_(cell_linked_list_.getParticleIndex()),
      dv_cell_offset_(cell_linked_list_.getCellOffset()),
      dv_current_cell_size_(DiscreteVariable<UnsignedInt>("CurrentCellSize", cell_offset_list_size_)),
      ex_policy_(ExecutionPolicy{}), kernel_implementation_(*this)
{
    particles_->addVariableToWrite<UnsignedInt>("ParticleIndex");
}
//=================================================================================================//
template <class ExecutionPolicy, typename CellLinkedListType>
UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType> &encloser)
    : mesh_(encloser.mesh_),
      cell_offset_list_size_(encloser.cell_offset_list_size_),
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      particle_index_(encloser.dv_particle_index_->DelegatedDataField(ex_policy)),
      cell_offset_(encloser.dv_cell_offset_->DelegatedDataField(ex_policy)),
      current_cell_size_(encloser.dv_current_cell_size_.DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename CellLinkedListType>
void UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::ComputingKernel::
    clearAllLists(UnsignedInt index_i)
{
    cell_offset_[index_i] = 0;
    current_cell_size_[index_i] = 0;
    particle_index_[index_i] = 0;
}
//=================================================================================================//
template <class ExecutionPolicy, typename CellLinkedListType>
void UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::ComputingKernel::
    incrementCellSize(UnsignedInt index_i)
{
    // Here, particle_index_ takes role of current_cell_size_list_.
    const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
    typename AtomicUnsignedIntRef<ExecutionPolicy>::type
        atomic_cell_size(particle_index_[linear_index]);
    ++atomic_cell_size;
}
//=================================================================================================//
template <class ExecutionPolicy, typename CellLinkedListType>
void UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::ComputingKernel::
    updateCellList(UnsignedInt index_i)
{
    // Here, particle_index_ takes its original role.
    const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
    typename AtomicUnsignedIntRef<ExecutionPolicy>::type
        atomic_current_cell_size(current_cell_size_[linear_index]);
    particle_index_[cell_offset_[linear_index] + atomic_current_cell_size++] = index_i;
}
//=================================================================================================//
template <class ExecutionPolicy, class CellLinkedListType>
void UpdateCellLinkedList<ExecutionPolicy, CellLinkedListType>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();

    particle_for(ex_policy_,
                 IndexRange(0, this->cell_offset_list_size_),
                 [=](size_t i)
                 { computing_kernel->clearAllLists(i); });

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementCellSize(i); });

    UnsignedInt *particle_index = this->dv_particle_index_->DelegatedDataField(ex_policy_);
    UnsignedInt *cell_offset = this->dv_cell_offset_->DelegatedDataField(ex_policy_);
    exclusive_scan(ex_policy_, particle_index, cell_offset,
                   this->cell_offset_list_size_,
                   typename PlusUnsignedInt<ExecutionPolicy>::type());

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateCellList(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
