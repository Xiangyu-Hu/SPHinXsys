#ifndef UPDATE_CELL_LINKED_LIST_HPP
#define UPDATE_CELL_LINKED_LIST_HPP

#include "update_cell_linked_list.h"

#include "adaptation.hpp"
#include "base_particles.hpp"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
template <class MeshType>
ParticleCellLinkedList<MeshType>::ParticleCellLinkedList(
    const MeshType &mesh, Vecd *pos,
    UnsignedInt *particle_id_list, UnsignedInt *particle_offset_list)
    : MeshType(mesh), pos_(pos), particle_id_list_(particle_id_list),
      particle_offset_list_(particle_offset_list) {}
//=================================================================================================//
template <class MeshType>
template <typename FunctionOnEach>
void ParticleCellLinkedList<MeshType>::
    forEachNeighbor(UnsignedInt index_i, const Vecd *source_pos,
                    const FunctionOnEach &function) const
{
    const Arrayi target_cell_index = this->CellIndexFromPosition(source_pos[index_i]);
    mesh_for_each(
        Arrayi::Zero().max(target_cell_index - Arrayi::Ones()),
        this->all_cells_.min(target_cell_index + 2 * Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            const UnsignedInt linear_index = this->transferMeshIndexTo1D(cell_index, this->all_cells_);
            // Since offset_cell_size_ has linear_cell_size_+1 elements, no boundary checks are needed.
            // offset_cell_size_[0] == 0 && offset_cell_size_[linear_cell_size_] == total_real_particles_
            for (UnsignedInt n = particle_offset_list_[linear_index]; n < particle_offset_list_[linear_index + 1]; ++n)
            {
                const UnsignedInt index_j = particle_id_list_[n];
                if ((source_pos[index_i] - pos_[index_j]).squaredNorm() < this->grid_spacing_)
                {
                    function(index_j);
                }
            }
        });
}
//=================================================================================================//
template <typename MeshType>
template <class ExecutionPolicy>
UpdateCellLinkedList<MeshType>::
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body)
    : LocalDynamics(real_body),
      mesh_(real_body.sph_adaptation_->createBackGroundMesh<MeshType>(real_body)),
      number_of_cells_plus_one_(mesh_.NumberOfCells() + 1),
      particle_id_list_size_(SMAX(particles_->ParticlesBound(), number_of_cells_plus_one_)),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleIDList", particle_id_list_size_)),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleOffsetList", number_of_cells_plus_one_)),
      current_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "CurrentCellSize", number_of_cells_plus_one_))
{
    particles_->addVariableToWrite<UnsignedInt>("ParticleIDList");
}
//=================================================================================================//
template <typename MeshType>
ParticleCellLinkedList<MeshType> UpdateCellLinkedList<MeshType>::getParticleCellLinkedList() const
{
    return ParticleCellLinkedList<MeshType>(mesh_, pos_, particle_id_list_, particle_offset_list_);
}
//=================================================================================================//
template <typename MeshType>
template <class T>
UpdateCellLinkedList<MeshType>::ComputingKernel<T>::
    ComputingKernel(UpdateCellLinkedList<MeshType> &update_cell_linked_list)
    : mesh_(update_cell_linked_list.mesh_),
      number_of_cells_plus_one_(update_cell_linked_list.number_of_cells_plus_one_),
      particle_id_list_size_(update_cell_linked_list.particle_id_list_size_),
      pos_(update_cell_linked_list.pos_),
      particle_id_list_(update_cell_linked_list.particle_id_list_),
      particle_offset_list_(update_cell_linked_list.particle_offset_list_),
      current_size_list_(update_cell_linked_list.current_size_list_) {}
//=================================================================================================//
template <typename MeshType>
template <class T>
void UpdateCellLinkedList<MeshType>::ComputingKernel<T>::clearAllLists(UnsignedInt index_i)
{
    particle_offset_list_[index_i] = 0;
    current_size_list_[index_i] = 0;
    particle_id_list_[index_i] = 0;
}
//=================================================================================================//
template <typename MeshType>
template <class T>
void UpdateCellLinkedList<MeshType>::ComputingKernel<T>::incrementCellSize(UnsignedInt index_i)
{
    // Here, particle_id_list_ takes role of current_cell_size_list_.
    const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
    typename AtomicUnsignedIntRef<T>::type atomic_cell_size(particle_id_list_[linear_index]);
    ++atomic_cell_size;
}
//=================================================================================================//
template <typename MeshType>
template <class T>
void UpdateCellLinkedList<MeshType>::ComputingKernel<T>::updateCellLists(UnsignedInt index_i)
{
    // Here, particle_id_list_ takes its original role.
    const UnsignedInt linear_index = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
    typename AtomicUnsignedIntRef<T>::type atomic_current_cell_size(current_size_list_[linear_index]);
    particle_id_list_[particle_offset_list_[linear_index] + atomic_current_cell_size++] = index_i;
}
//=================================================================================================//
template <class MeshType, class ExecutionPolicy>
UpdateCellLinkedList<MeshType, ExecutionPolicy>::UpdateCellLinkedList(RealBody &real_body)
    : UpdateCellLinkedList<MeshType>(ExecutionPolicy{}, real_body),
      BaseDynamics<void>(), kernel_implementation_(*this){};
//=================================================================================================//
template <class MeshType, class ExecutionPolicy>
void UpdateCellLinkedList<MeshType, ExecutionPolicy>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
    particle_for(ExecutionPolicy{},
                 IndexRange(0, this->number_of_cells_plus_one_),
                 [=](size_t i)
                 { computing_kernel->clearAllLists(i); });

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementCellSize(i); });

    exclusive_scan(ExecutionPolicy{}, this->particle_id_list_,
                   this->particle_id_list_ + this->number_of_cells_plus_one_,
                   this->particle_offset_list_,
                   typename PlusUnsignedInt<ExecutionPolicy>::type());

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateCellLists(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
