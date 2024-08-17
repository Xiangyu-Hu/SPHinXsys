#ifndef UPDATE_BODY_RELATION_HPP
#define UPDATE_BODY_RELATION_HPP

#include "update_body_relation.h"

namespace SPH
{
//=================================================================================================//
template <class ParticleCellLinkedListType>
template <class ExecutionPolicy>
Relation<Inner<ParticleCellLinkedListType>>::Relation(
    const ExecutionPolicy &ex_policy, RealBody &real_body,
    const ParticleCellLinkedListType &particle_cell_linked_list)
    : LocalDynamics(real_body),
      particle_cell_linked_list_(particle_cell_linked_list),
      real_particle_bound_plus_one_(particles_->RealParticlesBound() + 1),
      neighbor_id_list_size_(particles_->ParticlesBound() * NeighborSize<Dimensions>::value),
      pos_(particles_->getVariableDataByName<Vecd>(ex_policy, "Position")),
      neighbor_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(ex_policy, "NeighborIDList", neighbor_id_list_size_)),
      neighbor_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(ex_policy, "NeighborOffsetList", real_particle_bound_plus_one_)),
      neighbor_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(ex_policy, "NeighborSizeList", real_particle_bound_plus_one_))
{
    particles_->addVariableToWrite<UnsignedInt>("NeighborSizeList");
}
//=================================================================================================//
template <class ParticleCellLinkedListType>
template <class T>
Relation<Inner<ParticleCellLinkedListType>>::ComputingKernel<T>::
    ComputingKernel(Relation<Inner<ParticleCellLinkedListType>> &update_inner_relation)
    : particle_cell_linked_list_(update_inner_relation.particle_cell_linked_list_),
      neighborhood_(update_inner_relation.neighborhood_),
      real_particle_bound_plus_one_(update_inner_relation.real_particle_bound_plus_one_),
      neighbor_id_list_size_(update_inner_relation.neighbor_id_list_size_),
      pos_(update_inner_relation.pos_),
      neighbor_id_list_(update_inner_relation.neighbor_id_list_),
      neighbor_offset_list_(update_inner_relation.neighbor_offset_list_),
      neighbor_size_list_(update_inner_relation.neighbor_size_list_) {}
//=================================================================================================//
template <class ParticleCellLinkedListType>
template <class T>
void Relation<Inner<ParticleCellLinkedListType>>::ComputingKernel<T>::
    clearAllLists(UnsignedInt index_i)
{
    neighbor_offset_list_[index_i] = 0;
    neighbor_size_list_[index_i] = 0;
    neighbor_id_list_[index_i] = 0;
}
//=================================================================================================//
template <class ParticleCellLinkedListType>
template <class T>
void Relation<Inner<ParticleCellLinkedListType>>::ComputingKernel<T>::
    incrementNeighborSize(UnsignedInt index_i)
{
    // Here, neighbor_id_list_ takes role of neighbor_size_list_.
    particle_cell_linked_list_.forEachNeighbor(index_i, pos_, neighborhood_,
                                               [=](size_t j)
                                               { ++neighbor_id_list_[index_i]; });
}
//=================================================================================================//
template <class ParticleCellLinkedListType>
template <class T>
void Relation<Inner<ParticleCellLinkedListType>>::ComputingKernel<T>::
    updateNeighborIDList(UnsignedInt index_i)
{
    particle_cell_linked_list_.forEachNeighbor(
        index_i, pos_, neighborhood_,
        [=](size_t j)
        {
            ++neighbor_size_list_[index_i];
            neighbor_id_list_[neighbor_offset_list_[index_i] + neighbor_size_list_[index_i]] = j;
        });
}
//=================================================================================================//
template <typename... T, class ExecutionPolicy>
Relation<ExecutionPolicy, T...>::Relation(RealBody &real_body, T &&...parameters)
    : Relation<T...>(real_body, std::forward<T>(parameters)...),
      BaseDynamics<void>(), kernel_implementation_(*this) {}
//=================================================================================================//
template <typename... T, class ExecutionPolicy>
void Relation<ExecutionPolicy, T...>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
    particle_for(ExecutionPolicy{},
                 IndexRange(0, this->real_particle_bound_plus_one_),
                 [=](size_t i)
                 { computing_kernel->clearAllLists(i); });

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementNeighborSize(i); });

    exclusive_scan(ExecutionPolicy{}, this->neighbor_id_list_,
                   this->neighbor_id_list_ + this->real_particle_bound_plus_one_,
                   this->neighbor_offset_list_,
                   typename PlusUnsignedInt<ExecutionPolicy>::type());

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborIDList(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_BODY_RELATION_HPP
