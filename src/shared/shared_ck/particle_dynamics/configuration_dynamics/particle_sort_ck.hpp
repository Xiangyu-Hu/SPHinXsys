#ifndef PARTICLE_SORT_HPP
#define PARTICLE_SORT_HPP

#include "particle_sort_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
ParticleSortCK<ExecutionPolicy>::ParticleSortCK(RealBody &real_body)
    : LocalDynamics(real_body), BaseDynamics<void>(),
      ex_policy_(ExecutionPolicy{}),
      cell_linked_list_(real_body.getCellLinkedList()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_sequence_(particles_->registerDiscreteVariable<UnsignedInt>(
          "Sequence", particles_->ParticlesBound())),
      dv_index_permutation_(particles_->registerDiscreteVariable<UnsignedInt>(
          "IndexPermutation", particles_->ParticlesBound())),
      sort_method_(ExecutionPolicy{}, dv_sequence_, dv_index_permutation_),
      kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy>
ParticleSortCK<ExecutionPolicy>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, ParticleSortCK<ExecutionPolicy> &encloser)
    : mesh_(encloser.cell_linked_list_.getSortSequenceMesh()),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      sequence_(encloser.dv_sequence_->DelegatedData(ex_policy)),
      index_permutation_(encloser.dv_index_permutation_->DelegatedData(ex_policy)){}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleSortCK<ExecutionPolicy>::ComputingKernel::
    prepareSequence(UnsignedInt index_i)
{
    sequence_[index_i] = Mesh::transferMeshIndexToMortonOrder(mesh_.CellIndexFromPosition(pos_[index_i]));
    index_permutation_[index_i] = index_i;
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleSortCK<ExecutionPolicy>::exec(Real dt)
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();

    particle_for(ex_policy_, IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->prepareSequence(i); });

    sort_method_.sort(ex_policy_, total_real_particles);
    update_variables_to_sort_(particles_->EvolvingVariables(), ex_policy_,
                              0, total_real_particles, dv_index_permutation_);
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORT_HPP