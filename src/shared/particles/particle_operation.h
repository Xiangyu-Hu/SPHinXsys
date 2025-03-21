/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file particle_operation.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef PARTICLE_OPERATION_H
#define PARTICLE_OPERATION_H

#include "base_particles.hpp"
#include "sphinxsys_variable_array.h"

namespace SPH
{

struct CopyParticleStateCK
{
    template <typename DataType>
    void operator()(VariableAllocationPair<AllocatedDataArray<DataType>> &variable_allocation_pair,
                    size_t index, size_t another_index);
};

class SpawnRealParticle
{
    ParticleVariables &evolving_variables_;
    DiscreteVariableArrays copyable_states_;
    DiscreteVariable<UnsignedInt> *dv_original_id_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    UnsignedInt particles_bound_;

  public:
    SpawnRealParticle(BaseParticles *particles);

    class ComputingKernel // only run with sequenced policy for now
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        UnsignedInt operator()(UnsignedInt index_i)
        {
            AtomicRef<UnsignedInt> total_real_particles_ref(*total_real_particles_);
            UnsignedInt new_original_id = total_real_particles_ref.fetch_add(1);
            if (new_original_id < particles_bound_)
            {
                copy_particle_state_(copyable_state_data_arrays_, new_original_id, index_i);
                original_id_[new_original_id] = new_original_id;
            }
            else
            {
                return particles_bound_;
            }
            return new_original_id;
        };

      protected:
        UnsignedInt *total_real_particles_;
        UnsignedInt particles_bound_;
        UnsignedInt *original_id_;
        VariableDataArrays copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK> copy_particle_state_;
    };
};

class DespawnRealParticle
{
    ParticleVariables &evolving_variables_;
    DiscreteVariableArrays copyable_states_;
    DiscreteVariable<UnsignedInt> *dv_original_id_, *dv_sorted_id_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    UnsignedInt real_particles_bound_;

  public:
    DespawnRealParticle(BaseParticles *particles);

    class ComputingKernel // only run with sequenced policy for now
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        template <class IsDeletable>
        void operator()(UnsignedInt index_i, const IsDeletable &is_deletable)
        {
            AtomicRef<UnsignedInt> total_real_particles_ref(*total_real_particles_);
            UnsignedInt last_real_particle_index = total_real_particles_ref.fetch_sub(1) - 1;
            while (is_deletable(last_real_particle_index))
            {
                last_real_particle_index = total_real_particles_ref.fetch_sub(1) - 1;
            }

            if (index_i < last_real_particle_index)
            {
                copy_particle_state_(copyable_state_data_arrays_, index_i, last_real_particle_index);
                original_id_[index_i] = original_id_[last_real_particle_index];
                sorted_id_[original_id_[index_i]] = index_i;
            }
        };

      protected:
        UnsignedInt *total_real_particles_;
        UnsignedInt real_particles_bound_;
        UnsignedInt *original_id_, *sorted_id_;
        VariableDataArrays copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK> copy_particle_state_;
    };
};
} // namespace SPH
#endif // PARTICLE_OPERATION_H
