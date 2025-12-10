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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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
    void operator()(VariableAllocationSet<AllocatedDataArray<DataType>> &variable_allocation_pair,
                    size_t index, size_t another_index);
};

class SpawnRealParticle
{
    ParticleVariables &evolving_variables_;
    DiscreteVariableArrayAssemble copyable_states_;
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
            UnsignedInt last_real_particle_index = total_real_particles_ref.fetch_add(1);
            if (last_real_particle_index < particles_bound_)
            {
                copy_particle_state_(copyable_state_data_arrays_, last_real_particle_index, index_i);
                original_id_[last_real_particle_index] = last_real_particle_index; // reinitialize original id
            }
            return last_real_particle_index;
        };

      protected:
        UnsignedInt *total_real_particles_;
        UnsignedInt particles_bound_;
        UnsignedInt *original_id_;
        VariableDataArrayAssemble copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrayAssemble, CopyParticleStateCK> copy_particle_state_;
    };
};

class RemoveRealParticle
{
    ParticleVariables &evolving_variables_;
    DiscreteVariableArrayAssemble copyable_states_;
    DiscreteVariable<UnsignedInt> *dv_original_id_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;

  public:
    RemoveRealParticle(BaseParticles *particles);

    class ComputingKernel // only run with sequenced policy for now
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void operator()(UnsignedInt index_i, int *life_status)
        {
            AtomicRef<UnsignedInt> total_real_particles_ref(*total_real_particles_);
            UnsignedInt last_real_particle_index = total_real_particles_ref.fetch_sub(1) - 1;
            while (life_status[last_real_particle_index] == 1) // to delete
            {
                life_status[last_real_particle_index] = 0; // reset the life status
                last_real_particle_index = total_real_particles_ref.fetch_sub(1) - 1;
            }

            if (index_i < last_real_particle_index)
            {
                UnsignedInt old_original_id = original_id_[index_i];
                copy_particle_state_(copyable_state_data_arrays_, index_i, last_real_particle_index);
                life_status[index_i] = 0;                                 // reset the life status
                original_id_[last_real_particle_index] = old_original_id; // swap the original id
            }
        };

      protected:
        UnsignedInt *total_real_particles_;
        UnsignedInt *original_id_;
        VariableDataArrayAssemble copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrayAssemble, CopyParticleStateCK> copy_particle_state_;
    };
};
} // namespace SPH
#endif // PARTICLE_OPERATION_H
