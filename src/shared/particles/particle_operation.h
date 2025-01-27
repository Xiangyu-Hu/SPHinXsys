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

class CreateRealParticleFrom
{
    ParticleVariables &variables_to_sort_;
    DiscreteVariableArrays copyable_states_;
    DiscreteVariable<UnsignedInt> *dv_original_id_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    UnsignedInt real_particles_bound_;

  public:
    CreateRealParticleFrom(BaseParticles *particles);

    class ComputingKernel // only run with sequenced policy for now
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        UnsignedInt operator()(UnsignedInt index_i)
        {
            UnsignedInt new_original_id = *total_real_particles_;
            original_id_[new_original_id] = new_original_id;
            /** Buffer Particle state copied from real particle. */
            copy_particle_state_(copyable_state_data_arrays_, new_original_id, index_i);
            /** Realize the buffer particle by increasing the number of real particle by one.  */
            *total_real_particles_ += 1;
            return new_original_id;
        };

      protected:
        UnsignedInt *total_real_particles_;
        UnsignedInt real_particles_bound_;
        UnsignedInt *original_id_;
        VariableDataArrays copyable_state_data_arrays_;
        OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK> copy_particle_state_;
    };
};
} // namespace SPH
#endif // PARTICLE_OPERATION_H
