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
            UnsignedInt new_original_id = *total_real_particles_;
            if (new_original_id < particles_bound_)
            {
                /** Buffer Particle state copied from real particle. */
                copy_particle_state_(copyable_state_data_arrays_, new_original_id, index_i);
                /** Realize the buffer particle by increasing the number of real particle by one.  */
                *total_real_particles_ += 1;
                original_id_[new_original_id] = new_original_id;
            }
            // Use an unordered_map to count how many times each ID appears
            std::unordered_map<unsigned int, size_t> frequency;

            // Pass through the array once, counting occurrences
            for (size_t i = 0; i < *total_real_particles_; ++i)
            {
                frequency[original_id_[i]]++;
            }

            // Check which IDs appeared more than once
            for (const auto &entry : frequency)
            {
                if (entry.second > 1)
                {
                    std::cout << "Repeated ID: " << entry.first
                              << " occurs " << entry.second << " times.\n";
                }
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
    DiscreteVariable<UnsignedInt> *dv_original_id_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    UnsignedInt real_particles_bound_;

  public:
    DespawnRealParticle(BaseParticles *particles);

    class ComputingKernel // only run with sequenced policy for now
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        UnsignedInt operator()(UnsignedInt index_i)
        {
            UnsignedInt last_real_particle_index = *total_real_particles_ - 1;
            if (index_i < last_real_particle_index)
            {
                auto tid = std::this_thread::get_id();
                std::cout << "Thread ID: " << tid << std::endl;
                /** Buffer Particle state copied from real particle. */
                const UnsignedInt temp_original_id_index_i = original_id_[index_i];
                std::cout << "l109 :<<" << temp_original_id_index_i << ", original_id_[index_i]:" << original_id_[index_i] << ",original_id_[last_real_particle_index]:" << original_id_[last_real_particle_index] << std::endl;
                copy_particle_state_(copyable_state_data_arrays_, index_i, last_real_particle_index);
                // original_id_[index_i] = original_id_[last_real_particle_index]; // sorted_id_[original_id_[index_i]] = index_i;
                original_id_[last_real_particle_index] = original_id_[index_i]; // sorted_id_[original_id_[index_i]] = index_i;
                std::cout << "l111 :" << ", original_id_[index_i]:" << original_id_[index_i] << ",original_id_[last_real_particle_index]:" << original_id_[last_real_particle_index] << std::flush << std::endl;
                // update original and sorted_id as well
                // original_id_[index_i] = original_id_[last_real_particle_index];
                // original_id_[last_real_particle_index] = temp_original_id_index_i;
                std::cout << "l115 :" << ", original_id_[index_i]:" << original_id_[index_i] << ",original_id_[last_real_particle_index]:" << original_id_[last_real_particle_index] << std::flush << std::endl;
                std::cout << "lets print sth \n";

                // std::swap(original_id_[index_i], original_id_[last_real_particle_index]);
            }
            // MOVE IT OUT !! TODO: Move it out if its last particle need to be removed!
            *total_real_particles_ -= 1;
            // Use an unordered_map to count how many times each ID appears
            std::unordered_map<unsigned int, size_t> frequency;

            // Pass through the array once, counting occurrences
            for (size_t i = 0; i < *total_real_particles_; ++i)
            {
                frequency[original_id_[i]]++;
            }

            // Check which IDs appeared more than once
            for (const auto &entry : frequency)
            {
                if (entry.second > 1)
                {
                    std::cout << "Repeated ID: " << entry.first
                              << " occurs " << entry.second << " times.\n";
                }
            }
            return last_real_particle_index;
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
