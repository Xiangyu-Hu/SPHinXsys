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
 * @file particle_sorting.h
 * @brief Here gives the classes for particle sorting.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_SORTING_H
#define PARTICLE_SORTING_H

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.hpp"
#include "algorithm_primitive.h"

namespace SPH
{
class BaseParticles;

struct SwapParticleDataValue
{
    template <typename DataType>
    void operator()(DataContainerKeeper<AllocatedData<DataType>> &data_keeper, size_t index_a, size_t index_b) const
    {
        for (size_t i = 0; i != data_keeper.size(); ++i)
        {
            DataType *data_field = data_keeper[i];
            std::swap(data_field[index_a], data_field[index_b]);
        }
    };
};

/**
 * @class SwapSortableParticleData
 * @brief swap sortable particle data according to a sequence
 */
class SwapSortableParticleData
{
  protected:
    UnsignedInt *sequence_;
    ParticleData &evolving_variables_data_;
    OperationOnDataAssemble<ParticleData, SwapParticleDataValue> swap_particle_data_value_;

  public:
    explicit SwapSortableParticleData(BaseParticles *base_particles);
    ~SwapSortableParticleData() {};

    /** the operator overload for swapping particle data.
     *  the arguments are the same with std::iter_swap
     */
    void operator()(UnsignedInt *a, UnsignedInt *b);
};

class ParticleSequence : public LocalDynamics
{
  protected:
    Vecd *pos_;
    UnsignedInt *sequence_;
    BaseCellLinkedList &cell_linked_list_;

  public:
    explicit ParticleSequence(RealBody &real_body);
    virtual ~ParticleSequence() {};
    void update(size_t index_i, Real dt = 0.0);
};

template <class ExecutionPolicy>
class ParticleDataSort;

template <>
class ParticleDataSort<ParallelPolicy>
    : public LocalDynamics, public BaseDynamics<void>
{
  protected:
    UnsignedInt *sequence_;

    /** using pointer because it is constructed after particles. */
    SwapSortableParticleData swap_sortable_particle_data_;
    CompareSequence compare_;
    tbb::interface9::QuickSortRange<
        UnsignedInt *, CompareSequence, SwapSortableParticleData>
        quick_sort_particle_range_;
    tbb::interface9::QuickSortBody<
        UnsignedInt *, CompareSequence, SwapSortableParticleData>
        quick_sort_particle_body_;

  public:
    explicit ParticleDataSort(RealBody &real_body);
    virtual ~ParticleDataSort() {};
    virtual void exec(Real dt = 0.0) override;
};

class UpdateSortedID : public LocalDynamics
{
  protected:
    UnsignedInt *original_id_;
    UnsignedInt *sorted_id_;

  public:
    explicit UpdateSortedID(RealBody &real_body);
    virtual ~UpdateSortedID() {};
    void update(size_t index_i, Real dt = 0.0);
};

template <class ExecutionPolicy = ParallelPolicy>
class ParticleSorting : public BaseDynamics<void>
{
    SimpleDynamics<ParticleSequence, ExecutionPolicy> particle_sequence_;
    ParticleDataSort<ParallelPolicy> particle_data_sort_;
    SimpleDynamics<UpdateSortedID, ExecutionPolicy> update_sorted_id_;

  public:
    ParticleSorting(RealBody &real_body);
    virtual ~ParticleSorting() {};

    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // PARTICLE_SORTING_H
