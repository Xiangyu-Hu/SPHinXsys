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
 * @file particle_sort_ck.h
 * @brief Here gives the classes for particle sorting.
 * @author Xiangyu Hu
 */

#ifndef PARTICLE_SORT_H
#define PARTICLE_SORT_H

#include "particle_sorting.h"

/**
 * SPH implementation.
 */
namespace SPH
{
class UpdateSortableVariables
{
    typedef DataAssemble<UniquePtr, DiscreteVariable> TemporaryVariables;

    struct InitializeTemporaryVariables
    {
        template <typename DataType>
        void operator()(UniquePtr<DiscreteVariable<DataType>> &variable_ptr, UnsignedInt data_size);
    };

    BaseParticles *particles_;
    TemporaryVariables temp_variables_;
    OperationOnDataAssemble<TemporaryVariables, InitializeTemporaryVariables> initialize_temp_variables_;

  public:
    UpdateSortableVariables(BaseParticles *particles);

    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    ExecutionPolicy &ex_policy, BaseParticles *particles,
                    DiscreteVariable<UnsignedInt> *dv_index_permutation);
};

class QuickSort
{
    class SwapParticleIndex
    {
        UnsignedInt *sequence_;
        UnsignedInt *index_permutation_;

      public:
        SwapParticleIndex(UnsignedInt *sequence, UnsignedInt *index_permutation);
        ~SwapParticleIndex(){};

        void operator()(UnsignedInt *a, UnsignedInt *b);
    };

  public:
    template <class ExecutionPolicy>
    explicit QuickSort(const ExecutionPolicy &ex_policy,
                       DiscreteVariable<UnsignedInt> *dv_sequence,
                       DiscreteVariable<UnsignedInt> *dv_index_permutation);
    void sort(const ParallelPolicy &ex_policy, BaseParticles *particles);

  protected:
    UnsignedInt *sequence_;
    UnsignedInt *index_permutation_;
    SwapParticleIndex swap_particle_index_;
    CompareParticleSequence compare_;
    tbb::interface9::QuickSortParticleRange<
        UnsignedInt *, CompareParticleSequence, SwapParticleIndex>
        quick_sort_particle_range_;
    tbb::interface9::QuickSortParticleBody<
        UnsignedInt *, CompareParticleSequence, SwapParticleIndex>
        quick_sort_particle_body_;
};

template <class ExecutionPolicy, class SortMethodType>
class ParticleSortCK : public LocalDynamics, public BaseDynamics<void>
{
  public:
    explicit ParticleSortCK(RealBody &real_body);
    virtual ~ParticleSortCK(){};

    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleSortCK<ExecutionPolicy, SortMethodType> &encloser);
        void prepareSequence(UnsignedInt index_i);
        void updateSortedID(UnsignedInt index_i);

      protected:
        Mesh mesh_;

        Vecd *pos_;
        UnsignedInt *sequence_;
        UnsignedInt *index_permutation_;
        UnsignedInt *original_id_;
        UnsignedInt *sorted_id_;
    };

    virtual void exec(Real dt = 0.0) override;
    typedef ParticleSortCK<ExecutionPolicy, SortMethodType> LocalDynamicsType;
    using ComputingKernel = typename LocalDynamicsType::ComputingKernel;

  protected:
    ExecutionPolicy ex_policy_;
    CellLinkedList &cell_linked_list_;
    Mesh mesh_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
    DiscreteVariable<UnsignedInt> *dv_original_id_;
    DiscreteVariable<UnsignedInt> *dv_sorted_id_;
    OperationOnDataAssemble<ParticleVariables, UpdateSortableVariables>
        update_variables_to_sort_;
    SortMethodType sort_method_;
    Implementation<ExecutionPolicy, LocalDynamicsType, ComputingKernel> kernel_implementation_;
};
} // namespace SPH
#endif // PARTICLE_SORT_H
