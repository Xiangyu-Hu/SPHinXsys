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
 * @file particle_sort_ck.h
 * @brief Here gives the classes for particle sorting.
 * @author Xiangyu Hu
 */

#ifndef PARTICLE_SORT_H
#define PARTICLE_SORT_H

#include "base_configuration_dynamics.h"
#include "particle_sorting.h"

/**
 * SPH implementation.
 */
namespace SPH
{
template <class ExecutionPolicy>
class ParticleSortCK : public LocalDynamics, public BaseDynamics<void>
{
    using SortMethodType = typename SortMethod<ExecutionPolicy>::type;

  public:
    explicit ParticleSortCK(RealBody &real_body);
    virtual ~ParticleSortCK() {};

    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleSortCK<ExecutionPolicy> &encloser);
        void prepareSequence(UnsignedInt index_i);

      protected:
        Mesh mesh_;

        Vecd *pos_;
        UnsignedInt *sequence_;
        UnsignedInt *index_permutation_;
    };

    virtual void exec(Real dt = 0.0) override;
    typedef ParticleSortCK<ExecutionPolicy> LocalDynamicsType;

  protected:
    ExecutionPolicy ex_policy_;
    BaseCellLinkedList &cell_linked_list_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
    OperationOnDataAssemble<DiscreteVariables, UpdateSortableVariables<DiscreteVariable>> update_variables_to_sort_;
    SortMethodType sort_method_;
    Implementation<ExecutionPolicy, LocalDynamicsType, ComputingKernel> kernel_implementation_;
};
} // namespace SPH
#endif // PARTICLE_SORT_H
