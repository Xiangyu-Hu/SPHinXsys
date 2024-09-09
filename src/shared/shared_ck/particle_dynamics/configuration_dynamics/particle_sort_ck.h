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
    typedef DataAssemble<UniquePtrKeeper, DiscreteVariable> TemporaryVariables;

    struct InitializeTemporaryVariables
    {
        template <typename DataType>
        void operator()(UniquePtrKeeper<DiscreteVariable<DataType>> &variable_ptr_keeper, UnsignedInt data_size);
    };

    BaseParticles *particles_;
    TemporaryVariables temp_variables_;
    OperationOnDataAssemble<TemporaryVariables, InitializeTemporaryVariables> initialize_temp_variables_;

  public:
    UpdateSortableVariables(BaseParticles *particles);

    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    ExecutionPolicy &ex_policy, DiscreteVariable<UnsignedInt> *dv_index_permutation);
};

} // namespace SPH
#endif // PARTICLE_SORT_H
