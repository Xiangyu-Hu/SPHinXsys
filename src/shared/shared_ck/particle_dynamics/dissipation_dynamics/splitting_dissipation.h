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
 * @file splitting_dissipation.h
 * @brief TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef SPLITTING_DISSIPATION_H
#define SPLITTING_DISSIPATION_H

#include "base_dissipation.h"

namespace SPH
{
// Note: ProjectionDissipation method is for obtaining accurate solution for steady problems.
// This method is unconditionally stable, but has small conservative errors.
// Therefore, it is not used for transient system which requires conservation properties.
template <typename...>
class ProjectionDissipation;

template <typename DissipationType, typename... Parameters>
class ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>
    : public Dissipation<Base, DissipationType, Inner<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseDissipationType = Dissipation<Base, DissipationType, Inner<Parameters...>>;

  public:
    explicit ProjectionDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~ProjectionDissipation() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);
    };
};

// Note: PairwiseDissipation method is for obtaining solution for both for steady and transient problems.
// This method is unconditionally stable, conservative and able to obtaining converged solution.
// However the method is less accurate than explicit methods.
// Therefore, it is used when explicit time-step size is too small to be practical,
// or when the problem is too stiff for explicit methods.
template <typename...>
class PairwiseDissipation;

template <typename DissipationType, typename... Parameters>
class PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>
    : public Dissipation<Base, DissipationType, Inner<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseDissipationType = Dissipation<Base, DissipationType, Inner<Parameters...>>;
    using InverseVolumetricCapacity = typename DissipationType::InverseVolumetricCapacity;

  public:
    explicit PairwiseDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~PairwiseDissipation() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InverseVolumetricCapacity inverse_capacity_;
    };
};

template <class DissipationType, template <typename...> class BoundaryType, typename... Parameters>
class PairwiseDissipation<Contact<BoundaryType<DissipationType>, Parameters...>>
    : public Dissipation<Base, DissipationType, Contact<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseInteraction = Dissipation<Base, DissipationType, Interaction<Contact<Parameters...>>>;
    using BoundaryKernel = typename BoundaryType<DissipationType>::ComputingKernel;
    UniquePtrsKeeper<DiscreteVariableArray<DataType>> contact_transfer_array_ptrs_keeper_;
    UniquePtrsKeeper<BoundaryType<DissipationType>> boundary_ptrs_keeper_;

  public:
    PairwiseDissipation(Contact<Parameters...> &contact_relation, const std::string &variable_name);
    virtual ~PairwiseDissipation() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                       UnsignedInt contact_index);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataArray<DataType> *contact_transfer_;
        BoundaryKernel boundary_flux_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariableArray<DataType> *> contact_dv_transfer_array_;
    StdVec<BoundaryType<DissipationType> *> contact_boundary_method_;
};
} // namespace SPH
#endif // SPLITTING_DISSIPATION_H