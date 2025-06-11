/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'SourceTermsis) is an acronym from Smoothed Particle *
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

void atomic_add(Real &variable, const Real &value)
{
    AtomicRef<Real> atomic_variable(variable);
    atomic_variable.fetch_add(value);
};

void atomic_add(Vecd &variable, const Vecd &value)
{
    for (size_t i = 0; i < Dimensions; ++i)
    {
        AtomicRef<Real> atomic_variable(variable[i]);
        atomic_variable.fetch_add(value[i]);
    }
};

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
    explicit PairwiseDissipation(Inner<Parameters...> &inner_relation,
                                 const std::string &variable_name);
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

class ConstantSource
{
  public:
    ConstantSource(BaseParticles *particles, Real constant_value = 0.0)
        : constant_value_(constant_value) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : constant_value_(encloser.constant_value_){};

        // This function returns the constant value for all particles.
        Real operator()(size_t) const { return constant_value_; }

      private:
        Real constant_value_;
    };

  private:
    Real constant_value_;
};

template <class DissipationType, typename... Parameters>
class PairwiseDissipation<Contact<Dirichlet<DissipationType>, Parameters...>>
    : public Dissipation<Base, DissipationType, Contact<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseInteraction = Dissipation<Base, DissipationType, Contact<Parameters...>>;
    using InverseVolumetricCapacity = typename DissipationType::InverseVolumetricCapacity;

  public:
    PairwiseDissipation(Contact<Parameters...> &contact_relation, const std::string &variable_name);
    virtual ~PairwiseDissipation() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        InverseVolumetricCapacity inverse_capacity_;
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariable<DataType> *> contact_dv_variable_;
};
} // namespace SPH
#endif // SPLITTING_DISSIPATION_H