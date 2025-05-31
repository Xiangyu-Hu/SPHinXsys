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

#include "base_general_dynamics.h"
#include "interaction_ck.hpp"
#include "particle_dynamics_dissipation.h"

namespace SPH
{
template <typename... RelationTypes>
class Dissipation;

template <typename DissipationType, template <typename...> class RelationType, typename... Parameters>
class Dissipation<Base, DissipationType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using InterParticleDiffusionCoeff = typename DissipationType::InterParticleDiffusionCoeff;

  public:
    explicit Dissipation(RelationType<Parameters...> &relation, const std::string &variable_name);
    virtual ~Dissipation() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        InterParticleDiffusionCoeff dis_coeff_;
        Real *Vol_;
        DataType *variable_;
        DataType zero_error_;
    };

  protected:
    DissipationType &dissipation_model_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<DataType> *dv_variable_;
};

// Note: ConservativeDamping method is not for obtaining accurate solution,
// but for stabilization the system with a numerical damping coefficient.
// This method is unconditionally stable and has zero-order consistency.
// It is used to damp the high-frequency oscillations in the system.
template <typename...>
class ConservativeDamping;

template <typename DampingType, typename... Parameters>
class ConservativeDamping<Inner<Splitting, DampingType, Parameters...>>
    : public Dissipation<Base, DampingType, Inner<Parameters...>>
{
    using DataType = typename DampingType::DataType;
    using BaseDampingType = Dissipation<Base, DampingType, Inner<Parameters...>>;

  public:
    explicit ConservativeDamping(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~ConservativeDamping() {};

    class InteractKernel : public BaseDampingType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);
    };
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
    explicit PairwiseDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~PairwiseDissipation() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InverseVolumetricCapacity capacity_;
    };
};
} // namespace SPH
#endif // SPLITTING_DISSIPATION_H