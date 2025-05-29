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

template <typename DissipationType, typename SourceType, typename... Parameters>
class Dissipation<Inner<Splitting, DissipationType, SourceType, Parameters...>>
    : public Dissipation<Base, DissipationType, Inner<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseDissipationType = Dissipation<Base, DissipationType, Inner<Parameters...>>;
    using Sourcekernel = typename SourceType::ComputingKernel;

  public:
    explicit Dissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~Dissipation() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Sourcekernel source_;
    };

  protected:
    SourceType source_model_;
};

class NoSource
{
  public:
    NoSource(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, NoSource &encloser){};

        Real operator()(UnsignedInt index_i) { return 0.0; };
    };
};
} // namespace SPH
#endif // SPLITTING_DISSIPATION_H