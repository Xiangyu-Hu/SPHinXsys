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
 * @file base_dissipation.h
 * @brief TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef BASE_DISSIPATION_H
#define BASE_DISSIPATION_H

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
    explicit Dissipation(RelationType<Parameters...> &relation, const std::string &variable_name)
        : Interaction<RelationType<Parameters...>>(relation),
          dissipation_model_(DynamicCast<DissipationType>(this, this->particles_->getBaseMaterial())),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)) {};
    virtual ~Dissipation() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : Interaction<RelationType<Parameters...>>::InteractKernel(
                  ex_policy, encloser, std::forward<Args>(args)...),
              dis_coeff_(ex_policy, encloser.dissipation_model_),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
              zero_value_(ReduceReference<ReduceSum<DataType>>::value){};

      protected:
        InterParticleDiffusionCoeff dis_coeff_;
        Real *Vol_;
        DataType *variable_;
        DataType zero_value_;
    };

  protected:
    DissipationType &dissipation_model_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<DataType> *dv_variable_;
};
} // namespace SPH
#endif // BASE_DISSIPATION_H