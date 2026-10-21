/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPDX-License-Identifier: Apache-2.0
 * ------------------------------------------------------------------------- */
/**
 * @file base_dissipation.h
 * @brief TBD.
 * @details TBD.
 * @author	sphinxsim@einsimo.com
 */

#ifndef BASE_DISSIPATION_H
#define BASE_DISSIPATION_H

#include "interaction_ck.hpp"

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
              zero_value_(ReduceReference<ReduceSum<DataType>>::value) {};

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