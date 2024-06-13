#ifndef PARTICLE_DYNAMICS_DISSIPATION_HPP
#define PARTICLE_DYNAMICS_DISSIPATION_HPP

#include "particle_dynamics_dissipation.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType, typename DampingType, class DataDelegationType>
template <class BaseRelationType, typename... Args>
DampingBySplitting<DampingType, VariableType, DataDelegationType>::
    DampingBySplitting(BaseRelationType &base_relation, const std::string &variable_name, Args &&...args)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      damping_(this->particles_, std::forward<Args>(args)...),
      Vol_(*this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      mass_(*this->particles_->template getVariableByName<Real>("Mass")),
      variable_(*this->particles_->template getVariableByName<VariableType>(variable_name)) {}
//=================================================================================================//
template <typename VariableType, typename DampingType>
ErrorAndParameters<VariableType> DampingBySplitting<Inner<>, VariableType, DampingType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        // linear projection
        VariableType variable_derivative = (this->variable_[index_i] - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->damping_.DampingCoefficient(index_i, index_j) * inner_neighborhood.dW_ij_[n] *
                           this->Vol_[index_i] * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= this->mass_[index_i] * this->damping_.SpecificCapacity(index_i);
    return error_and_parameters;
}
//=================================================================================================//
template <typename VariableType, typename DampingType>
void DampingBySplitting<Inner<>, VariableType, DampingType>::updateStates(
    size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Real parameter_b = 2.0 * this->damping_.DampingCoefficient(index_i, index_j) * inner_neighborhood.dW_ij_[n] *
                           this->Vol_[index_i] * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        // predicted quantity at particle j
        VariableType variable_j = this->variable_[index_j] - parameter_k * parameter_b;
        VariableType variable_derivative = (this->variable_[index_i] - variable_j);

        // exchange in conservation form
        this->variable_[index_j] -= variable_derivative * parameter_b /
                                    this->mass_[index_j] / this->damping_.SpecificCapacity(index_j);
    }
}
//=================================================================================================//
template <typename VariableType, typename DampingType>
void DampingBySplitting<Inner<>, VariableType, DampingType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
}
//=================================================================================================//
template <class DampingAlgorithmType>
template <typename... Args>
DampingWithRandomChoice<DampingAlgorithmType>::
    DampingWithRandomChoice(Real random_ratio, Args &&...args)
    : DampingAlgorithmType(std::forward<Args>(args)...), random_ratio_(random_ratio) {}
//=================================================================================================//
template <class DampingAlgorithmType>
bool DampingWithRandomChoice<DampingAlgorithmType>::RandomChoice()
{
    return rand_uniform(0.0, 1.0) < random_ratio_ ? true : false;
}
//=================================================================================================//
template <class DampingAlgorithmType>
void DampingWithRandomChoice<DampingAlgorithmType>::exec(Real dt)
{
    if (RandomChoice())
        DampingAlgorithmType::exec(dt / random_ratio_);
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_DYNAMICS_DISSIPATION_HPP