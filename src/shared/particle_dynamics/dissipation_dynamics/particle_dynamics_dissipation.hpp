#ifndef PARTICLE_DYNAMICS_DISSIPATION_HPP
#define PARTICLE_DYNAMICS_DISSIPATION_HPP

#include "particle_dynamics_dissipation.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType, typename DampingRateType, class DataDelegationType>
template <class BaseRelationType, typename... Args>
Damping<Base, VariableType, DampingRateType, DataDelegationType>::
    Damping(BaseRelationType &base_relation, const std::string &variable_name, Args &&...args)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation), OperatorSplitting(),
      variable_name_(variable_name), damping_(this->particles_, std::forward<Args>(args)...),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      variable_(*this->particles_->template getVariableDataByName<VariableType>(variable_name)) {}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
ErrorAndParameters<VariableType> Damping<Inner<Projection>, VariableType, DampingRateType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        // linear projection
        VariableType variable_derivative = (this->variable_[index_i] - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->damping_.DampingRate(index_i, index_j) * inner_neighborhood.dW_ij_[n] *
                           this->Vol_[index_i] * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= this->damping_.Capacity(index_i);
    return error_and_parameters;
}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
void Damping<Inner<Projection>, VariableType, DampingRateType>::updateStates(
    size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Real parameter_b = 2.0 * this->damping_.DampingRate(index_i, index_j) * inner_neighborhood.dW_ij_[n] *
                           this->Vol_[index_i] * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        // predicted quantity at particle j
        VariableType variable_j = this->variable_[index_j] - parameter_k * parameter_b;
        VariableType variable_derivative = (this->variable_[index_i] - variable_j);

        // exchange in conservation form
        this->variable_[index_j] -= variable_derivative * parameter_b / this->damping_.Capacity(index_j);
    }
}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
void Damping<Inner<Projection>, VariableType, DampingRateType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
template <typename... Args>
Damping<Contact<Pairwise>, VariableType, DampingRateType>::Damping(Args &&...args)
    : Damping<Base, VariableType, DampingRateType, DataDelegateContact>(std::forward<Args>(args)...)
{
    for (auto &particles : this->contact_particles_)
    {
        contact_Vol_.push_back(particles->template getVariableDataByName<Real>("VolumetricMeasure"));
        contact_variable_.push_back(particles->template getVariableDataByName<VariableType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
void Damping<Inner<Pairwise>, VariableType, DampingRateType>::interaction(size_t index_i, Real dt)
{
    Real Vol_i = this->Vol_[index_i];
    Real capacity_i = this->damping_.Capacity(index_i);

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // forward sweep
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real capacity_j = this->damping_.Capacity(index_j);

        VariableType variable_derivative = (this->variable_[index_i] - this->variable_[index_j]);
        Real parameter_b = this->damping_.DampingRate(index_i, index_j) * inner_neighborhood.dW_ij_[n] *
                           Vol_i * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        VariableType increment = parameter_b * variable_derivative /
                                 (capacity_i * capacity_j - parameter_b * (capacity_i + capacity_j));
        this->variable_[index_i] += increment * capacity_j;
        this->variable_[index_j] -= increment * capacity_i;
    }

    for (size_t n = inner_neighborhood.current_size_; n != 0; --n) // backward sweep
    {
        size_t index_j = inner_neighborhood.j_[n - 1];
        Real capacity_j = this->damping_.Capacity(index_j);

        VariableType variable_derivative = (this->variable_[index_i] - this->variable_[index_j]);
        Real parameter_b = this->damping_.DampingRate(index_i, index_j) * inner_neighborhood.dW_ij_[n - 1] *
                           Vol_i * this->Vol_[index_j] * dt / inner_neighborhood.r_ij_[n - 1];
        VariableType increment = parameter_b * variable_derivative /
                                 (capacity_i * capacity_j - parameter_b * (capacity_i + capacity_j));

        this->variable_[index_i] += increment * capacity_j;
        this->variable_[index_j] -= increment * capacity_i;
    }
}
//=================================================================================================//
template <typename VariableType, typename DampingRateType>
void Damping<Contact<Pairwise, Wall>, VariableType, DampingRateType>::interaction(size_t index_i, Real dt)
{

    Real Vol_i = this->Vol_[index_i];
    Real capacity_i = this->damping_.Capacity(index_i);

    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<VariableType> &variable_k = *(this->contact_variable_[k]);
        StdLargeVec<Real> &Vol_k = *(this->contact_Vol_[k]);
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n) // forward sweep
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real parameter_b = this->damping_.DampingRate(index_i, index_j) * contact_neighborhood.dW_ij_[n] *
                               Vol_i * Vol_k[index_j] * dt / contact_neighborhood.r_ij_[n];

            // only update particle i
            this->variable_[index_i] += parameter_b * (this->variable_[index_i] - variable_k[index_j]) /
                                        (capacity_i - 2.0 * parameter_b);
        }

        for (size_t n = contact_neighborhood.current_size_; n != 0; --n) // backward sweep
        {
            size_t index_j = contact_neighborhood.j_[n - 1];
            Real parameter_b = this->damping_.DampingRate(index_i, index_j) * contact_neighborhood.dW_ij_[n - 1] *
                               Vol_i * Vol_k[index_j] * dt / contact_neighborhood.r_ij_[n - 1];

            // only update particle i
            this->variable_[index_i] += parameter_b * (this->variable_[index_i] - variable_k[index_j]) /
                                        (capacity_i - 2.0 * parameter_b);
        }
    }
}
//=================================================================================================//
template <class DampingAlgorithmType>
template <typename... Args>
DampingWithRandomChoice<DampingAlgorithmType>::
    DampingWithRandomChoice(Real random_ratio, Args &&...args)
    : DampingAlgorithmType(std::forward<Args>(args)...), random_ratio_(random_ratio)
{
}
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