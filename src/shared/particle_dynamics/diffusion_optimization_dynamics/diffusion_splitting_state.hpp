/**
 * @file 	diffusion_splitting_state.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_SPLITTING_STATE_HPP
#define DIFFUSION_SPLITTING_STATE_HPP

#include "diffusion_splitting_state.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType>
TemperatureSplittingByPDEInner<VariableType>::
    TemperatureSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : OptimizationBySplittingAlgorithmBase<VariableType>(inner_relation, variable_name){};
//=================================================================================================//
template <typename VariableType>
ErrorAndParameters<VariableType> TemperatureSplittingByPDEInner<VariableType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    VariableType &variable_i = this->variable_[index_i];
    ErrorAndParameters<VariableType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t &index_j = inner_neighborhood.j_[n];
        Real &r_ij_ = inner_neighborhood.r_ij_[n];
        Vecd &e_ij_ = inner_neighborhood.e_ij_[n];

        // linear projection
        VariableType variable_derivative = (variable_i - this->variable_[index_j]);
        Real diff_coff_ij = this->diffusion_.getInterParticleDiffusionCoeff(index_i, index_j, e_ij_);
        Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * dt / r_ij_;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= 1;
    error_and_parameters.error_ -= this->heat_source_[index_i] * dt;
    return error_and_parameters;
};
//=================================================================================================//
template <typename VariableType>
void TemperatureSplittingByPDEInner<VariableType>::
    updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t &index_j = inner_neighborhood.j_[n];
        Real &r_ij_ = inner_neighborhood.r_ij_[n];
        Vecd &e_ij_ = inner_neighborhood.e_ij_[n];

        Real diff_coff_ij = this->diffusion_.getInterParticleDiffusionCoeff(index_i, index_j, e_ij_);
        Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * dt / r_ij_;
        this->variable_[index_j] -= parameter_k * parameter_b;
    }
}
//=================================================================================================//
template <typename VariableType>
void TemperatureSplittingByPDEInner<VariableType>::
    interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStatesByError(index_i, dt, error_and_parameters);
    this->residual_T_local_[index_i] = error_and_parameters.error_;
}
//=================================================================================================//
template <typename VariableType>
TemperatureSplittingByPDEWithBoundary<VariableType>::
    TemperatureSplittingByPDEWithBoundary(BaseInnerRelation &inner_relation,
                                          BaseContactRelation &contact_relation, const std::string &variable_name)
    : TemperatureSplittingByPDEInner<VariableType>(inner_relation, variable_name),
      DataDelegateContactOnly(contact_relation)
{
    boundary_heat_flux_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        boundary_Vol_.push_back(this->contact_particles_[k]->template registerSharedVariable<Real>("VolumetricMeasure"));
        boundary_normal_vector_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        boundary_variable_.push_back(this->contact_particles_[k]->template registerSharedVariable<VariableType>(variable_name));
        boundary_heat_flux_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HeatFlux");
    }
};
//=================================================================================================//
template <typename VariableType>
ErrorAndParameters<VariableType> TemperatureSplittingByPDEWithBoundary<VariableType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters =
        TemperatureSplittingByPDEInner<VariableType>::computeErrorAndParameters(index_i, dt);

    VariableType &variable_i = this->variable_[index_i];
    /* contact interaction. */
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &Vol_k = *(this->boundary_Vol_[k]);
        StdLargeVec<Real> &heat_flux_k = *(this->boundary_heat_flux_[k]);
        StdLargeVec<Vecd> &normal_vector_k = *(this->boundary_normal_vector_[k]);
        StdLargeVec<VariableType> &variable_k = *(this->boundary_variable_[k]);

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t &index_j = contact_neighborhood.j_[n];

            if (variable_k[index_j] > 0.0)
            {
                // linear projection
                VariableType variable_derivative = 2 * (variable_i - variable_k[index_j]);
                Real diff_coff_ij = this->diffusion_.getDiffusionCoeffWithBoundary(index_i);
                Real parameter_b = 2.0 * diff_coff_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * dt / contact_neighborhood.r_ij_[n];

                error_and_parameters.error_ -= variable_derivative * parameter_b;
                error_and_parameters.a_ += parameter_b;
            }

            Vecd n_ij = this->normal_vector_[index_i] - normal_vector_k[index_j];
            error_and_parameters.error_ -= heat_flux_k[index_j] * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n].dot(n_ij) * dt;
        }
    }
    return error_and_parameters;
}
//=================================================================================================//
template <typename TemperatureSplittingType>
template <typename... Args>
UpdateTemperaturePDEResidual<TemperatureSplittingType>::
    UpdateTemperaturePDEResidual(Args &&...args)
    : TemperatureSplittingType(std::forward<Args>(args)...){};
//=================================================================================================//
template <typename TemperatureSplittingType>
void UpdateTemperaturePDEResidual<TemperatureSplittingType>::
    interaction(size_t index_i, Real dt)
{
    auto error_and_parameters = this->computeErrorAndParameters(index_i, dt);
    this->residual_T_global_[index_i] = error_and_parameters.error_;
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_SPLITTING_STATE_HPP
