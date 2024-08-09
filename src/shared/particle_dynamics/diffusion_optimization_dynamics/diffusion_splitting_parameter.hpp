/**
 * @file 	diffusion_splitting_parameter.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_SPLITTING_PARAMETER_HPP
#define DIFFUSION_SPLITTING_PARAMETER_HPP

#include "diffusion_splitting_parameter.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
ParameterSplittingByPDEInner<DataType>::
    ParameterSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &name)
    : OptimizationBySplittingAlgorithmBase<DataType>(inner_relation, name){};
//=================================================================================================//
template <typename DataType>
ErrorAndParameters<DataType> ParameterSplittingByPDEInner<DataType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    DataType &variable_i = this->variable_[index_i];
    ErrorAndParameters<DataType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        DataType variable_derivative = variable_i + this->variable_[index_j];
        Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
        Real parameter_b = phi_ij * dW_ijV_j * dt / r_ij;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= 1;
    error_and_parameters.error_ -= this->heat_source_[index_i] * dt;
    return error_and_parameters;
}
//=================================================================================================//
template <typename DataType>
void ParameterSplittingByPDEInner<DataType>::
    updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    DataType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->parameter_recovery_[index_i] = this->variable_[index_i];
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;
    if (this->variable_[index_i] < 0.1)
    {
        this->variable_[index_i] = 0.1;
    } // set lower bound

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
        Real parameter_b = phi_ij * dW_ijV_j * dt / r_ij;

        this->parameter_recovery_[index_j] = this->variable_[index_j];
        this->variable_[index_j] += parameter_k * parameter_b;
        if (this->variable_[index_j] < 0.1)
        {
            this->variable_[index_j] = 0.1;
        } // set lower bound
    }
}
//=================================================================================================//
template <typename DataType>
void ParameterSplittingByPDEInner<DataType>::interaction(size_t index_i, Real dt)
{
    /* With Pseudo Time Step */
    this->splitting_index_[index_i] = 0;
    ErrorAndParameters<DataType> error_and_parameters = computeErrorAndParameters(index_i, dt);
    error_and_parameters.error_ = error_and_parameters.error_ - this->residual_T_local_[index_i];
    updateStatesByError(index_i, dt, error_and_parameters);
    this->residual_k_local_[index_i] = error_and_parameters.error_;

    ErrorAndParameters<DataType> error_and_parameters_aft = computeErrorAndParameters(index_i, dt);
    error_and_parameters_aft.error_ = error_and_parameters_aft.error_ - this->residual_T_local_[index_i];
    this->residual_after_splitting_[index_i] = error_and_parameters_aft.error_;

    if (abs(this->residual_after_splitting_[index_i]) > abs(this->residual_k_local_[index_i]))
    {
        this->variable_[index_i] = this->parameter_recovery_[index_i];
        Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t &index_j = inner_neighborhood.j_[n];
            this->variable_[index_j] = this->parameter_recovery_[index_j];
        }

        ErrorAndParameters<DataType> error_and_parameters_reverse = computeErrorAndParameters(index_i, -dt);
        error_and_parameters_reverse.error_ = error_and_parameters_reverse.error_ + this->residual_T_local_[index_i];
        updateStatesByError(index_i, -dt, error_and_parameters_reverse);
        this->residual_k_local_[index_i] = error_and_parameters_reverse.error_;

        ErrorAndParameters<DataType> error_and_parameters_aft = computeErrorAndParameters(index_i, -dt);
        error_and_parameters_aft.error_ = error_and_parameters_aft.error_ + this->residual_T_local_[index_i];
        this->residual_after_splitting_[index_i] = error_and_parameters_aft.error_;

        if (abs(this->residual_after_splitting_[index_i]) > abs(this->residual_k_local_[index_i]))
        {
            this->splitting_index_[index_i] = 0;
            this->variable_[index_i] = this->parameter_recovery_[index_i];
            Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t &index_j = inner_neighborhood.j_[n];
                this->variable_[index_j] = this->parameter_recovery_[index_j];
            }
        }
        else if (abs(this->residual_after_splitting_[index_i]) <= abs(this->residual_k_local_[index_i]))
        {
            this->splitting_index_[index_i] = 2;
        }
    }
    else if (abs(this->residual_after_splitting_[index_i]) <= abs(this->residual_k_local_[index_i]))
    {
        this->splitting_index_[index_i] = 1;
    }
};
//=================================================================================================//
template <typename DataType>
ParameterSplittingByPDEWithBoundary<DataType>::
    ParameterSplittingByPDEWithBoundary(BaseInnerRelation &inner_relation,
                                        BaseContactRelation &contact_relation, const std::string &name)
    : ParameterSplittingByPDEInner<DataType>(inner_relation, name),
      DataDelegateContact(contact_relation)
{
    const std::string &species_name = this->diffusion_.DiffusionSpeciesName();
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        boundary_Vol_.push_back(this->contact_particles_[k]->template registerStateVariable<Real>("VolumetricMeasure"));
        boundary_normal_vector_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        boundary_species_.push_back(this->contact_particles_[k]->template registerStateVariable<Real>(species_name));
        boundary_heat_flux_.push_back(this->contact_particles_[k]->template registerStateVariable<Real>("HeatFlux"));
    }
}
//=================================================================================================//
template <typename DataType>
ErrorAndParameters<DataType> ParameterSplittingByPDEWithBoundary<DataType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<DataType> error_and_parameters =
        ParameterSplittingByPDEInner<DataType>::computeErrorAndParameters(index_i, dt);

    DataType &variable_i = this->variable_[index_i];
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *heat_flux_k = this->boundary_heat_flux_[k];
        Vecd *normal_vector_k = this->boundary_normal_vector_[k];
        Real *Vol_k = this->boundary_Vol_[k];
        Real *species_k = boundary_species_[k];

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t &index_j = contact_neighborhood.j_[n];

            if (species_k[index_j] > 0.0)
            {
                DataType variable_derivative = variable_i;
                Real phi_ij = 2 * (this->species_modified_[index_i] - species_k[index_j]);
                Real parameter_b = 2.0 * phi_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * dt / contact_neighborhood.r_ij_[n];

                error_and_parameters.error_ -= variable_derivative * parameter_b;
                error_and_parameters.a_ += parameter_b;
            }

            if (heat_flux_k[index_j] != 0.0)
            {
                Vecd n_ij = this->normal_vector_[index_i] - normal_vector_k[index_j];
                error_and_parameters.error_ -= heat_flux_k[index_j] * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n].dot(n_ij) * dt;
            }
        }
    }
    return error_and_parameters;
}
//=================================================================================================//
template <typename ParameterSplittingType>
template <typename... Args>
UpdateParameterPDEResidual<ParameterSplittingType>::UpdateParameterPDEResidual(Args &&...args)
    : ParameterSplittingType(std::forward<Args>(args)...){};
//=================================================================================================//
template <typename ParameterSplittingType>
void UpdateParameterPDEResidual<ParameterSplittingType>::interaction(size_t index_i, Real dt)
{
    auto error_and_parameters = this->computeErrorAndParameters(index_i, dt);
    error_and_parameters.error_ = error_and_parameters.error_ - this->residual_T_local_[index_i];
    this->residual_k_global_[index_i] = error_and_parameters.error_;
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_SPLITTING_PARAMETER_HPP
