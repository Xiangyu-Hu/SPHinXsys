#ifndef SPLITTING_DISSIPATION_HPP
#define SPLITTING_DISSIPATION_HPP

#include "splitting_dissipation.h"

namespace SPH
{
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    ProjectionDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDissipationType::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    // compute the error and parameters
    DataType error = this->zero_value_;
    Real parameter_a = 0.0;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    Real parameter_c = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        UnsignedInt neighbor_j = n - this->FirstNeighbor(index_i);

        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) *
                                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                                  this->vec_r_ij(index_i, index_j).norm();
        error -= (this->variable_[index_i] - this->variable_[index_j]) * parameter_b[neighbor_j];

        parameter_a += parameter_b[neighbor_j];
        parameter_c += parameter_b[neighbor_j] * parameter_b[neighbor_j];
    }
    parameter_a = (parameter_a - 1.0);

    // update the variable at particle i only seems more accurate
    DataType parameter_k = error / (parameter_c + parameter_a * parameter_a + TinyReal);
    this->variable_[index_i] += parameter_k * parameter_a;
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
ProjectionDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    ProjectionDissipation(Contact<Parameters...> &contact_relation, const std::string &variable_name)
    : BaseInteraction(contact_relation, variable_name)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_dv_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(variable_name));
    }
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ProjectionDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      inverse_capacity_(encloser.dissipation_model_),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.contact_dv_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
void ProjectionDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    DataType error = this->variable_[index_i];
    Real parameter_a = 1.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        Real parameter_b = 2.0 * this->dis_coeff_(index_i, index_j) *
                           this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                           dt * inverse_capacity_(index_i) /
                           this->vec_r_ij(index_i, index_j).norm();
        error -= parameter_b * contact_variable_[index_j];
        parameter_a -= parameter_b;
    }
    this->variable_[index_i] = error * parameter_a / (parameter_a * parameter_a + TinyReal);
}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    PairwiseDissipation(
        Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDissipationType::InteractKernel(ex_policy, encloser),
      inverse_capacity_(encloser.dissipation_model_) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real split_dt = dt * 0.5;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    Real Vol_i = this->Vol_[index_i];
    // forward sweep
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        UnsignedInt neighbor_j = n - this->FirstNeighbor(index_i);

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) * this->dW_ij(index_i, index_j) *
                                  split_dt / this->vec_r_ij(index_i, index_j).norm();

        DataType increment = 0.5 * parameter_b[neighbor_j] * pair_difference /
                             (1.0 - parameter_b[neighbor_j] * (Vol_i + this->Vol_[index_j]));
        atomic_add(this->variable_[index_i], increment * inverse_capacity_(index_i) * this->Vol_[index_j]);
        atomic_add(this->variable_[index_j], -increment * inverse_capacity_(index_j) * Vol_i);
    }
    // backward sweep
    for (UnsignedInt n = this->LastNeighbor(index_i); n != this->FirstNeighbor(index_i); --n)
    {
        UnsignedInt index_j = this->neighbor_index_[n - 1];
        UnsignedInt neighbor_j = n - 1 - this->FirstNeighbor(index_i);

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        DataType increment = 0.5 * parameter_b[neighbor_j] * pair_difference /
                             (1.0 - parameter_b[neighbor_j] * (Vol_i + this->Vol_[index_j]));
        atomic_add(this->variable_[index_i], increment * inverse_capacity_(index_i) * this->Vol_[index_j]);
        atomic_add(this->variable_[index_j], -increment * inverse_capacity_(index_j) * Vol_i);
    }
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
PairwiseDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    PairwiseDissipation(Contact<Parameters...> &contact_relation, const std::string &variable_name)
    : BaseInteraction(contact_relation, variable_name)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_dv_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(variable_name));
    }
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PairwiseDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      inverse_capacity_(encloser.dissipation_model_),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.contact_dv_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
void PairwiseDissipation<Contact<Splitting, Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    DataType error = this->variable_[index_i];
    Real parameter_a = 1.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        Real parameter_b = 2.0 * this->dis_coeff_(index_i, index_j) *
                           this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                           dt * inverse_capacity_(index_i) /
                           this->vec_r_ij(index_i, index_j).norm();
        error -= parameter_b * contact_variable_[index_j];
        parameter_a -= parameter_b;
    }
    this->variable_[index_i] = error * parameter_a / (parameter_a * parameter_a + TinyReal);
}
//=================================================================================================//
} // namespace SPH
#endif // SPLITTING_DISSIPATION_HPP
