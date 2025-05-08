#ifndef INTERPOLATION_DYNAMICS_HPP
#define INTERPOLATION_DYNAMICS_HPP

#include "interpolation_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, typename... Parameters>
Interpolation<Contact<Base, DataType, Parameters...>>::Interpolation(
    Contact<Parameters...> &pair_contact_relation, const std::string &variable_name)
    : Interaction<Contact<Parameters...>>(pair_contact_relation),
      dv_interpolated_quantities_(this->particles_->template registerStateVariableOnly<DataType>(variable_name))
{
    if (this->contact_particles_.size() > 1)
    {
        std::cout << "\n Error: Interpolation only works for single contact body!" << std::endl;
        exit(1);
    }
    // must be single contact body
    dv_contact_Vol_.push_back(this->contact_particles_[0]->template getVariableByName<Real>("VolumetricMeasure"));
    dv_contact_data_.push_back(this->contact_particles_[0]->template getVariableByName<DataType>(variable_name));
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interpolation<Contact<DataType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      zero_value_(ZeroData<DataType>::value),
      interpolated_quantities_(encloser.dv_interpolated_quantities_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_data_(encloser.dv_contact_data_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Interpolation<Contact<DataType, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    DataType interpolated_quantity(zero_value_);
    Real ttl_weight(0);

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real weight_j = this->W_ij(index_i, index_j) * contact_Vol_[index_j];

        interpolated_quantity += weight_j * contact_data_[index_j];
        ttl_weight += weight_j;
    }
    interpolated_quantities_[index_i] = interpolated_quantity / (ttl_weight + TinyReal);
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interpolation<Contact<DataType, RestoringCorrection, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      zero_prediction_(ZeroData<PredictVecd>::value),
      interpolated_quantities_(encloser.dv_interpolated_quantities_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_data_(encloser.dv_contact_data_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Interpolation<Contact<DataType, RestoringCorrection, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    PredictVecd prediction = zero_prediction_;
    RestoreMatd restoring_matrix = Eps * RestoreMatd::Identity();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real W_ijV_j = this->W_ij(index_i, index_j) * contact_Vol_[index_j];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];

        RestoreMatd approximation_matrix = RestoreMatd::Zero();
        approximation_matrix(0, 0) = W_ijV_j;
        approximation_matrix.block(0, 1, 1, Dimensions) = -W_ijV_j * r_ij.transpose();
        approximation_matrix.block(1, 0, Dimensions, 1) = dW_ijV_j * e_ij;
        approximation_matrix.block(1, 1, Dimensions, Dimensions) = -dW_ijV_j * r_ij * e_ij.transpose();

        RestoreVecd approximation_vector = approximation_matrix.col(0);
        prediction += scalarProduct(approximation_vector, Scalar<DataType>(contact_data_[index_j]));
        restoring_matrix += approximation_matrix;
    }

    RestoreVecd restoring_vector = restoring_matrix.inverse().row(0).transpose();
    interpolated_quantities_[index_i] = dotProduct(restoring_vector, prediction).get();
}
//=================================================================================================//
} // namespace SPH
#endif // INTERPOLATION_DYNAMICS_HPP
