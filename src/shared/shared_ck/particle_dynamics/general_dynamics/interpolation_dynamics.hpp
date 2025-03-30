#ifndef INTERPOLATION_DYNAMICS_HPP
#define INTERPOLATION_DYNAMICS_HPP

#include "interpolation_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, typename... Parameters>
Interpolation<Contact<DataType, Parameters...>>::Interpolation(
    Relation<Contact<Parameters...>> &pair_contact_relation, const std::string &variable_name)
    : BaseDynamicsType(pair_contact_relation),
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
      zero_value_(ZeroData<DataType>::value), zero_prediction_(ZeroData<PredictVec<DataType>>::value),
      interpolated_quantities_(encloser.dv_interpolated_quantities_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_data_(encloser.dv_contact_data_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Interpolation<Contact<DataType, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    DataType interpolated_quantity = zero_value_;
    PredictVec<DataType> prediction = zero_prediction_;
    RestoreMatd restoring_matrix = Eps * RestoreMatd::Identity();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real W_ijV_j = this->W_ij(index_i, index_j) * contact_Vol_[index_j];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];

        Real element1 = W_ijV_j;
        Vecd element2 = W_ijV_j * r_ij;
        Vecd element3 = dW_ijV_j * e_ij;
        Matd element4 = dW_ijV_j * r_ij * e_ij.transpose();

        prediction[0] += element1 * contact_data_[index_j];
        for (UnsignedInt i = 0; i < Dimensions; ++i)
        {
            prediction[i + 1] += element3[i] * data_k[index_j];
        }

        restoring_matrix(0, 0) += element1;
        restoring_matrix.block(0, 1, 1, Dimensions) -= element2.transpose();
        restoring_matrix.block(1, 0, Dimensions, 1) += element3;
        restoring_matrix.block(1, 1, Dimensions, Dimensions) -= element4;
    }

    Vecd first_row_components = restoring_matrix.inverse().row(0).transpose();
    interpolated_quantities_[index_i] = first_row_components.dot(prediction);
}
//=================================================================================================//
} // namespace SPH
#endif // INTERPOLATION_DYNAMICS_HPP
