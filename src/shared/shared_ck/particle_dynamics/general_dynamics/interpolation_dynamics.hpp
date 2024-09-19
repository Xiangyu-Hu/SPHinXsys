#ifndef INTERPOLATION_DYNAMICS_HPP
#define INTERPOLATION_DYNAMICS_HPP

#include "interpolation_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
Interpolation<Contact<DataType>>::Interpolation(Relation<Contact<>> &pair_contact_relation, const std::string &variable_name)
    : Interaction<Contact<>>(pair_contact_relation),
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
template <typename DataType>
template <class ExecutionPolicy, class EncloserType>
Interpolation<Contact<DataType>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : Interaction<Contact<>>::InteractKernel(ex_policy, encloser, contact_index),
      zero_value_(ZeroData<DataType>::value),
      interpolated_quantities_(encloser.dv_interpolated_quantities_->DelegatedDataField(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataField(ex_policy)),
      contact_data_(encloser.dv_contact_data_[contact_index]->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <typename DataType>
void Interpolation<Contact<DataType>>::InteractKernel::interact(size_t index_i, Real dt)
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
} // namespace SPH
#endif // INTERPOLATION_DYNAMICS_HPP
