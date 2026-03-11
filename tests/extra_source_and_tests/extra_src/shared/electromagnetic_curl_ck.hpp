#ifndef ELECTROMAGNETIC_CURL_CK_HPP
#define ELECTROMAGNETIC_CURL_CK_HPP

#include "electromagnetic_curl_ck.h"

namespace SPH
{
namespace electromagnetics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
AphiCurlCK<Base, RelationType<Parameters...>>::
    AphiCurlCK(RelationType<Parameters...> &relation, const std::string &variable_name)
    : BaseDynamicsType(relation, variable_name),
      dv_curl_(this->particles_->template registerStateVariable<AngularVecd>(
          variable_name + "Curl", ZeroData<AngularVecd>::value)) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
AphiCurlCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      curl_(encloser.dv_curl_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void AphiCurlCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    AngularVecd summation = ZeroData<AngularVecd>::value;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd difference = this->variable_[index_i] - this->variable_[index_j];
        summation += getCrossProduct(difference, corrected_gradW_ijV_j);
    }
    this->curl_[index_i] = summation;
}
//=================================================================================================//
template <typename... Parameters>
AphiCurlCK<Contact<Parameters...>>::AphiCurlCK(
    Contact<Parameters...> &contact_relation, const std::string &variable_name)
    : BaseDynamicsType(contact_relation, variable_name)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AphiCurlCK<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void AphiCurlCK<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    AngularVecd summation = ZeroData<AngularVecd>::value;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd difference = this->variable_[index_i] - contact_variable_[index_j];
        summation += getCrossProduct(difference, corrected_gradW_ijV_j);
    }
    this->curl_[index_i] += summation;
}
//=================================================================================================//
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_CURL_CK_HPP
