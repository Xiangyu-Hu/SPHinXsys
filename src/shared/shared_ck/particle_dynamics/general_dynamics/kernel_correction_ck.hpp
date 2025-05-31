#ifndef KERNEL_CORRECTION_CK_HPP
#define KERNEL_CORRECTION_CK_HPP

#include "kernel_correction_ck.h"

namespace SPH
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
LinearCorrectionMatrix<Base, RelationType<Parameters...>>::
    LinearCorrectionMatrix(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_B_(this->particles_->template registerStateVariableOnly<Matd>(
          "LinearCorrectionMatrix", IdentityMatrix<Matd>::value)) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
LinearCorrectionMatrix<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   LinearCorrectionMatrix<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>> &encloser)
    : LinearCorrectionMatrix<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Matd local_configuration = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] * this->e_ij(index_i, index_j);
        local_configuration -= this->vec_r_ij(index_i, index_j) * gradW_ij.transpose();
    }
    this->B_[index_i] = local_configuration;
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy,
                 LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>> &encloser)
    : LinearCorrectionMatrix<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      alpha_(encloser.alpha_) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    Real determinant = this->B_[index_i].determinant();
    Real det_sqr = SMAX(alpha_ - determinant, Real(0));
    Matd B_T = this->B_[index_i].transpose(); // for Tikhonov regularization
    Matd inverse = (B_T * this->B_[index_i] + SqrtEps * Matd::Identity()).inverse() * B_T;
    Real weight = determinant / (determinant + det_sqr);
    this->B_[index_i] = weight * inverse + (1.0 - weight) * Matd::Identity();
}
//=================================================================================================//
template <typename... Parameters>
LinearCorrectionMatrix<Contact<Parameters...>>::
    LinearCorrectionMatrix(Contact<Parameters...> &contact_relation)
    : LinearCorrectionMatrix<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
LinearCorrectionMatrix<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   LinearCorrectionMatrix<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : LinearCorrectionMatrix<Base, Contact<Parameters...>>::
          InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_k_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Matd local_configuration = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd gradW_ij = this->dW_ij(index_i, index_j) * contact_Vol_k_[index_j] * this->e_ij(index_i, index_j);
        local_configuration -= this->vec_r_ij(index_i, index_j) * gradW_ij.transpose();
    }
    this->B_[index_i] += local_configuration;
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_CORRECTION_CK_HPP
