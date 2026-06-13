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
      dv_B_(this->particles_->template registerStateVariable<Matd>(
          "LinearCorrectionMatrix", IdentityMatrix<Matd>::value)) {}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::LinearCorrectionMatrix(
    DynamicsIdentifier &identifier, Real alpha)
    : BaseInteraction(identifier), alpha_(alpha) {}
//=================================================================================================//
template <typename... Parameters>
template <typename BodyRelationType, typename FirstArg>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::LinearCorrectionMatrix(
    DynamicsArgs<BodyRelationType, FirstArg> parameters)
    : LinearCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      B_(encloser.dv_B_->DelegatedDataView(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    B_[index_i] = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd nablaW_ijV_j = this->nablaW_ij(index_i, index_j) * Vol_[index_j];
        B_[index_i] -= this->vec_r_ij(index_i, index_j) * nablaW_ijV_j.transpose();
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : alpha_(encloser.alpha_), B_(encloser.dv_B_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    Real determinant = B_[index_i].determinant();
    Real det_sqr = SMAX(alpha_ - determinant, Real(0));
    Matd inverse = inverseTikhonov(B_[index_i], SqrtEps);
    Real weight = determinant / (determinant + det_sqr);
    B_[index_i] = weight * inverse + (1.0 - weight) * Matd::Identity();
}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
LinearCorrectionMatrix<Contact<Parameters...>>::LinearCorrectionMatrix(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
LinearCorrectionMatrix<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      B_(encloser.dv_B_->DelegatedDataView(ex_policy)),
      contact_Vol_k_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCorrectionMatrix<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd nablaW_ijV_j = this->nablaW_ij(index_i, index_j) * contact_Vol_k_[index_j];
        B_[index_i] -= this->vec_r_ij(index_i, index_j) * nablaW_ijV_j.transpose();
    }
}
//=================================================================================================//
template <class DynamicsIdentifier, class ParticleScope>
template <typename... Args>
LinearCorrectionMatrixScope<DynamicsIdentifier, ParticleScope>::
    LinearCorrectionMatrixScope(DynamicsIdentifier &identifier, Args... args)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_B_(this->particles_->template getVariableByName<Matd>("LinearCorrectionMatrix")),
      within_scope_method_(this->particles_, args...) {}
//=================================================================================================//
template <class DynamicsIdentifier, class ParticleScope>
template <class ExecutionPolicy, class EncloserType>
LinearCorrectionMatrixScope<DynamicsIdentifier, ParticleScope>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : B_(encloser.dv_B_->DelegatedDataView(ex_policy)),
      within_scope_(ex_policy, encloser.within_scope_method_) {}
//=================================================================================================//
template <class DynamicsIdentifier, class ParticleScope>
void LinearCorrectionMatrixScope<DynamicsIdentifier, ParticleScope>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    if (!within_scope_(index_i))
    {
        B_[index_i] = Matd::Identity();
    }
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_CORRECTION_CK_HPP
