#ifndef HESSIAN_CORRECTION_CK_HPP
#define HESSIAN_CORRECTION_CK_HPP

#include "hessian_correction_ck.h"

namespace SPH
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
HessianCorrectionMatrix<Base, RelationType<Parameters...>>::
    HessianCorrectionMatrix(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_B_(this->particles_->template registerStateVariableOnly<Matd>(
          "LinearCorrectionMatrix", IdentityMatrix<Matd>::value)),
      dv_displacement_matrix_grad_(this->particles_->template registerStateVariableOnly<VecMatGrad>(
          "DisplacementMatrixGradient", ZeroData<VecMatGrad>::value)),
      dv_M_(this->particles_->template registerStateVariableOnly<MatTend>(
          "HessianCorrectionMatrix", IdentityMatrix<MatTend>::value)) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
HessianCorrectionMatrix<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   HessianCorrectionMatrix<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), B_(encloser.dv_B_->DelegatedData(ex_policy)),
      displacement_matrix_grad_(encloser.dv_displacement_matrix_grad_->DelegatedData(ex_policy)),
      M_(encloser.dv_M_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DisplacementMatrixGradient<Inner<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    VecMatGrad grad_displacement_matrix = VecMatGrad::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);

        grad_displacement_matrix -= (vectorizeTensorSquare(r_ij) * corrected_gradW_ij.transpose());
    }
    this->displacement_matrix_grad_[index_i] = grad_displacement_matrix;
}
//=================================================================================================//
template <typename... Parameters>
void DisplacementMatrixGradient<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    VecMatGrad grad_displacement_matrix = VecMatGrad::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);

        grad_displacement_matrix -= (vectorizeTensorSquare(r_ij) * corrected_gradW_ij.transpose());
    }
    this->displacement_matrix_grad_[index_i] += grad_displacement_matrix;
}
//=================================================================================================//
template <typename... Parameters>
void HessianCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    MatTend approximation_tensor = MatTend::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        VecMatd displacement_matrix = vectorizeTensorSquare(r_ij);
        MatTend displacement_tensor = displacement_matrix * displacement_matrix.transpose();
        VecMatd displacement_matrix_increment = this->displacement_matrix_grad_[index_i].dot(r_ij);

        approximation_tensor -= (displacement_tensor + displacement_matrix_increment) *
                                r_ij.dot(corrected_gradW_ij) / math::pow(r_ij.squaredNorm(), 2);
    }
    this->M_[index_i] = approximation_tensor;
}
//=================================================================================================//
template <typename... Parameters>
void HessianCorrectionMatrix<Inner<WithUpdate, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    Real determinant = this->M_[index_i].determinant();
    Real det_sqr = SMAX(alpha_ - determinant, Real(0));
    Matd M_T = this->M_[index_i].transpose(); // for Tikhonov regularization
    Matd inverse = (M_T * this->M_[index_i] + SqrtEps * Matd::Identity()).inverse() * M_T;
    Real weight = determinant / (determinant + det_sqr);
    this->M_[index_i] = weight * inverse + (1.0 - weight) * Matd::Identity();
}
//=================================================================================================//
template <typename... Parameters>
void HessianCorrectionMatrix<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    MatTend approximation_tensor = MatTend::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        VecMatd displacement_matrix = vectorizeTensorSquare(r_ij);
        MatTend displacement_tensor = displacement_matrix * displacement_matrix.transpose();
        VecMatd displacement_matrix_increment = this->displacement_matrix_grad_[index_i].dot(r_ij);

        approximation_tensor -= (displacement_tensor + displacement_matrix_increment) *
                                r_ij.dot(corrected_gradW_ij) / math::pow(r_ij.squaredNorm(), 2);
    }
    this->M_[index_i] += approximation_tensor;
}
//=================================================================================================//
} // namespace SPH
#endif // HESSIAN_CORRECTION_CK_HPP
