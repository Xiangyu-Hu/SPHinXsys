#ifndef HESSIAN_CORRECTION_CK_HPP
#define HESSIAN_CORRECTION_CK_HPP

#include "hessian_correction_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
Gradient<Base, DataType, RelationType<Parameters...>>::
    Gradient(DynamicsIdentifier &identifier, std::string &variable_name)
    : BaseDynamicsType(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      dv_gradient_(this->particles_->template registerStateVariableOnly<GradType>(
          variable_name + "Gradient", ZeroData<GradType>::value)) {}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
Gradient<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      gradient_(encloser.dv_gradient_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
HessianMatrix<Base, DataType, RelationType<Parameters...>>::
    HessianMatrix(DynamicsIdentifier &identifier, std::string &variable_name)
    : BaseDynamicsType(identifier, variable_name),
      dv_B_(this->particles_->template getVariableByName<Matd>("LinearCorrectionMatrix")),
      dv_M_(this->particles_->template getVariableByName<MatTend>("HessianCorrectionMatrix")),
      dv_hessian(this->particles_->template getVariableByName<HessianType>(variable_name + "Hessian")) {}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
HessianMatrix<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      B_(encloser.dv_B_->DelegatedData(ex_policy)), M_(encloser.dv_M_->DelegatedData(ex_policy)),
      hessian_(encloser.dv_hessian_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void HessianMatrix<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    HessianType summation = HessianType::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  (B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        VecMatd displacement_matrix = vectorizeSymMatrix(r_ij * r_ij.transpose());
        DataType corrected_difference = this->variable_[index_i] - this->variable_[index_j] +
                                        this->gradient_[index_i].dot(r_ij);
        summation += corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     displacement_matrix * transferToMatrix(corrected_difference).transpose();
    }
    hessian_[index_i] = summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void HessianMatrix<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    HessianType summation = HessianType::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  (B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        VecMatd displacement_matrix = vectorizeSymMatrix(r_ij * r_ij.transpose());
        DataType corrected_difference = this->variable_[index_i] - this->variable_[index_j] +
                                        this->gradient_[index_i].dot(r_ij);
        summation += corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     displacement_matrix * transferToMatrix(corrected_difference).transpose();
    }
    hessian_[index_i] += summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void SecondOrderGradient<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    GradType summation = GradType::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  B_[index_i] * this->e_ij(index_i, index_j);
        VecMatd displacement_matrix = vectorizeSymMatrix(r_ij * r_ij.transpose());

        summation -= corrected_gradW_ij;

                     displacement_matrix *
                     transferToMatrix(corrected_difference).transpose();
    }
    hessian_[index_i] = summation;
}
//=================================================================================================//
} // namespace SPH
#endif // HESSIAN_CORRECTION_CK_HPP
