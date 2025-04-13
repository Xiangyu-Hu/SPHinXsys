#ifndef GENERAL_GRADIENT_HPP
#define GENERAL_GRADIENT_HPP

#include "general_gradient.h"

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
      dv_gradient_(this->particles_->template registerStateVariableOnly<Grad<DataType>>(
          variable_name + "Gradient", ZeroData<Grad<DataType>>::value)),
      dv_B_(this->particles_->template getVariableByName<Matd>("LinearCorrectionMatrix")) {}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
Gradient<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      gradient_(encloser.dv_gradient_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void LinearGradient<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    ScalarVec<DataType, Dimensions> summation = ScalarVec<DataType, Dimensions>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        Scalar<DataType> difference(this->variable_[index_i] - this->variable_[index_j]);
        summation -= scalarProduct(corrected_gradW_ijV_j, difference);
    }
    this->gradient_[index_i] = ScalarVecVecToMatrix(summation);
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void LinearGradient<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    ScalarVec<DataType, Dimensions> summation = ScalarVec<DataType, Dimensions>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        Scalar<DataType> difference(this->variable_[index_i] - this->variable_[index_j]);
        summation -= scalarProduct(corrected_gradW_ijV_j, difference);
    }
    this->gradient_[index_i] += ScalarVecVecToMatrix(summation);
}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
Hessian<Base, DataType, RelationType<Parameters...>>::
    Hessian(DynamicsIdentifier &identifier, std::string &variable_name)
    : BaseDynamicsType(identifier, variable_name),
      dv_M_(this->particles_->template getVariableByName<MatTend>("HessianCorrectionMatrix")),
      dv_hessian_(this->particles_->template getVariableByName<Hess<DataType>>(variable_name + "Hessian")) {}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
Hessian<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      M_(encloser.dv_M_->DelegatedData(ex_policy)),
      hessian_(encloser.dv_hessian_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Hessian<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Hess<DataType> summation = Hess<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  (this->B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        DataType corrected_difference = this->variable_[index_i] - this->variable_[index_j] +
                                        this->gradient_[index_i].dot(r_ij);
        summation += corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     vectorizeTensorSquare(r_ij) * transferToMatrix(corrected_difference).transpose();
    }
    this->hessian_[index_i] = summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Hessian<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Hess<DataType> summation = Hess<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  (this->B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        DataType corrected_difference = this->variable_[index_i] - this->variable_[index_j] +
                                        this->gradient_[index_i].dot(r_ij);
        summation += corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     vectorizeTensorSquare(r_ij) * transferToMatrix(corrected_difference).transpose();
    }
    this->hessian_[index_i] += summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void SecondOrderGradient<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    auto summation = MatrixToScalarVecVec(Grad<DataType>::Zero());
    auto hessian_scalar = MatrixToScalarVecVec(this->hessian_[index_i]);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        auto difference = dotProduct(vectorizeTensorSquare(r_ij), hessian_scalar);
        summation -= 0.5 * scalarProduct(corrected_gradW_ij, difference);
    }
    this->gradient_[index_i] += ScalarVecVecToMatrix(summation);
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void SecondOrderGradient<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    auto summation = MatrixToScalarVecVec(Grad<DataType>::Zero());
    auto hessian_scalar = MatrixToScalarVecVec(this->hessian_[index_i]);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        auto difference = dotProduct(vectorizeTensorSquare(r_ij), hessian_scalar);
        summation -= 0.5 * scalarProduct(corrected_gradW_ij, difference);
    }
    this->gradient_[index_i] += ScalarVecVecToMatrix(summation);
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_GRADIENT_HPP
