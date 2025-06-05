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
    : BaseDynamicsType(identifier), variable_name_(variable_name),
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
    Grad<DataType> summation = Grad<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        DataType difference = this->variable_[index_i] - this->variable_[index_j];
        summation -= tensorProduct(corrected_gradW_ijV_j, difference);
    }
    this->gradient_[index_i] = summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <typename... Args>
LinearGradient<Contact<DataType, Parameters...>>::LinearGradient(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
LinearGradient<Contact<DataType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void LinearGradient<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Grad<DataType> summation = Grad<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd corrected_gradW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                     this->B_[index_i] * this->e_ij(index_i, index_j);
        DataType difference = this->variable_[index_i] - contact_variable_[index_j];
        summation -= tensorProduct(corrected_gradW_ijV_j, difference);
    }
    this->gradient_[index_i] += summation;
}
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <typename... Args>
Hessian<Base, DataType, RelationType<Parameters...>>::Hessian(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_M_(this->particles_->template getVariableByName<MatTend>("HessianCorrectionMatrix")),
      dv_hessian_(this->particles_->template registerStateVariableOnly<Hess<DataType>>(
          this->variable_name_ + "Hessian")) {}
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
        DataType corrected_difference = this->variable_[index_i] - this->variable_[index_j] -
                                        this->gradient_[index_i].dot(r_ij);
        summation += 2.0 * corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     tensorProduct(vectorizeTensorSquare(r_ij), corrected_difference);
    }
    this->hessian_[index_i] = this->M_[index_i] * summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <typename... Args>
Hessian<Contact<DataType, Parameters...>>::Hessian(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Hessian<Contact<DataType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
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
        DataType corrected_difference = this->variable_[index_i] - contact_variable_[index_j] -
                                        this->gradient_[index_i].dot(r_ij);
        summation += 2.0 * corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     tensorProduct(vectorizeTensorSquare(r_ij), corrected_difference);
    }
    this->hessian_[index_i] += this->M_[index_i] * summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void SecondOrderGradient<Inner<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Grad<DataType> summation = Grad<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        auto difference = this->variable_[index_i] - this->variable_[index_j] +
                          0.5 * vectorizeTensorSquare(r_ij).transpose() * this->hessian_[index_i];
        summation -= corrected_gradW_ij * difference;
    }
    this->gradient_[index_i] = summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <typename... Args>
SecondOrderGradient<Contact<DataType, Parameters...>>::SecondOrderGradient(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
SecondOrderGradient<Contact<DataType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void SecondOrderGradient<Contact<DataType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Grad<DataType> summation = Grad<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Vecd corrected_gradW_ij = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  this->B_[index_i] * this->e_ij(index_i, index_j);
        auto difference = this->variable_[index_i] - contact_variable_[index_j] +
                          0.5 * vectorizeTensorSquare(r_ij).transpose() * this->hessian_[index_i];
        summation -= corrected_gradW_ij * difference;
    }
    this->gradient_[index_i] += summation;
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_GRADIENT_HPP
