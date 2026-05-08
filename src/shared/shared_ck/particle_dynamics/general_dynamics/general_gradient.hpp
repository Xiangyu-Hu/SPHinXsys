#ifndef GENERAL_GRADIENT_HPP
#define GENERAL_GRADIENT_HPP

#include "general_gradient.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, template <typename...> class RelationType, typename... Parameters>
Gradient<Base, DataType, RelationType<Parameters...>>::
    Gradient(RelationType<Parameters...> &relation, const std::string &variable_name)
    : BaseDynamicsType(relation), variable_name_(variable_name),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      dv_gradient_(this->particles_->template registerStateVariable<Grad<DataType>>(
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
void Gradient<Inner<DataType, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
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
Gradient<Contact<DataType, Parameters...>>::Gradient(
    Contact<Parameters...> &contact_relation, const std::string &variable_name)
    : BaseDynamicsType(contact_relation, variable_name)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Gradient<Contact<DataType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, typename... Parameters>
void Gradient<Contact<DataType, Parameters...>>::
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
template <typename TransportType, template <typename...> class RelationType, typename... Parameters>
template <typename... Args>
Hessian<Base, TransportType, RelationType<Parameters...>>::Hessian(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      transport_(nullptr),
      dv_M_(this->particles_->template getVariableByName<MatTend>("HessianCorrectionMatrix")),
      dv_hessian_(this->particles_->template registerStateVariable<Hess<DataType>>(
          this->variable_name_ + "Hessian"))
{
    if constexpr (!std::is_same_v<DataType, TransportType>)
    {
        transport_ = &DynamicsCast<TransportType>(this, this->sph_body_->getBaseMaterial());
    }
}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Hessian<Inner<TransportType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser),
      inter_particle_coeff_(),
      M_(encloser.dv_M_->DelegatedData(ex_policy)),
      hessian_(encloser.dv_hessian_->DelegatedData(ex_policy))
{
    if constexpr (!std::is_same_v<DataType, TransportType>)
    {
        inter_particle_coeff_ =
            encloser.transport_->getInterParticleCoeff(ex_policy, *encloser.transport_);
    }
}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
void Hessian<Inner<TransportType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Hess<DataType> summation = Hess<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j] *
                                  (this->B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        DataType corrected_difference = inter_particle_coeff_(index_i, index_j) *
                                        (this->variable_[index_i] - this->variable_[index_j] -
                                         this->gradient_[index_i].dot(r_ij));
        summation += 2.0 * corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     tensorProduct(vectorizeTensorSquare(r_ij), corrected_difference);
    }
    this->hessian_[index_i] = this->M_[index_i] * summation;
}
//=================================================================================================//
template <typename TransportType, typename... OtherTransportType,
          template <typename...> class InterfaceType, typename... Parameters>
Hessian<Contact<InterfaceType<TransportType, OtherTransportType...>, Parameters...>>::
    Hessian(Contact<Parameters...> &contact_relation, std::string &variable_name,
            StdVec<InterfaceType<TransportType, OtherTransportType...> *> interface_models)
    : BaseDynamicsType(contact_relation, variable_name), interface_models_(interface_models)
{
    if (interface_models_.size() != this->contact_bodies_.size())
    {
        throw std::runtime_error(
            "Error: the number of interface models should be the same as that of contact bodies.");
    }

    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(this->variable_name_));
    }
}
//=================================================================================================//
template <typename TransportType, typename... OtherTransportType,
          template <typename...> class InterfaceType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Hessian<Contact<InterfaceType<TransportType, OtherTransportType...>, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      interface_coeff_(encloser.interface_models_[contact_index]->getInterfaceCoeff(ex_policy)),
      M_(encloser.dv_M_->DelegatedData(ex_policy)),
      hessian_(encloser.dv_hessian_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.dv_contact_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename TransportType, typename... OtherTransportType,
          template <typename...> class InterfaceType, typename... Parameters>
void Hessian<Contact<InterfaceType<TransportType, OtherTransportType...>, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Hess<DataType> summation = Hess<DataType>::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd r_ij = this->vec_r_ij(index_i, index_j);
        Real corrected_dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j] *
                                  (this->B_[index_i] * this->e_ij(index_i, index_j)).dot(r_ij);
        DataType corrected_difference = interface_coeff_(index_i, index_j) *
                                        (this->variable_[index_i] - contact_variable_[index_j] -
                                         this->gradient_[index_i].dot(r_ij));
        summation += 2.0 * corrected_dW_ijV_j / math::pow(r_ij.squaredNorm(), 2) *
                     tensorProduct(vectorizeTensorSquare(r_ij), corrected_difference);
    }
    hessian_[index_i] += M_[index_i] * summation;
}
//=================================================================================================//
template <typename DataType, typename... Parameters>
template <typename... Args>
Hessian<Contact<DataType, Parameters...>>::Hessian(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
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
      M_(encloser.dv_M_->DelegatedData(ex_policy)),
      hessian_(encloser.dv_hessian_->DelegatedData(ex_policy)),
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
    hessian_[index_i] += M_[index_i] * summation;
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
{}
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
template <typename... Parameters>
LinearCurl<Inner<WithUpdate, Parameters...>>::LinearCurl(
    Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDynamicsType(inner_relation, variable_name),
      dv_curl_(this->particles_->template registerStateVariable<Curl>(
          variable_name + "Curl", ZeroData<Curl>::value)) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
LinearCurl<Inner<WithUpdate, Parameters...>>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : gradient_(encloser.dv_gradient_->DelegatedData(ex_policy)),
      curl_(encloser.dv_curl_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void LinearCurl<Inner<WithUpdate, Parameters...>>::UpdateKernel::update(size_t index_i, Real dt)
{
    curl_[index_i] = skewSymmetric(gradient_[index_i]);
}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
DoubleCurl<Inner<TransportType, Parameters...>>::DoubleCurl(
    Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDynamicsType(inner_relation, variable_name),
      dv_double_curl_(this->particles_->template registerStateVariable<Vecd>(
          variable_name + "Curl", ZeroData<Vecd>::value)) {}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DoubleCurl<Inner<TransportType, Parameters...>>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : hessian_(encloser.dv_hessian_->DelegatedData(ex_policy)),
      double_curl_(encloser.dv_double_curl_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
void DoubleCurl<Inner<TransportType, Parameters...>>::UpdateKernel::update(size_t index_i, Real dt)
{
    double_curl_[index_i] = computeDoubleCurl(this->hessian_[index_i]);
}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
Vec2d DoubleCurl<Inner<TransportType, Parameters...>>::UpdateKernel::
    computeDoubleCurl(Hess<Vec2d> &hessian)
{
    double C1 = hessian(1, 1) - hessian(2, 0);  // ∂_xy f2 - ∂_yy f1
    double C2 = -hessian(0, 1) + hessian(1, 0); // -∂_xx f2 + ∂_xy f1
    return Vec2d(C1, C2);
}
//=================================================================================================//
template <typename TransportType, typename... Parameters>
Vec3d DoubleCurl<Inner<TransportType, Parameters...>>::UpdateKernel::
    computeDoubleCurl(Hess<Vec3d> &hessian)
{
    double C1 = hessian(1, 1) + hessian(2, 2) - hessian(3, 0) - hessian(5, 0); // ∂_xy f2 + ∂_xz f3 - ∂_yy f1 - ∂_zz f1
    double C2 = hessian(1, 0) + hessian(4, 2) - hessian(0, 1) - hessian(5, 1); // ∂_xy f1 + ∂_yz f3 - ∂_xx f2 - ∂_zz f2
    double C3 = hessian(2, 0) + hessian(4, 1) - hessian(0, 2) - hessian(3, 2); // ∂_xz f1 + ∂_yz f2 - ∂_xx f3 - ∂_yy f3
    return Vec3d(C1, C2, C3);
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_GRADIENT_HPP
