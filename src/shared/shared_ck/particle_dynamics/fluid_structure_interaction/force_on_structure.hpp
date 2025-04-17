#ifndef FORCE_ON_STRUCTURE_HPP
#define FORCE_ON_STRUCTURE_HPP

#include "force_on_structure.h"

namespace SPH
{
namespace FSI
{
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ContactRelationType>
ForceFromFluid<KernelCorrectionType, Parameters...>::
    ForceFromFluid(ContactRelationType &contact_relation, const std::string &force_name)
    : Interaction<Contact<Parameters...>>(contact_relation), ForcePriorCK(this->particles_, force_name),
      solid_(DynamicCast<Solid>(this, this->sph_body_.getBaseMaterial())),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_force_from_fluid_(this->dv_current_force_),
      dv_vel_ave_(solid_.AverageVelocityVariable(this->particles_))
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_kernel_correction_.push_back(KernelCorrectionType(this->contact_particles_[k]));
        contact_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_vel_.push_back(this->contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ForceFromFluid<KernelCorrectionType, Parameters...>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : Interaction<Contact<Parameters...>>::InteractKernel(ex_policy, encloser, contact_index),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      force_from_fluid_(encloser.dv_force_from_fluid_->DelegatedData(ex_policy)),
      vel_ave_(encloser.dv_vel_ave_->DelegatedData(ex_policy)),
      contact_correction_(ex_policy, encloser.contact_kernel_correction_[contact_index]),
      contact_Vol_(encloser.contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_vel_(encloser.contact_vel_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename ViscousForceType, typename... Parameters>
template <class ContactRelationType>
ViscousForceFromFluid<Contact<WithUpdate, ViscousForceType, Parameters...>>::
    ViscousForceFromFluid(ContactRelationType &contact_relation)
    : BaseForceFromFluid(contact_relation, "ViscousForceFromFluid")
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        ViscosityType *viscosity_model_k = DynamicCast<ViscosityType>(this, &this->contact_particles_[k]->getBaseMaterial());
        contact_viscosity_model_.push_back(viscosity_model_k);
        contact_smoothing_length_sq_.push_back(pow(this->contact_bodies_[k]->getSPHAdaptation().ReferenceSmoothingLength(), 2));
    }
}
//=================================================================================================//
template <typename ViscousForceType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ViscousForceFromFluid<Contact<WithUpdate, ViscousForceType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseForceFromFluid::InteractKernel(ex_policy, encloser, contact_index),
      viscosity_(ex_policy, *encloser.contact_viscosity_model_[contact_index]),
      smoothing_length_sq_(encloser.contact_smoothing_length_sq_[contact_index]) {}
//=================================================================================================//
template <typename ViscousForceType, typename... Parameters>
void ViscousForceFromFluid<Contact<WithUpdate, ViscousForceType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
        Vecd vel_derivative = (this->vel_ave_[index_i] - this->contact_vel_[index_j]) /
                              (vec_r_ij.squaredNorm() + 0.01 * smoothing_length_sq_);

        force += vec_r_ij.dot(this->contact_correction_(index_j) * e_ij) *
                 viscosity_(index_j) * vel_derivative *
                 this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j];
    }

    this->force_from_fluid_[index_i] = force * this->Vol_[index_i];
}
//=================================================================================================//
template <class AcousticStep2ndHalfType, typename... Parameters>
template <class ContactRelationType>
PressureForceFromFluid<Contact<WithUpdate, AcousticStep2ndHalfType, Parameters...>>::
    PressureForceFromFluid(ContactRelationType &contact_relation)
    : BaseForceFromFluid(contact_relation, "PressureForceFromFluid"),
      dv_acc_ave_(this->solid_.AverageAccelerationVariable(this->particles_)),
      dv_n_(this->particles_->template getVariableByName<Vecd>("NormalDirection"))
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        FluidType &contact_fluid_k = DynamicCast<FluidType>(this, this->contact_particles_[k]->getBaseMaterial());
        contact_riemann_solver_.push_back(RiemannSolverType(contact_fluid_k, contact_fluid_k));
        dv_contact_rho_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Density"));
        dv_contact_mass_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Mass"));
        dv_contact_p_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Pressure"));
        dv_contact_force_prior_.push_back(this->contact_particles_[k]->template getVariableByName<Vecd>("ForcePrior"));
    }
}
//=================================================================================================//
template <class AcousticStep2ndHalfType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PressureForceFromFluid<Contact<WithUpdate, AcousticStep2ndHalfType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseForceFromFluid::InteractKernel(ex_policy, encloser, contact_index),
      acc_ave_(encloser.dv_acc_ave_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      riemann_solver_(encloser.contact_riemann_solver_[contact_index]),
      contact_rho_(encloser.dv_contact_rho_[contact_index]->DelegatedData(ex_policy)),
      contact_mass_(encloser.dv_contact_mass_[contact_index]->DelegatedData(ex_policy)),
      contact_p_(encloser.dv_contact_p_[contact_index]->DelegatedData(ex_policy)),
      contact_force_prior_(encloser.dv_contact_force_prior_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class AcousticStep2ndHalfType, typename... Parameters>
void PressureForceFromFluid<Contact<WithUpdate, AcousticStep2ndHalfType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        Real face_wall_external_acceleration =
            (contact_force_prior_[index_j] / contact_mass_[index_j] - acc_ave_[index_i]).dot(e_ij);
        Real p_j_in_wall = contact_p_[index_j] + contact_rho_[index_j] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
        Real u_jump = 2.0 * (this->contact_vel_[index_j] - this->vel_ave_[index_i]).dot(n_[index_i]);
        force -= (riemann_solver_.DissipativePJump(u_jump) * n_[index_i] +
                  (p_j_in_wall + contact_p_[index_j]) * this->contact_correction_(index_j) * e_ij) *
                 this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j];
    }

    this->force_from_fluid_[index_i] = force * this->Vol_[index_i];
}
//=================================================================================================//
} // namespace FSI
} // namespace SPH
#endif // FORCE_ON_STRUCTURE_HPP