#ifndef EROSION_DYNAMICS_1ST_HPP
#define EROSION_DYNAMICS_1ST_HPP

#include "base_particles.hpp"
#include "erosion_dynamics_1st_ck.h"
namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
PlasticAcousticStepWithErosion<BaseInteractionType>::PlasticAcousticStepWithErosion(DynamicsIdentifier &identifier)
    : PlasticAcousticStep<BaseInteractionType>(identifier),
      dv_shear_vel_(this->particles_->template registerStateVariableOnly<Vecd>("ShearVelocity")),
      dv_friction_angle_(this->particles_->template registerStateVariableOnly<Real>("FrictionAngle",this->plastic_continuum_.getFrictionAngle())), 
      dv_cohesion_(this->particles_->template registerStateVariableOnly<Real>("Cohesion",this->plastic_continuum_.getCohesion())),
      dv_reduction_para_(this->particles_->template registerStateVariableOnly<Real>("ReductionParameter",Real(1.0))),
      dv_plastic_label_(this->particles_->template registerStateVariableOnly<int>("PlasticLabel",int(0))),
      dv_total_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("TotalStressTensor")),
      dv_viscous_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("ViscousStressTensor")),
      dv_shear_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("ShearStressTensor")),
      dv_indicator_(this->particles_->template getVariableByName<int>("Indicator"))
{
    this->particles_->template addEvolvingVariable<Vecd>("ShearVelocity");
    this->particles_->template addEvolvingVariable<Real>("FrictionAngle");
    this->particles_->template addEvolvingVariable<Real>("Cohesion");
    this->particles_->template addEvolvingVariable<Real>("ReductionParameter");
    this->particles_->template addEvolvingVariable<int>("PlasticLabel");
    this->particles_->template addEvolvingVariable<Mat3d>("TotalStressTensor");
    this->particles_->template addEvolvingVariable<Mat3d>("ViscousStressTensor");
    this->particles_->template addEvolvingVariable<Mat3d>("ShearStressTensor");
}

//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStepWithErosion1stHalf(Inner<Parameters...> &inner_relation)
    : PlasticAcousticStepWithErosion<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_), riemann_solver_(this->plastic_continuum_, this->plastic_continuum_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    p_[index_i] = -stress_tensor_3D_[index_i].trace() / 3;
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      viscous_stress_tensor_3D_(encloser.dv_viscous_stress_tensor_3D_->DelegatedData(ex_policy)),
      plastic_label_(encloser.dv_plastic_label_->DelegatedData(ex_policy)),
      shear_stress_tensor_3D_(encloser.dv_shear_stress_tensor_3D_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    Real rho_i = rho_[index_i];
    Matd stress_tensor_i = degradeToMatd(stress_tensor_3D_[index_i]);
    Matd viscous_stress_tensor_i = degradeToMatd(viscous_stress_tensor_3D_[index_i]);
    Matd shear_stress_tensor_i = degradeToMatd(shear_stress_tensor_3D_[index_i]);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd nablaW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j] * this->e_ij(index_i, index_j);
        Matd stress_tensor_j = degradeToMatd(stress_tensor_3D_[index_j]);
        force += mass_[index_i] * rho_[index_j] * ((stress_tensor_i + stress_tensor_j) / (rho_i * rho_[index_j])) * nablaW_ijV_j;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
        /*Viscous force*/
        Matd viscous_stress_tensor_j = degradeToMatd(viscous_stress_tensor_3D_[index_j]);
        force += mass_[index_i] * rho_[index_j] * ((viscous_stress_tensor_i + viscous_stress_tensor_j) / (rho_i * rho_[index_j])) * nablaW_ijV_j;
        /*Shear force*/
        Matd shear_stress_tensor_j = degradeToMatd(shear_stress_tensor_3D_[index_j]);
        if(plastic_label_[index_i] == 1)
            force += mass_[index_i] * rho_[index_j] * ((shear_stress_tensor_i + shear_stress_tensor_j) / (rho_i * rho_[index_j])) * nablaW_ijV_j;
    }
    force_[index_i] += force;
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}

//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStepWithErosion1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStepWithErosion1stHalf(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation), Interaction<Wall>(wall_contact_relation),
      correction_(this->particles_), riemann_solver_(this->plastic_continuum_, this->plastic_continuum_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStepWithErosion1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      viscous_stress_tensor_3D_(encloser.dv_viscous_stress_tensor_3D_->DelegatedData(ex_policy)),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedData(ex_policy)),
      wall_acc_ave_(encloser.dv_wall_acc_ave_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStepWithErosion1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);

    Matd stress_tensor_i = degradeToMatd(stress_tensor_3D_[index_i]);
    Matd viscous_stress_tensor_i = degradeToMatd(viscous_stress_tensor_3D_[index_i]);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        Real face_wall_external_acceleration = (force_prior_[index_i] / mass_[index_i] - wall_acc_ave_[index_j]).dot(-e_ij);
        Real p_in_wall = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
        force += 2 * mass_[index_i] * (stress_tensor_i + viscous_stress_tensor_i) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_in_wall) * dW_ijV_j;
    }
    force_[index_i] += force / rho_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
} // namespace continuum_dynamics
} // namespace SPH
#endif // EROSION_DYNAMICS_1ST_HPP