// SPDX-License-Identifier: Apache-2.0
#ifndef HJC_DYNAMICS_CK_HPP
#define HJC_DYNAMICS_CK_HPP

#include "hjc_dynamics_ck.h"
#include "adaptation.h"

namespace SPH
{
namespace solid_dynamics
{
template <typename... Parameters>
template <class DynamicsIdentifier>
HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::HJCIntegration1stHalfCK(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), BaseStructureIntegration1stHalf(this->particles_),
      material_(DynamicCast<HJCSolid>(this, this->sph_body_->getMatterMaterial())),
      h_ref_(this->sph_body_->getSPHAdaptation().ReferenceSmoothingLength()) {}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InitializeKernel::InitializeKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : constitute_(ex_policy, encloser.material_), rho0_(encloser.material_.ReferenceDensity()),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)),
      dF_dt_(encloser.dv_dF_dt_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)) {}

template <typename... Parameters>
void HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InitializeKernel::initialize(size_t i, Real dt)
{
    pos_[i] += vel_[i] * dt * Real(0.5);
    F_[i] += dF_dt_[i] * dt * Real(0.5);
    Matd stress = constitute_.UpdateStress(F_[i], i, dt);
    // The caller checks HJCIntegrationStatus before advancing the second half step.
    if (constitute_.Status(i) != 0)
        return;
    Real J = F_[i].determinant();
    rho_[i] = rho0_ / J;
    stress_on_particle_[i] = J * stress * F_[i].inverse().transpose() * B_[i].transpose();
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      damping_scale_(Real(0.5) * encloser.material_.ReferenceDensity() *
                     encloser.material_.ReferenceSoundSpeed() * encloser.h_ref_),
      Vol0_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)) {}

template <typename... Parameters>
void HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::interact(size_t i, Real dt)
{
    Vecd sum = Vecd::Zero();
    const Vecd zero = Vecd::Zero();
    Real inv_W0 = Real(1) / this->W0(i, zero);
    for (UnsignedInt n = this->FirstNeighbor(i); n != this->LastNeighbor(i); ++n)
    {
        UnsignedInt j = this->neighbor_index_[n];
        Real r = this->vec_r_ij(i, j).norm();
        Real dim_over_r = Real(Dimensions) / r;
        Real strain_rate = dim_over_r * dim_over_r * (pos_[i] - pos_[j]).dot(vel_[i] - vel_[j]);
        Real weight = this->W_ij(i, j) * inv_W0;
        Matd damping = Real(0.5) * (F_[i] + F_[j]) * damping_scale_ * strain_rate;
        sum += (stress_on_particle_[i] + stress_on_particle_[j] + Real(0.25) * weight * damping) *
               this->nablaW_ij(i, j) * Vol0_[j];
    }
    force_[i] = Vol0_[i] * sum;
}

inline HJCAcousticTimeStepCK::HJCAcousticTimeStepCK(SPHBody &body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMin<Real>>(body),
      h_ref_(body.getSPHAdaptation().ReferenceSmoothingLength()), acousticCFL_(acousticCFL),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_modulus_(particles_->getVariableByName<Real>("HJCAcousticModulus")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}

template <class ExecutionPolicy, class EncloserType>
HJCAcousticTimeStepCK::ReduceKernel::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : h_ref_(encloser.h_ref_), acousticCFL_(encloser.acousticCFL_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      modulus_(encloser.dv_modulus_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}

inline Real HJCAcousticTimeStepCK::ReduceKernel::reduce(size_t i, Real dt)
{
    Real acceleration = (force_[i] + force_prior_[i]).norm() / mass_[i];
    return acousticCFL_ * SMIN(math::sqrt(h_ref_ / (acceleration + TinyReal)),
        h_ref_ / (math::sqrt(modulus_[i] / rho_[i]) + vel_[i].norm()));
}
} // namespace solid_dynamics
} // namespace SPH
#endif
