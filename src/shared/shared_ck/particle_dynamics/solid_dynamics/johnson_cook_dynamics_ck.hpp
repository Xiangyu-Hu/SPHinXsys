#ifndef JOHNSON_COOK_DYNAMICS_CK_HPP
#define JOHNSON_COOK_DYNAMICS_CK_HPP

#include "adaptation.h"
#include "johnson_cook_dynamics_ck.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::JohnsonCookIntegration1stHalfCK(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), BaseStructureIntegration1stHalf(this->particles_),
      material_(DynamicCast<JohnsonCookSolid>(this, this->sph_body_->getMatterMaterial())),
      h_ref_(this->sph_body_->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_active_(this->particles_->template getVariableByName<int>("MaterialActive")) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InitializeKernel::InitializeKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : constitute_(ex_policy, encloser.material_), rho0_(encloser.material_.ReferenceDensity()),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)),
      dF_dt_(encloser.dv_dF_dt_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)),
      active_(encloser.dv_active_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InitializeKernel::initialize(size_t i, Real dt)
{
    pos_[i] += vel_[i] * dt * Real(0.5);
    if (!active_[i])
    {
        stress_on_particle_[i].setZero();
        return;
    }
    F_[i] += dF_dt_[i] * dt * Real(0.5);
    const Matd stress = constitute_.UpdateStress(F_[i], i, dt);
    if (constitute_.Status(i) != 0)
        return;
    const Real J = F_[i].determinant();
    rho_[i] = rho0_ / J;
    stress_on_particle_[i] = J * stress * F_[i].inverse().transpose() * B_[i].transpose();
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      damping_scale_(Real(0.5) * encloser.material_.ReferenceDensity() * encloser.material_.ReferenceSoundSpeed() * encloser.h_ref_),
      Vol0_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)),
      active_(encloser.dv_active_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::interact(size_t i, Real)
{
    Vecd sum = Vecd::Zero();
    if (active_[i])
    {
        const Vecd zero = Vecd::Zero();
        const Real inv_W0 = Real(1) / this->W0(i, zero);
        for (UnsignedInt n = this->FirstNeighbor(i); n != this->LastNeighbor(i); ++n)
        {
            const UnsignedInt j = this->neighbor_index_[n];
            if (!active_[j])
                continue;
            const Real inverse_distance = Real(Dimensions) / this->vec_r_ij(i, j).norm();
            const Real rate = inverse_distance * inverse_distance * (pos_[i] - pos_[j]).dot(vel_[i] - vel_[j]);
            const Real weight = this->W_ij(i, j) * inv_W0;
            const Matd damping = Real(0.5) * (F_[i] + F_[j]) * damping_scale_ * rate;
            sum += (stress_on_particle_[i] + stress_on_particle_[j] + Real(0.25) * weight * damping) *
                   this->nablaW_ij(i, j) * Vol0_[j];
        }
    }
    force_[i] = Vol0_[i] * sum;
}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
JohnsonCookIntegration2ndHalfCK<Inner<OneLevel, Parameters...>>::JohnsonCookIntegration2ndHalfCK(DynamicsIdentifier &identifier)
    : Base(identifier), dv_active_(this->particles_->template getVariableByName<int>("MaterialActive")) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
JohnsonCookIntegration2ndHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : Base::InteractKernel(ex_policy, encloser), active_(encloser.dv_active_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void JohnsonCookIntegration2ndHalfCK<Inner<OneLevel, Parameters...>>::InteractKernel::interact(size_t i, Real)
{
    Matd rate = Matd::Zero();
    if (active_[i])
        for (UnsignedInt n = this->FirstNeighbor(i); n != this->LastNeighbor(i); ++n)
        {
            const UnsignedInt j = this->neighbor_index_[n];
            if (active_[j])
                rate -= (this->vel_[i] - this->vel_[j]) * (this->nablaW_ij(i, j) * this->Vol0_[j]).transpose();
        }
    this->dF_dt_[i] = rate * this->B_[i];
}
//=================================================================================================//
inline JohnsonCookAcousticTimeStepCK::JohnsonCookAcousticTimeStepCK(SPHBody &body, Real cfl)
    : LocalDynamicsReduce<ReduceMin<Real>>(body),
      h_min_(body.getSPHAdaptation().MinimumSmoothingLength()),
      c0_(DynamicCast<JohnsonCookSolid>(this, body.getMatterMaterial()).ReferenceSoundSpeed()), cfl_(cfl),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
JohnsonCookAcousticTimeStepCK::ReduceKernel::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : h_min_(encloser.h_min_), c0_(encloser.c0_), cfl_(encloser.cfl_),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Real JohnsonCookAcousticTimeStepCK::ReduceKernel::reduce(size_t i, Real)
{
    const Real acceleration = ((force_[i] + force_prior_[i]) / mass_[i]).norm();
    return cfl_ * math::min(math::sqrt(h_min_ / (acceleration + TinyReal)), h_min_ / (c0_ + vel_[i].norm()));
}
} // namespace solid_dynamics
} // namespace SPH
#endif // JOHNSON_COOK_DYNAMICS_CK_HPP
