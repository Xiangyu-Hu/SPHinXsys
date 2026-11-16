// SPDX-License-Identifier: Apache-2.0
#include "hjc_solid.hpp"
#include "base_particles.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SPH
{
HJCSolid::HJCSolid(Real rho0, const HJCParameters &p)
    : ElasticSolid(rho0), parameters_(p), lock_compression_(0), transition_slope_(0),
      stress_(nullptr), previous_deformation_(nullptr),
      damage_(nullptr), plastic_strain_(nullptr), plastic_volume_(nullptr),
      maximum_compression_(nullptr), pressure_(nullptr), equivalent_stress_(nullptr), acoustic_modulus_(nullptr)
{
    const Real values[] = {rho0, p.shear_modulus, p.A, p.B, p.C, p.N, p.compressive_strength,
                           p.tensile_strength, p.reference_strain_rate, p.minimum_fracture_strain,
                           p.maximum_strength, p.crush_pressure, p.crush_strain, p.lock_pressure,
                           p.lock_strain, p.D1, p.D2, p.K1, p.K2, p.K3};
    for (Real value : values)
        if (!std::isfinite(value))
            throw std::invalid_argument("HJC parameters must be finite");
    if (!(rho0 > 0 && p.shear_modulus > 0 && p.A >= 0 && p.B >= 0 && p.C >= 0 && p.N > 0 &&
          p.compressive_strength > 0 && p.tensile_strength > 0 && p.reference_strain_rate > 0 &&
          p.minimum_fracture_strain > 0 && p.maximum_strength > 0 && p.crush_pressure > 0 &&
          p.crush_strain > 0 && p.lock_pressure > p.crush_pressure && p.lock_strain > 0 &&
          p.D1 > 0 && p.D2 >= 0 && p.K1 > 0 && p.K2 >= 0 && p.K3 >= 0))
        throw std::invalid_argument("Invalid HJC parameter range");
    G0_ = p.shear_modulus;
    K0_ = p.crush_pressure / p.crush_strain;
    E0_ = 9 * K0_ * G0_ / (3 * K0_ + G0_);
    nu_ = (3 * K0_ - 2 * G0_) / (2 * (3 * K0_ + G0_));
    setSoundSpeeds();
    if (!std::isfinite(K0_) || !std::isfinite(E0_) || !std::isfinite(c0_))
        throw std::invalid_argument("HJC elastic moduli must be finite");
    contact_stiffness_ = rho0 * c0_ * c0_;
    material_type_name_ = "HJCSolid";
    // Match the transition to the dense polynomial at P_lock, continuously.
    Real lo = p.lock_strain;
    Real hi = p.lock_strain + (1 + p.lock_strain) * p.lock_pressure / p.K1;
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        Real mid = (lo + hi) * 0.5;
        if (DensePressure(mid) < p.lock_pressure)
            lo = mid;
        else
            hi = mid;
    }
    lock_compression_ = (lo + hi) * 0.5;
    transition_slope_ = (p.lock_pressure - p.crush_pressure) / (lock_compression_ - p.crush_strain);
    if (!(lock_compression_ > p.crush_strain && transition_slope_ < K0_ &&
          transition_slope_ < p.lock_pressure / (lock_compression_ - p.lock_strain)))
        throw std::invalid_argument("HJC crushing branch must be softer than elastic unloading");
}

Real HJCSolid::DensePressure(Real compression) const
{
    return ConstitutiveModel().DensePressure(compression);
}

Real HJCSolid::Pressure(Real compression, HJCState &state) const
{
    return ConstitutiveModel().Pressure(compression, state);
}

Real HJCSolid::FractureStrain(Real pressure) const
{
    return ConstitutiveModel().FractureStrain(pressure);
}

Real HJCSolid::YieldStrength(Real pressure, Real damage, Real strain_rate) const
{
    return ConstitutiveModel().YieldStrength(pressure, damage, strain_rate);
}

void HJCSolid::Integrate(const Mat3d &increment, Real compression, Real dt, HJCState &state) const
{
    HJCIntegrationStatus status = ConstitutiveModel().Integrate(increment, compression, dt, state);
    if (status == HJCIntegrationStatus::invalid_increment)
        throw std::invalid_argument("HJC requires finite increments, positive J and a positive timestep");
    if (status != HJCIntegrationStatus::success)
        throw std::domain_error("HJC stress update produced a non-finite state");
}

Real HJCSolid::AcousticModulus(Real compression) const
{
    return ConstitutiveModel().AcousticModulus(compression);
}

void HJCSolid::initializeLocalParameters(BaseParticles *particles)
{
    ElasticSolid::initializeLocalParameters(particles);
    dv_stress_ = particles->registerStateVariable<Mat3d>("StressCauchy");
    dv_previous_deformation_ = particles->registerStateVariable<Matd>("HJCPreviousDeformation", IdentityMatrix<Matd>::value);
    dv_damage_ = particles->registerStateVariable<Real>("HJCDamage");
    dv_plastic_strain_ = particles->registerStateVariable<Real>("HJCPlasticStrain");
    dv_plastic_volume_ = particles->registerStateVariable<Real>("HJCPlasticVolume");
    dv_maximum_compression_ = particles->registerStateVariable<Real>("HJCMaximumCompression");
    dv_pressure_ = particles->registerStateVariable<Real>("Pressure");
    dv_equivalent_stress_ = particles->registerStateVariable<Real>("VonMisesStress");
    dv_acoustic_modulus_ = particles->registerStateVariable<Real>("HJCAcousticModulus", AcousticModulus(0));
    dv_status_ = particles->registerStateVariable<int>("HJCIntegrationStatus");
    stress_ = dv_stress_->Data();
    previous_deformation_ = dv_previous_deformation_->Data();
    damage_ = dv_damage_->Data();
    plastic_strain_ = dv_plastic_strain_->Data();
    plastic_volume_ = dv_plastic_volume_->Data();
    maximum_compression_ = dv_maximum_compression_->Data();
    pressure_ = dv_pressure_->Data();
    equivalent_stress_ = dv_equivalent_stress_->Data();
    acoustic_modulus_ = dv_acoustic_modulus_->Data();
    particles->addEvolvingVariable<Mat3d>("StressCauchy");
    particles->addEvolvingVariable<Matd>("HJCPreviousDeformation");
    for (const char *name : {"HJCDamage", "HJCPlasticStrain", "HJCPlasticVolume", "HJCMaximumCompression",
                             "Pressure", "VonMisesStress", "HJCAcousticModulus"})
        particles->addEvolvingVariable<Real>(name);
    particles->addEvolvingVariable<int>("HJCIntegrationStatus");
}

Matd HJCSolid::UpdateStress(const Matd &F, size_t i, Real dt)
{
    HJCState state{stress_[i], damage_[i], plastic_strain_[i], plastic_volume_[i], maximum_compression_[i]};
    HJCIntegrationStatus status = ConstitutiveModel().UpdateStress(F, previous_deformation_[i], dt, state);
    dv_status_->Data()[i] = static_cast<int>(status);
    if (status == HJCIntegrationStatus::invalid_deformation)
        throw std::domain_error("HJC requires positive finite deformation determinant");
    if (status == HJCIntegrationStatus::invalid_increment)
        throw std::invalid_argument("HJC requires finite increments, positive J and a positive timestep");
    if (status == HJCIntegrationStatus::singular_increment)
        throw std::domain_error("HJC incremental stretch is singular");
    if (status != HJCIntegrationStatus::success)
        throw std::domain_error("HJC stress update produced a non-finite state");
    if (dt == 0)
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
    Real pressure = -state.stress.trace() / Real(3);
    Real equivalent_stress = math::sqrt(Real(1.5) * (state.stress + pressure * Mat3d::Identity()).squaredNorm());
    Real modulus = AcousticModulus(math::max(Real(1) / F.determinant() - Real(1), state.maximum_compression));
    if (!std::isfinite(pressure) || !std::isfinite(equivalent_stress) || !std::isfinite(modulus))
    {
        dv_status_->Data()[i] = static_cast<int>(HJCIntegrationStatus::nonfinite_state);
        throw std::domain_error("HJC stress update produced a non-finite state");
    }
    stress_[i] = state.stress;
    previous_deformation_[i] = F;
    damage_[i] = state.damage;
    plastic_strain_[i] = state.plastic_strain;
    plastic_volume_[i] = state.plastic_volume;
    maximum_compression_[i] = state.maximum_compression;
    pressure_[i] = pressure;
    equivalent_stress_[i] = equivalent_stress;
    acoustic_modulus_[i] = modulus;
    return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
}

Matd HJCSolid::StressPK1(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Matd HJCSolid::StressPK2(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Matd HJCSolid::StressCauchy(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Real HJCSolid::VolumetricKirchhoff(Real) { throw std::logic_error("HJC pressure requires material history"); }
} // namespace SPH
