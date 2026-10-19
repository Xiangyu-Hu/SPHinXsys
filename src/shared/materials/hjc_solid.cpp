// SPDX-License-Identifier: Apache-2.0
#include "hjc_solid.h"
#include "base_particles.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SPH
{
HJCSolid::HJCSolid(Real rho0, const HJCParameters &p)
    : ElasticSolid(rho0), parameters_(p), stress_(nullptr), previous_deformation_(nullptr),
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
    const auto &p = parameters_;
    Real mu = (compression - p.lock_strain) / (1 + p.lock_strain);
    return mu < 0 ? p.K1 * mu : mu * (p.K1 + mu * (p.K2 + mu * p.K3));
}

Real HJCSolid::Pressure(Real compression, HJCState &state) const
{
    const auto &p = parameters_;
    state.maximum_compression = std::max(state.maximum_compression, compression);
    Real peak = std::min(state.maximum_compression, lock_compression_);
    Real fraction = std::clamp((peak - p.crush_strain) / (lock_compression_ - p.crush_strain), Real(0), Real(1));
    state.plastic_volume = p.lock_strain * fraction;
    if (state.maximum_compression >= lock_compression_)
        return DensePressure(compression);
    if (peak <= p.crush_strain)
        return K0_ * compression;
    Real peak_pressure = p.crush_pressure + transition_slope_ * (peak - p.crush_strain);
    Real unloading_modulus = peak_pressure / (peak - state.plastic_volume);
    return unloading_modulus * (compression - state.plastic_volume);
}

Real HJCSolid::FractureStrain(Real pressure) const
{
    const auto &p = parameters_;
    return std::max(p.minimum_fracture_strain,
                    p.D1 * std::pow(std::max(Real(0), (pressure + p.tensile_strength) / p.compressive_strength), p.D2));
}

Real HJCSolid::YieldStrength(Real pressure, Real damage, Real strain_rate) const
{
    const auto &p = parameters_;
    Real strength = pressure < 0
                        ? p.A * std::max(Real(0), 1 - damage + pressure / p.tensile_strength)
                        : p.A * (1 - damage) + p.B * std::pow(pressure / p.compressive_strength, p.N);
    // The reference rate defines the quasi-static floor.
    Real rate_factor = 1 + p.C * std::log(std::max(Real(1), strain_rate / p.reference_strain_rate));
    return p.compressive_strength * std::min(p.maximum_strength, strength * rate_factor);
}

void HJCSolid::Integrate(const Mat3d &increment, Real compression, Real dt, HJCState &state) const
{
    if (!(dt > 0 && std::isfinite(dt)) || !increment.allFinite() ||
        !(compression > -1 && std::isfinite(compression)))
        throw std::invalid_argument("HJC requires finite increments, positive J and a positive timestep");
    const auto &p = parameters_;
    Mat3d dev = increment - increment.trace() / 3 * Mat3d::Identity();
    Mat3d trial = state.stress - state.stress.trace() / 3 * Mat3d::Identity() + 2 * G0_ * dev;
    Real q = std::sqrt(1.5 * trial.squaredNorm());
    Real old_plastic_volume = state.plastic_volume;
    Real pressure = std::max(Pressure(compression, state), -p.tensile_strength * (1 - state.damage));
    Real rate = std::sqrt(Real(2) / 3 * dev.squaredNorm()) / dt;
    Real strength = YieldStrength(pressure, state.damage, rate);
    Real scale = q > strength ? strength / q : 1;
    Real plastic_increment = (1 - scale) * q / (3 * G0_);
    state.plastic_strain += plastic_increment;
    state.damage = std::min(Real(1), state.damage +
        (plastic_increment + state.plastic_volume - old_plastic_volume) / FractureStrain(pressure));
    // Damage is explicit in the deviatoric return; enforce the current tensile cutoff.
    pressure = std::max(pressure, -p.tensile_strength * (1 - state.damage));
    state.stress = scale * trial - pressure * Mat3d::Identity();
}

Real HJCSolid::AcousticModulus(Real compression) const
{
    const auto &p = parameters_;
    Real eta = std::max(Real(0), (compression - p.lock_strain) / (1 + p.lock_strain));
    Real dense_tangent = (p.K1 + 2 * p.K2 * eta + 3 * p.K3 * eta * eta) / (1 + p.lock_strain);
    return (1 + std::max(compression, Real(0))) * std::max(K0_, dense_tangent) + 4 * G0_ / 3;
}

void HJCSolid::initializeLocalParameters(BaseParticles *particles)
{
    ElasticSolid::initializeLocalParameters(particles);
    stress_ = particles->registerStateVariableData<Mat3d>("StressCauchy");
    previous_deformation_ = particles->registerStateVariableData<Matd>("HJCPreviousDeformation", IdentityMatrix<Matd>::value);
    damage_ = particles->registerStateVariableData<Real>("HJCDamage");
    plastic_strain_ = particles->registerStateVariableData<Real>("HJCPlasticStrain");
    plastic_volume_ = particles->registerStateVariableData<Real>("HJCPlasticVolume");
    maximum_compression_ = particles->registerStateVariableData<Real>("HJCMaximumCompression");
    pressure_ = particles->registerStateVariableData<Real>("Pressure");
    equivalent_stress_ = particles->registerStateVariableData<Real>("VonMisesStress");
    acoustic_modulus_ = particles->registerStateVariableData<Real>("HJCAcousticModulus", AcousticModulus(0));
    particles->addEvolvingVariable<Mat3d>("StressCauchy");
    particles->addEvolvingVariable<Matd>("HJCPreviousDeformation");
    for (const char *name : {"HJCDamage", "HJCPlasticStrain", "HJCPlasticVolume", "HJCMaximumCompression",
                                    "Pressure", "VonMisesStress", "HJCAcousticModulus"})
        particles->addEvolvingVariable<Real>(name);
}

Matd HJCSolid::UpdateStress(const Matd &F, size_t i, Real dt)
{
    Real J = F.determinant();
    if (!(J > 0) || !std::isfinite(J))
        throw std::domain_error("HJC requires positive finite deformation determinant");
    if (dt == 0)
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
    Mat3d delta = Mat3d::Identity();
    delta.template topLeftCorner<Dimensions, Dimensions>() = F * previous_deformation_[i].inverse();
    Eigen::SelfAdjointEigenSolver<Mat3d> solver(delta.transpose() * delta);
    if (solver.info() != Eigen::Success || solver.eigenvalues().minCoeff() <= 0)
        throw std::domain_error("HJC incremental stretch is singular");
    Mat3d U = solver.eigenvectors();
    Mat3d rotation = delta * U * solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
    Mat3d increment = U * (Real(0.5) * solver.eigenvalues().array().log()).matrix().asDiagonal() * U.transpose();
    HJCState state{stress_[i], damage_[i], plastic_strain_[i], plastic_volume_[i], maximum_compression_[i]};
    Integrate(increment, 1 / J - 1, dt, state);
    stress_[i] = rotation * state.stress * rotation.transpose();
    previous_deformation_[i] = F;
    damage_[i] = state.damage;
    plastic_strain_[i] = state.plastic_strain;
    plastic_volume_[i] = state.plastic_volume;
    maximum_compression_[i] = state.maximum_compression;
    pressure_[i] = -stress_[i].trace() / 3;
    equivalent_stress_[i] = std::sqrt(1.5 * (stress_[i] + pressure_[i] * Mat3d::Identity()).squaredNorm());
    acoustic_modulus_[i] = AcousticModulus(std::max(1 / J - 1, state.maximum_compression));
    return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
}

Matd HJCSolid::StressPK1(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Matd HJCSolid::StressPK2(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Matd HJCSolid::StressCauchy(Matd &, size_t) { throw std::logic_error("Use HJCIntegration1stHalf"); }
Real HJCSolid::VolumetricKirchhoff(Real) { throw std::logic_error("HJC pressure requires material history"); }
} // namespace SPH
