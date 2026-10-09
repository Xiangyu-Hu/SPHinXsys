// SPDX-License-Identifier: Apache-2.0
#ifndef HJC_SOLID_HPP
#define HJC_SOLID_HPP

#include "hjc_solid.h"
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace SPH
{
inline Real HJCConstitutiveModel::DensePressure(Real compression) const
{
    const auto &p = parameters_;
    Real mu = (compression - p.lock_strain) / (Real(1) + p.lock_strain);
    return mu < Real(0) ? p.K1 * mu : mu * (p.K1 + mu * (p.K2 + mu * p.K3));
}

inline Real HJCConstitutiveModel::Pressure(Real compression, HJCState &state) const
{
    const auto &p = parameters_;
    state.maximum_compression = math::max(state.maximum_compression, compression);
    Real peak = math::min(state.maximum_compression, lock_compression_);
    Real fraction = math::min(Real(1), math::max(Real(0),
                                                 (peak - p.crush_strain) / (lock_compression_ - p.crush_strain)));
    state.plastic_volume = p.lock_strain * fraction;
    if (state.maximum_compression >= lock_compression_)
        return DensePressure(compression);
    if (peak <= p.crush_strain)
        return K0_ * compression;
    Real peak_pressure = p.crush_pressure + transition_slope_ * (peak - p.crush_strain);
    Real unloading_modulus = peak_pressure / (peak - state.plastic_volume);
    return unloading_modulus * (compression - state.plastic_volume);
}

inline Real HJCConstitutiveModel::FractureStrain(Real pressure) const
{
    const auto &p = parameters_;
    return math::max(p.minimum_fracture_strain,
                     p.D1 * math::pow(math::max(Real(0), (pressure + p.tensile_strength) / p.compressive_strength), p.D2));
}

inline Real HJCConstitutiveModel::YieldStrength(Real pressure, Real damage, Real strain_rate) const
{
    const auto &p = parameters_;
    Real strength = pressure < Real(0)
                        ? p.A * math::max(Real(0), Real(1) - damage + pressure / p.tensile_strength)
                        : p.A * (Real(1) - damage) + p.B * math::pow(pressure / p.compressive_strength, p.N);
    // The reference rate defines the quasi-static floor.
    Real rate_factor = Real(1) + p.C * math::log(math::max(Real(1), strain_rate / p.reference_strain_rate));
    return p.compressive_strength * math::min(p.maximum_strength, strength * rate_factor);
}

inline HJCIntegrationStatus HJCConstitutiveModel::Integrate(
    const Mat3d &increment, Real compression, Real dt, HJCState &state) const
{
    if (!(dt > Real(0) && math::isfinite(dt)) || !increment.allFinite() ||
        !(compression > Real(-1) && math::isfinite(compression)))
        return HJCIntegrationStatus::invalid_increment;
    const auto &p = parameters_;
    HJCState next = state;
    Mat3d dev = increment - increment.trace() / Real(3) * Mat3d::Identity();
    Mat3d trial = state.stress - state.stress.trace() / Real(3) * Mat3d::Identity() + Real(2) * G0_ * dev;
    Real q = math::sqrt(Real(1.5) * trial.squaredNorm());
    Real old_plastic_volume = state.plastic_volume;
    Real pressure = math::max(Pressure(compression, next), -p.tensile_strength * (Real(1) - state.damage));
    Real rate = math::sqrt(Real(2) / Real(3) * dev.squaredNorm()) / dt;
    Real strength = YieldStrength(pressure, state.damage, rate);
    Real scale = q > strength ? strength / q : Real(1);
    Real plastic_increment = (Real(1) - scale) * q / (Real(3) * G0_);
    next.plastic_strain += plastic_increment;
    next.damage = math::min(Real(1), state.damage +
                                         (plastic_increment + next.plastic_volume - old_plastic_volume) / FractureStrain(pressure));
    // Damage is explicit in the deviatoric return; enforce the current tensile cutoff.
    pressure = math::max(pressure, -p.tensile_strength * (Real(1) - next.damage));
    next.stress = scale * trial - pressure * Mat3d::Identity();
    if (!next.stress.allFinite() || !math::isfinite(next.damage) || !math::isfinite(next.plastic_strain) ||
        !math::isfinite(next.plastic_volume) || !math::isfinite(next.maximum_compression))
        return HJCIntegrationStatus::nonfinite_state;
    state = next;
    return HJCIntegrationStatus::success;
}

inline Real HJCConstitutiveModel::AcousticModulus(Real compression) const
{
    const auto &p = parameters_;
    Real eta = math::max(Real(0), (compression - p.lock_strain) / (Real(1) + p.lock_strain));
    Real dense_tangent = (p.K1 + Real(2) * p.K2 * eta + Real(3) * p.K3 * eta * eta) / (Real(1) + p.lock_strain);
    return (Real(1) + math::max(compression, Real(0))) * math::max(K0_, dense_tangent) + Real(4) * G0_ / Real(3);
}

inline HJCIntegrationStatus HJCConstitutiveModel::UpdateStress(
    const Matd &F, const Matd &previous_deformation, Real dt, HJCState &state) const
{
    Real J = F.determinant();
    if (!(J > Real(0)) || !math::isfinite(J) || !F.allFinite())
        return HJCIntegrationStatus::invalid_deformation;
    if (dt == Real(0))
        return HJCIntegrationStatus::success;
    if (!(dt > Real(0)) || !math::isfinite(dt))
        return HJCIntegrationStatus::invalid_increment;
    Real previous_J = previous_deformation.determinant();
    if (!(previous_J > Real(0)) || !math::isfinite(previous_J) || !previous_deformation.allFinite())
        return HJCIntegrationStatus::singular_increment;
    Mat3d delta = Mat3d::Identity();
    delta.template topLeftCorner<Dimensions, Dimensions>() = F * previous_deformation.inverse();
    if (!delta.allFinite())
        return HJCIntegrationStatus::singular_increment;
    Eigen::SelfAdjointEigenSolver<Mat3d> solver(delta.transpose() * delta);
    if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite() || solver.eigenvalues().minCoeff() <= Real(0))
        return HJCIntegrationStatus::singular_increment;
    Mat3d U = solver.eigenvectors();
    Mat3d rotation = delta * U * solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
    Mat3d increment = U * (Real(0.5) * solver.eigenvalues().array().log()).matrix().asDiagonal() * U.transpose();
    HJCState next = state;
    HJCIntegrationStatus status = Integrate(increment, Real(1) / J - Real(1), dt, next);
    if (status != HJCIntegrationStatus::success)
        return status;
    next.stress = rotation * next.stress * rotation.transpose();
    if (!next.stress.allFinite())
        return HJCIntegrationStatus::nonfinite_state;
    state = next;
    return HJCIntegrationStatus::success;
}

template <typename ExecutionPolicy>
HJCSolid::ConstituteKernel::ConstituteKernel(const ExecutionPolicy &ex_policy, HJCSolid &encloser)
    : model_(encloser.ConstitutiveModel()),
      stress_(encloser.dv_stress_->DelegatedData(ex_policy)),
      previous_deformation_(encloser.dv_previous_deformation_->DelegatedData(ex_policy)),
      damage_(encloser.dv_damage_->DelegatedData(ex_policy)),
      plastic_strain_(encloser.dv_plastic_strain_->DelegatedData(ex_policy)),
      plastic_volume_(encloser.dv_plastic_volume_->DelegatedData(ex_policy)),
      maximum_compression_(encloser.dv_maximum_compression_->DelegatedData(ex_policy)),
      pressure_(encloser.dv_pressure_->DelegatedData(ex_policy)),
      equivalent_stress_(encloser.dv_equivalent_stress_->DelegatedData(ex_policy)),
      acoustic_modulus_(encloser.dv_acoustic_modulus_->DelegatedData(ex_policy)),
      status_(encloser.dv_status_->DelegatedData(ex_policy)) {}

inline Matd HJCSolid::ConstituteKernel::UpdateStress(const Matd &F, size_t i, Real dt)
{
    HJCState state{stress_[i], damage_[i], plastic_strain_[i], plastic_volume_[i], maximum_compression_[i]};
    HJCIntegrationStatus status = model_.UpdateStress(F, previous_deformation_[i], dt, state);
    status_[i] = static_cast<int>(status);
    if (status != HJCIntegrationStatus::success || dt == Real(0))
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
    Real pressure = -state.stress.trace() / Real(3);
    Real equivalent_stress = math::sqrt(Real(1.5) * (state.stress + pressure * Mat3d::Identity()).squaredNorm());
    Real modulus = model_.AcousticModulus(math::max(Real(1) / F.determinant() - Real(1), state.maximum_compression));
    if (!math::isfinite(pressure) || !math::isfinite(equivalent_stress) || !math::isfinite(modulus))
    {
        status_[i] = static_cast<int>(HJCIntegrationStatus::nonfinite_state);
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
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
} // namespace SPH
#endif // HJC_SOLID_HPP
