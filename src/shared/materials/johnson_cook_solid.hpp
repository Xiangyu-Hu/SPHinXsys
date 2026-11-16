#ifndef JOHNSON_COOK_SOLID_HPP
#define JOHNSON_COOK_SOLID_HPP
#include "base_particles.hpp"
#include "johnson_cook_solid.h"
#include <Eigen/Eigenvalues>

namespace SPH
{
//=================================================================================================//
inline Real JohnsonCookConstitutiveModel::NormalizedTemperature(Real temperature) const
{
    return math::min(Real(1), math::max(Real(0), (temperature - parameters_.reference_temperature) * inverse_temperature_interval_));
}
//=================================================================================================//
inline Real JohnsonCookConstitutiveModel::FractureStrain(Real mean_stress_ratio, Real strain_rate, Real temperature) const
{
    const auto &p = parameters_;
    return (p.failure_D1 + p.failure_D2 * math::exp(p.failure_D3 * mean_stress_ratio)) *
           (1 + p.failure_D4 * math::log(math::max(Real(1), strain_rate / p.reference_strain_rate))) *
           (1 + p.failure_D5 * NormalizedTemperature(temperature));
}
//=================================================================================================//
inline JohnsonCookIntegrationStatus JohnsonCookConstitutiveModel::Integrate(const Mat3d &strain_increment, Real J, Real dt, JohnsonCookState &state) const
{
    if (!(dt > 0 && J > 0 && math::isfinite(dt) && math::isfinite(J) && strain_increment.allFinite()))
        return JohnsonCookIntegrationStatus::invalid_increment;
    if (!state.active)
    {
        state.stress.setZero();
        state.equivalent_stress = state.mean_stress = 0;
        return JohnsonCookIntegrationStatus::success;
    }
    JohnsonCookState next = state;
    const auto &p = parameters_;
    const Mat3d identity = Mat3d::Identity();
    const Mat3d deviatoric_increment = strain_increment - strain_increment.trace() * identity / 3;
    Mat3d trial_stress = state.stress - state.stress.trace() * identity / 3 + 2 * G0_ * deviatoric_increment;
    const Mat3d trial_deviator = trial_stress - trial_stress.trace() * identity / 3;
    const Real trial_equivalent_stress = math::sqrt(Real(1.5) * trial_deviator.squaredNorm());
    const Real old_temperature = state.temperature;
    const Real rate_factor = 1 + p.strain_rate_coefficient *
                                     math::log(math::max(Real(1), state.plastic_strain_rate / p.reference_strain_rate));
    const Real thermal_factor = 1 - math::pow(NormalizedTemperature(old_temperature), p.thermal_softening_exponent);
    const Real factor = rate_factor * thermal_factor;
    const Real yield_stress = (p.yield_stress + p.hardening_coefficient * math::pow(state.plastic_strain, p.hardening_exponent)) * factor;
    Real hardening_tangent = 0;
    if (state.plastic_strain > 0)
        hardening_tangent = p.hardening_coefficient * p.hardening_exponent *
                            math::pow(state.plastic_strain, p.hardening_exponent - 1) * factor;
    if (p.hardening_exponent == 1)
        hardening_tangent = p.hardening_coefficient * factor;
    Real scale = math::min(Real(1), yield_stress / math::max(trial_equivalent_stress, Real(1e-15)));
    const Real plastic_increment = (1 - scale) * trial_equivalent_stress / math::max(3 * G0_ + hardening_tangent, Real(1e-15));
    const Real updated_yield_stress = yield_stress + plastic_increment * hardening_tangent;
    scale = math::min(Real(1), updated_yield_stress / math::max(trial_equivalent_stress, Real(1e-15)));
    const Real pressure = K0_ * (1 / J - 1);
    next.stress = scale * trial_stress - pressure * identity;
    next.plastic_strain += plastic_increment;
    next.plastic_strain_rate = plastic_increment / dt;
    next.mean_stress = -pressure;
    const Mat3d deviatoric_stress = next.stress - next.stress.trace() * identity / 3;
    next.equivalent_stress = math::sqrt(Real(1.5) * deviatoric_stress.squaredNorm());
    const Real plastic_work = math::max(yield_stress, updated_yield_stress) * plastic_increment;
    next.plastic_work_specific += plastic_work / rho0_;
    next.temperature += plastic_work * inverse_heat_capacity_;
    if (p.failure_enabled && plastic_increment > 0)
    {
        const Real failure_rate = math::sqrt(Real(2.0 / 3.0) * deviatoric_increment.squaredNorm()) / dt;
        const Real failure_strain = FractureStrain(next.mean_stress / math::max(next.equivalent_stress, Real(1e-20)), failure_rate, old_temperature);
        if (failure_strain > 0)
            next.damage = math::min(Real(1), next.damage + plastic_increment / failure_strain);
    }
    if (next.damage >= 1)
    {
        next.active = false;
        next.stress.setZero();
        next.equivalent_stress = next.mean_stress = 0;
    }
    if (!next.stress.allFinite() || !math::isfinite(next.plastic_strain) || !math::isfinite(next.plastic_strain_rate) ||
        !math::isfinite(next.temperature) || !math::isfinite(next.damage) || !math::isfinite(next.plastic_work_specific) ||
        !math::isfinite(next.equivalent_stress) || !math::isfinite(next.mean_stress))
        return JohnsonCookIntegrationStatus::nonfinite_state;
    state = next;
    return JohnsonCookIntegrationStatus::success;
}
//=================================================================================================//
inline JohnsonCookIntegrationStatus JohnsonCookConstitutiveModel::UpdateStress(
    const Matd &F, const Matd &previous_F, Real dt, JohnsonCookState &state) const
{
    if (!state.active || dt == 0)
        return JohnsonCookIntegrationStatus::success;
    if (!F.allFinite() || !(F.determinant() > 0))
        return JohnsonCookIntegrationStatus::invalid_deformation;
    if (!previous_F.allFinite() || !(previous_F.determinant() > 0))
        return JohnsonCookIntegrationStatus::singular_increment;
    Mat3d delta = Mat3d::Identity();
    delta.template topLeftCorner<Dimensions, Dimensions>() = F * previous_F.inverse();
    if (!delta.allFinite())
        return JohnsonCookIntegrationStatus::singular_increment;
    Eigen::SelfAdjointEigenSolver<Mat3d> solver(delta.transpose() * delta);
    if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite() || solver.eigenvalues().minCoeff() <= 0)
        return JohnsonCookIntegrationStatus::singular_increment;
    const Mat3d &U = solver.eigenvectors();
    const Mat3d rotation = delta * U * solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
    const Mat3d increment = U * (Real(0.5) * solver.eigenvalues().array().log()).matrix().asDiagonal() * U.transpose();
    JohnsonCookState next = state;
    const auto status = Integrate(increment, F.determinant(), dt, next);
    if (status != JohnsonCookIntegrationStatus::success)
        return status;
    next.stress = rotation * next.stress * rotation.transpose();
    if (!next.stress.allFinite())
        return JohnsonCookIntegrationStatus::nonfinite_state;
    state = next;
    return JohnsonCookIntegrationStatus::success;
}
//=================================================================================================//
template <class ExecutionPolicy>
JohnsonCookSolid::ConstituteKernel::ConstituteKernel(const ExecutionPolicy &ex_policy, JohnsonCookSolid &encloser)
    : model_(encloser.ConstitutiveModel()),
      stress_(encloser.particles_->getVariableByName<Mat3d>("StressCauchy")->DelegatedData(ex_policy)),
      previous_deformation_(encloser.particles_->getVariableByName<Matd>("JCPreviousDeformation")->DelegatedData(ex_policy)),
      previous_time_step_(encloser.particles_->getVariableByName<Real>("JCPreviousTimeStep")->DelegatedData(ex_policy)),
      plastic_strain_(encloser.particles_->getVariableByName<Real>("JCPlasticStrain")->DelegatedData(ex_policy)),
      plastic_strain_rate_(encloser.particles_->getVariableByName<Real>("JCPlasticStrainRate")->DelegatedData(ex_policy)),
      temperature_(encloser.particles_->getVariableByName<Real>("JCTemperature")->DelegatedData(ex_policy)),
      damage_(encloser.particles_->getVariableByName<Real>("JCDamage")->DelegatedData(ex_policy)),
      equivalent_stress_(encloser.particles_->getVariableByName<Real>("VonMisesStress")->DelegatedData(ex_policy)),
      mean_stress_(encloser.particles_->getVariableByName<Real>("MeanStress")->DelegatedData(ex_policy)),
      plastic_work_specific_(encloser.particles_->getVariableByName<Real>("PlasticWorkSpecific")->DelegatedData(ex_policy)),
      material_active_(encloser.particles_->getVariableByName<int>("MaterialActive")->DelegatedData(ex_policy)),
      status_(encloser.particles_->getVariableByName<int>("JCIntegrationStatus")->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Matd JohnsonCookSolid::ConstituteKernel::UpdateStress(const Matd &F, size_t i, Real dt)
{
    status_[i] = 0;
    if (!material_active_[i])
        return Matd::Zero();
    if (dt == 0)
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
    JohnsonCookState state;
    state.stress = stress_[i];
    state.plastic_strain = plastic_strain_[i];
    state.plastic_strain_rate = plastic_strain_rate_[i];
    state.temperature = temperature_[i];
    state.damage = damage_[i];
    state.plastic_work_specific = plastic_work_specific_[i];
    const auto status = model_.UpdateStress(F, previous_deformation_[i], Real(0.5) * (previous_time_step_[i] + dt), state);
    status_[i] = static_cast<int>(status);
    if (status != JohnsonCookIntegrationStatus::success)
        return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
    stress_[i] = state.stress;
    previous_deformation_[i] = F;
    previous_time_step_[i] = dt;
    plastic_strain_[i] = state.plastic_strain;
    plastic_strain_rate_[i] = state.plastic_strain_rate;
    temperature_[i] = state.temperature;
    damage_[i] = state.damage;
    equivalent_stress_[i] = state.equivalent_stress;
    mean_stress_[i] = state.mean_stress;
    plastic_work_specific_[i] = state.plastic_work_specific;
    material_active_[i] = state.active ? 1 : 0;
    return stress_[i].template topLeftCorner<Dimensions, Dimensions>();
}
} // namespace SPH
#endif // JOHNSON_COOK_SOLID_HPP
