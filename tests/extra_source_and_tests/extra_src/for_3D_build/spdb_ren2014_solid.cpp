#include "spdb_ren2014_solid.h"

#include "base_particles.hpp"

#include <algorithm>
#include <vector>

namespace SPH
{
//=================================================================================================//
Ren2014X20Cr13Solid::Ren2014X20Cr13Solid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real temperature_K,
                                         bool inelastic_enabled)
    : NeoHookeanSolid(rho0, youngs_modulus, poisson_ratio), dv_inverse_plastic_strain_(nullptr),
      dv_inelastic_eq_(nullptr), dv_range_flag_(nullptr), temperature_K_(temperature_K),
      inelastic_enabled_(inelastic_enabled), inverse_plastic_strain_(nullptr), inelastic_eq_(nullptr),
      range_flag_(nullptr)
{
    material_type_name_ = "Ren2014X20Cr13Solid";
}
//=================================================================================================//
void Ren2014X20Cr13Solid::initializeLocalParameters(BaseParticles *base_particles)
{
    NeoHookeanSolid::initializeLocalParameters(base_particles);
    dv_inverse_plastic_strain_ = base_particles->registerStateVariable<Matd>(
        "InversePlasticRightCauchyStrain",
        [&](size_t) -> Matd
        { return Matd::Identity(); });
    inverse_plastic_strain_ = dv_inverse_plastic_strain_->Data();
    dv_inelastic_eq_ = base_particles->registerStateVariable<Real>("EquivalentInelasticStrain");
    inelastic_eq_ = dv_inelastic_eq_->Data();
    dv_range_flag_ = base_particles->registerStateVariable<Real>("ArrheniusRangeFlag");
    range_flag_ = dv_range_flag_->Data();
    base_particles->addEvolvingVariable<Matd>(dv_inverse_plastic_strain_);
    base_particles->addEvolvingVariable<Real>(dv_inelastic_eq_);
    base_particles->addEvolvingVariable<Real>(dv_range_flag_);
}
//=================================================================================================//
Matd Ren2014X20Cr13Solid::ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt)
{
    return relax_ren2014_elastic_left_cauchy(deformation, index_i, dt, G0_, temperature_K_, inelastic_enabled_,
                                             inverse_plastic_strain_, inelastic_eq_, range_flag_);
}
//=================================================================================================//
Ren2014PhysicalCommit commit_ren2014_physical_increment(BaseParticles &particles, Real shear_modulus,
                                                        Real temperature_K, Real dt_physical, Real rate_scale)
{
    Ren2014PhysicalCommit report;
    const size_t n = particles.TotalRealParticles();
    report.particle_count = n;
    Matd *inverse_plastic = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *inelastic_eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *range_flag = particles.getVariableDataByName<Real>("ArrheniusRangeFlag");
    Matd *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    std::vector<Matd> saved_inverse(n);
    std::vector<Real> saved_eq(n);
    std::vector<Real> saved_flag(n);
    for (size_t i = 0; i < n; ++i)
    {
        saved_inverse[i] = inverse_plastic[i];
        saved_eq[i] = inelastic_eq[i];
        saved_flag[i] = range_flag[i];
    }

    Real sum_trial = 0.0;
    Real sum_relaxed = 0.0;
    Real sum_rate = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real J = deformation[i].determinant();
        const Real deformation_norm = (deformation[i] - Matd::Identity()).norm();
        const bool small_strain = deformation_norm <= 0.05 && ABS(J - Real(1.0)) <= Real(0.02);
        Real sigma_trial = 0.0;
        ren2014_cauchy_stress(deformation[i], inverse_plastic[i], shear_modulus, 0.0, &sigma_trial, nullptr);
        sum_trial += sigma_trial;
        if (!small_strain)
        {
            range_flag[i] = 1.0;
            ++report.large_strain;
            continue;
        }
        const Real eq_before = inelastic_eq[i];
        relax_ren2014_elastic_left_cauchy(deformation[i], i, dt_physical, shear_modulus, temperature_K, true,
                                          inverse_plastic, inelastic_eq, range_flag, rate_scale);
        Real sigma_relaxed = 0.0;
        ren2014_cauchy_stress(deformation[i], inverse_plastic[i], shear_modulus, 0.0, &sigma_relaxed, nullptr);
        sum_relaxed += sigma_relaxed;
        if (rate_scale > 0.0 && range_flag[i] > saved_flag[i] && !(inelastic_eq[i] > eq_before))
        {
            ++report.rejected;
            continue;
        }
        const Real rate = dt_physical > 0.0 ? (inelastic_eq[i] - eq_before) / dt_physical : 0.0;
        sum_rate += rate;
        ++report.updated;
    }

    report.mean_sigma_trial_Pa = n > 0 ? sum_trial / static_cast<Real>(n) : 0.0;
    report.mean_sigma_relaxed_Pa = n > 0 ? sum_relaxed / static_cast<Real>(n) : 0.0;
    report.mean_rate_per_s = n > 0 ? sum_rate / static_cast<Real>(n) : 0.0;
    if (report.large_strain > 0 || report.rejected > 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            inverse_plastic[i] = saved_inverse[i];
            inelastic_eq[i] = saved_eq[i];
            range_flag[i] = saved_flag[i];
        }
        report.accepted = 0;
        report.mean_inelastic = 0.0;
        report.max_inelastic = 0.0;
        return report;
    }

    Real sum_eq = 0.0;
    Real max_eq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum_eq += inelastic_eq[i];
        max_eq = std::max(max_eq, inelastic_eq[i]);
    }
    report.mean_inelastic = n > 0 ? sum_eq / static_cast<Real>(n) : 0.0;
    report.max_inelastic = max_eq;
    report.accepted = 1;
    return report;
}
//=================================================================================================//
} // namespace SPH
