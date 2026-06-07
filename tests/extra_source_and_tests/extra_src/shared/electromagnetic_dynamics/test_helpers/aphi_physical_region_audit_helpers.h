#ifndef APHI_PHYSICAL_REGION_AUDIT_HELPERS_H
#define APHI_PHYSICAL_REGION_AUDIT_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiBodyRegionAuditMetrics
{
    std::string body_name;
    size_t particle_count_total = 0;
    size_t particle_count_physical_box = 0;
    size_t particle_count_shell = 0;
    size_t particle_count_air = 0;
    size_t particle_count_conductor = 0;
    size_t particle_count_source = 0;
    Real sigma_min = 0.0;
    Real sigma_max = 0.0;
    Real nu_min = 0.0;
    Real nu_max = 0.0;
    Real joule_integral_total = 0.0;
    Real joule_integral_air = 0.0;
    Real joule_integral_conductor = 0.0;
    Real joule_integral_source = 0.0;
    Real joule_integral_shell = 0.0;
    Real max_J_conductor = 0.0;
    Real max_J_source = 0.0;
    Real max_Joule_conductor = 0.0;
    Real rhs_l2 = 0.0;
};

inline AphiBodyRegionAuditMetrics hostAuditBodyRegions(
    SPHBody &body, const std::string &body_name, const AphiVariableNames &names,
    const AphiJouleHeatingFieldNames &joule_names, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
    Real body_length, Real body_height, Real body_width, Real conductor_sigma)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, names.material.sigma);
    syncVariableToHost<Real>(particles, names.material.nu);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *sigma = particles.getVariableDataByName<Real>(names.material.sigma);
    const Real *nu = particles.getVariableDataByName<Real>(names.material.nu);

    AphiBodyRegionAuditMetrics metrics;
    metrics.body_name = body_name;
    metrics.particle_count_total = total_real_particles;
    metrics.rhs_l2 = hostBlockNorm(particles, names.rhs, total_real_particles);

    bool sigma_initialized = false;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!sigma_initialized)
        {
            metrics.sigma_min = sigma[i];
            metrics.sigma_max = sigma[i];
            metrics.nu_min = nu[i];
            metrics.nu_max = nu[i];
            sigma_initialized = true;
        }
        metrics.sigma_min = std::min(metrics.sigma_min, sigma[i]);
        metrics.sigma_max = std::max(metrics.sigma_max, sigma[i]);
        metrics.nu_min = std::min(metrics.nu_min, nu[i]);
        metrics.nu_max = std::max(metrics.nu_max, nu[i]);

        if (isInsidePassiveAirShellParticle(positions[i], body_length, body_height, body_width))
        {
            metrics.particle_count_shell += 1;
            continue;
        }
        metrics.particle_count_physical_box += 1;
        if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
        {
            metrics.particle_count_conductor += 1;
        }
        else if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Coil))
        {
            metrics.particle_count_source += 1;
        }
        else
        {
            metrics.particle_count_air += 1;
        }
    }

    const std::function<bool(const Vecd &)> in_physical_conductor = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, body_length, body_height, body_width,
                                             AphiBenchmarkMaterialRegion::Conductor);
    };
    const std::function<bool(const Vecd &)> in_physical_source = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, body_length, body_height, body_width,
                                             AphiBenchmarkMaterialRegion::Coil);
    };
    const std::function<bool(const Vecd &)> in_physical_air = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, body_length, body_height, body_width,
                                             AphiBenchmarkMaterialRegion::Air);
    };

    metrics.joule_integral_total = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, [](const Vecd &) { return true; }, joule_names.joule_heat_source);
    metrics.joule_integral_conductor = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_conductor, joule_names.joule_heat_source);
    metrics.joule_integral_source = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_source, joule_names.joule_heat_source);
    metrics.joule_integral_air = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_air, joule_names.joule_heat_source);
    const auto in_shell = [&](const Vecd &position) {
        return isInsidePassiveAirShellParticle(position, body_length, body_height, body_width);
    };
    metrics.joule_integral_shell = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_shell, joule_names.joule_heat_source);

    metrics.max_J_conductor =
        conductor_sigma * hostParticleRegionVecdMax(particles, joule_names.electric_field_a_real, positions,
                                                  total_real_particles, in_physical_conductor);
    metrics.max_J_source = layout.coil_material.sigma *
                           hostParticleRegionVecdMax(particles, joule_names.electric_field_a_real, positions,
                                                     total_real_particles, in_physical_source);
    metrics.max_Joule_conductor = hostParticleRegionScalarMax(particles, joule_names.joule_heat_source, positions,
                                                              total_real_particles, in_physical_conductor);
    return metrics;
}

inline void printBodyRegionAuditMetrics(const std::string &prefix, const AphiBodyRegionAuditMetrics &metrics)
{
    std::cout << prefix << " body=" << metrics.body_name << " particles=" << metrics.particle_count_total
              << " physical_box=" << metrics.particle_count_physical_box << " shell=" << metrics.particle_count_shell
              << " air=" << metrics.particle_count_air << " conductor=" << metrics.particle_count_conductor
              << " source=" << metrics.particle_count_source << " sigma_min=" << metrics.sigma_min
              << " sigma_max=" << metrics.sigma_max << " nu_min=" << metrics.nu_min << " nu_max=" << metrics.nu_max
              << " rhs_l2=" << metrics.rhs_l2 << " joule_total=" << metrics.joule_integral_total
              << " joule_air=" << metrics.joule_integral_air
              << " joule_conductor=" << metrics.joule_integral_conductor
              << " joule_source=" << metrics.joule_integral_source << " joule_shell=" << metrics.joule_integral_shell
              << " max_J_conductor=" << metrics.max_J_conductor
              << " max_Joule_conductor=" << metrics.max_Joule_conductor << std::endl;
}

struct AphiThreeBodyContactAuditSummary
{
    AphiBodyRegionAuditMetrics air{};
    AphiBodyRegionAuditMetrics coil{};
    AphiBodyRegionAuditMetrics plate{};
    Real total_source_rhs_l2 = 0.0;
    Real global_conductor_joule = 0.0;
    Real global_air_joule = 0.0;
    Real global_source_joule = 0.0;
    Real global_shell_joule = 0.0;
};

inline AphiThreeBodyContactAuditSummary hostThreeBodyContactAuditSummary(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
    const AphiJouleHeatingFieldNames &joule_names, Real body_length, Real body_height, Real body_width)
{
    const Real conductor_sigma = case_setup.layout.conductor_material.sigma;
    AphiThreeBodyContactAuditSummary summary;
    summary.air = hostAuditBodyRegions(case_setup.air_body, "air", names, joule_names, case_setup.layout, body_length,
                                       body_height, body_width, conductor_sigma);
    summary.coil = hostAuditBodyRegions(case_setup.coil_body, "coil", names, joule_names, case_setup.layout, body_length,
                                        body_height, body_width, conductor_sigma);
    summary.plate = hostAuditBodyRegions(case_setup.plate_body, "plate", names, joule_names, case_setup.layout,
                                         body_length, body_height, body_width, conductor_sigma);
    summary.total_source_rhs_l2 = summary.coil.rhs_l2;
    summary.global_conductor_joule = summary.plate.joule_integral_conductor;
    summary.global_air_joule = summary.air.joule_integral_air + summary.coil.joule_integral_air +
                               summary.plate.joule_integral_air;
    summary.global_source_joule = summary.coil.joule_integral_source;
    summary.global_shell_joule = summary.air.joule_integral_shell + summary.coil.joule_integral_shell +
                                 summary.plate.joule_integral_shell;
    return summary;
}

inline bool threeBodyContactAuditPassed(const AphiThreeBodyContactAuditSummary &summary,
                                        Real max_air_to_conductor_joule_ratio = 0.05)
{
    if (summary.air.particle_count_air == 0 || summary.coil.particle_count_source == 0 ||
        summary.plate.particle_count_conductor == 0)
    {
        return false;
    }
    if (summary.total_source_rhs_l2 <= 0.0 || summary.global_conductor_joule <= 1.0e-6)
    {
        return false;
    }
    const Real air_to_conductor = summary.global_air_joule / (summary.global_conductor_joule + TinyReal);
    if (air_to_conductor >= max_air_to_conductor_joule_ratio)
    {
        return false;
    }
    const Real max_source_joule = std::max(Real(1.0e-12), Real(1.0e-6) * summary.global_conductor_joule);
    if (summary.global_source_joule > max_source_joule)
    {
        return false;
    }
    return summary.coil.sigma_max <= 1.0e-14;
}

inline bool twoBodyContactAuditPassed(const AphiBodyRegionAuditMetrics &left_audit,
                                      const AphiBodyRegionAuditMetrics &right_audit,
                                      Real max_air_to_conductor_joule_ratio = 0.05)
{
    if (left_audit.rhs_l2 <= 0.0 || right_audit.particle_count_conductor == 0)
    {
        return false;
    }
    const Real conductor_joule = right_audit.joule_integral_conductor;
    if (conductor_joule <= 1.0e-6)
    {
        return false;
    }
    const Real air_joule = left_audit.joule_integral_air + right_audit.joule_integral_air;
    return air_joule / (conductor_joule + TinyReal) < max_air_to_conductor_joule_ratio;
}

inline void printThreeBodyContactAuditSummary(const std::string &prefix, const AphiThreeBodyContactAuditSummary &summary)
{
    printBodyRegionAuditMetrics(prefix, summary.air);
    printBodyRegionAuditMetrics(prefix, summary.coil);
    printBodyRegionAuditMetrics(prefix, summary.plate);
    std::cout << prefix << " total_source_rhs_l2=" << summary.total_source_rhs_l2
              << " global_conductor_joule=" << summary.global_conductor_joule
              << " global_air_joule=" << summary.global_air_joule
              << " global_source_joule=" << summary.global_source_joule
              << " global_shell_joule=" << summary.global_shell_joule << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PHYSICAL_REGION_AUDIT_HELPERS_H
