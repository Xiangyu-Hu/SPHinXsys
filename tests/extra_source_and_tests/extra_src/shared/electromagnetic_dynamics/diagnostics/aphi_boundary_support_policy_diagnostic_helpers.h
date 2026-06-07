#ifndef APHI_BOUNDARY_SUPPORT_POLICY_DIAGNOSTIC_HELPERS_H
#define APHI_BOUNDARY_SUPPORT_POLICY_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

enum class AphiBoundarySupportPolicy
{
    Baseline,
    EnlargedAir,
    PassiveAirShell,
    /** @deprecated use PassiveAirShell */
    DummyShell = PassiveAirShell,
};

inline const char *boundarySupportPolicyName(AphiBoundarySupportPolicy policy)
{
    switch (policy)
    {
    case AphiBoundarySupportPolicy::EnlargedAir:
        return "enlarged_air_domain_padding";
    case AphiBoundarySupportPolicy::PassiveAirShell:
        return "passive_air_shell";
    default:
        return "baseline";
    }
}

struct AphiBoundaryPolicyRow
{
    AphiBoundarySupportPolicy policy = AphiBoundarySupportPolicy::Baseline;
    Real boundary_width_scale = 1.0;
    Real passive_air_shell_width_scale = 0.0;
    AphiSourceDrivenEmSolveMetrics em{};
};

inline AphiBoundaryPolicyRow runBoundaryPolicyRow(int ac, char *av[], AphiBoundarySupportPolicy policy,
                                                  Real boundary_width_scale, Real passive_air_shell_width_scale = 0.0)
{
    AphiSourceDrivenEmSolveSpec spec;
    spec.write_vtp = false;
    spec.boundary_width_scale = boundary_width_scale;
    spec.passive_air_shell_width_scale = passive_air_shell_width_scale;

    AphiBoundaryPolicyRow row;
    row.policy = policy;
    row.boundary_width_scale = boundary_width_scale;
    row.passive_air_shell_width_scale = passive_air_shell_width_scale;
    row.em = runSourceDrivenEmSolve(ac, av, spec);
    return row;
}

inline bool boundarySupportPolicyDiagnosticPassed(const StdVec<AphiBoundaryPolicyRow> &rows)
{
    if (rows.empty())
    {
        return false;
    }
    auto findPolicyRow = [&](AphiBoundarySupportPolicy policy, Real width_scale) -> const AphiBoundaryPolicyRow * {
        for (const AphiBoundaryPolicyRow &row : rows)
        {
            if (row.policy == policy && std::abs(row.boundary_width_scale - width_scale) <= TinyReal)
            {
                return &row;
            }
        }
        return nullptr;
    };
    auto findPassiveRow = [&]() -> const AphiBoundaryPolicyRow * {
        for (const AphiBoundaryPolicyRow &row : rows)
        {
            if (row.policy == AphiBoundarySupportPolicy::PassiveAirShell)
            {
                return &row;
            }
        }
        return nullptr;
    };

    const AphiBoundaryPolicyRow *baseline_row = findPolicyRow(AphiBoundarySupportPolicy::Baseline, 1.0);
    const AphiBoundaryPolicyRow *enlarged_w2 = findPolicyRow(AphiBoundarySupportPolicy::EnlargedAir, 2.0);
    const AphiBoundaryPolicyRow *enlarged_w3 = findPolicyRow(AphiBoundarySupportPolicy::EnlargedAir, 3.0);
    const AphiBoundaryPolicyRow *passive_row = findPassiveRow();
    if (baseline_row == nullptr || enlarged_w2 == nullptr || enlarged_w3 == nullptr || passive_row == nullptr)
    {
        return false;
    }

    const AphiBoundaryPolicyRow &baseline = *baseline_row;
    if (!baseline.em.converged || !baseline.em.finite_field_check || baseline.em.conductor_Joule_integral <= 0.0 ||
        baseline.em.source_rhs_l2 <= 0.0 || baseline.em.particle_count_conductor == 0 ||
        baseline.em.particle_count_source == 0 || baseline.em.particle_count_air == 0)
    {
        return false;
    }
    const Real baseline_joule = baseline.em.conductor_Joule_integral;
    const Real air_ratio_baseline = baseline.em.air_Joule_integral / (baseline_joule + TinyReal);
    if (air_ratio_baseline > 0.05)
    {
        return false;
    }

    for (const AphiBoundaryPolicyRow *row_ptr : {baseline_row, enlarged_w2, enlarged_w3})
    {
        const AphiBoundaryPolicyRow &row = *row_ptr;
        if (!row.em.converged || !row.em.finite_field_check || row.em.conductor_Joule_integral <= 0.0)
        {
            return false;
        }
        const Real rel_change = std::abs(row.em.conductor_Joule_integral - baseline_joule) / (baseline_joule + TinyReal);
        if (rel_change > 0.02)
        {
            return false;
        }
    }

    const AphiBoundaryPolicyRow &passive_row_ref = *passive_row;
    if (!passive_row_ref.em.converged || !passive_row_ref.em.finite_field_check ||
        passive_row_ref.em.conductor_Joule_integral <= 0.0)
    {
        return false;
    }
    const Real passive_rel_change =
        std::abs(passive_row_ref.em.conductor_Joule_integral - baseline_joule) / (baseline_joule + TinyReal);
    if (passive_rel_change > 0.02)
    {
        return false;
    }
    if (passive_row_ref.em.shell_sigma_max > 1.0e-14)
    {
        return false;
    }
    const Real shell_joule_limit = std::max(1.0e-12, baseline_joule * 1.0e-6);
    if (passive_row_ref.em.shell_joule_integral > shell_joule_limit)
    {
        return false;
    }
    return true;
}

inline void printBoundaryPolicyRows(const char *test_name, const StdVec<AphiBoundaryPolicyRow> &rows, bool passed)
{
    for (const AphiBoundaryPolicyRow &row : rows)
    {
        const Real air_ratio = row.em.air_Joule_integral / (row.em.conductor_Joule_integral + TinyReal);
        std::cout << test_name << " policy=" << boundarySupportPolicyName(row.policy)
                  << " boundary_width_scale=" << row.boundary_width_scale
                  << " passive_air_shell_width_scale=" << row.passive_air_shell_width_scale
                  << " converged=" << (row.em.converged ? 1 : 0) << " gmres_iterations=" << row.em.num_iterations
                  << " final_residual=" << row.em.final_residual << " particles=" << row.em.particles
                  << " physical_box=" << row.em.particle_count_physical_box << " shell=" << row.em.particle_count_shell
                  << " conductor_particles=" << row.em.particle_count_conductor
                  << " source_particles=" << row.em.particle_count_source
                  << " air_particles=" << row.em.particle_count_air << " shell_sigma_max=" << row.em.shell_sigma_max
                  << " shell_joule=" << row.em.shell_joule_integral
                  << " conductor_Joule=" << row.em.conductor_Joule_integral << " air_Joule=" << row.em.air_Joule_integral
                  << " air_to_conductor_Joule=" << air_ratio << " source_rhs_l2=" << row.em.source_rhs_l2
                  << " max_abs_B=" << row.em.max_abs_B << " max_abs_E=" << row.em.max_abs_E
                  << " max_abs_J=" << row.em.max_abs_J << std::endl;
    }
    if (rows.size() >= 2)
    {
        const Real rel_change =
            std::abs(rows.back().em.conductor_Joule_integral - rows.front().em.conductor_Joule_integral) /
            (std::abs(rows.front().em.conductor_Joule_integral) + TinyReal);
        std::cout << test_name << " conductor_Joule_relative_change_last_vs_baseline=" << rel_change << std::endl;
    }
    std::cout << test_name << " default_boundary_policy=enlarged_air_domain_padding"
              << " passive_air_shell_status=diagnostic_only ghost_mirror_boundary=not_implemented"
              << " legacy_ghost_buffer_diva=decoupled_from_heating_solver passed=" << (passed ? 1 : 0) << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
