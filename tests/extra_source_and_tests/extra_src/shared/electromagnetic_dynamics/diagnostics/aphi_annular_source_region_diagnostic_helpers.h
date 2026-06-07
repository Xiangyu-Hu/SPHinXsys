#ifndef APHI_ANNULAR_SOURCE_REGION_DIAGNOSTIC_HELPERS_H
#define APHI_ANNULAR_SOURCE_REGION_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

#include <cmath>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiAnnularSourceRegionSpec
{
    Real dp = 0.1;
    Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    Vecd annulus_center = Vecd(0.5, 0.5, 0.5);
    Real annulus_inner_radius = 0.12;
    Real annulus_outer_radius = 0.18;
    Real annulus_z_half_height = 0.12;
    bool write_vtp = true;
};

/** TEAM7 coil-slot source diagnostic (sigma=0, RHS-only). Annulus fields are metadata only. */
inline benchmark::AphiTeam7LikeUnitBoxLayout buildAnnularSourceLayout(const AphiAnnularSourceRegionSpec &spec)
{
    benchmark::AphiTeam7LikeUnitBoxLayout layout = benchmark::buildTeam7LayoutForBox(
        spec.body_length, spec.body_height, spec.body_width);
    layout.coil_material.sigma = 0.0;
    (void)spec.annulus_center;
    (void)spec.annulus_inner_radius;
    (void)spec.annulus_outer_radius;
    (void)spec.annulus_z_half_height;
    return layout;
}

inline AphiSourceDrivenEmSolveMetrics runAnnularSourceRegionDiagnostic(int ac, char *av[],
                                                                       const AphiAnnularSourceRegionSpec &spec)
{
    AphiSourceDrivenEmSolveSpec solve_spec;
    solve_spec.dp = spec.dp;
    solve_spec.body_length = spec.body_length;
    solve_spec.body_height = spec.body_height;
    solve_spec.body_width = spec.body_width;
    solve_spec.write_vtp = spec.write_vtp;
    solve_spec.write_probe_csv = false;
    solve_spec.passive_air_shell_width_scale = 0.0;

    const benchmark::AphiTeam7LikeUnitBoxLayout layout = buildAnnularSourceLayout(spec);
    return runSourceDrivenEmSolveWithLayout(ac, av, solve_spec, layout);
}

inline bool annularSourceRegionDiagnosticPassed(const AphiSourceDrivenEmSolveMetrics &metrics)
{
    const Real air_ratio = metrics.air_Joule_integral / (metrics.conductor_Joule_integral + TinyReal);
    return metrics.converged && metrics.finite_field_check && metrics.particle_count_source > 0 &&
           metrics.particle_count_conductor > 0 && metrics.source_rhs_l2 > 0.0 && metrics.conductor_Joule_integral > 1.0e-6 &&
           metrics.max_abs_J > 1.0e-8 && air_ratio < 0.05;
}

inline void printAnnularSourceRegionDiagnosticMetrics(const char *test_name, const AphiSourceDrivenEmSolveMetrics &metrics,
                                                     bool passed)
{
    printSourceDrivenEmSolveMetrics(test_name, metrics, passed);
    std::cout << test_name
              << " source_region_model=team7_coil_slot_only annular_geometry_implemented=0"
              << " naming_note=legacy_annular_diagnostic" << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
