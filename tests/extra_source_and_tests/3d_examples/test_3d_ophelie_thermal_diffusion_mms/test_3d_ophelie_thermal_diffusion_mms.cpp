/**
 * @file test_3d_ophelie_thermal_diffusion_mms.cpp
 * @brief Uniform JouleHeat + isotropic diffusion + cold-wall Dirichlet on a box lattice.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
constexpr Real kRho = 2500.0;
constexpr Real kCp = 1200.0;
constexpr Real kK = 1.0;
constexpr Real kT0 = 300.0;
constexpr Real kQ = 5.0e4;
} // namespace

int main(int, char *[])
{
    const Real dp = 0.04;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const Real dt = 0.1;
    const size_t n_steps = 5;

    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));
    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    Real *joule_heat = particles.getVariableDataByName<Real>(glass_names.joule_heat);
    for (size_t i = 0; i < n; ++i)
    {
        joule_heat[i] = kQ;
    }
    syncVariableToDevice<Real>(particles, glass_names.joule_heat);

    registerOphelieJouleHeatTemperatureField(particles, kT0);
    registerOphelieThermalDiffusionAuxFields(particles, kK);

    OphelieThermalDiffusionOneWayOptions thermal_options;
    thermal_options.enable_diffusion = true;
    thermal_options.enable_cold_wall_dirichlet = true;
    thermal_options.boundary_width_factor = 1.5;

    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = halfsize;
    const size_t n_boundary = setupOphelieThermalDirichletBoundaryMask(particles, thermal_options, geom, dp);

    OphelieJouleHeatOneWayMaterialProps material;
    material.rho = kRho;
    material.cp = kCp;
    material.k = kK;
    material.t_initial = kT0;

    const OphelieThermalDiffusionOneWayStepResult thermal =
        applyOphelieJouleHeatDiffusionOneWaySteps<MainExecutionPolicy>(
            glass_body, *glass_inner, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n,
            n_steps, thermal_options);

    const Real energy_cap_rel_err =
        (thermal.total_thermal_energy_j - thermal.total_joule_energy_j) /
        (thermal.total_joule_energy_j + TinyReal);
    const bool passed = n > 0 && n_boundary > 0 && thermal.max_temperature > kT0 + TinyReal &&
                        thermal.boundary_dirichlet_compliance > Real(0.95) && energy_cap_rel_err <= Real(1.0e-6);

    std::cout << "test_3d_ophelie_thermal_diffusion_mms n=" << n << " n_boundary=" << n_boundary << " Q=" << kQ
              << " k=" << kK << " dt=" << dt << " steps=" << n_steps << " max_T=" << thermal.max_temperature
              << " boundary_compliance=" << thermal.boundary_dirichlet_compliance
              << " E_joule_J=" << thermal.total_joule_energy_j << " E_thermal_J=" << thermal.total_thermal_energy_j
              << " energy_cap_rel_err=" << energy_cap_rel_err << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
