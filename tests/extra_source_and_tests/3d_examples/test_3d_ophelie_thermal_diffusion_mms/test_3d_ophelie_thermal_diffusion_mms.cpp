/**
 * @file test_3d_ophelie_thermal_diffusion_mms.cpp
 * @brief One SYCL SPHSystem: host Stage-1 identities, insulated Q+diffusion, then Dirichlet k=4.
 *
 * Do not construct a second SPHSystem in this process. Tearing down a SYCL SPHSystem and
 * creating another often hangs the GPU runtime.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_material_laws.h"
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
constexpr Real kK = 4.0;
constexpr Real kT0 = 300.0;
constexpr Real kQ = 5.0e4;
} // namespace

int main(int, char *[])
{
    std::string stage1_msg;
    const bool stage1_ok = verifyFrenchThermalStage1HostAlgebra(stage1_msg);
    std::cout << "[ophelie][thermal-mms] thermal_stage1_host=" << stage1_msg << " passed=" << (stage1_ok ? 1 : 0)
              << std::endl
              << std::flush;
    if (!stage1_ok)
    {
        return 1;
    }

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

    OphelieJouleHeatOneWayMaterialProps material;
    material.rho = kRho;
    material.cp = kCp;
    material.k = kK;
    material.t_initial = kT0;

    OphelieThermalDiffusionOneWayOptions insulated_options;
    insulated_options.enable_diffusion = true;
    insulated_options.enable_cold_wall_dirichlet = false;
    insulated_options.enable_french_natural_bc = false;

    const OphelieThermalDiffusionOneWayStepResult insulated =
        applyOphelieJouleHeatDiffusionOneWaySteps<MainExecutionPolicy>(
            glass_body, *glass_inner, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n,
            n_steps, insulated_options);
    const Real insulated_energy_rel =
        (insulated.total_thermal_energy_j - insulated.total_joule_energy_j) /
        (insulated.total_joule_energy_j + TinyReal);
    std::cout << "[ophelie][thermal-mms] insulated n=" << n << " max_T=" << insulated.max_temperature
              << " E_joule=" << insulated.total_joule_energy_j << " E_thermal=" << insulated.total_thermal_energy_j
              << " energy_rel=" << insulated_energy_rel << std::endl
              << std::flush;

    StateDynamics<MainExecutionPolicy, ResetOphelieThermalStateToInitialCK> reset_thermal(glass_body, kT0);
    reset_thermal.exec();

    OphelieThermalDiffusionOneWayOptions dirichlet_options;
    dirichlet_options.enable_diffusion = true;
    dirichlet_options.enable_cold_wall_dirichlet = true;
    dirichlet_options.enable_french_natural_bc = false;
    dirichlet_options.boundary_width_factor = 1.5;
    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = halfsize;
    const size_t n_boundary = setupOphelieThermalDirichletBoundaryMask(particles, dirichlet_options, geom, dp);

    const OphelieThermalDiffusionOneWayStepResult dirichlet =
        applyOphelieJouleHeatDiffusionOneWaySteps<MainExecutionPolicy>(
            glass_body, *glass_inner, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n,
            n_steps, dirichlet_options);
    const Real dirichlet_energy_rel =
        (dirichlet.total_thermal_energy_j - dirichlet.total_joule_energy_j) /
        (dirichlet.total_joule_energy_j + TinyReal);

    const bool insulated_ok = n > 0 && std::abs(insulated_energy_rel) < Real(0.02);
    const bool dirichlet_ok = n_boundary > 0 && dirichlet.max_temperature > kT0 + TinyReal &&
                              dirichlet.boundary_dirichlet_compliance > Real(0.95) &&
                              dirichlet_energy_rel < Real(-0.01);

    // Doc 45 F05: nonuniform insulated diffusion — energy conserved, T_std decreases, k scales decay.
    // Index-based T pattern avoids Position host-pointer edge cases on this small body.
    reset_thermal.exec();
    for (size_t i = 0; i < n; ++i)
    {
        joule_heat[i] = Real(0);
    }
    syncVariableToDevice<Real>(particles, glass_names.joule_heat);
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    double t_mean0 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real tz = kT0 + Real(40) * std::sin(Real(0.15) * static_cast<Real>(i));
        temperature[i] = tz;
        delta_t[i] = tz - kT0;
        t_mean0 += static_cast<double>(tz);
    }
    t_mean0 /= static_cast<double>(n);
    double var0 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double d = static_cast<double>(temperature[i]) - t_mean0;
        var0 += d * d;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);

    auto run_nonuniform = [&](Real k_val, size_t steps) {
        material.k = k_val;
        registerOphelieThermalDiffusionAuxFields(particles, k_val);
        OphelieThermalDiffusionOneWayOptions opt = insulated_options;
        return applyOphelieJouleHeatDiffusionOneWaySteps<MainExecutionPolicy>(
            glass_body, *glass_inner, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n,
            steps, opt);
    };

    const Real e0 = hostThermalEnergyConstantCp(particles, kRho, kCp, kT0, kOphelieTemperatureField);
    const auto nu_k1 = run_nonuniform(Real(1), n_steps);
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    double t_mean1 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        t_mean1 += static_cast<double>(temperature[i]);
    }
    t_mean1 /= static_cast<double>(n);
    double var1 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double d = static_cast<double>(temperature[i]) - t_mean1;
        var1 += d * d;
    }
    const Real e1 = hostThermalEnergyConstantCp(particles, kRho, kCp, kT0, kOphelieTemperatureField);
    const Real e_rel = static_cast<Real>(std::abs(static_cast<double>(e1) - static_cast<double>(e0)) /
                                         (std::abs(static_cast<double>(e0)) + 1.0));

    reset_thermal.exec();
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    for (size_t i = 0; i < n; ++i)
    {
        const Real tz = kT0 + Real(40) * std::sin(Real(0.15) * static_cast<Real>(i));
        temperature[i] = tz;
        delta_t[i] = tz - kT0;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    (void)run_nonuniform(Real(4), n_steps);
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    double t_mean4 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        t_mean4 += static_cast<double>(temperature[i]);
    }
    t_mean4 /= static_cast<double>(n);
    double var4 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double d = static_cast<double>(temperature[i]) - t_mean4;
        var4 += d * d;
    }

    const bool nonuniform_ok = var0 > 1.0 && var1 < var0 && var4 < var1 && e_rel < Real(0.05);
    std::cout << "[ophelie][thermal-mms] nonuniform_diffusion var0=" << var0 << " var_k1=" << var1
              << " var_k4=" << var4 << " energy_rel=" << e_rel << " nonuniform_ok=" << (nonuniform_ok ? 1 : 0)
              << " nu_k1_maxT=" << nu_k1.max_temperature << std::endl
              << std::flush;

    const bool passed = stage1_ok && insulated_ok && dirichlet_ok && nonuniform_ok;

    std::cout << "test_3d_ophelie_thermal_diffusion_mms n=" << n << " n_boundary=" << n_boundary << " k=" << kK
              << " insulated_energy_rel=" << insulated_energy_rel << " dir_max_T=" << dirichlet.max_temperature
              << " dir_energy_rel=" << dirichlet_energy_rel
              << " boundary_compliance=" << dirichlet.boundary_dirichlet_compliance
              << " insulated_ok=" << (insulated_ok ? 1 : 0) << " dirichlet_ok=" << (dirichlet_ok ? 1 : 0)
              << " nonuniform_ok=" << (nonuniform_ok ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
