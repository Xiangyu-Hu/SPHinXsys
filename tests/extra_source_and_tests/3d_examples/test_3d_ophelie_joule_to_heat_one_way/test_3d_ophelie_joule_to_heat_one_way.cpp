/**
 * @file test_3d_ophelie_joule_to_heat_one_way.cpp
 * @brief One explicit step: DeltaT = Q * dt / (rho * cp) from uniform JouleHeat.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

#ifndef OPHELIE_TEST_RELOAD_DIR
#define OPHELIE_TEST_RELOAD_DIR "./reload"
#endif

namespace
{
constexpr const char *kTemperatureField = "Temperature";
constexpr Real kRho = 2500.0;
constexpr Real kCp = 1000.0;
constexpr Real kT0 = 300.0;
} // namespace

int main(int, char *[])
{
    OphelieParameters params;
    const Real q_uniform = 5.0e4;
    const Real dt = 0.1;

    const Real dp = 0.2;
    const BoundingBoxd system_bounds(Vecd(-0.2, -0.2, -0.2), Vecd(0.4, 0.4, 0.8));

    SPHSystem sph_system(system_bounds, dp);
    sph_system.setReloadParticles(true);
    sph_system.setRunParticleRelaxation(false);
    IO::getEnvironment().resetReloadFolder(OPHELIE_TEST_RELOAD_DIR, true);

    SolidBody glass_body(sph_system, makeShared<ComplexShape>("GlassBody"));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    hostRegisterOphelieTemperatureField(particles, kT0);

    Real *joule_heat = particles.getVariableDataByName<Real>(glass_names.joule_heat);
    for (size_t i = 0; i < n; ++i)
    {
        joule_heat[i] = q_uniform;
    }
    syncVariableToDevice<Real>(particles, glass_names.joule_heat);

    const Real delta_t_expected = ophelieJouleHeatOneWayDeltaTExpected(q_uniform, dt, kRho, kCp);
    OphelieJouleHeatOneWayMaterialProps material;
    material.rho = kRho;
    material.cp = kCp;
    material.t_initial = kT0;
    const OphelieJouleHeatOneWayStepResult thermal = applyOphelieJouleHeatOneWayTemperatureSteps<MainExecutionPolicy>(
        glass_body, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n, 1);
    const Real max_rel_err = thermal.max_per_particle_rel_err;

    const bool passed = n > 0 && max_rel_err < 1.0e-12;
    std::cout << "test_3d_ophelie_joule_to_heat_one_way n=" << n << " Q=" << q_uniform << " dt=" << dt
              << " rho=" << kRho << " cp=" << kCp << " DeltaT_expected=" << delta_t_expected
              << " max_rel_err=" << max_rel_err << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
