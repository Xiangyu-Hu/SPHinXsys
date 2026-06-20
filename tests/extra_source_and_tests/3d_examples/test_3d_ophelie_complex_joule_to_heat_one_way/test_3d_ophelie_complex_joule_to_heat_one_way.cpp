/**
 * @file test_3d_ophelie_complex_joule_to_heat_one_way.cpp
 * @brief One explicit thermal step from complex edge-flux Q_complex (no EM feedback).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "sphinxsys.h"

#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
constexpr Real kRho = 2500.0;
constexpr Real kCp = 1000.0;
constexpr Real kT0 = 300.0;
} // namespace

int main(int, char *[])
{
    const Real sigma = 16.0;
    const Real omega = Real(2) * Pi * 50.0e3;
    const Real dp = 0.04;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const Vecd a_uniform(0.002, 0.0, 0.0);
    const Real dt = 0.1;

    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));
    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieParameters params;
    params.frequency_ = omega / (Real(2) * Pi);
    params.sigma_glass_ = sigma;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    OphelieTestCliOptions cli_options;
    finalizeOphelieCurrentFormConfiguration(params, cli_options);

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, sigma);
    assign_sigma.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    hostRegisterOphelieTemperatureField(particles, kT0);

    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Real *phi_real = particles.getVariableDataByName<Real>(glass_names.phi_real);
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = a_uniform;
        a_coil_imag[i] = Vecd::Zero();
        phi_real[i] = Real(0);
        phi_imag[i] = Real(0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);

    (void)execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    syncOphelieJouleHeatPrimaryForThermalOneWay<MainExecutionPolicy>(glass_body, glass_names, params);

    const std::string &q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
    OphelieJouleHeatOneWayMaterialProps material;
    material.rho = kRho;
    material.cp = kCp;
    material.t_initial = kT0;
    const OphelieJouleHeatOneWayStepResult thermal = applyOphelieJouleHeatOneWayTemperatureSteps<MainExecutionPolicy>(
        glass_body, particles, glass_names.joule_heat, kOphelieTemperatureField, dt, material, n, 1);
    const Real max_rel_err = thermal.max_per_particle_rel_err;

    const Real p_recon = hostEdgeFluxReconPower(particles, glass_names, n, params);
    const bool passed = n > 0 && std::isfinite(p_recon) && p_recon > TinyReal && max_rel_err < Real(1.0e-6);

    std::cout << "test_3d_ophelie_complex_joule_to_heat_one_way n=" << n << " edge_flux_complex=1"
              << " q_field=" << q_field << " P_recon=" << p_recon << " dt=" << dt << " rho=" << kRho << " cp=" << kCp
              << " max_rel_err=" << max_rel_err << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
