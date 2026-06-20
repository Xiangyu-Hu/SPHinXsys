/**
 * @file test_3d_ophelie_edge_flux_scaling.cpp
 * @brief Edge-flux MMS box: I×k => P_recon×k² (imag + complex). Absolute target_P is diagnostic-only.
 */
#include "electromagnetic_ophelie.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

static Real runUniformInductionPower(const Vecd &a_base, const Real current_scale, const Real sigma, const Real omega,
                                     const Real dp, const Vecd &center, const Vecd &halfsize, const bool complex_mode)
{
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
    params.edge_flux_complex_ = complex_mode;
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
    const Vecd a_uniform = current_scale * a_base;
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    Real *phi_real = particles.getVariableDataByName<Real>(glass_names.phi_real);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = a_uniform;
        a_coil_imag[i] = Vecd::Zero();
        phi_imag[i] = -omega * a_uniform.dot(pos[i]);
        phi_real[i] = Real(0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);

    const OphelieEdgeFluxPowerMetrics power_metrics =
        execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    return power_metrics.p_total_recon;
}

static bool checkRatio(const char *label, Real ratio, Real expected, Real tol)
{
    const bool ok = std::abs(ratio - expected) <= tol;
    std::cout << "test_3d_ophelie_edge_flux_scaling " << label << " ratio=" << ratio << " expected=" << expected
              << " ok=" << (ok ? 1 : 0) << std::endl;
    return ok;
}

static bool runScalingSuite(const char *suite_label, const Vecd &a_base, const Real sigma, const Real omega,
                            const Real dp, const Vecd &center, const Vecd &halfsize, const bool complex_mode,
                            const Real tol)
{
    const Real p_ref = runUniformInductionPower(a_base, 1.0, sigma, omega, dp, center, halfsize, complex_mode);
    const Real p_half = runUniformInductionPower(a_base, 0.5, sigma, omega, dp, center, halfsize, complex_mode);
    const Real p_double = runUniformInductionPower(a_base, 2.0, sigma, omega, dp, center, halfsize, complex_mode);

    const bool i_squared_passed =
        p_ref > TinyReal && checkRatio((std::string(suite_label) + " P(0.5I)/P(1I)").c_str(), p_half / p_ref,
                                       Real(0.25), tol) &&
        checkRatio((std::string(suite_label) + " P(2I)/P(1I)").c_str(), p_double / p_ref, Real(4.0), tol);

    const Real targets[] = {12500.0, 50000.0, 200000.0};
    for (const Real target_p : targets)
    {
        const Real scale = std::sqrt(target_p / (p_ref + TinyReal));
        const Real p_scaled = runUniformInductionPower(a_base, scale, sigma, omega, dp, center, halfsize, complex_mode);
        const Real rel_err = std::abs(p_scaled - target_p) / target_p;
        std::cout << "test_3d_ophelie_edge_flux_scaling " << suite_label << " target_P_diagnostic=" << target_p
                  << " P_recon=" << p_scaled << " rel_err=" << rel_err << " (not gated on MMS box)" << std::endl;
    }

    std::cout << "test_3d_ophelie_edge_flux_scaling " << suite_label << " P_ref=" << p_ref
              << " i_squared_passed=" << (i_squared_passed ? 1 : 0) << std::endl;
    return i_squared_passed;
}

} // namespace

int main(int, char *[])
{
    const Real sigma = 16.0;
    const Real omega = Real(2) * Pi * 50.0e3;
    const Real dp = 0.04;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const Vecd a_base(0.002, 0.0, 0.0);
    const Real tol = 0.15;

    const bool imag_only_passed = runScalingSuite("imag_only", a_base, sigma, omega, dp, center, halfsize, false, tol);
    const bool complex_passed = runScalingSuite("complex", a_base, sigma, omega, dp, center, halfsize, true, tol);
    const bool passed = imag_only_passed && complex_passed;

    std::cout << "test_3d_ophelie_edge_flux_scaling edge_flux_scaling_passed=" << (imag_only_passed ? 1 : 0)
              << " complex_scaling_passed=" << (complex_passed ? 1 : 0)
              << " edge_flux_scaling_passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
