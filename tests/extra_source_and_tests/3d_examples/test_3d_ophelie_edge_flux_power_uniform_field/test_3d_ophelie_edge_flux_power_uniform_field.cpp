/**
 * @file test_3d_ophelie_edge_flux_power_uniform_field.cpp
 * @brief Uniform-field edge-flux power closure: P_recon vs P_exact; P_graph diagnostic only.
 */
#include "electromagnetic_ophelie.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

static Real hostTotalVolume(BaseParticles &particles, size_t n)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real total = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        total += vol[i];
    }
    return total;
}

static bool runUniformFieldPowerCase(const std::string &case_name, const Vecd &e_uniform, const Vecd &a_uniform,
                                     const Real omega, const Real sigma, const Real dp, const Vecd &center,
                                     const Vecd &halfsize, bool set_phi_from_e, bool edge_flux_complex,
                                     Real &p_recon_over_exact_out, Real &p_graph_over_exact_out)
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
    params.edge_flux_complex_ = edge_flux_complex;
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
        phi_imag[i] = (!edge_flux_complex && set_phi_from_e) ? -e_uniform.dot(pos[i]) : Real(0);
        phi_real[i] = (edge_flux_complex && set_phi_from_e) ? -e_uniform.dot(pos[i]) : Real(0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);

    const OphelieEdgeFluxPowerMetrics power_metrics =
        execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    const Real total_volume = hostTotalVolume(particles, n);
    const Real p_exact = Real(0.5) * sigma * e_uniform.squaredNorm() * total_volume;
    p_recon_over_exact_out = power_metrics.p_total_recon / (p_exact + TinyReal);
    p_graph_over_exact_out = power_metrics.p_graph_edge / (p_exact + TinyReal);

    const bool passed = n > 0 && p_recon_over_exact_out > Real(0.5) && p_recon_over_exact_out < Real(2.0);
    std::cout << "test_3d_ophelie_edge_flux_power_uniform_field case=" << case_name
              << " edge_flux_complex=" << (edge_flux_complex ? 1 : 0) << " n=" << n << " P_exact=" << p_exact
              << " P_recon=" << power_metrics.p_total_recon
              << " P_graph_edge=" << power_metrics.p_graph_edge << " P_recon/P_exact=" << p_recon_over_exact_out
              << " P_graph/P_exact=" << p_graph_over_exact_out << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed;
}

static bool runComplexInductionPowerCase(const std::string &case_name, const Vecd &a_uniform, bool a_on_imag_chain,
                                         const Real omega, const Real sigma, const Real dp, const Vecd &center,
                                         const Vecd &halfsize, Real &p_recon_over_exact_out)
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
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Real *phi_real = particles.getVariableDataByName<Real>(glass_names.phi_real);
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = a_on_imag_chain ? a_uniform : Vecd::Zero();
        a_coil_imag[i] = a_on_imag_chain ? Vecd::Zero() : a_uniform;
        phi_real[i] = Real(0);
        phi_imag[i] = Real(0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);

    const OphelieEdgeFluxPowerMetrics power_metrics =
        execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    const Real total_volume = hostTotalVolume(particles, n);
    const Real p_exact = Real(0.5) * sigma * omega * omega * a_uniform.squaredNorm() * total_volume;
    p_recon_over_exact_out = power_metrics.p_total_recon / (p_exact + TinyReal);

    const bool passed = n > 0 && p_recon_over_exact_out > Real(0.5) && p_recon_over_exact_out < Real(2.0);
    std::cout << "test_3d_ophelie_edge_flux_power_uniform_field case=" << case_name
              << " edge_flux_complex=1 n=" << n << " P_exact=" << p_exact
              << " P_recon=" << power_metrics.p_total_recon << " P_recon/P_exact=" << p_recon_over_exact_out
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed;
}

int main(int, char *[])
{
    const Real sigma = 16.0;
    const Real omega = Real(2) * Pi * 50.0e3;
    const Real dp = 0.04;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);

    Real p_recon_over_exact = 0.0;
    Real p_graph_over_exact = 0.0;
    const Vecd e_potential(120.0, 0.0, 0.0);
    const bool case1_passed =
        runUniformFieldPowerCase("potential_field", e_potential, Vecd::Zero(), omega, sigma, dp, center, halfsize, true,
                                 false, p_recon_over_exact, p_graph_over_exact);

    const Vecd a_induction(0.002, 0.0, 0.0);
    const Vecd e_induction = -omega * a_induction;
    const bool case2_passed =
        runUniformFieldPowerCase("induction_field", e_induction, a_induction, omega, sigma, dp, center, halfsize, false,
                                 false, p_recon_over_exact, p_graph_over_exact);

    const bool complex_case_a_passed =
        runUniformFieldPowerCase("complex_potential_real", e_potential, Vecd::Zero(), omega, sigma, dp, center, halfsize,
                                 true, true, p_recon_over_exact, p_graph_over_exact);

    const bool complex_case_b_passed =
        runComplexInductionPowerCase("complex_induction_imag_chain", a_induction, true, omega, sigma, dp, center,
                                     halfsize, p_recon_over_exact);

    const bool complex_case_c_passed =
        runComplexInductionPowerCase("complex_induction_real_chain", a_induction, false, omega, sigma, dp, center,
                                     halfsize, p_recon_over_exact);

    const bool passed = case1_passed && case2_passed && complex_case_a_passed && complex_case_b_passed &&
                        complex_case_c_passed;
    std::cout << "test_3d_ophelie_edge_flux_power_uniform_field summary passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
