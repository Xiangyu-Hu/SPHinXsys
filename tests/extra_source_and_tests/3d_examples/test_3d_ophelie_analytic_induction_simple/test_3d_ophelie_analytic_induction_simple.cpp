/**
 * @file test_3d_ophelie_analytic_induction_simple.cpp
 * @brief Batch E — weak-feedback uniform B_z on a cylinder: Q ∝ r² shell gate.
 *
 * Does NOT validate self-consistent A_ind. Prescribes rotational A for uniform B₀,
 * runs complex edge-flux (+φ), and compares Q shells to the analytic ratio 25.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_skin_depth_diagnostic.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

#if SPHINXSYS_USE_SYCL
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
}
#else
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
}
#endif

struct AnalyticLocalCli
{
    Real dp = 0.05;
    Real radius = 0.20;
    Real half_height = 0.20;
    Real b0 = 0.01;
    Real sigma = 16.0;
    Real frequency_hz = 300000.0;
};

inline void applyLocalCli(int ac, char *av[], AnalyticLocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--radius=", 9) == 0)
        {
            cli.radius = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--b0=", 5) == 0)
        {
            cli.b0 = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--sigma=", 8) == 0)
        {
            cli.sigma = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--frequency=", 12) == 0)
        {
            cli.frequency_hz = static_cast<Real>(std::atof(av[i] + 12));
        }
    }
}

} // namespace

int main(int ac, char *av[])
{
    AnalyticLocalCli local;
    applyLocalCli(ac, av, local);

    std::string analytic_msg;
    if (!verifyOphelieWeakFeedbackAnalyticIdentities(analytic_msg))
    {
        std::cerr << "analytic identity failed: " << analytic_msg << std::endl;
        return 1;
    }
    std::cout << "[ophelie][analytic-induction] host_identity=" << analytic_msg << " passed=1" << std::endl;

    OphelieParameters params;
    params.sigma_glass_ = local.sigma;
    params.frequency_ = local.frequency_hz;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = true;
    params.aind_one_way_feedback_resolve_ = false;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::SolverLocal;

    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);
    logOphelieFinalParams(params, cli_options);

    const Vecd center(0.0, 0.0, local.half_height);
    const Real pad = Real(2) * local.dp;
    const BoundingBoxd bounds(Vecd(-local.radius - pad, -local.radius - pad, -pad),
                              Vecd(local.radius + pad, local.radius + pad, Real(2) * local.half_height + pad));
    SPHSystem sph_system(bounds, local.dp);
    SolidBody glass_body(sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                                         "GlassBody", center, local.radius, local.half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames names;
    RegisterOphelieGlassFields register_glass(glass_body, names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, names, params.sigma_glass_);
    assign_sigma.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = ophelieRotationalAUniformBz(pos[i], center, local.b0);
        a_coil_imag[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_coil_imag);

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine(glass_body, names);
    // Zero induced fields (weak-feedback: do not solve A_ind).
    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    syncVariableToDevice<Vecd>(particles, names.a_ind_real);
    syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
    combine.exec();

    OphelieComplexEdgeFluxSolveReport solve_report;
    const Real p_recon = execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
        glass_body, *glass_inner, names, params, nullptr, local.dp, &solve_report);

    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
    pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    const Real omega = params.omega();
    const OphelieWeakFeedbackQr2ShellMetrics shells = hostOphelieWeakFeedbackQr2ShellMetrics(
        pos, q, vol, n, center, local.radius, omega, local.b0, params.sigma_glass_);
    logOphelieWeakFeedbackQr2ShellMetrics(shells, "sph_edge_flux");

    const FrenchCylindricalFrame frame = makeFrenchCylinderFrameFromGlass(center, local.radius, local.half_height,
                                                                         FrenchVerticalAxis::Z);
    const FrenchSkinDepthRadialReport skin = computeFrenchSkinDepthRadialReport(
        pos, q, vol, n, frame, params.frequency_, params.sigma_glass_, params.mu0_, local.dp);
    logFrenchSkinDepthRadialReport(skin, "skin-depth-context");

    const Real ratio_rel =
        std::abs(shells.q_outer_over_center - Real(25)) / Real(25);
    const bool ratio_ok = ratio_rel < Real(0.35); // coarse lattice; analytic asymptote is 25
    const bool corr_ok = shells.corr_q_r2 > Real(0.85);
    const bool power_ok = shells.rel_power_err < Real(0.25);
    const bool finite_ok = n > 0 && std::isfinite(p_recon) && std::isfinite(shells.q_outer_over_center);
    const bool passed = finite_ok && ratio_ok && corr_ok && power_ok;

    std::cout << "test_3d_ophelie_analytic_induction_simple n=" << n << " dp=" << local.dp << " R=" << local.radius
              << " B0=" << local.b0 << " sigma=" << params.sigma_glass_ << " f_hz=" << params.frequency_
              << " mode=weak_feedback_prescribed_A"
              << " phi_eq_res_vol=" << solve_report.phi_eq_res_vol_imag << " P_recon_W=" << p_recon
              << " Q_outer_over_center=" << shells.q_outer_over_center << " ratio_rel_vs_25=" << ratio_rel
              << " corr_Q_r2=" << shells.corr_q_r2 << " rel_P_vs_analytic=" << shells.rel_power_err
              << " ratio_ok=" << (ratio_ok ? 1 : 0) << " corr_ok=" << (corr_ok ? 1 : 0)
              << " power_ok=" << (power_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
