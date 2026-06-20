/**
 * @file test_3d_ophelie_phi_biot_rhs_solvability.cpp
 * @brief French-reduced Biot A: div(σA) RHS vs legacy-flux under DivSigmaGrad LHS (lattice or reload).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstring>
#include <filesystem>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

struct SolvabilityCli
{
    Real dp = 0.06;
    UnsignedInt gmres_outer = 80;
    UnsignedInt gmres_restart = 60;
    bool use_reload = false;
    bool gmres_device_ops_parity = false;
    std::string reload_dir;
    std::string csv_path = "./output/ophelie_phi_rhs_solvability.csv";
};

inline SolvabilityCli parseSolvabilityCli(int ac, char *av[])
{
    SolvabilityCli cli;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            cli.use_reload = true;
            cli.dp = 0.02;
            cli.gmres_outer = 120;
            cli.gmres_restart = 80;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--phi-gmres-max-outer-iter=", 25) == 0)
        {
            cli.gmres_outer = static_cast<UnsignedInt>(std::atoi(av[i] + 25));
        }
        else if (std::strncmp(av[i], "--phi-gmres-restart=", 18) == 0)
        {
            cli.gmres_restart = static_cast<UnsignedInt>(std::atoi(av[i] + 18));
        }
        else if (std::strncmp(av[i], "--csv=", 6) == 0)
        {
            cli.csv_path = std::string(av[i] + 6);
        }
        else if (std::strcmp(av[i], "--gmres-device-ops-parity") == 0)
        {
            cli.gmres_device_ops_parity = true;
        }
    }
    return cli;
}

inline bool ophelieReloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

/** Works when cwd is build/ or test .../bin/ (French relax writes build/reload). */
inline std::string resolveDefaultFrenchReloadFolder()
{
    const StdVec<std::string> candidates = {"./reload", "../../../../../reload"};
    for (const std::string &candidate : candidates)
    {
        if (ophelieReloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return "./reload";
}

inline void writeSolvabilityCsv(const std::string &path, const char *particle_source, size_t n, Real dp,
                                UnsignedInt gmres_outer, const OpheliePhiRhsAlignmentMetrics &alignment,
                                const OpheliePhiRhsSolveDiagnostic &div_diag,
                                const OpheliePhiRhsSolveDiagnostic &legacy_diag)
{
    std::ofstream csv(path, std::ios::app);
    if (!csv)
    {
        return;
    }
    csv.seekp(0, std::ios::end);
    const bool write_header = csv.tellp() == 0;
    if (write_header)
    {
        csv << "particle_source,n,dp,gmres_outer,rhs_div_vs_legacy_vol,rhs_cosine_div_legacy,"
               "rhs_div_vs_neg_legacy_vol,rhs_cosine_div_neg_legacy,"
               "div_eq_res,legacy_eq_res,div_cross_eq_res_legacy,legacy_cross_eq_res_div,"
               "div_divJ_red,legacy_divJ_red\n";
    }
    csv << particle_source << "," << n << "," << dp << "," << gmres_outer << "," << alignment.rhs_div_vs_legacy_vol
        << "," << alignment.rhs_cosine_div_legacy << "," << alignment.rhs_div_vs_neg_legacy_vol << ","
        << alignment.rhs_cosine_div_neg_legacy << "," << div_diag.phi_eq_res_vol << "," << legacy_diag.phi_eq_res_vol
        << "," << div_diag.cross_eq_res_other_rhs << "," << legacy_diag.cross_eq_res_other_rhs << ","
        << div_diag.div_j_l2_reduction << "," << legacy_diag.div_j_l2_reduction << "\n";
}

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

} // namespace

int main(int ac, char *av[])
{
    SolvabilityCli cli = parseSolvabilityCli(ac, av);
    if (cli.use_reload)
    {
        if (cli.reload_dir.empty())
        {
            cli.reload_dir = resolveDefaultFrenchReloadFolder();
        }
        if (!ophelieReloadXmlExists(cli.reload_dir))
        {
            std::cerr << "test_3d_ophelie_phi_biot_rhs_solvability: Reload.xml not found under \""
                      << cli.reload_dir << "\" (cwd=" << fs::current_path().string() << ")\n"
                      << "  Generate: cd build && "
                         "./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/"
                         "test_3d_ophelie_french_reduced --relax=1\n"
                      << "  Or pass:   --reload-dir=/path/to/build/reload\n";
            return 1;
        }
    }

    OphelieFrenchReducedCaseParams french;
    OphelieParameters params;
    french.dp = cli.dp;
    applyFrenchReducedDefaults(params, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
    params.phi_gmres_max_outer_iterations_ = cli.gmres_outer;
    params.phi_gmres_restart_dimension_ = cli.gmres_restart;
    params.phi_gauge_penalty_ = 1.0;
    params.enable_phi_correction_ = true;
    params.enable_power_scaling_ = false;

    const Real boundary_width = cli.use_reload ? 3.0 * french.dp : 2.0 * french.dp;
    const BoundingBoxd bounds = frenchReducedDomainBounds(french, boundary_width);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(cli.use_reload);
    if (cli.use_reload)
    {
        IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);
        std::cout << "[ophelie] reload folder: " << IO::getEnvironment().ReloadFolder() << std::endl;
    }

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);

    const char *particle_source = "lattice";
    if (cli.use_reload)
    {
        glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
        particle_source = "reload";
        std::cout << "[ophelie] loaded GlassBody from reload" << std::endl;
    }
    else
    {
        glass_body.generateParticles<BaseParticles, Lattice>();
    }

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    applyMultiloopFilamentBiotToGlass(particles, glass_names, french.coil, params.mu0_, params.softening_length_);
    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a(glass_body, glass_names);
    combine_a.exec();

    StdVec<Real> rhs_div(n, Real(0));
    StdVec<Real> rhs_legacy(n, Real(0));
    OphelieParameters div_params = params;
    div_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_params);
    hostReadScalarField(particles, glass_names.phi_rhs_imag, rhs_div.data(), n);
    OphelieParameters legacy_params = params;
    legacy_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, legacy_params);
    hostReadScalarField(particles, glass_names.phi_rhs_imag, rhs_legacy.data(), n);
    const OpheliePhiRhsAlignmentMetrics alignment =
        measurePhiRhsAlignmentMetrics(particles, rhs_div.data(), rhs_legacy.data(), n);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> level0(glass_body, glass_names, params);
    level0.exec();
    const OphelieDivJMetrics div_j_level0 =
        computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, french.glass_radius);
    const Real div_j_l2_level0 = div_j_level0.div_j_weighted_l2;

    OphelieParameters params_host = params;
    params_host.phi_gmres_use_device_vector_ops_ = false;
    OpheliePhiRhsSolveDiagnostic div_diag = runPhiRhsSolveDiagnostic<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params_host, OpheliePhiRhsOperatorKind::DivSigmaA, div_j_l2_level0,
        french.glass_radius);

    Real gmres_device_ops_parity_err = 0.0;
    bool gmres_device_ops_parity_ok = true;
    if (cli.gmres_device_ops_parity)
    {
        const Real eq_res_host = div_diag.phi_eq_res_vol;
        OphelieParameters params_device = params;
        params_device.phi_gmres_use_device_vector_ops_ = true;
        params_device.phi_gmres_use_device_krylov_storage_ = false;
        const OpheliePhiRhsSolveDiagnostic div_device = runPhiRhsSolveDiagnostic<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params_device, OpheliePhiRhsOperatorKind::DivSigmaA, div_j_l2_level0,
            french.glass_radius);
        div_diag = div_device;
        gmres_device_ops_parity_err =
            std::abs(eq_res_host - div_device.phi_eq_res_vol) / (std::abs(eq_res_host) + Real(1));
        gmres_device_ops_parity_ok = gmres_device_ops_parity_err < Real(1e-3);
        std::cout << "[ophelie] gmres_device_ops_parity eq_res_host=" << eq_res_host
                  << " eq_res_device=" << div_device.phi_eq_res_vol << " rel_err=" << gmres_device_ops_parity_err
                  << " ok=" << (gmres_device_ops_parity_ok ? 1 : 0) << std::endl;
    }

    const OpheliePhiRhsSolveDiagnostic legacy_diag = runPhiRhsSolveDiagnostic<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params, OpheliePhiRhsOperatorKind::LegacyFlux, div_j_l2_level0,
        french.glass_radius);

    const bool div_better_eq = div_diag.phi_eq_res_vol <= legacy_diag.phi_eq_res_vol;
    const bool div_better_divj = div_diag.div_j_l2_reduction >= legacy_diag.div_j_l2_reduction;
    const bool sign_flip_ok = alignment.rhs_cosine_div_neg_legacy > 0.95 &&
                              alignment.rhs_div_vs_neg_legacy_vol < 0.35;
    const bool passed =
        n > 0 && alignment.rhs_div_vs_legacy_vol > TinyReal && div_better_divj &&
        (!cli.gmres_device_ops_parity || gmres_device_ops_parity_ok);

    writeSolvabilityCsv(cli.csv_path, particle_source, n, french.dp, cli.gmres_outer, alignment, div_diag, legacy_diag);

    std::cout << "test_3d_ophelie_phi_biot_rhs_solvability particles=" << particle_source << " n=" << n
              << " dp=" << french.dp << " gmres_outer=" << cli.gmres_outer
              << " rhs_div_vs_legacy_vol=" << alignment.rhs_div_vs_legacy_vol
              << " rhs_cosine_div_legacy=" << alignment.rhs_cosine_div_legacy
              << " rhs_div_vs_neg_legacy_vol=" << alignment.rhs_div_vs_neg_legacy_vol
              << " rhs_cosine_div_neg_legacy=" << alignment.rhs_cosine_div_neg_legacy
              << " sign_flip_ok=" << (sign_flip_ok ? 1 : 0)
              << " div_sigma_a_eq_res=" << div_diag.phi_eq_res_vol
              << " legacy_flux_eq_res=" << legacy_diag.phi_eq_res_vol
              << " div_cross_eq_res_legacy=" << div_diag.cross_eq_res_other_rhs
              << " legacy_cross_eq_res_div=" << legacy_diag.cross_eq_res_other_rhs
              << " div_sigma_a_divJ_red=" << div_diag.div_j_l2_reduction
              << " legacy_flux_divJ_red=" << legacy_diag.div_j_l2_reduction
              << " divJ_better_with_div_rhs=" << (div_better_divj ? 1 : 0)
              << " eq_res_better_with_div_rhs=" << (div_better_eq ? 1 : 0)
              << " gmres_device_ops_parity=" << (cli.gmres_device_ops_parity ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
