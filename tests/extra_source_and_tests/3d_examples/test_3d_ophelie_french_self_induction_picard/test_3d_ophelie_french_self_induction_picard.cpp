/**
 * @file test_3d_ophelie_french_self_induction_picard.cpp
 * @brief Picard A_src = A_coil + K[J]: experimental self-induction on French reduced (not literature_passed).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

struct PicardLocalCli
{
    bool use_reload = false;
    std::string reload_dir;
    size_t max_iter = 5;
    Real relax = 0.3;
    Real j_tol = 0.05;
    Real phi_tol = 0.01;
};

inline bool reloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

inline std::string resolveDefaultFrenchReloadFolder()
{
    const StdVec<std::string> candidates = {"./reload", "../../../../../reload"};
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return "./reload";
}

inline void applyPicardLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, PicardLocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            cli.use_reload = true;
            french.dp = 0.02;
            cli.max_iter = 8;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
            cli.use_reload = true;
        }
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            french.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--self-induction-max-iter=", 26) == 0)
        {
            cli.max_iter = static_cast<size_t>(std::atoi(av[i] + 26));
        }
        else if (std::strncmp(av[i], "--self-induction-relax=", 23) == 0)
        {
            cli.relax = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strncmp(av[i], "--self-induction-tol=", 21) == 0)
        {
            cli.j_tol = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strncmp(av[i], "--self-induction-phi-tol=", 25) == 0)
        {
            cli.phi_tol = static_cast<Real>(std::atof(av[i] + 25));
        }
    }
    if (cli.use_reload && cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
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
    PicardLocalCli local_cli;
    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    applyPicardLocalCli(ac, av, french, local_cli);

    const StdVec<std::string> french_filtered = filterFrenchReducedCommandLine(ac, av, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    OphelieTestCliOptions cli_options;
    StdVec<char *> french_av;
    french_av.reserve(french_filtered.size());
    for (auto &argument : french_filtered)
    {
        french_av.push_back(const_cast<char *>(argument.c_str()));
    }
    const int french_ac = static_cast<int>(french_av.size());
    (void)filterOphelieTestCommandLine(french_ac, french_av.data(), params, cli_options);

    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
        local_cli.use_reload = true;
    }

    OphelieFrenchLiteratureProfile literature_profile;
    literature_profile.calibrate_coil_current = false;
    if (cli_options.literature_mode)
    {
        applyFrenchLiteratureMode(params, cli_options, literature_profile);
    }

    params.enable_phi_correction_ = true;
    params.enable_self_induction_ = true;
    params.self_induction_max_iterations_ = local_cli.max_iter;
    params.self_induction_relaxation_factor_ = local_cli.relax;
    params.self_induction_j_tolerance_ = local_cli.j_tol;
    params.self_induction_phi_eq_res_tolerance_ = local_cli.phi_tol;

    french.sigma_glass = params.sigma_glass_;
    french.frequency_hz = params.frequency_;
    french.target_joule_power = params.target_joule_power_;
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    logOphelieFinalParams(params, cli_options);

    if (local_cli.use_reload && !reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_self_induction_picard: Reload.xml not found under \"" << local_cli.reload_dir
                  << "\"\n";
        return 1;
    }

    const BoundingBoxd bounds =
        frenchReducedDomainBounds(french, local_cli.use_reload ? 3.0 * french.dp : 2.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setReloadParticles(local_cli.use_reload);
    if (local_cli.use_reload)
    {
        IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);
        std::cout << "[ophelie] reload folder: " << IO::getEnvironment().ReloadFolder() << std::endl;
    }

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);

    const char *particle_source = "lattice";
    if (local_cli.use_reload)
    {
        glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
        particle_source = "reload";
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
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    const OphelieFrenchSelfInductionPicardResult picard = runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params, french);

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const std::string &j_real_field = getOphelieAIndJRealFieldName(glass_names, params);
    const std::string &j_imag_field = getOphelieAIndJImagFieldName(glass_names, params);
    const Real max_j_real = hostVecdFieldMax(particles, j_real_field, n);
    const Real max_j_imag = hostVecdFieldMax(particles, j_imag_field, n);

    const bool complex_path = params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params);
    const Real phi_tol = ophelieSelfInductionPicardPhiEqResTolerance(params);
    const bool j_ok = picard.final_j_rel_change < local_cli.j_tol;
    const bool phi_ok = picard.phi_eq_res_vol < phi_tol;
    const bool converged = picard.picard_converged && j_ok && phi_ok;
    const bool passed =
        n > 0 && std::isfinite(picard.joule_power_w) && picard.a_ind_over_a_coil > TinyReal && converged;

    std::cout << "test_3d_ophelie_french_self_induction_picard particles=" << particle_source << " n=" << n
              << " dp=" << french.dp << " current_form=" << ophelieCurrentFormKindName(params.ophelie_current_form_)
              << " edge_flux_complex=" << (params.edge_flux_complex_ ? 1 : 0)
              << " complex_picard_path=" << (complex_path ? 1 : 0) << " max_iter=" << local_cli.max_iter
              << " relax=" << local_cli.relax << " j_tol=" << local_cli.j_tol << " phi_tol=" << phi_tol
              << " self_induction_iters=" << picard.self_induction_iterations
              << " final_J_rel=" << picard.final_j_rel_change << " j_ok=" << (j_ok ? 1 : 0)
              << " phi_eq_res_vol=" << picard.phi_eq_res_vol << " phi_ok=" << (phi_ok ? 1 : 0)
              << " picard_converged=" << (picard.picard_converged ? 1 : 0) << " converged=" << (converged ? 1 : 0)
              << " P_joule_W=" << picard.joule_power_w << " max_J_real=" << max_j_real << " max_J_imag=" << max_j_imag
              << " A_ind_over_A_coil=" << picard.a_ind_over_a_coil << " B_ind_over_B_coil=" << picard.b_ind_over_b_coil
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
