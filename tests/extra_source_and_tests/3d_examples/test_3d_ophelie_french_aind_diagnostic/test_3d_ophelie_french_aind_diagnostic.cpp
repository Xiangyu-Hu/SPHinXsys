/**
 * @file test_3d_ophelie_french_aind_diagnostic.cpp
 * @brief One-way A_ind = K[J_glass]: coil-only phi solve, induced A/B ratios (no feedback).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
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

inline void applyAIndDiagLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, bool &use_reload,
                                  std::string &reload_dir)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            use_reload = true;
            french.dp = 0.02;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            reload_dir = std::string(av[i] + 13);
            use_reload = true;
        }
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            french.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
    }
    if (use_reload && reload_dir.empty())
    {
        reload_dir = resolveDefaultFrenchReloadFolder();
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
    bool use_reload = false;
    std::string reload_dir;

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    applyAIndDiagLocalCli(ac, av, french, use_reload, reload_dir);

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
        reload_dir = cli_options.reload_dir;
        use_reload = true;
    }

    OphelieFrenchLiteratureProfile literature_profile;
    literature_profile.calibrate_coil_current = cli_options.literature_calibrate_current;
    if (cli_options.literature_mode)
    {
        applyFrenchLiteratureMode(params, cli_options, literature_profile);
    }
    if (cli_options.literature_mode && !params.enable_phi_correction_)
    {
        params.enable_phi_correction_ = true;
    }

    french.sigma_glass = params.sigma_glass_;
    french.frequency_hz = params.frequency_;
    french.target_joule_power = params.target_joule_power_;
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    logOphelieFinalParams(params, cli_options);

    if (use_reload && !reloadXmlExists(reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_aind_diagnostic: Reload.xml not found under \"" << reload_dir
                  << "\" (cwd=" << fs::current_path().string() << ")\n"
                  << "  Generate: cd build && "
                     "./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/"
                     "test_3d_ophelie_french_reduced --relax=1\n";
        return 1;
    }

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, use_reload ? 3.0 * french.dp : 2.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(use_reload);
    if (use_reload)
    {
        IO::getEnvironment().resetReloadFolder(reload_dir, true);
        std::cout << "[ophelie] reload folder: " << IO::getEnvironment().ReloadFolder() << std::endl;
    }

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);

    const char *particle_source = "lattice";
    if (use_reload)
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

    const OphelieAIndOneWayDiagnostic diag = runFrenchReducedAIndOneWayDiagnostic<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params, french);

    const bool finite = std::isfinite(diag.a_ind_over_a_coil) && std::isfinite(diag.b_ind_over_b_coil) &&
                        diag.a_coil_vol_norm > TinyReal;
    const bool passed = glass_body.getBaseParticles().TotalRealParticles() > 0 && finite;

    std::cout << "test_3d_ophelie_french_aind_diagnostic particles=" << particle_source
              << " n=" << glass_body.getBaseParticles().TotalRealParticles() << " dp=" << french.dp
              << " current_form=" << ophelieCurrentFormKindName(params.ophelie_current_form_)
              << " phi_eq_res_vol=" << diag.phi_eq_res_vol << " P_joule_W=" << diag.joule_power_w
              << " A_coil_vol_norm=" << diag.a_coil_vol_norm << " A_ind_vol_norm=" << diag.a_ind_vol_norm
              << " A_ind_over_A_coil=" << diag.a_ind_over_a_coil << " B_coil_vol_norm=" << diag.b_coil_vol_norm
              << " B_ind_vol_norm=" << diag.b_ind_vol_norm << " B_ind_over_B_coil=" << diag.b_ind_over_b_coil
              << " A_coil_real_norm=" << diag.a_coil_real_vol_norm << " A_coil_imag_norm=" << diag.a_coil_imag_vol_norm
              << " A_ind_real_norm=" << diag.a_ind_real_vol_norm << " A_ind_imag_norm=" << diag.a_ind_imag_vol_norm
              << " A_ind_real_over_A_src_real=" << diag.a_ind_real_over_a_src_real
              << " A_ind_imag_over_A_src_real=" << diag.a_ind_imag_over_a_src_real
              << " B_coil_real_norm=" << diag.b_coil_real_vol_norm << " B_coil_imag_norm=" << diag.b_coil_imag_vol_norm
              << " B_ind_real_norm=" << diag.b_ind_real_vol_norm << " B_ind_imag_norm=" << diag.b_ind_imag_vol_norm
              << " B_ind_real_over_B_src_real=" << diag.b_ind_real_over_b_src_real
              << " B_ind_imag_over_B_src_real=" << diag.b_ind_imag_over_b_src_real << " max_A_coil=" << diag.max_a_coil
              << " max_A_ind=" << diag.max_a_ind << " max_A_ind_real=" << diag.max_a_ind_real
              << " max_A_ind_imag=" << diag.max_a_ind_imag << " max_B_coil=" << diag.max_b_coil
              << " max_B_ind=" << diag.max_b_ind << " max_B_ind_real=" << diag.max_b_ind_real
              << " max_B_ind_imag=" << diag.max_b_ind_imag               << " max_J_real=" << diag.max_j_real << " max_J_imag=" << diag.max_j_imag
              << " P_complex_coil_only=" << diag.p_complex_coil_only << " P_complex_total_A=" << diag.p_complex_total_a
              << " feedback_resolve=" << (diag.feedback_resolve_done ? 1 : 0)
              << " edge_res_red_imag=" << diag.edge_res_red_imag << " edge_res_red_real=" << diag.edge_res_red_real
              << " max_J_real_feedback=" << diag.max_j_real_after_feedback
              << " max_J_imag_feedback=" << diag.max_j_imag_after_feedback << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
