/**
 * @file test_3d_ophelie_french_self_induction_picard.cpp
 * @brief Picard A_src = A_coil + K[J]: experimental self-induction on French reduced (not literature_passed).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_skin_depth_diagnostic.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_spatial_sigma_field.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

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
    std::string dump_rz_csv;
    std::string dump_em_particles_csv;
    std::string sigma_mode = "constant"; /**< constant | spatial-field */
    std::string sigma_field_file;
    std::string em_result_json;
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
        else if (std::strncmp(av[i], "--dump-rz-csv=", 14) == 0)
        {
            cli.dump_rz_csv = std::string(av[i] + 14);
        }
        else if (std::strncmp(av[i], "--dump-em-particles-csv=", 24) == 0)
        {
            cli.dump_em_particles_csv = std::string(av[i] + 24);
        }
        else if (std::strncmp(av[i], "--sigma-mode=", 13) == 0)
        {
            cli.sigma_mode = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--sigma-field-file=", 19) == 0)
        {
            cli.sigma_field_file = std::string(av[i] + 19);
        }
        else if (std::strncmp(av[i], "--em-result-json=", 17) == 0)
        {
            cli.em_result_json = std::string(av[i] + 17);
        }
    }
    if (cli.use_reload && cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
}

/** Cylindrical (r,z) probe dump for GetDP same-current comparison (doc 32). */
inline void dumpFrenchPicardRzCsv(const std::string &path, BaseParticles &particles, const OphelieGlassFieldNames &names,
                                  const OphelieFrenchReducedCaseParams &french)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.b_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_coil_imag);
    syncVariableToHost<Vecd>(particles, names.b_src_real);
    syncVariableToHost<Vecd>(particles, names.b_src_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *bc_r = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    const Vecd *bc_i = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
    const Vecd *bs_r = particles.getVariableDataByName<Vecd>(names.b_src_real);
    const Vecd *bs_i = particles.getVariableDataByName<Vecd>(names.b_src_imag);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
    const size_t n = particles.TotalRealParticles();
    const Real z0 = french.glass_center[2] - french.glass_half_height;

    std::ofstream out(path);
    out << "r_m,z_m,vol_m3,"
           "Br_coil_re_T,Br_coil_im_T,Bz_coil_re_T,Bz_coil_im_T,"
           "Br_re_T,Br_im_T,Bz_re_T,Bz_im_T,"
           "Jtheta_re_Apm2,Jtheta_im_Apm2,Q_Wpm3\n";
    out.setf(std::ios::scientific);
    out.precision(8);
    for (size_t k = 0; k < n; ++k)
    {
        const Real x = pos[k][0];
        const Real y = pos[k][1];
        const Real z = pos[k][2] - z0;
        const Real r = std::sqrt(x * x + y * y);
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const auto cyl_rz = [&](const Vecd &v) -> std::pair<Real, Real> {
            const Real br = inv_r * (v[0] * x + v[1] * y);
            return {br, v[2]};
        };
        const auto cyl_th = [&](const Vecd &v) -> Real { return inv_r * (-v[0] * y + v[1] * x); };
        const auto bc_re = cyl_rz(bc_r[k]);
        const auto bc_im = cyl_rz(bc_i[k]);
        const auto bs_re = cyl_rz(bs_r[k]);
        const auto bs_im = cyl_rz(bs_i[k]);
        out << r << "," << z << "," << vol[k] << "," << bc_re.first << "," << bc_im.first << "," << bc_re.second << ","
            << bc_im.second << "," << bs_re.first << "," << bs_im.first << "," << bs_re.second << "," << bs_im.second
            << "," << cyl_th(jr[k]) << "," << cyl_th(ji[k]) << "," << q[k] << "\n";
    }
    std::cout << "[ophelie] wrote rz csv n=" << n << " path=" << path << std::endl;
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
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    if (!local_cli.dump_em_particles_csv.empty())
    {
        if (!writeOphelieEmParticlesLatticeCsv(local_cli.dump_em_particles_csv, particles))
        {
            std::cerr << "[ophelie] failed to write --dump-em-particles-csv=" << local_cli.dump_em_particles_csv
                      << std::endl;
            return 2;
        }
        std::cout << "[ophelie] wrote em particles n=" << n << " path=" << local_cli.dump_em_particles_csv
                  << std::endl;
    }

    OphelieSpatialSigmaLoadReport sigma_load;
    sigma_load.match_mode = "constant";
    sigma_load.ok = true;
    sigma_load.n_particles = n;
    sigma_load.n_assigned = n;
    sigma_load.sigma_min = params.sigma_glass_;
    sigma_load.sigma_max = params.sigma_glass_;
    sigma_load.sigma_mean_vol = params.sigma_glass_;
    if (local_cli.sigma_mode == "constant")
    {
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names,
                                                                                    params.sigma_glass_);
        assign_sigma.exec();
    }
    else if (local_cli.sigma_mode == "spatial-field")
    {
        if (local_cli.sigma_field_file.empty())
        {
            std::cerr << "[ophelie][error] --sigma-mode=spatial-field requires --sigma-field-file=\n";
            return 2;
        }
        // Seed then overwrite; incomplete coverage hard-fails (no mean fallback).
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names,
                                                                                    params.sigma_glass_);
        assign_sigma.exec();
        sigma_load = loadAndAssignOphelieSpatialSigmaField(particles, glass_names, local_cli.sigma_field_file);
        if (!sigma_load.ok)
        {
            std::cerr << "[ophelie][error] spatial sigma load failed: " << sigma_load.error << std::endl;
            return 2;
        }
    }
    else
    {
        std::cerr << "[ophelie][error] --sigma-mode must be constant|spatial-field, got " << local_cli.sigma_mode
                  << std::endl;
        return 2;
    }

    const OphelieFrenchSelfInductionPicardResult picard = runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params, french);

    const std::string &j_real_field = getOphelieAIndJRealFieldName(glass_names, params);
    const std::string &j_imag_field = getOphelieAIndJImagFieldName(glass_names, params);
    const Real max_j_real = hostVecdFieldMax(particles, j_real_field, n);
    const Real max_j_imag = hostVecdFieldMax(particles, j_imag_field, n);

    if (!local_cli.dump_rz_csv.empty())
    {
        dumpFrenchPicardRzCsv(local_cli.dump_rz_csv, particles, glass_names, french);
    }

    {
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
        syncVariableToHost<Real>(particles, q_field);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *q = particles.getVariableDataByName<Real>(q_field);
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const FrenchCylindricalFrame skin_frame =
            makeFrenchCylinderFrameFromGlass(french.glass_center, french.glass_radius, french.glass_half_height,
                                             FrenchVerticalAxis::Z);
        const FrenchSkinDepthRadialReport skin = computeFrenchSkinDepthRadialReport(
            pos, q, vol, n, skin_frame, french.frequency_hz, sigma_load.sigma_mean_vol, params.mu0_, french.dp);
        logFrenchSkinDepthRadialReport(skin, "skin-depth");
    }

    const bool complex_path = params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params);
    const Real phi_tol = ophelieSelfInductionPicardPhiEqResTolerance(params);
    const bool j_ok = picard.fixed_point_residual < local_cli.j_tol;
    const bool phi_ok = picard.phi_eq_res_vol < phi_tol;
    const bool converged = picard.picard_converged && j_ok && phi_ok;
    const bool passed =
        n > 0 && std::isfinite(picard.joule_power_w) && picard.a_ind_over_a_coil > TinyReal && converged &&
        sigma_load.ok;

    std::cout << "test_3d_ophelie_french_self_induction_picard particles=" << particle_source << " n=" << n
              << " dp=" << french.dp << " current_form=" << ophelieCurrentFormKindName(params.ophelie_current_form_)
              << " edge_flux_complex=" << (params.edge_flux_complex_ ? 1 : 0)
              << " complex_picard_path=" << (complex_path ? 1 : 0) << " max_iter=" << local_cli.max_iter
              << " relax=" << local_cli.relax << " j_tol=" << local_cli.j_tol << " phi_tol=" << phi_tol
              << " sigma_mode=" << local_cli.sigma_mode << " sigma_match=" << sigma_load.match_mode
              << " sigma_min=" << sigma_load.sigma_min << " sigma_max=" << sigma_load.sigma_max
              << " sigma_mean_vol=" << sigma_load.sigma_mean_vol
              << " self_induction_iters=" << picard.self_induction_iterations
              << " fixed_point_res=" << picard.fixed_point_residual
              << " relaxed_step=" << picard.relaxed_step_change
              << " final_J_rel=" << picard.final_j_rel_change << " j_ok=" << (j_ok ? 1 : 0)
              << " phi_eq_res_vol=" << picard.phi_eq_res_vol << " phi_ok=" << (phi_ok ? 1 : 0)
              << " picard_converged=" << (picard.picard_converged ? 1 : 0) << " converged=" << (converged ? 1 : 0)
              << " P_joule_W=" << picard.joule_power_w << " max_J_real=" << max_j_real << " max_J_imag=" << max_j_imag
              << " A_ind_over_A_coil=" << picard.a_ind_over_a_coil << " B_ind_over_B_coil=" << picard.b_ind_over_b_coil
              << " passed=" << (passed ? 1 : 0) << std::endl;

    if (!local_cli.em_result_json.empty())
    {
        std::ofstream js(local_cli.em_result_json);
        js.setf(std::ios::scientific);
        js.precision(10);
        js << "{\n"
           << "  \"ok\": " << (passed ? "true" : "false") << ",\n"
           << "  \"converged\": " << (converged ? "true" : "false") << ",\n"
           << "  \"picard_converged\": " << (picard.picard_converged ? "true" : "false") << ",\n"
           << "  \"j_ok\": " << (j_ok ? "true" : "false") << ",\n"
           << "  \"phi_ok\": " << (phi_ok ? "true" : "false") << ",\n"
           << "  \"n_particles\": " << n << ",\n"
           << "  \"dp\": " << french.dp << ",\n"
           << "  \"sigma_mode\": \"" << local_cli.sigma_mode << "\",\n"
           << "  \"sigma_field_file\": \"" << local_cli.sigma_field_file << "\",\n"
           << "  \"sigma_match_mode\": \"" << sigma_load.match_mode << "\",\n"
           << "  \"sigma_min\": " << sigma_load.sigma_min << ",\n"
           << "  \"sigma_max\": " << sigma_load.sigma_max << ",\n"
           << "  \"sigma_mean_vol\": " << sigma_load.sigma_mean_vol << ",\n"
           << "  \"self_induction_iterations\": " << picard.self_induction_iterations << ",\n"
           << "  \"fixed_point_res\": " << picard.fixed_point_residual << ",\n"
           << "  \"phi_eq_res_vol\": " << picard.phi_eq_res_vol << ",\n"
           << "  \"P_joule_W\": " << picard.joule_power_w << ",\n"
           << "  \"dump_rz_csv\": \"" << local_cli.dump_rz_csv << "\",\n"
           << "  \"coil_excitation\": \"finite_section_prescribed_J\"\n"
           << "}\n";
        std::cout << "[ophelie] wrote em_result_json path=" << local_cli.em_result_json << std::endl;
    }
    return passed ? 0 : 1;
}
