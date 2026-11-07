/**
 * @file test_3d_ophelie_french_aind_diagnostic.cpp
 * @brief One-way A_ind = K[J_glass]: coil-only phi solve, induced A/B ratios (no feedback).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_skin_depth_diagnostic.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_spatial_sigma_field.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "io_environment.h"
#include "sphinxsys.h"

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
                                  std::string &reload_dir, std::string &dump_rz_csv, bool &probe_coil_point,
                                  Vecd &probe_xyz, std::string &probe_out_json, std::string &sigma_mode,
                                  std::string &sigma_field_file)
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
        else if (std::strncmp(av[i], "--dump-rz-csv=", 14) == 0)
        {
            dump_rz_csv = std::string(av[i] + 14);
        }
        else if (std::strncmp(av[i], "--sigma-mode=", 13) == 0)
        {
            sigma_mode = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--sigma-field-file=", 19) == 0)
        {
            sigma_field_file = std::string(av[i] + 19);
        }
        else if (std::strncmp(av[i], "--probe-coil-point=", 19) == 0)
        {
            // Doc 47 §8: evaluate production filament Biot at an explicit observer.
            probe_coil_point = true;
            double x = 0.0, y = 0.0, z = 0.0;
            if (std::sscanf(av[i] + 19, "%lf,%lf,%lf", &x, &y, &z) == 3)
            {
                probe_xyz = Vecd(static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z));
            }
        }
        else if (std::strncmp(av[i], "--probe-out-json=", 17) == 0)
        {
            probe_out_json = std::string(av[i] + 17);
        }
    }
    if (use_reload && reload_dir.empty())
    {
        reload_dir = resolveDefaultFrenchReloadFolder();
    }
}

/** Cylindrical (r,z) probe dump for GetDP same-current comparison (doc 23 §14.3).
 *  When coil_only_fields=true (feedback off), Br/Bz total columns reuse coil Biot so the
 *  CSV is a consistent coil-only mode dump (J/Q already come from A_coil-only recon). */
inline void dumpFrenchAIndRzCsv(const std::string &path, BaseParticles &particles, const OphelieGlassFieldNames &names,
                                const OphelieFrenchReducedCaseParams &french, bool coil_only_fields)
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

    Real z_world_min = 0.0;
    Real z_world_max = 0.0;
    if (n > 0)
    {
        z_world_min = pos[0][2];
        z_world_max = pos[0][2];
        for (size_t k = 1; k < n; ++k)
        {
            z_world_min = std::min(z_world_min, pos[k][2]);
            z_world_max = std::max(z_world_max, pos[k][2]);
        }
    }
    std::cout << "[ophelie] dump_rz_csv frame: glass_center=(" << french.glass_center[0] << ","
              << french.glass_center[1] << "," << french.glass_center[2] << ") half_h=" << french.glass_half_height
              << " z0=" << z0 << " pos_z_world=[" << z_world_min << "," << z_world_max << "]"
              << " coil_z=[" << french.coil.z_min << "," << french.coil.z_max << "] R=" << french.coil.loop_radius
              << " n_loops=" << french.coil.num_loops << std::endl;

    std::ofstream out(path);
    out << "x_m,y_m,z_world_m,r_m,z_m,vol_m3,"
           "Br_coil_re_T,Br_coil_im_T,Bz_coil_re_T,Bz_coil_im_T,"
           "Br_re_T,Br_im_T,Bz_re_T,Bz_im_T,"
           "Jtheta_re_Apm2,Jtheta_im_Apm2,Q_Wpm3\n";
    out.setf(std::ios::scientific);
    out.precision(8);
    for (size_t k = 0; k < n; ++k)
    {
        const Real x = pos[k][0];
        const Real y = pos[k][1];
        const Real z_world = pos[k][2];
        const Real z = z_world - z0; // GetDP glass uses z in [0,H]
        const Real r = std::sqrt(x * x + y * y);
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const auto cyl_rz = [&](const Vecd &v) -> std::pair<Real, Real> {
            const Real br = inv_r * (v[0] * x + v[1] * y);
            return {br, v[2]};
        };
        const auto cyl_th = [&](const Vecd &v) -> Real { return inv_r * (-v[0] * y + v[1] * x); };
        const auto bc_re = cyl_rz(bc_r[k]);
        const auto bc_im = cyl_rz(bc_i[k]);
        const auto bs_re = coil_only_fields ? bc_re : cyl_rz(bs_r[k]);
        const auto bs_im = coil_only_fields ? bc_im : cyl_rz(bs_i[k]);
        out << x << "," << y << "," << z_world << "," << r << "," << z << "," << vol[k] << "," << bc_re.first << ","
            << bc_im.first << "," << bc_re.second << "," << bc_im.second << "," << bs_re.first << "," << bs_im.first
            << "," << bs_re.second << "," << bs_im.second << "," << cyl_th(jr[k]) << "," << cyl_th(ji[k]) << "," << q[k]
            << "\n";
    }
    std::cout << "[ophelie] wrote rz csv n=" << n << " coil_only_fields=" << (coil_only_fields ? 1 : 0)
              << " path=" << path << std::endl;
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
    std::string dump_rz_csv;
    bool probe_coil_point = false;
    Vecd probe_xyz = Vecd(Real(1.0e-4), Real(0), Real(0.105));
    std::string probe_out_json;
    std::string sigma_mode = "constant";
    std::string sigma_field_file;

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    applyAIndDiagLocalCli(ac, av, french, use_reload, reload_dir, dump_rz_csv, probe_coil_point, probe_xyz,
                          probe_out_json, sigma_mode, sigma_field_file);

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

    // Doc 47 §8: production filament Biot at an explicit observer (no mesh / no phi).
    if (probe_coil_point)
    {
        StdVec<OphelieCurrentMomentSample> moments;
        buildMultiloopFilamentMoments(french.coil, moments);
        Vecd a_r = Vecd::Zero();
        Vecd a_i = Vecd::Zero();
        Vecd b_r = Vecd::Zero();
        Vecd b_i = Vecd::Zero();
        accumulateBiotSavartFromMoments(moments, probe_xyz, params.mu0_, params.softening_length_, a_r, a_i, b_r,
                                        b_i);
        const Real x = probe_xyz[0];
        const Real y = probe_xyz[1];
        const Real z = probe_xyz[2];
        const Real r = std::sqrt(x * x + y * y);
        const Real inv_r = (r > TinyReal) ? (Real(1) / r) : Real(0);
        const Real Br_re = inv_r * (b_r[0] * x + b_r[1] * y);
        const Real Br_im = inv_r * (b_i[0] * x + b_i[1] * y);
        const Real Bz_re = b_r[2];
        const Real Bz_im = b_i[2];
        std::cout.setf(std::ios::scientific);
        std::cout.precision(16);
        std::cout << "[ophelie][batchC] sph_production_coil_probe"
                  << " xyz=" << x << "," << y << "," << z << " r=" << r
                  << " Br_re=" << Br_re << " Br_im=" << Br_im << " Bz_re=" << Bz_re << " Bz_im=" << Bz_im
                  << " |Bz|=" << std::hypot(Bz_re, Bz_im)
                  << " n_moments=" << moments.size()
                  << " R=" << french.coil.loop_radius << " n_loops=" << french.coil.num_loops
                  << " I=" << french.coil.current_per_loop
                  << " z=[" << french.coil.z_min << "," << french.coil.z_max << "]"
                  << " segments=" << french.coil.segments_per_loop
                  << " cell_centered=" << (french.coil.use_cell_centered_loops ? 1 : 0)
                  << " softening=" << params.softening_length_ << " mu0=" << params.mu0_ << std::endl;
        if (!probe_out_json.empty())
        {
            std::ofstream out(probe_out_json);
            out.setf(std::ios::scientific);
            out.precision(16);
            out << "{\n"
                << "  \"source\": \"test_3d_ophelie_french_aind_diagnostic --probe-coil-point\",\n"
                << "  \"evaluator\": \"buildMultiloopFilamentMoments+accumulateBiotSavartFromMoments\",\n"
                << "  \"xyz_m\": [" << x << ", " << y << ", " << z << "],\n"
                << "  \"Br_re_T\": " << Br_re << ",\n"
                << "  \"Br_im_T\": " << Br_im << ",\n"
                << "  \"Bz_re_T\": " << Bz_re << ",\n"
                << "  \"Bz_im_T\": " << Bz_im << ",\n"
                << "  \"coil\": {\n"
                << "    \"R_m\": " << french.coil.loop_radius << ",\n"
                << "    \"n_loops\": " << french.coil.num_loops << ",\n"
                << "    \"I_per_loop_A\": " << french.coil.current_per_loop << ",\n"
                << "    \"z_min_m\": " << french.coil.z_min << ",\n"
                << "    \"z_max_m\": " << french.coil.z_max << ",\n"
                << "    \"segments_per_loop\": " << french.coil.segments_per_loop << ",\n"
                << "    \"cell_centered\": " << (french.coil.use_cell_centered_loops ? "true" : "false") << ",\n"
                << "    \"n_moments\": " << moments.size() << ",\n"
                << "    \"softening_length_m\": " << params.softening_length_ << ",\n"
                << "    \"mu0\": " << params.mu0_ << "\n"
                << "  }\n"
                << "}\n";
            std::cout << "[ophelie][batchC] wrote " << probe_out_json << std::endl;
        }
        return 0;
    }

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
    BaseParticles &particles = glass_body.getBaseParticles();
    OphelieSpatialSigmaLoadReport sigma_load{};
    if (sigma_mode == "spatial-field")
    {
        if (sigma_field_file.empty())
        {
            std::cerr << "test_3d_ophelie_french_aind_diagnostic: --sigma-mode=spatial-field requires "
                         "--sigma-field-file=\n";
            return 2;
        }
        sigma_load = loadAndAssignOphelieSpatialSigmaField(particles, glass_names, sigma_field_file);
        if (!sigma_load.ok)
        {
            std::cerr << "test_3d_ophelie_french_aind_diagnostic: spatial sigma load failed: " << sigma_load.error
                      << std::endl;
            return 2;
        }
        std::cout << "[ophelie] sigma_mode=spatial-field match=" << sigma_load.match_mode
                  << " sigma=[" << sigma_load.sigma_min << "," << sigma_load.sigma_max
                  << "] mean_vol=" << sigma_load.sigma_mean_vol << std::endl;
    }
    else
    {
        if (sigma_mode != "constant")
        {
            std::cerr << "test_3d_ophelie_french_aind_diagnostic: --sigma-mode must be constant|spatial-field\n";
            return 2;
        }
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names,
                                                                                    params.sigma_glass_);
        assign_sigma.exec();
        sigma_load.ok = true;
        sigma_load.match_mode = "constant";
        sigma_load.sigma_min = sigma_load.sigma_max = sigma_load.sigma_mean_vol = params.sigma_glass_;
    }

    const OphelieAIndOneWayDiagnostic diag = runFrenchReducedAIndOneWayDiagnostic<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, params, french);

    if (!dump_rz_csv.empty())
    {
        // feedback_resolve_done=false ⇒ J/Q are coil-only; keep B columns consistent (Doc 40 R05).
        dumpFrenchAIndRzCsv(dump_rz_csv, glass_body.getBaseParticles(), glass_names, french,
                            /*coil_only_fields=*/!diag.feedback_resolve_done);
    }

    {
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
        syncVariableToHost<Real>(particles, q_field);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *q = particles.getVariableDataByName<Real>(q_field);
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const size_t n = particles.TotalRealParticles();
        const FrenchCylindricalFrame skin_frame =
            makeFrenchCylinderFrameFromGlass(french.glass_center, french.glass_radius, french.glass_half_height,
                                             FrenchVerticalAxis::Z);
        // Doc156: skin-depth diagnostic uses volume-mean of particle Sigma when spatial-field.
        const Real sigma_for_skin =
            (sigma_mode == "spatial-field") ? sigma_load.sigma_mean_vol : params.sigma_glass_;
        const FrenchSkinDepthRadialReport skin = computeFrenchSkinDepthRadialReport(
            pos, q, vol, n, skin_frame, french.frequency_hz, sigma_for_skin, params.mu0_, french.dp);
        logFrenchSkinDepthRadialReport(skin, "skin-depth");
    }

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
