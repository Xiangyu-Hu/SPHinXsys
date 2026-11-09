/**
 * @file particle_relaxation.cpp
 * @brief GPU paper B1: body-fitted particle generation + relaxation (Stanford bunny).
 * @details Unmodified STL solid → LevelSet → Lattice → relax.
 *          Pipeline timed separately for the paper:
 *            (1) STL/shape setup (no Extrude / no geometry edit)
 *            (2) Level-set build / sign correct / write
 *            (3) Lattice particle generation
 *            (4) Relaxation iterations
 *          Domain bounds from this bunny's AABB + kernel pad (not triangle_prism).
 */
#include "gpu_paper_recording.h"
#include "sphinxsys.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;
using namespace sphinxsys_gpu_paper;

using Clock = std::chrono::steady_clock;
static double secsSince(Clock::time_point t0)
{
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

//----------------------------------------------------------------------
// Frozen Stanford bunny — native AABB from binary STL (sha256 a1b1badb...).
//----------------------------------------------------------------------
std::string stl_file_name = "StanfordBunny_fixed.stl";
const Real stl_xmin = -43.148468017578125;
const Real stl_xmax = 43.148468017578125;
const Real stl_ymin = -33.433746337890625;
const Real stl_ymax = 33.433746337890625;
const Real stl_zmin = 0.8808363676071167;
const Real stl_zmax = 83.74813079833984;
const Real stl_Lx = stl_xmax - stl_xmin; // ≈86.297 (longest → L_char)
const Real stl_Ly = stl_ymax - stl_ymin; // ≈66.867
const Real stl_Lz = stl_zmax - stl_zmin; // ≈82.867

Real L_char = 1.0;     // scaled max extent (X)
Real l_over_dp = 80.0; // L1 default (B1: L_char/dp=20 不足以成形)
Real particle_spacing = L_char / l_over_dp;
Real scaling = L_char / stl_Lx;
Vecd stl_center(0.5 * (stl_xmin + stl_xmax),
                0.5 * (stl_ymin + stl_ymax),
                0.5 * (stl_zmin + stl_zmax));
Vecd translation = -scaling * stl_center;

// Half-extents after scale+center-to-origin
Real hx = 0.5 * L_char;
Real hy = 0.5 * stl_Ly * scaling;
Real hz = 0.5 * stl_Lz * scaling;
// Computational-domain pad only (kernel / level-set margin); geometry itself is untouched.
Real domain_pad = 4.0 * particle_spacing;
Vec3d domain_lower_bound(-(hx + domain_pad), -(hy + domain_pad), -(hz + domain_pad));
Vec3d domain_upper_bound(+(hx + domain_pad), +(hy + domain_pad), +(hz + domain_pad));

int relax_iterations = 1000;

void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        l_over_dp = Real(std::atof(v));
    if (l_over_dp < Real(20.0))
        l_over_dp = Real(80.0);
    if (const char *v = std::getenv("GPU_PAPER_RELAX_ITERS"))
        relax_iterations = std::atoi(v);
    if (relax_iterations < 1)
        relax_iterations = 1000;

    particle_spacing = L_char / l_over_dp;
    scaling = L_char / stl_Lx;
    translation = -scaling * stl_center;
    hx = 0.5 * L_char;
    hy = 0.5 * stl_Ly * scaling;
    hz = 0.5 * stl_Lz * scaling;
    domain_pad = Real(4.0) * particle_spacing;
    domain_lower_bound = Vec3d(-(hx + domain_pad), -(hy + domain_pad), -(hz + domain_pad));
    domain_upper_bound = Vec3d(+(hx + domain_pad), +(hy + domain_pad), +(hz + domain_pad));
}

/** Unmodified STL solid — body-fitted volume particles. */
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stl_file_name, translation, scaling);
    }
};

int main(int ac, char *av[])
{
    applyResolutionFromEnv();

    std::cout << std::fixed << std::setprecision(6)
              << "[B1] bunny native AABB=[" << stl_xmin << "," << stl_ymin << "," << stl_zmin
              << "]..[" << stl_xmax << "," << stl_ymax << "," << stl_zmax << "]\n"
              << "[B1] scaling=" << scaling << " translation=(" << translation[0] << ","
              << translation[1] << "," << translation[2] << ")\n"
              << "[B1] scaled half-extents hx,hy,hz=(" << hx << "," << hy << "," << hz << ")\n"
              << "[B1] geometry=unmodified_STL_solid domain_pad=" << domain_pad << "\n"
              << "[B1] domain=[" << domain_lower_bound[0] << "," << domain_lower_bound[1] << ","
              << domain_lower_bound[2] << "]..[" << domain_upper_bound[0] << ","
              << domain_upper_bound[1] << "," << domain_upper_bound[2] << "]\n"
              << "[B1] L_char/dp=" << l_over_dp << " dp=" << particle_spacing
              << " relax_iters=" << relax_iterations << "\n";

    auto wall_all0 = Clock::now();

    // ---- (1) STL / shape / system setup ----
    auto t_init0 = Clock::now();
    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    RealBody input_body(sph_system, makeShared<SolidBodyFromMesh>("StanfordBunny"));
    double wall_shape_s = secsSince(t_init0);

    // ---- (2) Level-set build + sign correct + write (预处理核心) ----
    auto t_ls0 = Clock::now();
    LevelSetShape &level_set_shape = input_body.defineBodyLevelSetShape(par_ck)
                                         .correctLevelSetSign()
                                         .writeLevelSet();
    double wall_levelset_s = secsSince(t_ls0);

    // ---- (3) Lattice particle generation ----
    auto t_pg0 = Clock::now();
    input_body.generateParticles<BaseParticles, Lattice>();
    double wall_particle_gen_s = secsSince(t_pg0);
    const int n_particles = int(input_body.getBaseParticles().TotalRealParticles());
    std::cout << "[B1] N_particles after Lattice = " << n_particles
              << "  wall_shape=" << wall_shape_s
              << "s  wall_levelset=" << wall_levelset_s
              << "s  wall_particle_gen=" << wall_particle_gen_s << "s\n";

    NearShapeSurface near_body_surface(input_body);
    Inner<> input_body_inner(input_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    auto &host_methods = sph_solver.getHostMethodContainer();

    auto &input_body_cell_linked_list = main_methods.addCellLinkedListDynamics(input_body);
    auto &input_body_update_inner_relation = main_methods.addRelationDynamics(input_body_inner);
    auto &random_input_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(input_body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(input_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(input_body, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(input_body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(input_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(input_body);

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "B1_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "benchmark";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/B1_stanford_bunny/schedules/"
              "B1_benchmark_I1000_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "B1", mode);
    recorder.setMeta("case_name", "3d_gpu_paper_b1_particle_relaxation");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("L_char", std::to_string(L_char));
    recorder.setMeta("L_over_dp", std::to_string(l_over_dp));
    recorder.setMeta("particle_spacing", std::to_string(particle_spacing));
    recorder.setMeta("stl_file", stl_file_name);
    recorder.setMeta("stl_sha256", "a1b1badb60b1739c90cf34c49d81f3ece36fcffba988330b72ad05664811812d");
    recorder.setMeta("stl_scaling", std::to_string(scaling));
    recorder.setMeta("domain_pad", std::to_string(domain_pad));
    recorder.setMeta("geometry_mode", "unmodified_stl_solid_body_fitted");
    recorder.setMeta("domain_lower",
                     std::to_string(domain_lower_bound[0]) + "," + std::to_string(domain_lower_bound[1]) + "," +
                         std::to_string(domain_lower_bound[2]));
    recorder.setMeta("domain_upper",
                     std::to_string(domain_upper_bound[0]) + "," + std::to_string(domain_upper_bound[1]) + "," +
                         std::to_string(domain_upper_bound[2]));
    recorder.setMeta("scaled_half_extents",
                     std::to_string(hx) + "," + std::to_string(hy) + "," + std::to_string(hz));
    recorder.setMeta("relax_iterations", std::to_string(relax_iterations));
    recorder.setMeta("N_particles_lattice", std::to_string(n_particles));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_ck_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("api", "SPHSolver_RelaxationCK_LevelSet_BunnySolid");
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");
    recorder.setMeta("wall_shape_s", std::to_string(wall_shape_s));
    recorder.setMeta("wall_levelset_s", std::to_string(wall_levelset_s));
    recorder.setMeta("wall_particle_gen_s", std::to_string(wall_particle_gen_s));
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected", "ck_device_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected", "");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected", "");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "");
    recorder.registerMissing({"surface_distance_error", "not_collected",
                              "expensive surface-distance quality metric not wired yet",
                              "case_b1", "optional post-process vs level-set"});
    recorder.registerMissing({"neighbor_distribution", "not_collected",
                              "nearest-neighbor histogram not wired yet",
                              "case_b1", "optional end-of-run quality"});

    std::vector<double> probe_iters = loadProbeTimesFromBoundSchedule(schedule_path);
    if (probe_iters.empty())
    {
        for (int j = 0; j <= 100; ++j)
            probe_iters.push_back(double(j) * double(relax_iterations) / 100.0);
    }
    std::vector<double> field_iters = loadFieldTimesFromBoundSchedule(schedule_path);
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    auto sample_probes = [&](int ite, double requested, Real residual)
    {
        ProbeSample ps;
        ps.physical_time_s = double(ite);
        ps.requested_time_s = requested;
        ps.probe_id = "relax";
        ps.body = "StanfordBunny";
        ps.quantity = "relaxation_step";
        ps.value = double(residual);
        ps.unit = "1";
        recorder.addProbe(ps);

        ProbeSample pn = ps;
        pn.quantity = "N_active";
        pn.value = double(input_body.getBaseParticles().TotalRealParticles());
        recorder.addProbe(pn);

        recorder.addDiagnostic(double(ite), "StanfordBunny", "relaxation_step", double(residual), "1", "measured");
        recorder.addDiagnostic(double(ite), "StanfordBunny", "N_active",
                               double(input_body.getBaseParticles().TotalRealParticles()), "1", "measured");
    };

    // ---- relax prep (randomize + first configuration) ----
    // Stage prints: L4 (N~6e7) previously segfaulted here before any relax_step log.
    auto t_prep0 = Clock::now();
    std::cout << "[B1] prep: randomize ..." << std::endl << std::flush;
    random_input_body_particles.exec();
    std::cout << "[B1] prep: cell_linked_list ..." << std::endl << std::flush;
    input_body_cell_linked_list.exec();
    std::cout << "[B1] prep: inner_relation ..." << std::endl << std::flush;
    input_body_update_inner_relation.exec();
    std::cout << "[B1] prep: levelset_bounding ..." << std::endl << std::flush;
    level_set_bounding.exec();
    double wall_relax_prep_s = secsSince(t_prep0);
    std::cout << "[B1] prep done in " << wall_relax_prep_s << "s" << std::endl << std::flush;

    double wall_vtp_s = 0.0;
    double wall_diag_s = 0.0;
    if (enable_vtp)
    {
        auto tv = Clock::now();
        body_state_recorder.writeToFile(0);
        wall_vtp_s += secsSince(tv);
    }
    {
        auto td = Clock::now();
        double req0 = probe_iters.empty() ? 0.0 : probe_iters[0];
        sample_probes(0, req0, Real(0.0));
        if (!probe_iters.empty())
            probe_cursor = 1;
        if (enable_vtp && !field_iters.empty() && field_iters[0] <= 0.0)
            field_cursor = 1;
        wall_diag_s += secsSince(td);
    }

    // ---- (4) Relaxation loop (solver window) ----
    auto t_relax0 = Clock::now();
    int ite_p = 0;
    Real last_residual = Real(0.0);
    while (ite_p < relax_iterations)
    {
        input_body_cell_linked_list.exec();
        input_body_update_inner_relation.exec();

        relaxation_residual.exec();
        last_residual = relaxation_scaling.exec();
        update_particle_position.exec(last_residual);
        level_set_bounding.exec();

        ite_p += 1;
        double requested = 0.0;
        {
            auto td = Clock::now();
            if (dueForSample(double(ite_p), probe_iters, probe_cursor, &requested))
                sample_probes(ite_p, requested, last_residual);
            wall_diag_s += secsSince(td);
        }
        if (enable_vtp && dueForSample(double(ite_p), field_iters, field_cursor, &requested))
        {
            auto tv = Clock::now();
            body_state_recorder.writeToFile(ite_p);
            wall_vtp_s += secsSince(tv);
        }

        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9)
                      << "Relaxation steps N = " << ite_p
                      << "	residual_step = " << last_residual
                      << "	N_particles = " << input_body.getBaseParticles().TotalRealParticles()
                      << "\n";
        }
    }
    double wall_relax_loop_s = secsSince(t_relax0);
    std::cout << "The physics relaxation process finish !" << std::endl;
    sample_probes(ite_p, double(relax_iterations), last_residual);

    // Final VTP: always for modest N; large N only if GPU_PAPER_ENABLE_VTP=1 (L4~6e7 would be huge).
    const int vtp_auto_n_max = 2000000;
    if (enable_vtp || n_particles <= vtp_auto_n_max)
    {
        auto tv = Clock::now();
        body_state_recorder.writeToFile(ite_p);
        wall_vtp_s += secsSince(tv);
        std::cout << "[B1] wrote final particle VTP at iteration " << ite_p << "\n";
    }
    else
    {
        std::cout << "[B1] skip final VTP (N=" << n_particles << " > " << vtp_auto_n_max
                  << "; set GPU_PAPER_ENABLE_VTP=1 to force)\n";
    }

    // Copy particle / level-set outputs into run_dir for the paper archive.
    try
    {
        namespace fs = std::filesystem;
        IO::initEnvironment();
        const fs::path out_dir = IO::getEnvironment().OutputFolder();
        const fs::path dst_dir = fs::path(run_dir) / "particles";
        fs::create_directories(dst_dir);
        int n_copied = 0;
        if (fs::exists(out_dir))
        {
            for (const auto &entry : fs::directory_iterator(out_dir))
            {
                if (!entry.is_regular_file())
                    continue;
                const std::string name = entry.path().filename().string();
                if (name.find("StanfordBunny") == std::string::npos &&
                    name.find("LevelSet_StanfordBunny") == std::string::npos &&
                    name.find("ShapeSPHSystemDomain") == std::string::npos)
                    continue;
                fs::copy_file(entry.path(), dst_dir / entry.path().filename(),
                              fs::copy_options::overwrite_existing);
                ++n_copied;
            }
        }
        std::cout << "[B1] archived " << n_copied << " files -> " << dst_dir.string() << "\n";
        recorder.setMeta("particle_archive_dir", dst_dir.string());
        recorder.setMeta("particle_archive_n", std::to_string(n_copied));
    }
    catch (const std::exception &ex)
    {
        std::cout << "[B1] particle archive skipped: " << ex.what() << "\n";
        recorder.registerMissing({"particle_vtp_archive", "not_collected", ex.what(),
                                  "case_b1", "copy output/*.vtp into run_dir/particles"});
    }

    double wall_process = secsSince(wall_all0);
    // Duration V2: 预处理另报；wall_solver = 纯松弛循环
    double wall_init = wall_shape_s;
    double wall_prewindow = wall_levelset_s; // level-set 预处理
    double wall_warmup = wall_particle_gen_s + wall_relax_prep_s;
    double wall_solver = wall_relax_loop_s;
    double wall_solver_inclusive = wall_solver + wall_diag_s + wall_vtp_s;

    recorder.stages().add("shape_stl", wall_shape_s);
    recorder.stages().add("levelset_build_correct_write", wall_levelset_s);
    recorder.stages().add("particle_generation_lattice", wall_particle_gen_s);
    recorder.stages().add("relax_prep_randomize_config", wall_relax_prep_s);
    recorder.stages().add("relax_loop", wall_relax_loop_s);
    recorder.stages().add("diagnostics_probes", wall_diag_s);
    recorder.stages().add("vtp_or_io", wall_vtp_s);

    recorder.setMeta("wall_relax_prep_s", std::to_string(wall_relax_prep_s));
    recorder.setMeta("wall_relax_loop_s", std::to_string(wall_relax_loop_s));

    recorder.addStepBlock("relaxation", ite_p, 0, 0, 0, 0, 0, 0, "relax_sycl_ck", "measured");
    recorder.addBodyStats(double(ite_p), "StanfordBunny",
                          int(input_body.getBaseParticles().TotalRealParticles()),
                          int(input_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    std::cout << std::fixed << std::setprecision(6)
              << "[B1] timing summary (s):\n"
              << "  shape/STL         = " << wall_shape_s << "\n"
              << "  levelset          = " << wall_levelset_s << "\n"
              << "  particle_gen      = " << wall_particle_gen_s << "\n"
              << "  relax_prep        = " << wall_relax_prep_s << "\n"
              << "  relax_loop        = " << wall_relax_loop_s << "\n"
              << "  diagnostics       = " << wall_diag_s << "\n"
              << "  vtp               = " << wall_vtp_s << "\n"
              << "  wall_process      = " << wall_process << "\n";

    recorder.finalize("ok", 0.0, double(ite_p), wall_process,
                      wall_init, wall_warmup, wall_prewindow,
                      wall_solver, wall_solver_inclusive, wall_diag_s,
                      wall_vtp_s, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return 0;
}
