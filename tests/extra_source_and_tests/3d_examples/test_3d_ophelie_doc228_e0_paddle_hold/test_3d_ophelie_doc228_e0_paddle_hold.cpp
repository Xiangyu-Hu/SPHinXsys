/**
 * @file test_3d_ophelie_doc228_e0_paddle_hold.cpp
 * @brief Doc228 E0 — client stirring paddle solid hold (ω=0, no fluid, no gravity).
 *
 * Requires solver-frame paddle STL from client stirring relax data/
 *   (client_stir_paddle_solver_m.stl; CMake copies into this case bin/data).
 * Regenerate: test_3d_ophelie_client_stirring_glass_relax/preprocess_client_paddle_stl.py
 *
 * Clamp: provisional top-of-shaft band (highest 50 mm in Z). Documented approximate
 * until drive-section CAD confirmation (Doc228 §2.3).
 */
#include "electromagnetic_ophelie_client_stainless_steel.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{
struct E0Cli
{
    std::string stl_path = "./data/client_stir_paddle_solver_m.stl";
    Real dp = Real(0.004); // coarse; Doc228 may mark under-resolved blade layers
    int solid_steps = 100;
    Real clamp_top_band_m = Real(0.05);
    std::string out_dir = ".";
    std::string summary_json = "doc228_e0_summary.json";
};

class ClampTopBand : public BodyPartByParticle
{
  public:
    ClampTopBand(SPHBody &body, Real z_min) : BodyPartByParticle(body), z_min_(z_min)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ClampTopBand::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    }

  private:
    Real z_min_;
    bool tagManually(size_t index_i) { return pos_[index_i][2] >= z_min_; }
};

void applyCli(int ac, char *av[], E0Cli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--stl=", 6) == 0)
            cli.stl_path = av[i] + 6;
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        else if (std::strncmp(av[i], "--solid-steps=", 14) == 0)
            cli.solid_steps = std::atoi(av[i] + 14);
        else if (std::strncmp(av[i], "--clamp-top-band=", 17) == 0)
            cli.clamp_top_band_m = static_cast<Real>(std::atof(av[i] + 17));
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            cli.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            cli.summary_json = av[i] + 15;
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc228 E0 paddle hold: --stl= --dp= --solid-steps= --clamp-top-band= "
                         "--out-dir= --summary-json=\n";
            std::exit(0);
        }
    }
}

class ClientPaddleShape : public ComplexShape
{
  public:
    ClientPaddleShape(const std::string &name, const std::string &stl_path)
        : ComplexShape(name)
    {
        // STL already in solver metres (+Z up); no extra scale/translation.
        add<TriangleMeshShapeSTL>(stl_path, Vecd::Zero(), Real(1), "ClientPaddle");
    }
};
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_doc228_e0_paddle_hold requires SYCL build.\n";
    return 1;
#else
    E0Cli cli;
    applyCli(ac, av, cli);
    fs::create_directories(cli.out_dir);

    if (!fs::exists(cli.stl_path))
    {
        std::cerr << "[ophelie][doc228][E0] FATAL: missing STL \"" << cli.stl_path
                  << "\". From stirring relax: python3 preprocess_client_paddle_stl.py "
                     "(then rebuild so CMake copies data/).\n";
        return 2;
    }

    const OphelieClientStainless304Development steel = makeClientStainless304Development();
    writeClientStainless304Json(steel, (fs::path(cli.out_dir) / "stainless_304_dev_preset.json").string());

    auto mesh = makeShared<ClientPaddleShape>("PaddleBody", cli.stl_path);
    const Real pad = Real(4) * cli.dp;
    BoundingBoxd bb = mesh->getBounds().expand(pad);

    SPHSystem system(bb, cli.dp);
    SolidBody paddle(system, mesh);
    paddle.defineAdaptationRatios(1.3, 1.0);
    paddle.defineBodyLevelSetShape(par_ck, 2).correctLevelSetSign().cleanLevelSet();
    // Doc228: real stainless E — do not secretly soften.
    paddle.defineMatterMaterial<NeoHookeanSolid>(steel.rho_kg_m3, steel.youngs_modulus_pa, steel.poisson);
    paddle.generateParticles<BaseParticles, Lattice>();

    const size_t n = paddle.getBaseParticles().TotalRealParticles();
    if (n == 0)
    {
        std::cerr << "[ophelie][doc228][E0] FATAL: zero paddle particles (check STL / dp)\n";
        return 2;
    }

    // Clamp top band of current AABB (provisional drive zone).
    const Real z_max = mesh->getBounds().upper_[2];
    const Real z_clamp = z_max - cli.clamp_top_band_m;
    ClampTopBand clamp_part(paddle, z_clamp);
    const size_t n_clamp = clamp_part.body_part_particles_.size();
    if (n_clamp == 0 || n_clamp >= n)
    {
        std::cerr << "[ophelie][doc228][E0] FATAL: clamp region empty or whole body "
                  << "n_clamp=" << n_clamp << " n=" << n << " z_clamp=" << z_clamp << "\n";
        return 2;
    }

    InnerRelation inner(paddle);
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected(inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_first(inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_second(inner);
    SimpleDynamics<FixBodyPartConstraint> fix_clamp(clamp_part);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_size(paddle);
    SimpleDynamics<VonMisesStress> von_mises(paddle);

    BodyStatesRecordingToVtp vtp(system);
    vtp.addDerivedVariableRecording<SimpleDynamics<Displacement>>(paddle);
    vtp.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(paddle);

    corrected.exec();
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    const Real dt0 = dt_size.exec();

    std::cout << "[ophelie][doc228][E0] paddle hold stl=" << cli.stl_path << " dp=" << cli.dp
              << " n=" << n << " n_clamp=" << n_clamp << " clamp_frac=" << (Real(n_clamp) / Real(n))
              << " z_clamp=" << z_clamp << " E=" << steel.youngs_modulus_pa << " dt0=" << dt0
              << " steps=" << cli.solid_steps << "\n";
    std::cout << "[ophelie][doc228][E0] note: clamp is provisional top-" << cli.clamp_top_band_m
              << " m band; verify drive section before E2\n";

    vtp.writeToFile(0);

    Vecd *pos = paddle.getBaseParticles().getVariableDataByName<Vecd>("Position");
    StdVec<Vecd> pos0(n);
    for (size_t i = 0; i < n; ++i)
        pos0[i] = pos[i];

    Real max_disp = Real(0);
    Real max_von = Real(0);
    bool finite_ok = true;
    const TickCount t1 = TickCount::now();

    for (int ite = 0; ite < cli.solid_steps; ++ite)
    {
        const Real dt = dt_size.exec();
        stress_first.exec(dt);
        fix_clamp.exec();
        stress_second.exec(dt);
        physical_time += dt;
        if ((ite + 1) % 20 == 0)
            std::cout << "[ophelie][doc228][E0] N=" << (ite + 1) << " t=" << physical_time << " dt=" << dt << "\n";
    }

    von_mises.exec();
    Real *von = paddle.getBaseParticles().getVariableDataByName<Real>("VonMisesStress");
    for (size_t i = 0; i < n; ++i)
    {
        const Real d = (pos[i] - pos0[i]).norm();
        max_disp = std::max(max_disp, d);
        max_von = std::max(max_von, std::abs(von[i]));
        if (!std::isfinite(static_cast<double>(pos[i].norm())) || !std::isfinite(static_cast<double>(von[i])))
            finite_ok = false;
    }

    vtp.writeToFile(cli.solid_steps);
    const Real wall_s = (TickCount::now() - t1).seconds();
    const Real disp_over_dp = (cli.dp > TinyReal) ? max_disp / cli.dp : Real(0);
    const Real sigma_over_E =
        (steel.youngs_modulus_pa > TinyReal) ? max_von / steel.youngs_modulus_pa : Real(0);

    // Doc228 E0: no NaN, no spontaneous large motion; residual OK if << dp and sigma/E tiny.
    const bool clamp_ok = (n_clamp > 0) && (n_clamp < n);
    const bool quiet_ok = (disp_over_dp < Real(0.5)) && (sigma_over_E < Real(1.0e-3));
    const bool pass = finite_ok && clamp_ok && quiet_ok && (cli.solid_steps >= 100);

    std::cout << "[ophelie][doc228][E0] max_disp=" << max_disp << " disp/dp=" << disp_over_dp
              << " max_vonMises=" << max_von << " sigma/E=" << sigma_over_E
              << " finite=" << (finite_ok ? 1 : 0) << " wall_s=" << wall_s
              << " PADDLE_HOLD_PASS=" << (pass ? 1 : 0) << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc228_e0_paddle_hold\",\n";
        js << "  \"tag\": \"PADDLE_HOLD\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"PADDLE_HOLD_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"stl\": \"" << cli.stl_path << "\",\n";
        js << "  \"dp_m\": " << cli.dp << ",\n";
        js << "  \"n_particles\": " << n << ",\n";
        js << "  \"n_clamp\": " << n_clamp << ",\n";
        js << "  \"clamp_top_band_m\": " << cli.clamp_top_band_m << ",\n";
        js << "  \"clamp_status\": \"provisional_top_band_unconfirmed_drive_section\",\n";
        js << "  \"E_Pa\": " << steel.youngs_modulus_pa << ",\n";
        js << "  \"rho\": " << steel.rho_kg_m3 << ",\n";
        js << "  \"poisson\": " << steel.poisson << ",\n";
        js << "  \"solid_steps\": " << cli.solid_steps << ",\n";
        js << "  \"physical_time\": " << physical_time << ",\n";
        js << "  \"max_disp_m\": " << max_disp << ",\n";
        js << "  \"disp_over_dp\": " << disp_over_dp << ",\n";
        js << "  \"max_von_Mises\": " << max_von << ",\n";
        js << "  \"sigma_over_E\": " << sigma_over_E << ",\n";
        js << "  \"finite_ok\": " << (finite_ok ? 1 : 0) << ",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"note\": \"E0 isothermal hold; no fluid/EM/gravity; stainless E not softened\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][doc228][E0] summary_json=" << summary_path.string() << "\n";

    // Soft warn if blade thickness under-resolved (Doc228 ~3–4 layers).
    if (cli.dp > Real(0.003))
    {
        std::cout << "[ophelie][doc228][E0] WARN: coarse dp=" << cli.dp
                  << " may be BLOCKED_UNDERRESOLVED_SOLID for thin blades; refine for E1 stress claims\n";
    }
    return pass ? 0 : 1;
#endif
}
