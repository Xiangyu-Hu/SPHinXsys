/**
 * @file test_3d_ophelie_doc228_e1_paddle_load.cpp
 * @brief Doc228 E1 — client stirring paddle non-zero load response (no fluid).
 *
 * Two known total-force levels on the free (non-clamp) region; mass-weighted into
 * ForcePrior. Stainless E not softened. Checks resolvable elastic offset,
 * stress update, clamp stillness, and approximate force–offset proportionality.
 *
 * STL: client_stir_paddle_solver_m.stl from stirring relax data/ (CMake copy).
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
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{
Real physicalViscosity(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = Real(0.4))
{
    return shape_constant / Real(4) * std::sqrt(rho * youngs_modulus) * length_scale;
}

struct E1Cli
{
    std::string stl_path = "./data/client_stir_paddle_solver_m.stl";
    Real dp = Real(0.004);
    int solid_steps = 5000;
    Real clamp_top_band_m = Real(0.05);
    Real force1_n = Real(1.5e4); // total force on free region (level 1)
    Real force2_n = Real(3.0e4); // level 2 (= 2× default); keep sigma/E ≲ 1e-2
    Real force_ramp_end = Real(1.0e-4);
    Real ratio_rel_tol = Real(0.45); // transient window; not quasi-static analytic
    Real min_e_off_over_dp = Real(1.0e-4);
    Real max_sigma_over_E = Real(1.5e-2); // Doc228 small-strain gate (not strength rating)
    std::string out_dir = ".";
    std::string summary_json = "doc228_e1_summary.json";
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

class FreeLoadBand : public BodyPartByParticle
{
  public:
    FreeLoadBand(SPHBody &body, Real z_max) : BodyPartByParticle(body), z_max_(z_max)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FreeLoadBand::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    }

  private:
    Real z_max_;
    bool tagManually(size_t index_i) { return pos_[index_i][2] < z_max_; }
};

void applyCli(int ac, char *av[], E1Cli &cli)
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
        else if (std::strncmp(av[i], "--force1=", 9) == 0)
            cli.force1_n = static_cast<Real>(std::atof(av[i] + 9));
        else if (std::strncmp(av[i], "--force2=", 9) == 0)
            cli.force2_n = static_cast<Real>(std::atof(av[i] + 9));
        else if (std::strncmp(av[i], "--force-ramp=", 13) == 0)
            cli.force_ramp_end = static_cast<Real>(std::atof(av[i] + 13));
        else if (std::strncmp(av[i], "--ratio-tol=", 12) == 0)
            cli.ratio_rel_tol = static_cast<Real>(std::atof(av[i] + 12));
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            cli.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            cli.summary_json = av[i] + 15;
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc228 E1 paddle load: --stl= --dp= --solid-steps= --clamp-top-band= "
                         "--force1= --force2= --force-ramp= --ratio-tol= --out-dir= --summary-json=\n";
            std::exit(0);
        }
    }
}

class ClientPaddleShape : public ComplexShape
{
  public:
    ClientPaddleShape(const std::string &name, const std::string &stl_path) : ComplexShape(name)
    {
        add<TriangleMeshShapeSTL>(stl_path, Vecd::Zero(), Real(1), "ClientPaddle");
    }
};

struct LoadMetrics
{
    Real force_n = Real(0);
    Real max_e_off_free = Real(0);
    Real max_e_off_clamp = Real(0);
    Real max_von = Real(0);
    Real physical_time = Real(0);
    bool finite_ok = true;
};

void zeroField(Vecd *v, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        v[i] = Vecd::Zero();
}

void resetSolidState(BaseParticles &particles, const StdVec<Vecd> &pos0, Real rho0)
{
    const size_t n = particles.TotalRealParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *force = particles.getVariableDataByName<Vecd>("Force");
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Matd *F = particles.getVariableDataByName<Matd>("DeformationGradient");
    Matd *dF = particles.getVariableDataByName<Matd>("DeformationRate");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Matd *stress = particles.getVariableDataByName<Matd>("StressPK1OnParticle");

    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = pos0[i];
        vel[i] = Vecd::Zero();
        force[i] = Vecd::Zero();
        fp[i] = Vecd::Zero();
        F[i] = Matd::Identity();
        dF[i] = Matd::Zero();
        rho[i] = rho0;
        stress[i] = Matd::Zero();
    }
}

void applyMassWeightedForce(BaseParticles &particles, const StdVec<size_t> &load_ids, const Vecd &F_total,
                            Real ramp)
{
    const size_t n = particles.TotalRealParticles();
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    zeroField(fp, n);

    Real M = Real(0);
    for (size_t id : load_ids)
        M += mass[id];
    if (!(M > TinyReal) || !(ramp > TinyReal))
        return;

    const Vecd F = F_total * ramp;
    for (size_t id : load_ids)
        fp[id] = F * (mass[id] / M);
}

LoadMetrics runLoadLevel(SPHSystem &system, SolidBody &paddle, InnerRelation &inner,
                         const StdVec<size_t> &load_ids, const StdVec<size_t> &clamp_ids,
                         const StdVec<Vecd> &pos0, Real rho0, Real force_n, const Vecd &force_dir,
                         const E1Cli &cli,
                         Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> &stress_first,
                         Dynamics1Level<solid_dynamics::Integration2ndHalf> &stress_second,
                         SimpleDynamics<FixBodyPartConstraint> &fix_clamp,
                         DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> &damping,
                         ReduceDynamics<solid_dynamics::AcousticTimeStep> &dt_size,
                         SimpleDynamics<VonMisesStress> &von_mises, BodyStatesRecordingToVtp &vtp,
                         int vtp_tag)
{
    resetSolidState(paddle.getBaseParticles(), pos0, rho0);
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    physical_time = Real(0);

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected(inner);
    corrected.exec();
    vtp.writeToFile(vtp_tag);

    const Vecd F_total = force_dir * force_n;
    const size_t n = paddle.getBaseParticles().TotalRealParticles();
    Vecd *pos = paddle.getBaseParticles().getVariableDataByName<Vecd>("Position");

    for (int ite = 0; ite < cli.solid_steps; ++ite)
    {
        const Real ramp =
            (cli.force_ramp_end > TinyReal) ? std::min(Real(1), physical_time / cli.force_ramp_end) : Real(1);
        applyMassWeightedForce(paddle.getBaseParticles(), load_ids, F_total, ramp);

        const Real dt = dt_size.exec();
        stress_first.exec(dt);
        fix_clamp.exec();
        damping.exec(dt);
        fix_clamp.exec();
        stress_second.exec(dt);
        physical_time += dt;

        if ((ite + 1) % 1000 == 0)
        {
            std::cout << "[ophelie][doc228][E1] F=" << force_n << " N=" << (ite + 1)
                      << " t=" << physical_time << " dt=" << dt << " ramp=" << ramp << "\n";
        }
    }

    von_mises.exec();
    Real *von = paddle.getBaseParticles().getVariableDataByName<Real>("VonMisesStress");

    LoadMetrics m;
    m.force_n = force_n;
    m.physical_time = physical_time;
    for (size_t id : load_ids)
    {
        const Real d = (pos[id] - pos0[id]).norm();
        m.max_e_off_free = std::max(m.max_e_off_free, d);
        m.max_von = std::max(m.max_von, std::abs(von[id]));
        if (!std::isfinite(static_cast<double>(pos[id].norm())) || !std::isfinite(static_cast<double>(von[id])))
            m.finite_ok = false;
    }
    for (size_t id : clamp_ids)
    {
        const Real d = (pos[id] - pos0[id]).norm();
        m.max_e_off_clamp = std::max(m.max_e_off_clamp, d);
        if (!std::isfinite(static_cast<double>(pos[id].norm())))
            m.finite_ok = false;
    }
    for (size_t i = 0; i < n; ++i)
        m.max_von = std::max(m.max_von, std::abs(von[i]));

    vtp.writeToFile(vtp_tag + 1);
    return m;
}
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_doc228_e1_paddle_load requires SYCL build.\n";
    return 1;
#else
    E1Cli cli;
    applyCli(ac, av, cli);
    fs::create_directories(cli.out_dir);

    if (!fs::exists(cli.stl_path))
    {
        std::cerr << "[ophelie][doc228][E1] FATAL: missing STL \"" << cli.stl_path
                  << "\". Rebuild after stirring relax data/ has client_stir_paddle_solver_m.stl.\n";
        return 2;
    }
    if (!(cli.force1_n > TinyReal) || !(cli.force2_n > TinyReal) || !(cli.force2_n > cli.force1_n))
    {
        std::cerr << "[ophelie][doc228][E1] FATAL: need 0 < force1 < force2\n";
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
    paddle.defineMatterMaterial<NeoHookeanSolid>(steel.rho_kg_m3, steel.youngs_modulus_pa, steel.poisson);
    paddle.generateParticles<BaseParticles, Lattice>();

    const size_t n = paddle.getBaseParticles().TotalRealParticles();
    if (n == 0)
    {
        std::cerr << "[ophelie][doc228][E1] FATAL: zero paddle particles\n";
        return 2;
    }

    const Real z_max = mesh->getBounds().upper_[2];
    const Real z_clamp = z_max - cli.clamp_top_band_m;
    ClampTopBand clamp_part(paddle, z_clamp);
    FreeLoadBand load_part(paddle, z_clamp);
    const size_t n_clamp = clamp_part.body_part_particles_.size();
    const size_t n_load = load_part.body_part_particles_.size();
    if (n_clamp == 0 || n_clamp >= n || n_load == 0)
    {
        std::cerr << "[ophelie][doc228][E1] FATAL: bad clamp/load split n_clamp=" << n_clamp
                  << " n_load=" << n_load << " n=" << n << "\n";
        return 2;
    }

    InnerRelation inner(paddle);
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected(inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_first(inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_second(inner);
    SimpleDynamics<FixBodyPartConstraint> fix_clamp(clamp_part);
    const Real visc = physicalViscosity(steel.rho_kg_m3, steel.youngs_modulus_pa, cli.dp);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping(
        Real(0.2), inner, "Velocity", visc);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_size(paddle);
    SimpleDynamics<VonMisesStress> von_mises(paddle);

    corrected.exec();

    BodyStatesRecordingToVtp vtp(system);
    vtp.addDerivedVariableRecording<SimpleDynamics<Displacement>>(paddle);
    vtp.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(paddle);

    Vecd *pos = paddle.getBaseParticles().getVariableDataByName<Vecd>("Position");
    StdVec<Vecd> pos0(n);
    for (size_t i = 0; i < n; ++i)
        pos0[i] = pos[i];

    const StdVec<size_t> load_ids = load_part.body_part_particles_;
    const StdVec<size_t> clamp_ids = clamp_part.body_part_particles_;
    const Vecd force_dir = Vecd(0, 1, 0); // lateral +Y

    std::cout << "[ophelie][doc228][E1] paddle load stl=" << cli.stl_path << " dp=" << cli.dp
              << " n=" << n << " n_clamp=" << n_clamp << " n_load=" << n_load
              << " E=" << steel.youngs_modulus_pa << " F1=" << cli.force1_n << " F2=" << cli.force2_n
              << " steps=" << cli.solid_steps << " visc=" << visc << "\n";
    std::cout << "[ophelie][doc228][E1] note: stainless E not softened; two load levels on free band "
                 "z<"
              << z_clamp << "\n";

    const TickCount t1 = TickCount::now();
    const LoadMetrics m1 =
        runLoadLevel(system, paddle, inner, load_ids, clamp_ids, pos0, steel.rho_kg_m3, cli.force1_n,
                     force_dir, cli, stress_first, stress_second, fix_clamp, damping, dt_size, von_mises, vtp,
                     10);
    const LoadMetrics m2 =
        runLoadLevel(system, paddle, inner, load_ids, clamp_ids, pos0, steel.rho_kg_m3, cli.force2_n,
                     force_dir, cli, stress_first, stress_second, fix_clamp, damping, dt_size, von_mises, vtp,
                     20);
    const Real wall_s = (TickCount::now() - t1).seconds();

    const Real force_ratio = cli.force2_n / cli.force1_n;
    const Real e_off_ratio =
        (m1.max_e_off_free > TinyReal) ? (m2.max_e_off_free / m1.max_e_off_free) : Real(0);
    const Real ratio_rel_err =
        (force_ratio > TinyReal) ? std::abs(e_off_ratio - force_ratio) / force_ratio : Real(1);
    const Real sigma_over_E_1 =
        (steel.youngs_modulus_pa > TinyReal) ? m1.max_von / steel.youngs_modulus_pa : Real(0);
    const Real sigma_over_E_2 =
        (steel.youngs_modulus_pa > TinyReal) ? m2.max_von / steel.youngs_modulus_pa : Real(0);
    const Real e_off1_over_dp = (cli.dp > TinyReal) ? m1.max_e_off_free / cli.dp : Real(0);
    const Real e_off2_over_dp = (cli.dp > TinyReal) ? m2.max_e_off_free / cli.dp : Real(0);

    const bool finite_ok = m1.finite_ok && m2.finite_ok;
    const bool response_ok =
        (e_off1_over_dp >= cli.min_e_off_over_dp) && (e_off2_over_dp >= cli.min_e_off_over_dp) &&
        (m1.max_von > TinyReal) && (m2.max_von > TinyReal);
    const bool clamp_ok =
        (m1.max_e_off_clamp < Real(0.05) * cli.dp) && (m2.max_e_off_clamp < Real(0.05) * cli.dp);
    const bool small_strain_ok =
        (sigma_over_E_1 < cli.max_sigma_over_E) && (sigma_over_E_2 < cli.max_sigma_over_E);
    const bool proportional_ok = (ratio_rel_err <= cli.ratio_rel_tol);
    const bool pass = finite_ok && response_ok && clamp_ok && small_strain_ok && proportional_ok &&
                      (cli.solid_steps >= 100);

    std::cout << "[ophelie][doc228][E1] F1 e_off_free=" << m1.max_e_off_free
              << " e_off/dp=" << e_off1_over_dp << " von=" << m1.max_von << " sigma/E=" << sigma_over_E_1
              << " clamp_e_off=" << m1.max_e_off_clamp << "\n";
    std::cout << "[ophelie][doc228][E1] F2 e_off_free=" << m2.max_e_off_free
              << " e_off/dp=" << e_off2_over_dp << " von=" << m2.max_von << " sigma/E=" << sigma_over_E_2
              << " clamp_e_off=" << m2.max_e_off_clamp << "\n";
    std::cout << "[ophelie][doc228][E1] force_ratio=" << force_ratio << " e_off_ratio=" << e_off_ratio
              << " ratio_rel_err=" << ratio_rel_err << " wall_s=" << wall_s
              << " LOAD_RESPONSE_PASS=" << (pass ? 1 : 0)
              << " ELASTIC_PADDLE_ACTIVE=" << (pass ? 1 : 0) << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc228_e1_paddle_load\",\n";
        js << "  \"tag\": \"LOAD_RESPONSE\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"LOAD_RESPONSE_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"ELASTIC_PADDLE_ACTIVE\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"operator_path\": \"NeoHookeanSolid+Integration1stHalfPK2RightCauchy+"
              "massWeightedForcePrior+FixBodyPartConstraint\",\n";
        js << "  \"stl\": \"" << cli.stl_path << "\",\n";
        js << "  \"dp_m\": " << cli.dp << ",\n";
        js << "  \"n_particles\": " << n << ",\n";
        js << "  \"n_clamp\": " << n_clamp << ",\n";
        js << "  \"n_load\": " << n_load << ",\n";
        js << "  \"clamp_top_band_m\": " << cli.clamp_top_band_m << ",\n";
        js << "  \"clamp_status\": \"provisional_top_band_unconfirmed_drive_section\",\n";
        js << "  \"E_Pa\": " << steel.youngs_modulus_pa << ",\n";
        js << "  \"rho\": " << steel.rho_kg_m3 << ",\n";
        js << "  \"poisson\": " << steel.poisson << ",\n";
        js << "  \"force_direction\": [0,1,0],\n";
        js << "  \"force1_N\": " << m1.force_n << ",\n";
        js << "  \"force2_N\": " << m2.force_n << ",\n";
        js << "  \"force_ratio\": " << force_ratio << ",\n";
        js << "  \"solid_steps_per_level\": " << cli.solid_steps << ",\n";
        js << "  \"physical_time_level1\": " << m1.physical_time << ",\n";
        js << "  \"physical_time_level2\": " << m2.physical_time << ",\n";
        js << "  \"max_e_off_free_1_m\": " << m1.max_e_off_free << ",\n";
        js << "  \"max_e_off_free_2_m\": " << m2.max_e_off_free << ",\n";
        js << "  \"e_off_ratio\": " << e_off_ratio << ",\n";
        js << "  \"ratio_rel_err\": " << ratio_rel_err << ",\n";
        js << "  \"ratio_rel_tol\": " << cli.ratio_rel_tol << ",\n";
        js << "  \"max_e_off_clamp_1_m\": " << m1.max_e_off_clamp << ",\n";
        js << "  \"max_e_off_clamp_2_m\": " << m2.max_e_off_clamp << ",\n";
        js << "  \"max_von_1\": " << m1.max_von << ",\n";
        js << "  \"max_von_2\": " << m2.max_von << ",\n";
        js << "  \"sigma_over_E_1\": " << sigma_over_E_1 << ",\n";
        js << "  \"sigma_over_E_2\": " << sigma_over_E_2 << ",\n";
        js << "  \"finite_ok\": " << (finite_ok ? 1 : 0) << ",\n";
        js << "  \"response_ok\": " << (response_ok ? 1 : 0) << ",\n";
        js << "  \"clamp_ok\": " << (clamp_ok ? 1 : 0) << ",\n";
        js << "  \"small_strain_ok\": " << (small_strain_ok ? 1 : 0) << ",\n";
        js << "  \"proportional_ok\": " << (proportional_ok ? 1 : 0) << ",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"note\": \"E1 isothermal two-level load; no fluid/EM/gravity; E not softened; "
              "not a strength rating\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][doc228][E1] summary_json=" << summary_path.string() << "\n";

    if (cli.dp > Real(0.003))
    {
        std::cout << "[ophelie][doc228][E1] WARN: coarse dp=" << cli.dp
                  << "; stress magnitude not claimed accurate — gate is free-zone response + "
                     "proportionality\n";
    }
    return pass ? 0 : 1;
#endif
}
