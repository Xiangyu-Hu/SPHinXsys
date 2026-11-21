/**
 * @file test_3d_ophelie_doc228_e2_paddle_spin.cpp
 * @brief Doc228 E2 — client paddle prescribed clamp rotation (no fluid).
 *
 * Smooth ω ramp at client 10 rpm; clamp top-band driven; free zone integrated.
 * Short smoke window (default 0.05 s). Full 2π revolution is a separate later run.
 * Stainless E not softened. Tag: DRIVEN_ELASTIC_ROTATION_SHORT_PASS.
 */
#include "electromagnetic_ophelie_client_stirring_geometry.h"
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

struct E2Cli
{
    std::string stl_path = "./data/client_stir_paddle_solver_m.stl";
    Real dp = Real(0.004);
    Real end_time = Real(0.05); // Doc228 short smoke
    Real omega_ramp_s = Real(0.05);
    Real clamp_top_band_m = Real(0.05);
    Real max_wall_hours = Real(2.0);
    std::string out_dir = ".";
    std::string summary_json = "doc228_e2_summary.json";
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

inline Vecd rotateAboutAxis(const Vecd &rel, const Vecd &axis, Real c, Real s)
{
    return rel * c + axis.cross(rel) * s + axis * (axis.dot(rel) * (Real(1) - c));
}

inline void omegaThetaRamp(Real t, Real tramp, Real omega_max, Real &omega, Real &theta)
{
    if (tramp <= TinyReal)
    {
        omega = omega_max;
        theta = omega_max * t;
        return;
    }
    if (t < tramp)
    {
        omega = omega_max * (t / tramp);
        theta = omega_max * (t * t) / (Real(2) * tramp);
    }
    else
    {
        omega = omega_max;
        theta = omega_max * (Real(0.5) * tramp + (t - tramp));
    }
}

/** Host-only clamp drive (no FSI / no AverageVelocity). */
inline void applyClampSpinDriveHost(BaseParticles &particles, const StdVec<size_t> &clamp_ids,
                                    const StdVec<Vecd> &pos0, const Vecd &center, const Vecd &axis,
                                    Real theta, Real omega)
{
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd e = axis.normalized();
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    for (size_t id : clamp_ids)
    {
        const Vecd rotated = rotateAboutAxis(pos0[id] - center, e, c, s);
        pos[id] = center + rotated;
        vel[id] = e.cross(rotated) * omega;
    }
}

void applyCli(int ac, char *av[], E2Cli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--stl=", 6) == 0)
            cli.stl_path = av[i] + 6;
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        else if (std::strncmp(av[i], "--end-time=", 11) == 0)
            cli.end_time = static_cast<Real>(std::atof(av[i] + 11));
        else if (std::strncmp(av[i], "--omega-ramp=", 13) == 0)
            cli.omega_ramp_s = static_cast<Real>(std::atof(av[i] + 13));
        else if (std::strncmp(av[i], "--clamp-top-band=", 17) == 0)
            cli.clamp_top_band_m = static_cast<Real>(std::atof(av[i] + 17));
        else if (std::strncmp(av[i], "--max-wall-hours=", 17) == 0)
            cli.max_wall_hours = static_cast<Real>(std::atof(av[i] + 17));
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            cli.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            cli.summary_json = av[i] + 15;
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc228 E2 paddle spin: --stl= --dp= --end-time= --omega-ramp= "
                         "--clamp-top-band= --max-wall-hours= --out-dir= --summary-json=\n";
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
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_doc228_e2_paddle_spin requires SYCL build.\n";
    return 1;
#else
    E2Cli cli;
    applyCli(ac, av, cli);
    fs::create_directories(cli.out_dir);

    if (!fs::exists(cli.stl_path))
    {
        std::cerr << "[ophelie][doc228][E2] FATAL: missing STL \"" << cli.stl_path << "\"\n";
        return 2;
    }

    OphelieClientStirringCaseParams g_case;
    applyClientStirringDevelopmentDefaults(g_case);
    const OphelieClientStainless304Development steel = g_case.steel;
    writeClientStainless304Json(steel, (fs::path(cli.out_dir) / "stainless_304_dev_preset.json").string());

    const Real omega_max = g_case.french_stirring.rotation_speed_rad_s; // 10 rpm
    const Vecd axis = g_case.french_stirring.rotation_axis;
    Vecd center = g_case.french_stirring.rotation_center;

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
        std::cerr << "[ophelie][doc228][E2] FATAL: zero paddle particles\n";
        return 2;
    }

    const Real z_max = mesh->getBounds().upper_[2];
    const Real z_clamp = z_max - cli.clamp_top_band_m;
    ClampTopBand clamp_part(paddle, z_clamp);
    const StdVec<size_t> clamp_ids = clamp_part.body_part_particles_;
    const size_t n_clamp = clamp_ids.size();
    if (n_clamp == 0 || n_clamp >= n)
    {
        std::cerr << "[ophelie][doc228][E2] FATAL: clamp empty or whole body n_clamp=" << n_clamp
                  << " n=" << n << "\n";
        return 2;
    }

    InnerRelation inner(paddle);
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected(inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_first(inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_second(inner);
    const Real visc = physicalViscosity(steel.rho_kg_m3, steel.youngs_modulus_pa, cli.dp);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping(
        Real(0.2), inner, "Velocity", visc);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_size(paddle);
    SimpleDynamics<VonMisesStress> von_mises(paddle);

    corrected.exec();

    Vecd *pos = paddle.getBaseParticles().getVariableDataByName<Vecd>("Position");
    StdVec<Vecd> pos0(n);
    for (size_t i = 0; i < n; ++i)
        pos0[i] = pos[i];

    // Prefer clamp-centroid XY if client eccentric axis is far from shaft particles.
    {
        Vecd clamp_c = Vecd::Zero();
        for (size_t id : clamp_ids)
            clamp_c += pos0[id];
        clamp_c /= Real(n_clamp);
        const Real dx = clamp_c[0] - center[0];
        const Real dy = clamp_c[1] - center[1];
        if (std::sqrt(dx * dx + dy * dy) > Real(2) * cli.dp)
        {
            std::cout << "[ophelie][doc228][E2] WARN: client rotation_center XY=(" << center[0] << ","
                      << center[1] << ") far from clamp centroid (" << clamp_c[0] << "," << clamp_c[1]
                      << "); using clamp centroid XY for drive axis\n";
            center[0] = clamp_c[0];
            center[1] = clamp_c[1];
        }
        center[2] = clamp_c[2];
    }

    BodyStatesRecordingToVtp vtp(system);
    vtp.addDerivedVariableRecording<SimpleDynamics<Displacement>>(paddle);
    vtp.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(paddle);
    vtp.writeToFile(0);

    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    const Real wall_budget_s = cli.max_wall_hours * Real(3600);
    const TickCount t1 = TickCount::now();

    std::cout << "[ophelie][doc228][E2] paddle spin stl=" << cli.stl_path << " dp=" << cli.dp
              << " n=" << n << " n_clamp=" << n_clamp << " E=" << steel.youngs_modulus_pa
              << " rpm=" << g_case.client.stirrer_rpm << " omega_max=" << omega_max
              << " ramp=" << cli.omega_ramp_s << " end=" << cli.end_time << " center=(" << center[0]
              << "," << center[1] << "," << center[2] << ")\n";
    std::cout << "[ophelie][doc228][E2] note: short smoke only; full 2π revolution not claimed\n";

    std::ofstream csv((fs::path(cli.out_dir) / "e2_drive.csv").string());
    csv << "step,physical_time_s,omega,theta,dt,max_e_off,clamp_err,tip_speed\n";

    int ite = 0;
    Real last_omega = 0;
    Real last_theta = 0;
    Real max_e_off = 0;
    Real max_clamp_err = 0;
    Real max_tip_speed = 0;
    Real max_von = 0;
    bool finite_ok = true;
    bool budget_exhausted = false;
    bool sign_ok = false;

    const Vecd e_axis = axis.normalized();

    while (physical_time < cli.end_time - TinyReal)
    {
        const Real wall_s = (TickCount::now() - t1).seconds();
        if (wall_s > wall_budget_s)
        {
            budget_exhausted = true;
            break;
        }

        Real dt = dt_size.exec();
        if (!(dt > TinyReal) || !std::isfinite(static_cast<double>(dt)))
        {
            finite_ok = false;
            break;
        }

        Real omega = 0;
        Real theta = 0;
        omegaThetaRamp(physical_time + Real(0.5) * dt, cli.omega_ramp_s, omega_max, omega, theta);

        stress_first.exec(dt);
        applyClampSpinDriveHost(paddle.getBaseParticles(), clamp_ids, pos0, center, axis, theta, omega);
        damping.exec(dt);
        applyClampSpinDriveHost(paddle.getBaseParticles(), clamp_ids, pos0, center, axis, theta, omega);
        stress_second.exec(dt);
        applyClampSpinDriveHost(paddle.getBaseParticles(), clamp_ids, pos0, center, axis, theta, omega);

        physical_time += dt;
        ++ite;
        last_omega = omega;
        last_theta = theta;

        if (ite % 2000 == 0 || physical_time >= cli.end_time - TinyReal)
        {
            const Real c = std::cos(theta);
            const Real s = std::sin(theta);
            Real e_off = 0;
            Real clamp_err = 0;
            Real tip_speed = 0;
            Vecd *vel = paddle.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
            for (size_t i = 0; i < n; ++i)
            {
                if (!std::isfinite(static_cast<double>(pos[i].norm())))
                {
                    finite_ok = false;
                    break;
                }
                const Vecd rigid = center + rotateAboutAxis(pos0[i] - center, e_axis, c, s);
                const Real d = (pos[i] - rigid).norm();
                tip_speed = std::max(tip_speed, vel[i].norm());
                // classify clamp vs free by membership
                bool is_clamp = false;
                // linear scan rare path — clamp set small; use z for speed
                if (pos0[i][2] >= z_clamp)
                {
                    is_clamp = true;
                    clamp_err = std::max(clamp_err, d);
                }
                else
                    e_off = std::max(e_off, d);
                (void)is_clamp;
            }
            max_e_off = std::max(max_e_off, e_off);
            max_clamp_err = std::max(max_clamp_err, clamp_err);
            max_tip_speed = std::max(max_tip_speed, tip_speed);

            // Sign check: free particle with largest |r_xy| should have v · (e × r_xy_hat) > 0 when ω>0
            if (!sign_ok && omega > TinyReal)
            {
                Real best_r = 0;
                size_t best_i = 0;
                for (size_t i = 0; i < n; ++i)
                {
                    if (pos0[i][2] >= z_clamp)
                        continue;
                    const Vecd rel = pos[i] - center;
                    const Real r_xy = (rel - e_axis * e_axis.dot(rel)).norm();
                    if (r_xy > best_r)
                    {
                        best_r = r_xy;
                        best_i = i;
                    }
                }
                if (best_r > TinyReal)
                {
                    const Vecd rel = pos[best_i] - center;
                    const Vecd r_xy = rel - e_axis * e_axis.dot(rel);
                    const Vecd v_expect_dir = e_axis.cross(r_xy);
                    sign_ok = vel[best_i].dot(v_expect_dir) > Real(0);
                }
            }

            csv << ite << "," << physical_time << "," << omega << "," << theta << "," << dt << "," << e_off
                << "," << clamp_err << "," << tip_speed << "\n";
            std::cout << "[ophelie][doc228][E2] N=" << ite << " t=" << physical_time << " omega=" << omega
                      << " theta=" << theta << " e_off=" << e_off << " clamp_err=" << clamp_err
                      << " tip_u=" << tip_speed << "\n";
        }
    }

    von_mises.exec();
    Real *von = paddle.getBaseParticles().getVariableDataByName<Real>("VonMisesStress");
    for (size_t i = 0; i < n; ++i)
        max_von = std::max(max_von, std::abs(von[i]));

    vtp.writeToFile(ite);
    const Real wall_s = (TickCount::now() - t1).seconds();

    const bool time_ok = (physical_time >= Real(0.9) * cli.end_time) && !budget_exhausted;
    const bool clamp_ok = (max_clamp_err < Real(0.05) * cli.dp);
    const bool motion_ok = (last_theta > Real(0.01)) && (last_omega > TinyReal) && (max_tip_speed > TinyReal);
    const bool pass = finite_ok && time_ok && clamp_ok && motion_ok && sign_ok && (ite > 100);

    std::cout << "[ophelie][doc228][E2] t=" << physical_time << " theta=" << last_theta
              << " omega=" << last_omega << " max_e_off=" << max_e_off
              << " max_clamp_err=" << max_clamp_err << " tip_u=" << max_tip_speed
              << " von=" << max_von << " steps=" << ite << " wall_s=" << wall_s
              << " sign_ok=" << (sign_ok ? 1 : 0)
              << " DRIVEN_ELASTIC_ROTATION_SHORT_PASS=" << (pass ? 1 : 0) << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc228_e2_paddle_spin\",\n";
        js << "  \"tag\": \"DRIVEN_ELASTIC_ROTATION_SHORT\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"DRIVEN_ELASTIC_ROTATION_SHORT_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"full_revolution_pass\": 0,\n";
        js << "  \"full_revolution_note\": \"not claimed; short smoke only\",\n";
        js << "  \"stl\": \"" << cli.stl_path << "\",\n";
        js << "  \"dp_m\": " << cli.dp << ",\n";
        js << "  \"n_particles\": " << n << ",\n";
        js << "  \"n_clamp\": " << n_clamp << ",\n";
        js << "  \"clamp_top_band_m\": " << cli.clamp_top_band_m << ",\n";
        js << "  \"clamp_status\": \"provisional_top_band_unconfirmed_drive_section\",\n";
        js << "  \"E_Pa\": " << steel.youngs_modulus_pa << ",\n";
        js << "  \"stirrer_rpm\": " << g_case.client.stirrer_rpm << ",\n";
        js << "  \"omega_max_rad_s\": " << omega_max << ",\n";
        js << "  \"omega_ramp_s\": " << cli.omega_ramp_s << ",\n";
        js << "  \"end_time_s\": " << cli.end_time << ",\n";
        js << "  \"physical_time_s\": " << physical_time << ",\n";
        js << "  \"theta_rad\": " << last_theta << ",\n";
        js << "  \"omega_final_rad_s\": " << last_omega << ",\n";
        js << "  \"rotation_center\": [" << center[0] << ", " << center[1] << ", " << center[2] << "],\n";
        js << "  \"rotation_axis\": [" << axis[0] << ", " << axis[1] << ", " << axis[2] << "],\n";
        js << "  \"max_elastic_offset_m\": " << max_e_off << ",\n";
        js << "  \"max_clamp_err_m\": " << max_clamp_err << ",\n";
        js << "  \"max_tip_speed_m_s\": " << max_tip_speed << ",\n";
        js << "  \"max_von_Mises\": " << max_von << ",\n";
        js << "  \"solid_steps\": " << ite << ",\n";
        js << "  \"sign_ok\": " << (sign_ok ? 1 : 0) << ",\n";
        js << "  \"finite_ok\": " << (finite_ok ? 1 : 0) << ",\n";
        js << "  \"budget_exhausted\": " << (budget_exhausted ? 1 : 0) << ",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"note\": \"E2 no fluid; clamp-driven 10 rpm ramp; stainless E not softened\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][doc228][E2] summary_json=" << summary_path.string() << "\n";
    return pass ? 0 : 1;
#endif
}
