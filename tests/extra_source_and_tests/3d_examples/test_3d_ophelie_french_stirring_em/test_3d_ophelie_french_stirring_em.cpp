/**
 * @file test_3d_ophelie_french_stirring_em.cpp
 * @brief French cold-crucible glass melter: induction heating + mechanical stirring.
 *
 * Phase A solves the EM problem once on the reloaded glass particles and deposits the
 * resulting Joule power density onto a fixed Eulerian grid. Phase B runs WCSPH with a
 * Simbody-driven paddle, Boussinesq buoyancy and the French Robin/radiation thermal BC,
 * resampling Q from that grid every advection step.
 *
 * The Eulerian handoff (instead of freezing Q on particles as the natural-convection case
 * does) is what makes the stirred run correct: the inductor is fixed in the laboratory frame
 * while the melt is transported past it.
 *
 * Prerequisites:
 *   1) Run test_3d_ophelie_french_stirring_glass_relax --bodies=all
 *      (GlassBody + Rotor + WallBoundary reload).
 *   2) Point --reload-dir= to that reload folder.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_stirring_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "io_environment.h"
#include "sphinxsys.h"

// Eulerian Q grid (CIC deposit, device trilinear sample) is shared with the RH200 stirring case
// rather than duplicated; it lives next to its rh200_fake_joule_heat.h dependency.
#include "../test_3d_ophelie_rh200_glass_em_stirring/rh200_joule_heat_grid.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

inline std::string resolveDefaultReloadFolder()
{
    const StdVec<std::string> candidates = {
        "./reload",
        "../reload",
        "../../test_3d_ophelie_french_stirring_glass_relax/bin/reload",
        "../../../test_3d_ophelie_french_stirring_glass_relax/bin/reload",
    };
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return candidates.front();
}

struct LocalCli
{
    std::string reload_dir;
    Real end_time = 1200.0;     // 200 paddle revolutions at 10 rpm
    Real max_wall_hours = 10.0; // stop cleanly when the compute budget is spent
    // Upper bound set by stability, not accuracy: the reloaded glass and the separately reloaded
    // rotor/wall packings overlap slightly, and the resulting pressure spikes grow with c0^2.
    // c0 = 10 blows up within 50 acoustic steps; c0 = 8 runs. The residual 3% hydrostatic
    // compression is absorbed by tracking the free surface (see refreshMeltExtentFromParticles).
    Real c0 = 8.0;
    Real q_grid_spacing = 0.0;  // 0 -> dp
    Real t_min = 300.0;         // coolant temperature floor (see README, skull layer)
    Real t_max = 2500.0;
    Real h_side = -1.0;         // <0 -> literature value from OphelieFrenchStirringPhysicsParams
    Real h_bottom = -1.0;
    Real h_free = -1.0;
    bool balance_heat_loss = false;
    int bc_retag_every = 50;    // advection steps between Robin/radiation face re-tagging
    int screen_every = 20;
    int csv_every = 20;
    bool state_recording = true;
    Real state_record_interval = 10.0;
    bool enable_self_induction = true;
    bool enable_boussinesq = true;
};

inline bool isLocalCliOption(const char *arg)
{
    return std::strncmp(arg, "--reload-dir=", 13) == 0 || std::strncmp(arg, "--end-time=", 11) == 0 ||
           std::strncmp(arg, "--max-wall-hours=", 17) == 0 || std::strncmp(arg, "--c0=", 5) == 0 ||
           std::strncmp(arg, "--q-grid-spacing=", 17) == 0 || std::strncmp(arg, "--t-min=", 8) == 0 ||
           std::strncmp(arg, "--t-max=", 8) == 0 || std::strncmp(arg, "--h-side=", 9) == 0 ||
           std::strncmp(arg, "--h-bottom=", 11) == 0 || std::strncmp(arg, "--h-free=", 9) == 0 ||
           std::strcmp(arg, "--balance-heat-loss") == 0 || std::strncmp(arg, "--bc-retag-every=", 17) == 0 ||
           std::strncmp(arg, "--screen-every=", 15) == 0 || std::strncmp(arg, "--state-record-interval=", 24) == 0 ||
           std::strcmp(arg, "--no-state-recording") == 0 || std::strcmp(arg, "--no-self-induction") == 0 ||
           std::strcmp(arg, "--no-boussinesq") == 0;
}

inline void applyLocalCli(int ac, char *av[], LocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--end-time=", 11) == 0)
        {
            cli.end_time = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--max-wall-hours=", 17) == 0)
        {
            cli.max_wall_hours = static_cast<Real>(std::atof(av[i] + 17));
        }
        else if (std::strncmp(av[i], "--c0=", 5) == 0)
        {
            cli.c0 = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--q-grid-spacing=", 17) == 0)
        {
            cli.q_grid_spacing = static_cast<Real>(std::atof(av[i] + 17));
        }
        else if (std::strncmp(av[i], "--t-min=", 8) == 0)
        {
            cli.t_min = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--t-max=", 8) == 0)
        {
            cli.t_max = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--h-side=", 9) == 0)
        {
            cli.h_side = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--h-bottom=", 11) == 0)
        {
            cli.h_bottom = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--h-free=", 9) == 0)
        {
            cli.h_free = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strcmp(av[i], "--balance-heat-loss") == 0)
        {
            cli.balance_heat_loss = true;
        }
        else if (std::strncmp(av[i], "--bc-retag-every=", 17) == 0)
        {
            cli.bc_retag_every = std::atoi(av[i] + 17);
        }
        else if (std::strncmp(av[i], "--screen-every=", 15) == 0)
        {
            cli.screen_every = std::atoi(av[i] + 15);
        }
        else if (std::strncmp(av[i], "--state-record-interval=", 24) == 0)
        {
            cli.state_record_interval = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strcmp(av[i], "--no-state-recording") == 0)
        {
            cli.state_recording = false;
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.enable_self_induction = false;
        }
        else if (std::strcmp(av[i], "--no-boussinesq") == 0)
        {
            cli.enable_boussinesq = false;
        }
    }
    if (cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultReloadFolder();
    }
}

/**
 * Track the free surface so the Robin/radiation face tagging keeps up with the melt.
 *
 * setupOphelieThermalFrenchNaturalBoundaryFaces tags the top face by distance to the nominal
 * plane z = glass_center_z + glass_half_height. The melt settles a few dp under gravity (WCSPH
 * is compressible) and sloshes while stirred, so a fixed plane loses most of the radiating
 * particles within a second of physical time — and radiation carries ~2/3 of the heat budget.
 */
inline void refreshMeltExtentFromParticles(BaseParticles &particles, OphelieFrenchReducedCaseParams &french,
                                           Real z_floor, Real top_quantile = Real(0.995))
{
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        return;
    }
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    StdVec<Real> z_values(n);
    for (size_t i = 0; i < n; ++i)
    {
        z_values[i] = pos[i][2];
    }
    const size_t k = static_cast<size_t>(top_quantile * static_cast<Real>(n - 1));
    std::nth_element(z_values.begin(), z_values.begin() + k, z_values.end());
    const Real z_top = z_values[k];
    if (!(z_top > z_floor))
    {
        return;
    }
    french.glass_half_height = Real(0.5) * (z_top - z_floor);
    french.glass_center[2] = Real(0.5) * (z_top + z_floor);
}

/** Explicit Joule heating that keeps OphelieThermalDeltaT synchronized with Temperature. */
class ApplyGridJouleHeatToTemperatureCK : public LocalDynamics
{
  public:
    ApplyGridJouleHeatToTemperatureCK(SPHBody &sph_body, Real rho, Real cp, Real t_initial)
        : LocalDynamics(sph_body), inv_rho_cp_(Real(1) / (rho * cp + TinyReal)), t_initial_(t_initial),
          dv_q_(particles_->template getVariableByName<Real>("JouleHeat")),
          dv_temperature_(particles_->template getVariableByName<Real>(kOphelieTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : inv_rho_cp_(encloser.inv_rho_cp_), t_initial_(encloser.t_initial_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            temperature_[index_i] += q_[index_i] * dt * inv_rho_cp_;
            delta_t_[index_i] = temperature_[index_i] - t_initial_;
        }

      protected:
        Real inv_rho_cp_;
        Real t_initial_;
        Real *q_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real inv_rho_cp_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

/**
 * Bound the melt temperature.
 *
 * The Robin BC extracts ~170 kW against 60 kW of induction, because the real cold crucible
 * grows an insulating solidified skull that this model does not carry. Without a floor the
 * boundary shell would cool past 0 K within a few seconds of physical time.
 */
class ClampTemperatureCK : public LocalDynamics
{
  public:
    ClampTemperatureCK(SPHBody &sph_body, Real t_min, Real t_max, Real t_initial)
        : LocalDynamics(sph_body), t_min_(t_min), t_max_(t_max), t_initial_(t_initial),
          dv_temperature_(particles_->template getVariableByName<Real>(kOphelieTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_min_(encloser.t_min_), t_max_(encloser.t_max_), t_initial_(encloser.t_initial_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            (void)dt;
            Real t = temperature_[index_i];
            t = t < t_min_ ? t_min_ : t;
            t = t > t_max_ ? t_max_ : t;
            temperature_[index_i] = t;
            delta_t_[index_i] = t - t_initial_;
        }

      protected:
        Real t_min_;
        Real t_max_;
        Real t_initial_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real t_min_;
    Real t_max_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

struct RotorLoadSnapshot
{
    Real tip_speed = 0;
    Real r_tip = 0;
    Real f_visc = 0;
    Real f_pres = 0;
    Real torque_z = 0;
};

/**
 * Rest pose kept on the host so paddle kinematics never depend on a device Position
 * buffer. The previous GPU spin kernel left Rotor_*.vtp with all-nan coordinates:
 * device Position was either never copied back, or InitialPosition on device was
 * uninitialized, and ParaView then renders an empty black view.
 */
struct HostSpinPose
{
    StdVec<Vecd> pos0;
    StdVec<Vecd> n0;
    size_t n_nonfinite = 0;
    Real r_tip = 0;
};

inline Vecd hostRotateAboutAxis(const Vecd &rel, const Vecd &axis, Real c, Real s)
{
    return rel * c + axis.cross(rel) * s + axis * (axis.dot(rel) * (Real(1) - c));
}

inline HostSpinPose captureHostSpinPose(BaseParticles &particles, const Vecd &center, const Vecd &axis,
                                        bool with_normals)
{
    // Read the host reload/STL coordinates only. A device->host sync here would poison
    // the rest pose if an earlier kernel allocated an uninitialized device Position.
    const size_t n = particles.TotalRealParticles();
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    HostSpinPose pose;
    pose.pos0.assign(pos, pos + n);
    const Vecd e_z = axis.normalized();
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(pose.pos0[i][0]) || !std::isfinite(pose.pos0[i][1]) || !std::isfinite(pose.pos0[i][2]))
        {
            ++pose.n_nonfinite;
            continue;
        }
        const Vecd rel = pose.pos0[i] - center;
        pose.r_tip = std::max(pose.r_tip, (rel - e_z * e_z.dot(rel)).norm());
    }
    if (with_normals)
    {
        const Vecd *n_dir = particles.getVariableDataByName<Vecd>("NormalDirection");
        pose.n0.assign(n_dir, n_dir + n);
    }
    return pose;
}

inline void applyHostSpinPose(BaseParticles &particles, const HostSpinPose &pose, const Vecd &center,
                              const Vecd &axis, Real omega, Real physical_time, bool with_normals)
{
    const size_t n = particles.TotalRealParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *n_dir = with_normals ? particles.getVariableDataByName<Vecd>("NormalDirection") : nullptr;
    const Vecd e = axis.normalized();
    const Real theta = omega * physical_time;
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd rotated = hostRotateAboutAxis(pose.pos0[i] - center, e, c, s);
        pos[i] = center + rotated;
        vel[i] = e.cross(rotated) * omega;
        if (n_dir != nullptr)
        {
            n_dir[i] = hostRotateAboutAxis(pose.n0[i], e, c, s);
        }
    }
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    if (n_dir != nullptr)
    {
        syncVariableToDevice<Vecd>(particles, "NormalDirection");
    }
}

inline RotorLoadSnapshot hostRotorLoadSnapshot(BaseParticles &particles, const Vecd &center, const Vecd &axis,
                                               Real omega)
{
    RotorLoadSnapshot snap;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ViscousForceFromFluid");
    syncVariableToHost<Vecd>(particles, "PressureForceFromFluid");
    const size_t n = particles.TotalRealParticles();
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *f_visc = particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    const Vecd *f_pres = particles.getVariableDataByName<Vecd>("PressureForceFromFluid");
    const Vecd e_z = axis.normalized();
    Vecd f_visc_sum = Vecd::Zero();
    Vecd f_pres_sum = Vecd::Zero();
    Real torque = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd rel = pos[i] - center;
        if (!std::isfinite(rel[0]) || !std::isfinite(rel[1]) || !std::isfinite(rel[2]))
        {
            continue;
        }
        const Real r_xy = (rel - e_z * e_z.dot(rel)).norm();
        snap.r_tip = std::max(snap.r_tip, r_xy);
        if (std::isfinite(vel[i][0]))
        {
            snap.tip_speed = std::max(snap.tip_speed, vel[i].norm());
        }
        const Vecd f = f_visc[i] + f_pres[i];
        if (!std::isfinite(f[0]) || !std::isfinite(f[1]) || !std::isfinite(f[2]))
        {
            continue;
        }
        torque += e_z.dot(rel.cross(f));
        f_visc_sum += f_visc[i];
        f_pres_sum += f_pres[i];
    }
    if (snap.tip_speed < TinyReal && snap.r_tip > TinyReal)
    {
        snap.tip_speed = std::abs(omega) * snap.r_tip;
    }
    snap.f_visc = f_visc_sum.norm();
    snap.f_pres = f_pres_sum.norm();
    snap.torque_z = torque;
    return snap;
}

inline void hostMapRotorLoadToProxy(BaseParticles &rotor, BaseParticles &proxy)
{
    syncVariableToHost<Vecd>(rotor, "Position");
    syncVariableToHost<Vecd>(rotor, "ViscousForceFromFluid");
    syncVariableToHost<Vecd>(rotor, "PressureForceFromFluid");
    syncVariableToHost<Vecd>(proxy, "Position");
    const size_t n_rotor = rotor.TotalRealParticles();
    const size_t n_proxy = proxy.TotalRealParticles();
    const Vecd *r_pos = rotor.getVariableDataByName<Vecd>("Position");
    const Vecd *f_visc = rotor.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    const Vecd *f_pres = rotor.getVariableDataByName<Vecd>("PressureForceFromFluid");
    const Vecd *p_pos = proxy.getVariableDataByName<Vecd>("Position");
    Vecd *p_visc = proxy.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    Vecd *p_pres = proxy.getVariableDataByName<Vecd>("PressureForceFromFluid");
    for (size_t j = 0; j < n_proxy; ++j)
    {
        Real best = MaxReal;
        size_t best_i = 0;
        for (size_t i = 0; i < n_rotor; ++i)
        {
            if (!std::isfinite(r_pos[i][0]))
            {
                continue;
            }
            const Real d2 = (p_pos[j] - r_pos[i]).squaredNorm();
            if (d2 < best)
            {
                best = d2;
                best_i = i;
            }
        }
        p_visc[j] = f_visc[best_i];
        p_pres[j] = f_pres[best_i];
    }
    syncVariableToDevice<Vecd>(proxy, "ViscousForceFromFluid");
    syncVariableToDevice<Vecd>(proxy, "PressureForceFromFluid");
}

} // namespace

int main(int ac, char *av[])
{
    const OphelieFrenchStirringPhysicsParams physics;
    OphelieParameters params;
    OphelieFrenchStirringCaseParams stirring;
    applyFrenchReducedDefaults(params, stirring.french);
    applyFrenchReducedStirringDefaults(stirring, physics);

    LocalCli local_cli;
    applyLocalCli(ac, av, local_cli);
    // Resolves the CAD placement and pushes the melt bbox into french.glass_* / coil stack.
    const StdVec<std::string> filtered = filterFrenchStirringCommandLine(ac, av, stirring);
    OphelieFrenchReducedCaseParams &french = stirring.french;

    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_self_induction_ = local_cli.enable_self_induction;
    params.enable_power_scaling_ = false;
    if (params.self_induction_phi_eq_res_tolerance_ > Real(1.0e-3))
    {
        params.self_induction_phi_eq_res_tolerance_ = Real(2.0e-4);
    }
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);

    StdVec<std::string> sph_arguments;
    for (size_t i = 0; i < filtered.size(); ++i)
    {
        if (i == 0 || !isLocalCliOption(filtered[i].c_str()))
        {
            sph_arguments.push_back(filtered[i]);
        }
    }
    StdVec<char *> sph_av;
    for (auto &argument : sph_arguments)
    {
        sph_av.push_back(const_cast<char *>(argument.c_str()));
    }

    if (!reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_stirring_em: Reload.xml required under \"" << local_cli.reload_dir
                  << "\"\nRun test_3d_ophelie_french_stirring_glass_relax --bodies=all first.\n";
        return 1;
    }

    printFrenchStirringCaseSummary(stirring);

    const Real rho0 = physics.rho0_glass;
    const Real mu = physics.mu_glass;
    const Real cp = physics.cp_glass;
    const Real k_th = physics.k_glass;
    const Real beta = local_cli.enable_boussinesq ? physics.beta_glass : Real(0);
    const Real t0 = physics.t_initial;
    const Real target_power = french.target_joule_power;

    //----------------------------------------------------------------------
    // Phase A — EM on the reloaded glass particles, then CIC-deposit Q onto a fixed grid.
    // The EM SPHSystem stays alive for the whole process: tearing down a SYCL SPHSystem and
    // creating another one often hangs the GPU runtime.
    //----------------------------------------------------------------------
    rh200::Rh200ScalarEulerianGrid q_grid;
    Real p_joule = 0.0;
    Real phi_eq = 0.0;
    Real q_grid_scale = 1.0;
    Real q_sample_rel_l2 = 0.0;
    UniquePtr<SPHSystem> em_system;
    UniquePtr<SolidBody> glass_em;
    UniquePtr<Inner<>> glass_em_inner;
    UniquePtr<RegisterOphelieGlassFields> register_glass_em;
    {
        const BoundingBoxd bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
        em_system = makeUnique<SPHSystem>(bounds, french.dp);
        em_system->setReloadParticles(true);
        IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

        // Analytic cylinder shape: the particles come from Reload (they already carry the
        // paddle cavity), and the EM phi boundary treats the melt as a z-cylinder anyway.
        glass_em = makeUnique<SolidBody>(
            *em_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                            "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height));
        glass_em->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        glass_em->defineMatterMaterial<Solid>();
        glass_em->generateParticles<BaseParticles, Reload>(glass_em->Name());
        em_system->initializeSystemCellLinkedLists();
        em_system->initializeSystemConfigurations();

        OphelieGlassFieldNames glass_names;
        register_glass_em = makeUnique<RegisterOphelieGlassFields>(*glass_em, glass_names);
        glass_em_inner = makeUnique<Inner<>>(*glass_em);
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(*glass_em, glass_names,
                                                                                  params.sigma_glass_);
        assign_sigma.exec();

        OphelieFrenchEmJouleHeatOneWayResult probe;
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_names, params, french, probe);
        (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_names, params, french, probe);
        p_joule = probe.joule_power_w;
        phi_eq = probe.phi_eq_res_vol;

        BaseParticles &particles = glass_em->getBaseParticles();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
        syncVariableToHost<Real>(particles, q_field);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *q = particles.getVariableDataByName<Real>(q_field);
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const size_t n_em = particles.TotalRealParticles();

        StdVec<Vecd> em_pos(pos, pos + n_em);
        StdVec<Real> em_q(q, q + n_em);
        StdVec<Real> em_vol(vol, vol + n_em);

        const Real spacing = local_cli.q_grid_spacing > TinyReal ? local_cli.q_grid_spacing : french.dp;
        const Vecd &gc = french.glass_center;
        const BoundingBoxd melt_bounds(
            Vecd(gc[0] - french.glass_radius, gc[1] - french.glass_radius, gc[2] - french.glass_half_height),
            Vecd(gc[0] + french.glass_radius, gc[1] + french.glass_radius, gc[2] + french.glass_half_height));
        q_grid.spec_ = rh200::makeRh200JouleHeatGridSpecFromBounds(melt_bounds, spacing, Real(2));
        q_grid.resetAccumulators();
        q_grid.depositScalarCloudInCell(em_pos, em_q, em_vol);
        q_grid.finalizeFromAccumulators();

        // Interpolation is not conservative, so renormalize the grid to the calibrated power.
        const Real p_sampled = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
        q_grid_scale = target_power / (p_sampled + TinyReal);
        q_grid.scaleField(q_grid_scale);

        rh200::Rh200EmParticleDepositionHost em_host;
        em_host.position = em_pos;
        em_host.volumetric_measure = em_vol;
        em_host.joule_heat = em_q;
        Real rel_max = 0.0;
        rh200::hostCompareSampleToEmField(q_grid, em_host, q_sample_rel_l2, rel_max);

        std::cout << "[ophelie][stirring-em] EM handoff: n=" << n_em << " P_joule_W=" << p_joule
                  << " phi_eq_res_vol=" << phi_eq << " self_induction=" << (probe.used_self_induction ? 1 : 0)
                  << "\n[ophelie][stirring-em] Q grid: " << q_grid.spec_.nx_ << "x" << q_grid.spec_.ny_ << "x"
                  << q_grid.spec_.nz_ << " h=" << q_grid.spec_.spacing_ << " P_sampled_W=" << p_sampled
                  << " scale=" << q_grid_scale << " rel_l2=" << q_sample_rel_l2 << " rel_max=" << rel_max << std::endl;
    }

    //----------------------------------------------------------------------
    // Phase B — WCSPH melt + Simbody paddle + Boussinesq + French thermal BC.
    //----------------------------------------------------------------------
    const BoundingBoxd flow_bounds = frenchStirringFlowDomainBounds(stirring);
    SPHSystem sph_system(flow_bounds, french.dp);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

    const Real u_ref = frenchStirringReferenceSpeed(stirring);
    const Real wall_thickness = frenchStirringWallThickness(stirring);

    SolidBody rotor(sph_system, makeShared<OphelieFrenchStirringRotorShape>("Rotor", stirring));
    rotor.defineAdaptationRatios(1.3, 1.0);
    rotor.defineMatterMaterial<Solid>(stirring.rho_rotor);
    rotor.generateParticles<BaseParticles, Reload>(rotor.Name());

    SolidBody wall(sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                                   "WallBoundary", french, wall_thickness, stirring.wall_mesh_resolution));
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    (void)wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall.generateParticles<BaseParticles, Reload>(wall.Name());

    FluidBody glass(sph_system, makeShared<OphelieFrenchStirringGlassFluidShape>("GlassBody", stirring));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, local_cli.c0);
    glass.addMaterialProperty<Viscosity>(mu);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());

    // Report the reload contents before anything consumes them. A truncated relax (or one run
    // without --bodies=all) otherwise surfaces much later as a NaN Simbody inertia matrix, because
    // SolidBodyPartForSimbody::setMassProperties divides by a zero body-part volume.
    const size_t n_rotor = rotor.getBaseParticles().TotalRealParticles();
    const size_t n_wall = wall.getBaseParticles().TotalRealParticles();
    const size_t n_glass = glass.getBaseParticles().TotalRealParticles();
    std::cout << "[ophelie][stirring-em] reload: n_glass=" << n_glass << " n_rotor=" << n_rotor
              << " n_wall=" << n_wall << std::endl;
    if (n_rotor == 0 || n_wall == 0 || n_glass == 0)
    {
        std::cerr << "[ophelie][stirring-em] empty body in \"" << local_cli.reload_dir
                  << "/Reload.xml\". Rerun test_3d_ophelie_french_stirring_glass_relax --bodies=all "
                     "and let it finish.\n";
        return 1;
    }

    TriangleMeshShapeSTL rotor_surface(stirring.rotor_stl_path, stirring.rotor_translation, stirring.geometry_scale);
    ObserverBody rotor_proxy(sph_system, "RotorProxy");
    rotor_proxy.defineAdaptationRatios(2.0);
    rotor_proxy.generateParticles<ObserverParticles>(rotor_surface);

    Inner<> glass_inner(glass);
    Contact<> glass_solid_contact(glass, {&wall, &rotor});
    Contact<> rotor_fluid_contact(rotor, {&glass});

    BaseParticles &glass_particles = glass.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(glass_particles, t0);
    registerOphelieThermalDiffusionAuxFields(glass_particles, k_th);
    glass_particles.registerStateVariable<Real>(rh200::kJouleHeatField, Real(0));
    (void)glass_particles.getVariableByName<Real>(rh200::kJouleHeatField)->DelegatedData(MainExecutionPolicy{});

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = true;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = false;
    thermal_bc.boundary_width_factor = params.phi_boundary_distance_factor_;
    thermal_bc.h_side = local_cli.h_side > Real(0) ? local_cli.h_side : physics.h_side;
    thermal_bc.h_bottom = local_cli.h_bottom > Real(0) ? local_cli.h_bottom : physics.h_bottom;
    thermal_bc.h_free = local_cli.h_free > Real(0) ? local_cli.h_free : physics.h_free;
    thermal_bc.emissivity = physics.emissivity;
    thermal_bc.t_cool = physics.t_cool;
    thermal_bc.t_ambient = physics.t_ambient;
    thermal_bc.t_rad_ambient = physics.t_rad_ambient;
    // Thermal tagging follows the settling melt, so it needs its own copy of the melt extent;
    // the EM handoff above is already done and must keep the nominal CAD cylinder.
    OphelieFrenchReducedCaseParams french_thermal = french;
    const Real melt_floor_z = french.glass_center[2] - french.glass_half_height;
    refreshMeltExtentFromParticles(glass_particles, french_thermal, melt_floor_z);
    const size_t n_bc_particles =
        setupOphelieThermalFrenchNaturalBoundaryFaces(glass_particles, thermal_bc, french_thermal);

    // The literature coefficients extract ~145 kW at 1473 K against 60 kW of induction, because the
    // real cold crucible grows an insulating solidified skull on the walls and a crust on the free
    // surface, neither of which this model carries. --balance-heat-loss scales every boundary
    // coefficient by one common factor so the melt starts at the paper's steady 1473 K operating
    // point while keeping the literature's relative split between the faces.
    Real bc_scale = 1.0;
    {
        Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                      side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
        std::cout << "[ophelie][stirring-em] thermal BC at T0: side=" << side_w << " W bottom=" << bottom_w
                  << " W free_conv=" << free_conv_w << " W free_rad=" << free_rad_w << " W total=" << total_loss_w
                  << " W vs P_joule=" << target_power << " W" << std::endl;
        if (local_cli.balance_heat_loss && total_loss_w > TinyReal)
        {
            bc_scale = target_power / total_loss_w;
            thermal_bc.h_side *= bc_scale;
            thermal_bc.h_bottom *= bc_scale;
            thermal_bc.h_free *= bc_scale;
            thermal_bc.emissivity *= bc_scale;
            std::cout << "[ophelie][stirring-em] balanced boundary coefficients by " << bc_scale
                      << ": h_side=" << thermal_bc.h_side << " h_bottom=" << thermal_bc.h_bottom
                      << " h_free=" << thermal_bc.h_free << " emissivity=" << thermal_bc.emissivity << std::endl;
        }
    }

    rh200::Rh200JouleHeatGridDevice q_grid_device;
    q_grid_device.upload(q_grid);
    rh200::writeRh200ScalarGridVtiAscii("./output/StirringJouleHeatGrid.vti", q_grid.spec_,
                                        {{"JouleHeat_Wpm3", &q_grid}});

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall).exec();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(rotor).exec();

    auto &update_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(glass);
    auto &update_wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall);
    auto &update_rotor_cell_linked_list = main_methods.addCellLinkedListDynamics(rotor);
    auto &update_glass_relations = main_methods.addRelationDynamics(glass_inner, glass_solid_contact);
    auto &update_rotor_fluid_contact = main_methods.addRelationDynamics(rotor_fluid_contact);
    auto &fluid_particle_sorting = main_methods.addSortDynamics(glass);

    Gravity gravity(Vecd(0.0, 0.0, -physics.gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(glass, gravity);
    auto &boussinesq_force =
        main_methods.addStateDynamics<BoussinesqBuoyancyForceCK>(glass, gravity, beta, t0, kOphelieTemperatureField);
    auto &fluid_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(glass);
    auto &fluid_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(glass);

    auto &acoustic_step_1st_half =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(glass_solid_contact);
    auto &acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK,
                                                    NoKernelCorrectionCK>(glass_inner);
    auto &acoustic_step_2nd_half_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_solid_contact);
    acoustic_step_2nd_half.addPostContactInteraction(acoustic_step_2nd_half_wall);
    auto &density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
            .addPostContactInteraction(glass_solid_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(glass);

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(
            glass_inner);
    auto &fluid_viscous_force_from_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(
            glass_solid_contact);
    fluid_viscous_force.addPostContactInteraction(fluid_viscous_force_from_wall);
    auto &viscous_force_from_fluid_on_rotor =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::ViscousForceFromFluid, std::remove_reference_t<decltype(fluid_viscous_force_from_wall)>>(
            rotor_fluid_contact);
    auto &pressure_force_from_fluid_on_rotor =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::PressureForceFromFluid, std::remove_reference_t<decltype(acoustic_step_2nd_half_wall)>>(
            rotor_fluid_contact);

    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(glass, u_ref);
    auto &acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(glass, 0.55);
    auto &u_max_reduce = main_methods.addReduceDynamics<OphelieVelocityMaxReduceCK>(glass);
    auto &t_mean_reduce = main_methods.addReduceDynamics<OphelieTemperatureMeanReduceCK>(glass);
    auto &t_max_reduce =
        main_methods.addReduceDynamics<OphelieThermalMaxTemperatureReduceCK>(glass, kOphelieTemperatureField);

    StateDynamics<MainExecutionPolicy, rh200::SampleJouleHeatFromGridCK> sample_joule_heat(glass, q_grid_device);
    auto &out_of_grid_probe = main_methods.addReduceDynamics<rh200::OutOfGridJouleSampleReduceCK>(glass, q_grid_device);
    StateDynamics<MainExecutionPolicy, ApplyGridJouleHeatToTemperatureCK> apply_q(glass, rho0, cp, t0);
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc);
    StateDynamics<MainExecutionPolicy, ClampTemperatureCK> clamp_temperature(glass, local_cli.t_min, local_cli.t_max,
                                                                             t0);

    auto &glass_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass);
    glass_state_recorder.addToWrite<Vecd>(glass, "Velocity");
    glass_state_recorder.addToWrite<Real>(glass, "Pressure");
    glass_state_recorder.addToWrite<Real>(glass, kOphelieTemperatureField);
    glass_state_recorder.addToWrite<Real>(glass, rh200::kJouleHeatField);

    auto &rotor_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(rotor);
    rotor_state_recorder.addToWrite<Vecd>(rotor, "Velocity");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "ViscousForceFromFluid");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "PressureForceFromFluid");

    rotor.getBaseParticles().registerStateVariable<Vecd>("Velocity");
    rotor_proxy.getBaseParticles().registerStateVariable<Vecd>("Velocity");
    rotor_proxy.getBaseParticles().registerStateVariable<Vecd>("ViscousForceFromFluid");
    rotor_proxy.getBaseParticles().registerStateVariable<Vecd>("PressureForceFromFluid");
    BodyStatesRecordingToTriangleMeshVtpCK<MainExecutionPolicy> write_rotor_surface(rotor_proxy, rotor_surface);
    write_rotor_surface.addToWrite<Vecd>(rotor_proxy, "Velocity");
    write_rotor_surface.addToWrite<Vecd>(rotor_proxy, "ViscousForceFromFluid");
    write_rotor_surface.addToWrite<Vecd>(rotor_proxy, "PressureForceFromFluid");

    // Snapshot the reload/STL rest pose on the host *before* any device kernel touches
    // Position. The GPU spin path wrote nan into Rotor_*.vtp because device Position
    // never came back as finite coordinates.
    const HostSpinPose rotor_rest =
        captureHostSpinPose(rotor.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis, true);
    const HostSpinPose proxy_rest =
        captureHostSpinPose(rotor_proxy.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis, false);
    std::cout << "[ophelie][stirring-em] rotor rest pose: n=" << rotor_rest.pos0.size()
              << " n_nan=" << rotor_rest.n_nonfinite << " r_tip=" << rotor_rest.r_tip
              << " m; proxy n=" << proxy_rest.pos0.size() << " n_nan=" << proxy_rest.n_nonfinite
              << std::endl;
    if (rotor_rest.n_nonfinite > 0 || rotor_rest.r_tip < TinyReal)
    {
        std::cerr << "[ophelie][stirring-em] Rotor reload positions are non-finite or collapsed. "
                     "Rerun test_3d_ophelie_french_stirring_glass_relax --bodies=all.\n";
        return 1;
    }

    auto apply_rotor_kinematics = [&](Real physical_time) {
        applyHostSpinPose(rotor.getBaseParticles(), rotor_rest, stirring.rotation_center, stirring.rotation_axis,
                          stirring.rotation_speed_rad_s, physical_time, true);
        applyHostSpinPose(rotor_proxy.getBaseParticles(), proxy_rest, stirring.rotation_center,
                          stirring.rotation_axis, stirring.rotation_speed_rad_s, physical_time, false);
        rotor.setNewlyUpdated();
        rotor_proxy.setNewlyUpdated();
    };

    auto write_rotor_visualization = [&](size_t ite, Real physical_time) {
        apply_rotor_kinematics(physical_time);
        hostMapRotorLoadToProxy(rotor.getBaseParticles(), rotor_proxy.getBaseParticles());
        rotor.setNewlyUpdated();
        rotor_proxy.setNewlyUpdated();
        rotor_state_recorder.writeToFile(ite);
        write_rotor_surface.writeToFile(ite);
    };

    apply_rotor_kinematics(Real(0));

    update_glass_cell_linked_list.exec();
    update_wall_cell_linked_list.exec();
    update_rotor_cell_linked_list.exec();
    update_glass_relations.exec();
    update_rotor_fluid_contact.exec();
    density_regularization.exec();
    fluid_advection_step_setup.exec();
    sample_joule_heat.exec();
    constant_gravity.exec();
    boussinesq_force.exec();

    const Real n_out_of_grid = out_of_grid_probe.exec();

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording = time_stepper.addTriggerByInterval(std::max(local_cli.state_record_interval, Real(1.0e-6)));

    fs::create_directories("./output");
    std::ofstream monitor("./output/french_stirring_em_monitor.csv");
    monitor << "advection_step,physical_time_s,revolutions,acoustic_dt_s,u_max_mps,t_mean_K,t_max_K,"
               "rotor_tip_mps,rotor_r_tip_m,f_visc_N,f_pres_N,torque_z_Nm,wall_clock_s\n";
    monitor << std::setprecision(10);

    if (local_cli.state_recording)
    {
        glass_state_recorder.writeToFile(0);
        write_rotor_visualization(0, Real(0));
    }

    const Real revolution_time = Real(2.0 * M_PI) / std::max(stirring.rotation_speed_rad_s, TinyReal);
    const RotorLoadSnapshot rotor0 =
        hostRotorLoadSnapshot(rotor.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis,
                              stirring.rotation_speed_rad_s);
    const Real omega = stirring.rotation_speed_rad_s;
    const Real rpm = omega * Real(60.0) / Real(2.0 * M_PI);
    const Real u_tip_expected = omega * rotor0.r_tip;
    const Real re_impeller = rho0 * (omega / Real(2.0 * M_PI)) * (Real(2) * rotor0.r_tip) * (Real(2) * rotor0.r_tip) / mu;
    std::cout << "[ophelie][stirring-em] flow start: end_time=" << local_cli.end_time << " s ("
              << local_cli.end_time / revolution_time << " revolutions) budget=" << local_cli.max_wall_hours << " h\n"
              << "  n_glass=" << glass_particles.TotalRealParticles()
              << " n_rotor=" << rotor.getBaseParticles().TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " n_bc_shell=" << n_bc_particles
              << " out_of_grid=" << n_out_of_grid << "\n"
              << "  rho=" << rho0 << " mu=" << mu << " cp=" << cp << " k=" << k_th << " beta=" << beta
              << " c0=" << local_cli.c0 << " U_ref=" << u_ref << " T0=" << t0 << " T_floor=" << local_cli.t_min
              << "\n  paddle: rpm=" << rpm << " omega=" << omega << " rad/s r_tip=" << rotor0.r_tip
              << " m U_tip=" << u_tip_expected << " m/s Re_imp=" << re_impeller
              << " (creeping: swirl is weak, look at GlassBody Velocity glyphs)\n"
              << std::endl;

    const auto wall_clock_start = std::chrono::steady_clock::now();
    const double wall_clock_budget_s = static_cast<double>(local_cli.max_wall_hours) * 3600.0;
    auto elapsed_s = [&wall_clock_start]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_clock_start).count();
    };

    size_t advection_steps = 0;
    size_t acoustic_steps = 0;
    bool budget_exhausted = false;
    bool diverged = false;
    Real acoustic_dt = 0.0;
    while (!time_stepper.isEndTime(local_cli.end_time))
    {
        acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        ++acoustic_steps;
        // A non-finite or collapsed step means the flow blew up; without this the loop would
        // spin forever because comparisons against a NaN physical time are always false.
        if (!std::isfinite(acoustic_dt) || acoustic_dt < Real(1.0e-9) ||
            !std::isfinite(time_stepper.getPhysicalTime()))
        {
            std::cerr << "[ophelie][stirring-em] diverged at acoustic step " << acoustic_steps
                      << " t=" << time_stepper.getPhysicalTime() << " dt=" << acoustic_dt
                      << " U_max=" << u_max_reduce.exec() << std::endl;
            diverged = true;
            break;
        }
        acoustic_step_1st_half.exec(acoustic_dt);
        pressure_force_from_fluid_on_rotor.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);

        apply_q.exec(acoustic_dt);
        apply_natural_bc.exec(acoustic_dt);
        clamp_temperature.exec(acoustic_dt);
        boussinesq_force.exec();

        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            apply_rotor_kinematics(time_stepper.getPhysicalTime());

            if (local_cli.state_recording && state_recording())
            {
                glass_state_recorder.writeToFile(advection_steps);
                write_rotor_visualization(advection_steps, time_stepper.getPhysicalTime());
            }

            if (advection_steps % 50 == 0)
            {
                fluid_particle_sorting.exec();
            }
            update_glass_cell_linked_list.exec();
            update_rotor_cell_linked_list.exec();
            update_glass_relations.exec();
            update_rotor_fluid_contact.exec();
            density_regularization.exec();
            fluid_advection_step_setup.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
            boussinesq_force.exec();
            viscous_force_from_fluid_on_rotor.exec();

            // The melt has moved through the laboratory-fixed inductor field.
            sample_joule_heat.exec();
            if (local_cli.bc_retag_every > 0 && advection_steps % static_cast<size_t>(local_cli.bc_retag_every) == 0)
            {
                refreshMeltExtentFromParticles(glass_particles, french_thermal, melt_floor_z);
                (void)setupOphelieThermalFrenchNaturalBoundaryFaces(glass_particles, thermal_bc, french_thermal);
            }

            const bool write_csv =
                local_cli.csv_every > 0 && advection_steps % static_cast<size_t>(local_cli.csv_every) == 0;
            const bool write_screen =
                local_cli.screen_every > 0 && advection_steps % static_cast<size_t>(local_cli.screen_every) == 0;
            if (write_csv || write_screen)
            {
                const Real t_now = time_stepper.getPhysicalTime();
                const Real u_max_now = u_max_reduce.exec();
                const Real t_mean_now = t_mean_reduce.exec();
                const Real t_max_now = t_max_reduce.exec();
                const RotorLoadSnapshot rotor_load = hostRotorLoadSnapshot(
                    rotor.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis,
                    stirring.rotation_speed_rad_s);
                if (write_csv)
                {
                    monitor << advection_steps << "," << t_now << "," << t_now / revolution_time << "," << acoustic_dt
                            << "," << u_max_now << "," << t_mean_now << "," << t_max_now << "," << rotor_load.tip_speed
                            << "," << rotor_load.r_tip << "," << rotor_load.f_visc << "," << rotor_load.f_pres
                            << "," << rotor_load.torque_z << "," << elapsed_s() << "\n";
                    monitor.flush();
                }
                if (write_screen)
                {
                    std::cout << std::fixed << std::setprecision(6) << "[ophelie][stirring-em] N=" << advection_steps
                              << " n_ac=" << acoustic_steps << " t=" << t_now << " dt=" << std::scientific
                              << acoustic_dt << std::fixed << " rev=" << t_now / revolution_time
                              << " U_max=" << u_max_now << " U_tip=" << rotor_load.tip_speed
                              << " Fv=" << rotor_load.f_visc << " Fp=" << rotor_load.f_pres
                              << " Tz=" << rotor_load.torque_z << " T_mean=" << t_mean_now
                              << " T_max=" << t_max_now << " wall=" << elapsed_s() / 3600.0 << " h" << std::endl;
                }
            }

            if (wall_clock_budget_s > 0.0 && elapsed_s() > wall_clock_budget_s)
            {
                budget_exhausted = true;
                break;
            }
        }
    }

    if (local_cli.state_recording)
    {
        glass_state_recorder.writeToFile(advection_steps + 1);
        write_rotor_visualization(advection_steps + 1, time_stepper.getPhysicalTime());
    }

    const Real u_max = u_max_reduce.exec();
    const Real t_mean = t_mean_reduce.exec();
    const Real t_max = t_max_reduce.exec();
    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                  side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
    monitor.close();

    const Real power_rel_err = std::abs(p_joule - target_power) / (target_power + TinyReal);
    const bool em_ok = std::isfinite(p_joule) && power_rel_err < Real(1.0e-2);
    const bool grid_ok = std::isfinite(q_sample_rel_l2) && q_sample_rel_l2 < Real(0.2) && n_out_of_grid < Real(1);
    const bool flow_ok = !diverged && std::isfinite(u_max) && u_max > TinyReal && u_max <= local_cli.c0;
    const bool thermal_ok = std::isfinite(t_mean) && std::isfinite(t_max) && total_loss_w > TinyReal;
    const bool passed = em_ok && grid_ok && flow_ok && thermal_ok;

    std::cout << "test_3d_ophelie_french_stirring_em"
              << " P_joule_W=" << p_joule << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << phi_eq
              << " q_grid_rel_l2=" << q_sample_rel_l2 << " out_of_grid=" << n_out_of_grid
              << " physical_time_s=" << time_stepper.getPhysicalTime()
              << " revolutions=" << time_stepper.getPhysicalTime() / revolution_time
              << " advection_steps=" << advection_steps << " U_max=" << u_max << " T_mean=" << t_mean
              << " T_max=" << t_max << " wall_loss_side_W=" << side_w << " wall_loss_bottom_W=" << bottom_w
              << " free_conv_loss_W=" << free_conv_w << " free_rad_loss_W=" << free_rad_w
              << " total_heat_loss_W=" << total_loss_w << " bc_scale=" << bc_scale
              << " h_side_used=" << thermal_bc.h_side << " wall_clock_h=" << elapsed_s() / 3600.0
              << " acoustic_steps=" << acoustic_steps << " budget_exhausted=" << (budget_exhausted ? 1 : 0)
              << " diverged=" << (diverged ? 1 : 0) << " em_ok=" << (em_ok ? 1 : 0)
              << " grid_ok=" << (grid_ok ? 1 : 0) << " flow_ok=" << (flow_ok ? 1 : 0)
              << " thermal_ok=" << (thermal_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
