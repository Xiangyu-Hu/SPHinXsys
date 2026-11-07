/**
 * @file test_3d_ophelie_french_stirring_em.cpp
 * @brief French cold-crucible glass melter: induction heating + mechanical stirring.
 *
 * Phase A solves the EM problem on the reloaded glass particles and deposits Joule
 * power density onto a laboratory-fixed Eulerian grid. Phase B runs WCSPH with a
 * host-spun paddle, Boussinesq buoyancy and the French Robin/radiation thermal BC,
 * resampling Q from that grid every advection step.
 *
 * Default is one-shot EM. `--thermo-em-coupling=periodic` reuses the natural-convection
 * controller: azimuthal σ(T) → EM → Euler Q, on physical time. Periodic mode also
 * turns SPH conduction on unless `--no-thermal-diffusion`.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_stirring_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_french_skull_resistance.h"
#include "electromagnetic_ophelie_french_skull_shell.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "electromagnetic_ophelie_french_skin_depth_diagnostic.h"
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
    /** Doc 40 R03: scale SPH face weights to geometric A (default on for literature). */
    bool boundary_area_normalize = true;
    int bc_retag_every = 50;    // advection steps between Robin/radiation face re-tagging
    int screen_every = 20;
    int csv_every = 20;
    bool state_recording = true;
    Real state_record_interval = 10.0;
    bool enable_self_induction = true;
    bool enable_boussinesq = true;
    bool q_conservative_remap = false;
    bool thermal_diffusion_explicit = false;
    /** If true, process exit code requires validation_passed (doc 23), not only smoke. */
    bool require_validation = false;
    /** Doc 69 P1: write output/french_f07_state_audit.csv with same-state Mass/ρ/V snapshots. */
    bool f07_state_audit = false;
    FrenchLiteratureParameters literature = makeFrenchStirredLiteratureParameters();
    FrenchSkullResistanceParams skull;
    FrenchFreeCrustParams free_crust;
    FrenchSkullShellParams skull_shell;
};

inline bool isLocalCliOption(const char *arg)
{
    return std::strncmp(arg, "--reload-dir=", 13) == 0 || std::strncmp(arg, "--end-time=", 11) == 0 ||
           std::strncmp(arg, "--max-wall-hours=", 17) == 0 || std::strncmp(arg, "--c0=", 5) == 0 ||
           std::strncmp(arg, "--q-grid-spacing=", 17) == 0 || std::strncmp(arg, "--t-min=", 8) == 0 ||
           std::strncmp(arg, "--t-max=", 8) == 0 || std::strncmp(arg, "--h-side=", 9) == 0 ||
           std::strncmp(arg, "--h-bottom=", 11) == 0 || std::strncmp(arg, "--h-free=", 9) == 0 ||
           std::strcmp(arg, "--balance-heat-loss") == 0 || std::strcmp(arg, "--boundary-area-normalize") == 0 ||
           std::strcmp(arg, "--no-boundary-area-normalize") == 0 ||
           std::strncmp(arg, "--bc-retag-every=", 17) == 0 ||
           std::strncmp(arg, "--screen-every=", 15) == 0 || std::strncmp(arg, "--csv-every=", 12) == 0 ||
           std::strncmp(arg, "--state-record-interval=", 24) == 0 ||
           std::strcmp(arg, "--no-state-recording") == 0 || std::strcmp(arg, "--no-self-induction") == 0 ||
           std::strcmp(arg, "--no-boussinesq") == 0 || std::strcmp(arg, "--q-conservative-remap") == 0 ||
           std::strcmp(arg, "--require-validation") == 0 || std::strcmp(arg, "--f07-state-audit") == 0 ||
           isFrenchLiteratureCouplingCliOption(arg) ||
           isFrenchSkullCliOption(arg) || isFrenchSkullShellCliOption(arg);
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
        else if (std::strcmp(av[i], "--boundary-area-normalize") == 0)
        {
            cli.boundary_area_normalize = true;
        }
        else if (std::strcmp(av[i], "--no-boundary-area-normalize") == 0)
        {
            cli.boundary_area_normalize = false;
        }
        else if (tryApplyFrenchSkullCli(av[i], cli.skull))
        {
        }
        else if (tryApplyFrenchFreeCrustCli(av[i], cli.free_crust))
        {
        }
        else if (tryApplyFrenchSkullShellCli(av[i], cli.skull_shell))
        {
        }
        else if (std::strncmp(av[i], "--bc-retag-every=", 17) == 0)
        {
            cli.bc_retag_every = std::atoi(av[i] + 17);
        }
        else if (std::strncmp(av[i], "--screen-every=", 15) == 0)
        {
            cli.screen_every = std::atoi(av[i] + 15);
        }
        else if (std::strncmp(av[i], "--csv-every=", 12) == 0)
        {
            // Doc 47 F02 genuine cadence: must vary with screen (or set 0 to disable CSV).
            // Shared energy/host block runs under write_csv || write_screen.
            cli.csv_every = std::atoi(av[i] + 12);
        }
        else if (std::strncmp(av[i], "--state-record-interval=", 24) == 0)
        {
            cli.state_record_interval = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strcmp(av[i], "--no-state-recording") == 0)
        {
            cli.state_recording = false;
        }
        else if (std::strcmp(av[i], "--require-validation") == 0)
        {
            cli.require_validation = true;
        }
        else if (std::strcmp(av[i], "--f07-state-audit") == 0)
        {
            cli.f07_state_audit = true;
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.enable_self_induction = false;
            cli.literature.enable_aind = false;
        }
        else if (tryApplyFrenchLiteratureCli(av[i], cli.literature, cli.thermal_diffusion_explicit))
        {
            if (std::strncmp(av[i], "--aind=", 7) == 0)
            {
                cli.enable_self_induction = cli.literature.enable_aind;
            }
            else if (std::strcmp(av[i], "--q-conservative-remap") == 0)
            {
                cli.q_conservative_remap = true;
            }
            else if (std::strcmp(av[i], "--no-boussinesq") == 0)
            {
                cli.enable_boussinesq = false;
            }
        }
        else if (std::strcmp(av[i], "--no-boussinesq") == 0)
        {
            cli.enable_boussinesq = false;
            cli.literature.enable_boussinesq = false;
        }
        else if (std::strcmp(av[i], "--q-conservative-remap") == 0)
        {
            cli.q_conservative_remap = true;
            cli.literature.q_conservative_remap = true;
        }
    }
    if (cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultReloadFolder();
    }
    finalizeFrenchPeriodicLiteratureDefaults(cli.literature, cli.thermal_diffusion_explicit);
    if (cli.literature.coupling == ThermoEMCouplingMode::Periodic)
    {
        cli.enable_self_induction = cli.literature.enable_aind;
    }
    cli.q_conservative_remap = cli.q_conservative_remap || cli.literature.q_conservative_remap;
    cli.enable_boussinesq = cli.enable_boussinesq && cli.literature.enable_boussinesq;
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

/** ViscousForceCK caches μ in InteractKernel; rebuild after mean-μ(T) refresh (Batch D / P1-01). */
inline void invalidateViscosityDependentCkKernels(Inner<> &glass_inner, Contact<> &glass_solid_contact,
                                                  Contact<> *rotor_fluid_contact = nullptr)
{
    glass_inner.resetComputingKernelUpdated();
    const StdVec<SPHBody *> glass_contacts = glass_solid_contact.getContactBodies();
    for (size_t k = 0; k < glass_contacts.size(); ++k)
    {
        glass_solid_contact.resetComputingKernelUpdated(k);
    }
    if (rotor_fluid_contact != nullptr)
    {
        const StdVec<SPHBody *> rotor_contacts = rotor_fluid_contact->getContactBodies();
        for (size_t k = 0; k < rotor_contacts.size(); ++k)
        {
            rotor_fluid_contact->resetComputingKernelUpdated(k);
        }
    }
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
            // Doc 45 F01: ΔT is authoritative; derive absolute T after the increment.
            const Real heating = q_[index_i] * dt * inv_rho_cp_;
            delta_t_[index_i] += heating;
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
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
            // Doc 45 F01: inactive clamp must not rebuild ΔT from rounded absolute T.
            Real t = t_initial_ + delta_t_[index_i];
            if (t < t_min_)
            {
                t = t_min_;
                delta_t_[index_i] = t - t_initial_;
            }
            else if (t > t_max_)
            {
                t = t_max_;
                delta_t_[index_i] = t - t_initial_;
            }
            temperature_[index_i] = t;
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

inline Real hostThermalEnergy(BaseParticles &particles, Real rho, Real cp, Real t_ref)
{
    return hostThermalEnergyConstantCp(particles, rho, cp, t_ref, kOphelieTemperatureField);
}

inline bool depositEmQAndRecord(BaseParticles &em_particles, const OphelieGlassFieldNames &names,
                                const OphelieParameters &params, const FrenchCylindricalFrame &frame,
                                const FrenchLiteratureParameters &lit, rh200::Rh200ScalarEulerianGrid &q_grid,
                                ThermoEMCouplingState &state, Real p_em_recon)
{
    const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
    syncVariableToHost<Real>(em_particles, q_field);
    syncVariableToHost<Vecd>(em_particles, "Position");
    syncVariableToHost<Real>(em_particles, "VolumetricMeasure");
    const Real *q = em_particles.getVariableDataByName<Real>(q_field);
    const Vecd *pos = em_particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = em_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = em_particles.TotalRealParticles();

    state.q_rz.frame = frame;
    azimuthalDeposit(state.q_rz, pos, q, vol, n);
    fillEmptyAxisymmetricBins(state.q_rz, Real(0));

    StdVec<Vecd> em_pos(pos, pos + n);
    StdVec<Real> em_q(q, q + n);
    StdVec<Real> em_vol(vol, vol + n);
    q_grid.resetAccumulators();
    q_grid.depositScalarCloudInCell(em_pos, em_q, em_vol);
    q_grid.finalizeFromAccumulators();

    FrenchPowerMapBook book;
    book.schema_version = 2;
    book.p_em_recon = p_em_recon;
    book.p_axisymmetric_grid = integrateAxisymmetricPower(state.q_rz);
    book.p_euler_sample = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
    book.p_em_particles = hostParticlePower(q, vol, n);
    book.p_sph_particles = book.p_em_particles; // historical alias — EM body, not fluid
    book.p_euler_before_remap = book.p_euler_sample;
    std::cout << "[ophelie][q-map] deposit_path=3D_EM_q_direct_to_euler"
                 " P_axi_diagnostic=" << book.p_axisymmetric_grid << std::endl;
    if (lit.q_conservative_remap && book.p_euler_sample > TinyReal && p_em_recon > TinyReal)
    {
        book.remap_scale = p_em_recon / book.p_euler_sample;
        q_grid.scaleField(book.remap_scale);
        book.conservative_remap_applied = true;
        book.p_euler_sample = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
        book.p_euler_after_remap = book.p_euler_sample;
        std::cout << "[ophelie][q-map] conservative remapping (not EM calibration): P_before="
                  << book.p_euler_before_remap << " scale=" << book.remap_scale << " P_after="
                  << book.p_euler_after_remap << " domain=EM_particles" << std::endl;
    }
    logFrenchPowerMapBook(book);
    state.last_power = book;
    state.reconstructed_glass_power_w = p_em_recon;
    return true;
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
    stirring.french.target_joule_power = local_cli.literature.target_glass_absorbed_power_w;
    // Resolves the CAD placement and pushes the melt bbox into french.glass_* / coil stack.
    const StdVec<std::string> filtered = filterFrenchStirringCommandLine(ac, av, stirring);
    OphelieFrenchReducedCaseParams &french = stirring.french;

    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    {
        OphelieTestCliOptions em_cli;
        em_cli.ophelie_current_form_user_set = true;
        finalizeOphelieCurrentFormConfiguration(params, em_cli);
    }
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
    logFrenchLiteratureParameters(local_cli.literature);
    std::cout << "[ophelie][stirring-em] cadence csv_every=" << local_cli.csv_every
              << " screen_every=" << local_cli.screen_every
              << " (energy/host block: write_csv||write_screen; use --csv-every= with --screen-every= for F02)"
              << std::endl;

    const OphelieFrenchGlassMaterialLaws &laws = local_cli.literature.material;
    const Real rho0 = physics.rho0_glass; // OPEN_PARAMETER vs Table-1 2750
    const Real t0 = physics.t_initial;
    const Real mu = laws.mu.evaluate(t0);
    const Real cp = laws.cp.evaluate(t0);
    const Real k_th = laws.k.evaluate(t0);
    const Real beta = local_cli.enable_boussinesq ? physics.beta_glass : Real(0);
    const Real target_power = french.target_joule_power;
    const OphelieMaterialLawCheckResult material_check = verifyFrenchMaterialLawsAtReferenceTemperatures(laws);
    std::cout << "[ophelie][stirring-em] material_check=" << material_check.message
              << " mu_1473=" << material_check.mu_1473 << " cp_1473=" << material_check.cp_1473
              << " k_1473=" << material_check.k_1473 << " sigma_1473=" << material_check.sigma_1473
              << " passed=" << (material_check.passed ? 1 : 0) << std::endl;
    {
        std::string enthalpy_msg;
        const bool enthalpy_ok = verifyOphelieEnthalpyEnergyAccounting(enthalpy_msg);
        std::cout << "[ophelie][stirring-em] enthalpy_check=" << enthalpy_msg
                  << " passed=" << (enthalpy_ok ? 1 : 0) << std::endl;
        if (!enthalpy_ok)
        {
            std::cerr << "[ophelie][stirring-em] enthalpy accounting regression failed" << std::endl;
            return 1;
        }
        std::string applied_msg;
        const bool applied_ok = verifyFrenchAppliedEnergyResidualIdentity(applied_msg);
        std::cout << "[ophelie][stirring-em] applied_energy_check=" << applied_msg
                  << " passed=" << (applied_ok ? 1 : 0) << std::endl;
        if (!applied_ok)
        {
            std::cerr << "[ophelie][stirring-em] applied energy residual identity failed" << std::endl;
            return 1;
        }
        std::string gate_msg;
        const bool gate_ok = verifyFrenchCouplingGateNegativeTests(gate_msg);
        std::cout << "[ophelie][stirring-em] gate_negative_check=" << gate_msg
                  << " passed=" << (gate_ok ? 1 : 0) << std::endl;
        if (!gate_ok)
        {
            std::cerr << "[ophelie][stirring-em] coupling gate negative tests failed" << std::endl;
            return 1;
        }
    }
    {
        std::string stage1_msg;
        const bool stage1_ok = verifyFrenchThermalStage1HostAlgebra(stage1_msg);
        std::cout << "[ophelie][stirring-em] thermal_stage1_host=" << stage1_msg
                  << " passed=" << (stage1_ok ? 1 : 0) << std::endl;
        if (!stage1_ok)
        {
            std::cerr << "[ophelie][stirring-em] thermal stage-1 host identities failed" << std::endl;
            return 1;
        }
        std::string visc_msg;
        const bool visc_ok = verifyFrenchViscosityFunctorCachesByValue(visc_msg);
        std::cout << "[ophelie][stirring-em] viscosity_ck_cache_check=" << visc_msg
                  << " passed=" << (visc_ok ? 1 : 0) << std::endl;
        if (!visc_ok)
        {
            std::cerr << "[ophelie][stirring-em] viscosity CK cache regression failed" << std::endl;
            return 1;
        }
    }

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
    OphelieGlassFieldNames glass_em_names;
    ThermoEMCouplingState coupling_state;
    FrenchCylindricalFrame cyl_frame = makeFrenchCylinderFrameFromGlass(
        french.glass_center, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
    const bool periodic_coupling = local_cli.literature.coupling == ThermoEMCouplingMode::Periodic;
    fs::create_directories("./output");
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
        glass_em_names = glass_names;
        register_glass_em = makeUnique<RegisterOphelieGlassFields>(*glass_em, glass_em_names);
        glass_em_inner = makeUnique<Inner<>>(*glass_em);
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(*glass_em, glass_em_names,
                                                                                  params.sigma_glass_);
        assign_sigma.exec();

        OphelieFrenchEmJouleHeatOneWayResult probe;
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
        (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
        p_joule = probe.joule_power_w;
        phi_eq = probe.phi_eq_res_vol;

        BaseParticles &particles = glass_em->getBaseParticles();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_em_names, params);
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

        const Real p_sampled = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
        q_grid_scale = 1.0;
        std::cout << "[ophelie][stirring-em] Q map before remap: P_em=" << p_joule << " P_sampled=" << p_sampled
                  << std::endl;
        if (local_cli.q_conservative_remap && p_sampled > TinyReal && p_joule > TinyReal)
        {
            q_grid_scale = p_joule / (p_sampled + TinyReal);
            q_grid.scaleField(q_grid_scale);
            const Real p_after = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
            std::cout << "[ophelie][stirring-em] conservative remapping (not EM calibration): P_before=" << p_sampled
                      << " scale=" << q_grid_scale << " P_after=" << p_after << std::endl;
        }

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

        {
            const FrenchCylindricalFrame skin_frame = makeFrenchCylinderFrameFromGlass(
                gc, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
            const FrenchSkinDepthRadialReport skin = computeFrenchSkinDepthRadialReport(
                pos, q, vol, n_em, skin_frame, french.frequency_hz, params.sigma_glass_, params.mu0_, french.dp);
            logFrenchSkinDepthRadialReport(skin, "skin-depth");
        }

        if (periodic_coupling)
        {
            cyl_frame = makeFrenchCylinderFrameFromGlass(gc, french.glass_radius, french.glass_half_height,
                                                         local_cli.literature.vertical_axis);
            coupling_state.generator_power_w = local_cli.literature.generator_power_w;
            coupling_state.target_glass_absorbed_power_w = local_cli.literature.target_glass_absorbed_power_w;
            coupling_state.coil_current_per_loop = french.coil.current_per_loop;
            coupling_state.last_phi_eq_res = phi_eq;
            coupling_state.sigma_rz.frame = cyl_frame;
            coupling_state.q_rz.frame = cyl_frame;
            azimuthalDeposit(coupling_state.q_rz, pos, q, vol, n_em);
            fillEmptyAxisymmetricBins(coupling_state.q_rz, Real(0));
            FrenchPowerMapBook book;
            book.schema_version = 2;
            book.p_em_recon = p_joule;
            book.p_axisymmetric_grid = integrateAxisymmetricPower(coupling_state.q_rz);
            book.p_euler_sample = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
            book.p_em_particles = hostParticlePower(q, vol, n_em);
            book.p_sph_particles = book.p_em_particles;
            book.remap_scale = q_grid_scale;
            book.conservative_remap_applied = local_cli.q_conservative_remap;
            logFrenchPowerMapBook(book);
            coupling_state.last_power = book;
            coupling_state.reconstructed_glass_power_w = p_joule;
            syncVariableToHost<Real>(particles, glass_em_names.sigma);
            const Real *sigma0 = particles.getVariableDataByName<Real>(glass_em_names.sigma);
            azimuthalDeposit(coupling_state.sigma_rz, pos, sigma0, vol, n_em);
            fillEmptyAxisymmetricBins(coupling_state.sigma_rz, params.sigma_glass_);
            coupling_state.sigma_em.assign(sigma0, sigma0 + n_em);
            sigmaFieldStats(coupling_state.sigma_em, vol, n_em, coupling_state.last_sigma_min,
                            coupling_state.last_sigma_max, coupling_state.last_sigma_mean);
            coupling_state.last_update_time = 0.0;
            armFrenchPeriodicRefreshSchedule(coupling_state, local_cli.literature.em_update_interval_s);
            coupling_state.update_count = 1;
            writeAxisymmetricCsv("output/french_sigma_rz.csv", coupling_state.sigma_rz, "sigma_Sm");
            writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
            writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
            appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                      local_cli.literature.em_control);
            writeFrenchEnergyBudgetCsvHeader("output/french_energy_budget.csv");
            std::cout << "[ophelie][thermo-em] initial coupling book written; mode="
                      << emControlModeName(local_cli.literature.em_control) << std::endl;
        }
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
    Viscosity &glass_viscosity = glass.getMaterialProperty<Viscosity>();

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
    registerOphelieFrenchThermalRheologyFields(glass_particles, mu, cp, k_th);
    glass_particles.registerStateVariable<Real>(rh200::kJouleHeatField, Real(0));
    (void)glass_particles.getVariableByName<Real>(rh200::kJouleHeatField)->DelegatedData(MainExecutionPolicy{});

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = true;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = local_cli.literature.enable_thermal_diffusion;
    thermal_bc.boundary_width_factor = params.phi_boundary_distance_factor_;
    thermal_bc.renormalize_boundary_area_to_geometry = local_cli.boundary_area_normalize;
    thermal_bc.h_side = local_cli.h_side > Real(0) ? local_cli.h_side : physics.h_side;
    thermal_bc.h_bottom = local_cli.h_bottom > Real(0) ? local_cli.h_bottom : physics.h_bottom;
    thermal_bc.h_free = local_cli.h_free > Real(0) ? local_cli.h_free : physics.h_free;
    thermal_bc.emissivity = physics.emissivity;
    thermal_bc.t_cool = physics.t_cool;
    thermal_bc.t_ambient = physics.t_ambient;
    thermal_bc.t_rad_ambient = physics.t_rad_ambient;

    {
        std::string skull_msg;
        if (!verifyFrenchSkullResistanceRegression(skull_msg))
        {
            std::cerr << "[ophelie][skull-resistance] regression failed: " << skull_msg << std::endl;
            return 1;
        }
        std::cout << "[ophelie][skull-resistance] " << skull_msg << std::endl;
    }
    if (local_cli.skull.enable && local_cli.balance_heat_loss)
    {
        std::cerr << "[ophelie][stirring-em] refuse --skull-resistance together with --balance-heat-loss "
                     "(legacy engineering scale). Use skull alone.\n";
        return 1;
    }
    if (local_cli.boundary_area_normalize && local_cli.balance_heat_loss)
    {
        std::cerr << "[ophelie][stirring-em] refuse --boundary-area-normalize together with --balance-heat-loss "
                     "(area normalize keeps h; balance-heat-loss scales h).\n";
        return 1;
    }
    {
        std::string shell_err;
        if (!assertFrenchSkullShellExclusive(local_cli.skull_shell, local_cli.skull, local_cli.balance_heat_loss,
                                             shell_err))
        {
            std::cerr << "[ophelie][stirring-em] " << shell_err << std::endl;
            return 1;
        }
        if (local_cli.skull_shell.enable)
        {
            std::cerr << "[ophelie][stirring-em] --skull-shell not wired for stirring yet "
                         "(natural parked after Round 3C no-go; see docs/ophelie/19–20). "
                         "Use --skull-resistance + --free-crust-resistance.\n";
            return 1;
        }
    }
    {
        std::string skull_err;
        if (!applyFrenchSkullResistanceToThermalBc(thermal_bc, local_cli.skull, skull_err))
        {
            std::cerr << "[ophelie][stirring-em] skull setup failed: " << skull_err << std::endl;
            return 1;
        }
    }
    {
        std::string crust_err;
        if (!applyFrenchFreeCrustToThermalBc(thermal_bc, local_cli.free_crust.enable, local_cli.free_crust.delta_m,
                                            local_cli.free_crust.k_crust, crust_err))
        {
            std::cerr << "[ophelie][stirring-em] free-crust setup failed: " << crust_err << std::endl;
            return 1;
        }
    }
    Real area_wall_wb = Real(0);
    {
        const Real melt_h = Real(2) * french.glass_half_height;
        const FrenchZeroDHeatBudget budget0d =
            estimateFrenchZeroDHeatBudget(french.glass_radius, melt_h, t0, thermal_bc, target_power);
        const char *tag = "bare_walls";
        if (local_cli.skull.enable && local_cli.free_crust.enable)
        {
            tag = "with_skull_and_free_crust";
        }
        else if (local_cli.skull.enable)
        {
            tag = "with_skull";
        }
        else if (local_cli.free_crust.enable)
        {
            tag = "with_free_crust";
        }
        logFrenchZeroDHeatBudget(budget0d, tag);
        area_wall_wb = budget0d.area_side + budget0d.area_bottom;
    }

    // Thermal tagging follows the settling melt, so it needs its own copy of the melt extent;
    // the EM handoff above is already done and must keep the nominal CAD cylinder.
    OphelieFrenchReducedCaseParams french_thermal = french;
    const Real melt_floor_z = french.glass_center[2] - french.glass_half_height;
    refreshMeltExtentFromParticles(glass_particles, french_thermal, melt_floor_z);
    const size_t n_bc_particles =
        setupOphelieThermalFrenchNaturalBoundaryFaces(glass_particles, thermal_bc, french_thermal);

    // Legacy engineering knob only. Preferred French path is --skull-resistance (Level-1 series model).
    Real bc_scale = 1.0;
    {
        Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                      side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
        std::cout << "[ophelie][stirring-em] thermal BC at T0: side=" << side_w << " W bottom=" << bottom_w
                  << " W free_conv=" << free_conv_w << " W free_rad=" << free_rad_w << " W total=" << total_loss_w
                  << " W vs P_joule=" << target_power << " W" << std::endl;
        {
            // P0-06: compare discrete shell areas to CAD and to melt-extent geometry (post-settling).
            const OphelieFrenchBoundaryAreaAudit audit_cad =
                hostOphelieThermalFrenchNaturalBoundaryAreaAudit(glass_particles, thermal_bc, french,
                                                                 thermal_bc.boundary_shell_thickness);
            logOphelieFrenchBoundaryAreaAudit(audit_cad, "vs_CAD_cylinder");
            const OphelieFrenchBoundaryAreaAudit audit_melt =
                hostOphelieThermalFrenchNaturalBoundaryAreaAudit(glass_particles, thermal_bc, french_thermal,
                                                                 thermal_bc.boundary_shell_thickness);
            logOphelieFrenchBoundaryAreaAudit(audit_melt, "vs_melt_extent");
            const FrenchZeroDHeatBudget budget_melt = estimateFrenchZeroDHeatBudget(
                french_thermal.glass_radius, Real(2) * french_thermal.glass_half_height, t0, thermal_bc, target_power);
            logFrenchZeroDHeatBudget(budget_melt, "0d_on_melt_extent");
            std::string nu_msg;
            const bool nu_ok =
                verifyFrenchNonuniformBoundaryFluxOnParticles(glass_particles, thermal_bc, french_thermal, t0, nu_msg);
            std::cout << "[ophelie][stirring-em] nonuniform_flux_check=" << nu_msg
                      << " passed=" << (nu_ok ? 1 : 0) << std::endl;
            if (!nu_ok)
            {
                std::cerr << "[ophelie][stirring-em] nonuniform boundary flux check failed" << std::endl;
                return 1;
            }
        }
        if (local_cli.balance_heat_loss && total_loss_w > TinyReal)
        {
            bc_scale = target_power / total_loss_w;
            thermal_bc.h_side *= bc_scale;
            thermal_bc.h_bottom *= bc_scale;
            thermal_bc.h_free *= bc_scale;
            thermal_bc.emissivity *= bc_scale;
            std::cout << "[ophelie][stirring-em] WARNING legacy --balance-heat-loss scale=" << bc_scale
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
    auto &t_min_reduce =
        main_methods.addReduceDynamics<OphelieThermalMinTemperatureReduceCK>(glass, kOphelieTemperatureField);

    StateDynamics<MainExecutionPolicy, rh200::SampleJouleHeatFromGridCK> sample_joule_heat(glass, q_grid_device);
    auto &out_of_grid_probe = main_methods.addReduceDynamics<rh200::OutOfGridJouleSampleReduceCK>(glass, q_grid_device);
    StateDynamics<MainExecutionPolicy, ApplyGridJouleHeatToTemperatureCK> apply_q_const(glass, rho0, cp, t0);
    StateDynamics<MainExecutionPolicy, ApplyJouleHeatWithParticleCpCK> apply_q_variable(glass, rho0, t0);
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc);
    // Doc 45 F02: operator-consistent boundary power (device reduce), not screen-cadence p_loss_live.
    ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK> boundary_loss_power(
        glass, t0, thermal_bc.boundary_shell_thickness, thermal_bc);
    StateDynamics<MainExecutionPolicy, ClampTemperatureCK> clamp_temperature(glass, local_cli.t_min, local_cli.t_max,
                                                                             t0);
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>>> laplace_temperature;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>> apply_diffusion;
    if (thermal_bc.enable_diffusion)
    {
        laplace_temperature = makeUnique<InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>>>(
            glass_inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
            thermal_bc.pair_weight_regularization);
        apply_diffusion = makeUnique<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>>(glass, rho0, cp,
                                                                                                   t0);
    }

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
    const Real n_out_of_grid_init_sample = out_of_grid_probe.exec();
    recordFrenchFluidHeatSample(coupling_state, glass_particles, n_out_of_grid_init_sample);
    if (periodic_coupling)
    {
        // Rewrite init CSV row with fluid power filled in.
        writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
        appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                  local_cli.literature.em_control);
    }
    constant_gravity.exec();
    boussinesq_force.exec();

    const Real n_out_of_grid_init = n_out_of_grid_init_sample;
    Real n_out_of_grid = n_out_of_grid_init;

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording = time_stepper.addTriggerByInterval(std::max(local_cli.state_record_interval, Real(1.0e-6)));

    fs::create_directories("./output");
    std::ofstream monitor("./output/french_stirring_em_monitor.csv");
    monitor << "advection_step,physical_time_s,revolutions,acoustic_dt_s,u_max_mps,t_mean_K,t_max_K,"
               "rotor_tip_mps,rotor_r_tip_m,f_visc_N,f_pres_N,torque_z_Nm,P_mech_W,Np_sph,Nu_sph,"
               "wall_clock_s,t_min_K,P_joule,em_updates,"
               "U_rms,U_theta_rms,U_z_rms,T_std,T_outer,T_center,T_bottom,T_top\n";
    monitor << std::setprecision(10);
    writeFrenchSpatialStatsCsvHeader("output/french_spatial_stats.csv");

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
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0)
              << " coupling=" << thermoEMCouplingModeName(local_cli.literature.coupling) << "\n"
              << "  paddle: rpm=" << rpm << " omega=" << omega << " rad/s r_tip=" << rotor0.r_tip
              << " m U_tip=" << u_tip_expected << " m/s Re_imp=" << re_impeller
              << " (creeping swirl: use U_rms/U_th on screen and output/french_spatial_stats.csv; VTP optional)\n"
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
    OphelieThermalRheologyStats rheology_stats;
    const Real mu_init = hostUpdateFrenchThermalRheologyFromTemperature(
        glass_particles, laws, local_cli.literature.constant_thermal_rheology, t0, rheology_stats);
    glass_viscosity.setReferenceViscosity(mu_init);
    invalidateViscosityDependentCkKernels(glass_inner, glass_solid_contact, &rotor_fluid_contact);
    if (!local_cli.literature.constant_thermal_rheology)
    {
        std::cout << "[ophelie][stirring-em] viscosity_mode=mean_mu_T_volume_average mu_ref=" << mu_init
                  << " (particle OphelieDynamicViscosity written; hydro force uses body mean μ after CK invalidate)"
                  << std::endl;
    }
    Real prev_energy = local_cli.literature.constant_thermal_rheology
                           ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                           : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
    Real prev_energy_time = 0.0;
    coupling_state.ledger.setH0(prev_energy);
    if (local_cli.f07_state_audit)
    {
        writeFrenchF07StateAuditCsvHeader("output/french_f07_state_audit.csv");
        std::cout << "[ophelie][stirring-em] f07_state_audit=1 -> output/french_f07_state_audit.csv\n";
    }
    auto f07_snap = [&](const char *stage) {
        if (!local_cli.f07_state_audit)
        {
            return;
        }
        const OphelieThermalMassVolumeAudit mv = hostThermalMassVolumeAudit(glass_particles, rho0);
        const Real h_rho0v = local_cli.literature.constant_thermal_rheology
                                 ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                 : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
        const Real h_mass = local_cli.literature.constant_thermal_rheology
                                ? hostThermalEnergyConstantCpFromMass(glass_particles, cp, t0)
                                : Real(0);
        appendFrenchF07StateAuditCsv("output/french_f07_state_audit.csv", stage, advection_steps, acoustic_steps,
                                     time_stepper.getPhysicalTime(), mv, h_rho0v, h_mass,
                                     static_cast<Real>(coupling_state.ledger.e_q_applied),
                                     static_cast<Real>(coupling_state.ledger.e_boundary_out));
    };
    {
        Real side0 = 0, bottom0 = 0, free_conv0 = 0, free_rad0 = 0, total0 = 0;
        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                      side0, bottom0, free_conv0, free_rad0, total0);
        coupling_state.ledger.noteLiveLoss(side0, bottom0, free_conv0 + free_rad0);
    }
    const Real glass_h = frenchReducedGlassHeight(french);
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
        if (local_cli.f07_state_audit && local_cli.csv_every > 0 &&
            advection_steps > 0 && advection_steps % static_cast<size_t>(local_cli.csv_every) == 0 &&
            acoustic_steps % 20 == 0)
        {
            // Sparse: WCSPH acoustic update of Density without neighbor summation.
            f07_snap("post_acoustic_2nd_half");
        }

        if (local_cli.literature.constant_thermal_rheology)
        {
            apply_q_const.exec(acoustic_dt);
        }
        else
        {
            apply_q_variable.exec(acoustic_dt);
        }
        if (!local_cli.literature.zero_q)
        {
            coupling_state.ledger.accumulateQ(coupling_state.last_power.p_heat_fluid, acoustic_dt);
        }
        if (!local_cli.literature.constant_thermal_rheology)
        {
            const Real mu_now = hostUpdateFrenchThermalRheologyFromTemperature(
                glass_particles, laws, false, t0, rheology_stats);
            glass_viscosity.setReferenceViscosity(mu_now);
            invalidateViscosityDependentCkKernels(glass_inner, glass_solid_contact, &rotor_fluid_contact);
        }
        if (laplace_temperature && apply_diffusion)
        {
            laplace_temperature->exec();
            apply_diffusion->exec(acoustic_dt);
        }
        // Doc 45 F02: integrate the flux the BC operator is about to apply (pre-update state).
        {
            const Real p_boundary_now = boundary_loss_power.exec();
            coupling_state.ledger.p_loss_live = p_boundary_now;
            coupling_state.ledger.accumulateBoundary(p_boundary_now, acoustic_dt);
        }
        apply_natural_bc.exec(acoustic_dt);
        clamp_temperature.exec(acoustic_dt);
        if (local_cli.f07_state_audit && local_cli.csv_every > 0 &&
            advection_steps > 0 && advection_steps % static_cast<size_t>(local_cli.csv_every) == 0 &&
            acoustic_steps % 20 == 0)
        {
            f07_snap("post_thermal_q_diff_bc");
        }
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
            // Doc 69/70: DensitySummationCK alone leaves Mass≠ρV; AdvectionStepSetup restores V=m/ρ.
            if (local_cli.f07_state_audit && local_cli.csv_every > 0 &&
                advection_steps % static_cast<size_t>(local_cli.csv_every) == 0)
            {
                f07_snap("post_density_regularization");
            }
            fluid_advection_step_setup.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
            boussinesq_force.exec();
            viscous_force_from_fluid_on_rotor.exec();

            // The melt has moved through the laboratory-fixed inductor field.
            sample_joule_heat.exec();
            {
                const Real oog_sample = out_of_grid_probe.exec();
                recordFrenchFluidHeatSample(coupling_state, glass_particles, oog_sample);
                n_out_of_grid = oog_sample;
            }
            if (periodic_coupling)
            {
                const Real t_now = time_stepper.getPhysicalTime();
                if (shouldUpdateEM(t_now, local_cli.literature.em_update_interval_s, coupling_state.next_due_time))
                {
                    syncVariableToHost<Real>(glass_particles, kOphelieTemperatureField);
                    syncVariableToHost<Vecd>(glass_particles, "Position");
                    syncVariableToHost<Real>(glass_particles, "VolumetricMeasure");
                    const Real *t_host = glass_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
                    const Vecd *pos_f = glass_particles.getVariableDataByName<Vecd>("Position");
                    const Real *vol_f = glass_particles.getVariableDataByName<Real>("VolumetricMeasure");
                    const size_t n_f = glass_particles.TotalRealParticles();
                    BaseParticles &em_particles = glass_em->getBaseParticles();
                    syncVariableToHost<Vecd>(em_particles, "Position");
                    syncVariableToHost<Real>(em_particles, "VolumetricMeasure");
                    const Vecd *pos_em = em_particles.getVariableDataByName<Vecd>("Position");
                    const Real *vol_em = em_particles.getVariableDataByName<Real>("VolumetricMeasure");
                    const size_t n_em = em_particles.TotalRealParticles();
                    StdVec<Real> sigma_on_em;
                    std::string sigma_error;
                    if (!mapFluidTemperatureToUnderRelaxedEmSigma(local_cli.literature, cyl_frame, pos_f, vol_f,
                                                                  t_host, n_f, pos_em, vol_em, n_em, coupling_state,
                                                                  sigma_on_em, sigma_error))
                    {
                        std::cerr << "[ophelie][thermo-em] abort: " << sigma_error << std::endl;
                        return 1;
                    }
                    hostAssignSigmaField(em_particles, glass_em_names.sigma, sigma_on_em);
                    OphelieFrenchEmJouleHeatOneWayResult probe;
                    runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                        *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
                    if (local_cli.literature.em_control == EMControlMode::FixedAbsorbedPower)
                    {
                        (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
                        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                            *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
                    }
                    p_joule = probe.joule_power_w;
                    phi_eq = probe.phi_eq_res_vol;
                    coupling_state.coil_current_per_loop = french.coil.current_per_loop;
                    coupling_state.last_phi_eq_res = phi_eq;
                    depositEmQAndRecord(em_particles, glass_em_names, params, cyl_frame, local_cli.literature, q_grid,
                                        coupling_state, p_joule);
                    q_grid_device.upload(q_grid);
                    sample_joule_heat.exec();
                    const Real oog = out_of_grid_probe.exec();
                    recordFrenchFluidHeatSample(coupling_state, glass_particles, oog);
                    const bool refresh_solver_ok =
                        std::isfinite(phi_eq) && phi_eq < Real(1.0e-2) && std::isfinite(p_joule) && p_joule > TinyReal;
                    noteFrenchPeriodicRefresh(coupling_state, t_now, local_cli.literature.em_update_interval_s,
                                              refresh_solver_ok);
                    ++coupling_state.update_count;
                    ++coupling_state.periodic_refresh_count;
                    n_out_of_grid = oog;
                    writeAxisymmetricCsv("output/french_sigma_rz.csv", coupling_state.sigma_rz, "sigma_Sm");
                    writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
                    appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, t_now,
                                              local_cli.literature.em_control);
                    std::cout << "[ophelie][thermo-em] update=" << coupling_state.update_count
                              << " periodic_refresh=" << coupling_state.periodic_refresh_count << " t=" << t_now
                              << " next_due=" << coupling_state.next_due_time
                              << " max_field_age=" << coupling_state.max_field_age_s
                              << " max_lateness=" << coupling_state.max_nominal_lateness_s
                              << " reconstructed_glass_power=" << p_joule
                              << " P_heat_fluid=" << coupling_state.last_power.p_heat_fluid
                              << " coil_current=" << coupling_state.coil_current_per_loop
                              << " phi_eq_res_vol=" << phi_eq << " out_of_grid=" << oog << std::endl;
                }
            }
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
                    const Real t_min_now = t_min_reduce.exec();
                    const RotorLoadSnapshot rotor_load = hostRotorLoadSnapshot(
                        rotor.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis,
                        stirring.rotation_speed_rad_s);
                    const FrenchMeltSpatialStats spatial =
                        hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
                    appendFrenchSpatialStatsCsv("output/french_spatial_stats.csv", t_now, spatial);
                    writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame,
                                               kOphelieTemperatureField);
                    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
                    Real residual = 0.0;
                    const Real dT_now = t_max_now - t_min_now;
                    const Real p_mech = frenchStirringMechanicalPowerW(rotor_load.torque_z, omega);
                    const Real np_sph =
                        frenchStirringPowerNumberSph(p_mech, rho0, omega, rotor_load.r_tip);
                    if (periodic_coupling)
                    {
                        const Real energy = local_cli.literature.constant_thermal_rheology
                                                ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                                : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
                        const Real dt_e = t_now - prev_energy_time;
                        const Real dE_dt = dt_e > TinyReal ? (energy - prev_energy) / dt_e : Real(0);
                        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc,
                                                                      thermal_bc.boundary_shell_thickness, side_w,
                                                                      bottom_w, free_conv_w, free_rad_w, total_loss_w);
                        coupling_state.ledger.noteLiveLoss(side_w, bottom_w, free_conv_w + free_rad_w);
                        const Real u_buoy = frenchBuoyancySpeedScale(9.81, beta, dT_now, glass_h);
                        const Real p_applied = coupling_state.last_power.p_heat_fluid;
                        residual = frenchEnergyResidual(dE_dt, p_applied, side_w, free_conv_w + free_rad_w, bottom_w);
                        FrenchEnergyBudgetRowExtras extras;
                        extras.p_em_recon = p_joule;
                        extras.p_heat_fluid = p_applied;
                        extras.e_q_applied = static_cast<Real>(coupling_state.ledger.e_q_applied);
                        extras.e_boundary_out = static_cast<Real>(coupling_state.ledger.e_boundary_out);
                        extras.residual_H = coupling_state.ledger.residualH(energy);
                        extras.residual_H_rel = coupling_state.ledger.residualHRel(energy);
                        extras.residual_H_rel_vs_net = coupling_state.ledger.residualHRelVsNet(energy);
                        {
                            const OphelieThermalMassVolumeAudit mv = hostThermalMassVolumeAudit(glass_particles, rho0);
                            extras.sum_mass = static_cast<Real>(mv.sum_mass);
                            extras.sum_vol = static_cast<Real>(mv.sum_vol);
                            extras.sum_rho0_vol = static_cast<Real>(mv.sum_rho0_vol);
                            extras.sum_rho_vol = static_cast<Real>(mv.sum_rho_vol);
                            extras.rho_min = static_cast<Real>(mv.rho_min);
                            extras.rho_max = static_cast<Real>(mv.rho_max);
                            extras.rho_mean = static_cast<Real>(mv.rho_mean);
                            extras.max_abs_mass_minus_rho_vol = static_cast<Real>(mv.max_abs_mass_minus_rho_vol);
                            extras.rms_mass_minus_rho_vol = static_cast<Real>(mv.rms_mass_minus_rho_vol);
                            extras.rel_l2_mass_vs_rho_vol = static_cast<Real>(mv.rel_l2_mass_vs_rho_vol);
                            if (local_cli.literature.constant_thermal_rheology)
                            {
                                extras.thermal_energy_mass =
                                    hostThermalEnergyConstantCpFromMass(glass_particles, cp, t0);
                            }
                            // Same particle state as energy-CSV Mass−ρV columns (Doc 69 P1).
                            f07_snap("csv_energy_block");
                        }
                        appendFrenchEnergyBudgetCsv("output/french_energy_budget.csv", t_now, energy, dE_dt, p_applied,
                                                    side_w, free_conv_w + free_rad_w, bottom_w, t_min_now, t_mean_now,
                                                    t_max_now, u_max_now, u_buoy,
                                                    thermal_bc.enable_diffusion ? 1 : 0, &rheology_stats,
                                                    local_cli.literature.constant_thermal_rheology ? 1 : 0, &extras);
                        prev_energy = energy;
                        prev_energy_time = t_now;
                    }
                    else
                    {
                        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc,
                                                                      thermal_bc.boundary_shell_thickness, side_w,
                                                                      bottom_w, free_conv_w, free_rad_w, total_loss_w);
                        coupling_state.ledger.noteLiveLoss(side_w, bottom_w, free_conv_w + free_rad_w);
                    }
                    const Real nu_sph =
                        frenchStirringNusseltSph(side_w + bottom_w, area_wall_wb, dT_now, glass_h, k_th);
                    if (write_csv)
                    {
                        monitor << advection_steps << "," << t_now << "," << t_now / revolution_time << "," << acoustic_dt
                                << "," << u_max_now << "," << t_mean_now << "," << t_max_now << "," << rotor_load.tip_speed
                                << "," << rotor_load.r_tip << "," << rotor_load.f_visc << "," << rotor_load.f_pres
                                << "," << rotor_load.torque_z << "," << p_mech << "," << np_sph << "," << nu_sph << ","
                                << elapsed_s() << "," << t_min_now << "," << p_joule << ","
                                << coupling_state.update_count << "," << spatial.u_rms << "," << spatial.u_theta_rms
                                << "," << spatial.u_z_rms << "," << spatial.t_std << "," << spatial.t_outer << ","
                                << spatial.t_center << "," << spatial.t_bottom << "," << spatial.t_top << "\n";
                        monitor.flush();
                    }
                    if (write_screen)
                    {
                        std::cout << std::fixed << std::setprecision(6) << "[ophelie][stirring-em] N=" << advection_steps
                                  << " n_ac=" << acoustic_steps << " t=" << t_now << " dt=" << std::scientific
                                  << acoustic_dt << std::fixed << " rev=" << t_now / revolution_time
                                  << " U_max=" << u_max_now << " U_rms=" << spatial.u_rms << " U_th=" << spatial.u_theta_rms
                                  << " U_z=" << spatial.u_z_rms << " U_tip=" << rotor_load.tip_speed
                                  << " Fv=" << rotor_load.f_visc << " Fp=" << rotor_load.f_pres
                                  << " Tz=" << rotor_load.torque_z << " P_mech=" << p_mech << " Np_sph=" << np_sph
                                  << " Nu_sph=" << nu_sph << " T_min=" << t_min_now << " T_mean=" << t_mean_now
                                  << " T_max=" << t_max_now << " T_std=" << spatial.t_std << " T_out=" << spatial.t_outer
                                  << " T_ctr=" << spatial.t_center << " em_updates=" << coupling_state.update_count
                                  << " wall=" << elapsed_s() / 3600.0 << " h" << std::endl;
                        (void)residual;
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
    const Real t_min = t_min_reduce.exec();
    const FrenchMeltSpatialStats spatial_end =
        hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
    writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame, kOphelieTemperatureField);
    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                  side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
    const RotorLoadSnapshot rotor_end =
        hostRotorLoadSnapshot(rotor.getBaseParticles(), stirring.rotation_center, stirring.rotation_axis,
                              stirring.rotation_speed_rad_s);
    const Real p_mech_end = frenchStirringMechanicalPowerW(rotor_end.torque_z, omega);
    const Real np_end = frenchStirringPowerNumberSph(p_mech_end, rho0, omega, rotor_end.r_tip);
    const Real nu_end =
        frenchStirringNusseltSph(side_w + bottom_w, area_wall_wb, t_max - t_min, glass_h, k_th);
    monitor.close();

    const Real power_rel_err = std::abs(p_joule - target_power) / (target_power + TinyReal);
    const Real physical_time = time_stepper.getPhysicalTime();
    const Real energy_end = local_cli.literature.constant_thermal_rheology
                                ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
    FrenchCouplingGateInputs gate_in;
    gate_in.periodic_coupling = periodic_coupling;
    gate_in.zero_q = local_cli.literature.zero_q;
    gate_in.constant_thermal_rheology = local_cli.literature.constant_thermal_rheology;
    gate_in.require_u_max_positive = true;
    gate_in.energy_residual_available = coupling_state.ledger.h0_set && physical_time > TinyReal;
    gate_in.budget_exhausted = budget_exhausted;
    gate_in.diverged = diverged;
    gate_in.fluid_power_recorded = coupling_state.last_power.fluid_power_recorded;
    gate_in.ledger_invalid_accumulate = coupling_state.ledger.invalid_accumulate;
    gate_in.out_of_grid_measured = true;
    gate_in.end_time = local_cli.end_time;
    gate_in.em_interval = local_cli.literature.em_update_interval_s;
    gate_in.physical_time = physical_time;
    gate_in.p_em = p_joule;
    gate_in.p_heat_fluid = coupling_state.last_power.p_heat_fluid;
    gate_in.phi_eq = phi_eq;
    gate_in.power_rel_err = power_rel_err;
    gate_in.q_sample_rel_l2 = q_sample_rel_l2;
    gate_in.n_out_of_grid_init = n_out_of_grid_init;
    gate_in.n_out_of_grid_now = n_out_of_grid;
    gate_in.max_out_of_grid = std::max(coupling_state.max_out_of_grid, coupling_state.ledger.max_out_of_grid_sampling);
    gate_in.u_max = u_max;
    gate_in.c0 = local_cli.c0;
    gate_in.t_mean = t_mean;
    gate_in.t_max = t_max;
    gate_in.cumulative_energy_residual_rel = coupling_state.ledger.residualHRel(energy_end);
    gate_in.update_count = coupling_state.update_count;
    gate_in.periodic_refresh_count = coupling_state.periodic_refresh_count;
    gate_in.freshness_measured = coupling_state.freshness_measured;
    gate_in.max_field_age_s =
        periodic_coupling
            ? std::max(coupling_state.max_field_age_s, frenchCurrentFieldAge(coupling_state, physical_time))
            : Real(-1);
    gate_in.max_nominal_lateness_s =
        periodic_coupling ? coupling_state.max_nominal_lateness_s : Real(-1);
    gate_in.refresh_solver_fail_count = coupling_state.refresh_solver_fail_count;
    FrenchCouplingRunGates gates;
    evaluateFrenchCouplingRunGates(gate_in, gates);

    const bool passed = local_cli.require_validation ? gates.validation_passed : gates.smoke_passed;
    f07_snap("endpoint_as_left_by_time_loop");
    const OphelieThermalMassVolumeAudit mv_end = hostThermalMassVolumeAudit(glass_particles, rho0);
    const Real energy_mass_end = local_cli.literature.constant_thermal_rheology
                                     ? hostThermalEnergyConstantCpFromMass(glass_particles, cp, t0)
                                     : Real(0);

    std::cout << "test_3d_ophelie_french_stirring_em"
              << " P_joule_W=" << p_joule << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << phi_eq
              << " q_grid_rel_l2=" << q_sample_rel_l2 << " out_of_grid=" << n_out_of_grid
              << " out_of_grid_max=" << coupling_state.max_out_of_grid
              << " P_heat_fluid=" << coupling_state.last_power.p_heat_fluid
              << " P_em_particles=" << coupling_state.last_power.p_em_particles
              << " E_Q_applied=" << coupling_state.ledger.e_q_applied
              << " E_boundary_out=" << coupling_state.ledger.e_boundary_out
              << " energy_residual_H=" << coupling_state.ledger.residualH(energy_end)
              << " energy_residual_H_rel=" << gate_in.cumulative_energy_residual_rel
              << " energy_residual_H_rel_vs_net=" << coupling_state.ledger.residualHRelVsNet(energy_end)
              << " sum_mass=" << mv_end.sum_mass << " sum_vol=" << mv_end.sum_vol
              << " sum_rho0_vol=" << mv_end.sum_rho0_vol << " sum_rho_vol=" << mv_end.sum_rho_vol
              << " rho_min=" << mv_end.rho_min << " rho_max=" << mv_end.rho_max
              << " rho_mean=" << mv_end.rho_mean
              << " mass_vs_rho0V_rel="
              << ((mv_end.sum_mass > 0.0)
                      ? (mv_end.sum_mass - mv_end.sum_rho0_vol) / mv_end.sum_mass
                      : 0.0)
              << " mass_vs_rhoV_rel="
              << ((mv_end.sum_mass > 0.0) ? (mv_end.sum_mass - mv_end.sum_rho_vol) / mv_end.sum_mass : 0.0)
              << " max_abs_mass_minus_rho_vol=" << mv_end.max_abs_mass_minus_rho_vol
              << " rms_mass_minus_rho_vol=" << mv_end.rms_mass_minus_rho_vol
              << " rel_l2_mass_vs_rho_vol=" << mv_end.rel_l2_mass_vs_rho_vol
              << " thermal_energy_rho0V=" << energy_end << " thermal_energy_mass=" << energy_mass_end
              << " fluid_samples=" << coupling_state.ledger.fluid_sample_count
              << " physical_time_s=" << physical_time
              << " revolutions=" << physical_time / revolution_time
              << " advection_steps=" << advection_steps << " U_max=" << u_max << " U_rms=" << spatial_end.u_rms
              << " U_th=" << spatial_end.u_theta_rms << " U_z=" << spatial_end.u_z_rms
              << " Tz=" << rotor_end.torque_z << " P_mech_W=" << p_mech_end << " Np_sph=" << np_end
              << " Nu_sph=" << nu_end << " T_min=" << t_min
              << " T_mean=" << t_mean << " T_max=" << t_max << " T_std=" << spatial_end.t_std
              << " T_out=" << spatial_end.t_outer << " T_ctr=" << spatial_end.t_center
              << " T_bot=" << spatial_end.t_bottom << " T_top=" << spatial_end.t_top << " wall_loss_side_W=" << side_w
              << " wall_loss_bottom_W=" << bottom_w << " free_conv_loss_W=" << free_conv_w
              << " free_rad_loss_W=" << free_rad_w << " total_heat_loss_W=" << total_loss_w << " bc_scale=" << bc_scale
              << " h_side_used=" << thermal_bc.h_side << " wall_clock_h=" << elapsed_s() / 3600.0
              << " acoustic_steps=" << acoustic_steps << " budget_exhausted=" << (budget_exhausted ? 1 : 0)
              << " diverged=" << (diverged ? 1 : 0) << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0)
              << " coupling=" << thermoEMCouplingModeName(local_cli.literature.coupling)
              << " coupling_updates=" << coupling_state.update_count
              << " periodic_refreshes=" << coupling_state.periodic_refresh_count
              << " expected_periodic_refreshes=" << gates.expected_periodic_refreshes
              << " startup_ok=" << (gates.startup_ok ? 1 : 0)
              << " completed_requested_time=" << (gates.completed_requested_time ? 1 : 0)
              << " periodic_refresh_ok=" << (gates.periodic_refresh_ok ? 1 : 0)
              << " freshness_ok=" << (gates.freshness_ok ? 1 : 0)
              << " max_field_age_s=" << gate_in.max_field_age_s
              << " max_nominal_lateness_s=" << gate_in.max_nominal_lateness_s
              << " refresh_solver_fail_count=" << gate_in.refresh_solver_fail_count
              << " em_solver_ok=" << (gates.em_solver_ok ? 1 : 0) << " q_transfer_ok=" << (gates.q_transfer_ok ? 1 : 0)
              << " flow_sanity_ok=" << (gates.flow_sanity_ok ? 1 : 0)
              << " thermal_energy_ok=" << (gates.thermal_energy_ok ? 1 : 0)
              << " smoke_passed=" << (gates.smoke_passed ? 1 : 0)
              << " validation_exercised=" << (gates.validation_exercised ? 1 : 0)
              << " validation_passed=" << (gates.validation_passed ? 1 : 0)
              << " require_validation=" << (local_cli.require_validation ? 1 : 0)
              << " fail_reasons=" << (gates.fail_reasons.empty() ? "none" : gates.fail_reasons)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    // Diagnostic only (after all production prints): DensitySummation to test whether
    // endpoint Mass−ρV mismatch is cleared by re-summation (Doc 69). Must not precede mv_end.
    if (local_cli.f07_state_audit)
    {
        density_regularization.exec();
        f07_snap("endpoint_after_forced_density_regularization");
    }
    return passed ? 0 : 1;
}
