/**
 * @file test_3d_ophelie_client_stirring_f1_heated_elastic.cpp
 * @brief Doc228/Doc248 C2 — client stirring heated elastic FSI.
 *
 * Builds on F0 (elastic_driven NeoHookean paddle + Reload FSI). Prescribed
 * volumetric heat (NOT EM-computed). Constant frozen μ; no Boussinesq; no solid
 * thermal strain. Rotor VTP writes VonMises stress/strain + Green–Lagrange + F.
 *
 * Stability note: stainless E=2e11 → ~500 solid substeps / fluid step and tip
 * runaway near t≈0.3 s (Doc248 C2). Long runs (≥2–5 s) use --youngs-modulus=
 * soft ASSUMED E + longer --omega-ramp= (see 03_run_C2_stir.sh). Soft E is
 * NOT a stainless strength certificate.
 *
 * Prerequisite: test_3d_ophelie_client_stirring_glass_relax --bodies=all → Reload.xml
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_stirring_geometry.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_fixed_space_probe.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_init_trz_qrz.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;

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
        "../../test_3d_ophelie_client_stirring_glass_relax/bin/reload",
        "../../../test_3d_ophelie_client_stirring_glass_relax/bin/reload",
    };
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
            return candidate;
    }
    return candidates.front();
}

struct F1Cli
{
    std::string reload_dir = resolveDefaultReloadFolder();
    std::string rotor_stl = "./data/client_stir_paddle_solver_m.stl";
    /** Short heated smoke; wall ~8× F0@0.05s → budget ~5 min. */
    Real end_time = Real(0.2);
    Real c0 = Real(8);
    /** Longer ramp reduces early tip runaway under elastic FSI. */
    Real omega_ramp_s = Real(2.0);
    Real clamp_top_band_m = Real(0.05);
    /**
     * Development prescribed absorbed power [W] for visible ΔT in short window.
     * Not generator input; not EM-computed. Override with --heat-power=.
     */
    Real heat_power_w = Real(1.0e5);
    /**
     * ≤0 → use stainless preset 2e11 Pa. Soft ASSUMED values (e.g. 2e8) needed
     * for multi-second elastic FSI; label clearly in summary.
     */
    Real youngs_modulus_pa = Real(-1);
    /** ≤0 → keep client stirrer_rpm default (10). */
    Real rpm = Real(-1);
    /** StructureIntegration1stHalf numerical damping (library default 0.125). */
    Real solid_numerical_damping = Real(0.5);
    /** Early abort if elastic offset vs rigid drive exceeds this [m]. */
    Real max_elastic_offset_m = Real(0.05);
    /** Doc248 ASSUMED: remapped French Fig.8a T / Fig.10 Q (empty = uniform T0 + uniform Q). */
    std::string init_T_rz_csv;
    std::string init_Q_rz_csv;
    Real state_record_interval = Real(0.05);
    int screen_every = 20;
    int csv_every = 20;
    Real max_wall_hours = Real(2.0);
    bool state_recording = true;
    std::string out_dir = ".";
    std::string summary_json = "doc228_f1_summary.json";
};

/** T += Q dt/(ρ cp); ΔT authoritative (same contract as french ApplyGridJouleHeat). */
class ApplyPrescribedVolumetricHeatCK : public LocalDynamics
{
  public:
    ApplyPrescribedVolumetricHeatCK(SPHBody &sph_body, Real rho, Real cp, Real t_initial)
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

inline Real hostMeltVolume(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real V = 0;
    for (size_t i = 0; i < n; ++i)
        V += vol[i];
    return V;
}

inline void hostAssignUniformJouleHeat(BaseParticles &particles, Real q_w_m3)
{
    particles.registerStateVariable<Real>("JouleHeat", q_w_m3);
    const size_t n = particles.TotalRealParticles();
    Real *q = particles.getVariableDataByName<Real>("JouleHeat");
    for (size_t i = 0; i < n; ++i)
        q[i] = q_w_m3;
    syncVariableToDevice<Real>(particles, "JouleHeat");
}

inline Real hostVolumeMeanTemperature(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_Tv = 0;
    Real sum_v = 0;
    for (size_t i = 0; i < n; ++i)
    {
        sum_Tv += T[i] * vol[i];
        sum_v += vol[i];
    }
    return sum_v > TinyReal ? sum_Tv / sum_v : Real(0);
}

inline Real hostMaxRealField(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Real>(particles, name);
    const size_t n = particles.TotalRealParticles();
    const Real *v = particles.getVariableDataByName<Real>(name);
    Real mx = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(v[i]))
            continue;
        mx = std::max(mx, std::abs(v[i]));
    }
    return mx;
}

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

/** Smooth ramp: ω(t)=ω_max·min(1,t/tramp); θ=∫ω. */
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

inline void applyClampSpinDrive(BaseParticles &particles, const StdVec<size_t> &clamp_ids,
                                const StdVec<Vecd> &pos0, const StdVec<Vecd> &n0, const Vecd &center,
                                const Vecd &axis, Real theta, Real omega, bool with_normals)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "AverageVelocity");
    if (with_normals)
        syncVariableToHost<Vecd>(particles, "NormalDirection");

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *vel_ave = particles.getVariableDataByName<Vecd>("AverageVelocity");
    Vecd *n_dir = with_normals ? particles.getVariableDataByName<Vecd>("NormalDirection") : nullptr;
    const Vecd e = axis.normalized();
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    for (size_t k = 0; k < clamp_ids.size(); ++k)
    {
        const size_t i = clamp_ids[k];
        const Vecd rotated = rotateAboutAxis(pos0[i] - center, e, c, s);
        pos[i] = center + rotated;
        vel[i] = e.cross(rotated) * omega;
        vel_ave[i] = vel[i];
        if (n_dir != nullptr)
            n_dir[i] = rotateAboutAxis(n0[i], e, c, s);
    }
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "AverageVelocity");
    if (with_normals)
        syncVariableToDevice<Vecd>(particles, "NormalDirection");
}

/** Fluid Wall/FSI reads ElasticSolid::AverageVelocity; solid integrator writes Velocity. */
inline void copyVelocityToAverageVelocity(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "AverageVelocity");
    const size_t n = particles.TotalRealParticles();
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *vel_ave = particles.getVariableDataByName<Vecd>("AverageVelocity");
    for (size_t i = 0; i < n; ++i)
        vel_ave[i] = vel[i];
    syncVariableToDevice<Vecd>(particles, "AverageVelocity");
}

struct RotorLoadSnap
{
    Real r_tip = 0;
    Real tip_speed = 0;
    Real f_visc = 0;
    Real f_pres = 0;
    Real torque_z = 0;
    Real max_elastic_offset = 0;
};

inline void hostSampleFluidAtProbe(BaseParticles &glass, const Vecd &probe, Real &T_K, Real &u_mps, Real &p_Pa)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Vecd>(glass, "Velocity");
    syncVariableToHost<Real>(glass, "Pressure");
    syncVariableToHost<Real>(glass, kOphelieTemperatureField);
    const size_t n = glass.TotalRealParticles();
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = glass.getVariableDataByName<Vecd>("Velocity");
    const Real *pres = glass.getVariableDataByName<Real>("Pressure");
    const Real *temp = glass.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real best = MaxReal;
    size_t best_i = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(pos[i][0]))
            continue;
        const Real d2 = (pos[i] - probe).squaredNorm();
        if (d2 < best)
        {
            best = d2;
            best_i = i;
        }
    }
    T_K = temp[best_i];
    u_mps = vel[best_i].norm();
    p_Pa = pres[best_i];
}

inline RotorLoadSnap hostRotorSnap(BaseParticles &particles, const StdVec<Vecd> &pos0, const Vecd &center,
                                   const Vecd &axis, Real theta)
{
    RotorLoadSnap snap;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ViscousForceFromFluid");
    syncVariableToHost<Vecd>(particles, "PressureForceFromFluid");
    const size_t n = particles.TotalRealParticles();
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *f_visc = particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    const Vecd *f_pres = particles.getVariableDataByName<Vecd>("PressureForceFromFluid");
    const Vecd e = axis.normalized();
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    Vecd f_visc_sum = Vecd::Zero();
    Vecd f_pres_sum = Vecd::Zero();
    Real torque = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(pos[i][0]))
            continue;
        const Vecd rel = pos[i] - center;
        const Real r_xy = (rel - e * e.dot(rel)).norm();
        snap.r_tip = std::max(snap.r_tip, r_xy);
        snap.tip_speed = std::max(snap.tip_speed, vel[i].norm());
        const Vecd rigid = center + rotateAboutAxis(pos0[i] - center, e, c, s);
        snap.max_elastic_offset = std::max(snap.max_elastic_offset, (pos[i] - rigid).norm());
        const Vecd f = f_visc[i] + f_pres[i];
        if (!std::isfinite(f[0]))
            continue;
        torque += e.dot(rel.cross(f));
        f_visc_sum += f_visc[i];
        f_pres_sum += f_pres[i];
    }
    snap.f_visc = f_visc_sum.norm();
    snap.f_pres = f_pres_sum.norm();
    snap.torque_z = torque;
    return snap;
}

void applyCli(int ac, char *av[], F1Cli &cli, StdVec<char *> &sph_av)
{
    sph_av.push_back(av[0]);
    for (int i = 1; i < ac; ++i)
    {
        const char *arg = av[i];
        if (std::strncmp(arg, "--reload-dir=", 13) == 0)
            cli.reload_dir = arg + 13;
        else if (std::strncmp(arg, "--rotor-stl=", 12) == 0)
            cli.rotor_stl = arg + 12;
        else if (std::strncmp(arg, "--end-time=", 11) == 0)
            cli.end_time = static_cast<Real>(std::atof(arg + 11));
        else if (std::strncmp(arg, "--c0=", 5) == 0)
            cli.c0 = static_cast<Real>(std::atof(arg + 5));
        else if (std::strncmp(arg, "--omega-ramp=", 13) == 0)
            cli.omega_ramp_s = static_cast<Real>(std::atof(arg + 13));
        else if (std::strncmp(arg, "--clamp-top-band=", 17) == 0)
            cli.clamp_top_band_m = static_cast<Real>(std::atof(arg + 17));
        else if (std::strncmp(arg, "--heat-power=", 13) == 0)
            cli.heat_power_w = static_cast<Real>(std::atof(arg + 13));
        else if (std::strncmp(arg, "--youngs-modulus=", 17) == 0)
            cli.youngs_modulus_pa = static_cast<Real>(std::atof(arg + 17));
        else if (std::strncmp(arg, "--rpm=", 6) == 0)
            cli.rpm = static_cast<Real>(std::atof(arg + 6));
        else if (std::strncmp(arg, "--solid-damping=", 16) == 0)
            cli.solid_numerical_damping = static_cast<Real>(std::atof(arg + 16));
        else if (std::strncmp(arg, "--max-elastic-offset=", 21) == 0)
            cli.max_elastic_offset_m = static_cast<Real>(std::atof(arg + 21));
        else if (std::strncmp(arg, "--init-T-rz-csv=", 16) == 0)
            cli.init_T_rz_csv = arg + 16;
        else if (std::strncmp(arg, "--init-Q-rz-csv=", 16) == 0)
            cli.init_Q_rz_csv = arg + 16;
        else if (std::strncmp(arg, "--state-record-interval=", 24) == 0)
            cli.state_record_interval = static_cast<Real>(std::atof(arg + 24));
        else if (std::strncmp(arg, "--screen-every=", 15) == 0)
            cli.screen_every = std::atoi(arg + 15);
        else if (std::strncmp(arg, "--csv-every=", 12) == 0)
            cli.csv_every = std::atoi(arg + 12);
        else if (std::strncmp(arg, "--max-wall-hours=", 17) == 0)
            cli.max_wall_hours = static_cast<Real>(std::atof(arg + 17));
        else if (std::strncmp(arg, "--out-dir=", 10) == 0)
            cli.out_dir = arg + 10;
        else if (std::strncmp(arg, "--summary-json=", 15) == 0)
            cli.summary_json = arg + 15;
        else if (std::strcmp(arg, "--no-state-recording") == 0)
            cli.state_recording = false;
        else if (std::strcmp(arg, "--help") == 0 || std::strcmp(arg, "-h") == 0)
        {
            std::cout << "Doc248 C2 / Doc228 F1 heated elastic FSI (prescribed volumetric heat):\n"
                         "  --reload-dir= --rotor-stl= --end-time= --c0= --omega-ramp=\n"
                         "  --heat-power= --clamp-top-band= --out-dir= --summary-json=\n"
                         "  --youngs-modulus= (ASSUMED soft E for multi-s FSI; ≤0 = stainless 2e11)\n"
                         "  --rpm= --solid-damping= --max-elastic-offset=\n"
                         "  --init-T-rz-csv= --init-Q-rz-csv=  (ASSUMED Fig.8a/Fig.10 remap)\n"
                         "  --max-wall-hours= --no-state-recording\n";
            std::exit(0);
        }
        else
            sph_av.push_back(const_cast<char *>(arg));
    }
}
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_client_stirring_f1_heated_elastic requires SYCL build.\n";
    return 1;
#else
    F1Cli cli;
    StdVec<char *> sph_av;
    applyCli(ac, av, cli, sph_av);
    fs::create_directories(cli.out_dir);
    fs::create_directories("./output");

    OphelieClientStirringCaseParams g_case;
    g_case.rotor_stl_solver_m = cli.rotor_stl;
    if (cli.rpm > TinyReal)
        g_case.client.stirrer_rpm = cli.rpm;
    applyClientStirringDevelopmentDefaults(g_case);
    // Reload particles use dp≈0.005 (Vol=1.25e-7); keep system resolution consistent.
    if (std::abs(g_case.french_stirring.french.dp - Real(0.005)) > Real(1.0e-6))
    {
        std::cout << "[ophelie][client-f1] note: default client dp=" << g_case.french_stirring.french.dp
                  << "; override to 0.005 to match existing Reload unless you regenerated.\n";
        g_case.client.particle_spacing_m = Real(0.005);
        applyClientStirringDevelopmentDefaults(g_case);
    }
    const bool soft_E_ASSUMED = cli.youngs_modulus_pa > TinyReal;
    if (soft_E_ASSUMED)
    {
        g_case.steel.youngs_modulus_pa = cli.youngs_modulus_pa;
        g_case.steel.source_status = "ASSUMED_soft_E_for_elastic_FSI_longevity";
        g_case.steel.material_id = "soft_neoHookean_paddle_dev_NOT_stainless_certificate";
    }

    if (!reloadXmlExists(cli.reload_dir))
    {
        std::cerr << "[ophelie][client-f1] FATAL: missing " << cli.reload_dir
                  << "/Reload.xml. Run test_3d_ophelie_client_stirring_glass_relax --bodies=all.\n";
        return 2;
    }
    if (!fs::exists(cli.rotor_stl))
    {
        std::cerr << "[ophelie][client-f1] FATAL: missing rotor STL \"" << cli.rotor_stl << "\"\n";
        return 2;
    }
    if (cli.heat_power_w <= TinyReal)
    {
        std::cerr << "[ophelie][client-f1] FATAL: --heat-power must be > 0 (prescribed absorbed W)\n";
        return 2;
    }

    const OphelieFrenchGlassMaterialLaws laws = buildClientCaseMaterialLaws(g_case.client);
    const Real rho0 = laws.rho.evaluate(g_case.client.t_initial_k);
    const Real mu = laws.mu.evaluate(g_case.client.mu_freeze_t_k);
    const Real cp = laws.cp.evaluate(g_case.client.t_initial_k);
    const Real omega_max = g_case.french_stirring.rotation_speed_rad_s;
    const Vecd axis = g_case.french_stirring.rotation_axis;
    const Vecd center = g_case.french_stirring.rotation_center;
    const Real dp = g_case.french_stirring.french.dp;

    writeClientStainless304Json(g_case.steel, (fs::path(cli.out_dir) / "stainless_304_dev_preset.json").string());

    const BoundingBoxd flow_bounds = clientStirringRelaxDomainBounds(g_case);
    SPHSystem sph_system(flow_bounds, dp);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    SolidBody rotor(sph_system, makeShared<OphelieClientStirringRotorShape>("Rotor", g_case));
    rotor.defineAdaptationRatios(1.3, 1.0);
    rotor.defineMatterMaterial<NeoHookeanSolid>(g_case.steel.rho_kg_m3, g_case.steel.youngs_modulus_pa,
                                                 g_case.steel.poisson);
    rotor.generateParticles<BaseParticles, Reload>(rotor.Name());
    // NeoHookean: structure integration uses "Velocity"; Wall/FSI AverageVelocityVariable
    // registers "AverageVelocity". Both must exist before CK Wall / solid dynamics construct.
    rotor.getBaseParticles().registerStateVariable<Vecd>("Velocity");
    rotor.getBaseParticles().registerStateVariable<Vecd>("AverageVelocity");
    rotor.getBaseParticles().registerStateVariable<Vecd>("AverageAcceleration");
    rotor.getBaseParticles().addEvolvingVariable<Vecd>("Velocity");
    rotor.getBaseParticles().addEvolvingVariable<Vecd>("AverageVelocity");

    OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
    melt.melt_height_m = g_case.melt_height_m;
    melt.glass_half_height = Real(0.5) * melt.melt_height_m;
    melt.glass_center_z_m = melt.inner_bottom_z_m + melt.glass_half_height;
    melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);
    SolidBody wall(sph_system, makeShared<OphelieClientNcCrucibleWallShape>(
                                   "WallBoundary", melt, Real(4) * dp, 20));
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    (void)wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall.generateParticles<BaseParticles, Reload>(wall.Name());

    FluidBody glass(sph_system, makeShared<OphelieClientStirringGlassFluidShape>("GlassBody", g_case));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, cli.c0);
    glass.addMaterialProperty<Viscosity>(mu);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());
    const Real t0 = g_case.client.t_initial_k;
    BaseParticles &glass_particles = glass.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(glass_particles, t0);
    glass_particles.addVariableToWrite<Real>(kOphelieTemperatureField);

    const Vecd melt_center = g_case.french_stirring.french.glass_center;
    const Real melt_half_h = g_case.french_stirring.french.glass_half_height;
    const Real melt_R = g_case.melt_radius_m;
    if (!cli.init_T_rz_csv.empty())
    {
        Doc248RzTemperatureGrid t_grid;
        std::string terr;
        if (!hostLoadDoc248RzTemperatureGrid(cli.init_T_rz_csv, t_grid, terr) || !t_grid.ok)
        {
            std::cerr << "[ophelie][client-f1] FATAL: " << terr << " path=" << cli.init_T_rz_csv << "\n";
            return 2;
        }
        hostApplyInitTRzAbsolute(glass_particles, melt_center, t0, t_grid, "client-f1/init-T-rz");
    }

    const Real melt_volume = hostMeltVolume(glass_particles);
    Real p_heat_realized = cli.heat_power_w;
    Real q_uniform = cli.heat_power_w / std::max(melt_volume, TinyReal);
    if (!cli.init_Q_rz_csv.empty())
    {
        Doc247RzQGrid q_grid;
        std::string qerr;
        if (!hostLoadDoc247RzQGrid(cli.init_Q_rz_csv, q_grid, qerr) || !q_grid.ok)
        {
            std::cerr << "[ophelie][client-f1] FATAL: " << qerr << " path=" << cli.init_Q_rz_csv << "\n";
            return 2;
        }
        // register then overwrite with shaped field
        hostAssignUniformJouleHeat(glass_particles, Real(0));
        p_heat_realized = hostAssignPrescribedShapeJouleHeat(
            glass_particles, "JouleHeat", q_grid, melt_center, melt_R, melt_half_h, cli.heat_power_w,
            "client-f1/init-Q-rz");
        q_uniform = p_heat_realized / std::max(melt_volume, TinyReal); // mean Q after shape scale
    }
    else
    {
        hostAssignUniformJouleHeat(glass_particles, q_uniform);
    }
    glass_particles.addVariableToWrite<Real>("JouleHeat");

    // Legacy CAD tip (AABB) + Doc248 midradius virtual probe (r=0.5R, depth 0.12 m).
    const OphelieClientProbeDefinition thermo_probe_def = clientStirringThermocoupleProbe();
    const Vecd thermo_probe_pos = clientStirringThermocoupleProbePosition();
    const OphelieClientProbeDefinition mid_probe_def = makeClientMidradiusDepth12cmProbe(
        g_case.melt_radius_m, g_case.melt_bottom_z_m, g_case.melt_height_m, Real(1673.15));
    const Vecd mid_probe_pos = Vecd(mid_probe_def.solver_m.x, mid_probe_def.solver_m.y, mid_probe_def.solver_m.z);
    ObserverBody thermocouple_probe(sph_system, "ThermocoupleProbe");
    thermocouple_probe.getSPHAdaptation().resetAdaptationRatios(1.15, 2.0);
    thermocouple_probe.generateParticles<ObserverParticles>(StdVec<Vecd>{thermo_probe_pos});

    const size_t n_rotor = rotor.getBaseParticles().TotalRealParticles();
    const size_t n_wall = wall.getBaseParticles().TotalRealParticles();
    const size_t n_glass = glass_particles.TotalRealParticles();
    std::cout << "[ophelie][client-f1] reload: n_glass=" << n_glass << " n_rotor=" << n_rotor
              << " n_wall=" << n_wall << " dir=" << cli.reload_dir << "\n";
    if (n_rotor == 0 || n_wall == 0 || n_glass == 0)
    {
        std::cerr << "[ophelie][client-f1] empty body in Reload.xml\n";
        return 1;
    }

    // Clamp: provisional top band (same convention as E0).
    Real z_max = -std::numeric_limits<Real>::max();
    {
        const Vecd *pos = rotor.getBaseParticles().getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < n_rotor; ++i)
            z_max = std::max(z_max, pos[i][2]);
    }
    const Real z_clamp = z_max - cli.clamp_top_band_m;
    ClampTopBand clamp_part(rotor, z_clamp);
    const StdVec<size_t> clamp_ids = clamp_part.body_part_particles_;
    const size_t n_clamp = clamp_ids.size();
    if (n_clamp == 0 || n_clamp >= n_rotor)
    {
        std::cerr << "[ophelie][client-f1] FATAL: clamp empty or whole body n_clamp=" << n_clamp
                  << " n=" << n_rotor << " z_clamp=" << z_clamp << "\n";
        return 2;
    }

    Inner<> glass_inner(glass);
    Contact<> glass_solid_contact(glass, {&wall, &rotor});
    Contact<> rotor_fluid_contact(rotor, {&glass});
    Contact<> thermocouple_probe_contact(thermocouple_probe, {&glass});
    auto &rotor_inner = sph_system.addInnerRelation(rotor, ConfigType::Lagrangian);

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
    auto &update_thermocouple_probe_contact = main_methods.addRelationDynamics(thermocouple_probe_contact);
    auto &update_rotor_inner = main_methods.addRelationDynamics(rotor_inner);
    auto &fluid_particle_sorting = main_methods.addSortDynamics(glass);

    Gravity gravity(Vecd(0.0, 0.0, -Real(9.81)));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(glass, gravity);
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

    const Real u_ref = std::max(omega_max * Real(0.15), Real(0.05));
    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(glass, u_ref);
    auto &acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(glass, 0.55);
    auto &u_max_reduce = main_methods.addReduceDynamics<OphelieVelocityMaxReduceCK>(glass);

    auto &rotor_linear_correction =
        main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(rotor_inner);
    auto &solid_1st =
        main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid,
                                                    NoKernelCorrectionCK>(rotor_inner, cli.solid_numerical_damping);
    auto &solid_2nd =
        main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(rotor_inner);
    auto &solid_acoustic_dt = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(rotor, 0.5);
    auto &update_rotor_normal = main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(rotor);

    auto &glass_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass);
    glass_state_recorder.addToWrite<Vecd>(glass, "Velocity");
    glass_state_recorder.addToWrite<Real>(glass, "Pressure");
    glass_state_recorder.addToWrite<Real>(glass, kOphelieTemperatureField);
    glass_state_recorder.addToWrite<Real>(glass, "JouleHeat");
    // Host derived fields (non-CK LocalDynamics): construct first so names exist, then VTP write list.
    SimpleDynamics<Displacement> displacement_host(rotor);
    SimpleDynamics<VonMisesStress> von_mises_host(rotor);
    SimpleDynamics<VonMisesStrain> von_strain_host(rotor);
    SimpleDynamics<GreenLagrangeStrain> green_lagrange_host(rotor);

    auto &rotor_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(rotor);
    rotor_state_recorder.addToWrite<Vecd>(rotor, "Velocity");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "ViscousForceFromFluid");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "PressureForceFromFluid");
    rotor_state_recorder.addToWrite<Matd>(rotor, "DeformationGradient");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "Displacement");
    rotor_state_recorder.addToWrite<Real>(rotor, "VonMisesStress");
    rotor_state_recorder.addToWrite<Real>(rotor, "VonMisesStrain");
    rotor_state_recorder.addToWrite<Matd>(rotor, "GreenLagrangeStrain");

    auto refresh_rotor_stress_strain = [&]() {
        // Solid CK updates live on device; host SimpleDynamics needs synced F/ρ/x.
        syncVariableToHost<Matd>(rotor.getBaseParticles(), "DeformationGradient");
        syncVariableToHost<Real>(rotor.getBaseParticles(), "Density");
        syncVariableToHost<Vecd>(rotor.getBaseParticles(), "Position");
        displacement_host.exec();
        von_mises_host.exec();
        von_strain_host.exec();
        green_lagrange_host.exec();
    };

    StateDynamics<execution::MainExecutionPolicy, ApplyPrescribedVolumetricHeatCK> apply_prescribed_heat(
        glass, rho0, cp, t0);

    // Same pattern as dambreak / heat_transfer probes: Observer samples neighboring fluid.
    auto &write_probe_temperature =
        main_methods.addObserveRecorder<Real>(thermocouple_probe_contact, kOphelieTemperatureField);
    auto &write_probe_velocity =
        main_methods.addObserveRecorder<Vecd>(thermocouple_probe_contact, "Velocity");
    auto &write_probe_pressure =
        main_methods.addObserveRecorder<Real>(thermocouple_probe_contact, "Pressure");

    // Rest pose for clamp drive + elastic offset (host reload coordinates).
    StdVec<Vecd> pos0(n_rotor);
    StdVec<Vecd> n0(n_rotor);
    {
        const Vecd *pos = rotor.getBaseParticles().getVariableDataByName<Vecd>("Position");
        const Vecd *n_dir = rotor.getBaseParticles().getVariableDataByName<Vecd>("NormalDirection");
        for (size_t i = 0; i < n_rotor; ++i)
        {
            pos0[i] = pos[i];
            n0[i] = n_dir[i];
        }
    }

    update_rotor_cell_linked_list.exec();
    update_rotor_inner.exec();
    rotor_linear_correction.exec();

    update_glass_cell_linked_list.exec();
    update_wall_cell_linked_list.exec();
    update_glass_relations.exec();
    update_rotor_fluid_contact.exec();
    update_thermocouple_probe_contact.exec();
    density_regularization.exec();
    fluid_advection_step_setup.exec();
    constant_gravity.exec();
    fluid_viscous_force.exec();
    viscous_force_from_fluid_on_rotor.exec();

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording = time_stepper.addTriggerByInterval(std::max(cli.state_record_interval, Real(1.0e-6)));

    std::ofstream monitor("./output/client_f1_heated_elastic_monitor.csv");
    monitor << "advection_step,physical_time_s,omega,theta,acoustic_dt_s,solid_dt_s,solid_substeps,"
               "u_max_mps,r_tip_m,tip_mps,f_visc_N,f_pres_N,torque_z_Nm,elastic_offset_m,"
               "von_mises_max_Pa,von_mises_strain_max,"
               "T_mean_K,probe_T_K,probe_u_mps,probe_p_Pa,"
               "T_glass_virtual_midradius_K,midradius_n_neighbors,midradius_valid,"
               "E_heat_J,wall_clock_s\n";
    monitor << std::setprecision(10);
    {
        std::ofstream map(fs::path(cli.out_dir) / "geometry_and_probe_mapping.json");
        map << std::setprecision(10) << "{\n"
            << "  \"case\": \"doc248_C2_client_stir_f1\",\n"
            << "  \"melt_radius_m\": " << g_case.melt_radius_m << ",\n"
            << "  \"melt_height_m\": " << g_case.melt_height_m << ",\n"
            << "  \"melt_bottom_z_m\": " << g_case.melt_bottom_z_m << ",\n"
            << "  \"paddle_axis_x_m\": " << g_case.paddle_axis_x_m << ",\n"
            << "  \"legacy_thermocouple_tip_m\": [" << thermo_probe_pos[0] << "," << thermo_probe_pos[1] << ","
            << thermo_probe_pos[2] << "],\n"
            << "  \"client_midradius_probe_m\": [" << mid_probe_pos[0] << "," << mid_probe_pos[1] << ","
            << mid_probe_pos[2] << "],\n"
            << "  \"client_midradius_status\": \"" << mid_probe_def.status << "\",\n"
            << "  \"note\": \"midradius = 0.5R, depth 0.12 m below free surface, +X azimuth ASSUMED\"\n"
            << "}\n";
    }

    write_probe_temperature.writeToFile(0);
    write_probe_velocity.writeToFile(0);
    write_probe_pressure.writeToFile(0);

    if (cli.state_recording)
    {
        refresh_rotor_stress_strain();
        glass_state_recorder.writeToFile(0);
        rotor_state_recorder.writeToFile(0);
    }

    const Real dt_s0 = solid_acoustic_dt.exec();
    const Real dT_expected =
        cli.heat_power_w * cli.end_time / std::max(rho0 * cp * melt_volume, TinyReal);
    std::cout << "[ophelie][client-f1] heated elastic FSI start\n"
              << "  mechanics=" << clientPaddleMechanicsName(g_case.paddle_mechanics)
              << " E=" << g_case.steel.youngs_modulus_pa
              << (soft_E_ASSUMED ? " (ASSUMED_soft_E)" : " (stainless_dev)")
              << " rho_s=" << g_case.steel.rho_kg_m3
              << " solid_damp=" << cli.solid_numerical_damping
              << " n_clamp=" << n_clamp << " clamp_frac=" << (Real(n_clamp) / Real(n_rotor))
              << " z_clamp=" << z_clamp << "\n"
              << "  rho_f=" << rho0 << " mu=" << mu << " cp=" << cp << " c0=" << cli.c0 << " dp=" << dp
              << " rpm=" << g_case.client.stirrer_rpm
              << " omega_max=" << omega_max << " ramp=" << cli.omega_ramp_s << " s end=" << cli.end_time
              << " solid_dt0=" << dt_s0 << "\n"
              << "  heat_mode=prescribed_uniform_volumetric P=" << cli.heat_power_w
              << " W V_melt=" << melt_volume << " m3 Q=" << q_uniform << " W/m3"
              << " dT_expected~" << dT_expected << " K (adiabatic estimate)\n"
              << "  rotor VTP: Velocity, FSI forces, DeformationGradient, Displacement,"
                 " VonMisesStress, VonMisesStrain, GreenLagrangeStrain\n"
              << "  thermocouple probe id=" << thermo_probe_def.probe_id
              << " stl=" << clientStirringThermocoupleStlName()
              << " status=" << thermo_probe_def.status
              << " pos_m=[" << thermo_probe_pos[0] << "," << thermo_probe_pos[1] << "," << thermo_probe_pos[2]
              << "] T0=" << t0 << " K\n"
              << "  doc248 midradius id=" << mid_probe_def.probe_id << " status=" << mid_probe_def.status
              << " pos_m=[" << mid_probe_pos[0] << "," << mid_probe_pos[1] << "," << mid_probe_pos[2] << "]\n"
              << "  volumetric heat on; no Boussinesq; no solid thermal strain; gravity on\n";

    const auto wall_clock_start = std::chrono::steady_clock::now();
    const double wall_clock_budget_s = static_cast<double>(cli.max_wall_hours) * 3600.0;
    auto elapsed_s = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_clock_start).count();
    };

    size_t advection_steps = 0;
    size_t acoustic_steps = 0;
    size_t solid_steps_total = 0;
    bool budget_exhausted = false;
    bool diverged = false;
    Real last_omega = 0;
    Real last_theta = 0;
    Real last_elastic = 0;
    Real last_torque = 0;
    Real last_u_max = 0;
    Real last_von_mises = 0;
    Real last_von_strain = 0;
    Real last_T_mean = t0;
    Real last_probe_T = t0;
    Real E_heat_applied = 0;
    RotorLoadSnap last_load{};

    auto run_solid_subcycle = [&](Real dt_fluid, Real physical_time) {
        Real dt_s = solid_acoustic_dt.exec();
        Real summed = 0;
        size_t n_sub = 0;
        while (summed < dt_fluid - TinyReal)
        {
            Real dt_use = dt_s;
            if (summed + dt_use > dt_fluid)
                dt_use = dt_fluid - summed;
            Real omega = 0;
            Real theta = 0;
            omegaThetaRamp(physical_time + summed + Real(0.5) * dt_use, cli.omega_ramp_s, omega_max, omega, theta);
            solid_1st.exec(dt_use);
            applyClampSpinDrive(rotor.getBaseParticles(), clamp_ids, pos0, n0, center, axis, theta, omega, true);
            solid_2nd.exec(dt_use);
            applyClampSpinDrive(rotor.getBaseParticles(), clamp_ids, pos0, n0, center, axis, theta, omega, true);
            summed += dt_use;
            ++n_sub;
            ++solid_steps_total;
            last_omega = omega;
            last_theta = theta;
            if (n_sub > 0 && (n_sub % 200 == 0))
                dt_s = solid_acoustic_dt.exec();
        }
        // Free particles: publish Velocity to AverageVelocity for fluid Wall / FSI.
        copyVelocityToAverageVelocity(rotor.getBaseParticles());
        rotor.setNewlyUpdated();
        return n_sub;
    };

    while (!time_stepper.isEndTime(cli.end_time))
    {
        if (elapsed_s() > wall_clock_budget_s)
        {
            budget_exhausted = true;
            break;
        }

        const Real dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        ++acoustic_steps;
        if (!std::isfinite(dt) || dt < Real(1.0e-9) || !std::isfinite(time_stepper.getPhysicalTime()))
        {
            std::cerr << "[ophelie][client-f1] diverged at acoustic step " << acoustic_steps
                      << " t=" << time_stepper.getPhysicalTime() << " dt=" << dt << "\n";
            diverged = true;
            break;
        }

        acoustic_step_1st_half.exec(dt);
        pressure_force_from_fluid_on_rotor.exec(dt);
        acoustic_step_2nd_half.exec(dt);
        apply_prescribed_heat.exec(dt);
        E_heat_applied += cli.heat_power_w * dt;

        // Solid subcycle uses time at start of this acoustic interval.
        const Real t_start = time_stepper.getPhysicalTime() - dt;
        const size_t n_sub = run_solid_subcycle(dt, t_start);

        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            update_rotor_normal.exec();

            if (cli.state_recording && state_recording())
            {
                refresh_rotor_stress_strain();
                glass_state_recorder.writeToFile(advection_steps);
                rotor_state_recorder.writeToFile(advection_steps);
            }

            if (advection_steps % 50 == 0)
                fluid_particle_sorting.exec();

            update_glass_cell_linked_list.exec();
            update_rotor_cell_linked_list.exec();
            update_glass_relations.exec();
            update_rotor_fluid_contact.exec();
            update_thermocouple_probe_contact.exec();
            density_regularization.exec();
            fluid_advection_step_setup.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
            viscous_force_from_fluid_on_rotor.exec();

            const bool write_csv = cli.csv_every > 0 && advection_steps % static_cast<size_t>(cli.csv_every) == 0;
            const bool write_screen =
                cli.screen_every > 0 && advection_steps % static_cast<size_t>(cli.screen_every) == 0;
            if (write_csv || write_screen)
            {
                const Real t = time_stepper.getPhysicalTime();
                last_u_max = u_max_reduce.exec();
                last_load = hostRotorSnap(rotor.getBaseParticles(), pos0, center, axis, last_theta);
                last_elastic = last_load.max_elastic_offset;
                last_torque = last_load.torque_z;
                refresh_rotor_stress_strain();
                last_von_mises = hostMaxRealField(rotor.getBaseParticles(), "VonMisesStress");
                last_von_strain = hostMaxRealField(rotor.getBaseParticles(), "VonMisesStrain");
                Real probe_T = t0;
                Real probe_u = 0;
                Real probe_p = 0;
                hostSampleFluidAtProbe(glass_particles, thermo_probe_pos, probe_T, probe_u, probe_p);
                const Real h_probe = Real(1.3) * dp;
                const auto mid_s =
                    hostSampleFixedSpaceTemperature(glass_particles, mid_probe_pos, h_probe, kOphelieTemperatureField);
                last_T_mean = hostVolumeMeanTemperature(glass_particles);
                last_probe_T = probe_T;
                write_probe_temperature.writeToFile(advection_steps);
                write_probe_velocity.writeToFile(advection_steps);
                write_probe_pressure.writeToFile(advection_steps);
                const Real wall_now = static_cast<Real>(elapsed_s());
                if (write_csv)
                {
                    monitor << advection_steps << "," << t << "," << last_omega << "," << last_theta << "," << dt
                            << "," << solid_acoustic_dt.exec() << "," << n_sub << "," << last_u_max << ","
                            << last_load.r_tip << "," << last_load.tip_speed << "," << last_load.f_visc << ","
                            << last_load.f_pres << "," << last_load.torque_z << "," << last_elastic << ","
                            << last_von_mises << "," << last_von_strain << ","
                            << last_T_mean << "," << probe_T << "," << probe_u << "," << probe_p << ","
                            << mid_s.T << "," << mid_s.n_neighbors << "," << (mid_s.valid ? 1 : 0) << ","
                            << E_heat_applied << "," << wall_now << "\n";
                    monitor.flush();
                }
                if (write_screen)
                {
                    std::cout << "[ophelie][client-f1] adv=" << advection_steps << " t=" << t
                              << " dt=" << dt << " solid_sub=" << n_sub << " u_max=" << last_u_max
                              << " tip=" << last_load.tip_speed << " Tq=" << last_torque
                              << " e_off=" << last_elastic << " VM=" << last_von_mises
                              << " T_mean=" << last_T_mean
                              << " probe_T=" << probe_T << " mid_T=" << mid_s.T
                              << " E_heat=" << E_heat_applied << " wall_s=" << wall_now << "\n";
                }
                const Real tip_rigid = last_omega * std::max(last_load.r_tip, TinyReal);
                const bool tip_runaway =
                    last_load.tip_speed > std::max(Real(5) * tip_rigid, Real(2)) && last_load.tip_speed > Real(1);
                if (!std::isfinite(last_u_max) || !std::isfinite(last_elastic) || !std::isfinite(last_von_mises) ||
                    last_u_max > Real(50) || last_elastic > cli.max_elastic_offset_m || tip_runaway)
                {
                    std::cerr << "[ophelie][client-f1] diverge trip t=" << t << " u_max=" << last_u_max
                              << " e_off=" << last_elastic << " tip=" << last_load.tip_speed
                              << " tip_rigid~" << tip_rigid << " VM=" << last_von_mises << "\n";
                    diverged = true;
                    break;
                }
            }
        }
    }

    const Real wall_s = static_cast<Real>(elapsed_s());
    const Real t_final = time_stepper.getPhysicalTime();
    last_load = hostRotorSnap(rotor.getBaseParticles(), pos0, center, axis, last_theta);
    last_u_max = u_max_reduce.exec();
    last_elastic = last_load.max_elastic_offset;
    last_torque = last_load.torque_z;
    refresh_rotor_stress_strain();
    last_von_mises = hostMaxRealField(rotor.getBaseParticles(), "VonMisesStress");
    last_von_strain = hostMaxRealField(rotor.getBaseParticles(), "VonMisesStrain");
    last_T_mean = hostVolumeMeanTemperature(glass_particles);
    {
        Real probe_u = 0;
        Real probe_p = 0;
        hostSampleFluidAtProbe(glass_particles, thermo_probe_pos, last_probe_T, probe_u, probe_p);
        (void)probe_u;
        (void)probe_p;
    }
    const Real dT_mean = last_T_mean - t0;
    const Real E_thermal =
        hostOphelieJouleHeatOneWayTotalThermalEnergy(glass_particles, kOphelieTemperatureField, t0, rho0, cp,
                                                     n_glass, false);
    const Real E_expected = cli.heat_power_w * t_final;
    const Real energy_rel =
        E_expected > TinyReal ? std::abs(E_thermal - E_expected) / E_expected : Real(1);

    const bool finite_ok = std::isfinite(last_u_max) && std::isfinite(last_elastic) && std::isfinite(last_torque) &&
                           std::isfinite(last_T_mean) && std::isfinite(last_probe_T);
    const bool clamp_ok = (n_clamp > 0) && (n_clamp < n_rotor);
    const bool coupled_ok = (last_load.f_pres > TinyReal || last_load.f_visc > TinyReal || last_omega > TinyReal);
    const bool time_ok = t_final >= Real(0.8) * cli.end_time || budget_exhausted;
    // Remapped Fig.8a IC already embeds ΔT vs t0 → energy_rel gate is meaningless then.
    const bool heat_ok = cli.init_T_rz_csv.empty()
                             ? ((dT_mean > Real(0.05)) && (energy_rel < Real(0.15)))
                             : (std::isfinite(last_T_mean) && E_heat_applied > TinyReal && !diverged);
    const bool pass =
        finite_ok && clamp_ok && coupled_ok && heat_ok && !diverged && solid_steps_total > 0 && time_ok;

    std::cout << "[ophelie][client-f1] done t=" << t_final << " wall_s=" << wall_s
              << " solid_steps=" << solid_steps_total << " u_max=" << last_u_max
              << " tip=" << last_load.tip_speed << " Tq=" << last_torque << " e_off=" << last_elastic
              << " VM=" << last_von_mises << " VMeps=" << last_von_strain
              << " T_mean=" << last_T_mean << " dT=" << dT_mean << " probe_T=" << last_probe_T
              << " E_th=" << E_thermal << " E_exp=" << E_expected << " energy_rel=" << energy_rel
              << " budget_exhausted=" << (budget_exhausted ? 1 : 0) << " diverged=" << (diverged ? 1 : 0)
              << " HEATED_ELASTIC_INTEGRATION_PASS=" << (pass ? 1 : 0) << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"client_stirring_f1_heated_elastic\",\n";
        js << "  \"tag\": \"HEATED_ELASTIC_INTEGRATION\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"HEATED_ELASTIC_INTEGRATION_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"mechanics_mode\": \"elastic_driven\",\n";
        js << "  \"heat_source_mode\": \""
           << (cli.init_Q_rz_csv.empty() ? "prescribed_uniform_volumetric" : "prescribed_shaped_rz") << "\",\n";
        js << "  \"heat_source_note\": \"NOT EM-computed; ASSUMED prescribed absorbed power"
           << (cli.init_Q_rz_csv.empty() ? "" : " with Fig.10-shaped Q remap") << "\",\n";
        js << "  \"init_T_rz_csv\": \"" << cli.init_T_rz_csv << "\",\n";
        js << "  \"init_Q_rz_csv\": \"" << cli.init_Q_rz_csv << "\",\n";
        js << "  \"init_T_ASSUMED\": " << (cli.init_T_rz_csv.empty() ? 0 : 1) << ",\n";
        js << "  \"reload_dir\": \"" << cli.reload_dir << "\",\n";
        js << "  \"n_glass\": " << n_glass << ",\n";
        js << "  \"n_rotor\": " << n_rotor << ",\n";
        js << "  \"n_wall\": " << n_wall << ",\n";
        js << "  \"n_clamp\": " << n_clamp << ",\n";
        js << "  \"clamp_top_band_m\": " << cli.clamp_top_band_m << ",\n";
        js << "  \"clamp_status\": \"provisional_top_band_unconfirmed_drive_section\",\n";
        js << "  \"E_Pa\": " << g_case.steel.youngs_modulus_pa << ",\n";
        js << "  \"E_ASSUMED_soft\": " << (soft_E_ASSUMED ? 1 : 0) << ",\n";
        js << "  \"solid_numerical_damping\": " << cli.solid_numerical_damping << ",\n";
        js << "  \"rho_solid\": " << g_case.steel.rho_kg_m3 << ",\n";
        js << "  \"rho_fluid\": " << rho0 << ",\n";
        js << "  \"mu_fluid\": " << mu << ",\n";
        js << "  \"cp_fluid\": " << cp << ",\n";
        js << "  \"c0\": " << cli.c0 << ",\n";
        js << "  \"dp_m\": " << dp << ",\n";
        js << "  \"rpm\": " << g_case.client.stirrer_rpm << ",\n";
        js << "  \"melt_volume_m3\": " << melt_volume << ",\n";
        js << "  \"heat_power_W\": " << cli.heat_power_w << ",\n";
        js << "  \"Q_uniform_Wpm3\": " << q_uniform << ",\n";
        js << "  \"omega_max_rad_s\": " << omega_max << ",\n";
        js << "  \"omega_ramp_s\": " << cli.omega_ramp_s << ",\n";
        js << "  \"end_time_s\": " << cli.end_time << ",\n";
        js << "  \"physical_time_s\": " << t_final << ",\n";
        js << "  \"solid_steps_total\": " << solid_steps_total << ",\n";
        js << "  \"u_max_mps\": " << last_u_max << ",\n";
        js << "  \"r_tip_m\": " << last_load.r_tip << ",\n";
        js << "  \"tip_speed_mps\": " << last_load.tip_speed << ",\n";
        js << "  \"f_visc_N\": " << last_load.f_visc << ",\n";
        js << "  \"f_pres_N\": " << last_load.f_pres << ",\n";
        js << "  \"torque_z_Nm\": " << last_torque << ",\n";
        js << "  \"max_elastic_offset_m\": " << last_elastic << ",\n";
        js << "  \"von_mises_max_Pa\": " << last_von_mises << ",\n";
        js << "  \"von_mises_strain_max\": " << last_von_strain << ",\n";
        js << "  \"rotor_vtp_fields\": ["
              "\"Velocity\",\"ViscousForceFromFluid\",\"PressureForceFromFluid\","
              "\"DeformationGradient\",\"Displacement\",\"VonMisesStress\","
              "\"VonMisesStrain\",\"GreenLagrangeStrain\"],\n";
        js << "  \"T0_K\": " << t0 << ",\n";
        js << "  \"T_mean_K\": " << last_T_mean << ",\n";
        js << "  \"dT_mean_K\": " << dT_mean << ",\n";
        js << "  \"probe_T_K\": " << last_probe_T << ",\n";
        js << "  \"E_heat_applied_J\": " << E_heat_applied << ",\n";
        js << "  \"E_thermal_J\": " << E_thermal << ",\n";
        js << "  \"E_expected_J\": " << E_expected << ",\n";
        js << "  \"energy_rel_err\": " << energy_rel << ",\n";
        js << "  \"finite_ok\": " << (finite_ok ? 1 : 0) << ",\n";
        js << "  \"heat_ok\": " << (heat_ok ? 1 : 0) << ",\n";
        js << "  \"budget_exhausted\": " << (budget_exhausted ? 1 : 0) << ",\n";
        js << "  \"diverged\": " << (diverged ? 1 : 0) << ",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"em_enabled\": 0,\n";
        js << "  \"boussinesq_enabled\": 0,\n";
        js << "  \"solid_thermal_strain_enabled\": 0,\n";
        js << "  \"thermocouple_probe\": {\n";
        js << "    \"id\": \"" << thermo_probe_def.probe_id << "\",\n";
        js << "    \"stl\": \"" << clientStirringThermocoupleStlName() << "\",\n";
        js << "    \"status\": \"" << thermo_probe_def.status << "\",\n";
        js << "    \"solver_m\": [" << thermo_probe_pos[0] << ", " << thermo_probe_pos[1] << ", "
           << thermo_probe_pos[2] << "],\n";
        js << "    \"target_temperature_k\": " << thermo_probe_def.target_temperature_k << ",\n";
        js << "    \"target_meaning\": \"" << thermo_probe_def.target_meaning << "\",\n";
        js << "    \"probe_T_K\": " << last_probe_T << "\n";
        js << "  },\n";
        js << "  \"note\": \"F1/C2 prescribed volumetric heat + elastic FSI; not EM Q; "
           << (soft_E_ASSUMED
                   ? "ASSUMED soft E for multi-second FSI longevity (NOT stainless certificate)"
                   : "stainless E=2e11 — expect early tip runaway without soft E")
           << "\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][client-f1] summary_json=" << summary_path.string() << "\n";
    return pass ? 0 : 1;
#endif
}
