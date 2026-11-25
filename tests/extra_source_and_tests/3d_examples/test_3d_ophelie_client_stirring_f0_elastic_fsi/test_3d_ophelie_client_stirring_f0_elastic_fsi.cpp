/**
 * @file test_3d_ophelie_client_stirring_f0_elastic_fsi.cpp
 * @brief Doc228 F0 — client stirring isothermal elastic FSI short run.
 *
 * Structure mirrors french_stirring_em Phase B (GlassBody / Rotor / WallBoundary Reload),
 * but Rotor is NeoHookean stainless 304 with clamp-driven rotation (elastic_driven).
 * No EM, no thermal source, no Boussinesq. Gravity kept for melt hydrostatics.
 *
 * Prerequisite: test_3d_ophelie_client_stirring_glass_relax --bodies=all → Reload.xml
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_stirring_geometry.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_material_laws.h"
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
#include <sstream>
#include <string>
#include <vector>

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

/** Doc241: extreme finite sentinels must not count as valid physics. */
inline bool isPhysicallyPlausibleScalar(Real v, Real abs_cap = Real(1.0e12))
{
    return std::isfinite(v) && std::abs(v) < abs_cap;
}

struct ReloadDpAudit
{
    bool ok = false;
    Real dp_from_median_cbrt_vol = Real(-1);
    Real vol_median = Real(-1);
    size_t n_samples = 0;
    size_t n_rotor = 0;
    size_t n_glass = 0;
    size_t n_wall = 0;
    std::string fail_reason;
};

inline ReloadDpAudit inferReloadDpFromXml(const std::string &reload_dir, size_t max_samples = 4000)
{
    ReloadDpAudit a;
    const fs::path xml = fs::path(reload_dir) / "Reload.xml";
    std::ifstream in(xml);
    if (!in)
    {
        a.fail_reason = "cannot_open_Reload.xml";
        return a;
    }
    std::string line;
    std::string body;
    std::vector<Real> vols;
    vols.reserve(max_samples);
    auto flush_body_count = [&]() {
        if (body == "Rotor")
            a.n_rotor += 0; // counted per particle below
    };
    (void)flush_body_count;
    while (std::getline(in, line))
    {
        const auto bn = line.find("body name=\"");
        if (bn != std::string::npos)
        {
            const auto s = bn + 11;
            const auto e = line.find('"', s);
            if (e != std::string::npos)
                body = line.substr(s, e - s);
            continue;
        }
        const auto vm = line.find("VolumetricMeasure=\"");
        if (vm == std::string::npos)
            continue;
        if (body == "Rotor")
            ++a.n_rotor;
        else if (body == "GlassBody")
            ++a.n_glass;
        else if (body == "WallBoundary")
            ++a.n_wall;
        if (vols.size() >= max_samples)
            continue;
        const auto s = vm + 19;
        const auto e = line.find('"', s);
        if (e == std::string::npos)
            continue;
        const Real v = static_cast<Real>(std::atof(line.substr(s, e - s).c_str()));
        if (std::isfinite(v) && v > TinyReal)
            vols.push_back(v);
    }
    a.n_samples = vols.size();
    if (vols.empty())
    {
        a.fail_reason = "no_positive_VolumetricMeasure";
        return a;
    }
    std::nth_element(vols.begin(), vols.begin() + static_cast<std::ptrdiff_t>(vols.size() / 2), vols.end());
    a.vol_median = vols[vols.size() / 2];
    a.dp_from_median_cbrt_vol = std::cbrt(a.vol_median);
    a.ok = std::isfinite(a.dp_from_median_cbrt_vol) && a.dp_from_median_cbrt_vol > TinyReal;
    if (!a.ok)
        a.fail_reason = "invalid_cbrt_vol";
    return a;
}

struct FieldValidityAudit
{
    bool ok = false;
    size_t n_checked = 0;
    size_t n_nonfinite = 0;
    size_t n_sentinel_abs = 0;
    size_t n_nonpositive_density = 0;
    size_t n_nonpositive_vol = 0;
    Real u_max = Real(-1);
    Real rho_min = Real(-1);
    Real rho_max = Real(-1);
    std::string fail_reason;
};

inline FieldValidityAudit auditFluidFields(BaseParticles &glass)
{
    FieldValidityAudit a;
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Vecd>(glass, "Velocity");
    syncVariableToHost<Real>(glass, "Density");
    syncVariableToHost<Real>(glass, "VolumetricMeasure");
    const size_t n = glass.TotalRealParticles();
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = glass.getVariableDataByName<Vecd>("Velocity");
    const Real *rho = glass.getVariableDataByName<Real>("Density");
    const Real *vol = glass.getVariableDataByName<Real>("VolumetricMeasure");
    a.n_checked = n;
    Real umax = 0;
    Real rmin = MaxReal;
    Real rmax = -MaxReal;
    for (size_t i = 0; i < n; ++i)
    {
        const bool pos_ok = isPhysicallyPlausibleScalar(pos[i][0]) && isPhysicallyPlausibleScalar(pos[i][1]) &&
                            isPhysicallyPlausibleScalar(pos[i][2]);
        const bool vel_ok = isPhysicallyPlausibleScalar(vel[i][0]) && isPhysicallyPlausibleScalar(vel[i][1]) &&
                            isPhysicallyPlausibleScalar(vel[i][2]);
        if (!pos_ok || !vel_ok)
        {
            if (!std::isfinite(pos[i][0]) || !std::isfinite(vel[i][0]))
                ++a.n_nonfinite;
            else
                ++a.n_sentinel_abs;
            continue;
        }
        umax = std::max(umax, vel[i].norm());
        if (!isPhysicallyPlausibleScalar(rho[i]) || rho[i] <= TinyReal)
        {
            ++a.n_nonpositive_density;
            continue;
        }
        if (!isPhysicallyPlausibleScalar(vol[i]) || vol[i] <= TinyReal)
        {
            ++a.n_nonpositive_vol;
            continue;
        }
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
    }
    a.u_max = umax;
    a.rho_min = (rmin < MaxReal * Real(0.5)) ? rmin : Real(-1);
    a.rho_max = (rmax > -MaxReal * Real(0.5)) ? rmax : Real(-1);
    if (n == 0)
        a.fail_reason = "zero_fluid_particles";
    else if (a.n_nonfinite > 0)
        a.fail_reason = "nonfinite_position_or_velocity";
    else if (a.n_sentinel_abs > 0)
        a.fail_reason = "sentinel_abs_cap_exceeded";
    else if (a.n_nonpositive_density > 0)
        a.fail_reason = "nonpositive_density";
    else if (a.n_nonpositive_vol > 0)
        a.fail_reason = "nonpositive_volume";
    else if (!(a.u_max >= Real(0) && isPhysicallyPlausibleScalar(a.u_max, Real(1.0e3))))
        a.fail_reason = "invalid_u_max";
    else if (!(a.rho_min > TinyReal && a.rho_min <= a.rho_max))
        a.fail_reason = "invalid_density_range";
    else
        a.ok = true;
    return a;
}

struct F0Cli
{
    std::string reload_dir = resolveDefaultReloadFolder();
    std::string rotor_stl = "./data/client_stir_paddle_solver_m.stl";
    Real end_time = Real(0.05);
    Real c0 = Real(8);
    Real omega_ramp_s = Real(0.05);
    Real clamp_top_band_m = Real(0.05);
    Real state_record_interval = Real(0.01);
    int screen_every = 20;
    int csv_every = 20;
    Real max_wall_hours = Real(2.0);
    bool state_recording = true;
    /** If >=0, force this system dp (metres). Else infer from Reload Vol^(1/3). */
    Real dp_override_m = Real(-1);
    /** Relative tolerance between system dp and reload-inferred dp (Doc241). */
    Real dp_match_tol = Real(0.15);
    bool allow_dp_mismatch = false;
    /** Doc242 B1: scale commanded omega (0 = parked paddle). */
    Real omega_scale = Real(1);
    /** Doc242 B2: omit fluid→solid pressure/viscous load (solid/drive diagnostic only). */
    bool omit_fluid_load_on_solid = false;
    /** Absolute overshoot/undershoot allowed when declaring reached_end_time (seconds). */
    Real end_time_abs_tol_s = Real(1.0e-4);
    /**
     * Doc242: apply WCSPH hydrostatic IC (ρ,p) from z_fs=max(z) before time marching.
     * Default ON — t=0 VTP showed P≡0 with gravity ON → free-surface spike.
     */
    bool hydrostatic_ic = true;
    /**
     * Doc242: after each dens-sum advection, re-seed hydrostatic ρ/p (diagnostic hold).
     * B0' showed IC kills early spike but dens-sum/acoustics re-break FS by ~0.011s.
     * Keep DIAGNOSTIC only — not a production fix (Doc243).
     */
    bool hold_hydrostatic_ic = false;
    /**
     * Doc243 experiment (FAILED as default): skip runtime DensitySummation under hydrostatic_ic.
     * Continuity-only ρ drifted worse (rho_rel 0.56→0.74). Kept as opt-in A/B only.
     * Production path remains dens-sum; true fix = background_hydro / compatible dens-sum (OPEN).
     */
    bool skip_runtime_dens_sum = false;
    /**
     * Doc244 S2: French free-surface p' hydro (not closed-box P2b copy).
     * BackgroundPressure p_h=ρ0 g (z_fs−z), Density=ρ0, PerturbationMomentum=1,
     * AcousticDensMode=1 with a_ref=0, omit GravityForceCK (g in ∇p_h).
     * Addresses S1 bottom-wall first-step Force kick from total-p + g double-count.
     */
    bool doc244_pprime_hydro = false;
    /**
     * Doc244 S1: treat rotor as a fully fixed contact wall (no elastic integrate, vel≡0).
     * Not an elastic-FSI validation path.
     */
    bool diag_fixed_paddle = false;
    /**
     * Doc244 S3: after solid subcycle, zero rotor Velocity/AverageVelocity (keep Position).
     * Localizes fluid kick from wall-velocity BC vs geometry/deformation. Diagnostic only.
     */
    bool diag_parked_zero_velocity = false;
    /**
     * Doc244 S3: publish commanded rigid-rotation AverageVelocity for fluid Wall/FSI
     * (ω×R(θ)(x0−c)), keep elastic Position/Velocity. Suppresses free-particle chatter
     * on the fluid BC channel diagnosed by parked-zero-velocity.
     */
    bool doc244_rigid_ave_velocity = false;
    /**
     * Doc244 S3: after solid, snap all rotor Position/Velocity/AverageVelocity to rigid
     * drive field (full-body clamp kinematics). Stronger than rigid-ave alone when
     * deformed pose still kicks fluid. Elastic integrate still runs inside the subcycle.
     */
    bool doc244_rigid_rotor_pose = false;
    /**
     * Doc244 S3 / French Doc222 E2a-ish: with pprime, freeze Density→ρ0 only below
     * z_fs−band; keep dens-sum Density in free-surface band (vel-zero was whack-a-mole).
     * FAILED as default (first-step u~5, rho_rel~0.29) — keep as A/B only.
     */
    bool doc244_pprime_fs_keep_density = false;
    Real doc244_pprime_fs_band_m = Real(0.05);
    /**
     * Doc244 S3 / French atmospheric FS with p': in z>=z_fs−band, set
     * Density so p_eos=−p_h ⇒ p_total≈0; bulk still Density=ρ0.
     * FAILED (oscillatory u~8) — A/B only.
     */
    bool doc244_pprime_fs_atmosphere = false;
    /**
     * Doc244 S3 / French Doc105 P2r: after AcousticStep1stHalf, scale pressure Force
     * by ρ0/ρ_sum on free-surface-band incomplete particles (ρ_sum<ρ0) and correct vel.
     * Not a density rewrite. Opt-in. REJECTED (≈S3b3 even with dens-sum refresh).
     */
    bool doc244_fs_pack_scale_force = false;
    Real doc244_fs_pack_band_m = Real(0.05);
    /**
     * Doc244 S3: umax-trip showed hot FS particles have dens_sum≫ρ0 (clump), not incomplete.
     * With global Density≡ρ0, EOS cannot repel clumps. Opt-in: Density=dens_sum where
     * dens_sum > ratio·ρ0 (and z>=z_cut if band>0), else Density=ρ0.
     * ratio=1 / 1.5 FAILED (wall-elevated dens_sum → first-step u~5). Prefer ratio=2 + FS band.
     */
    bool doc244_pprime_fs_overpack_density = false;
    Real doc244_pprime_fs_overpack_ratio = Real(2.0);
    /** If >0 with overpack: only apply above z_fs−band. 0 = whole domain (FAILED for ratio≤1.5). */
    Real doc244_pprime_fs_overpack_band_m = Real(0.05);
    /** Dump glass particle CSV at S0 / after first acoustic 1st+2nd half (Doc244 first-step). */
    std::string diag_first_step_csv_dir;
    bool diag_first_step_exit = false;
    /**
     * Doc244 S3: first time acoustic u_max exceeds threshold, dump full glass CSV once
     * (locate 1–4 m/s spike particles). Empty dir = off.
     */
    std::string diag_umax_csv_dir;
    Real diag_umax_threshold = Real(0.5);
    /**
     * Doc244 S3: after each acoustic step, zero glass Velocity in free-surface band
     * z >= z_fs − band (default 3·dp). Confirms umax-trip localization (FS near-axis).
     * Diagnostic only — not production.
     */
    bool diag_fs_band_vel_zero = false;
    Real diag_fs_band_depth_m = Real(-1); // <0 ⇒ 3·dp
    std::string out_dir = ".";
    std::string summary_json = "doc228_f0_summary.json";
};

inline Real hostGlassZMax(BaseParticles &glass)
{
    syncVariableToHost<Vecd>(glass, "Position");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const size_t n = glass.TotalRealParticles();
    Real zmax = -MaxReal;
    for (size_t i = 0; i < n; ++i)
        zmax = std::max(zmax, pos[i][2]);
    return zmax;
}

/** Doc244: glass particle snapshot for first-acoustic-step stages (OriginalID-stable). */
inline void hostDumpGlassFirstStepCsv(BaseParticles &glass, const std::string &path, Real rho0, Real c0,
                                      const char *stage_tag)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Vecd>(glass, "Velocity");
    syncVariableToHost<Vecd>(glass, "Force");
    syncVariableToHost<Vecd>(glass, "ForcePrior");
    syncVariableToHost<Real>(glass, "Density");
    syncVariableToHost<Real>(glass, "Mass");
    syncVariableToHost<Real>(glass, "Pressure");
    syncVariableToHost<UnsignedInt>(glass, "OriginalID");
    syncVariableToHost<Real>(glass, "DensitySummation");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = glass.getVariableDataByName<Vecd>("Velocity");
    const Vecd *f = glass.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = glass.getVariableDataByName<Vecd>("ForcePrior");
    const Real *rho = glass.getVariableDataByName<Real>("Density");
    const Real *mass = glass.getVariableDataByName<Real>("Mass");
    const Real *pres = glass.getVariableDataByName<Real>("Pressure");
    const UnsignedInt *oid = glass.getVariableDataByName<UnsignedInt>("OriginalID");
    const Real *rho_sum = glass.getVariableDataByName<Real>("DensitySummation");
    const size_t n = glass.TotalRealParticles();
    const Real c2 = c0 * c0;
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "stage,OriginalID,x,y,z,rho,rho_sum,p_field,p_eos,mass,ux,uy,uz,u,"
           "Fx,Fy,Fz,FPx,FPy,FPz,ax,ay,az,a\n";
    for (size_t i = 0; i < n; ++i)
    {
        const Real m = mass[i] + TinyReal;
        const Real u = vel[i].norm();
        const Vecd a = (f[i] + fp[i]) / m;
        out << stage_tag << "," << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << ","
            << rho[i] << "," << rho_sum[i] << "," << pres[i] << "," << (c2 * (rho[i] - rho0)) << "," << mass[i]
            << "," << vel[i][0] << "," << vel[i][1] << "," << vel[i][2] << "," << u << "," << f[i][0] << ","
            << f[i][1] << "," << f[i][2] << "," << fp[i][0] << "," << fp[i][1] << "," << fp[i][2] << ","
            << a[0] << "," << a[1] << "," << a[2] << "," << a.norm() << "\n";
    }
}

inline void hostFreezeRotorAsWall(BaseParticles &rotor, const StdVec<Vecd> &pos0)
{
    syncVariableToHost<Vecd>(rotor, "Position");
    syncVariableToHost<Vecd>(rotor, "Velocity");
    syncVariableToHost<Vecd>(rotor, "AverageVelocity");
    Vecd *pos = rotor.getVariableDataByName<Vecd>("Position");
    Vecd *vel = rotor.getVariableDataByName<Vecd>("Velocity");
    Vecd *avg = rotor.getVariableDataByName<Vecd>("AverageVelocity");
    const size_t n = rotor.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = pos0[i];
        vel[i] = Vecd::Zero();
        avg[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(rotor, "Position");
    syncVariableToDevice<Vecd>(rotor, "Velocity");
    syncVariableToDevice<Vecd>(rotor, "AverageVelocity");
}

/**
 * Doc105 P2r / Doc244: Shepard packing scale of pressure Force on FS-band incomplete
 * particles (ρ_sum<ρ0 and z>z_cut). Force *= ρ0/ρ_sum; correct the 1st-half vel kick.
 */
inline size_t hostScaleIncompleteFsForceByPacking(BaseParticles &particles, Real rho0, Real dt, Real z_cut)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    size_t n_scale = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal < rho0 && rho_sum[i] > Real(0.2) * rho0 && pos[i][2] > z_cut)
        {
            const Real s = rho0 / (rho_sum[i] + TinyReal);
            const Real inv_m = Real(1) / (mass[i] + TinyReal);
            vel[i] += f[i] * (s - Real(1)) * inv_m * dt;
            f[i] *= s;
            ++n_scale;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "Force");
    return n_scale;
}

/** Zero glass Velocity for particles with z >= z_cut (free-surface band diagnostic). */
inline void hostZeroGlassVelocityAboveZ(BaseParticles &glass, Real z_cut)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Vecd>(glass, "Velocity");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    Vecd *vel = glass.getVariableDataByName<Vecd>("Velocity");
    const size_t n = glass.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] >= z_cut)
            vel[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(glass, "Velocity");
}

/** Keep deformed Position; zero Velocity + AverageVelocity seen by fluid Wall/FSI. */
inline void hostZeroRotorVelocityKeepPosition(BaseParticles &rotor)
{
    syncVariableToHost<Vecd>(rotor, "Velocity");
    syncVariableToHost<Vecd>(rotor, "AverageVelocity");
    Vecd *vel = rotor.getVariableDataByName<Vecd>("Velocity");
    Vecd *avg = rotor.getVariableDataByName<Vecd>("AverageVelocity");
    const size_t n = rotor.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        vel[i] = Vecd::Zero();
        avg[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(rotor, "Velocity");
    syncVariableToDevice<Vecd>(rotor, "AverageVelocity");
}

/** WCSPH hydrostatic: ρ=ρ0(1+g·depth/c²), p=c²(ρ−ρ0); free surface ASSUMED at z_fs. */
inline void hostApplyHydrostaticIc(BaseParticles &glass, Real rho0, Real c0, Real gravity_g, Real z_fs)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Real>(glass, "Density");
    syncVariableToHost<Real>(glass, "Pressure");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    Real *rho = glass.getVariableDataByName<Real>("Density");
    Real *pres = glass.getVariableDataByName<Real>("Pressure");
    const size_t n = glass.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    const Real c2 = c0 * c0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real depth = std::max(Real(0), z_fs - pos[i][2]);
        rho[i] = rho0 * (Real(1) + gravity_g * depth * inv_c2);
        pres[i] = c2 * (rho[i] - rho0);
    }
    syncVariableToDevice<Real>(glass, "Density");
    syncVariableToDevice<Real>(glass, "Pressure");
}

/** Doc244 S2 / French doc220: Density→ρ0 so EOS p'≈0; hydro lives in BackgroundPressure. */
inline void hostFreezeDensityToRho0(BaseParticles &glass, Real rho0)
{
    syncVariableToHost<Real>(glass, "Density");
    Real *rho = glass.getVariableDataByName<Real>("Density");
    const size_t n = glass.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        rho[i] = rho0;
    syncVariableToDevice<Real>(glass, "Density");
}

/** Freeze Density→ρ0 only for z < z_cut; leave free-surface band for dens-sum/EOS. */
inline void hostFreezeDensityToRho0BelowZ(BaseParticles &glass, Real rho0, Real z_cut)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Real>(glass, "Density");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    Real *rho = glass.getVariableDataByName<Real>("Density");
    const size_t n = glass.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] < z_cut)
            rho[i] = rho0;
    }
    syncVariableToDevice<Real>(glass, "Density");
}

/**
 * Doc244: Density=dens_sum where dens_sum > ratio·ρ0 and (band off or z>=z_cut);
 * else Density=ρ0. Requires fresh DensitySummation. Returns number written as overpack.
 */
inline size_t hostApplyPprimeOverpackDensity(BaseParticles &glass, Real rho0, Real ratio, Real z_cut,
                                             bool use_z_cut)
{
    syncVariableToHost<Real>(glass, "Density");
    syncVariableToHost<Real>(glass, "DensitySummation");
    syncVariableToHost<Vecd>(glass, "Position");
    Real *rho = glass.getVariableDataByName<Real>("Density");
    const Real *rho_sum = glass.getVariableDataByName<Real>("DensitySummation");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const size_t n = glass.TotalRealParticles();
    const Real thresh = ratio * rho0;
    size_t n_over = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const bool in_band = (!use_z_cut) || (pos[i][2] >= z_cut);
        if (in_band && rho_sum[i] > thresh)
        {
            rho[i] = rho_sum[i];
            ++n_over;
        }
        else
            rho[i] = rho0;
    }
    syncVariableToDevice<Real>(glass, "Density");
    return n_over;
}

/**
 * Bulk: Density=ρ0. FS band z>=z_cut: Density = ρ0 − p_h/c² so p_eos≈−p_h, p_total≈0.
 * Requires BackgroundPressure already written.
 */
inline void hostApplyPprimeFsAtmosphereDensity(BaseParticles &glass, Real rho0, Real c0, Real z_cut)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Real>(glass, "Density");
    syncVariableToHost<Real>(glass, "BackgroundPressure");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    Real *rho = glass.getVariableDataByName<Real>("Density");
    const Real *p_bg = glass.getVariableDataByName<Real>("BackgroundPressure");
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    const size_t n = glass.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] < z_cut)
            rho[i] = rho0;
        else
            rho[i] = std::max(Real(0.1) * rho0, rho0 - p_bg[i] * inv_c2);
    }
    syncVariableToDevice<Real>(glass, "Density");
}

/** Free-surface p_h = ρ0 g (z_fs − z), signed depth (French Doc222 E1 default). */
inline void hostApplyBackgroundHydroPressure(BaseParticles &glass, Real rho0, Real gravity_g, Real z_fs)
{
    syncVariableToHost<Vecd>(glass, "Position");
    syncVariableToHost<Real>(glass, "BackgroundPressure");
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    Real *p_bg = glass.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = glass.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        p_bg[i] = rho0 * gravity_g * (z_fs - pos[i][2]);
    syncVariableToDevice<Real>(glass, "BackgroundPressure");
}

inline void hostSetParticleRealFlag(BaseParticles &particles, const std::string &name, Real value)
{
    syncVariableToHost<Real>(particles, name);
    Real *v = particles.getVariableDataByName<Real>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = value;
    syncVariableToDevice<Real>(particles, name);
}

inline void hostSetReferenceBodyAcceleration(BaseParticles &particles, const Vecd &a_ref)
{
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    Vecd *a = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        a[i] = a_ref;
    syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
}

inline void hostZeroVecField(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Vecd>(particles, name);
    Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Vecd::Zero();
    syncVariableToDevice<Vecd>(particles, name);
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

/** Fluid Wall reads AverageVelocity: publish rigid drive field from (θ,ω,pos0). */
inline void hostPublishRigidAverageVelocity(BaseParticles &rotor, const StdVec<Vecd> &pos0,
                                            const Vecd &center, const Vecd &axis, Real theta, Real omega)
{
    syncVariableToHost<Vecd>(rotor, "AverageVelocity");
    Vecd *avg = rotor.getVariableDataByName<Vecd>("AverageVelocity");
    const Vecd e = axis.normalized();
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    const size_t n = rotor.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd rotated = rotateAboutAxis(pos0[i] - center, e, c, s);
        avg[i] = e.cross(rotated) * omega;
    }
    syncVariableToDevice<Vecd>(rotor, "AverageVelocity");
}

/** Snap rotor kinematics to rigid drive; keep NormalDirection consistent with θ. */
inline void hostPublishRigidRotorPose(BaseParticles &rotor, const StdVec<Vecd> &pos0,
                                      const StdVec<Vecd> &n0, const Vecd &center, const Vecd &axis,
                                      Real theta, Real omega)
{
    syncVariableToHost<Vecd>(rotor, "Position");
    syncVariableToHost<Vecd>(rotor, "Velocity");
    syncVariableToHost<Vecd>(rotor, "AverageVelocity");
    syncVariableToHost<Vecd>(rotor, "NormalDirection");
    Vecd *pos = rotor.getVariableDataByName<Vecd>("Position");
    Vecd *vel = rotor.getVariableDataByName<Vecd>("Velocity");
    Vecd *avg = rotor.getVariableDataByName<Vecd>("AverageVelocity");
    Vecd *n_dir = rotor.getVariableDataByName<Vecd>("NormalDirection");
    const Vecd e = axis.normalized();
    const Real c = std::cos(theta);
    const Real s = std::sin(theta);
    const size_t n = rotor.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd rotated = rotateAboutAxis(pos0[i] - center, e, c, s);
        pos[i] = center + rotated;
        vel[i] = e.cross(rotated) * omega;
        avg[i] = vel[i];
        n_dir[i] = rotateAboutAxis(n0[i], e, c, s);
    }
    syncVariableToDevice<Vecd>(rotor, "Position");
    syncVariableToDevice<Vecd>(rotor, "Velocity");
    syncVariableToDevice<Vecd>(rotor, "AverageVelocity");
    syncVariableToDevice<Vecd>(rotor, "NormalDirection");
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
    syncVariableToHost<Real>(glass, "Temperature");
    const size_t n = glass.TotalRealParticles();
    const Vecd *pos = glass.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = glass.getVariableDataByName<Vecd>("Velocity");
    const Real *pres = glass.getVariableDataByName<Real>("Pressure");
    const Real *temp = glass.getVariableDataByName<Real>("Temperature");
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

void applyCli(int ac, char *av[], F0Cli &cli, StdVec<char *> &sph_av)
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
        else if (std::strncmp(arg, "--dp=", 5) == 0)
            cli.dp_override_m = static_cast<Real>(std::atof(arg + 5));
        else if (std::strncmp(arg, "--dp-match-tol=", 15) == 0)
            cli.dp_match_tol = static_cast<Real>(std::atof(arg + 15));
        else if (std::strcmp(arg, "--allow-dp-mismatch") == 0)
            cli.allow_dp_mismatch = true;
        else if (std::strcmp(arg, "--no-state-recording") == 0)
            cli.state_recording = false;
        else if (std::strncmp(arg, "--omega-scale=", 14) == 0)
            cli.omega_scale = static_cast<Real>(std::atof(arg + 14));
        else if (std::strcmp(arg, "--diag-omit-fluid-load-on-solid") == 0)
            cli.omit_fluid_load_on_solid = true;
        else if (std::strncmp(arg, "--end-time-abs-tol=", 19) == 0)
            cli.end_time_abs_tol_s = static_cast<Real>(std::atof(arg + 19));
        else if (std::strcmp(arg, "--no-hydrostatic-ic") == 0)
            cli.hydrostatic_ic = false;
        else if (std::strcmp(arg, "--hydrostatic-ic") == 0)
            cli.hydrostatic_ic = true;
        else if (std::strcmp(arg, "--hold-hydrostatic-ic") == 0)
            cli.hold_hydrostatic_ic = true;
        else if (std::strcmp(arg, "--skip-runtime-dens-sum") == 0)
            cli.skip_runtime_dens_sum = true;
        else if (std::strcmp(arg, "--diag-fixed-paddle") == 0)
            cli.diag_fixed_paddle = true;
        else if (std::strcmp(arg, "--diag-parked-zero-velocity") == 0)
            cli.diag_parked_zero_velocity = true;
        else if (std::strcmp(arg, "--doc244-rigid-ave-velocity") == 0)
            cli.doc244_rigid_ave_velocity = true;
        else if (std::strcmp(arg, "--doc244-rigid-rotor-pose") == 0)
            cli.doc244_rigid_rotor_pose = true;
        else if (std::strncmp(arg, "--diag-first-step-csv-dir=", 26) == 0)
            cli.diag_first_step_csv_dir = arg + 26;
        else if (std::strcmp(arg, "--diag-first-step-exit") == 0)
            cli.diag_first_step_exit = true;
        else if (std::strncmp(arg, "--diag-umax-csv-dir=", 20) == 0)
            cli.diag_umax_csv_dir = arg + 20;
        else if (std::strncmp(arg, "--diag-umax-threshold=", 22) == 0)
            cli.diag_umax_threshold = static_cast<Real>(std::atof(arg + 22));
        else if (std::strcmp(arg, "--diag-fs-band-vel-zero") == 0)
            cli.diag_fs_band_vel_zero = true;
        else if (std::strncmp(arg, "--diag-fs-band-depth=", 21) == 0)
            cli.diag_fs_band_depth_m = static_cast<Real>(std::atof(arg + 21));
        else if (std::strcmp(arg, "--doc244-pprime-hydro") == 0)
            cli.doc244_pprime_hydro = true;
        else if (std::strcmp(arg, "--doc244-pprime-fs-keep-density") == 0)
            cli.doc244_pprime_fs_keep_density = true;
        else if (std::strncmp(arg, "--doc244-pprime-fs-band=", 24) == 0)
            cli.doc244_pprime_fs_band_m = static_cast<Real>(std::atof(arg + 24));
        else if (std::strcmp(arg, "--doc244-pprime-fs-atmosphere") == 0)
            cli.doc244_pprime_fs_atmosphere = true;
        else if (std::strcmp(arg, "--doc244-fs-pack-scale-force") == 0)
            cli.doc244_fs_pack_scale_force = true;
        else if (std::strncmp(arg, "--doc244-fs-pack-band=", 22) == 0)
            cli.doc244_fs_pack_band_m = static_cast<Real>(std::atof(arg + 22));
        else if (std::strcmp(arg, "--doc244-pprime-fs-overpack-density") == 0)
            cli.doc244_pprime_fs_overpack_density = true;
        else if (std::strncmp(arg, "--doc244-pprime-fs-overpack-ratio=", 34) == 0)
            cli.doc244_pprime_fs_overpack_ratio = static_cast<Real>(std::atof(arg + 34));
        else if (std::strncmp(arg, "--doc244-pprime-fs-overpack-band=", 33) == 0)
            cli.doc244_pprime_fs_overpack_band_m = static_cast<Real>(std::atof(arg + 33));
        else if (std::strcmp(arg, "--help") == 0 || std::strcmp(arg, "-h") == 0)
        {
            std::cout << "Doc242/243/244 F0 isothermal elastic FSI:\n"
                         "  --reload-dir= --rotor-stl= --end-time= --c0= --omega-ramp=\n"
                         "  --clamp-top-band= --out-dir= --summary-json= --max-wall-hours=\n"
                         "  --dp=  (optional; else inferred from Reload Vol^(1/3))\n"
                         "  --dp-match-tol=0.15  --allow-dp-mismatch\n"
                         "  --omega-scale=1  (0 = B1 parked paddle; still elastic)\n"
                         "  --diag-fixed-paddle  (Doc244: rotor frozen as wall)\n"
                         "  --diag-parked-zero-velocity  (S3: zero rotor vel after solid; keep Position)\n"
                         "  --doc244-rigid-ave-velocity  (S3: fluid Wall uses rigid ω×r AveVel)\n"
                         "  --doc244-rigid-rotor-pose  (S3: snap rotor pose/vel to rigid drive)\n"
                         "  --diag-first-step-csv-dir=  --diag-first-step-exit\n"
                         "  --diag-umax-csv-dir=  --diag-umax-threshold=0.5  (dump once when u_max trips)\n"
                         "  --diag-fs-band-vel-zero  [--diag-fs-band-depth=]  (S3 FS band vel=0)\n"
                         "  --doc244-pprime-hydro  (S2: p'+BackgroundPressure, omit GravityForceCK)\n"
                         "  --doc244-pprime-fs-keep-density [--doc244-pprime-fs-band=0.05]  (A/B FAILED)\n"
                         "  --doc244-pprime-fs-atmosphere [--doc244-pprime-fs-band=0.05]  (A/B FAILED)\n"
                         "  --doc244-fs-pack-scale-force [--doc244-fs-pack-band=0.05]  (REJECTED)\n"
                         "  --doc244-pprime-fs-overpack-density\n"
                         "     [--doc244-pprime-fs-overpack-ratio=2] [--doc244-pprime-fs-overpack-band=0.05]\n"
                         "     (A/B: Density=dens_sum if dens_sum>ratio*rho0 in FS band; ratio≤1.5 FAILED)\n"
                         "  --diag-omit-fluid-load-on-solid  (B2 solid/drive only)\n"
                         "  --hydrostatic-ic / --no-hydrostatic-ic  (default ON)\n"
                         "  --hold-hydrostatic-ic  (DIAGNOSTIC: re-seed rho/p after dens-sum)\n"
                         "  --skip-runtime-dens-sum  (A/B only; Doc243 FAILED as default)\n"
                         "  --end-time-abs-tol=1e-4\n"
                         "  --no-state-recording\n";
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
    std::cerr << "test_3d_ophelie_client_stirring_f0_elastic_fsi requires SYCL build.\n";
    return 1;
#else
    F0Cli cli;
    StdVec<char *> sph_av;
    applyCli(ac, av, cli, sph_av);
    fs::create_directories(cli.out_dir);
    fs::create_directories("./output");

    if (!reloadXmlExists(cli.reload_dir))
    {
        std::cerr << "[ophelie][client-f0] FATAL: missing " << cli.reload_dir
                  << "/Reload.xml. Run test_3d_ophelie_client_stirring_glass_relax --bodies=all --dp=0.005.\n";
        return 2;
    }
    if (!fs::exists(cli.rotor_stl))
    {
        std::cerr << "[ophelie][client-f0] FATAL: missing rotor STL \"" << cli.rotor_stl << "\"\n";
        return 2;
    }

    const ReloadDpAudit reload_dp = inferReloadDpFromXml(cli.reload_dir);
    if (!reload_dp.ok)
    {
        std::cerr << "[ophelie][client-f0] FATAL: reload dp audit failed: " << reload_dp.fail_reason << "\n";
        return 2;
    }

    OphelieClientStirringCaseParams g_case;
    g_case.rotor_stl_solver_m = cli.rotor_stl;
    // Doc241: never silently force 0.005. Prefer Reload-inferred spacing; optional --dp=.
    const Real dp_requested =
        (cli.dp_override_m > TinyReal) ? cli.dp_override_m : reload_dp.dp_from_median_cbrt_vol;
    g_case.client.particle_spacing_m = dp_requested;
    applyClientStirringDevelopmentDefaults(g_case);

    const Real dp = g_case.french_stirring.french.dp;
    const Real dp_rel = std::abs(dp - reload_dp.dp_from_median_cbrt_vol) /
                        std::max(reload_dp.dp_from_median_cbrt_vol, TinyReal);
    std::cout << "[ophelie][client-f0] reload_audit n_rotor=" << reload_dp.n_rotor
              << " n_glass=" << reload_dp.n_glass << " n_wall=" << reload_dp.n_wall
              << " vol_median=" << reload_dp.vol_median
              << " dp_reload≈" << reload_dp.dp_from_median_cbrt_vol
              << " dp_system=" << dp << " dp_rel_err=" << dp_rel
              << (cli.dp_override_m > TinyReal ? " (cli --dp)" : " (inferred)") << "\n";
    if (dp_rel > cli.dp_match_tol && !cli.allow_dp_mismatch)
    {
        std::cerr << "[ophelie][client-f0] FATAL: system dp vs Reload Vol^(1/3) mismatch > "
                  << cli.dp_match_tol
                  << ". Regenerate Reload with matching --dp, or pass matching --dp=/--allow-dp-mismatch.\n"
                  << "  Hint: historical PASS used ~4332 rotor particles at dp≈0.005; current coarse "
                  << "Reload has " << reload_dp.n_rotor << " rotor particles at dp≈"
                  << reload_dp.dp_from_median_cbrt_vol << ".\n";
        return 3;
    }
    if (dp_rel > cli.dp_match_tol)
    {
        std::cerr << "[ophelie][client-f0] WARN: allowing dp mismatch (--allow-dp-mismatch); "
                     "results are not a valid resolution-consistent baseline.\n";
    }

    {
        std::ofstream man(fs::path(cli.out_dir) / "reload_manifest.json");
        man << std::setprecision(10);
        man << "{\n";
        man << "  \"reload_dir\": \"" << cli.reload_dir << "\",\n";
        man << "  \"n_rotor_xml\": " << reload_dp.n_rotor << ",\n";
        man << "  \"n_glass_xml\": " << reload_dp.n_glass << ",\n";
        man << "  \"n_wall_xml\": " << reload_dp.n_wall << ",\n";
        man << "  \"vol_median\": " << reload_dp.vol_median << ",\n";
        man << "  \"dp_from_median_cbrt_vol\": " << reload_dp.dp_from_median_cbrt_vol << ",\n";
        man << "  \"dp_system\": " << dp << ",\n";
        man << "  \"dp_rel_err\": " << dp_rel << ",\n";
        man << "  \"dp_override_cli\": " << cli.dp_override_m << ",\n";
        man << "  \"geometry_resolution_consistent\": " << ((dp_rel <= cli.dp_match_tol) ? 1 : 0) << "\n";
        man << "}\n";
    }

    const OphelieFrenchGlassMaterialLaws laws = buildClientCaseMaterialLaws(g_case.client);
    const Real rho0 = laws.rho.evaluate(g_case.client.t_initial_k);
    const Real mu = laws.mu.evaluate(g_case.client.mu_freeze_t_k);
    const Real omega_max = g_case.french_stirring.rotation_speed_rad_s * cli.omega_scale;
    const Vecd axis = g_case.french_stirring.rotation_axis;
    const Vecd center = g_case.french_stirring.rotation_center;

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
    // Isothermal F0: register Temperature at T0 so the client thermocouple Observer can sample
    // the same field F1/heated runs will evolve (flat T until heat/EM is restored).
    const Real t0 = g_case.client.t_initial_k;
    glass.getBaseParticles().registerStateVariable<Real>("Temperature", t0);
    glass.getBaseParticles().addEvolvingVariable<Real>("Temperature");
    glass.getBaseParticles().addVariableToWrite<Real>("Temperature");

    // Client thermocouple tip (thermal-couple-1.STL AABB lower end → solver m); probe-only.
    const OphelieClientProbeDefinition thermo_probe_def = clientStirringThermocoupleProbe();
    const Vecd thermo_probe_pos = clientStirringThermocoupleProbePosition();
    ObserverBody thermocouple_probe(sph_system, "ThermocoupleProbe");
    thermocouple_probe.getSPHAdaptation().resetAdaptationRatios(1.15, 2.0);
    thermocouple_probe.generateParticles<ObserverParticles>(StdVec<Vecd>{thermo_probe_pos});

    const size_t n_rotor = rotor.getBaseParticles().TotalRealParticles();
    const size_t n_wall = wall.getBaseParticles().TotalRealParticles();
    const size_t n_glass = glass.getBaseParticles().TotalRealParticles();
    std::cout << "[ophelie][client-f0] reload: n_glass=" << n_glass << " n_rotor=" << n_rotor
              << " n_wall=" << n_wall << " dir=" << cli.reload_dir << "\n";
    if (n_rotor == 0 || n_wall == 0 || n_glass == 0)
    {
        std::cerr << "[ophelie][client-f0] empty body in Reload.xml\n";
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
        std::cerr << "[ophelie][client-f0] FATAL: clamp empty or whole body n_clamp=" << n_clamp
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
                                                    NoKernelCorrectionCK>(rotor_inner);
    auto &solid_2nd =
        main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(rotor_inner);
    auto &solid_acoustic_dt = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(rotor, 0.5);
    auto &update_rotor_normal = main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(rotor);

    auto &glass_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass);
    glass_state_recorder.addToWrite<Vecd>(glass, "Velocity");
    glass_state_recorder.addToWrite<Real>(glass, "Pressure");
    glass_state_recorder.addToWrite<Real>(glass, "Density");
    glass_state_recorder.addToWrite<Real>(glass, "Temperature");
    auto &rotor_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(rotor);
    rotor_state_recorder.addToWrite<Vecd>(rotor, "Velocity");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "ViscousForceFromFluid");
    rotor_state_recorder.addToWrite<Vecd>(rotor, "PressureForceFromFluid");

    // Same pattern as dambreak / heat_transfer probes: Observer samples neighboring fluid.
    auto &write_probe_temperature =
        main_methods.addObserveRecorder<Real>(thermocouple_probe_contact, "Temperature");
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

    const Real z_fs = hostGlassZMax(glass.getBaseParticles());
    const Real gravity_g = Real(9.81);
    size_t last_overpack_n = 0;
    auto apply_doc244_pprime_hydro_state = [&]() {
        // AcousticStep ctor already registered BackgroundPressure / PM / AcousticDensMode / a_ref.
        hostApplyBackgroundHydroPressure(glass.getBaseParticles(), rho0, gravity_g, z_fs);
        const Real band =
            (cli.doc244_pprime_fs_band_m > TinyReal) ? cli.doc244_pprime_fs_band_m : Real(0.05);
        const Real z_cut = z_fs - band;
        if (cli.doc244_pprime_fs_atmosphere)
            hostApplyPprimeFsAtmosphereDensity(glass.getBaseParticles(), rho0, cli.c0, z_cut);
        else if (cli.doc244_pprime_fs_keep_density)
            hostFreezeDensityToRho0BelowZ(glass.getBaseParticles(), rho0, z_cut);
        else if (cli.doc244_pprime_fs_overpack_density)
        {
            const Real ratio = (cli.doc244_pprime_fs_overpack_ratio > TinyReal)
                                   ? cli.doc244_pprime_fs_overpack_ratio
                                   : Real(2.0);
            const Real oband = cli.doc244_pprime_fs_overpack_band_m;
            const bool use_z = oband > TinyReal;
            const Real z_over = use_z ? (z_fs - oband) : Real(0);
            last_overpack_n =
                hostApplyPprimeOverpackDensity(glass.getBaseParticles(), rho0, ratio, z_over, use_z);
        }
        else
            hostFreezeDensityToRho0(glass.getBaseParticles(), rho0);
        hostSetParticleRealFlag(glass.getBaseParticles(), "PerturbationMomentum", Real(1));
        hostSetParticleRealFlag(glass.getBaseParticles(), "AcousticDensMode", Real(1));
        hostSetReferenceBodyAcceleration(glass.getBaseParticles(), Vecd::Zero());
        hostZeroVecField(glass.getBaseParticles(), "GravityForceCK");
        hostZeroVecField(glass.getBaseParticles(), "PreviousGravityForceCK");
    };
    if (cli.doc244_pprime_hydro)
    {
        apply_doc244_pprime_hydro_state();
        std::cout << "[ophelie][client-f0][doc244] pprime_hydro=1 z_fs=" << z_fs
                  << " BackgroundPressure=p_h PerturbationMomentum=1 "
                     "AcousticDensMode=1 a_ref=0 omit_GravityForceCK=1";
        if (cli.doc244_pprime_fs_atmosphere)
            std::cout << " fs_atmosphere=1 band=" << cli.doc244_pprime_fs_band_m
                      << " (FS: p_eos=-p_h; bulk Density=rho0)";
        else if (cli.doc244_pprime_fs_keep_density)
            std::cout << " fs_keep_density=1 band=" << cli.doc244_pprime_fs_band_m
                      << " (A/B; bulk Density=rho0; FS dens-sum)";
        else if (cli.doc244_pprime_fs_overpack_density)
            std::cout << " fs_overpack_density=1 ratio=" << cli.doc244_pprime_fs_overpack_ratio
                      << " band=" << cli.doc244_pprime_fs_overpack_band_m
                      << " (Density=dens_sum if >ratio*rho0 in band; else rho0)";
        else
            std::cout << " Density=rho0";
        std::cout << " (French free-surface path; not closed-box P2b)\n";
    }
    else if (cli.hydrostatic_ic)
    {
        hostApplyHydrostaticIc(glass.getBaseParticles(), rho0, cli.c0, gravity_g, z_fs);
        std::cout << "[ophelie][client-f0] hydrostatic_ic=1 z_fs=" << z_fs
                  << " ASSUMED free-surface at particle z_max; p=c0^2(rho-rho0)\n";
    }
    else
    {
        std::cout << "[ophelie][client-f0] hydrostatic_ic=0 (legacy Doc241: P starts at 0 with gravity ON)\n";
    }

    density_regularization.exec();
    // Re-assert hydro state after dens-sum (dens-sum may pull rho toward kernel sum).
    if (cli.doc244_pprime_hydro)
        apply_doc244_pprime_hydro_state();
    else if (cli.hydrostatic_ic)
        hostApplyHydrostaticIc(glass.getBaseParticles(), rho0, cli.c0, gravity_g, z_fs);
    // Doc243: default keeps runtime dens-sum. --skip-runtime-dens-sum is FAILED A/B only.
    const bool runtime_dens_sum = cli.hold_hydrostatic_ic || !cli.skip_runtime_dens_sum;
    const char *dens_sum_mode = "runtime_dens_sum";
    if (cli.hold_hydrostatic_ic)
        dens_sum_mode = "runtime_then_hold_DIAGNOSTIC";
    else if (cli.skip_runtime_dens_sum)
        dens_sum_mode = "skip_runtime_A_B_FAILED_default";
    std::cout << "[ophelie][client-f0][doc243] dens_sum_mode=" << dens_sum_mode
              << " runtime_dens_sum=" << (runtime_dens_sum ? 1 : 0)
              << " (skip-as-default FAILED: rho_rel worsened; need background_hydro)\n";
    fluid_advection_step_setup.exec();
    if (!cli.doc244_pprime_hydro)
        constant_gravity.exec();
    else
    {
        // g lives in BackgroundPressure; do not put GravityForceCK into ForcePrior.
        hostZeroVecField(glass.getBaseParticles(), "GravityForceCK");
        hostZeroVecField(glass.getBaseParticles(), "PreviousGravityForceCK");
    }
    fluid_viscous_force.exec();
    viscous_force_from_fluid_on_rotor.exec();

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording = time_stepper.addTriggerByInterval(std::max(cli.state_record_interval, Real(1.0e-6)));

    // Prefer --out-dir so smoke artifacts are not buried under cwd ./output/.
    const fs::path monitor_path = fs::path(cli.out_dir) / "client_f0_elastic_fsi_monitor.csv";
    std::ofstream monitor(monitor_path.string());
    monitor << "advection_step,physical_time_s,omega,theta,acoustic_dt_s,solid_dt_used_s,solid_dt_candidate_s,"
               "solid_substeps,solid_covered_frac,u_max_mps,r_tip_m,tip_mps,f_visc_N,f_pres_N,torque_z_Nm,"
               "elastic_offset_m,probe_T_K,probe_u_mps,probe_p_Pa,wall_clock_s\n";
    const fs::path acoustic_mon_path = fs::path(cli.out_dir) / "client_f0_acoustic_umax.csv";
    std::ofstream acoustic_mon(acoustic_mon_path.string());
    acoustic_mon << "acoustic_step,physical_time_s,u_max_mps,omega,wall_clock_s\n";

    if (cli.state_recording)
    {
        glass_state_recorder.writeToFile(0);
        rotor_state_recorder.writeToFile(0);
    }
    {
        const Real u0 = u_max_reduce.exec();
        acoustic_mon << "0,0," << u0 << ",0,0\n";
        acoustic_mon.flush();
        std::cout << "[ophelie][client-f0] IC u_max=" << u0 << " (post hydrostatic_ic="
                  << (cli.hydrostatic_ic ? 1 : 0) << ")\n";
    }
    if (cli.diag_fixed_paddle)
    {
        hostFreezeRotorAsWall(rotor.getBaseParticles(), pos0);
        rotor.setNewlyUpdated();
        std::cout << "[ophelie][client-f0][doc244] diag_fixed_paddle=1 rotor frozen as wall contact\n";
    }
    if (cli.diag_parked_zero_velocity)
        std::cout << "[ophelie][client-f0][doc244] diag_parked_zero_velocity=1 "
                     "(zero rotor vel after each solid subcycle; keep Position)\n";
    if (cli.doc244_rigid_ave_velocity)
        std::cout << "[ophelie][client-f0][doc244] rigid_ave_velocity=1 "
                     "(fluid Wall AverageVelocity = rigid ω×R(θ)(x0-c))\n";
    if (cli.doc244_rigid_rotor_pose)
        std::cout << "[ophelie][client-f0][doc244] rigid_rotor_pose=1 "
                     "(snap rotor Position/Velocity/AveVel/Normal to rigid drive after solid)\n";
    if (cli.doc244_fs_pack_scale_force)
        std::cout << "[ophelie][client-f0][doc244] fs_pack_scale_force=1 band=" << cli.doc244_fs_pack_band_m
                  << " (dens-sum refresh each ac; Force*=rho0/rho_sum on FS incomplete after 1st half)\n";
    if (cli.diag_fs_band_vel_zero)
    {
        const Real band =
            (cli.diag_fs_band_depth_m > TinyReal) ? cli.diag_fs_band_depth_m : (Real(3) * dp);
        std::cout << "[ophelie][client-f0][doc244] diag_fs_band_vel_zero=1 band=" << band
                  << " z_cut=" << (z_fs - band) << "\n";
    }
    if (!cli.diag_first_step_csv_dir.empty())
    {
        fs::create_directories(cli.diag_first_step_csv_dir);
        hostDumpGlassFirstStepCsv(glass.getBaseParticles(),
                                  (fs::path(cli.diag_first_step_csv_dir) / "S0_post_ic.csv").string(), rho0, cli.c0,
                                  "S0_post_ic");
        std::cout << "[ophelie][client-f0][doc244] wrote S0_post_ic.csv\n";
    }
    monitor << std::setprecision(10);
    monitor.flush();

    write_probe_temperature.writeToFile(0);
    write_probe_velocity.writeToFile(0);
    write_probe_pressure.writeToFile(0);

    if (cli.state_recording)
    {
        glass_state_recorder.writeToFile(0);
        rotor_state_recorder.writeToFile(0);
    }

    const Real dt_s0 = solid_acoustic_dt.exec();
    std::cout << "[ophelie][client-f0] isothermal elastic FSI start\n"
              << "  mechanics=" << clientPaddleMechanicsName(g_case.paddle_mechanics)
              << " E=" << g_case.steel.youngs_modulus_pa << " rho_s=" << g_case.steel.rho_kg_m3
              << " n_clamp=" << n_clamp << " clamp_frac=" << (Real(n_clamp) / Real(n_rotor))
              << " z_clamp=" << z_clamp << "\n"
              << "  rho_f=" << rho0 << " mu=" << mu << " c0=" << cli.c0 << " dp=" << dp
              << " omega_max=" << omega_max << " ramp=" << cli.omega_ramp_s << " s end=" << cli.end_time
              << " solid_dt0=" << dt_s0 << "\n"
              << "  thermocouple probe id=" << thermo_probe_def.probe_id
              << " stl=" << clientStirringThermocoupleStlName()
              << " status=" << thermo_probe_def.status
              << " pos_m=[" << thermo_probe_pos[0] << "," << thermo_probe_pos[1] << "," << thermo_probe_pos[2]
              << "] T0=" << t0 << " K (isothermal until F1)\n"
              << "  no EM / no Q / no Boussinesq; gravity on; clamp provisional top band\n"
              << "  note: first [adv=] line waits for first advection trigger; "
                 "acoustic heartbeats print sooner (Doc240 smoke)\n";
    std::cout.flush();

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
    bool solid_time_uncovered = false;
    Real last_omega = 0;
    Real last_theta = 0;
    Real last_elastic = 0;
    Real last_torque = 0;
    Real last_u_max = 0;
    Real last_solid_covered_frac = Real(1);
    Real last_solid_dt_candidate = dt_s0;
    RotorLoadSnap last_load{};
    bool diag_umax_dumped = false;

    // Freeze solid dt at start-of-run CFL (dt_s0). Adaptive solid_acoustic_dt collapses after
    // clamp/drive on coarse client reload (Doc240). Doc242: still report live candidate each step;
    // uncovered solid time after kSolidSubMax is FAIL (not WARN-only).
    constexpr size_t kSolidSubMax = 5000;
    struct SolidSubResult
    {
        size_t n_sub = 0;
        Real covered = 0;
        Real covered_frac = 0;
        bool capped = false;
    };
    auto run_solid_subcycle = [&](Real dt_fluid, Real physical_time, size_t ac_idx) {
        SolidSubResult r;
        if (cli.diag_fixed_paddle)
        {
            hostFreezeRotorAsWall(rotor.getBaseParticles(), pos0);
            rotor.setNewlyUpdated();
            r.covered = dt_fluid;
            r.covered_frac = Real(1);
            last_omega = 0;
            last_theta = 0;
            last_solid_covered_frac = Real(1);
            if (cli.screen_every > 0 && (ac_idx <= 3 || ac_idx % 10 == 0))
                std::cout << "[ophelie][client-f0][doc244] fixed_paddle skip solid integrate ac=" << ac_idx
                          << " dt_ac=" << dt_fluid << "\n";
            return r;
        }
        const Real dt_s = dt_s0;
        const size_t n_est =
            dt_s > TinyReal ? static_cast<size_t>(std::ceil(static_cast<double>(dt_fluid / dt_s))) : 0;
        if (cli.screen_every > 0 && (ac_idx <= 3 || ac_idx % 10 == 0))
        {
            std::cout << "[ophelie][client-f0] solid_begin ac=" << ac_idx
                      << " dt_ac=" << dt_fluid << " dt_s_frozen=" << dt_s
                      << " dt_s_candidate=" << solid_acoustic_dt.exec()
                      << " n_sub_est=" << n_est << " n_sub_max=" << kSolidSubMax
                      << " wall_s=" << elapsed_s() << "\n";
            std::cout.flush();
        }
        while (r.covered < dt_fluid - TinyReal)
        {
            Real dt_use = dt_s;
            if (r.covered + dt_use > dt_fluid)
                dt_use = dt_fluid - r.covered;
            if (dt_use <= TinyReal)
                break;
            Real omega = 0;
            Real theta = 0;
            omegaThetaRamp(physical_time + r.covered + Real(0.5) * dt_use, cli.omega_ramp_s, omega_max, omega, theta);
            solid_1st.exec(dt_use);
            applyClampSpinDrive(rotor.getBaseParticles(), clamp_ids, pos0, n0, center, axis, theta, omega, true);
            solid_2nd.exec(dt_use);
            applyClampSpinDrive(rotor.getBaseParticles(), clamp_ids, pos0, n0, center, axis, theta, omega, true);
            r.covered += dt_use;
            ++r.n_sub;
            ++solid_steps_total;
            last_omega = omega;
            last_theta = theta;
            const bool print_sub =
                cli.screen_every > 0 && (ac_idx <= 3 || ac_idx % 10 == 0) &&
                (r.n_sub == 1 || r.n_sub == n_est || r.n_sub % 200 == 0 || r.n_sub >= kSolidSubMax);
            if (print_sub)
            {
                std::cout << "[ophelie][client-f0] solid_sub ac=" << ac_idx
                          << " n_sub=" << r.n_sub << "/" << n_est
                          << " solid_steps_tot=" << solid_steps_total
                          << " wall_s=" << elapsed_s() << "\n";
                std::cout.flush();
            }
            if (r.n_sub >= kSolidSubMax)
            {
                r.capped = true;
                break;
            }
        }
        r.covered_frac = dt_fluid > TinyReal ? r.covered / dt_fluid : Real(1);
        last_solid_covered_frac = r.covered_frac;
        last_solid_dt_candidate = solid_acoustic_dt.exec();
        if (r.capped || r.covered_frac < Real(0.999))
        {
            solid_time_uncovered = true;
            std::cerr << "[ophelie][client-f0] FAIL solid subcycle uncovered "
                      << " (ac=" << ac_idx << " covered=" << r.covered << "/" << dt_fluid
                      << " frac=" << r.covered_frac << " capped=" << (r.capped ? 1 : 0)
                      << " wall_s=" << elapsed_s() << ")\n";
        }
        // Fluid Wall/FSI reads AverageVelocity (+ Position). Default: copy elastic Velocity.
        if (cli.doc244_rigid_rotor_pose)
            hostPublishRigidRotorPose(rotor.getBaseParticles(), pos0, n0, center, axis, last_theta,
                                      last_omega);
        else if (cli.doc244_rigid_ave_velocity)
            hostPublishRigidAverageVelocity(rotor.getBaseParticles(), pos0, center, axis, last_theta,
                                            last_omega);
        else
            copyVelocityToAverageVelocity(rotor.getBaseParticles());
        if (cli.diag_parked_zero_velocity)
        {
            hostZeroRotorVelocityKeepPosition(rotor.getBaseParticles());
            if (cli.screen_every > 0 && (ac_idx <= 3 || ac_idx % 10 == 0))
                std::cout << "[ophelie][client-f0][doc244] parked_zero_velocity after solid ac=" << ac_idx
                          << "\n";
        }
        rotor.setNewlyUpdated();
        return r;
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
            std::cerr << "[ophelie][client-f0] diverged at acoustic step " << acoustic_steps
                      << " t=" << time_stepper.getPhysicalTime() << " dt=" << dt << "\n";
            diverged = true;
            break;
        }

        if (cli.screen_every > 0 && (acoustic_steps <= 3 || acoustic_steps % 10 == 0))
        {
            std::cout << "[ophelie][client-f0] fluid_ac ac=" << acoustic_steps
                      << " t=" << time_stepper.getPhysicalTime() << " dt_ac=" << dt
                      << " wall_s=" << elapsed_s() << "\n";
            std::cout.flush();
        }
        // Pack-scale / overpack need a fresh DensitySummation; runtime dens-sum only on advection.
        // Refresh then re-apply p' Density (freeze or overpack) before 1st half.
        if (cli.doc244_fs_pack_scale_force || cli.doc244_pprime_fs_overpack_density)
        {
            density_regularization.exec();
            if (cli.doc244_pprime_hydro)
                apply_doc244_pprime_hydro_state();
            if (cli.doc244_pprime_fs_overpack_density && cli.screen_every > 0 &&
                (acoustic_steps <= 3 || acoustic_steps % 10 == 0))
                std::cout << "[ophelie][client-f0][doc244] fs_overpack n=" << last_overpack_n
                          << " ratio=" << cli.doc244_pprime_fs_overpack_ratio
                          << " band=" << cli.doc244_pprime_fs_overpack_band_m
                          << " ac=" << acoustic_steps << "\n";
        }
        acoustic_step_1st_half.exec(dt);
        if (cli.doc244_fs_pack_scale_force)
        {
            const Real band =
                (cli.doc244_fs_pack_band_m > TinyReal) ? cli.doc244_fs_pack_band_m : Real(0.05);
            const size_t n_sc =
                hostScaleIncompleteFsForceByPacking(glass.getBaseParticles(), rho0, dt, z_fs - band);
            if (cli.screen_every > 0 && (acoustic_steps <= 3 || acoustic_steps % 10 == 0))
                std::cout << "[ophelie][client-f0][doc244] fs_pack_scale n=" << n_sc
                          << " z_cut=" << (z_fs - band) << " ac=" << acoustic_steps << "\n";
        }
        if (acoustic_steps == 1 && !cli.diag_first_step_csv_dir.empty())
        {
            hostDumpGlassFirstStepCsv(glass.getBaseParticles(),
                                      (fs::path(cli.diag_first_step_csv_dir) / "S2S3_after_acoustic_1st_half.csv").string(),
                                      rho0, cli.c0, "S2S3_after_acoustic_1st_half");
            std::cout << "[ophelie][client-f0][doc244] wrote S2S3_after_acoustic_1st_half.csv "
                         "(inner+wall/paddle contact fused in 1st half)\n";
        }
        if (!cli.omit_fluid_load_on_solid)
            pressure_force_from_fluid_on_rotor.exec(dt);
        acoustic_step_2nd_half.exec(dt);
        // Continuity in 2nd half can drift ρ; French p' path keeps Density=ρ0.
        if (cli.doc244_pprime_hydro)
            apply_doc244_pprime_hydro_state();
        if (acoustic_steps == 1 && !cli.diag_first_step_csv_dir.empty())
        {
            hostDumpGlassFirstStepCsv(glass.getBaseParticles(),
                                      (fs::path(cli.diag_first_step_csv_dir) / "S5_after_acoustic_2nd_half.csv").string(),
                                      rho0, cli.c0, "S5_after_acoustic_2nd_half");
            std::cout << "[ophelie][client-f0][doc244] wrote S5_after_acoustic_2nd_half.csv\n";
        }

        // Solid subcycle uses time at start of this acoustic interval.
        const Real t_start = time_stepper.getPhysicalTime() - dt;
        const SolidSubResult solid_r = run_solid_subcycle(dt, t_start, acoustic_steps);
        const size_t n_sub = solid_r.n_sub;
        if (solid_time_uncovered)
        {
            diverged = true;
            break;
        }
        if (acoustic_steps == 1 && !cli.diag_first_step_csv_dir.empty())
        {
            hostDumpGlassFirstStepCsv(glass.getBaseParticles(),
                                      (fs::path(cli.diag_first_step_csv_dir) / "S_after_solid_subcycle.csv").string(),
                                      rho0, cli.c0, "S_after_solid_subcycle");
            {
                std::ofstream note(fs::path(cli.diag_first_step_csv_dir) / "STAGE_NOTES.md");
                note << "# Doc244 first acoustic step stages\n\n"
                     << "- S0: after hydrostatic IC (+ dens-sum reassert)\n"
                     << "- S2S3: after AcousticStep1stHalf (inner pressure + wall/paddle contact fused)\n"
                     << "- S5: after AcousticStep2ndHalf (continuity)\n"
                     << "- S_after_solid: after solid subcycle (fixed-paddle = freeze only)\n"
                     << "- Gravity/viscous ForcePrior updates mainly on advection steps (later)\n"
                     << "- diag_fixed_paddle=" << (cli.diag_fixed_paddle ? 1 : 0) << "\n";
            }
            if (cli.diag_first_step_exit)
            {
                std::cout << "[ophelie][client-f0][doc244] diag_first_step_exit after ac=1\n";
                break;
            }
        }
        if (cli.diag_fs_band_vel_zero)
        {
            const Real band =
                (cli.diag_fs_band_depth_m > TinyReal) ? cli.diag_fs_band_depth_m : (Real(3) * dp);
            hostZeroGlassVelocityAboveZ(glass.getBaseParticles(), z_fs - band);
            if (cli.screen_every > 0 && (acoustic_steps <= 3 || acoustic_steps % 10 == 0))
                std::cout << "[ophelie][client-f0][doc244] fs_band_vel_zero z>=" << (z_fs - band)
                          << " ac=" << acoustic_steps << "\n";
        }

        // Heartbeat before first advection: otherwise smoke looks "stuck" for minutes
        // while solid_dt~1e-6 subcycles inside acoustic steps (GPU still busy).
        if (cli.screen_every > 0 && (acoustic_steps <= 3 || acoustic_steps % 10 == 0))
        {
            std::cout << "[ophelie][client-f0] ac=" << acoustic_steps
                      << " t=" << time_stepper.getPhysicalTime() << " dt_ac=" << dt
                      << " solid_sub=" << n_sub << " solid_steps_tot=" << solid_steps_total
                      << " omega=" << last_omega << " wall_s=" << elapsed_s() << "\n";
            std::cout.flush();
        }
        // Doc242: sample u_max on acoustic steps (B1 had no advection → empty CSV).
        if (cli.csv_every > 0 && (acoustic_steps <= 5 || acoustic_steps % 5 == 0))
        {
            const Real u_ac = u_max_reduce.exec();
            acoustic_mon << acoustic_steps << "," << time_stepper.getPhysicalTime() << "," << u_ac << ","
                         << last_omega << "," << elapsed_s() << "\n";
            acoustic_mon.flush();
        }
        // Doc244: one-shot particle dump when u_max first exceeds threshold (spike localization).
        if (!diag_umax_dumped && !cli.diag_umax_csv_dir.empty())
        {
            const Real u_trip = u_max_reduce.exec();
            if (u_trip > cli.diag_umax_threshold)
            {
                // Fresh dens-sum so rho_sum in CSV matches current neighbors (advection-only otherwise).
                density_regularization.exec();
                if (cli.doc244_pprime_hydro)
                    apply_doc244_pprime_hydro_state();
                fs::create_directories(cli.diag_umax_csv_dir);
                const std::string tag = "UmaxTrip_ac" + std::to_string(acoustic_steps);
                hostDumpGlassFirstStepCsv(
                    glass.getBaseParticles(),
                    (fs::path(cli.diag_umax_csv_dir) / (tag + ".csv")).string(), rho0, cli.c0,
                    tag.c_str());
                // Log dens_sum of the hottest particle (OriginalID + rho_sum) for incomplete-support check.
                {
                    syncVariableToHost<Vecd>(glass.getBaseParticles(), "Velocity");
                    syncVariableToHost<Real>(glass.getBaseParticles(), "DensitySummation");
                    syncVariableToHost<UnsignedInt>(glass.getBaseParticles(), "OriginalID");
                    syncVariableToHost<Vecd>(glass.getBaseParticles(), "Position");
                    const Vecd *vel_h = glass.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
                    const Real *rs_h = glass.getBaseParticles().getVariableDataByName<Real>("DensitySummation");
                    const UnsignedInt *oid_h =
                        glass.getBaseParticles().getVariableDataByName<UnsignedInt>("OriginalID");
                    const Vecd *pos_h = glass.getBaseParticles().getVariableDataByName<Vecd>("Position");
                    size_t i_max = 0;
                    Real u_best = Real(0);
                    const size_t n_g = glass.getBaseParticles().TotalRealParticles();
                    for (size_t i = 0; i < n_g; ++i)
                    {
                        const Real ui = vel_h[i].norm();
                        if (ui > u_best)
                        {
                            u_best = ui;
                            i_max = i;
                        }
                    }
                    const Real rxy =
                        std::sqrt(pos_h[i_max][0] * pos_h[i_max][0] + pos_h[i_max][1] * pos_h[i_max][1]);
                    std::cout << "[ophelie][client-f0][doc244] umax_trip dump " << tag << " u_max=" << u_trip
                              << " thresh=" << cli.diag_umax_threshold << " omega=" << last_omega
                              << " t=" << time_stepper.getPhysicalTime()
                              << " hot_oid=" << oid_h[i_max] << " z=" << pos_h[i_max][2] << " rxy=" << rxy
                              << " rho_sum=" << rs_h[i_max] << " rho0=" << rho0 << "\n";
                }
                diag_umax_dumped = true;
                std::cout.flush();
            }
        }

        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            update_rotor_normal.exec();

            if (cli.state_recording && state_recording())
            {
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
            if (runtime_dens_sum)
            {
                density_regularization.exec();
                if (cli.doc244_pprime_hydro)
                    apply_doc244_pprime_hydro_state();
                else if (cli.hold_hydrostatic_ic)
                    hostApplyHydrostaticIc(glass.getBaseParticles(), rho0, cli.c0, gravity_g, z_fs);
            }
            fluid_advection_step_setup.exec();
            fluid_viscous_force.exec();
            if (!cli.doc244_pprime_hydro)
                constant_gravity.exec();
            if (!cli.omit_fluid_load_on_solid)
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
                Real probe_T = t0;
                Real probe_u = 0;
                Real probe_p = 0;
                hostSampleFluidAtProbe(glass.getBaseParticles(), thermo_probe_pos, probe_T, probe_u, probe_p);
                write_probe_temperature.writeToFile(advection_steps);
                write_probe_velocity.writeToFile(advection_steps);
                write_probe_pressure.writeToFile(advection_steps);
                const Real wall_now = static_cast<Real>(elapsed_s());
                if (write_csv)
                {
                    monitor << advection_steps << "," << t << "," << last_omega << "," << last_theta << "," << dt
                            << "," << dt_s0 << "," << last_solid_dt_candidate << "," << n_sub << ","
                            << last_solid_covered_frac << "," << last_u_max << ","
                            << last_load.r_tip << "," << last_load.tip_speed << "," << last_load.f_visc << ","
                            << last_load.f_pres << "," << last_load.torque_z << "," << last_elastic << ","
                            << probe_T << "," << probe_u << "," << probe_p << "," << wall_now << "\n";
                    monitor.flush();
                }
                if (write_screen)
                {
                    std::cout << "[ophelie][client-f0] adv=" << advection_steps << " t=" << t
                              << " dt=" << dt << " solid_sub=" << n_sub << " u_max=" << last_u_max
                              << " tip=" << last_load.tip_speed << " Tq=" << last_torque
                              << " e_off=" << last_elastic << " probe_T=" << probe_T
                              << " probe_u=" << probe_u << " wall_s=" << wall_now << "\n";
                }
                if (!std::isfinite(last_u_max) || !std::isfinite(last_elastic) || last_u_max > Real(50))
                {
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
    const FieldValidityAudit fields = auditFluidFields(glass.getBaseParticles());
    if (fields.ok && isPhysicallyPlausibleScalar(last_u_max, Real(1.0e3)) && last_u_max >= Real(0))
        last_u_max = std::max(last_u_max, fields.u_max);

    const bool scalars_ok = isPhysicallyPlausibleScalar(last_u_max, Real(1.0e3)) && last_u_max >= Real(0) &&
                            isPhysicallyPlausibleScalar(last_elastic) && last_elastic >= Real(0) &&
                            isPhysicallyPlausibleScalar(last_torque);
    const bool finite_ok = fields.ok && scalars_ok;
    const bool clamp_ok = (n_clamp > 0) && (n_clamp < n_rotor);
    const bool drive_executed = (cli.omega_scale > TinyReal) ? (last_omega > TinyReal) : true;
    const bool fluid_load_sample_valid =
        cli.omit_fluid_load_on_solid ? true
                                    : (std::isfinite(last_load.f_pres) && std::isfinite(last_load.f_visc));
    const bool fluid_contact_load_seen =
        cli.omit_fluid_load_on_solid ? false : (last_load.f_pres > TinyReal || last_load.f_visc > TinyReal);
    // Doc242: reached_end_time requires actual end (small abs tol), NOT 0.8*end_time.
    const Real t_deficit = cli.end_time - t_final;
    const bool time_ok = (!budget_exhausted) && (t_deficit <= cli.end_time_abs_tol_s);
    const bool geometry_ok = dp_rel <= cli.dp_match_tol;
    const bool reached_end = time_ok && !diverged && !budget_exhausted && !solid_time_uncovered;
    const Real tip_ref = std::max(last_load.tip_speed, std::max(std::abs(omega_max) * Real(0.05), Real(1.0e-3)));
    const Real rho_rel_max =
        (rho0 > TinyReal && fields.rho_max > TinyReal) ? (fields.rho_max / rho0 - Real(1)) : Real(-1);
    // Doc242 engineering screen (not contract): U not >> tip; density vs rho0 within ~5%.
    const bool numerical_quality_ok =
        finite_ok && (last_u_max <= Real(20) * tip_ref) && (rho_rel_max >= Real(0)) && (rho_rel_max <= Real(0.05));
    std::string fail_reason;
    if (!geometry_ok)
        fail_reason = "geometry_resolution_inconsistent";
    else if (solid_time_uncovered)
        fail_reason = "solid_time_uncovered";
    else if (diverged)
        fail_reason = "diverged";
    else if (budget_exhausted)
        fail_reason = "wall_budget_exhausted";
    else if (!finite_ok)
        fail_reason = fields.ok ? "invalid_summary_scalars" : fields.fail_reason;
    else if (!clamp_ok)
        fail_reason = "clamp_coverage_invalid";
    else if (!drive_executed)
        fail_reason = "drive_not_executed";
    else if (!fluid_load_sample_valid)
        fail_reason = "fluid_load_sample_invalid";
    else if (solid_steps_total == 0 && !cli.diag_fixed_paddle)
        fail_reason = "zero_solid_steps";
    else if (!time_ok && !cli.diag_first_step_exit)
        fail_reason = "did_not_reach_end_time";
    // workflow pass ≠ numerical quality / FSI physics validation
    // Doc244: fixed-paddle / first-step-exit are diagnostic — not ELASTIC_FSI validation.
    const bool workflow_ok = fail_reason.empty();
    const bool pass = workflow_ok; // legacy alias; quality reported separately

    std::cout << "[ophelie][client-f0] done t=" << t_final << " wall_s=" << wall_s
              << " solid_steps=" << solid_steps_total << " u_max=" << last_u_max
              << " tip=" << last_load.tip_speed << " Tq=" << last_torque << " e_off=" << last_elastic
              << " budget_exhausted=" << (budget_exhausted ? 1 : 0) << " diverged=" << (diverged ? 1 : 0)
              << " solid_time_uncovered=" << (solid_time_uncovered ? 1 : 0)
              << " finite_fields_ok=" << (finite_ok ? 1 : 0)
              << " numerical_quality_ok=" << (numerical_quality_ok ? 1 : 0)
              << " drive_executed=" << (drive_executed ? 1 : 0)
              << " fluid_contact_load_seen=" << (fluid_contact_load_seen ? 1 : 0)
              << " ELASTIC_FSI_SHORT_PASS=" << (pass ? 1 : 0);
    if (!fail_reason.empty())
        std::cout << " fail_reason=" << fail_reason;
    std::cout << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"client_stirring_f0_elastic_fsi\",\n";
        js << "  \"tag\": \"ELASTIC_FSI_SHORT\",\n";
        js << "  \"identity\": \"CLIENT_STIRRING_ELASTIC_BASELINE\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"ELASTIC_FSI_SHORT_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"workflow_completed\": " << (workflow_ok ? 1 : 0) << ",\n";
        js << "  \"reached_end_time\": " << (reached_end ? 1 : 0) << ",\n";
        js << "  \"requested_end_time_s\": " << cli.end_time << ",\n";
        js << "  \"actual_end_time_s\": " << t_final << ",\n";
        js << "  \"end_time_deficit_s\": " << t_deficit << ",\n";
        js << "  \"em_executed\": 0,\n";
        js << "  \"thermal_advanced\": 0,\n";
        js << "  \"flow_advanced\": " << (acoustic_steps > 0 ? 1 : 0) << ",\n";
        js << "  \"solid_advanced\": " << (solid_steps_total > 0 ? 1 : 0) << ",\n";
        js << "  \"drive_executed\": " << (drive_executed ? 1 : 0) << ",\n";
        js << "  \"fluid_contact_evaluated\": " << (cli.omit_fluid_load_on_solid ? 0 : 1) << ",\n";
        js << "  \"fluid_contact_load_seen\": " << (fluid_contact_load_seen ? 1 : 0) << ",\n";
        js << "  \"load_sample_valid\": " << (fluid_load_sample_valid ? 1 : 0) << ",\n";
        js << "  \"omit_fluid_load_on_solid\": " << (cli.omit_fluid_load_on_solid ? 1 : 0) << ",\n";
        js << "  \"diag_fixed_paddle\": " << (cli.diag_fixed_paddle ? 1 : 0) << ",\n";
        js << "  \"diag_parked_zero_velocity\": " << (cli.diag_parked_zero_velocity ? 1 : 0) << ",\n";
        js << "  \"doc244_rigid_ave_velocity\": " << (cli.doc244_rigid_ave_velocity ? 1 : 0) << ",\n";
        js << "  \"doc244_rigid_rotor_pose\": " << (cli.doc244_rigid_rotor_pose ? 1 : 0) << ",\n";
        js << "  \"diag_fs_band_vel_zero\": " << (cli.diag_fs_band_vel_zero ? 1 : 0) << ",\n";
        js << "  \"diag_first_step_csv_dir\": \"" << cli.diag_first_step_csv_dir << "\",\n";
        js << "  \"diag_first_step_exit\": " << (cli.diag_first_step_exit ? 1 : 0) << ",\n";
        js << "  \"doc244_pprime_hydro\": " << (cli.doc244_pprime_hydro ? 1 : 0) << ",\n";
        js << "  \"doc244_pprime_fs_keep_density\": " << (cli.doc244_pprime_fs_keep_density ? 1 : 0) << ",\n";
        js << "  \"doc244_pprime_fs_atmosphere\": " << (cli.doc244_pprime_fs_atmosphere ? 1 : 0) << ",\n";
        js << "  \"doc244_pprime_fs_overpack_density\": " << (cli.doc244_pprime_fs_overpack_density ? 1 : 0)
           << ",\n";
        js << "  \"doc244_pprime_fs_overpack_ratio\": " << cli.doc244_pprime_fs_overpack_ratio << ",\n";
        js << "  \"doc244_pprime_fs_overpack_band_m\": " << cli.doc244_pprime_fs_overpack_band_m << ",\n";
        js << "  \"doc244_fs_pack_scale_force\": " << (cli.doc244_fs_pack_scale_force ? 1 : 0) << ",\n";
        js << "  \"doc244_fs_pack_band_m\": " << cli.doc244_fs_pack_band_m << ",\n";
        js << "  \"doc244_pprime_fs_band_m\": " << cli.doc244_pprime_fs_band_m << ",\n";
        js << "  \"omega_scale\": " << cli.omega_scale << ",\n";
        js << "  \"hydrostatic_ic\": " << (cli.hydrostatic_ic ? 1 : 0) << ",\n";
        js << "  \"hold_hydrostatic_ic\": " << (cli.hold_hydrostatic_ic ? 1 : 0) << ",\n";
        js << "  \"skip_runtime_dens_sum\": " << (cli.skip_runtime_dens_sum ? 1 : 0) << ",\n";
        js << "  \"runtime_dens_sum\": " << (runtime_dens_sum ? 1 : 0) << ",\n";
        js << "  \"dens_sum_mode\": \"" << dens_sum_mode << "\",\n";
        js << "  \"z_fs_assumed_m\": " << z_fs << ",\n";
        js << "  \"french_pressure_path\": \""
           << (cli.doc244_pprime_hydro ? "p_prime_plus_BackgroundPressure_no_GravityForceCK"
                                      : "total_p_plus_GravityForceCK")
           << "\",\n";
        js << "  \"mechanics_mode\": \""
           << (cli.omit_fluid_load_on_solid ? "elastic_drive_diag_no_fluid_load" : "elastic_driven") << "\",\n";
        js << "  \"reload_dir\": \"" << cli.reload_dir << "\",\n";
        js << "  \"n_glass\": " << n_glass << ",\n";
        js << "  \"n_rotor\": " << n_rotor << ",\n";
        js << "  \"n_wall\": " << n_wall << ",\n";
        js << "  \"n_clamp\": " << n_clamp << ",\n";
        js << "  \"clamp_top_band_m\": " << cli.clamp_top_band_m << ",\n";
        js << "  \"clamp_status\": \"provisional_top_band_unconfirmed_drive_section\",\n";
        js << "  \"E_Pa\": " << g_case.steel.youngs_modulus_pa << ",\n";
        js << "  \"rho_solid\": " << g_case.steel.rho_kg_m3 << ",\n";
        js << "  \"rho_fluid\": " << rho0 << ",\n";
        js << "  \"mu_fluid\": " << mu << ",\n";
        js << "  \"c0\": " << cli.c0 << ",\n";
        js << "  \"dp_m\": " << dp << ",\n";
        js << "  \"dp_reload_inferred_m\": " << reload_dp.dp_from_median_cbrt_vol << ",\n";
        js << "  \"dp_rel_err\": " << dp_rel << ",\n";
        js << "  \"geometry_resolution_consistent\": " << (geometry_ok ? 1 : 0) << ",\n";
        js << "  \"omega_max_rad_s\": " << omega_max << ",\n";
        js << "  \"omega_ramp_s\": " << cli.omega_ramp_s << ",\n";
        js << "  \"end_time_s\": " << cli.end_time << ",\n";
        js << "  \"physical_time_s\": " << t_final << ",\n";
        js << "  \"solid_steps_total\": " << solid_steps_total << ",\n";
        js << "  \"solid_dt0_s\": " << dt_s0 << ",\n";
        js << "  \"solid_dt_candidate_last_s\": " << last_solid_dt_candidate << ",\n";
        js << "  \"solid_time_uncovered\": " << (solid_time_uncovered ? 1 : 0) << ",\n";
        js << "  \"solid_covered_frac_last\": " << last_solid_covered_frac << ",\n";
        js << "  \"u_max_mps\": " << last_u_max << ",\n";
        js << "  \"r_tip_m\": " << last_load.r_tip << ",\n";
        js << "  \"tip_speed_mps\": " << last_load.tip_speed << ",\n";
        js << "  \"f_visc_N\": " << last_load.f_visc << ",\n";
        js << "  \"f_pres_N\": " << last_load.f_pres << ",\n";
        js << "  \"torque_z_Nm\": " << last_torque << ",\n";
        js << "  \"max_elastic_offset_m\": " << last_elastic << ",\n";
        js << "  \"finite_ok\": " << (finite_ok ? 1 : 0) << ",\n";
        js << "  \"finite_fields_ok\": " << (fields.ok ? 1 : 0) << ",\n";
        js << "  \"numerical_quality_ok\": " << (numerical_quality_ok ? 1 : 0) << ",\n";
        js << "  \"validation_passed\": null,\n";
        js << "  \"rho_min\": " << fields.rho_min << ",\n";
        js << "  \"rho_max\": " << fields.rho_max << ",\n";
        js << "  \"rho_rel_max\": " << rho_rel_max << ",\n";
        js << "  \"budget_exhausted\": " << (budget_exhausted ? 1 : 0) << ",\n";
        js << "  \"diverged\": " << (diverged ? 1 : 0) << ",\n";
        js << "  \"fail_reason\": \"" << fail_reason << "\",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"em_enabled\": 0,\n";
        js << "  \"boussinesq_enabled\": 0,\n";
        js << "  \"thermocouple_probe\": {\n";
        js << "    \"id\": \"" << thermo_probe_def.probe_id << "\",\n";
        js << "    \"stl\": \"" << clientStirringThermocoupleStlName() << "\",\n";
        js << "    \"status\": \"" << thermo_probe_def.status << "\",\n";
        js << "    \"solver_m\": [" << thermo_probe_pos[0] << ", " << thermo_probe_pos[1] << ", "
           << thermo_probe_pos[2] << "],\n";
        js << "    \"target_temperature_k\": " << thermo_probe_def.target_temperature_k << ",\n";
        js << "    \"target_meaning\": \"" << thermo_probe_def.target_meaning << "\",\n";
        js << "    \"note\": \"Observer samples glass Temperature/Velocity/Pressure; F0 Temperature held at T0\"\n";
        js << "  },\n";
        js << "  \"note\": \"Doc242: workflow≠numerical_quality; solid uncovered FAIL; drive/load flags split; dp from Reload\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][client-f0] summary_json=" << summary_path.string() << "\n";
    return pass ? 0 : 1;
#endif
}
