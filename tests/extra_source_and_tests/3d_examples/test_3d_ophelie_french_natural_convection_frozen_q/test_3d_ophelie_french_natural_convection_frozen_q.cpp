/**
 * @file test_3d_ophelie_french_natural_convection_frozen_q.cpp
 * @brief Stage 4.1: frozen Q_50kW (default) or periodic thermo–EM coupling.
 *
 * Side/bottom no-slip wall with a short open rim above the melt (not a lid).
 * Hydrostatic density is applied, then a gravity-only settle, before Q / Boussinesq.
 * Periodic coupling defaults SPH conduction on (`--no-thermal-diffusion` to disable).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_french_skull_resistance.h"
#include "electromagnetic_ophelie_french_skull_shell.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "io_environment.h"
#include "interaction_ck.h"
#include "sphinxsys.h"

#include "../test_3d_ophelie_rh200_glass_em_stirring/rh200_joule_heat_grid.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

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
    const StdVec<std::string> candidates = {
        "./reload",
        "../reload",
        "../../../../../reload",
        "../test_3d_ophelie_french_natural_glass_relax/bin/reload",
        "../../test_3d_ophelie_french_natural_glass_relax/bin/reload",
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
    Real end_time = 0.02;
    Real u_ref = 0.0;              // 0 -> sqrt(gH); was 0.05 and delayed the first advection
    Real c0 = 0.0;                 // 0 -> 10*sqrt(gH); was 5 and too soft for hydrostatics
    Real beta = 1.0e-5;
    Real gravity_g = 9.81;
    Real wall_thickness_factor = 4.0; // still-water / stirring dummy-wall width; 2 dp leaks through the floor
    Real wall_rim_height = -1.0;   // <0 -> 6*dp splash guard above the melt (still open top)
    Real hydro_settle_time = 0.15; // gravity-only settle before Q / Boussinesq
    Real escaped_frac_max = 0.01;
    Real state_record_interval = 0.0; // 0 -> no VTP
    /** Start time-averaged meridional (r,z) sampling after this physical time (post-settle clock). */
    Real meridional_average_start = 0.5;
    /** P1-03: if >0, prescribe hot-bottom linear T: T=T0+dT*(1-ξ), ξ=(z-zb)/H (independent of EM). */
    Real prescribe_t_hot_bottom_dT = 0.0;
    /** Doc77 §3: if >0, prescribe radial hot-core T=T0+dT*(1-r/R); overrides hot-bottom if both set. */
    Real prescribe_t_radial_dT = 0.0;
    bool freeze_prescribed_temperature = false;
    bool prescribe_radial_z_taper = false;
    bool require_validation = false;
    bool thermal_diffusion_explicit = false;
    bool boundary_area_normalize = true;
    bool enable_self_induction = true;
    bool enable_target_power = true;
    /**
     * Doc95: quieter WCSPH defaults from Doc91/92 (free-surface: dens-sum normalize +
     * analytic p_h + LinearCorrectionCK). Off by default for regression; enable for lit-β SNR.
     */
    bool dens_sum_normalize = false;
    /**
     * Doc101/102 opt-in: apply dens_sum only for z <= z_free - surface_depth (keep FreeSurface
     * density in the surface band). Default false = historical full-domain overwrite.
     */
    bool dens_sum_interior_only = false;
    /**
     * Doc105 P2b opt-in: after dens update, for incomplete-support particles (rho_sum < rho0)
     * set Density so EOS pressure cancels BackgroundPressure: p_eos = −p_h ⇒ p_total ≈ 0
     * (atmospheric free-surface condition with background decomposition).
     * Requires BackgroundPressure field (use with --background-hydro-p). Default off.
     */
    bool dens_fs_atm_consistent = false;
    bool background_hydro_p = false;
    bool kernel_correction = false;
    /**
     * Doc105 P2d opt-in: AcousticDensMode for wall dens dissipation.
     * 0 = legacy (default): dens jump on Adami momentum pressure (−C_ext).
     * 1 = dens on p' with ReferenceBodyAcceleration; French keeps total-p momentum
     *     (PerturbationMomentum=0) + GravityForceCK, so a_ref must be −ĝ (g in ForcePrior).
     * Do NOT copy closed-box P2b (PM=1, GravityForceCK off, a_ref=0).
     */
    int acoustic_dens_mode = 0;
    /** Doc97 P4: scale acoustic + advection CFL together (1=default; 0.5 => ~dt/2). */
    Real dt_scale = 1.0;
    Real acoustic_cfl = 0.55;
    Real advection_cfl = 0.25;
    /** Doc97 P4: dens_sum/sigma0 vs FreeSurface audit CSV. */
    bool dens_sum_audit = false;
    int dens_sum_audit_interval = 50; // advection steps
    /** Doc97 P4: post-reload position jitter amplitude as fraction of dp (0=off). */
    Real layout_jitter_dp_frac = 0.0;
    unsigned layout_seed = 1;
    /**
     * Doc103 opt-in: every N advection steps, refresh sigma0=z-band particles from current rho_sum
     * (0=off). Intended with interior_only to reduce stale sigma0 after migration.
     */
    int dens_sum_surface_recapture_interval = 0;
    /**
     * Doc105 P2: after settle, one dens→ForcePrior→acoustic stage dump with zone stats; then exit.
     * Uniform T / zero-buoyancy recommended. Does not migrate closed-box P2b to free surface.
     */
    bool doc105_p2_hydro_stage_probe = false;
    /**
     * Doc105 P2c: after settle, run quiet keep for this many seconds (zero-buoyancy recommended),
     * dump zone residuals periodically, score pollution budget, then exit (skip EM/long thermo).
     */
    Real doc105_p2_quiet_keep = 0.0;
    /**
     * Doc105 P2f diagnostic: after dens/bgp(/fsatm), force Density=rho0 so p_eos=0 and
     * p_total=p_h only. Isolates SPH discrete ∇p_h error from dens_sum packing error.
     * Opt-in; never a production default. Applied only in stage/keep probes (not settle).
     */
    bool hydro_rho_freeze_rho0 = false;
    /**
     * Doc105 P2h: apply dens_sum_normalize only during hydro settle (layout compaction),
     * then disable for measurement/probe so FreeSurface path is used without extra normalize.
     * Opt-in; never a production default.
     */
    bool settle_dens_sum_normalize = false;
    /**
     * Doc105 P2i: after FreeSurface dens (and only when extra normalize is NOT writing),
     * if DensitySummation >= rho0 (complete/over-complete support), force Density=rho0.
     * Fixes one-sided FreeSurface max(ρ_sum,ρ0) overpack p_eos after settle packing.
     * Opt-in; never a production default. Distinct from probe-only hydro_rho_freeze_rho0.
     */
    bool dens_fs_cap_overpack = false;
    /**
     * Doc105 P2l/P2m diagnostic: after each acoustic step in quiet-keep, zero Velocity on
     * particles with DensitySummation < rho0. Bounds whether U_rms floor is carried by
     * incomplete-support motion. Opt-in; never a production default.
     */
    bool incomplete_vel_zero = false;
    /**
     * Doc105 P2n: after AcousticStep1stHalf, for ρ_sum<ρ0 undo vel kick and set
     * Force=−ForcePrior (a_net=0) before 2nd half. Opt-in diagnostic; never default.
     */
    bool incomplete_force_cancel_live = false;
    /**
     * Doc105 P2r: after AcousticStep1stHalf, for surface-band incomplete (ρ_sum<ρ0),
     * Shepard-like packing scale Force *= ρ0/ρ_sum and correct the 1st-half vel kick.
     * Distinct from force-cancel. Opt-in; never a production default.
     */
    bool fs_incomplete_pack_scale_force = false;
    /**
     * Doc105 P2o diagnostic: after each acoustic step in quiet-keep, zero Velocity on
     * wall-side / wall-corner / wall-surface band particles (geometric). Bounds whether
     * U_rms floor is carried by the side band. Opt-in; never a production default.
     */
    bool side_band_vel_zero = false;
    /**
     * Doc105 P2q diagnostic: after each acoustic step in quiet-keep, zero Velocity on
     * WallBottom (+ WallCorner) particles. Bounds residual U after side-band freeze.
     * Opt-in; never a production default.
     */
    bool bottom_band_vel_zero = false;
    Real target_power_w = 50000.0;
    std::string summary_json;
    FrenchLiteratureParameters literature = makeFrenchNaturalLiteratureParameters();
    FrenchSkullResistanceParams skull;
    FrenchFreeCrustParams free_crust;
    FrenchSkullShellParams skull_shell;
};

inline void applyLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, LocalCli &cli)
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
        else if (std::strncmp(av[i], "--u-ref=", 8) == 0)
        {
            cli.u_ref = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--c0=", 5) == 0)
        {
            cli.c0 = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--beta=", 7) == 0)
        {
            cli.beta = static_cast<Real>(std::atof(av[i] + 7));
        }
        else if (std::strncmp(av[i], "--target-power=", 15) == 0)
        {
            cli.enable_target_power = true;
            cli.target_power_w = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.enable_self_induction = false;
            cli.literature.enable_aind = false;
        }
        else if (tryApplyFrenchLiteratureCli(av[i], cli.literature, cli.thermal_diffusion_explicit))
        {
            if (std::strncmp(av[i], "--target-absorbed-power=", 24) == 0)
            {
                cli.enable_target_power = true;
                cli.target_power_w = cli.literature.target_glass_absorbed_power_w;
            }
            else if (std::strncmp(av[i], "--aind=", 7) == 0)
            {
                cli.enable_self_induction = cli.literature.enable_aind;
            }
            else if (std::strcmp(av[i], "--no-boussinesq") == 0 || std::strcmp(av[i], "--zero-buoyancy") == 0)
            {
                // literature.enable_boussinesq already cleared
            }
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
        else if (std::strncmp(av[i], "--wall-thickness-factor=", 24) == 0)
        {
            cli.wall_thickness_factor = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--wall-rim=", 11) == 0)
        {
            cli.wall_rim_height = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--hydro-settle-time=", 20) == 0)
        {
            cli.hydro_settle_time = static_cast<Real>(std::atof(av[i] + 20));
        }
        else if (std::strncmp(av[i], "--escaped-frac-max=", 19) == 0)
        {
            cli.escaped_frac_max = static_cast<Real>(std::atof(av[i] + 19));
        }
        else if (std::strncmp(av[i], "--state-record-interval=", 24) == 0)
        {
            cli.state_record_interval = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--meridional-average-start=", 27) == 0)
        {
            cli.meridional_average_start = static_cast<Real>(std::atof(av[i] + 27));
        }
        else if (std::strncmp(av[i], "--prescribe-t-hot-bottom=", 25) == 0)
        {
            cli.prescribe_t_hot_bottom_dT = static_cast<Real>(std::atof(av[i] + 25));
        }
        else if (std::strncmp(av[i], "--prescribe-t-radial=", 21) == 0)
        {
            cli.prescribe_t_radial_dT = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strcmp(av[i], "--freeze-prescribed-t") == 0)
        {
            cli.freeze_prescribed_temperature = true;
        }
        else if (std::strcmp(av[i], "--prescribe-radial-z-taper") == 0)
        {
            cli.prescribe_radial_z_taper = true;
        }
        else if (std::strcmp(av[i], "--require-validation") == 0)
        {
            cli.require_validation = true;
        }
        else if (std::strcmp(av[i], "--boundary-area-normalize") == 0)
        {
            cli.boundary_area_normalize = true;
        }
        else if (std::strcmp(av[i], "--no-boundary-area-normalize") == 0)
        {
            cli.boundary_area_normalize = false;
        }
        else if (std::strcmp(av[i], "--dens-sum-normalize") == 0)
        {
            cli.dens_sum_normalize = true;
        }
        else if (std::strcmp(av[i], "--no-dens-sum-normalize") == 0 || std::strcmp(av[i], "--nodens") == 0 ||
                 std::strcmp(av[i], "--no-extra-normalize") == 0)
        {
            // Doc105: nodens / no-extra-normalize = disable extra dens_sum_normalize only
            // (DensitySummationCK + FreeSurface regularization still run).
            cli.dens_sum_normalize = false;
            cli.dens_sum_interior_only = false;
        }
        else if (std::strcmp(av[i], "--doc105-p2-hydro-stage-probe") == 0)
        {
            cli.doc105_p2_hydro_stage_probe = true;
        }
        else if (std::strncmp(av[i], "--doc105-p2-quiet-keep=", 23) == 0)
        {
            cli.doc105_p2_quiet_keep = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strcmp(av[i], "--hydro-rho-freeze-rho0") == 0)
        {
            cli.hydro_rho_freeze_rho0 = true;
        }
        else if (std::strcmp(av[i], "--no-hydro-rho-freeze-rho0") == 0)
        {
            cli.hydro_rho_freeze_rho0 = false;
        }
        else if (std::strcmp(av[i], "--settle-dens-sum-normalize") == 0)
        {
            cli.settle_dens_sum_normalize = true;
        }
        else if (std::strcmp(av[i], "--no-settle-dens-sum-normalize") == 0)
        {
            cli.settle_dens_sum_normalize = false;
        }
        else if (std::strcmp(av[i], "--dens-fs-cap-overpack") == 0)
        {
            cli.dens_fs_cap_overpack = true;
        }
        else if (std::strcmp(av[i], "--no-dens-fs-cap-overpack") == 0)
        {
            cli.dens_fs_cap_overpack = false;
        }
        else if (std::strcmp(av[i], "--incomplete-vel-zero") == 0)
        {
            cli.incomplete_vel_zero = true;
        }
        else if (std::strcmp(av[i], "--no-incomplete-vel-zero") == 0)
        {
            cli.incomplete_vel_zero = false;
        }
        else if (std::strcmp(av[i], "--incomplete-force-cancel-live") == 0)
        {
            cli.incomplete_force_cancel_live = true;
        }
        else if (std::strcmp(av[i], "--no-incomplete-force-cancel-live") == 0)
        {
            cli.incomplete_force_cancel_live = false;
        }
        else if (std::strcmp(av[i], "--fs-incomplete-pack-scale-force") == 0)
        {
            cli.fs_incomplete_pack_scale_force = true;
        }
        else if (std::strcmp(av[i], "--no-fs-incomplete-pack-scale-force") == 0)
        {
            cli.fs_incomplete_pack_scale_force = false;
        }
        else if (std::strcmp(av[i], "--side-band-vel-zero") == 0)
        {
            cli.side_band_vel_zero = true;
        }
        else if (std::strcmp(av[i], "--no-side-band-vel-zero") == 0)
        {
            cli.side_band_vel_zero = false;
        }
        else if (std::strcmp(av[i], "--bottom-band-vel-zero") == 0)
        {
            cli.bottom_band_vel_zero = true;
        }
        else if (std::strcmp(av[i], "--no-bottom-band-vel-zero") == 0)
        {
            cli.bottom_band_vel_zero = false;
        }
        else if (std::strcmp(av[i], "--dens-sum-interior-only") == 0)
        {
            cli.dens_sum_normalize = true;
            cli.dens_sum_interior_only = true;
        }
        else if (std::strcmp(av[i], "--no-dens-sum-interior-only") == 0)
        {
            cli.dens_sum_interior_only = false;
        }
        else if (std::strcmp(av[i], "--dens-fs-atm-consistent") == 0)
        {
            cli.dens_fs_atm_consistent = true;
            cli.background_hydro_p = true; // needs p_h for p_eos = -p_h
        }
        else if (std::strcmp(av[i], "--no-dens-fs-atm-consistent") == 0)
        {
            cli.dens_fs_atm_consistent = false;
        }
        else if (std::strcmp(av[i], "--background-hydro-p") == 0)
        {
            cli.background_hydro_p = true;
        }
        else if (std::strcmp(av[i], "--no-background-hydro-p") == 0)
        {
            cli.background_hydro_p = false;
        }
        else if (std::strcmp(av[i], "--kernel-correction") == 0)
        {
            cli.kernel_correction = true;
        }
        else if (std::strcmp(av[i], "--no-kernel-correction") == 0)
        {
            cli.kernel_correction = false;
        }
        else if (std::strncmp(av[i], "--acoustic-dens-mode=", 21) == 0)
        {
            cli.acoustic_dens_mode = std::atoi(av[i] + 21);
            if (cli.acoustic_dens_mode != 0 && cli.acoustic_dens_mode != 1)
            {
                std::cerr << "[ophelie] --acoustic-dens-mode must be 0 or 1\n";
                std::exit(1);
            }
        }
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
        {
            cli.summary_json = std::string(av[i] + 15);
        }
        else if (std::strncmp(av[i], "--dt-scale=", 11) == 0)
        {
            cli.dt_scale = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--acoustic-cfl=", 15) == 0)
        {
            cli.acoustic_cfl = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strncmp(av[i], "--advection-cfl=", 16) == 0)
        {
            cli.advection_cfl = static_cast<Real>(std::atof(av[i] + 16));
        }
        else if (std::strcmp(av[i], "--dens-sum-audit") == 0)
        {
            cli.dens_sum_audit = true;
        }
        else if (std::strncmp(av[i], "--dens-sum-audit-interval=", 26) == 0)
        {
            cli.dens_sum_audit = true;
            cli.dens_sum_audit_interval = std::max(1, std::atoi(av[i] + 26));
        }
        else if (std::strncmp(av[i], "--layout-jitter-dp-frac=", 24) == 0)
        {
            cli.layout_jitter_dp_frac = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--layout-seed=", 14) == 0)
        {
            cli.layout_seed = static_cast<unsigned>(std::strtoul(av[i] + 14, nullptr, 10));
        }
        else if (std::strncmp(av[i], "--dens-sum-surface-recapture-interval=", 38) == 0)
        {
            cli.dens_sum_surface_recapture_interval = std::max(0, std::atoi(av[i] + 38));
        }
    }
    if (!(cli.dt_scale > TinyReal))
        cli.dt_scale = Real(1);
    if (!(cli.acoustic_cfl > TinyReal))
        cli.acoustic_cfl = Real(0.55);
    if (!(cli.advection_cfl > TinyReal))
        cli.advection_cfl = Real(0.25);
    if (cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
    if (cli.enable_target_power)
    {
        cli.literature.target_glass_absorbed_power_w = cli.target_power_w;
    }
    if (cli.prescribe_t_hot_bottom_dT > TinyReal || cli.prescribe_t_radial_dT > TinyReal)
    {
        // Independent buoyancy path: no EM heating; optional freeze keeps ΔT fixed.
        cli.literature.zero_q = true;
        if (cli.freeze_prescribed_temperature)
        {
            cli.literature.enable_thermal_diffusion = false;
            cli.thermal_diffusion_explicit = false;
        }
    }
    finalizeFrenchPeriodicLiteratureDefaults(cli.literature, cli.thermal_diffusion_explicit);
}

inline void finalizeNaturalHydroDefaults(LocalCli &cli, const OphelieFrenchReducedCaseParams &french)
{
    const Real glass_h = frenchReducedGlassHeight(french);
    const Real hydro_speed = std::sqrt(std::max(cli.gravity_g * glass_h, TinyReal));
    if (cli.c0 <= TinyReal)
    {
        cli.c0 = Real(10) * hydro_speed;
    }
    if (cli.u_ref <= TinyReal)
    {
        cli.u_ref = hydro_speed;
    }
    if (cli.wall_rim_height < Real(0))
    {
        cli.wall_rim_height = Real(6) * french.dp;
    }
}

struct GlassEscapeCensus
{
    size_t escaped = 0;
    size_t total = 0;
    Real r_max = 0.0;
    Real z_min = 0.0;
    Real z_max = 0.0;
};

inline GlassEscapeCensus hostCountEscapedGlassParticles(BaseParticles &particles, const Vecd &center, Real radius,
                                                        Real z_bottom, Real z_top, Real dp)
{
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    GlassEscapeCensus census;
    census.total = particles.TotalRealParticles();
    if (census.total == 0)
    {
        return census;
    }
    census.z_min = pos[0][2];
    census.z_max = pos[0][2];
    const Real r_limit = radius + Real(2) * dp;
    const Real z_lo = z_bottom - Real(2) * dp;
    const Real z_hi = z_top + Real(2) * dp;
    for (size_t i = 0; i < census.total; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        census.r_max = std::max(census.r_max, r);
        census.z_min = std::min(census.z_min, pos[i][2]);
        census.z_max = std::max(census.z_max, pos[i][2]);
        if (r > r_limit || pos[i][2] < z_lo || pos[i][2] > z_hi)
        {
            ++census.escaped;
        }
    }
    return census;
}

inline void hostApplyHydrostaticDensity(BaseParticles &particles, Real rho0, Real c0, Real gravity_g, Real z_top)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real depth = std::max(Real(0), z_top - pos[i][2]);
        rho[i] = rho0 * (Real(1) + gravity_g * depth * inv_c2);
    }
    syncVariableToDevice<Real>(particles, "Density");
}

inline void hostZeroVelocity(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        vel[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
}

inline void hostZeroVecFieldFrench(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Vecd>(particles, name);
    Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Vecd::Zero();
    syncVariableToDevice<Vecd>(particles, name);
}

/** Doc105 P2d: wall dens mode (fields registered by AcousticStep1stHalf). */
inline void hostSetAcousticDensModeFrench(BaseParticles &particles, Real mode)
{
    syncVariableToHost<Real>(particles, "AcousticDensMode");
    Real *m = particles.getVariableDataByName<Real>("AcousticDensMode");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        m[i] = mode;
    syncVariableToDevice<Real>(particles, "AcousticDensMode");
}

inline void hostSetReferenceBodyAccelerationFrench(BaseParticles &particles, const Vecd &a_ref)
{
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    Vecd *a = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        a[i] = a_ref;
    syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
}

inline void hostZeroRealFieldFrench(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Real>(particles, name);
    Real *v = particles.getVariableDataByName<Real>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Real(0);
    syncVariableToDevice<Real>(particles, name);
}

/** Doc105 P2/P2c zone tags for French free-surface hydrostatic residuals. */
enum class FrenchHydroZone : int
{
    Interior = 0,
    Surface = 1,
    WallSide = 2,
    WallBottom = 3,
    WallCorner = 4,     // side ∩ bottom
    WallSurface = 5,    // (side|bottom|corner) ∩ surface
};

struct FrenchZoneHydroStats
{
    size_t n = 0;
    size_t n_incomplete = 0; // DensitySummation < rho0
    Real uz_rms = 0;
    Real ur_rms = 0; // cylindrical radial velocity
    Real utheta_rms = 0;
    Real u_rms = 0; // Doc126: sqrt(mean(ux^2+uy^2+uz^2)), particle-count weighted
    Real a_prior_z_rms = 0;
    Real a_force_z_rms = 0;
    Real a_force_r_rms = 0; // Doc126: radial Force/m RMS (component budget)
    Real a_net_z_rms = 0;
    Real a_net_r_rms = 0; // radial accel residual
    Real drho_dt_rms = 0;
    Real rho_mean = 0;
    Real packing_mean = 0; // mean(rho_sum/rho0)
    Real packing_rms = 0;
    Real p_eos_rms = 0;
    Real p_eos_mean = 0; // Doc126 EOS identity: |mean| bound for rms
    Real p_bg_rms = 0;
    Real p_total_rms = 0;
    Real dp_h_rms = 0; // rms(p_total − ρ0 g (z_free−z)_+)
    Real mass_vol_mean = 0; // mean(mass / (ρ0 Vol)) — ~1 if mass=ρ0·Vol
};

struct FrenchHydroStageSnapshot
{
    FrenchZoneHydroStats all{};
    FrenchZoneHydroStats interior{};
    FrenchZoneHydroStats surface{};
    FrenchZoneHydroStats wall_side{};
    FrenchZoneHydroStats wall_bottom{};
    FrenchZoneHydroStats wall_corner{};
    FrenchZoneHydroStats wall_surface{};
    // Aggregate wall = side+bottom+corner (no surface) for Doc108 compatibility
    FrenchZoneHydroStats wall{};
    // Doc105 P2g: interior split by clearance from FS/walls (deep_band, typically 4·dp).
    FrenchZoneHydroStats deep_interior{};     // clearance ≥ deep_band from side/bot/FS
    FrenchZoneHydroStats near_fs_interior{};  // Interior zone but within deep_band of FS
    // Doc105 P2l: surface split by DensitySummation completeness.
    FrenchZoneHydroStats surface_incomplete{}; // Surface ∩ (ρ_sum < ρ0)
    FrenchZoneHydroStats surface_complete{};   // Surface ∩ (ρ_sum ≥ ρ0)
    // Doc105 P2o: wall-side split by support completeness.
    FrenchZoneHydroStats wall_side_incomplete{};
    FrenchZoneHydroStats wall_side_complete{};
    // Doc105 P2q: wall-bottom split by support completeness.
    FrenchZoneHydroStats wall_bottom_incomplete{};
    FrenchZoneHydroStats wall_bottom_complete{};
};

inline FrenchHydroZone hostClassifyFrenchHydroZone(const Vecd &pos, const Vecd &center, Real radius,
                                                   Real z_bottom, Real z_free, Real band)
{
    const Real dx = pos[0] - center[0];
    const Real dy = pos[1] - center[1];
    const Real r = std::sqrt(dx * dx + dy * dy);
    const bool near_surf = pos[2] > z_free - band;
    const bool near_side = r > radius - band;
    const bool near_bot = pos[2] < z_bottom + band;
    if (near_surf && (near_side || near_bot))
        return FrenchHydroZone::WallSurface;
    if (near_surf)
        return FrenchHydroZone::Surface;
    if (near_side && near_bot)
        return FrenchHydroZone::WallCorner;
    if (near_side)
        return FrenchHydroZone::WallSide;
    if (near_bot)
        return FrenchHydroZone::WallBottom;
    return FrenchHydroZone::Interior;
}

inline void hostAccumulateZone(FrenchZoneHydroStats &z, Real uz, Real ur, Real utheta, Real u_mag2, Real a_p,
                               Real a_f_z, Real a_f_r, Real a_n_z, Real a_n_r, Real drho, Real rho, Real packing,
                               Real p_eos, Real p_bg, Real p_tot, Real dp_h, bool incomplete, Real mass_vol)
{
    ++z.n;
    if (incomplete)
        ++z.n_incomplete;
    z.uz_rms += uz * uz;
    z.ur_rms += ur * ur;
    z.utheta_rms += utheta * utheta;
    z.u_rms += u_mag2;
    z.a_prior_z_rms += a_p * a_p;
    z.a_force_z_rms += a_f_z * a_f_z;
    z.a_force_r_rms += a_f_r * a_f_r;
    z.a_net_z_rms += a_n_z * a_n_z;
    z.a_net_r_rms += a_n_r * a_n_r;
    z.drho_dt_rms += drho * drho;
    z.rho_mean += rho;
    z.packing_mean += packing;
    z.packing_rms += packing * packing;
    z.p_eos_rms += p_eos * p_eos;
    z.p_eos_mean += p_eos;
    z.p_bg_rms += p_bg * p_bg;
    z.p_total_rms += p_tot * p_tot;
    z.dp_h_rms += dp_h * dp_h;
    z.mass_vol_mean += mass_vol;
}

inline void hostFinalizeZone(FrenchZoneHydroStats &z)
{
    if (z.n == 0)
        return;
    const Real inv = Real(1) / static_cast<Real>(z.n);
    z.uz_rms = std::sqrt(z.uz_rms * inv);
    z.ur_rms = std::sqrt(z.ur_rms * inv);
    z.utheta_rms = std::sqrt(z.utheta_rms * inv);
    z.u_rms = std::sqrt(z.u_rms * inv);
    z.a_prior_z_rms = std::sqrt(z.a_prior_z_rms * inv);
    z.a_force_z_rms = std::sqrt(z.a_force_z_rms * inv);
    z.a_force_r_rms = std::sqrt(z.a_force_r_rms * inv);
    z.a_net_z_rms = std::sqrt(z.a_net_z_rms * inv);
    z.a_net_r_rms = std::sqrt(z.a_net_r_rms * inv);
    z.drho_dt_rms = std::sqrt(z.drho_dt_rms * inv);
    z.rho_mean *= inv;
    z.packing_mean *= inv;
    z.packing_rms = std::sqrt(z.packing_rms * inv);
    z.p_eos_mean *= inv;
    z.p_eos_rms = std::sqrt(z.p_eos_rms * inv);
    z.p_bg_rms = std::sqrt(z.p_bg_rms * inv);
    z.p_total_rms = std::sqrt(z.p_total_rms * inv);
    z.dp_h_rms = std::sqrt(z.dp_h_rms * inv);
    z.mass_vol_mean *= inv;
}

inline FrenchHydroStageSnapshot hostFrenchHydroStageSnapshot(BaseParticles &particles, Real rho0, Real c0,
                                                            Real gravity_g, const Vecd &center, Real radius,
                                                            Real z_bottom, Real z_free, Real band,
                                                            Real deep_band)
{
    FrenchHydroStageSnapshot s;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real c2 = c0 * c0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real m = mass[i] + TinyReal;
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real rxy = std::sqrt(dx * dx + dy * dy) + TinyReal;
        const Real er_x = dx / rxy;
        const Real er_y = dy / rxy;
        const Real ur = vel[i][0] * er_x + vel[i][1] * er_y;
        const Real utheta = -vel[i][0] * er_y + vel[i][1] * er_x;
        const Real u_mag2 = vel[i].squaredNorm();
        const Real a_p = fp[i][2] / m;
        const Real a_f_z = f[i][2] / m;
        const Real a_f_r = (f[i][0] * er_x + f[i][1] * er_y) / m;
        const Real a_n_z = a_p + a_f_z;
        const Real a_n_r = ((fp[i][0] + f[i][0]) * er_x + (fp[i][1] + f[i][1]) * er_y) / m;
        const Real p_eos = c2 * (rho[i] - rho0);
        const Real p_tot = p_eos + p_bg[i];
        const Real depth = SMAX(z_free - pos[i][2], Real(0));
        const Real p_h_analytic = rho0 * gravity_g * depth;
        const Real dp_h = p_tot - p_h_analytic;
        const Real packing = rho_sum[i] / (rho0 + TinyReal);
        const bool incomplete = rho_sum[i] + TinyReal < rho0;
        const Real mass_vol = mass[i] / (rho0 * vol[i] + TinyReal);
        const FrenchHydroZone zone =
            hostClassifyFrenchHydroZone(pos[i], center, radius, z_bottom, z_free, band);
        auto acc = [&](FrenchZoneHydroStats &z) {
            hostAccumulateZone(z, vel[i][2], ur, utheta, u_mag2, a_p, a_f_z, a_f_r, a_n_z, a_n_r, drho[i], rho[i],
                               packing, p_eos, p_bg[i], p_tot, dp_h, incomplete, mass_vol);
        };
        acc(s.all);
        if (zone == FrenchHydroZone::Interior)
        {
            acc(s.interior);
            const Real r = rxy;
            const bool deep = (r < radius - deep_band) && (pos[i][2] > z_bottom + deep_band) &&
                              (pos[i][2] < z_free - deep_band);
            if (deep)
                acc(s.deep_interior);
            else if (pos[i][2] > z_free - deep_band)
                acc(s.near_fs_interior);
        }
        else if (zone == FrenchHydroZone::Surface)
        {
            acc(s.surface);
            if (incomplete)
                acc(s.surface_incomplete);
            else
                acc(s.surface_complete);
        }
        else if (zone == FrenchHydroZone::WallSide)
        {
            acc(s.wall_side);
            acc(s.wall);
            if (incomplete)
                acc(s.wall_side_incomplete);
            else
                acc(s.wall_side_complete);
        }
        else if (zone == FrenchHydroZone::WallBottom)
        {
            acc(s.wall_bottom);
            acc(s.wall);
            if (incomplete)
                acc(s.wall_bottom_incomplete);
            else
                acc(s.wall_bottom_complete);
        }
        else if (zone == FrenchHydroZone::WallCorner)
        {
            acc(s.wall_corner);
            acc(s.wall);
        }
        else
            acc(s.wall_surface);
    }
    hostFinalizeZone(s.all);
    hostFinalizeZone(s.interior);
    hostFinalizeZone(s.surface);
    hostFinalizeZone(s.wall_side);
    hostFinalizeZone(s.wall_bottom);
    hostFinalizeZone(s.wall_corner);
    hostFinalizeZone(s.wall_surface);
    hostFinalizeZone(s.wall);
    hostFinalizeZone(s.deep_interior);
    hostFinalizeZone(s.near_fs_interior);
    hostFinalizeZone(s.surface_incomplete);
    hostFinalizeZone(s.surface_complete);
    hostFinalizeZone(s.wall_side_incomplete);
    hostFinalizeZone(s.wall_side_complete);
    hostFinalizeZone(s.wall_bottom_incomplete);
    hostFinalizeZone(s.wall_bottom_complete);
    return s;
}

/**
 * Doc105 P2l diagnostic: for DensitySummation < ρ0, set Force = −ForcePrior so a_net=0.
 * Used only to bound how much surface residual is carried by incomplete-support particles.
 * Does not write Velocity; call on a force-snap after acoustic_1st (stage/keep end).
 */
inline size_t hostCancelNetForceIncompleteSupport(BaseParticles &particles, Real rho0)
{
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    size_t n_cancel = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal < rho0)
        {
            f[i] = -fp[i];
            ++n_cancel;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Force");
    return n_cancel;
}

/**
 * Doc105 P2o: zero Velocity on geometric side band (WallSide / WallCorner / WallSurface).
 */
inline size_t hostZeroVelocitySideBand(BaseParticles &particles, const Vecd &center, Real radius,
                                       Real z_bottom, Real z_free, Real band)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    size_t n_zero = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const FrenchHydroZone zone =
            hostClassifyFrenchHydroZone(pos[i], center, radius, z_bottom, z_free, band);
        if (zone == FrenchHydroZone::WallSide || zone == FrenchHydroZone::WallCorner ||
            zone == FrenchHydroZone::WallSurface)
        {
            vel[i] = Vecd::Zero();
            ++n_zero;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return n_zero;
}

/**
 * Doc105 P2q: zero Velocity on geometric bottom band (WallBottom / WallCorner).
 * Corner overlaps side0 when both are enabled; that is intentional for isolation trials.
 */
inline size_t hostZeroVelocityBottomBand(BaseParticles &particles, const Vecd &center, Real radius,
                                         Real z_bottom, Real z_free, Real band)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    size_t n_zero = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const FrenchHydroZone zone =
            hostClassifyFrenchHydroZone(pos[i], center, radius, z_bottom, z_free, band);
        if (zone == FrenchHydroZone::WallBottom || zone == FrenchHydroZone::WallCorner)
        {
            vel[i] = Vecd::Zero();
            ++n_zero;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return n_zero;
}

/** Doc105 P2m: zero Velocity on incomplete-support particles (diagnostic). */
inline size_t hostZeroVelocityIncompleteSupport(BaseParticles &particles, Real rho0)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    size_t n_zero = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal < rho0)
        {
            vel[i] = Vecd::Zero();
            ++n_zero;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return n_zero;
}

/**
 * Doc105 P2n: after AcousticStep1stHalf Update (vel += (Fprior+F)/m*dt), undo that kick on
 * incomplete particles in the free-surface band only (z > z_free - surface_depth), and set
 * Force=−ForcePrior. Restricting to the FS band avoids wrecking near-wall incomplete support.
 */
inline size_t hostUndoIncompleteAcoustic1stKick(BaseParticles &particles, Real rho0, Real dt, Real z_free,
                                                Real surface_depth)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real z_surf_lo = z_free - surface_depth;
    const size_t n = particles.TotalRealParticles();
    size_t n_fix = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal < rho0 && pos[i][2] > z_surf_lo)
        {
            const Real inv_m = Real(1) / (mass[i] + TinyReal);
            vel[i] -= (fp[i] + f[i]) * inv_m * dt;
            f[i] = -fp[i];
            ++n_fix;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "Force");
    return n_fix;
}

/**
 * Doc105 P2r: Shepard-like packing scale of pressure Force on surface-band incomplete particles.
 * After AcousticStep1stHalf Update, Force *= ρ0/ρ_sum and correct vel by the Force delta.
 * If correct_velocity=false (stage snap), only Force is scaled.
 */
inline size_t hostScaleIncompleteFsForceByPacking(BaseParticles &particles, Real rho0, Real dt, Real z_free,
                                                  Real surface_depth, bool correct_velocity)
{
    if (correct_velocity)
        syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Vecd *vel = correct_velocity ? particles.getVariableDataByName<Vecd>("Velocity") : nullptr;
    Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real z_surf_lo = z_free - surface_depth;
    const size_t n = particles.TotalRealParticles();
    size_t n_scale = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal < rho0 && pos[i][2] > z_surf_lo)
        {
            const Real s = rho0 / (rho_sum[i] + TinyReal);
            if (correct_velocity && vel != nullptr)
            {
                const Real inv_m = Real(1) / (mass[i] + TinyReal);
                vel[i] += f[i] * (s - Real(1)) * inv_m * dt;
            }
            f[i] *= s;
            ++n_scale;
        }
    }
    if (correct_velocity)
        syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "Force");
    return n_scale;
}

inline void hostWriteFrenchZoneJson(std::ostream &js, const char *name, const FrenchZoneHydroStats &z)
{
    // Doc126: particle-count weighted RMS; U_rms is 3D magnitude (not Uz alias).
    const Real eos_floor = std::abs(z.p_eos_mean);
    const bool eos_identity_ok = (z.n == 0) || (z.p_eos_rms + TinyReal >= eos_floor);
    js << "  \"" << name << "\": {\"n\": " << z.n << ", \"n_incomplete\": " << z.n_incomplete
       << ", \"uz_rms\": " << z.uz_rms << ", \"ur_rms\": " << z.ur_rms << ", \"utheta_rms\": " << z.utheta_rms
       << ", \"u_rms\": " << z.u_rms << ", \"weighting\": \"particle_count\""
       << ", \"a_prior_z_rms\": " << z.a_prior_z_rms << ", \"a_force_z_rms\": " << z.a_force_z_rms
       << ", \"a_force_r_rms\": " << z.a_force_r_rms
       << ", \"a_net_z_rms\": " << z.a_net_z_rms << ", \"a_net_r_rms\": " << z.a_net_r_rms
       << ", \"drho_dt_rms\": " << z.drho_dt_rms << ", \"rho_mean\": " << z.rho_mean
       << ", \"packing_mean\": " << z.packing_mean << ", \"packing_rms\": " << z.packing_rms
       << ", \"p_eos_rms\": " << z.p_eos_rms << ", \"p_eos_mean\": " << z.p_eos_mean
       << ", \"eos_identity_ok\": " << (eos_identity_ok ? 1 : 0)
       << ", \"p_bg_rms\": " << z.p_bg_rms
       << ", \"p_total_rms\": " << z.p_total_rms << ", \"dp_h_rms\": " << z.dp_h_rms
       << ", \"mass_vol_mean\": " << z.mass_vol_mean << "}";
}

/** Doc126 P0: snapshot evolving fields so acoustic probes can restore state. */
struct FrenchHydroFieldSnapshot
{
    std::vector<Vecd> velocity;
    std::vector<Vecd> force;
    std::vector<Vecd> force_prior;
    std::vector<Vecd> displacement;
    std::vector<Real> density;
    std::vector<Real> density_change_rate;
    std::vector<Real> pressure;
};

inline void hostCaptureFrenchHydroFields(BaseParticles &particles, FrenchHydroFieldSnapshot &s)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Displacement");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    syncVariableToHost<Real>(particles, "Pressure");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *dpos = particles.getVariableDataByName<Vecd>("Displacement");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const Real *p = particles.getVariableDataByName<Real>("Pressure");
    s.velocity.assign(vel, vel + n);
    s.force.assign(f, f + n);
    s.force_prior.assign(fp, fp + n);
    s.displacement.assign(dpos, dpos + n);
    s.density.assign(rho, rho + n);
    s.density_change_rate.assign(drho, drho + n);
    s.pressure.assign(p, p + n);
}

inline void hostRestoreFrenchHydroFields(BaseParticles &particles, const FrenchHydroFieldSnapshot &s)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Displacement");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    syncVariableToHost<Real>(particles, "Pressure");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Vecd *dpos = particles.getVariableDataByName<Vecd>("Displacement");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    Real *p = particles.getVariableDataByName<Real>("Pressure");
    for (size_t i = 0; i < n; ++i)
    {
        vel[i] = s.velocity[i];
        f[i] = s.force[i];
        fp[i] = s.force_prior[i];
        dpos[i] = s.displacement[i];
        rho[i] = s.density[i];
        drho[i] = s.density_change_rate[i];
        p[i] = s.pressure[i];
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "Force");
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    syncVariableToDevice<Vecd>(particles, "Displacement");
    syncVariableToDevice<Real>(particles, "Density");
    syncVariableToDevice<Real>(particles, "DensityChangeRate");
    syncVariableToDevice<Real>(particles, "Pressure");
}

/** Doc105 P2f: force Density=ρ0 (p_eos=0); keeps BackgroundPressure / mass unchanged. */
inline void hostFreezeDensityToRho0(BaseParticles &particles, Real rho0)
{
    syncVariableToHost<Real>(particles, "Density");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        rho[i] = rho0;
    syncVariableToDevice<Real>(particles, "Density");
}

/**
 * Doc105 P2i: FreeSurface is one-sided (ρ=max(ρ_sum,ρ0)). After packing settle, ρ_sum>ρ0
 * would keep over-density. Cap those particles back to ρ0; incomplete (ρ_sum<ρ0) already ρ0.
 * Returns number of particles capped.
 */
inline size_t hostCapOverpackDensityToRho0(BaseParticles &particles, Real rho0)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    size_t n_capped = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (rho_sum[i] + TinyReal >= rho0)
        {
            rho[i] = rho0;
            ++n_capped;
        }
    }
    syncVariableToDevice<Real>(particles, "Density");
    return n_capped;
}

/** Doc105 dens-write audit: Δρ and inferred Δp=c0²Δρ when extra normalize overwrites Density. */
struct FrenchDensWriteAudit
{
    size_t n_written = 0;
    size_t n_skipped_surface = 0;
    Real max_abs_drho = 0;
    Real rms_drho = 0;
    Real max_abs_dp = 0;
    Real rms_dp = 0;
    Real mass_sum = 0;
};

inline FrenchDensWriteAudit hostApplyDensSumNormalizeFieldAudited(BaseParticles &particles, Real rho0, Real c0,
                                                                  bool interior_only, Real z_free,
                                                                  Real surface_depth)
{
    FrenchDensWriteAudit a;
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, "Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_surf_lo = z_free - surface_depth;
    const Real c2 = c0 * c0;
    double s_dr2 = 0.0, s_dp2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        a.mass_sum += mass[i];
        if (interior_only && pos[i][2] > z_surf_lo)
        {
            ++a.n_skipped_surface;
            continue;
        }
        const Real rho_new = rho0 * (rho_sum[i] / (sigma0[i] + TinyReal));
        const Real drho = rho_new - rho[i];
        const Real dp = c2 * drho;
        rho[i] = rho_new;
        ++a.n_written;
        a.max_abs_drho = std::max(a.max_abs_drho, std::abs(drho));
        a.max_abs_dp = std::max(a.max_abs_dp, std::abs(dp));
        s_dr2 += static_cast<double>(drho) * static_cast<double>(drho);
        s_dp2 += static_cast<double>(dp) * static_cast<double>(dp);
    }
    if (a.n_written > 0)
    {
        a.rms_drho = static_cast<Real>(std::sqrt(s_dr2 / static_cast<double>(a.n_written)));
        a.rms_dp = static_cast<Real>(std::sqrt(s_dp2 / static_cast<double>(a.n_written)));
    }
    syncVariableToDevice<Real>(particles, "Density");
    return a;
}

/**
 * Doc105 P2b: atmospheric free-surface density with background pressure.
 * Only incomplete-support particles IN THE FREE-SURFACE BAND
 *   (rho_sum < rho0 && z > z_free - surface_depth)
 * get rho = rho0 - p_h/c0² so p_eos + p_h ≈ 0.
 * Do NOT apply on wall-near incomplete support (rho_sum<rho0 alone) — that destroys hydrostatics.
 */
struct FrenchFsAtmAudit
{
    size_t n_candidates = 0; // rho_sum < rho0 in surface band
    size_t n_written = 0;
    size_t n_skipped_non_surface = 0; // rho_sum < rho0 but below surface band
    Real max_abs_drho = 0;
    Real rms_drho = 0;
    Real max_abs_p_total_after = 0;
    Real rms_p_total_after = 0;
};

inline FrenchFsAtmAudit hostApplyFrenchFsAtmConsistentDensity(BaseParticles &particles, Real rho0, Real c0,
                                                             Real z_free, Real surface_depth)
{
    FrenchFsAtmAudit a;
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    syncVariableToHost<Vecd>(particles, "Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real c2 = c0 * c0 + TinyReal;
    const Real rho_lo = Real(0.5) * rho0;
    const Real z_surf_lo = z_free - surface_depth;
    double s_dr2 = 0.0, s_pt2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        if (!(rho_sum[i] < rho0))
            continue;
        if (pos[i][2] <= z_surf_lo)
        {
            ++a.n_skipped_non_surface;
            continue;
        }
        ++a.n_candidates;
        Real rho_new = rho0 - p_bg[i] / c2;
        rho_new = std::min(rho0, std::max(rho_lo, rho_new));
        const Real drho = rho_new - rho[i];
        rho[i] = rho_new;
        ++a.n_written;
        const Real p_tot = c2 * (rho_new - rho0) + p_bg[i];
        a.max_abs_drho = std::max(a.max_abs_drho, std::abs(drho));
        a.max_abs_p_total_after = std::max(a.max_abs_p_total_after, std::abs(p_tot));
        s_dr2 += static_cast<double>(drho) * static_cast<double>(drho);
        s_pt2 += static_cast<double>(p_tot) * static_cast<double>(p_tot);
    }
    if (a.n_written > 0)
    {
        a.rms_drho = static_cast<Real>(std::sqrt(s_dr2 / static_cast<double>(a.n_written)));
        a.rms_p_total_after = static_cast<Real>(std::sqrt(s_pt2 / static_cast<double>(a.n_written)));
    }
    syncVariableToDevice<Real>(particles, "Density");
    return a;
}

/** Explicit heating that keeps OphelieThermalDeltaT synchronized with Temperature. */
class ApplyFrozenQToTemperatureCK : public LocalDynamics
{
  public:
    ApplyFrozenQToTemperatureCK(SPHBody &sph_body, Real rho, Real cp, Real t_initial)
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
            // Doc 45 F01: ΔT-authoritative heating (same contract as stirring production path).
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

inline void hostSampleEulerQOntoParticles(BaseParticles &particles, const rh200::Rh200ScalarEulerianGrid &grid)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "JouleHeat");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *q = particles.getVariableDataByName<Real>("JouleHeat");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        q[i] = grid.sampleHost(pos[i]);
    }
    syncVariableToDevice<Real>(particles, "JouleHeat");
}

inline Real hostThermalEnergy(BaseParticles &particles, Real rho, Real cp, Real t_ref)
{
    return hostThermalEnergyConstantCp(particles, rho, cp, t_ref, kOphelieTemperatureField);
}

/** ViscousForceCK caches μ; rebuild after mean-μ(T) refresh (Batch D / P1-01). */
inline void invalidateViscosityDependentCkKernels(Inner<> &glass_inner, Contact<> &glass_wall_contact)
{
    glass_inner.resetComputingKernelUpdated();
    const StdVec<SPHBody *> contacts = glass_wall_contact.getContactBodies();
    for (size_t k = 0; k < contacts.size(); ++k)
    {
        glass_wall_contact.resetComputingKernelUpdated(k);
    }
}

inline void hostCaptureDensSumSigma0Field(BaseParticles &particles, double &sigma0_mean)
{
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "DensSumZ0");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *z0 = particles.getVariableDataByName<Real>("DensSumZ0");
    const size_t n = particles.TotalRealParticles();
    sigma0_mean = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sigma0[i] = std::max(rho_sum[i], TinyReal);
        z0[i] = pos[i][2];
        sigma0_mean += static_cast<double>(sigma0[i]);
    }
    sigma0_mean /= static_cast<double>(std::max<size_t>(n, 1));
    syncVariableToDevice<Real>(particles, "DensSumSigma0");
    syncVariableToDevice<Real>(particles, "DensSumZ0");
}

/**
 * dens_sum normalize. If interior_only, particles with z > z_free - surface_depth keep the
 * Density already written by FreeSurface regularization (no overwrite).
 * Returns number of particles rewritten.
 */
inline size_t hostApplyDensSumNormalizeField(BaseParticles &particles, Real rho0, bool interior_only,
                                             Real z_free, Real surface_depth)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    syncVariableToHost<Vecd>(particles, "Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_surf_lo = z_free - surface_depth;
    size_t n_written = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (interior_only && pos[i][2] > z_surf_lo)
            continue; // keep FreeSurface rho
        rho[i] = rho0 * (rho_sum[i] / (sigma0[i] + TinyReal));
        ++n_written;
    }
    syncVariableToDevice<Real>(particles, "Density");
    return n_written;
}

/** Doc103: refresh sigma0 for surface-band particles from current rho_sum. Returns count updated. */
inline size_t hostRecaptureSurfaceSigma0(BaseParticles &particles, Real z_free, Real surface_depth)
{
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_surf_lo = z_free - surface_depth;
    size_t n_updated = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] <= z_surf_lo)
            continue;
        sigma0[i] = std::max(rho_sum[i], TinyReal);
        ++n_updated;
    }
    syncVariableToDevice<Real>(particles, "DensSumSigma0");
    return n_updated;
}

/**
 * Doc97 P4 dens_sum audit: FreeSurface regularization sets rho = max(rho_sum, rho0);
 * dens_sum_normalize then overwrites rho = rho0 * rho_sum / sigma0 (sigma0 frozen at first dens-sum).
 * DensSumSigma0 is an evolving particle field (moves with sorting), so surface particles may keep
 * interior packing references and vice versa.
 */
struct FrenchDensSumAuditSnapshot
{
    size_t n_total = 0;
    size_t n_fs_clamp = 0;          // rho_sum < rho0 (FreeSurface would clamp)
    size_t n_surface_band = 0;      // z > z_free - surface_depth
    size_t n_interior = 0;
    size_t n_migrated_to_surface = 0; // z0 deep, now surface
    size_t n_override_active = 0;   // |rho_fs - rho_norm| > tol among fs_clamp
    Real mean_abs_drho_override = 0;
    Real max_abs_drho_override = 0;
    Real mass_sum = 0;
    Real mean_ratio_surf = 0;  // mean(rho_sum/sigma0) surface
    Real mean_ratio_int = 0;
    Real mean_sigma0_surf = 0;
    Real mean_sigma0_int = 0;
    Real mean_rho_surf = 0;
    Real mean_rho_int = 0;
};

inline FrenchDensSumAuditSnapshot hostFrenchDensSumAudit(BaseParticles &particles, Real rho0, Real z_free,
                                                        Real surface_depth)
{
    FrenchDensSumAuditSnapshot s;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    syncVariableToHost<Real>(particles, "DensSumZ0");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const Real *z0 = particles.getVariableDataByName<Real>("DensSumZ0");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    s.n_total = particles.TotalRealParticles();
    const Real z_surf_lo = z_free - surface_depth;
    const Real z_deep_hi = z_free - Real(3) * surface_depth;
    const Real tol = Real(1.0e-3) * rho0;
    double sum_abs = 0, sum_rs = 0, sum_ri = 0, sum_s0s = 0, sum_s0i = 0, sum_rhos = 0, sum_rhoi = 0;
    size_t ns = 0, ni = 0;
    for (size_t i = 0; i < s.n_total; ++i)
    {
        s.mass_sum = static_cast<Real>(static_cast<double>(s.mass_sum) + static_cast<double>(mass[i]));
        const bool fs_clamp = rho_sum[i] + TinyReal < rho0;
        if (fs_clamp)
            ++s.n_fs_clamp;
        const Real rho_fs = std::max(rho_sum[i], rho0); // FreeSurface Regularization intent
        // Actual Density vs FreeSurface intent (captures full-domain overwrite or interior_only keep).
        const Real adr = std::abs(rho[i] - rho_fs);
        if (fs_clamp && adr > tol)
        {
            ++s.n_override_active;
            sum_abs += static_cast<double>(adr);
            s.max_abs_drho_override = std::max(s.max_abs_drho_override, adr);
        }
        const bool surface = pos[i][2] > z_surf_lo;
        const bool was_deep = z0[i] < z_deep_hi;
        if (surface && was_deep)
            ++s.n_migrated_to_surface;
        const Real ratio = rho_sum[i] / (sigma0[i] + TinyReal);
        if (surface)
        {
            ++s.n_surface_band;
            ++ns;
            sum_rs += static_cast<double>(ratio);
            sum_s0s += static_cast<double>(sigma0[i]);
            sum_rhos += static_cast<double>(rho[i]);
        }
        else
        {
            ++s.n_interior;
            ++ni;
            sum_ri += static_cast<double>(ratio);
            sum_s0i += static_cast<double>(sigma0[i]);
            sum_rhoi += static_cast<double>(rho[i]);
        }
    }
    if (s.n_override_active > 0)
        s.mean_abs_drho_override = static_cast<Real>(sum_abs / static_cast<double>(s.n_override_active));
    if (ns > 0)
    {
        s.mean_ratio_surf = static_cast<Real>(sum_rs / static_cast<double>(ns));
        s.mean_sigma0_surf = static_cast<Real>(sum_s0s / static_cast<double>(ns));
        s.mean_rho_surf = static_cast<Real>(sum_rhos / static_cast<double>(ns));
    }
    if (ni > 0)
    {
        s.mean_ratio_int = static_cast<Real>(sum_ri / static_cast<double>(ni));
        s.mean_sigma0_int = static_cast<Real>(sum_s0i / static_cast<double>(ni));
        s.mean_rho_int = static_cast<Real>(sum_rhoi / static_cast<double>(ni));
    }
    return s;
}

inline void appendFrenchDensSumAuditCsv(const std::string &path, Real time, size_t adv_step,
                                        const FrenchDensSumAuditSnapshot &s, bool write_header)
{
    std::ofstream out(path, write_header ? std::ios::out : std::ios::app);
    if (write_header)
        out << "time_s,advection_step,n_total,n_fs_clamp,n_surface_band,n_interior,n_migrated_to_surface,"
               "n_override_active,mean_abs_drho_override,max_abs_drho_override,mass_sum,"
               "mean_ratio_surf,mean_ratio_int,mean_sigma0_surf,mean_sigma0_int,mean_rho_surf,mean_rho_int\n";
    out << time << "," << adv_step << "," << s.n_total << "," << s.n_fs_clamp << "," << s.n_surface_band << ","
        << s.n_interior << "," << s.n_migrated_to_surface << "," << s.n_override_active << ","
        << s.mean_abs_drho_override << "," << s.max_abs_drho_override << "," << s.mass_sum << ","
        << s.mean_ratio_surf << "," << s.mean_ratio_int << "," << s.mean_sigma0_surf << "," << s.mean_sigma0_int
        << "," << s.mean_rho_surf << "," << s.mean_rho_int << "\n";
}

inline void hostApplyLayoutJitter(BaseParticles &particles, Real amplitude, unsigned seed)
{
    if (!(amplitude > TinyReal))
        return;
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    for (size_t i = 0; i < n; ++i)
    {
        pos[i][0] += static_cast<Real>(uni(rng) * static_cast<double>(amplitude));
        pos[i][1] += static_cast<Real>(uni(rng) * static_cast<double>(amplitude));
        pos[i][2] += static_cast<Real>(uni(rng) * static_cast<double>(amplitude));
    }
    syncVariableToDevice<Vecd>(particles, "Position");
    std::cout << "[ophelie][doc97-p4] layout_jitter amplitude=" << amplitude << " seed=" << seed
              << " n=" << n << std::endl;
}

/** Free-surface melt: p_h = ρ0 g (z_free − z), clamped at free surface. */
inline void hostApplyFrenchBackgroundHydroPressure(BaseParticles &particles, Real rho0, Real g, Real z_free)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real depth = std::max(Real(0), z_free - pos[i][2]);
        p_bg[i] = rho0 * g * depth;
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
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
    LocalCli local_cli;
    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    // Natural geometry defaults (override reduced CAD placeholders).
    french.glass_radius = 0.25;
    french.glass_half_height = 0.5 * 0.185;
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    french.dp = 0.015;
    french.coil.loop_radius = 0.285;
    french.coil.num_loops = 7;
    french.coil.segments_per_loop = 64;
    french.coil.z_min = -0.0225;
    french.coil.z_max = 0.2075;
    french.auto_coil_z = false;
    params.frequency_ = 282000.0;
    params.sigma_glass_ = 16.0;

    applyLocalCli(ac, av, french, local_cli);
    if (local_cli.literature.coupling == ThermoEMCouplingMode::Periodic)
    {
        local_cli.enable_self_induction = local_cli.literature.enable_aind;
        local_cli.literature.target_glass_absorbed_power_w = local_cli.target_power_w;
    }
    const StdVec<std::string> french_filtered = filterFrenchReducedCommandLine(ac, av, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    OphelieTestCliOptions cli_options;
    StdVec<char *> french_av;
    for (auto &argument : french_filtered)
    {
        french_av.push_back(const_cast<char *>(argument.c_str()));
    }
    (void)filterOphelieTestCommandLine(static_cast<int>(french_av.size()), french_av.data(), params, cli_options);
    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
    }

    // Doc212: complex EdgeFlux must be set BEFORE finalize, otherwise filterOphelieTestCommandLine
    // leaves phi_rhs=DivΣA / lhs unset and final_params lie about the production profile.
    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    cli_options.ophelie_current_form_user_set = true;
    finalizeOphelieCurrentFormConfiguration(params, cli_options);
    params.enable_self_induction_ = local_cli.enable_self_induction;
    if (local_cli.enable_target_power)
    {
        params.target_joule_power_ = local_cli.target_power_w;
        params.enable_power_scaling_ = false;
    }
    if (params.self_induction_phi_eq_res_tolerance_ > Real(1.0e-3))
    {
        params.self_induction_phi_eq_res_tolerance_ = Real(2.0e-4);
    }
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    logOphelieFinalParams(params, cli_options);

    if (!reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_natural_convection_frozen_q: Reload.xml required under \""
                  << local_cli.reload_dir << "\"\n";
        return 1;
    }

    //----------------------------------------------------------------------
    // Phase A — EM on SolidBody; freeze Q_50kW particle field.
    // Keep EM SPHSystem alive for process lifetime: tearing down a SYCL
    // SPHSystem then creating another often hangs the GPU runtime.
    //----------------------------------------------------------------------
    StdVec<Real> frozen_q;
    Real p_joule = 0.0;
    Real phi_eq = 0.0;
    UniquePtr<SPHSystem> em_system;
    UniquePtr<SolidBody> glass_em;
    UniquePtr<Inner<>> glass_em_inner;
    UniquePtr<RegisterOphelieGlassFields> register_glass_em;
    OphelieGlassFieldNames glass_em_names;
    ThermoEMCouplingState coupling_state;
    rh200::Rh200ScalarEulerianGrid q_grid;
    const bool periodic_coupling = local_cli.literature.coupling == ThermoEMCouplingMode::Periodic;
    fs::create_directories("output");
    logFrenchLiteratureParameters(local_cli.literature);
    {
        const BoundingBoxd bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
        em_system = makeUnique<SPHSystem>(bounds, french.dp);
        em_system->setReloadParticles(true);
        IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

        glass_em = makeUnique<SolidBody>(
            *em_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                            "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height));
        glass_em->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        glass_em->defineMatterMaterial<Solid>();
        (void)defineOphelieSolidLevelSet(*glass_em);
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
        if (local_cli.enable_target_power)
        {
            (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
            runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
        }
        p_joule = probe.joule_power_w;
        phi_eq = probe.phi_eq_res_vol;
        coupling_state.calibrated_current_per_loop = french.coil.current_per_loop;
        coupling_state.coil_current_per_loop = french.coil.current_per_loop;
        coupling_state.current_frozen = local_cli.literature.em_control == EMControlMode::FixedCoilCurrent;
        coupling_state.generator_power_w = local_cli.literature.generator_power_w;
        coupling_state.target_glass_absorbed_power_w = local_cli.target_power_w;
        coupling_state.reconstructed_glass_power_w = p_joule;
        coupling_state.last_phi_eq_res = phi_eq;

        BaseParticles &particles = glass_em->getBaseParticles();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_em_names, params);
        syncVariableToHost<Real>(particles, q_field);
        const Real *q = particles.getVariableDataByName<Real>(q_field);
        const size_t n = particles.TotalRealParticles();
        frozen_q.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            frozen_q[i] = q[i];
        }
        std::cout << "[ophelie][stage4.1] EM freeze: n=" << n << " P_joule_W=" << p_joule
                  << " phi_eq_res_vol=" << phi_eq << " self_induction=" << (probe.used_self_induction ? 1 : 0)
                  << " generator_power=" << coupling_state.generator_power_w
                  << " target_glass_absorbed_power=" << coupling_state.target_glass_absorbed_power_w
                  << " reconstructed_glass_power=" << p_joule
                  << " coil_current=" << coupling_state.coil_current_per_loop << std::endl;

        if (periodic_coupling)
        {
            const Vecd &gc_em = french.glass_center;
            const FrenchCylindricalFrame frame = makeFrenchCylinderFrameFromGlass(
                gc_em, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
            const BoundingBoxd melt_bounds(
                Vecd(gc_em[0] - french.glass_radius, gc_em[1] - french.glass_radius,
                     gc_em[2] - french.glass_half_height),
                Vecd(gc_em[0] + french.glass_radius, gc_em[1] + french.glass_radius,
                     gc_em[2] + french.glass_half_height));
            q_grid.spec_ = rh200::makeRh200JouleHeatGridSpecFromBounds(melt_bounds, french.dp, Real(2));
            coupling_state.sigma_rz.frame = frame;
            coupling_state.q_rz.frame = frame;
            depositEmQAndRecord(particles, glass_em_names, params, frame, local_cli.literature, q_grid, coupling_state,
                                p_joule);
            syncVariableToHost<Real>(particles, glass_em_names.sigma);
            const Real *sigma0 = particles.getVariableDataByName<Real>(glass_em_names.sigma);
            const Vecd *pos0 = particles.getVariableDataByName<Vecd>("Position");
            const Real *vol0 = particles.getVariableDataByName<Real>("VolumetricMeasure");
            azimuthalDeposit(coupling_state.sigma_rz, pos0, sigma0, vol0, n);
            fillEmptyAxisymmetricBins(coupling_state.sigma_rz, params.sigma_glass_);
            coupling_state.sigma_em.assign(sigma0, sigma0 + n);
            sigmaFieldStats(coupling_state.sigma_em, vol0, n, coupling_state.last_sigma_min,
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
    // Phase B — FluidBody WCSPH + Boussinesq + Natural thermal BC.
    // Domain is glass±wall only (not coil): Lattice scans the full SPHSystem
    // AABB, so a coil-sized box makes wall generation extremely slow.
    //----------------------------------------------------------------------
    finalizeNaturalHydroDefaults(local_cli, french);
    {
        std::string shell_err;
        if (!assertFrenchSkullShellExclusive(local_cli.skull_shell, local_cli.skull, false, shell_err))
        {
            std::cerr << "[ophelie][stage4.1] " << shell_err << std::endl;
            return 1;
        }
    }
    Real wall_thickness = local_cli.wall_thickness_factor * french.dp;
    if (local_cli.skull_shell.enable)
    {
        // Hydro containment needs O(few·dp). At natural dp=0.015, a 5 mm thermal δ cannot be
        // spatially resolved; keep hydro wall ≥ factor·dp and treat δ as the thermal OPEN_PARAMETER.
        const Real hydro_min = local_cli.wall_thickness_factor * french.dp;
        wall_thickness = std::max(hydro_min, local_cli.skull_shell.delta_m);
        logFrenchSkullShellParams(local_cli.skull_shell, wall_thickness);
        if (wall_thickness > local_cli.skull_shell.delta_m + TinyReal)
        {
            std::cout << "[ophelie][skull-shell] WARNING hydro wall_thickness=" << wall_thickness
                      << " m > thermal delta=" << local_cli.skull_shell.delta_m
                      << " m (dp=" << french.dp << "); resolved shell needs finer dp.\n";
        }
        // Contact path is wired below; no Level-1 melt Robin bridge.
    }
    const Real rim_height = local_cli.wall_rim_height;
    const Real flow_margin = Real(4.0) * french.dp + wall_thickness;
    const Vecd &gc = french.glass_center;
    const Real flow_r = french.glass_radius + wall_thickness + flow_margin;
    const Real flow_down = french.glass_half_height + wall_thickness + flow_margin;
    const Real flow_up = french.glass_half_height + rim_height + wall_thickness + flow_margin;
    const BoundingBoxd flow_bounds(Vecd(gc[0] - flow_r, gc[1] - flow_r, gc[2] - flow_down),
                                   Vecd(gc[0] + flow_r, gc[1] + flow_r, gc[2] + flow_up));
    SPHSystem sph_system(flow_bounds, french.dp);
    sph_system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

    const OphelieFrenchGlassMaterialLaws laws = local_cli.literature.material;
    const Real t0 = local_cli.literature.t_initial_k;
    const Real rho0 = laws.rho.evaluate(t0);
    const Real mu = laws.mu.evaluate(t0);
    const Real cp = laws.cp.evaluate(t0);
    const Real k_th = laws.k.evaluate(t0);
    const Real beta = local_cli.literature.enable_boussinesq
                          ? (local_cli.beta > TinyReal ? local_cli.beta : laws.beta.evaluate(t0))
                          : Real(0);

    FluidBody glass(sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                                    "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, local_cli.c0);
    glass.addMaterialProperty<Viscosity>(mu);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());
    Viscosity &glass_viscosity = glass.getMaterialProperty<Viscosity>();

    // SPHinXsys cylinder wall path: TriangleMeshShapeCylinder + LevelSet + NormalFromBodyShapeCK
    // (same pattern as taylor_bar_sycl / french_natural_glass_relax).
    const int wall_mesh_resolution =
        french.glass_mesh_resolution > 0 ? french.glass_mesh_resolution : 20;
    SolidBody wall(sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                                   "WallBoundary", french, wall_thickness, wall_mesh_resolution, rim_height));
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    (void)wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> glass_inner(glass);
    Contact<> glass_wall_contact(glass, {&wall});
    Inner<> wall_inner(wall);
    Contact<> wall_glass_contact(wall, {&glass});

    BaseParticles &glass_particles = glass.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(glass_particles, t0);
    registerOphelieThermalDiffusionAuxFields(glass_particles, k_th);
    registerOphelieFrenchThermalRheologyFields(glass_particles, mu, cp, k_th);
    glass_particles.registerStateVariable<Real>("DensSumSigma0");
    glass_particles.addEvolvingVariable<Real>("DensSumSigma0");
    glass_particles.registerStateVariable<Real>("DensSumZ0");
    glass_particles.addEvolvingVariable<Real>("DensSumZ0");
    if (local_cli.layout_jitter_dp_frac > TinyReal)
    {
        hostApplyLayoutJitter(glass_particles, local_cli.layout_jitter_dp_frac * french.dp, local_cli.layout_seed);
    }
    const OphelieMaterialLawCheckResult material_check = verifyFrenchMaterialLawsAtReferenceTemperatures(laws);
    std::cout << "[ophelie][stage4.1] material_check=" << material_check.message
              << " mu_1473=" << material_check.mu_1473 << " cp_1473=" << material_check.cp_1473
              << " k_1473=" << material_check.k_1473 << " sigma_1473=" << material_check.sigma_1473
              << " passed=" << (material_check.passed ? 1 : 0) << std::endl;
    {
        std::string enthalpy_msg;
        const bool enthalpy_ok = verifyOphelieEnthalpyEnergyAccounting(enthalpy_msg);
        std::cout << "[ophelie][stage4.1] enthalpy_check=" << enthalpy_msg
                  << " passed=" << (enthalpy_ok ? 1 : 0) << std::endl;
        if (!enthalpy_ok)
        {
            std::cerr << "[ophelie][stage4.1] enthalpy accounting regression failed" << std::endl;
            return 1;
        }
        std::string applied_msg;
        const bool applied_ok = verifyFrenchAppliedEnergyResidualIdentity(applied_msg);
        std::cout << "[ophelie][stage4.1] applied_energy_check=" << applied_msg
                  << " passed=" << (applied_ok ? 1 : 0) << std::endl;
        if (!applied_ok)
        {
            std::cerr << "[ophelie][stage4.1] applied energy residual identity failed" << std::endl;
            return 1;
        }
        std::string stage1_msg;
        const bool stage1_ok = verifyFrenchThermalStage1HostAlgebra(stage1_msg);
        std::cout << "[ophelie][stage4.1] thermal_stage1_host=" << stage1_msg
                  << " passed=" << (stage1_ok ? 1 : 0) << std::endl;
        if (!stage1_ok)
        {
            std::cerr << "[ophelie][stage4.1] thermal stage-1 host identities failed" << std::endl;
            return 1;
        }
    }
    {
        std::string visc_msg;
        const bool visc_ok = verifyFrenchViscosityFunctorCachesByValue(visc_msg);
        std::cout << "[ophelie][stage4.1] viscosity_ck_cache_check=" << visc_msg
                  << " passed=" << (visc_ok ? 1 : 0) << std::endl;
        if (!visc_ok)
        {
            std::cerr << "[ophelie][stage4.1] viscosity CK cache regression failed" << std::endl;
            return 1;
        }
    }
    glass_particles.registerStateVariable<Real>("JouleHeat", Real(0));
    {
        Real *q_host = glass_particles.getVariableDataByName<Real>("JouleHeat");
        const size_t n = glass_particles.TotalRealParticles();
        if (n != frozen_q.size())
        {
            std::cerr << "[ophelie][stage4.1] particle count mismatch EM n=" << frozen_q.size() << " fluid n=" << n
                      << std::endl;
            return 1;
        }
        for (size_t i = 0; i < n; ++i)
        {
            q_host[i] = frozen_q[i];
        }
        syncVariableToDevice<Real>(glass_particles, "JouleHeat");
        if (periodic_coupling)
        {
            hostSampleEulerQOntoParticles(glass_particles, q_grid);
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
            writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
            appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                      local_cli.literature.em_control);
        }
        else
        {
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
        }
    }

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = true;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = local_cli.literature.enable_thermal_diffusion;
    thermal_bc.boundary_width_factor = params.phi_boundary_distance_factor_;
    thermal_bc.renormalize_boundary_area_to_geometry = local_cli.boundary_area_normalize;
    thermal_bc.h_side = Real(300);
    thermal_bc.h_bottom = Real(35);
    thermal_bc.h_free = Real(20);
    thermal_bc.emissivity = Real(0.8);
    thermal_bc.t_cool = Real(300);
    thermal_bc.t_ambient = Real(300);
    thermal_bc.t_rad_ambient = Real(300);
    {
        std::string skull_msg;
        if (!verifyFrenchSkullResistanceRegression(skull_msg))
        {
            std::cerr << "[ophelie][skull-resistance] regression failed: " << skull_msg << std::endl;
            return 1;
        }
        std::cout << "[ophelie][skull-resistance] " << skull_msg << std::endl;
        std::string skull_err;
        if (local_cli.skull_shell.enable)
        {
            // Level-2: melt side/bottom heat leaves via Contact to WallBoundary; zero Robin there.
            thermal_bc.h_side = Real(0);
            thermal_bc.h_bottom = Real(0);
            thermal_bc.t_cool = local_cli.skull_shell.t_cool;
            std::cout << "[ophelie][skull-shell] melt side/bottom Robin disabled; Contact + outer cool active\n";
        }
        else if (!applyFrenchSkullResistanceToThermalBc(thermal_bc, local_cli.skull, skull_err))
        {
            std::cerr << "[ophelie][stage4.1] skull setup failed: " << skull_err << std::endl;
            return 1;
        }
        {
            std::string crust_err;
            if (!applyFrenchFreeCrustToThermalBc(thermal_bc, local_cli.free_crust.enable, local_cli.free_crust.delta_m,
                                                local_cli.free_crust.k_crust, crust_err))
            {
                std::cerr << "[ophelie][stage4.1] free-crust setup failed: " << crust_err << std::endl;
                return 1;
            }
        }
        const Real melt_h = Real(2) * french.glass_half_height;
        const FrenchZeroDHeatBudget budget0d = estimateFrenchZeroDHeatBudget(
            french.glass_radius, melt_h, t0, thermal_bc, local_cli.target_power_w);
        const char *tag = "bare_walls";
        if (local_cli.skull_shell.enable && local_cli.free_crust.enable)
        {
            tag = "skull_shell_contact_free_crust";
        }
        else if (local_cli.skull_shell.enable)
        {
            tag = "skull_shell_contact";
        }
        else if (local_cli.skull.enable && local_cli.free_crust.enable)
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
    }
    {
        const Real melt_h = Real(2) * french.glass_half_height;
        // Proxy ΔT: 20 K until paper figure ΔT is wired; OPEN_PARAMETER for Round 3C.
        const Real dT_proxy = Real(20);
        const Real mu_hot = laws.mu.evaluate(t0 + dT_proxy);
        logFrenchNaturalRaPrEstimate(french.glass_radius, melt_h, dT_proxy, rho0, mu, cp, k_th, beta, local_cli.gravity_g,
                                     mu_hot, Real(0.02));
    }
    (void)setupOphelieThermalFrenchNaturalBoundaryFaces(glass_particles, thermal_bc, french);
    {
        const OphelieFrenchBoundaryAreaAudit audit =
            hostOphelieThermalFrenchNaturalBoundaryAreaAudit(glass_particles, thermal_bc, french,
                                                             thermal_bc.boundary_shell_thickness);
        logOphelieFrenchBoundaryAreaAudit(audit, "natural_vs_CAD");
    }

    // Level-2: register wall temperature / conductivity for Contact conduction.
    BaseParticles &wall_particles = wall.getBaseParticles();
    Real k_wall_eff = local_cli.skull_shell.k_skull;
    size_t n_wall_cool = 0;
    if (local_cli.skull_shell.enable)
    {
        registerOphelieJouleHeatTemperatureField(wall_particles, t0);
        k_wall_eff = frenchSkullShellEffectiveWallConductivity(local_cli.skull_shell.k_skull, wall_thickness,
                                                               local_cli.skull_shell.delta_m);
        registerOphelieThermalDiffusionAuxFields(wall_particles, k_wall_eff);
        hostInitializeFrenchSkullShellWallTemperatureProfile(wall_particles, french, wall_thickness, t0,
                                                             local_cli.skull_shell.t_cool, t0);
        n_wall_cool = hostApplyFrenchSkullShellOuterCoolDirichlet(wall_particles, french, wall_thickness,
                                                                  local_cli.skull_shell.t_cool, t0);
        std::cout << "[ophelie][skull-shell] wall particles n=" << wall_particles.TotalRealParticles()
                  << " k_eff=" << k_wall_eff << " (matches R=delta/k_skull) outer_cool_n=" << n_wall_cool
                  << " (linear T profile melt→cool)\n";
    }

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall).exec();

    auto &update_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(glass);
    auto &update_wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall);
    auto &update_glass_relations = main_methods.addRelationDynamics(glass_inner, glass_wall_contact);
    auto &update_wall_relations = main_methods.addRelationDynamics(wall_inner, wall_glass_contact);
    auto &fluid_particle_sorting = main_methods.addSortDynamics(glass);

    Gravity gravity(Vecd(0.0, 0.0, -local_cli.gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(glass, gravity);
    auto &boussinesq_force =
        main_methods.addStateDynamics<BoussinesqBuoyancyForceCK>(glass, gravity, beta, t0, kOphelieTemperatureField);
    auto &fluid_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(glass);
    auto &fluid_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(glass);

    auto &acoustic_step_1st_half =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(glass_wall_contact);
    auto &acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK,
                                                    NoKernelCorrectionCK>(glass_inner);
    auto &acoustic_step_2nd_half_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_wall_contact);
    acoustic_step_2nd_half.addPostContactInteraction(acoustic_step_2nd_half_wall);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(glass_inner, 0.5), glass_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        acoustic_1st_corr(glass_inner, glass_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        acoustic_2nd_corr(glass_inner, glass_wall_contact);
    // Doc105 P2e: separate fluid–fluid vs wall-contact acoustic_1st for force residual split.
    auto &acoustic_1st_inner_only =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_inner);
    auto &acoustic_1st_wall_only =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep1stHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep1stHalf<
                              Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>>>
        acoustic_1st_inner_corr_only(glass_inner);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep1stHalf<
                              Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>>>
        acoustic_1st_wall_corr_only(glass_wall_contact);
    auto &density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
            .addPostContactInteraction(glass_wall_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(glass);

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(
            glass_inner);
    auto &fluid_viscous_force_from_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(
            glass_wall_contact);
    fluid_viscous_force.addPostContactInteraction(fluid_viscous_force_from_wall);

    const Real acoustic_cfl_eff = local_cli.acoustic_cfl * local_cli.dt_scale;
    const Real advection_cfl_eff = local_cli.advection_cfl * local_cli.dt_scale;
    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(
        glass, local_cli.u_ref, advection_cfl_eff);
    auto &acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(glass, acoustic_cfl_eff);
    auto &u_max_reduce = main_methods.addReduceDynamics<OphelieVelocityMaxReduceCK>(glass);
    auto &t_mean_reduce = main_methods.addReduceDynamics<OphelieTemperatureMeanReduceCK>(glass);
    auto &t_max_reduce =
        main_methods.addReduceDynamics<OphelieThermalMaxTemperatureReduceCK>(glass, kOphelieTemperatureField);
    auto &t_min_reduce =
        main_methods.addReduceDynamics<OphelieThermalMinTemperatureReduceCK>(glass, kOphelieTemperatureField);

    StateDynamics<MainExecutionPolicy, ApplyFrozenQToTemperatureCK> apply_q_const(glass, rho0, cp, t0);
    StateDynamics<MainExecutionPolicy, ApplyJouleHeatWithParticleCpCK> apply_q_variable(glass, rho0, t0);
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc);
    ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK> boundary_loss_power(
        glass, t0, thermal_bc.boundary_shell_thickness, thermal_bc);
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>>> laplace_temperature;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>> apply_diffusion;
    // Level-2 shell: RH200-style Inner+Contact diffusion (replaces Ophelie melt Laplace when enabled).
    using SkullShellThermalRelaxation = RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                              Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                              Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>>>;
    UniquePtr<IsotropicDiffusion> glass_heat_diffusion;
    UniquePtr<IsotropicDiffusion> heat_diffusion_glass_wall;
    UniquePtr<IsotropicDiffusion> wall_heat_diffusion;
    UniquePtr<IsotropicDiffusion> heat_diffusion_wall_glass;
    UniquePtr<SkullShellThermalRelaxation> glass_shell_thermal;
    UniquePtr<SkullShellThermalRelaxation> wall_shell_thermal;
    const bool use_skull_shell_thermal = local_cli.skull_shell.enable;
    if (use_skull_shell_thermal)
    {
        const Real cv_glass = rho0 * cp;
        const Real cv_wall = local_cli.skull_shell.rho_skull * local_cli.skull_shell.cp_skull;
        const Real k_contact = Real(2) * k_th * k_wall_eff / (k_th + k_wall_eff + TinyReal);
        glass_heat_diffusion = makeUnique<IsotropicDiffusion>(kOphelieTemperatureField, k_th, cv_glass);
        heat_diffusion_glass_wall = makeUnique<IsotropicDiffusion>(kOphelieTemperatureField, k_contact, cv_glass);
        wall_heat_diffusion = makeUnique<IsotropicDiffusion>(kOphelieTemperatureField, k_wall_eff, cv_wall);
        heat_diffusion_wall_glass = makeUnique<IsotropicDiffusion>(kOphelieTemperatureField, k_contact, cv_wall);
        glass_shell_thermal = makeUnique<SkullShellThermalRelaxation>(
            DynamicsArgs(glass_inner, glass_heat_diffusion.get()),
            DynamicsArgs(glass_wall_contact, heat_diffusion_glass_wall.get()));
        wall_shell_thermal = makeUnique<SkullShellThermalRelaxation>(
            DynamicsArgs(wall_inner, wall_heat_diffusion.get()),
            DynamicsArgs(wall_glass_contact, heat_diffusion_wall_glass.get()));
        std::cout << "[ophelie][skull-shell] Contact thermal on: k_melt=" << k_th << " k_wall_eff=" << k_wall_eff
                  << " k_contact=" << k_contact << " cv_wall=" << cv_wall << "\n";
    }
    else if (thermal_bc.enable_diffusion)
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
    glass_state_recorder.addToWrite<Real>(glass, "JouleHeat");
    const bool write_vtp = local_cli.state_record_interval > TinyReal;

    const FrenchCylindricalFrame cyl_frame = makeFrenchCylinderFrameFromGlass(
        french.glass_center, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
    OphelieThermalRheologyStats rheology_stats;
    const Real mu_init = hostUpdateFrenchThermalRheologyFromTemperature(
        glass_particles, laws, local_cli.literature.constant_thermal_rheology, t0, rheology_stats);
    glass_viscosity.setReferenceViscosity(mu_init);
    invalidateViscosityDependentCkKernels(glass_inner, glass_wall_contact);
    if (!local_cli.literature.constant_thermal_rheology)
    {
        std::cout << "[ophelie][stage4.1] viscosity_mode=mean_mu_T_volume_average mu_ref=" << mu_init
                  << std::endl;
    }
    Real prev_energy = local_cli.literature.constant_thermal_rheology
                           ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                           : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
    Real prev_energy_time = 0.0;
    coupling_state.ledger.setH0(prev_energy);
    if (local_cli.freeze_prescribed_temperature)
    {
        coupling_state.ledger.noteLiveLoss(Real(0), Real(0), Real(0));
        std::cout << "[ophelie][stage4.1] freeze_prescribed_t=1 disables Q, diffusion, Robin/radiation "
                     "and skull-shell thermal updates (doc 40 R10 / doc 41)\n";
    }
    else
    {
        Real side0 = 0, bottom0 = 0, free_conv0 = 0, free_rad0 = 0, total0 = 0;
        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                      side0, bottom0, free_conv0, free_rad0, total0);
        coupling_state.ledger.noteLiveLoss(side0, bottom0, free_conv0 + free_rad0);
    }

    update_wall_cell_linked_list.exec();
    update_glass_cell_linked_list.exec();
    update_glass_relations.exec();
    update_wall_relations.exec();
    density_regularization.exec();
    const Real z_bottom = gc[2] - french.glass_half_height;
    const Real z_free = gc[2] + french.glass_half_height;
    const Real z_rim = z_free + rim_height;
    hostApplyHydrostaticDensity(glass_particles, rho0, local_cli.c0, local_cli.gravity_g, z_free);
    fluid_advection_step_setup.exec();
    constant_gravity.exec();

    double dens_sigma0_mean = 0.0;
    bool dens_sigma0_ready = false;
    bool dens_audit_header = true;
    FrenchDensSumAuditSnapshot dens_audit_last{};
    const Real dens_audit_surface_depth = Real(2) * french.dp;
    const std::string dens_audit_path = "output/french_dens_sum_audit.csv";
    size_t dens_sum_n_written_last = 0;
    size_t dens_sum_n_recaptured_last = 0;
    FrenchDensWriteAudit dens_write_last{};
    FrenchFsAtmAudit fs_atm_last{};
    size_t dens_fs_cap_n_last = 0;
    size_t advection_steps = 0;
    bool in_hydro_settle = true;
    auto apply_quiet_hydro_after_dens = [&]() {
        const bool do_dens_norm =
            local_cli.dens_sum_normalize || (in_hydro_settle && local_cli.settle_dens_sum_normalize);
        if (do_dens_norm)
        {
            if (!dens_sigma0_ready)
            {
                hostCaptureDensSumSigma0Field(glass_particles, dens_sigma0_mean);
                dens_sigma0_ready = true;
                std::cout << "[ophelie][doc95] dens_sum_normalize sigma0_mean=" << dens_sigma0_mean
                          << " packing=" << (dens_sigma0_mean / (rho0 + TinyReal))
                          << " interior_only=" << (local_cli.dens_sum_interior_only ? 1 : 0)
                          << " settle_only="
                          << ((in_hydro_settle && local_cli.settle_dens_sum_normalize &&
                               !local_cli.dens_sum_normalize)
                                  ? 1
                                  : 0)
                          << " surface_recapture_interval="
                          << local_cli.dens_sum_surface_recapture_interval
                          << " surface_depth=" << dens_audit_surface_depth << std::endl;
                std::cout << "[ophelie][doc97-p4] dens_sum_semantics: FreeSurface sets rho=max(rho_sum,rho0); "
                             "extra_normalize then overwrites rho=rho0*rho_sum/sigma0"
                          << (local_cli.dens_sum_interior_only
                                  ? " except surface band (interior_only keeps FreeSurface rho)"
                                  : " (full domain)")
                          << "; DensSumSigma0 is evolving (moves with particle sorting)\n";
                std::cout << "[ophelie][doc105-p2] no_extra_normalize="
                          << (do_dens_norm ? 0 : 1)
                          << " dens_fs_atm_consistent=" << (local_cli.dens_fs_atm_consistent ? 1 : 0)
                          << " (DensitySummation+FreeSurface still always run when this path executes)\n";
            }
            else if (local_cli.dens_sum_surface_recapture_interval > 0 &&
                     advection_steps % static_cast<size_t>(local_cli.dens_sum_surface_recapture_interval) == 0)
            {
                dens_sum_n_recaptured_last =
                    hostRecaptureSurfaceSigma0(glass_particles, z_free, dens_audit_surface_depth);
            }
            dens_write_last = hostApplyDensSumNormalizeFieldAudited(
                glass_particles, rho0, local_cli.c0, local_cli.dens_sum_interior_only, z_free,
                dens_audit_surface_depth);
            dens_sum_n_written_last = dens_write_last.n_written;
        }
        // Doc105 P2i: cap overpack only when FreeSurface path is active (no extra normalize write).
        dens_fs_cap_n_last = 0;
        if (local_cli.dens_fs_cap_overpack && !do_dens_norm)
        {
            dens_fs_cap_n_last = hostCapOverpackDensityToRho0(glass_particles, rho0);
        }
        if (local_cli.background_hydro_p)
            hostApplyFrenchBackgroundHydroPressure(glass_particles, rho0, local_cli.gravity_g, z_free);
        // Doc105 P2b: after p_h write, restore atm condition on incomplete-support particles.
        if (local_cli.dens_fs_atm_consistent)
            fs_atm_last = hostApplyFrenchFsAtmConsistentDensity(glass_particles, rho0, local_cli.c0, z_free,
                                                                dens_audit_surface_depth);
        // Note: hydro_rho_freeze_rho0 is applied only in stage/keep probes (not settle),
        // so layout matches the unfrozen baseline.
    };
    auto maybe_dens_sum_audit = [&](Real t_now, size_t adv_step) {
        if (!(local_cli.dens_sum_audit && local_cli.dens_sum_normalize && dens_sigma0_ready))
            return;
        if (adv_step % static_cast<size_t>(local_cli.dens_sum_audit_interval) != 0)
            return;
        dens_audit_last = hostFrenchDensSumAudit(glass_particles, rho0, z_free, dens_audit_surface_depth);
        appendFrenchDensSumAuditCsv(dens_audit_path, t_now, adv_step, dens_audit_last, dens_audit_header);
        dens_audit_header = false;
    };
    size_t incomplete_force_cancel_n_last = 0;
    size_t fs_pack_scale_n_last = 0;
    auto run_acoustic = [&](Real dt_ac) {
        // Doc105 P2n/P2r: never rewrite Force during settle — packing depends on gravity kick.
        const bool do_incomplete_fcancel =
            local_cli.incomplete_force_cancel_live && !in_hydro_settle;
        const bool do_fs_pack_scale =
            local_cli.fs_incomplete_pack_scale_force && !in_hydro_settle;
        if (local_cli.kernel_correction)
        {
            acoustic_1st_corr.exec(dt_ac);
            if (do_incomplete_fcancel)
                incomplete_force_cancel_n_last = hostUndoIncompleteAcoustic1stKick(
                    glass_particles, rho0, dt_ac, z_free, dens_audit_surface_depth);
            else if (do_fs_pack_scale)
                fs_pack_scale_n_last = hostScaleIncompleteFsForceByPacking(
                    glass_particles, rho0, dt_ac, z_free, dens_audit_surface_depth, true);
            acoustic_2nd_corr.exec(dt_ac);
        }
        else
        {
            acoustic_step_1st_half.exec(dt_ac);
            if (do_incomplete_fcancel)
                incomplete_force_cancel_n_last = hostUndoIncompleteAcoustic1stKick(
                    glass_particles, rho0, dt_ac, z_free, dens_audit_surface_depth);
            else if (do_fs_pack_scale)
                fs_pack_scale_n_last = hostScaleIncompleteFsForceByPacking(
                    glass_particles, rho0, dt_ac, z_free, dens_audit_surface_depth, true);
            acoustic_step_2nd_half.exec(dt_ac);
        }
    };
    if (local_cli.kernel_correction)
        fluid_linear_correction_matrix.exec();
    if (local_cli.background_hydro_p)
        hostApplyFrenchBackgroundHydroPressure(glass_particles, rho0, local_cli.gravity_g, z_free);

    // Doc105 P2d: French wall dens mode. Momentum stays total-p (PerturbationMomentum default 0).
    // mode 1 ⇒ a_ref = −gẑ because GravityForceCK keeps g in ForcePrior (≠ closed-box P2b).
    {
        hostSetAcousticDensModeFrench(glass_particles, static_cast<Real>(local_cli.acoustic_dens_mode));
        const Vecd a_ref = (local_cli.acoustic_dens_mode == 1) ? Vecd(0.0, 0.0, -local_cli.gravity_g) : Vecd::Zero();
        hostSetReferenceBodyAccelerationFrench(glass_particles, a_ref);
        std::cout << "[ophelie][doc105-p2d] acoustic_dens_mode=" << local_cli.acoustic_dens_mode
                  << " reference_body_acc_z=" << a_ref[2]
                  << " perturbation_momentum=0 french_pressure_path=total_p_plus_GravityForceCK"
                  << " closed_box_p2b_copied=0\n";
    }

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording =
        time_stepper.addTriggerByInterval(std::max(local_cli.state_record_interval, Real(1.0e-6)));

    const Real glass_h = frenchReducedGlassHeight(french);
    const Real alpha_th = k_th / (rho0 * cp + TinyReal);
    const Real h_min_ts = glass.getSPHAdaptation().MinimumSmoothingLength();
    const Real nu_ts = mu / (rho0 + TinyReal);
    const Real u_visc_ref = nu_ts / (h_min_ts + TinyReal);
    const Real dt_ac_est = acoustic_cfl_eff * h_min_ts / (local_cli.c0 + TinyReal);
    const Real dt_adv_inertial = advection_cfl_eff * h_min_ts / (local_cli.u_ref + TinyReal);
    const Real dt_adv_viscous =
        advection_cfl_eff * h_min_ts / (std::max(u_visc_ref, local_cli.u_ref) + TinyReal);
    std::cout << "[ophelie][stage4.1] hydro: c0=" << local_cli.c0 << " u_ref=" << local_cli.u_ref
              << " rim=" << rim_height << " wall_thickness=" << wall_thickness
              << " settle_s=" << local_cli.hydro_settle_time
              << " gH/c0^2=" << (local_cli.gravity_g * glass_h) / (local_cli.c0 * local_cli.c0 + TinyReal)
              << " k=" << k_th << " alpha=" << alpha_th
              << " dens_sum_normalize=" << (local_cli.dens_sum_normalize ? 1 : 0)
              << " settle_dens_sum_normalize=" << (local_cli.settle_dens_sum_normalize ? 1 : 0)
              << " background_hydro_p=" << (local_cli.background_hydro_p ? 1 : 0)
              << " kernel_correction=" << (local_cli.kernel_correction ? 1 : 0)
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0)
              << " dt_scale=" << local_cli.dt_scale << " acoustic_cfl=" << acoustic_cfl_eff
              << " advection_cfl=" << advection_cfl_eff << std::endl;
    std::cout << "[ophelie][doc97-p4] timestep_audit h_min=" << h_min_ts << " nu=" << nu_ts
              << " u_visc_ref=" << u_visc_ref << " dt_acoustic_est=" << dt_ac_est
              << " dt_adv_inertial=" << dt_adv_inertial << " dt_adv_with_viscous_ref=" << dt_adv_viscous
              << " viscous_tightens_advection=" << (u_visc_ref > local_cli.u_ref + TinyReal ? 1 : 0)
              << " (French still uses AdvectionTimeStepCK; AdvectionViscousTimeStepCK not enabled)\n";

    std::ofstream monitor("output/french_natural_convection_monitor.csv");
    monitor << "advection_step,physical_time_s,U_max,T_min,T_mean,T_max,dT,U_buoy,escaped,escaped_frac,"
               "P_joule,dE_dt,energy_residual,energy_residual_rel,total_heat_loss,em_updates,"
               "U_rms,U_theta_rms,U_z_rms,T_std,T_outer,T_center,T_bottom,T_top\n";
    writeFrenchSpatialStatsCsvHeader("output/french_spatial_stats.csv");

    while (time_stepper.getPhysicalTime() + TinyReal < local_cli.hydro_settle_time)
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        run_acoustic(acoustic_dt);
        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            if (advection_steps % 50 == 0)
            {
                fluid_particle_sorting.exec();
            }
            update_glass_cell_linked_list.exec();
            update_glass_relations.exec();
            density_regularization.exec();
            apply_quiet_hydro_after_dens();
            maybe_dens_sum_audit(time_stepper.getPhysicalTime(), advection_steps);
            fluid_advection_step_setup.exec();
            if (local_cli.kernel_correction)
                fluid_linear_correction_matrix.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
        }
    }
    if (local_cli.hydro_settle_time > TinyReal)
    {
        const Real u_settle = u_max_reduce.exec();
        const GlassEscapeCensus settle_escape = hostCountEscapedGlassParticles(
            glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
        std::cout << "[ophelie][stage4.1] hydro settle done: t=" << time_stepper.getPhysicalTime()
                  << " U_max=" << u_settle << " escaped=" << settle_escape.escaped << "/" << settle_escape.total
                  << " r_max=" << settle_escape.r_max << " z=[" << settle_escape.z_min << "," << settle_escape.z_max
                  << "]" << std::endl;
        hostZeroVelocity(glass_particles);
        time_stepper.setPhysicalTime(0.0);
        advection_steps = 0;
        // Doc105 P2h: settle-only normalize ends here; measurement uses FreeSurface path.
        if (local_cli.settle_dens_sum_normalize && !local_cli.dens_sum_normalize)
        {
            dens_sigma0_ready = false;
            dens_sigma0_mean = 0.0;
            dens_sum_n_written_last = 0;
            dens_sum_n_recaptured_last = 0;
            dens_write_last = {};
            std::cout << "[ophelie][doc105-p2h] settle_dens_sum_normalize off for probe/keep "
                         "(sigma0 cleared; FreeSurface dens path)\n";
        }
        in_hydro_settle = false;
    }
    else
    {
        in_hydro_settle = false;
    }

    // Doc105 P2: one hydrostatic stage probe with zone residuals (French free-surface path).
    // Pressure convention stays total-p + GravityForceCK (PerturbationMomentum default 0).
    // Not closed-box P2b. Prefer --zero-buoyancy --end-time=0 for a pure quiet keep probe.
    if (local_cli.doc105_p2_hydro_stage_probe)
    {
        const Real zone_band = Real(2) * french.dp;
        const Real deep_band = Real(4) * french.dp;
        const Vecd glass_c = french.glass_center;
        auto snap = [&](const char *stage) {
            const FrenchHydroStageSnapshot st = hostFrenchHydroStageSnapshot(
                glass_particles, rho0, local_cli.c0, local_cli.gravity_g, glass_c, french.glass_radius, z_bottom,
                z_free, zone_band, deep_band);
            std::cout << "[ophelie][doc105-p2][" << stage << "] n_all=" << st.all.n
                      << " n_int=" << st.interior.n << " n_deep=" << st.deep_interior.n
                      << " n_nfs=" << st.near_fs_interior.n << " n_surf=" << st.surface.n
                      << " n_side=" << st.wall_side.n << " n_bot=" << st.wall_bottom.n
                      << " n_corner=" << st.wall_corner.n << " n_ws=" << st.wall_surface.n
                      << " a_net_int=" << st.interior.a_net_z_rms << " a_net_deep=" << st.deep_interior.a_net_z_rms
                      << " a_net_nfs=" << st.near_fs_interior.a_net_z_rms
                      << " a_net_surf=" << st.surface.a_net_z_rms
                      << " a_net_side=" << st.wall_side.a_net_z_rms
                      << " a_net_r_side=" << st.wall_side.a_net_r_rms
                      << " a_net_bot=" << st.wall_bottom.a_net_z_rms << " a_net_ws=" << st.wall_surface.a_net_z_rms
                      << " a_f_int=" << st.interior.a_force_z_rms << " a_f_surf=" << st.surface.a_force_z_rms
                      << " a_f_side=" << st.wall_side.a_force_z_rms << " a_f_bot=" << st.wall_bottom.a_force_z_rms
                      << " pack_int=" << st.interior.packing_mean << " pack_deep=" << st.deep_interior.packing_mean
                      << " pack_surf=" << st.surface.packing_mean
                      << " n_inc_surf=" << st.surface.n_incomplete << " dp_h_int=" << st.interior.dp_h_rms
                      << " dp_h_surf=" << st.surface.dp_h_rms << " p_eos_int=" << st.interior.p_eos_rms
                      << " drho_all=" << st.all.drho_dt_rms << std::endl
                      << std::flush;
            return st;
        };

        hostZeroVecFieldFrench(glass_particles, "Velocity");
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroVecFieldFrench(glass_particles, "ForcePrior");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        FrenchHydroStageSnapshot s0 = snap("00_init_zero");

        density_regularization.exec();
        apply_quiet_hydro_after_dens();
        // Doc105 P2f: freeze ρ=ρ0 only for probe (after settle), same layout as baseline.
        if (local_cli.hydro_rho_freeze_rho0)
        {
            hostFreezeDensityToRho0(glass_particles, rho0);
            std::cout << "[ophelie][doc105-p2f] hydro_rho_freeze_rho0 applied after dens (probe only)\n";
        }
        fluid_advection_step_setup.exec();
        if (local_cli.kernel_correction)
            fluid_linear_correction_matrix.exec();
        FrenchHydroStageSnapshot s_dens = snap("01_after_dens_bgp");

        hostZeroVecFieldFrench(glass_particles, "ForcePrior");
        hostZeroVecFieldFrench(glass_particles, "GravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousGravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "BoussinesqBuoyancyForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
        constant_gravity.exec();
        // Intentionally no boussinesq here — probe isolates gravity/pressure path.
        FrenchHydroStageSnapshot s_prior = snap("02_after_gravity");

        const Real dt_ac = acoustic_time_step.exec();

        // Doc126 P0: capture pre-acoustic state so ff/wall/full budgets share the same input
        // and destructive Update kicks can be restored (stage_04 is NOT a single normal step).
        FrenchHydroFieldSnapshot probe_base;
        hostCaptureFrenchHydroFields(glass_particles, probe_base);

        // Doc105 P2e: fluid–fluid only (Inner OneLevel). Force after Update still holds ∇p_ff.
        hostRestoreFrenchHydroFields(glass_particles, probe_base);
        hostZeroVecFieldFrench(glass_particles, "Velocity");
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        if (local_cli.kernel_correction)
            acoustic_1st_inner_corr_only.exec(dt_ac);
        else
            acoustic_1st_inner_only.exec(dt_ac);
        FrenchHydroStageSnapshot s_ff = snap("03a_ff_only");

        // Wall-contact only from the same pre-acoustic snapshot.
        hostRestoreFrenchHydroFields(glass_particles, probe_base);
        hostZeroVecFieldFrench(glass_particles, "Velocity");
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        if (local_cli.kernel_correction)
            acoustic_1st_wall_corr_only.exec(dt_ac);
        else
            acoustic_1st_wall_only.exec(dt_ac);
        FrenchHydroStageSnapshot s_wc = snap("03b_wall_only");

        // Full acoustic_1st (ff + wall) — canonical residual used for pass/fail.
        hostRestoreFrenchHydroFields(glass_particles, probe_base);
        hostZeroVecFieldFrench(glass_particles, "Velocity");
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        if (local_cli.kernel_correction)
            acoustic_1st_corr.exec(dt_ac);
        else
            acoustic_step_1st_half.exec(dt_ac);
        FrenchHydroStageSnapshot s_ac = snap("03_after_acoustic_1st");
        FrenchHydroFieldSnapshot after_ac1;
        hostCaptureFrenchHydroFields(glass_particles, after_ac1);

        // Doc105 P2l: bound incomplete-support share of residual (Force-only; no vel rewrite).
        const size_t n_incomplete_cancel = hostCancelNetForceIncompleteSupport(glass_particles, rho0);
        FrenchHydroStageSnapshot s_ac_inc_cancel = snap("03c_incomplete_force_cancel");
        std::cout << "[ophelie][doc105-p2l] incomplete_force_cancel n=" << n_incomplete_cancel
                  << " a_surf=" << s_ac_inc_cancel.surface.a_net_z_rms
                  << " a_surf_inc=" << s_ac_inc_cancel.surface_incomplete.a_net_z_rms
                  << " a_surf_cmp=" << s_ac_inc_cancel.surface_complete.a_net_z_rms
                  << " a_deep=" << s_ac_inc_cancel.deep_interior.a_net_z_rms << std::endl
                  << std::flush;
        // Doc105 P2r: Shepard packing scale snap from restored acoustic_1st state.
        hostRestoreFrenchHydroFields(glass_particles, after_ac1);
        const size_t n_fs_pack_scale_stage = hostScaleIncompleteFsForceByPacking(
            glass_particles, rho0, dt_ac, z_free, dens_audit_surface_depth, false);
        FrenchHydroStageSnapshot s_ac_pack_scale = snap("03d_fs_incomplete_pack_scale");
        std::cout << "[ophelie][doc105-p2r] fs_incomplete_pack_scale n=" << n_fs_pack_scale_stage
                  << " a_surf=" << s_ac_pack_scale.surface.a_net_z_rms
                  << " a_surf_inc=" << s_ac_pack_scale.surface_incomplete.a_net_z_rms
                  << " a_deep=" << s_ac_pack_scale.deep_interior.a_net_z_rms << std::endl
                  << std::flush;
        // stage_04: continue from acoustic_1st (not a single normal time step — Doc126).
        hostRestoreFrenchHydroFields(glass_particles, after_ac1);
        if (local_cli.kernel_correction)
            acoustic_2nd_corr.exec(dt_ac);
        else
            acoustic_step_2nd_half.exec(dt_ac);
        FrenchHydroStageSnapshot s_ac2 = snap("04_after_acoustic_2nd");
        // Restore production state to pre-probe (diagnostics must not pollute subsequent keep).
        hostRestoreFrenchHydroFields(glass_particles, probe_base);

        const Real a_ref = std::abs(local_cli.gravity_g) + TinyReal;
        const Real u_hydro = std::sqrt(a_ref * glass_h + TinyReal);
        // Pre-declared Doc105 P2c pollution budgets (mechanism screen vs signal-safe).
        const Real budget_a_net_mech = Real(0.25) * a_ref;   // coarse screen
        const Real budget_a_net_quiet = Real(0.05) * a_ref;  // quiet target
        const Real budget_u_rms_mech = Real(5.0e-4);         // absolute, lit-β order
        const Real budget_u_rms_signal = Real(1.0e-4);       // << typical dUz~1e-4
        const bool dens_quiet = s_ac.all.drho_dt_rms < Real(1.0e-3) * rho0 * local_cli.c0 / (glass_h + TinyReal);
        const bool interior_mom_quiet = s_ac.interior.a_net_z_rms < budget_a_net_quiet;
        const bool interior_mom_mech = s_ac.interior.a_net_z_rms < budget_a_net_mech;
        const bool side_radial_mech = s_ac.wall_side.a_net_r_rms < budget_a_net_mech;
        const bool provisional_inner_pass = dens_quiet && interior_mom_quiet;
        const Real ff_dom_int =
            s_ff.interior.a_force_z_rms / (s_wc.interior.a_force_z_rms + TinyReal);
        const Real ff_dom_surf =
            s_ff.surface.a_force_z_rms / (s_wc.surface.a_force_z_rms + TinyReal);
        const Real wall_dom_side =
            s_wc.wall_side.a_force_z_rms / (s_ff.wall_side.a_force_z_rms + TinyReal);
        // Doc105 P2q: bottom vertical wall/ff force ratio (P2e style; RMS not additive).
        const Real wall_dom_bot =
            s_wc.wall_bottom.a_force_z_rms / (s_ff.wall_bottom.a_force_z_rms + TinyReal);
        const Real wall_dom_corner =
            s_wc.wall_corner.a_force_z_rms / (s_ff.wall_corner.a_force_z_rms + TinyReal);

        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"french_doc105_p2_hydro_stage\",\n"
               << "  \"doc105_p2\": 1,\n"
               << "  \"doc105_p2c\": 1,\n"
               << "  \"doc105_p2e\": 1,\n"
               << "  \"doc105_p2f\": 1,\n"
               << "  \"doc105_p2g\": 1,\n"
               << "  \"doc105_p2h\": 1,\n"
               << "  \"doc105_p2i\": 1,\n"
               << "  \"doc105_p2l\": 1,\n"
               << "  \"doc105_p2q\": 1,\n"
               << "  \"doc105_p2r\": 1,\n"
               << "  \"fs_incomplete_pack_scale_force\": "
               << (local_cli.fs_incomplete_pack_scale_force ? 1 : 0) << ",\n"
               << "  \"hydro_rho_freeze_rho0\": " << (local_cli.hydro_rho_freeze_rho0 ? 1 : 0) << ",\n"
               << "  \"dp\": " << french.dp << ",\n"
               << "  \"c0\": " << local_cli.c0 << ",\n"
               << "  \"g\": " << local_cli.gravity_g << ",\n"
               << "  \"dt_ac\": " << dt_ac << ",\n"
               << "  \"zone_band\": " << zone_band << ",\n"
               << "  \"deep_band\": " << deep_band << ",\n"
               << "  \"budget_a_net_mech\": " << budget_a_net_mech << ",\n"
               << "  \"budget_a_net_quiet\": " << budget_a_net_quiet << ",\n"
               << "  \"budget_u_rms_mech\": " << budget_u_rms_mech << ",\n"
               << "  \"budget_u_rms_signal\": " << budget_u_rms_signal << ",\n"
               << "  \"u_hydro\": " << u_hydro << ",\n"
               << "  \"dens_sum_normalize\": " << (local_cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"settle_dens_sum_normalize\": " << (local_cli.settle_dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"dens_fs_cap_overpack\": " << (local_cli.dens_fs_cap_overpack ? 1 : 0) << ",\n"
               << "  \"dens_fs_cap_n\": " << dens_fs_cap_n_last << ",\n"
               << "  \"no_extra_normalize\": "
               << ((local_cli.dens_sum_normalize || local_cli.settle_dens_sum_normalize) ? 0 : 1) << ",\n"
               << "  \"probe_extra_normalize\": " << (local_cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"dens_sum_interior_only\": " << (local_cli.dens_sum_interior_only ? 1 : 0) << ",\n"
               << "  \"dens_fs_atm_consistent\": " << (local_cli.dens_fs_atm_consistent ? 1 : 0) << ",\n"
               << "  \"background_hydro_p\": " << (local_cli.background_hydro_p ? 1 : 0) << ",\n"
               << "  \"kernel_correction\": " << (local_cli.kernel_correction ? 1 : 0) << ",\n"
               << "  \"acoustic_dens_mode\": " << local_cli.acoustic_dens_mode << ",\n"
               << "  \"reference_body_acc_z\": "
               << (local_cli.acoustic_dens_mode == 1 ? -local_cli.gravity_g : Real(0)) << ",\n"
               << "  \"french_pressure_path\": \"total_p_plus_GravityForceCK\",\n"
               << "  \"perturbation_momentum\": 0,\n"
               << "  \"closed_box_p2b_copied\": 0,\n"
               << "  \"free_surface_regularization\": \"rho=max(rho_sum,rho0)\",\n"
               << "  \"extra_normalize\": \"rho=rho0*rho_sum/sigma0_overwrites_FreeSurface\",\n"
               << "  \"fs_atm_rule\": \"surface_band && rho_sum<rho0 => rho=rho0-p_h/c0^2 so p_total~0\",\n"
               << "  \"settle_reset\": {\"velocity_zeroed\": 1, \"physical_time_zeroed\": 1, "
                  "\"ForcePrior_history_not_cleared_at_settle\": 1, "
                  "\"rho_drho_Displacement_not_cleared_at_settle\": 1},\n"
               << "  \"dens_write_n\": " << dens_write_last.n_written << ",\n"
               << "  \"dens_write_n_skipped_surface\": " << dens_write_last.n_skipped_surface << ",\n"
               << "  \"dens_write_max_abs_drho\": " << dens_write_last.max_abs_drho << ",\n"
               << "  \"dens_write_rms_drho\": " << dens_write_last.rms_drho << ",\n"
               << "  \"dens_write_max_abs_dp\": " << dens_write_last.max_abs_dp << ",\n"
               << "  \"dens_write_rms_dp\": " << dens_write_last.rms_dp << ",\n"
               << "  \"dens_write_mass_sum\": " << dens_write_last.mass_sum << ",\n"
               << "  \"fs_atm_n_candidates\": " << fs_atm_last.n_candidates << ",\n"
               << "  \"fs_atm_n_written\": " << fs_atm_last.n_written << ",\n"
               << "  \"fs_atm_n_skipped_non_surface\": " << fs_atm_last.n_skipped_non_surface << ",\n"
               << "  \"fs_atm_rms_drho\": " << fs_atm_last.rms_drho << ",\n"
               << "  \"fs_atm_rms_p_total_after\": " << fs_atm_last.rms_p_total_after << ",\n"
               << "  \"fs_atm_max_abs_p_total_after\": " << fs_atm_last.max_abs_p_total_after << ",\n"
               << "  \"a_prior_z_rms\": " << s_prior.all.a_prior_z_rms << ",\n"
               << "  \"provisional_inner_pass\": " << (provisional_inner_pass ? 1 : 0) << ",\n"
               << "  \"interior_mom_mech_pass\": " << (interior_mom_mech ? 1 : 0) << ",\n"
               << "  \"side_radial_mech_pass\": " << (side_radial_mech ? 1 : 0) << ",\n"
               << "  \"dens_quiet\": " << (dens_quiet ? 1 : 0) << ",\n"
               << "  \"interior_mom_quiet\": " << (interior_mom_quiet ? 1 : 0) << ",\n"
               << "  \"p2e_ff_over_wall_a_force_interior\": " << ff_dom_int << ",\n"
               << "  \"p2e_ff_over_wall_a_force_surface\": " << ff_dom_surf << ",\n"
               << "  \"p2e_wall_over_ff_a_force_side\": " << wall_dom_side << ",\n"
               << "  \"p2e_wall_over_ff_a_force_bottom\": " << wall_dom_bot << ",\n"
               << "  \"p2e_wall_over_ff_a_force_corner\": " << wall_dom_corner << ",\n"
               << "  \"deep_interior_mom_quiet\": "
               << (s_ac.deep_interior.a_net_z_rms < budget_a_net_quiet ? 1 : 0) << ",\n"
               << "  \"deep_interior_mom_mech\": "
               << (s_ac.deep_interior.a_net_z_rms < budget_a_net_mech ? 1 : 0) << ",\n";
            hostWriteFrenchZoneJson(js, "stage_00_all", s0.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_01_all", s_dens.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_02_all", s_prior.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_all", s_ff.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_interior", s_ff.interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_surface", s_ff.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_wall_side", s_ff.wall_side);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_wall_bottom", s_ff.wall_bottom);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_wall_corner", s_ff.wall_corner);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03a_ff_wall_surface", s_ff.wall_surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_all", s_wc.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_interior", s_wc.interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_surface", s_wc.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_wall_side", s_wc.wall_side);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_wall_bottom", s_wc.wall_bottom);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_wall_corner", s_wc.wall_corner);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03b_wall_wall_surface", s_wc.wall_surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_all", s_ac.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_interior", s_ac.interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_deep_interior", s_ac.deep_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_near_fs_interior", s_ac.near_fs_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_surface", s_ac.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_surface_incomplete", s_ac.surface_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_surface_complete", s_ac.surface_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall", s_ac.wall);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_side", s_ac.wall_side);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_side_incomplete", s_ac.wall_side_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_side_complete", s_ac.wall_side_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_bottom", s_ac.wall_bottom);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_bottom_incomplete", s_ac.wall_bottom_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_bottom_complete", s_ac.wall_bottom_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_corner", s_ac.wall_corner);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03_wall_surface", s_ac.wall_surface);
            js << ",\n";
            js << "  \"incomplete_force_cancel_n\": " << n_incomplete_cancel << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03c_surface", s_ac_inc_cancel.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03c_surface_incomplete", s_ac_inc_cancel.surface_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03c_surface_complete", s_ac_inc_cancel.surface_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03c_deep_interior", s_ac_inc_cancel.deep_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03c_interior", s_ac_inc_cancel.interior);
            js << ",\n";
            js << "  \"fs_incomplete_pack_scale_n_stage\": " << n_fs_pack_scale_stage << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03d_surface", s_ac_pack_scale.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03d_surface_incomplete", s_ac_pack_scale.surface_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03d_deep_interior", s_ac_pack_scale.deep_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_03d_interior", s_ac_pack_scale.interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "stage_04_all", s_ac2.all);
            js << ",\n"
               << "  \"stage_04_not_single_normal_step\": 1,\n"
               << "  \"stage_04_note\": \"acoustic_1st then acoustic_2nd from restored ac1; not integration acceptance\",\n"
               << "  \"doc126_p0_snapshot_restore\": 1,\n"
               << "  \"eos_identity_ok_stage03\": "
               << ((s_ac.all.p_eos_rms + TinyReal >= std::abs(s_ac.all.p_eos_mean)) ? 1 : 0) << ",\n"
               << "  \"component_budget_note\": \"use a_force_r_rms/a_net_r_rms vs a_*_z_rms; RMS not subtractable\",\n"
               << "  \"passed\": " << (provisional_inner_pass ? 1 : 0) << ",\n"
               << "  \"note\": \"P2l/P2q/P2r: surface/bottom incomplete; 03c=cancel; "
                  "03d=Force*=rho0/rho_sum on FS incomplete; P2e/P2q force split; Doc126 P0 restore\"\n"
               << "}\n";
            std::cout << "[ophelie][doc105-p2] summary_json=" << local_cli.summary_json
                      << " provisional_inner_pass=" << (provisional_inner_pass ? 1 : 0)
                      << " mech_int=" << (interior_mom_mech ? 1 : 0)
                      << " mech_side_r=" << (side_radial_mech ? 1 : 0)
                      << " a_net_deep=" << s_ac.deep_interior.a_net_z_rms
                      << " a_net_nfs=" << s_ac.near_fs_interior.a_net_z_rms
                      << " a_surf_inc=" << s_ac.surface_incomplete.a_net_z_rms
                      << " a_surf_cmp=" << s_ac.surface_complete.a_net_z_rms
                      << " p2e_ff/wall_int=" << ff_dom_int << " p2e_ff/wall_surf=" << ff_dom_surf
                      << " p2e_wall/ff_side=" << wall_dom_side << std::endl;
        }
        return provisional_inner_pass ? 0 : 2;
    }

    // Doc105 P2c quiet keep with pollution budget (zero-buoyancy / freeze recommended).
    // Mid-keep CSV tracks velocity/drho only (single acoustic path — no diagnostic double-step).
    // End-of-keep: one clean acoustic_1st snap for zone a_net (same convention as stage probe).
    if (local_cli.doc105_p2_quiet_keep > TinyReal)
    {
        const Real zone_band = Real(2) * french.dp;
        const Real deep_band = Real(4) * french.dp;
        const Vecd glass_c = french.glass_center;
        const Real a_ref = std::abs(local_cli.gravity_g) + TinyReal;
        const Real budget_a_net_mech = Real(0.25) * a_ref;
        const Real budget_a_net_quiet = Real(0.05) * a_ref;
        const Real budget_u_rms_mech = Real(5.0e-4);
        const Real budget_u_rms_signal = Real(1.0e-4);
        const Real t_keep = local_cli.doc105_p2_quiet_keep;
        std::ofstream zcsv("output/french_doc105_p2c_quiet_keep_zones.csv");
        zcsv << "t_s,adv_step,Uz_rms_legacy_col,U_max,uz_int,uz_deep,uz_nfs,uz_surf,uz_surf_inc,uz_side,ur_side,"
                "uz_side_inc,uz_corner,uz_ws,uz_bot,pack_side,n_inc_side,pack_surf,n_inc_surf,drho_all\n";
        Real u_rms_end = 0, u_max_end = 0;
        Real uz_rms_end = 0, ur_rms_end = 0, utheta_rms_end = 0;
        FrenchHydroStageSnapshot mid{};
        size_t incomplete_vel_zero_n_last = 0;
        size_t side_band_vel_zero_n_last = 0;
        size_t bottom_band_vel_zero_n_last = 0;
        while (time_stepper.getPhysicalTime() + TinyReal < t_keep)
        {
            const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
            run_acoustic(acoustic_dt);
            if (local_cli.incomplete_vel_zero)
                incomplete_vel_zero_n_last = hostZeroVelocityIncompleteSupport(glass_particles, rho0);
            if (local_cli.side_band_vel_zero)
                side_band_vel_zero_n_last = hostZeroVelocitySideBand(
                    glass_particles, glass_c, french.glass_radius, z_bottom, z_free, zone_band);
            if (local_cli.bottom_band_vel_zero)
                bottom_band_vel_zero_n_last = hostZeroVelocityBottomBand(
                    glass_particles, glass_c, french.glass_radius, z_bottom, z_free, zone_band);
            if (advection_step(advection_time_step))
            {
                ++advection_steps;
                fluid_update_particle_position.exec();
                if (advection_steps % 50 == 0)
                    fluid_particle_sorting.exec();
                update_glass_cell_linked_list.exec();
                update_glass_relations.exec();
                density_regularization.exec();
                apply_quiet_hydro_after_dens();
                fluid_advection_step_setup.exec();
                if (local_cli.kernel_correction)
                    fluid_linear_correction_matrix.exec();
                fluid_viscous_force.exec();
                constant_gravity.exec();
                // No boussinesq in quiet-keep probe path.
                if (advection_steps % 20 == 0)
                {
                    // Velocity/drho only — Force after acoustic_2nd is not a_net-comparable.
                    mid = hostFrenchHydroStageSnapshot(glass_particles, rho0, local_cli.c0, local_cli.gravity_g,
                                                       glass_c, french.glass_radius, z_bottom, z_free, zone_band,
                                                       deep_band);
                    u_max_end = u_max_reduce.exec();
                    // Doc126: U_rms is 3D magnitude; legacy CSV column was Uz_rms.
                    u_rms_end = mid.all.u_rms;
                    zcsv << time_stepper.getPhysicalTime() << "," << advection_steps << "," << mid.all.uz_rms << ","
                         << u_max_end << "," << mid.interior.uz_rms << "," << mid.deep_interior.uz_rms << ","
                         << mid.near_fs_interior.uz_rms << "," << mid.surface.uz_rms << ","
                         << mid.surface_incomplete.uz_rms << "," << mid.wall_side.uz_rms << ","
                         << mid.wall_side.ur_rms << "," << mid.wall_side_incomplete.uz_rms << ","
                         << mid.wall_corner.uz_rms << "," << mid.wall_surface.uz_rms << ","
                         << mid.wall_bottom.uz_rms << "," << mid.wall_side.packing_mean << ","
                         << mid.wall_side.n_incomplete << "," << mid.surface.packing_mean << ","
                         << mid.surface.n_incomplete << "," << mid.all.drho_dt_rms << "\n";
                    std::cout << "[ophelie][doc105-p2c-keep] t=" << time_stepper.getPhysicalTime()
                              << " U_max=" << u_max_end << " U_rms=" << mid.all.u_rms
                              << " Uz_rms=" << mid.all.uz_rms << " Ur_rms=" << mid.all.ur_rms
                              << " uz_int=" << mid.interior.uz_rms
                              << " uz_deep=" << mid.deep_interior.uz_rms
                              << " uz_side=" << mid.wall_side.uz_rms
                              << " ur_side=" << mid.wall_side.ur_rms
                              << " uz_corner=" << mid.wall_corner.uz_rms
                              << " drho_all=" << mid.all.drho_dt_rms << std::endl
                              << std::flush;
                }
            }
        }

        // Doc131 A1: formal U/K from true production_end (no extra acoustic). Force residual is a
        // separate diagnostic_force_evaluation that snapshot/restores and does not advance PhysicalTime.
        FrenchHydroStageSnapshot production_end = hostFrenchHydroStageSnapshot(
            glass_particles, rho0, local_cli.c0, local_cli.gravity_g, glass_c, french.glass_radius, z_bottom, z_free,
            zone_band, deep_band);
        u_max_end = u_max_reduce.exec();
        u_rms_end = production_end.all.u_rms;
        uz_rms_end = production_end.all.uz_rms;
        ur_rms_end = production_end.all.ur_rms;
        utheta_rms_end = production_end.all.utheta_rms;
        Real K_end = 0, mass_sum_end = 0;
        Vecd U_mean_end = Vecd::Zero();
        {
            syncVariableToHost<Vecd>(glass_particles, "Velocity");
            syncVariableToHost<Real>(glass_particles, "Mass");
            const Vecd *vel_e = glass_particles.getVariableDataByName<Vecd>("Velocity");
            const Real *mass_e = glass_particles.getVariableDataByName<Real>("Mass");
            const size_t n_e = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n_e; ++i)
            {
                mass_sum_end += mass_e[i];
                U_mean_end += mass_e[i] * vel_e[i];
                K_end += Real(0.5) * mass_e[i] * vel_e[i].squaredNorm();
            }
            if (mass_sum_end > TinyReal)
                U_mean_end /= mass_sum_end;
        }
        const Real t_production_end = time_stepper.getPhysicalTime();
        FrenchHydroFieldSnapshot keep_prod_snap;
        hostCaptureFrenchHydroFields(glass_particles, keep_prod_snap);

        // Diagnostic force budget (acoustic_1st Interact+Update) from production snapshot; then restore.
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        density_regularization.exec();
        apply_quiet_hydro_after_dens();
        if (local_cli.hydro_rho_freeze_rho0)
            hostFreezeDensityToRho0(glass_particles, rho0);
        fluid_advection_step_setup.exec();
        if (local_cli.kernel_correction)
            fluid_linear_correction_matrix.exec();
        hostZeroVecFieldFrench(glass_particles, "ForcePrior");
        hostZeroVecFieldFrench(glass_particles, "GravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousGravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "BoussinesqBuoyancyForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
        constant_gravity.exec();
        hostZeroVecFieldFrench(glass_particles, "Force");
        const Real dt_ac_end = acoustic_time_step.exec();
        if (local_cli.kernel_correction)
            acoustic_1st_corr.exec(dt_ac_end);
        else
            acoustic_step_1st_half.exec(dt_ac_end);
        FrenchHydroStageSnapshot last = hostFrenchHydroStageSnapshot(
            glass_particles, rho0, local_cli.c0, local_cli.gravity_g, glass_c, french.glass_radius, z_bottom, z_free,
            zone_band, deep_band);
        hostRestoreFrenchHydroFields(glass_particles, keep_prod_snap);
        FrenchHydroStageSnapshot after_restore = hostFrenchHydroStageSnapshot(
            glass_particles, rho0, local_cli.c0, local_cli.gravity_g, glass_c, french.glass_radius, z_bottom, z_free,
            zone_band, deep_band);
        const bool probe_restore_ok =
            std::abs(after_restore.all.u_rms - production_end.all.u_rms) <=
                Real(1.0e-6) * (production_end.all.u_rms + Real(1)) &&
            std::abs(after_restore.all.rho_mean - production_end.all.rho_mean) <=
                Real(1.0e-5) * (std::abs(production_end.all.rho_mean) + Real(1));
        const Real t_after_diag = time_stepper.getPhysicalTime();
        const bool physical_time_unchanged =
            std::abs(t_after_diag - t_production_end) <= TinyReal * (std::abs(t_production_end) + Real(1));

        const bool mech_u = u_rms_end < budget_u_rms_mech;
        const bool signal_u = u_rms_end < budget_u_rms_signal;
        // Zone a_net thresholds use diagnostic_force_evaluation state (not production velocity).
        const bool interior_mom_quiet = last.interior.a_net_z_rms < budget_a_net_quiet;
        const bool interior_mom_mech = last.interior.a_net_z_rms < budget_a_net_mech;
        const bool side_radial_mech = last.wall_side.a_net_r_rms < budget_a_net_mech;
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            const bool deep_mom_quiet = last.deep_interior.a_net_z_rms < budget_a_net_quiet;
            const bool deep_mom_mech = last.deep_interior.a_net_z_rms < budget_a_net_mech;
            js << "{\n"
               << "  \"case\": \"french_doc105_p2c_quiet_keep\",\n"
               << "  \"doc105_p2c\": 1,\n"
               << "  \"doc105_p2i\": 1,\n"
               << "  \"doc105_p2j\": 1,\n"
               << "  \"doc105_p2l\": 1,\n"
               << "  \"doc105_p2o\": 1,\n"
               << "  \"doc105_p2q\": 1,\n"
               << "  \"doc105_p2r\": 1,\n"
               << "  \"keep_time_s\": " << t_keep << ",\n"
               << "  \"dens_sum_normalize\": " << (local_cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"settle_dens_sum_normalize\": " << (local_cli.settle_dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"dens_fs_cap_overpack\": " << (local_cli.dens_fs_cap_overpack ? 1 : 0) << ",\n"
               << "  \"dens_fs_cap_n\": " << dens_fs_cap_n_last << ",\n"
               << "  \"incomplete_vel_zero\": " << (local_cli.incomplete_vel_zero ? 1 : 0) << ",\n"
               << "  \"incomplete_vel_zero_n\": " << incomplete_vel_zero_n_last << ",\n"
               << "  \"side_band_vel_zero\": " << (local_cli.side_band_vel_zero ? 1 : 0) << ",\n"
               << "  \"side_band_vel_zero_n\": " << side_band_vel_zero_n_last << ",\n"
               << "  \"bottom_band_vel_zero\": " << (local_cli.bottom_band_vel_zero ? 1 : 0) << ",\n"
               << "  \"bottom_band_vel_zero_n\": " << bottom_band_vel_zero_n_last << ",\n"
               << "  \"incomplete_force_cancel_live\": "
               << (local_cli.incomplete_force_cancel_live ? 1 : 0) << ",\n"
               << "  \"incomplete_force_cancel_n\": " << incomplete_force_cancel_n_last << ",\n"
               << "  \"fs_incomplete_pack_scale_force\": "
               << (local_cli.fs_incomplete_pack_scale_force ? 1 : 0) << ",\n"
               << "  \"fs_incomplete_pack_scale_n\": " << fs_pack_scale_n_last << ",\n"
               << "  \"hydro_rho_freeze_rho0\": " << (local_cli.hydro_rho_freeze_rho0 ? 1 : 0) << ",\n"
               << "  \"no_extra_normalize\": "
               << ((local_cli.dens_sum_normalize || local_cli.settle_dens_sum_normalize) ? 0 : 1) << ",\n"
               << "  \"probe_extra_normalize\": " << (local_cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"dens_sum_interior_only\": " << (local_cli.dens_sum_interior_only ? 1 : 0) << ",\n"
               << "  \"dens_fs_atm_consistent\": " << (local_cli.dens_fs_atm_consistent ? 1 : 0) << ",\n"
               << "  \"background_hydro_p\": " << (local_cli.background_hydro_p ? 1 : 0) << ",\n"
               << "  \"kernel_correction\": " << (local_cli.kernel_correction ? 1 : 0) << ",\n"
               << "  \"acoustic_dens_mode\": " << local_cli.acoustic_dens_mode << ",\n"
               << "  \"reference_body_acc_z\": "
               << (local_cli.acoustic_dens_mode == 1 ? -local_cli.gravity_g : Real(0)) << ",\n"
               << "  \"budget_a_net_mech\": " << budget_a_net_mech << ",\n"
               << "  \"budget_a_net_quiet\": " << budget_a_net_quiet << ",\n"
               << "  \"budget_u_rms_mech\": " << budget_u_rms_mech << ",\n"
               << "  \"budget_u_rms_signal\": " << budget_u_rms_signal << ",\n"
               << "  \"doc131_a1_production_end\": 1,\n"
               << "  \"measurement_phase_U_K\": \"production_end\",\n"
               << "  \"measurement_phase_a_net\": \"diagnostic_force_evaluation\",\n"
               << "  \"physical_time_production_end\": " << t_production_end << ",\n"
               << "  \"physical_time_unchanged_by_diag\": " << (physical_time_unchanged ? 1 : 0) << ",\n"
               << "  \"probe_restore_ok\": " << (probe_restore_ok ? 1 : 0) << ",\n"
               << "  \"U_rms_end\": " << u_rms_end << ",\n"
               << "  \"Uz_rms_end\": " << uz_rms_end << ",\n"
               << "  \"Ur_rms_end\": " << ur_rms_end << ",\n"
               << "  \"Utheta_rms_end\": " << utheta_rms_end << ",\n"
               << "  \"U_rms_end_legacy_alias\": \"was_Uz_rms_before_Doc126\",\n"
               << "  \"U_rms_weighting\": \"particle_count\",\n"
               << "  \"K_end\": " << K_end << ",\n"
               << "  \"U_mean_x\": " << U_mean_end[0] << ",\n"
               << "  \"U_mean_y\": " << U_mean_end[1] << ",\n"
               << "  \"U_mean_z\": " << U_mean_end[2] << ",\n"
               << "  \"U_mean_weighting\": \"mass\",\n"
               << "  \"U_max_end\": " << u_max_end << ",\n"
               << "  \"eos_identity_ok_production\": "
               << ((production_end.all.p_eos_rms + TinyReal >= std::abs(production_end.all.p_eos_mean)) ? 1 : 0)
               << ",\n"
               << "  \"eos_identity_ok_all\": "
               << ((production_end.all.p_eos_rms + TinyReal >= std::abs(production_end.all.p_eos_mean)) ? 1 : 0)
               << ",\n"
               << "  \"doc126_p0\": 1,\n"
               << "  \"doc131_a1\": 1,\n"
               << "  \"mech_u_pass\": " << (mech_u ? 1 : 0) << ",\n"
               << "  \"signal_u_pass\": " << (signal_u ? 1 : 0) << ",\n"
               << "  \"interior_mom_mech_pass\": " << (interior_mom_mech ? 1 : 0) << ",\n"
               << "  \"interior_mom_quiet_pass\": " << (interior_mom_quiet ? 1 : 0) << ",\n"
               << "  \"deep_interior_mom_quiet_pass\": " << (deep_mom_quiet ? 1 : 0) << ",\n"
               << "  \"deep_interior_mom_mech_pass\": " << (deep_mom_mech ? 1 : 0) << ",\n"
               << "  \"side_radial_mech_pass\": " << (side_radial_mech ? 1 : 0) << ",\n";
            hostWriteFrenchZoneJson(js, "production_end_all", production_end.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "diagnostic_force_all", last.all);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_interior", last.interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_deep_interior", last.deep_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_near_fs_interior", last.near_fs_interior);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_surface", last.surface);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_surface_incomplete", last.surface_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_surface_complete", last.surface_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_side", last.wall_side);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_side_incomplete", last.wall_side_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_side_complete", last.wall_side_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_bottom", last.wall_bottom);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_bottom_incomplete", last.wall_bottom_incomplete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_bottom_complete", last.wall_bottom_complete);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_corner", last.wall_corner);
            js << ",\n";
            hostWriteFrenchZoneJson(js, "end_wall_surface", last.wall_surface);
            js << ",\n"
               << "  \"zone_csv\": \"output/french_doc105_p2c_quiet_keep_zones.csv\",\n"
               << "  \"passed\": " << (mech_u ? 1 : 0) << ",\n"
               << "  \"note\": \"Mid-keep CSV=velocity/drho only; end a_net from one acoustic_1st snap; "
                  "signal_u stricter lit-β floor\"\n"
               << "}\n";
        }
        std::cout << "[ophelie][doc105-p2c-keep] done U_rms=" << u_rms_end << " mech_u=" << (mech_u ? 1 : 0)
                  << " signal_u=" << (signal_u ? 1 : 0) << " a_net_int=" << last.interior.a_net_z_rms
                  << " a_net_deep=" << last.deep_interior.a_net_z_rms
                  << " a_net_r_side=" << last.wall_side.a_net_r_rms << std::endl;
        return mech_u ? 0 : 2;
    }

    if (local_cli.literature.enable_boussinesq)
    {
        boussinesq_force.exec();
    }
    auto apply_prescribed_temperature = [&]() {
        if (local_cli.prescribe_t_radial_dT > TinyReal)
        {
            hostPrescribeRadialHotCoreTemperature(glass_particles, french.glass_center, french.glass_radius,
                                                  french.glass_half_height, t0, local_cli.prescribe_t_radial_dT,
                                                  local_cli.prescribe_radial_z_taper, kOphelieTemperatureField,
                                                  kOphelieThermalDeltaTField);
        }
        else if (local_cli.prescribe_t_hot_bottom_dT > TinyReal)
        {
            hostPrescribeHotBottomLinearTemperature(glass_particles, french.glass_center, french.glass_half_height, t0,
                                                    local_cli.prescribe_t_hot_bottom_dT, kOphelieTemperatureField,
                                                    kOphelieThermalDeltaTField);
        }
    };
    const bool has_prescribed_t =
        local_cli.prescribe_t_radial_dT > TinyReal || local_cli.prescribe_t_hot_bottom_dT > TinyReal;
    if (has_prescribed_t)
    {
        apply_prescribed_temperature();
        if (local_cli.literature.enable_boussinesq)
        {
            boussinesq_force.exec();
            const OphelieBoussinesqForceAudit audit = hostAuditBoussinesqForcePrior(
                glass_particles, beta, t0, -local_cli.gravity_g, kOphelieTemperatureField);
            const char *tag = local_cli.prescribe_t_radial_dT > TinyReal ? "prescribed_radial" : "prescribed_hot_bottom";
            logOphelieBoussinesqForceAudit(audit, tag);
            // Radial drive: net Fz may be small; require force field matches analytic (rel_err).
            const bool force_ok = audit.rel_err < Real(1.0e-4);
            std::cout << "[ophelie][p1-03] prescribe_mode=" << tag
                      << " dT=" << (local_cli.prescribe_t_radial_dT > TinyReal ? local_cli.prescribe_t_radial_dT
                                                                               : local_cli.prescribe_t_hot_bottom_dT)
                      << " freeze=" << (local_cli.freeze_prescribed_temperature ? 1 : 0)
                      << " force_sign_up=" << (audit.f_z_sum > 0 ? 1 : 0) << " force_ok=" << (force_ok ? 1 : 0)
                      << " Fz_sum=" << audit.f_z_sum << std::endl;
        }
    }
    if (write_vtp)
    {
        glass_state_recorder.writeToFile(0);
    }

    std::cout << "[ophelie][stage4.1] flow start: end_time=" << local_cli.end_time
              << " n_glass=" << glass_particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " beta=" << beta << " mu=" << mu
              << " k=" << k_th << " c0=" << local_cli.c0 << " T0=" << t0
              << " boussinesq=" << (local_cli.literature.enable_boussinesq ? 1 : 0)
              << " meridional_avg_start=" << local_cli.meridional_average_start
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0) << std::endl;

    FrenchMeridionalFlowAccumulator meridional;
    meridional.reset(cyl_frame);
    bool meridional_ts_header = true;
    const std::string meridional_ts_path = "output/french_meridional_timeseries.csv";
    const Real r_axis_hi_ts = Real(0.08);
    const Real r_rim_lo_ts = Real(0.18);

    while (!time_stepper.isEndTime(local_cli.end_time))
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        run_acoustic(acoustic_dt);

        if (!local_cli.literature.zero_q)
        {
            if (local_cli.literature.constant_thermal_rheology)
            {
                apply_q_const.exec(acoustic_dt);
            }
            else
            {
                apply_q_variable.exec(acoustic_dt);
            }
            coupling_state.ledger.accumulateQ(coupling_state.last_power.p_heat_fluid, acoustic_dt);
        }
        if (!local_cli.literature.constant_thermal_rheology)
        {
            const Real mu_now = hostUpdateFrenchThermalRheologyFromTemperature(
                glass_particles, laws, false, t0, rheology_stats);
            glass_viscosity.setReferenceViscosity(mu_now);
            invalidateViscosityDependentCkKernels(glass_inner, glass_wall_contact);
        }
        if (!local_cli.freeze_prescribed_temperature)
        {
            if (use_skull_shell_thermal)
            {
                glass_shell_thermal->exec(acoustic_dt);
                wall_shell_thermal->exec(acoustic_dt);
            }
            else if (laplace_temperature && apply_diffusion)
            {
                laplace_temperature->exec();
                apply_diffusion->exec(acoustic_dt);
            }
            {
                const Real p_boundary_now = boundary_loss_power.exec();
                coupling_state.ledger.p_loss_live = p_boundary_now;
                coupling_state.ledger.accumulateBoundary(p_boundary_now, acoustic_dt);
            }
            apply_natural_bc.exec(acoustic_dt);
        }
        if (local_cli.literature.enable_boussinesq)
        {
            boussinesq_force.exec();
        }

        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            // Keep geometric freeze-T after motion (Doc77-C style); else residual advection scrambles drive.
            if (local_cli.freeze_prescribed_temperature && has_prescribed_t)
            {
                apply_prescribed_temperature();
            }
            if (advection_steps % 50 == 0)
            {
                fluid_particle_sorting.exec();
            }
            update_glass_cell_linked_list.exec();
            update_wall_cell_linked_list.exec();
            update_glass_relations.exec();
            update_wall_relations.exec();
            if (use_skull_shell_thermal)
            {
                (void)hostApplyFrenchSkullShellOuterCoolDirichlet(wall_particles, french, wall_thickness,
                                                                  local_cli.skull_shell.t_cool, t0);
            }
            density_regularization.exec();
            apply_quiet_hydro_after_dens();
            maybe_dens_sum_audit(time_stepper.getPhysicalTime(), advection_steps);
            fluid_advection_step_setup.exec();
            if (local_cli.kernel_correction)
                fluid_linear_correction_matrix.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
            if (local_cli.literature.enable_boussinesq)
            {
                boussinesq_force.exec();
            }
            if (periodic_coupling)
            {
                hostSampleEulerQOntoParticles(glass_particles, q_grid);
                recordFrenchFluidHeatPower(coupling_state, glass_particles);
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
                    hostSampleEulerQOntoParticles(glass_particles, q_grid);
                    recordFrenchFluidHeatPower(coupling_state, glass_particles);
                    const bool refresh_solver_ok =
                        std::isfinite(phi_eq) && phi_eq < Real(1.0e-2) && std::isfinite(p_joule) && p_joule > TinyReal;
                    noteFrenchPeriodicRefresh(coupling_state, t_now, local_cli.literature.em_update_interval_s,
                                              refresh_solver_ok);
                    ++coupling_state.update_count;
                    ++coupling_state.periodic_refresh_count;
                    writeAxisymmetricCsv("output/french_sigma_rz.csv", coupling_state.sigma_rz, "sigma_Sm");
                    writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
                    appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, t_now,
                                              local_cli.literature.em_control);
                    std::cout << "[ophelie][thermo-em] update=" << coupling_state.update_count
                              << " periodic_refresh=" << coupling_state.periodic_refresh_count << " t=" << t_now
                              << " next_due=" << coupling_state.next_due_time
                              << " max_field_age=" << coupling_state.max_field_age_s
                              << " reconstructed_glass_power=" << p_joule
                              << " P_heat_fluid=" << coupling_state.last_power.p_heat_fluid
                              << " coil_current=" << coupling_state.coil_current_per_loop
                              << " phi_eq_res_vol=" << phi_eq << std::endl;
                }
            }
            if (write_vtp && state_recording())
            {
                glass_state_recorder.writeToFile(advection_steps);
            }
            if (advection_steps % 20 == 0)
            {
                const Real t_now = time_stepper.getPhysicalTime();
                const Real energy = local_cli.literature.constant_thermal_rheology
                                        ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                        : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
                const Real dt_e = t_now - prev_energy_time;
                const Real dE_dt = dt_e > TinyReal ? (energy - prev_energy) / dt_e : Real(0);
                Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
                hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc,
                                                              thermal_bc.boundary_shell_thickness, side_w, bottom_w,
                                                              free_conv_w, free_rad_w, total_loss_w);
                Real contact_melt_loss_w = 0, contact_wall_gain_w = 0;
                if (use_skull_shell_thermal)
                {
                    hostFrenchSkullShellContactHeatPowers(glass_particles, wall_particles, glass.Name(), wall.Name(),
                                                          contact_melt_loss_w, contact_wall_gain_w);
                    // Contact replaces side/bottom Robin in the energy residual.
                    side_w = contact_melt_loss_w;
                    bottom_w = Real(0);
                    total_loss_w = contact_melt_loss_w + free_conv_w + free_rad_w;
                }
                const Real u_now = u_max_reduce.exec();
                const Real t_mean_now = t_mean_reduce.exec();
                const Real t_max_now = t_max_reduce.exec();
                const Real t_min_now = t_min_reduce.exec();
                const Real dT_now = t_max_now - t_min_now;
                const Real u_buoy = frenchBuoyancySpeedScale(local_cli.gravity_g, beta, dT_now, glass_h);
                const Real residual =
                    frenchEnergyResidual(dE_dt, coupling_state.last_power.p_heat_fluid, side_w,
                                         free_conv_w + free_rad_w, bottom_w);
                const Real residual_rel =
                    residual / (std::abs(coupling_state.last_power.p_heat_fluid) + TinyReal);
                if (!local_cli.freeze_prescribed_temperature)
                {
                    coupling_state.ledger.noteLiveLoss(side_w, bottom_w, free_conv_w + free_rad_w);
                }
                const FrenchMeltSpatialStats spatial =
                    hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
                appendFrenchSpatialStatsCsv("output/french_spatial_stats.csv", t_now, spatial);
                writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame,
                                           kOphelieTemperatureField);
                if (t_now + TinyReal >= local_cli.meridional_average_start)
                {
                    hostAccumulateFrenchMeridionalFlow(meridional, glass_particles, t_now);
                }
                // Doc97 P3: instantaneous band scalars at monitor cadence (true-Δt window weights later).
                {
                    const FrenchMeridionalInstantBands bands =
                        hostFrenchMeridionalInstantBands(glass_particles, cyl_frame, r_axis_hi_ts, r_rim_lo_ts);
                    appendFrenchMeridionalTimeseriesCsv(meridional_ts_path, t_now, bands, meridional_ts_header);
                    meridional_ts_header = false;
                }
                if (periodic_coupling)
                {
                    FrenchEnergyBudgetRowExtras extras;
                    extras.p_em_recon = p_joule;
                    extras.p_heat_fluid = coupling_state.last_power.p_heat_fluid;
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
                    }
                    appendFrenchEnergyBudgetCsv("output/french_energy_budget.csv", t_now, energy, dE_dt,
                                                coupling_state.last_power.p_heat_fluid, side_w,
                                                free_conv_w + free_rad_w, bottom_w, t_min_now, t_mean_now,
                                                t_max_now, u_now, u_buoy,
                                                thermal_bc.enable_diffusion ? 1 : 0, &rheology_stats,
                                                local_cli.literature.constant_thermal_rheology ? 1 : 0, &extras);
                }
                prev_energy = energy;
                prev_energy_time = t_now;
                const GlassEscapeCensus mid_escape = hostCountEscapedGlassParticles(
                    glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
                const Real mid_frac =
                    mid_escape.total > 0 ? static_cast<Real>(mid_escape.escaped) / static_cast<Real>(mid_escape.total)
                                         : Real(0);
                monitor << advection_steps << "," << t_now << "," << u_now << "," << t_min_now << "," << t_mean_now
                        << "," << t_max_now << "," << dT_now << "," << u_buoy << "," << mid_escape.escaped << ","
                        << mid_frac << "," << p_joule << "," << dE_dt << "," << residual << "," << residual_rel << ","
                        << total_loss_w << "," << coupling_state.update_count << "," << spatial.u_rms << ","
                        << spatial.u_theta_rms << "," << spatial.u_z_rms << "," << spatial.t_std << ","
                        << spatial.t_outer << "," << spatial.t_center << "," << spatial.t_bottom << ","
                        << spatial.t_top << "\n";
                monitor.flush();
                std::cout << "[ophelie][stage4.1] N=" << advection_steps << " t=" << t_now << " U_max=" << u_now
                          << " U_rms=" << spatial.u_rms << " U_th=" << spatial.u_theta_rms << " U_z=" << spatial.u_z_rms
                          << " U_buoy=" << u_buoy << " T_min=" << t_min_now << " T_mean=" << t_mean_now
                          << " T_max=" << t_max_now << " dT=" << dT_now << " T_std=" << spatial.t_std
                          << " T_out=" << spatial.t_outer << " T_ctr=" << spatial.t_center
                          << " escaped=" << mid_escape.escaped << "/" << mid_escape.total
                          << " residual=" << residual << " residual_rel=" << residual_rel
                          << " em_updates=" << coupling_state.update_count;
                if (use_skull_shell_thermal)
                {
                    std::cout << " contact_melt_loss_W=" << contact_melt_loss_w
                              << " contact_wall_gain_W=" << contact_wall_gain_w;
                }
                std::cout << std::endl;
            }
        }
    }

    const Real u_max = u_max_reduce.exec();
    const Real t_mean = t_mean_reduce.exec();
    const Real t_max = t_max_reduce.exec();
    const Real t_min = t_min_reduce.exec();
    const Real dT = t_max - t_min;
    const Real u_buoy = frenchBuoyancySpeedScale(local_cli.gravity_g, beta, dT, glass_h);
    const FrenchMeltSpatialStats spatial_end =
        hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
    writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame, kOphelieTemperatureField);
    const GlassEscapeCensus escape = hostCountEscapedGlassParticles(glass_particles, gc, french.glass_radius, z_bottom,
                                                                   z_rim, french.dp);
    const Real escaped_frac =
        escape.total > 0 ? static_cast<Real>(escape.escaped) / static_cast<Real>(escape.total) : Real(0);
    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                  side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
    Real contact_melt_loss_w = 0, contact_wall_gain_w = 0;
    if (use_skull_shell_thermal)
    {
        hostFrenchSkullShellContactHeatPowers(glass_particles, wall_particles, glass.Name(), wall.Name(),
                                              contact_melt_loss_w, contact_wall_gain_w);
        side_w = contact_melt_loss_w;
        bottom_w = Real(0);
        total_loss_w = contact_melt_loss_w + free_conv_w + free_rad_w;
        std::cout << "[ophelie][skull-shell] contact audit: melt_loss_W=" << contact_melt_loss_w
                  << " wall_gain_W=" << contact_wall_gain_w
                  << " free_W=" << (free_conv_w + free_rad_w) << " total_W=" << total_loss_w << " vs P=" << p_joule
                  << "\n";
    }
    if (write_vtp)
    {
        glass_state_recorder.writeToFile(advection_steps + 1);
    }
    monitor.close();

    // Ensure at least one terminal sample for Doc97 P3 window post-process / short runs.
    {
        const Real t_end_ts = time_stepper.getPhysicalTime();
        const FrenchMeridionalInstantBands bands =
            hostFrenchMeridionalInstantBands(glass_particles, cyl_frame, r_axis_hi_ts, r_rim_lo_ts);
        appendFrenchMeridionalTimeseriesCsv(meridional_ts_path, t_end_ts, bands, meridional_ts_header);
        meridional_ts_header = false;
    }
    const char *merid_tag = local_cli.literature.enable_boussinesq ? "boussinesq_on" : "zero_buoyancy";
    meridional.writeCsv("output/french_meridional_ur_uz_mean.csv");
    meridional.logSummary(merid_tag);
    const Real R_melt = french.glass_radius;
    const Real uz_axis = meridional.meanUzRadiusBand(Real(0), Real(0.08));
    const Real uz_rim = meridional.meanUzRadiusBand(Real(0.18), R_melt + TinyReal);
    const Real u_merid = meridional.meanMeridionalSpeed();
    const Real circ_proxy = meridional.meanCircProxy();
    const Real u_rms_avg = meridional.meanInstUrms();
    {
        // Refresh Ra estimate with measured end-of-run ΔT for the detectability note.
        const Real mu_end = local_cli.literature.constant_thermal_rheology ? mu : rheology_stats.mu_mean;
        const Real mu_hot_end = laws.mu.evaluate(t_mean + Real(0.5) * dT);
        logFrenchNaturalRaPrEstimate(french.glass_radius, glass_h, dT, rho0, mu_end, cp, k_th, beta,
                                     local_cli.gravity_g, mu_hot_end, u_rms_avg > TinyReal ? u_rms_avg : Real(0.02));
        std::cout << "[ophelie][natural-detect] boussinesq=" << (local_cli.literature.enable_boussinesq ? 1 : 0)
                  << " U_merid_mean=" << u_merid << " circ_proxy=" << circ_proxy << " U_rms_avg=" << u_rms_avg
                  << " Uz_axis_rlt0.08=" << uz_axis << " Uz_rim_rgt0.18=" << uz_rim
                  << " dens_sum_normalize=" << (local_cli.dens_sum_normalize ? 1 : 0)
                  << " background_hydro_p=" << (local_cli.background_hydro_p ? 1 : 0)
                  << " kernel_correction=" << (local_cli.kernel_correction ? 1 : 0)
                  << " compare_zero_buoyancy_run_with_--zero-buoyancy_same_end-time\n";
    }

    const Real power_rel_err =
        std::abs(p_joule - local_cli.target_power_w) / (local_cli.target_power_w + TinyReal);
    const Real physical_time = time_stepper.getPhysicalTime();
    const Real energy_end = local_cli.literature.constant_thermal_rheology
                                ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
    const bool periodic_coupling_end = periodic_coupling;
    FrenchCouplingGateInputs gate_in;
    gate_in.periodic_coupling = periodic_coupling_end;
    gate_in.zero_q = local_cli.literature.zero_q;
    gate_in.constant_thermal_rheology = local_cli.literature.constant_thermal_rheology;
    gate_in.require_u_max_positive = false;
    gate_in.energy_residual_available = coupling_state.ledger.h0_set && physical_time > TinyReal;
    gate_in.budget_exhausted = false;
    gate_in.diverged = false;
    gate_in.fluid_power_recorded = coupling_state.last_power.fluid_power_recorded;
    gate_in.ledger_invalid_accumulate = coupling_state.ledger.invalid_accumulate;
    // Doc 45: natural currently does not measure out-of-grid; do not treat hard-coded zeros as a pass.
    gate_in.out_of_grid_measured = false;
    gate_in.end_time = local_cli.end_time;
    gate_in.em_interval = local_cli.literature.em_update_interval_s;
    gate_in.physical_time = physical_time;
    gate_in.p_em = p_joule;
    gate_in.p_heat_fluid = coupling_state.last_power.p_heat_fluid;
    gate_in.phi_eq = phi_eq;
    gate_in.power_rel_err = power_rel_err;
    gate_in.q_sample_rel_l2 = Real(-1);
    gate_in.n_out_of_grid_init = Real(0);
    gate_in.n_out_of_grid_now = Real(0);
    gate_in.max_out_of_grid = Real(0);
    gate_in.u_max = u_max;
    gate_in.c0 = local_cli.c0;
    gate_in.t_mean = t_mean;
    gate_in.t_max = t_max;
    gate_in.escaped_frac = escaped_frac;
    gate_in.escaped_frac_max = local_cli.escaped_frac_max;
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

    std::cout << "test_3d_ophelie_french_natural_convection_frozen_q"
              << " P_joule_W=" << p_joule << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << phi_eq
              << " end_time=" << local_cli.end_time << " U_max=" << u_max << " U_rms=" << spatial_end.u_rms
              << " U_th=" << spatial_end.u_theta_rms << " U_z=" << spatial_end.u_z_rms << " U_buoy=" << u_buoy
              << " T_min=" << t_min << " T_mean=" << t_mean << " T_max=" << t_max << " dT=" << dT
              << " T_std=" << spatial_end.t_std << " T_out=" << spatial_end.t_outer << " T_ctr=" << spatial_end.t_center
              << " escaped=" << escape.escaped << "/" << escape.total << " escaped_frac=" << escaped_frac
              << " r_max=" << escape.r_max << " z=[" << escape.z_min << "," << escape.z_max << "]"
              << " wall_loss_side_W=" << side_w << " wall_loss_bottom_W=" << bottom_w
              << " free_conv_loss_W=" << free_conv_w << " free_rad_loss_W=" << free_rad_w
              << " total_heat_loss_W=" << total_loss_w;
    if (use_skull_shell_thermal)
    {
        std::cout << " contact_melt_loss_W=" << contact_melt_loss_w
                  << " contact_wall_gain_W=" << contact_wall_gain_w;
    }
    std::cout << " thermal_diffusion=" << (thermal_bc.enable_diffusion || use_skull_shell_thermal ? 1 : 0)
              << " P_heat_fluid=" << coupling_state.last_power.p_heat_fluid
              << " E_Q_applied=" << coupling_state.ledger.e_q_applied
              << " E_boundary_out=" << coupling_state.ledger.e_boundary_out
              << " energy_residual_H=" << coupling_state.ledger.residualH(energy_end)
              << " energy_residual_H_rel=" << gate_in.cumulative_energy_residual_rel
              << " freeze_prescribed_t=" << (local_cli.freeze_prescribed_temperature ? 1 : 0)
              << " startup_ok=" << (gates.startup_ok ? 1 : 0)
              << " completed_requested_time=" << (gates.completed_requested_time ? 1 : 0)
              << " periodic_refresh_ok=" << (gates.periodic_refresh_ok ? 1 : 0)
              << " em_solver_ok=" << (gates.em_solver_ok ? 1 : 0) << " q_transfer_ok=" << (gates.q_transfer_ok ? 1 : 0)
              << " flow_sanity_ok=" << (gates.flow_sanity_ok ? 1 : 0)
              << " thermal_energy_ok=" << (gates.thermal_energy_ok ? 1 : 0)
              << " smoke_passed=" << (gates.smoke_passed ? 1 : 0)
              << " validation_exercised=" << (gates.validation_exercised ? 1 : 0)
              << " validation_passed=" << (gates.validation_passed ? 1 : 0)
              << " require_validation=" << (local_cli.require_validation ? 1 : 0)
              << " fail_reasons=" << (gates.fail_reasons.empty() ? "none" : gates.fail_reasons)
              << " coupling_updates=" << coupling_state.update_count
              << " periodic_refreshes=" << coupling_state.periodic_refresh_count
              << " coupling=" << thermoEMCouplingModeName(local_cli.literature.coupling)
              << " em_control=" << emControlModeName(local_cli.literature.em_control)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    const std::vector<std::pair<Real, Real>> default_windows = {
        {Real(5), Real(10)}, {Real(10), Real(15)}, {Real(15), Real(20)}};
    const auto meridional_windows = frenchMeridionalWindowsFromTimeseriesCsv(meridional_ts_path, default_windows);
    if (!local_cli.summary_json.empty())
    {
        const Real requested_beta = local_cli.beta;
        const Real effective_beta = beta;
        std::ofstream js(local_cli.summary_json);
        js << "{\n"
           << "  \"case\": \"french_natural_frozen_q\",\n"
           << "  \"doc95\": 1,\n"
           << "  \"doc97_p3\": 1,\n"
           << "  \"dp\": " << french.dp << ",\n"
           << "  \"end_time\": " << local_cli.end_time << ",\n"
           << "  \"beta\": " << beta << ",\n"
           << "  \"requested_beta\": " << requested_beta << ",\n"
           << "  \"effective_beta\": " << effective_beta << ",\n"
           << "  \"boussinesq\": " << (local_cli.literature.enable_boussinesq ? 1 : 0) << ",\n"
           << "  \"freeze_prescribed_t\": " << (local_cli.freeze_prescribed_temperature ? 1 : 0) << ",\n"
           << "  \"prescribe_t_radial_dT\": " << local_cli.prescribe_t_radial_dT << ",\n"
           << "  \"c0\": " << local_cli.c0 << ",\n"
           << "  \"hydro_settle_time\": " << local_cli.hydro_settle_time << ",\n"
           << "  \"meridional_average_start\": " << local_cli.meridional_average_start << ",\n"
           << "  \"dens_sum_normalize\": " << (local_cli.dens_sum_normalize ? 1 : 0) << ",\n"
           << "  \"no_extra_normalize\": " << (local_cli.dens_sum_normalize ? 0 : 1) << ",\n"
           << "  \"dens_sum_interior_only\": " << (local_cli.dens_sum_interior_only ? 1 : 0) << ",\n"
           << "  \"dens_fs_atm_consistent\": " << (local_cli.dens_fs_atm_consistent ? 1 : 0) << ",\n"
           << "  \"dens_sum_surface_depth\": " << dens_audit_surface_depth << ",\n"
           << "  \"dens_sum_n_written_last\": " << dens_sum_n_written_last << ",\n"
           << "  \"dens_write_max_abs_drho\": " << dens_write_last.max_abs_drho << ",\n"
           << "  \"dens_write_rms_drho\": " << dens_write_last.rms_drho << ",\n"
           << "  \"dens_write_max_abs_dp\": " << dens_write_last.max_abs_dp << ",\n"
           << "  \"dens_write_rms_dp\": " << dens_write_last.rms_dp << ",\n"
           << "  \"fs_atm_n_written\": " << fs_atm_last.n_written << ",\n"
           << "  \"fs_atm_rms_p_total_after\": " << fs_atm_last.rms_p_total_after << ",\n"
           << "  \"dens_sum_surface_recapture_interval\": "
           << local_cli.dens_sum_surface_recapture_interval << ",\n"
           << "  \"dens_sum_n_recaptured_last\": " << dens_sum_n_recaptured_last << ",\n"
           << "  \"background_hydro_p\": " << (local_cli.background_hydro_p ? 1 : 0) << ",\n"
           << "  \"french_pressure_path\": \"total_p_plus_GravityForceCK\",\n"
           << "  \"perturbation_momentum\": 0,\n"
           << "  \"kernel_correction\": " << (local_cli.kernel_correction ? 1 : 0) << ",\n"
           << "  \"dt_scale\": " << local_cli.dt_scale << ",\n"
           << "  \"acoustic_cfl\": " << acoustic_cfl_eff << ",\n"
           << "  \"advection_cfl\": " << advection_cfl_eff << ",\n"
           << "  \"dt_acoustic_est\": " << dt_ac_est << ",\n"
           << "  \"viscous_tightens_advection\": " << (u_visc_ref > local_cli.u_ref + TinyReal ? 1 : 0)
           << ",\n"
           << "  \"layout_jitter_dp_frac\": " << local_cli.layout_jitter_dp_frac << ",\n"
           << "  \"layout_seed\": " << local_cli.layout_seed << ",\n"
           << "  \"sigma0_mean\": " << dens_sigma0_mean << ",\n"
           << "  \"dens_sum_audit\": " << (local_cli.dens_sum_audit ? 1 : 0) << ",\n"
           << "  \"dens_sum_n_fs_clamp\": " << dens_audit_last.n_fs_clamp << ",\n"
           << "  \"dens_sum_n_override\": " << dens_audit_last.n_override_active << ",\n"
           << "  \"dens_sum_n_migrated_to_surface\": " << dens_audit_last.n_migrated_to_surface << ",\n"
           << "  \"dens_sum_mean_abs_drho_override\": " << dens_audit_last.mean_abs_drho_override << ",\n"
           << "  \"dens_sum_max_abs_drho_override\": " << dens_audit_last.max_abs_drho_override << ",\n"
           << "  \"dens_sum_mass_sum\": " << dens_audit_last.mass_sum << ",\n"
           << "  \"dens_sum_mean_ratio_surf\": " << dens_audit_last.mean_ratio_surf << ",\n"
           << "  \"dens_sum_mean_ratio_int\": " << dens_audit_last.mean_ratio_int << ",\n"
           << "  \"U_merid_mean\": " << u_merid << ",\n"
           << "  \"circ_proxy\": " << circ_proxy << ",\n"
           << "  \"U_rms_avg\": " << u_rms_avg << ",\n"
           << "  \"Uz_axis\": " << uz_axis << ",\n"
           << "  \"Uz_rim\": " << uz_rim << ",\n"
           << "  \"U_max\": " << u_max << ",\n"
           << "  \"U_rms_end\": " << spatial_end.u_rms << ",\n"
           << "  \"dT\": " << dT << ",\n"
           << "  \"T_ctr\": " << spatial_end.t_center << ",\n"
           << "  \"T_out\": " << spatial_end.t_outer << ",\n"
           << "  \"escaped_frac\": " << escaped_frac << ",\n"
           << "  \"smoke_passed\": " << (gates.smoke_passed ? 1 : 0) << ",\n"
           << "  \"passed\": " << (passed ? 1 : 0) << ",\n"
           << "  \"hydrostatic_pass\": null,\n"
           << "  \"meridional_timeseries\": \"" << meridional_ts_path << "\",\n";
        writeFrenchMeridionalWindowsJsonArray(js, meridional_windows);
        js << ",\n"
           << "  \"note\": \"Doc97 P3: lit-β SNR + non-overlap windows; not u_max gate; signal_verified OPEN until repeats\"\n"
           << "}\n";
        std::cout << "[ophelie] summary_json=" << local_cli.summary_json << std::endl;
    }
    return passed ? 0 : 1;
}
