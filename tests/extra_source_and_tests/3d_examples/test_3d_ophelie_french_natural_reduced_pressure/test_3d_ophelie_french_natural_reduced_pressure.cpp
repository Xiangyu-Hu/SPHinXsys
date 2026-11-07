/**
 * @file test_3d_ophelie_french_natural_reduced_pressure.cpp
 * @brief Doc126/Doc131: French cylinder melt, closed flat-lid domain, reduced-pressure (p') WCSPH
 *        with production BoussinesqBuoyancyForceCK, plus the Doc131 thermal completion
 *        (power contract, boundary-area audit, conduction + Robin/radiation boundary, T0-T5 tests).
 *
 * Physics contract (flow, unchanged from Doc126):
 *  - Melt: cylinder R=0.25 m, H=0.185 m (center (0,0,H/2), bottom at z=0).
 *  - Reduced pressure: PerturbationMomentum=1, AcousticDensMode=1, BackgroundPressure p_h = rho0 g (H - z)
 *    rewritten each advection step (host, base only). Default ReferenceBodyAcceleration=0 (g not in ForcePrior).
 *    Optional --cref-buoyancy: a_ref := F_boussinesq/m each advection step so AcousticDensMode=1 wall dens
 *    jump C_ext-C_ref cancels the buoyancy part (Doc133/134 D1 mean-density drift fix).
 *  - NO GravityForceCK. BoussinesqBuoyancyForceCK is kept (Gravity object exists for it only).
 *  - Closed domain: continuity-only density by default (no FreeSurface regularization, no density
 *    summation cap, NO analytic density rewrite per step). Optional --density-summation uses the
 *    Internal DensityRegularization path (diagnostic; it relaxes p' towards 0).
 *  - Acoustic pressure contact: wall_closed (side + bottom + flat lid, impermeable).
 *  - Viscous contact: wall_noslip (open cup, rim=0; side + bottom no-slip). The lid has no viscous
 *    wall contact => top is flat free-slip (top_bc = flat_free_slip_viscous_omit).
 *
 * Doc131 thermal contract:
 *  - Power: fixed_discrete_absorbed_power (default; V_i^Q = m_i/rho0, q_i = P s_i / sum(s V^Q)) or the
 *    diagnostic fixed_volumetric_density (q = P / V_geom).
 *  - Conduction + boundary: the production CK kernels OpheliePairwiseLaplaceCK + ApplyLaplaceDiffusionDtCK
 *    + ApplyOphelieThermalFrenchNaturalBcCK (same kernels as the natural/stirring cases), on the fluid body.
 *    Face area weights come from the host audit helper (B3): shell thickness dp, A_i = m_i/(rho0 dp),
 *    exclusive primary-face assignment, optional normalisation to pi R^2 / 2 pi R H (--area-normalize).
 *    Side/bottom exchange with t_cool, top with t_ambient (+ radiation); heat flows both ways (no clamp).
 *  - Modes: R0 R1 R2 R3 (flow mechanism) | baseline (thermal smoke baseline) | T0..T5 (thermal tests).
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_export.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_implicit_viscosity_split.h"
#include "electromagnetic_ophelie_viscosity_operator_diag.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

extern char **environ;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
constexpr double kStefanBoltzmann = 5.670374419e-8;
constexpr double kPiD = 3.14159265358979323846;
constexpr const char *kDoc131QField = "OphelieThermalSourceQ";
constexpr const char *kPowerDiscrete = "fixed_discrete_absorbed_power";
constexpr const char *kPowerVolumetric = "fixed_volumetric_density";

//----------------------------------------------------------------------------------------------------
// CLI
//----------------------------------------------------------------------------------------------------
struct Cli
{
    std::string mode = "R0"; /**< R0 | R1 | R2 | R3 | baseline | T0..T5 */
    Real dp = 0.025;
    Real settle_time = 0.5;
    Real end_time = 1.0; /**< observation window after settle */
    bool end_time_set = false;
    Real rho0 = 2750.0;
    Real gravity_g = 9.81;
    Real beta = 3.0e-5;
    Real t_ref = 1473.0;
    Real mu = 4.0;
    Real c0 = 20.0;
    Real u_ref = 0.5;
    Real a_strat = 50.0;   /**< K/m (R1, T2) */
    Real dT_radial = 10.0; /**< K (R2, T2) */
    Real perturb_u = 1.0e-3; /**< m/s (R3) */
    Real glass_radius = 0.25;
    Real glass_height = 0.185;
    Real q_power = 50000.0; /**< W */
    Real cp = 1150.0;       /**< J/(kg K); CEP contract value (Doc131 A3). 1300 is a historical value only. */
    Real k = 4.0;           /**< W/(m K), registered reference (Doc131 A3) */
    Real h_side = 300.0;    /**< W/(m^2 K), CODE_PROVISIONAL */
    Real h_bottom = 35.0;   /**< W/(m^2 K), CODE_PROVISIONAL */
    Real h_top = 20.0;      /**< W/(m^2 K), CODE_PROVISIONAL */
    Real t_cool = 300.0;    /**< K, side/bottom coolant sink, CODE_PROVISIONAL */
    Real emissivity = 0.8;  /**< top radiation */
    Real t_ambient = 300.0; /**< K, top convection ambient */
    Real t_rad_ambient = 300.0; /**< K, top radiation ambient; Doc139 B uses 423 vs conv 303 */
    bool t_rad_ambient_set = false;
    bool variable_properties = false; /**< Doc139 B1: μ/cp/k(T) particle fields */
    bool local_particle_viscosity = true; /**< Doc143: pairwise μ from OphelieDynamicViscosity (not body-mean) */
    bool mean_mu_viscous = false;     /**< legacy: force ViscousForceCK to use volume-mean μ */
    /** Doc175: explicit=ViscousForceCK+dt_visc; implicit_split=InteractionSplit DampingPairwise (CPU). */
    std::string viscosity_mode = "explicit";
    bool viscosity_operator_diag = false; /**< Doc177: fixed-geometry a_exp vs Δv/dt_imp then exit */
    bool viscosity_diag_inner_only = false; /**< Doc177: skip wall damping/force in operator-diag */
    Real viscosity_diag_dt = 1.0e-4;      /**< Doc177 operator-diag probe dt */
    Real viscosity_diag_shear = 0.05;     /**< m/s amplitude; Ux = A*sin(2π(z-z_b)/H) */
    /** Doc177: acoustic=every dt_ac (default); advection=once per Dt like explicit ForcePrior hold. */
    std::string implicit_visc_schedule = "acoustic";
    Real mu_cap = 1.0e5;              /**< Pa·s clamp for variable μ (thesis calc upper bound) */
    Real marangoni_dgamma_dt = 0.0;   /**< Doc139 B2: dγ_surface/dT [N/(m K)]; 0 = off */
    bool marangoni_verify = false;    /**< Doc144: host axisymmetric direction+∫ check (auto-on if dγ/dT≠0) */
    bool marangoni_verify_set = false;
    std::string literature_id;        /**< e.g. cep_natural_A1 / thesis_static_B0 */
    std::string material_curve_dir;   /**< directory with thesis_III8/9/10_*_digitized.csv */
    Real wall_thickness_factor = 4.0;
    Real thermal_dt = 0.0;  /**< frozen-flow thermal tests only; 0 = mode default */
    Real t4_amp = 10.0;     /**< K, T4 cosine amplitude */
    int mesh_resolution = 20;
    bool kernel_correction = false;
    bool use_density_summation = false; /**< continuity-only by default */
    bool thaw_T_set = false;
    bool thaw_T = false;
    bool top_loss = true;
    bool r2_keep_mean = false;
    bool r2_reverse = false;          /**< Doc131 D2: negate the radial dT (cold wall) */
    bool r2_fixed_spatial_T = false;  /**< Doc131 D2: re-apply T(r) from the CURRENT position every advection step (diagnostic) */
    Real acoustic_cfl = 0.6;          /**< Doc131 D5: AcousticTimeStepCK CFL (default = library default) */
    Real advection_cfl = 0.25;        /**< Doc131 D5: AdvectionTimeStepCK CFL (default = library default) */
    bool area_normalize = false;
    /** Doc207 B/C: legacy-shell (default) | wall-temperature (series-gap candidate; not production default). */
    std::string thermal_wall_model = "legacy-shell";
    bool freeze_flow = false;
    bool legacy_rad_clamp = false;
    bool legacy_thermal_sort = false; /**< diagnostic: do NOT make thermal fields follow ParticleSortCK (pre-Doc131 behaviour) */
    std::string power_mode = kPowerDiscrete;
    std::string q_shape = "uniform";   /**< uniform | radial_skin | em_rz_csv | em_particles_csv */
    std::string q_from_rz_csv;         /**< Doc134 E: prior EM Q(r,z) CSV; sets q_shape=em_rz_csv */
    std::string q_from_em_particles_csv; /**< Doc156: EM particle Q cloud (r_m,z_m,Q_Wpm3,vol_m3) */
    Real em_q_support_radius = 0.03;   /**< kernel support for em_particles_csv [m]; ~1.5*dp_em */
    Real state_record_interval = 0.0;  /**< Doc133: VTP output interval [s of physical time]; 0 = off */
    /** Doc183: lightweight particle CSV at listed PhysicalTime values (and after EM updates if listed). */
    std::vector<Real> field_export_times;
    std::string field_export_dir; /**< default: <output_dir>/field_exports */
    /** Doc151/Doc156: file-driven T→σ→L1 Picard→Q reload (subprocess). off|subprocess */
    std::string em_update_mode = "off";
    Real em_update_interval = 1.0;     /**< physical-time interval when --em-update-times unset */
    long max_em_updates = 2;           /**< Doc151 short test: at most two refreshes */
    std::string em_update_helper;      /**< python helper path */
    std::string em_solver_bin;         /**< Picard binary path */
    Real em_diag_T_perturb_K = 0.0;    /**< formal Doc156 must stay 0; diagnostic only if allow */
    Real em_dp = 0.02;                 /**< Doc156: fixed EM lattice dp (init=updates) */
    Real em_j_tol = 0.01;
    Real em_I_A = 237.341;
    std::vector<Real> em_update_times; /**< Doc156 absolute PhysicalTime schedule, e.g. 0.8,1.5 */
    bool em_end_flush = false;         /**< Doc156: default off; Doc151 used end flush */
    std::string em_particles_csv;      /**< fixed EM OriginalID lattice for spatial σ */
    std::string em_sigma_reduction = "spatial_field"; /**< spatial_field | volume_mean */
    std::string em_q_delivery = "rz_table"; /**< rz_table | em_particles — init+updates must match */
    std::string em_schedule = "absolute"; /**< Doc164: absolute | periodic */
    Real em_update_origin = 0.0;          /**< Doc164 periodic: next = origin + n*interval */
    std::string material_mode; /**< cep_anchor_constant_control | cep_thesis_NUM_shape_anchored | legacy_CES_shape_control */
    /** Doc161 P3 diagnostic: wall-cold/core-hot T init once (default off: T_wall<=0). */
    Real diag_init_T_wall = 0.0;
    Real diag_init_T_core = 0.0;
    Real diag_init_T_delta = 0.04; /**< blend thickness [m] */
    std::string diag_q_sample_rz; /**< Doc164 B1: write C++ samples at listed r,z CSV */
    std::string init_T_rz_csv;       /**< Doc164 E: r_m,z_m,T_K CSV initial temperature (fresh only) */
    std::string init_T_particle_csv; /**< Doc192: OriginalID,T_K (optional x,y,z) particle-wise T (fresh only) */
    std::string init_T_source_note;  /**< provenance label for init T CSV */
    /** Doc203 F1: re-apply T from init-T-rz-csv at current (r,z) each advection (fixed spatial field). */
    bool fixed_spatial_T = false;
    /** Doc203 F1: after Shepard sample, subtract particle-mean ΔT (R2-like zero net buoyancy). */
    bool fixed_spatial_T_demean = false;
    /**
     * Doc203 F1: keep T=Tref through settle; activate fixed-space T (CSV and/or R2) only after settle.
     * Isolates hydrostatic/acoustic packing from the circulation driver.
     */
    bool fixed_T_after_settle = false;
    /**
     * Doc203 F1: after settle (and before activating deferred T), zero all fluid velocities
     * and reset Density→ρ0 (DensityChangeRate=0). Needed because packing leaves
     * δρ/ρ~O(1%) which re-excites U~c0·δρ/ρ even after velocity clear — swamping √(βΔTgH).
     */
    bool zero_velocity_after_settle = false;
    /**
     * Doc203 F1 diagnostic: with --zero-velocity-after-settle, also set Mass=ρ0·Vol.
     * Breaks mass-vs-init ledger (mass_ok may fail) but removes Mass–ρ–Vol mismatch that
     * re-feeds the acoustic floor (§1e U~0.07 vs §1f keepM U~0.18).
     */
    bool quiet_reset_mass = false;
    std::string em_replay_dir;       /**< Doc169: freeze dir with uN/particles_rz.csv (+meta); skip Picard */
    bool em_initial_consistent = false; /**< Doc171: solve/reload Q at t=0 from current T (u0); not a scheduled event */
    long restart_step = 0;           /**< Doc164: RestartIO step to load (0=fresh) */
    Real checkpoint_interval = 0.0;  /**< Doc164: write RestartIO every Δt (0=off; EM events still write if >0 or force) */
    bool write_em_checkpoints = false; /**< Doc164: write restart + em_state.json at EM events / interval */
    bool drift_probe = false;          /**< Doc133 D1: per-substep host accumulation of d(rho)/dt (read-only) */
    bool probe_freeze_velocity = false; /**< Doc133 D1 probe: zero Velocity after the acoustic 1st half (isolates dissipative terms) */
    bool cref_buoyancy = false;        /**< Doc134: production ReferenceBodyAcceleration := Boussinesq a_b each advection step */
    bool probe_cref_buoyancy = false;  /**< alias of --cref-buoyancy that also enables --drift-probe */
    /** Doc186 U2: init ρ from uniform-buoyancy balance p'=ρ0 g β(T−Tref)(ζ−H); does not change Tref. */
    bool uniform_buoyancy_balance_init = false;
    /**
     * Doc200 B1: stratified p' balances force gβ(T−Tmean) instead of (T−Tref):
     * p'(ζ)=ρ0 g β [½ a (ζ²−H²)+(Tref−Tmean)(ζ−H)], p'(H)=0. Opt-in; default keeps Doc198/199.
     */
    bool doc200_strat_init_vs_tmean = false;
    /**
     * Doc201: with --thaw-T, default drops strat_p_prime_owns_buoyancy (ForcePrior back on).
     * Opt-in keep owns under thaw for short-case diagnostic only.
     */
    bool doc201_keep_strat_owns_when_thaw = false;
    /**
     * Doc201 B1 diagnostic ablation (NOT production): skip one update family while keeping positions live.
     * off | visc | buoyancy | wall | acoustic
     */
    std::string doc201_b1_ablate = "off";
    /**
     * Doc201 B2: live acoustic top shear — positions move; wall_closed pressure lid + viscous omit on top.
     * IC: ux=U0 sin(πζ/(2H)), uy=uz=0. Analytic decay assumes infinite plate; metrics use core r-band.
     */
    bool doc201_b2_live_top_shear = false;
    Real doc201_b2_U0 = 0.1;
    /**
     * Doc202 H1: omit acoustic wall pairs on side/bottom (AcousticWallContactOmit=1);
     * lid keeps free-slip impermeability. Diagnostic only — not production default.
     */
    bool doc202_b2_ablate_side_acoustic = false;
    /**
     * Doc202 H1b: omit acoustic wall pairs on lid only; side/bottom keep legacy noslip Riemann.
     * Mutually exclusive with ablate-side. Diagnostic only.
     */
    bool doc202_b2_ablate_lid_acoustic = false;
    /** Doc202 H1c-a: omit lid AcousticStep1stHalf pairs only (pressure force ∥ n). */
    bool doc202_b2_ablate_lid_acoustic_1st = false;
    /** Doc202 H1c-b: omit lid AcousticStep2ndHalf pairs only (n-reflect + DissipativePJump·n). */
    bool doc202_b2_ablate_lid_acoustic_2nd = false;
    /**
     * Doc202 H2a: paired dens-channel silence on lid — skip 2nd-half continuity dens;
     * keep 1st force∥n + 2nd DissipativePJump·n. Complements existing 1st-half UJump omit on FS.
     */
    bool doc202_b2_lid_continuity_silent = false;
    /** Doc202 H2b: omit DissipativePJump on lid; keep continuity dens + 1st force∥n. */
    bool doc202_b2_lid_omit_pjump = false;
    /** AcousticDensMode: 0=legacy, 1=p' dens (default), 2=Doc186 wall mom omits C_ext, 3=Doc201 omit p dens jump. */
    long acoustic_dens_mode = 1;
    /**
     * Doc186: absorb mass-mean buoyancy into BackgroundPressure; Boussinesq ForcePrior uses T_mean
     * (EOS/material Tref unchanged). Removes uniform-offset bulk drift without zeroing Tref.
     */
    bool absorb_mean_buoyancy_bg = false;
    /** Doc189 P0: same Boussinesq object, change force Tref, no flow advance. */
    bool doc189_dynamic_force_repro = false;
    /** Doc189 P3: prescribed T=Tc+A[f-<f>], reapplied from current position; energy gate N/A. */
    bool doc189_mechanism = false;
    Real doc189_mechanism_amp = 0.0;
    bool density_increment_audit = false; /**< Doc136 P0: half-step rho write budget vs actual store (read-only) */
    long density_audit_stride = 1;     /**< sample every N acoustic substeps (1 = every substep; short runs only) */
    long density_audit_max_rows = 400; /**< cap detailed CSV rows (global+IDs); 0 = unlimited */
    std::string t3_face = "top";
    std::string output_dir;
    std::string summary_json;
    std::string out_csv = "output/doc126_french_natural_reduced_pressure.csv";
    bool out_csv_set = false;
};

inline bool isThermalTestMode(const std::string &m)
{
    return m == "T0" || m == "T1" || m == "T2" || m == "T3" || m == "T4" || m == "T5";
}

/** Strict double parse: non-empty, whole string consumed, no leading space, finite, within Real range. */
inline bool parseRealStrict(const std::string &s, double &out)
{
    if (s.empty() || std::isspace(static_cast<unsigned char>(s[0])))
        return false;
    errno = 0;
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0' || errno == ERANGE)
        return false;
    if (!std::isfinite(v))
        return false;
    if (std::abs(v) > static_cast<double>(std::numeric_limits<Real>::max()))
        return false;
    out = v;
    return true;
}

inline bool parseIntStrict(const std::string &s, long &out)
{
    if (s.empty() || std::isspace(static_cast<unsigned char>(s[0])))
        return false;
    errno = 0;
    char *end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0' || errno == ERANGE)
        return false;
    out = v;
    return true;
}

inline std::string prefixOutputDir(const std::string &dir, const std::string &path)
{
    if (path.empty() || dir.empty())
        return path;
    if (std::filesystem::path(path).is_absolute())
    {
        std::cerr << "[doc131] note: absolute path '" << path << "' is not prefixed by --output-dir" << std::endl;
        return path;
    }
    return (std::filesystem::path(dir) / path).string();
}

bool validateCli(const Cli &c)
{
    auto fail = [](const std::string &msg) {
        std::cerr << "[doc131] invalid config: " << msg << std::endl;
        return false;
    };
    if (!(c.dp > 0))
        return fail("--dp must be > 0");
    if (!(c.c0 > 0))
        return fail("--c0 must be > 0");
    if (!(c.cp > 0))
        return fail("--cp must be > 0");
    if (!(c.rho0 > 0))
        return fail("--rho0 must be > 0");
    if (!(c.glass_radius > 0))
        return fail("--glass-radius (R) must be > 0");
    if (!(c.glass_height > 0))
        return fail("--glass-height (H) must be > 0");
    if (!(c.k > 0))
        return fail("--k must be > 0");
    if (!(c.u_ref > 0))
        return fail("--u-ref must be > 0");
    if (!(c.t_ref > 0) || !(c.t_cool > 0) || !(c.t_ambient > 0) || !(c.t_rad_ambient > 0))
        return fail("--t-ref, --t-cool, --t-ambient, --t-rad-ambient must be > 0 K");
    if (!(c.mu_cap > 0))
        return fail("--mu-cap must be > 0");
    if (c.settle_time < 0 || c.end_time < 0)
        return fail("--settle-time and --end-time must be >= 0");
    if (c.thermal_dt < 0)
        return fail("--thermal-dt must be >= 0 (0 = mode default)");
    if (c.emissivity < 0 || c.emissivity > 1)
        return fail("--emissivity must be in [0,1]");
    if (!c.q_from_rz_csv.empty() && c.q_shape != "uniform" && c.q_shape != "em_rz_csv")
        return fail("--q-from-rz-csv conflicts with --q-shape (use uniform/em_rz_csv or omit q-shape)");
    if (!c.q_from_em_particles_csv.empty() && c.q_shape != "uniform" && c.q_shape != "em_particles_csv")
        return fail("--q-from-em-particles-csv conflicts with --q-shape");
    if (c.q_shape != "uniform" && c.q_shape != "radial_skin" && c.q_shape != "em_rz_csv" &&
        c.q_shape != "em_particles_csv")
        return fail("--q-shape must be uniform|radial_skin|em_rz_csv|em_particles_csv");
    if ((c.q_shape == "em_rz_csv") != !c.q_from_rz_csv.empty())
        return fail("--q-shape=em_rz_csv and --q-from-rz-csv=path must be used together (or omit both)");
    if ((c.q_shape == "em_particles_csv") != !c.q_from_em_particles_csv.empty())
        return fail("--q-shape=em_particles_csv and --q-from-em-particles-csv=path must be used together");
    if (!c.q_from_rz_csv.empty() && !c.q_from_em_particles_csv.empty())
        return fail("use only one of --q-from-rz-csv or --q-from-em-particles-csv");
    if (c.q_shape == "em_particles_csv" && !(c.em_q_support_radius > 0))
        return fail("--em-q-support-radius must be > 0 for em_particles_csv");
    if (c.em_q_delivery != "rz_table" && c.em_q_delivery != "em_particles")
        return fail("--em-q-delivery must be rz_table|em_particles");
    if (c.em_schedule != "absolute" && c.em_schedule != "periodic")
        return fail("--em-schedule must be absolute|periodic");
    if (c.em_schedule == "periodic")
    {
        if (!c.em_update_times.empty())
            return fail("--em-schedule=periodic conflicts with --em-update-times (use one)");
        if (!(c.em_update_interval > 0))
            return fail("--em-schedule=periodic requires --em-update-interval>0");
        if (c.em_update_origin < 0)
            return fail("--em-update-origin must be >= 0");
    }
    if (c.em_schedule == "absolute" && c.em_update_times.empty() && c.em_update_mode == "subprocess")
    {
        // legacy interval cadence still allowed under absolute only if times empty → treat as periodic opt-in reminder
    }
    if (c.q_shape != "uniform" && c.power_mode != kPowerDiscrete)
        return fail("--q-shape!=uniform requires --power-mode=fixed_discrete_absorbed_power");
    if (c.state_record_interval < 0)
        return fail("--state-record-interval must be >= 0 (0 = off)");
    if (c.mu < 0 || c.gravity_g < 0 || c.beta < 0 || c.q_power < 0)
        return fail("--mu, --g, --beta, --q-power must be >= 0");
    if (c.h_side < 0 || c.h_bottom < 0 || c.h_top < 0)
        return fail("--h-side/--h-bottom/--h-top must be >= 0");
    if (!(c.acoustic_cfl > 0) || c.acoustic_cfl > 1.0)
        return fail("--acoustic-cfl must be in (0, 1]");
    if (!(c.advection_cfl > 0) || c.advection_cfl > 1.0)
        return fail("--advection-cfl must be in (0, 1]");
    if (!(c.wall_thickness_factor > 0))
        return fail("--wall-thickness-factor must be > 0");
    if (c.mesh_resolution < 4)
        return fail("--mesh-resolution must be an integer >= 4");
    if (c.density_audit_stride < 1)
        return fail("--density-audit-stride must be >= 1");
    if (c.density_audit_max_rows < 0)
        return fail("--density-audit-max-rows must be >= 0 (0 = unlimited)");
    if (c.power_mode != kPowerDiscrete && c.power_mode != kPowerVolumetric)
        return fail("--power-mode must be fixed_discrete_absorbed_power|fixed_volumetric_density");
    if (c.t3_face != "top" && c.t3_face != "side" && c.t3_face != "bottom")
        return fail("--t3-face must be top|side|bottom");
    if (c.em_update_mode != "off" && c.em_update_mode != "subprocess")
        return fail("--em-update-mode must be off|subprocess");
    if (c.viscosity_mode != "explicit" && c.viscosity_mode != "implicit_split")
        return fail("--viscosity-mode must be explicit|implicit_split");
    if (c.implicit_visc_schedule != "acoustic" && c.implicit_visc_schedule != "advection")
        return fail("--implicit-visc-schedule must be acoustic|advection");
    if (c.viscosity_mode == "implicit_split" && c.variable_properties &&
        (c.mean_mu_viscous || !c.local_particle_viscosity))
        return fail("--viscosity-mode=implicit_split with --variable-properties requires particle μ "
                    "(not --mean-mu-viscous)");
    if (c.viscosity_mode == "implicit_split" && !c.variable_properties && c.mean_mu_viscous)
        return fail("--viscosity-mode=implicit_split without --variable-properties uses a constant "
                    "--mu particle field");
    if (c.em_update_mode == "subprocess")
    {
        if (c.em_update_helper.empty())
            return fail("--em-update-helper required when --em-update-mode=subprocess");
        if (c.em_solver_bin.empty())
            return fail("--em-solver-bin required when --em-update-mode=subprocess");
        if (c.max_em_updates < 0)
            return fail("--max-em-updates must be >= 0");
        if (c.em_q_delivery == "em_particles")
        {
            if (c.q_shape != "em_particles_csv" || c.q_from_em_particles_csv.empty())
                return fail("em_q_delivery=em_particles requires --q-from-em-particles-csv");
        }
        else if (c.q_shape != "em_rz_csv" || c.q_from_rz_csv.empty())
        {
            return fail("subprocess EM update requires --q-from-rz-csv (or --em-q-delivery=em_particles)");
        }
        if (!c.em_update_times.empty())
        {
            if (static_cast<long>(c.em_update_times.size()) != c.max_em_updates && c.max_em_updates > 0)
            {
                // Allow max_em_updates to follow schedule size.
            }
            for (size_t i = 0; i < c.em_update_times.size(); ++i)
            {
                if (!(c.em_update_times[i] > 0))
                    return fail("--em-update-times entries must be > 0");
                if (i > 0 && !(c.em_update_times[i] > c.em_update_times[i - 1]))
                    return fail("--em-update-times must be strictly increasing");
            }
        }
        else if (!(c.em_update_interval > 0))
        {
            return fail("--em-update-interval must be > 0 when --em-update-times unset");
        }
        if (c.em_sigma_reduction != "spatial_field" && c.em_sigma_reduction != "volume_mean")
            return fail("--em-sigma-reduction must be spatial_field|volume_mean");
        if (c.em_sigma_reduction == "spatial_field" && c.em_particles_csv.empty())
            return fail("--em-particles-csv required when --em-sigma-reduction=spatial_field");
        if (!(c.em_dp > 0))
            return fail("--em-dp must be > 0");
    }
    if (!c.material_mode.empty() && c.material_mode != "cep_anchor_constant_control" &&
        c.material_mode != "cep_thesis_NUM_shape_anchored" && c.material_mode != "legacy_CES_shape_control")
        return fail("--material-mode must be cep_anchor_constant_control|cep_thesis_NUM_shape_anchored|legacy_CES_shape_control");
    if (c.material_mode == "cep_thesis_NUM_shape_anchored")
    {
        if (!c.variable_properties)
            return fail("--material-mode=cep_thesis_NUM_shape_anchored requires --variable-properties");
        if (c.material_curve_dir.empty())
            return fail("--material-mode=cep_thesis_NUM_shape_anchored requires --material-curve-dir (no CES fallback)");
    }
    return true;
}

bool parseCli(int argc, char *argv[], Cli &cli)
{
    struct RealKey
    {
        const char *name;
        Real *dst;
    };
    const RealKey real_keys[] = {{"--dp", &cli.dp},
                                 {"--settle-time", &cli.settle_time},
                                 {"--end-time", &cli.end_time},
                                 {"--rho0", &cli.rho0},
                                 {"--g", &cli.gravity_g},
                                 {"--beta", &cli.beta},
                                 {"--t-ref", &cli.t_ref},
                                 {"--mu", &cli.mu},
                                 {"--c0", &cli.c0},
                                 {"--u-ref", &cli.u_ref},
                                 {"--a-strat", &cli.a_strat},
                                 {"--dT-radial", &cli.dT_radial},
                                 {"--perturb-u", &cli.perturb_u},
                                 {"--glass-radius", &cli.glass_radius},
                                 {"--glass-height", &cli.glass_height},
                                 {"--q-power", &cli.q_power},
                                 {"--cp", &cli.cp},
                                 {"--k", &cli.k},
                                 {"--h-side", &cli.h_side},
                                 {"--h-bottom", &cli.h_bottom},
                                 {"--h-top", &cli.h_top},
                                 {"--t-cool", &cli.t_cool},
                                 {"--emissivity", &cli.emissivity},
                                 {"--t-ambient", &cli.t_ambient},
                                 {"--t-rad-ambient", &cli.t_rad_ambient},
                                 {"--mu-cap", &cli.mu_cap},
                                 {"--marangoni-dgamma-dt", &cli.marangoni_dgamma_dt},
                                 {"--wall-thickness-factor", &cli.wall_thickness_factor},
                                 {"--thermal-dt", &cli.thermal_dt},
                                 {"--t4-amp", &cli.t4_amp},
                                 {"--acoustic-cfl", &cli.acoustic_cfl},
                                 {"--advection-cfl", &cli.advection_cfl},
                                 {"--state-record-interval", &cli.state_record_interval},
                                 {"--em-update-interval", &cli.em_update_interval},
                                 {"--em-diag-T-perturb-K", &cli.em_diag_T_perturb_K},
                                 {"--em-update-origin", &cli.em_update_origin},
                                 {"--diag-init-T-wall", &cli.diag_init_T_wall},
                                 {"--diag-init-T-core", &cli.diag_init_T_core},
                                 {"--diag-init-T-delta", &cli.diag_init_T_delta},
                                 {"--checkpoint-interval", &cli.checkpoint_interval},
                                 {"--em-dp", &cli.em_dp},
                                 {"--em-j-tol", &cli.em_j_tol},
                                 {"--em-I-A", &cli.em_I_A},
                                 {"--em-q-support-radius", &cli.em_q_support_radius},
                                 {"--viscosity-diag-dt", &cli.viscosity_diag_dt},
                                 {"--viscosity-diag-shear", &cli.viscosity_diag_shear}};

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        bool matched = false;
        for (const RealKey &rk : real_keys)
        {
            const std::string p = std::string(rk.name) + "=";
            if (arg.rfind(p, 0) == 0)
            {
                double v = 0.0;
                if (!parseRealStrict(arg.substr(p.size()), v))
                {
                    std::cerr << "[doc131] invalid numeric value for " << rk.name << ": '" << arg.substr(p.size())
                              << "' (must be a finite real number)" << std::endl;
                    return false;
                }
                *rk.dst = static_cast<Real>(v);
                if (arg.rfind("--end-time=", 0) == 0)
                    cli.end_time_set = true;
                if (arg.rfind("--t-rad-ambient=", 0) == 0)
                    cli.t_rad_ambient_set = true;
                matched = true;
                break;
            }
        }
        if (matched)
            continue;

        auto take_s = [&](const std::string &key, std::string &dst) {
            const std::string p = key + "=";
            if (arg.rfind(p, 0) == 0)
            {
                dst = arg.substr(p.size());
                matched = true;
            }
        };
        take_s("--mode", cli.mode);
        take_s("--power-mode", cli.power_mode);
        take_s("--q-shape", cli.q_shape);
        take_s("--q-from-rz-csv", cli.q_from_rz_csv);
        take_s("--q-from-em-particles-csv", cli.q_from_em_particles_csv);
        take_s("--em-q-delivery", cli.em_q_delivery);
        take_s("--t3-face", cli.t3_face);
        take_s("--output-dir", cli.output_dir);
        take_s("--summary-json", cli.summary_json);
        take_s("--literature-id", cli.literature_id);
        take_s("--material-curve-dir", cli.material_curve_dir);
        take_s("--em-update-mode", cli.em_update_mode);
        take_s("--em-update-helper", cli.em_update_helper);
        take_s("--em-solver-bin", cli.em_solver_bin);
        take_s("--em-particles-csv", cli.em_particles_csv);
        take_s("--em-sigma-reduction", cli.em_sigma_reduction);
        take_s("--material-mode", cli.material_mode);
        take_s("--em-schedule", cli.em_schedule);
        take_s("--diag-q-sample-rz", cli.diag_q_sample_rz);
        take_s("--init-T-rz-csv", cli.init_T_rz_csv);
        take_s("--init-T-particle-csv", cli.init_T_particle_csv);
        take_s("--init-T-source-note", cli.init_T_source_note);
        if (arg == "--fixed-spatial-T")
        {
            cli.fixed_spatial_T = true;
            continue;
        }
        if (arg == "--fixed-spatial-T-demean")
        {
            cli.fixed_spatial_T_demean = true;
            continue;
        }
        if (arg == "--fixed-T-after-settle")
        {
            cli.fixed_T_after_settle = true;
            continue;
        }
        if (arg == "--zero-velocity-after-settle")
        {
            cli.zero_velocity_after_settle = true;
            continue;
        }
        if (arg == "--quiet-reset-mass")
        {
            cli.quiet_reset_mass = true;
            continue;
        }
        take_s("--em-replay-dir", cli.em_replay_dir);
        take_s("--viscosity-mode", cli.viscosity_mode);
        take_s("--implicit-visc-schedule", cli.implicit_visc_schedule);
        if (arg == "--viscosity-operator-diag")
        {
            cli.viscosity_operator_diag = true;
            continue;
        }
        if (arg == "--viscosity-diag-inner-only")
        {
            cli.viscosity_diag_inner_only = true;
            continue;
        }
        if (arg == "--em-initial-consistent")
        {
            cli.em_initial_consistent = true;
            continue;
        }
        if (arg.rfind("--restart-step=", 0) == 0)
        {
            cli.restart_step = std::stol(arg.substr(std::strlen("--restart-step=")));
            continue;
        }
        if (arg == "--write-em-checkpoints")
        {
            cli.write_em_checkpoints = true;
            continue;
        }
        if (arg.rfind("--field-export-times=", 0) == 0)
        {
            cli.field_export_times.clear();
            const std::string body = arg.substr(std::strlen("--field-export-times="));
            std::stringstream ss(body);
            std::string tok;
            while (std::getline(ss, tok, ','))
            {
                if (tok.empty())
                    continue;
                double v = 0;
                if (!parseRealStrict(tok, v))
                {
                    std::cerr << "[doc183] bad --field-export-times token: " << tok << std::endl;
                    return false;
                }
                cli.field_export_times.push_back(static_cast<Real>(v));
            }
            std::sort(cli.field_export_times.begin(), cli.field_export_times.end());
            continue;
        }
        take_s("--field-export-dir", cli.field_export_dir);
        if (arg.rfind("--em-update-times=", 0) == 0)
        {
            cli.em_update_times.clear();
            const std::string body = arg.substr(18);
            std::stringstream ss(body);
            std::string tok;
            while (std::getline(ss, tok, ','))
            {
                double v = 0.0;
                if (!parseRealStrict(tok, v) || !(v > 0))
                {
                    std::cerr << "[doc156] invalid --em-update-times token: '" << tok << "'" << std::endl;
                    return false;
                }
                cli.em_update_times.push_back(static_cast<Real>(v));
            }
            if (cli.em_update_times.empty())
            {
                std::cerr << "[doc156] --em-update-times empty" << std::endl;
                return false;
            }
            cli.max_em_updates = static_cast<long>(cli.em_update_times.size());
            matched = true;
        }
        if (arg == "--em-end-flush")
        {
            cli.em_end_flush = true;
            matched = true;
        }
        if (arg == "--no-em-end-flush")
        {
            cli.em_end_flush = false;
            matched = true;
        }
        if (arg.rfind("--out-csv=", 0) == 0)
        {
            cli.out_csv = arg.substr(10);
            cli.out_csv_set = true;
            matched = true;
        }
        if (matched)
        {
            if ((arg.rfind("--mode=", 0) == 0 && cli.mode.empty()) ||
                (arg.rfind("--out-csv=", 0) == 0 && cli.out_csv.empty()) ||
                (arg.rfind("--output-dir=", 0) == 0 && cli.output_dir.empty()))
            {
                std::cerr << "[doc131] empty value for " << arg << std::endl;
                return false;
            }
            continue;
        }

        if (arg.rfind("--mesh-resolution=", 0) == 0)
        {
            long v = 0;
            if (!parseIntStrict(arg.substr(18), v) || v < 4 || v > 1000000)
            {
                std::cerr << "[doc131] invalid --mesh-resolution (integer >= 4 required)" << std::endl;
                return false;
            }
            cli.mesh_resolution = static_cast<int>(v);
            continue;
        }
        if (arg.rfind("--max-em-updates=", 0) == 0)
        {
            long v = 0;
            if (!parseIntStrict(arg.substr(17), v) || v < 0)
            {
                std::cerr << "[doc151] invalid --max-em-updates (integer >= 0)" << std::endl;
                return false;
            }
            cli.max_em_updates = v;
            continue;
        }
        if (arg.rfind("--density-audit-stride=", 0) == 0)
        {
            long v = 0;
            if (!parseIntStrict(arg.substr(23), v) || v < 1)
            {
                std::cerr << "[doc136] invalid --density-audit-stride (integer >= 1)" << std::endl;
                return false;
            }
            cli.density_audit_stride = v;
            continue;
        }
        if (arg.rfind("--density-audit-max-rows=", 0) == 0)
        {
            long v = 0;
            if (!parseIntStrict(arg.substr(25), v) || v < 0)
            {
                std::cerr << "[doc136] invalid --density-audit-max-rows (integer >= 0)" << std::endl;
                return false;
            }
            cli.density_audit_max_rows = v;
            continue;
        }
        if (arg == "--kernel-correction")
        {
            cli.kernel_correction = true;
            continue;
        }
        if (arg == "--density-summation")
        {
            cli.use_density_summation = true;
            continue;
        }
        if (arg == "--thaw-T")
        {
            cli.thaw_T_set = true;
            cli.thaw_T = true;
            continue;
        }
        if (arg == "--no-thaw-T")
        {
            cli.thaw_T_set = true;
            cli.thaw_T = false;
            continue;
        }
        if (arg == "--no-top-loss")
        {
            cli.top_loss = false;
            continue;
        }
        if (arg == "--r2-keep-mean")
        {
            cli.r2_keep_mean = true;
            continue;
        }
        if (arg == "--r2-reverse")
        {
            cli.r2_reverse = true;
            continue;
        }
        if (arg == "--r2-fixed-spatial-T")
        {
            cli.r2_fixed_spatial_T = true;
            continue;
        }
        if (arg == "--area-normalize")
        {
            cli.area_normalize = true;
            continue;
        }
        if (arg.rfind("--thermal-wall-model=", 0) == 0)
        {
            cli.thermal_wall_model = arg.substr(std::string("--thermal-wall-model=").size());
            continue;
        }
        if (arg == "--freeze-flow")
        {
            cli.freeze_flow = true;
            continue;
        }
        if (arg == "--legacy-rad-clamp")
        {
            cli.legacy_rad_clamp = true;
            continue;
        }
        if (arg == "--legacy-thermal-sort")
        {
            cli.legacy_thermal_sort = true;
            continue;
        }
        if (arg == "--drift-probe")
        {
            cli.drift_probe = true;
            continue;
        }
        if (arg == "--probe-freeze-velocity")
        {
            cli.probe_freeze_velocity = true;
            cli.drift_probe = true;
            continue;
        }
        if (arg == "--cref-buoyancy")
        {
            cli.cref_buoyancy = true;
            continue;
        }
        if (arg == "--uniform-buoyancy-balance-init")
        {
            cli.uniform_buoyancy_balance_init = true;
            continue;
        }
        if (arg == "--doc200-strat-init-vs-tmean")
        {
            cli.doc200_strat_init_vs_tmean = true;
            continue;
        }
        if (arg == "--doc201-keep-strat-owns-when-thaw")
        {
            cli.doc201_keep_strat_owns_when_thaw = true;
            continue;
        }
        if (arg.rfind("--doc201-b1-ablate=", 0) == 0)
        {
            cli.doc201_b1_ablate = arg.substr(std::strlen("--doc201-b1-ablate="));
            if (cli.doc201_b1_ablate != "off" && cli.doc201_b1_ablate != "visc" &&
                cli.doc201_b1_ablate != "buoyancy" && cli.doc201_b1_ablate != "wall" &&
                cli.doc201_b1_ablate != "acoustic")
            {
                std::cerr << "[doc201][b1] --doc201-b1-ablate must be off|visc|buoyancy|wall|acoustic"
                          << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--doc201-b2-live-top-shear")
        {
            cli.doc201_b2_live_top_shear = true;
            continue;
        }
        if (arg == "--doc202-b2-ablate-side-acoustic")
        {
            cli.doc202_b2_ablate_side_acoustic = true;
            continue;
        }
        if (arg == "--doc202-b2-ablate-lid-acoustic")
        {
            cli.doc202_b2_ablate_lid_acoustic = true;
            continue;
        }
        if (arg == "--doc202-b2-ablate-lid-acoustic-1st")
        {
            cli.doc202_b2_ablate_lid_acoustic_1st = true;
            continue;
        }
        if (arg == "--doc202-b2-ablate-lid-acoustic-2nd")
        {
            cli.doc202_b2_ablate_lid_acoustic_2nd = true;
            continue;
        }
        if (arg == "--doc202-b2-lid-continuity-silent")
        {
            cli.doc202_b2_lid_continuity_silent = true;
            continue;
        }
        if (arg == "--doc202-b2-lid-omit-pjump")
        {
            cli.doc202_b2_lid_omit_pjump = true;
            continue;
        }
        if (arg.rfind("--doc201-b2-U0=", 0) == 0)
        {
            double v = 0.0;
            if (!parseRealStrict(arg.substr(std::strlen("--doc201-b2-U0=")), v) || !(v > 0.0))
            {
                std::cerr << "[doc201][b2] --doc201-b2-U0 must be a positive Real" << std::endl;
                return false;
            }
            cli.doc201_b2_U0 = static_cast<Real>(v);
            continue;
        }
        if (arg.rfind("--acoustic-dens-mode=", 0) == 0)
        {
            long v = 0;
            if (!parseIntStrict(arg.substr(21), v) || v < 0 || v > 3)
            {
                std::cerr << "[doc186/doc201] invalid --acoustic-dens-mode (integer 0|1|2|3)" << std::endl;
                return false;
            }
            cli.acoustic_dens_mode = v;
            continue;
        }
        if (arg == "--absorb-mean-buoyancy-bg")
        {
            cli.absorb_mean_buoyancy_bg = true;
            continue;
        }
        if (arg == "--doc189-dynamic-force-repro")
        {
            cli.doc189_dynamic_force_repro = true;
            continue;
        }
        if (arg.rfind("--doc189-mechanism-amp=", 0) == 0)
        {
            double v = 0.0;
            if (!parseRealStrict(arg.substr(23), v) || !std::isfinite(v))
            {
                std::cerr << "[doc189] invalid --doc189-mechanism-amp (Kelvin)" << std::endl;
                return false;
            }
            cli.doc189_mechanism = true;
            cli.doc189_mechanism_amp = static_cast<Real>(v);
            cli.thaw_T = false;
            cli.thaw_T_set = true;
            continue;
        }
        if (arg == "--probe-cref-buoyancy")
        {
            cli.cref_buoyancy = true;
            cli.probe_cref_buoyancy = true;
            cli.drift_probe = true;
            continue;
        }
        if (arg == "--density-increment-audit")
        {
            cli.density_increment_audit = true;
            continue;
        }
        if (arg == "--variable-properties")
        {
            cli.variable_properties = true;
            continue;
        }
        if (arg == "--local-particle-viscosity")
        {
            cli.local_particle_viscosity = true;
            cli.mean_mu_viscous = false;
            continue;
        }
        if (arg == "--mean-mu-viscous")
        {
            cli.mean_mu_viscous = true;
            cli.local_particle_viscosity = false;
            continue;
        }
        if (arg == "--marangoni-verify")
        {
            cli.marangoni_verify = true;
            cli.marangoni_verify_set = true;
            continue;
        }
        if (arg == "--no-marangoni-verify")
        {
            cli.marangoni_verify = false;
            cli.marangoni_verify_set = true;
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            std::cout
                << "Doc126/Doc131/Doc139 French natural convection, reduced pressure (closed flat lid) + thermal model\n"
                << "  --mode=R0|R1|R2|R3|baseline|T0|T1|T2|T3|T4|T5 (default R0)\n"
                << "  --dp= --settle-time= --end-time= --rho0= --g= --beta= --t-ref= --mu= --c0= --u-ref=\n"
                << "  --a-strat= (R1/T2 K/m)  --dT-radial= (R2/T2 K)  --perturb-u= (R3 m/s)  --t4-amp= (T4 K)\n"
                << "  --glass-radius= --glass-height= --mesh-resolution= --wall-thickness-factor=\n"
                << "  --q-power= (W)  --cp=1150 --k=4 --emissivity= --t-ambient=300 --t-rad-ambient= --t-cool=300\n"
                << "  --h-side=300 --h-bottom=35 --h-top=20  --variable-properties --mu-cap=1e5 --literature-id=\n"
                << "  --material-curve-dir=path (thesis III.8=k / III.9=cp / III.10=μ digitized CSVs)\n"
                << "  --local-particle-viscosity (default with --variable-properties) | --mean-mu-viscous (legacy)\n"
                << "  --viscosity-mode=explicit|implicit_split (Doc175; default explicit; implicit_split=CPU "
                   "DampingPairwise; constant --mu allowed without --variable-properties)\n"
                << "  --viscosity-operator-diag (Doc177: compare ViscousForceCK vs DampingPairwise then exit)\n"
                << "  --viscosity-diag-dt= --viscosity-diag-shear= (Doc177 operator-diag knobs)\n"
                << "  --marangoni-dgamma-dt=0 (Doc139 B2; thesis −1e-4 N/(m K); 0=off)\n"
                << "  --marangoni-verify | --no-marangoni-verify (axisymmetric ∇_s T bin-FD direction+area integral)\n"
                << "  --no-top-loss\n"
                << "  --power-mode=fixed_discrete_absorbed_power|fixed_volumetric_density\n"
                << "  --q-shape=uniform|radial_skin (effective_source_study) | em_rz_csv (via --q-from-rz-csv)\n"
                << "  --q-from-rz-csv=path (prior EM Q(r,z) CSV: ir,iz,r_m,z_m,Q_Wpm3; discrete-renormalised to P)\n"
                << "  --state-record-interval= (s; VTP of fluid to ./output; 0=off)\n"
                << "  --field-export-times=t0,t1,... (Doc183 CSV dumps; e.g. 0,0.7,1.4,2.1,2.3)\n"
                << "  --field-export-dir=path (default <output-dir>/field_exports)\n"
                << "  --cref-buoyancy (candidate: a_ref:=F_bouss/m; buoyancy-aware wall dens dissipation)\n"
                << "  --uniform-buoyancy-balance-init (Doc186: ρ from p'=ρ0 g β(T-Tref)(ζ-H); keep Tref)\n"
                << "  --doc200-strat-init-vs-tmean (Doc200: R1 p' balances gβ(T-Tmean); opt-in)\n"
                << "  --doc201-keep-strat-owns-when-thaw (Doc201: keep ForcePrior=0 under --thaw-T; diagnostic)\n"
                << "  --doc201-b1-ablate=off|visc|buoyancy|wall|acoustic (Doc201 diagnostic; skip one update)\n"
                << "  --doc201-b2-live-top-shear (Doc201: live positions + acoustic lid + top viscous omit shear)\n"
                << "  --doc201-b2-U0= (Doc201 B2 shear amplitude m/s; default 0.1)\n"
                << "  --doc202-b2-ablate-side-acoustic (Doc202 H1: omit side/bottom acoustic wall pairs; lid keeps FS)\n"
                << "  --doc202-b2-ablate-lid-acoustic (Doc202 H1b: omit lid acoustic pairs; side/bottom keep noslip)\n"
                << "  --doc202-b2-ablate-lid-acoustic-1st (Doc202 H1c-a: omit lid 1st-half only)\n"
                << "  --doc202-b2-ablate-lid-acoustic-2nd (Doc202 H1c-b: omit lid 2nd-half only)\n"
                << "  --doc202-b2-lid-continuity-silent (Doc202 H2a: lid skip 2nd continuity dens; keep PJump)\n"
                << "  --doc202-b2-lid-omit-pjump (Doc202 H2b: lid omit DissipativePJump; keep continuity)\n"
                << "  --acoustic-dens-mode=0|1|2|3 (default 1; 2=wall mom omits C_ext; 3=omit p dens jump / Doc201 strat)\n"
                << "  --absorb-mean-buoyancy-bg (Doc186: mean buoyancy → BackgroundPressure; force Tref:=T_mean)\n"
                << "  --doc189-dynamic-force-repro (force-only same-object Tref sequence; no flow)\n"
                << "  --doc189-mechanism-amp=A (K; prescribed radial T, reapplied; energy gate N/A)\n"
                << "  --drift-probe | --probe-freeze-velocity | --probe-cref-buoyancy (D1 probes; last aliases --cref-buoyancy)\n"
                << "  --density-increment-audit (Doc136 P0: half-step rho write budget; OriginalID samples)\n"
                << "  --density-audit-stride=N --density-audit-max-rows=N (audit sampling; defaults 1 / 400)\n"
                << "  --area-normalize (normalise B3 face areas to pi R^2 / 2 pi R H)\n"
                << "  --thermal-wall-model=legacy-shell|wall-temperature (Doc207 B; default legacy-shell)\n"
                << "  --freeze-flow (forced for T0-T4)  --t3-face=top|side|bottom  --thermal-dt= (frozen T tests)\n"
                << "  --legacy-thermal-sort (diagnostic: thermal fields NOT permuted by ParticleSortCK, pre-Doc131)\n"
                << "  --legacy-rad-clamp (baseline: old host top-only radiation with max(T,Tamb) clamp, no conduction)\n"
                << "  --thaw-T | --no-thaw-T (baseline: default thaw)\n"
                << "  --kernel-correction  --density-summation (Internal regularization; diagnostic)\n"
                << "  --r2-keep-mean (R2: keep raw dT*r/R instead of zero-mean)\n"
                << "  --r2-reverse (R2: negate dT -> cold wall)  --r2-fixed-spatial-T (R2 diagnostic: re-apply T(r) at the\n"
                << "      current position every advection step; default = Lagrangian carry of the initial T)\n"
                << "  --fixed-spatial-T (Doc203 F1: with --init-T-rz-csv, re-sample T(r,z) from CSV each advection)\n"
                << "  --fixed-spatial-T-demean (F1: subtract particle-mean ΔT after CSV sample; R2-like zero net buoyancy)\n"
                << "  --fixed-T-after-settle (F1: T=Tref through settle; activate fixed-space T after settle)\n"
                << "  --zero-velocity-after-settle (F1: clear acoustic/pseudo U after settle before observe)\n"
                << "  --quiet-reset-mass (F1 diagnostic: with zero-velocity-after-settle, also Mass=rho0*Vol)\n"
                << "  --acoustic-cfl=0.6 --advection-cfl=0.25 (time-step CFL numbers, (0,1])\n"
                << "  --output-dir= (prefix + create dir for summary-json/out-csv)  --summary-json= --out-csv=\n";
            return false;
        }
        std::cerr << "[doc131] unknown arg: " << arg << std::endl;
        return false;
    }

    const bool mode_ok = cli.mode == "R0" || cli.mode == "R1" || cli.mode == "R2" || cli.mode == "R3" ||
                         cli.mode == "baseline" || isThermalTestMode(cli.mode);
    if (!mode_ok)
    {
        std::cerr << "[doc131] unknown --mode=" << cli.mode << " (R0|R1|R2|R3|baseline|T0|T1|T2|T3|T4|T5)"
                  << std::endl;
        return false;
    }
    if (!cli.q_from_rz_csv.empty() && (cli.q_shape == "uniform" || cli.q_shape == "em_rz_csv"))
        cli.q_shape = "em_rz_csv";
    if (!cli.q_from_em_particles_csv.empty() &&
        (cli.q_shape == "uniform" || cli.q_shape == "em_particles_csv"))
        cli.q_shape = "em_particles_csv";
    if (cli.q_shape == "em_particles_csv")
        cli.em_q_delivery = "em_particles";
    else if (cli.q_shape == "em_rz_csv" && cli.em_q_delivery == "em_particles")
    {
        std::cerr << "[doc156] --em-q-delivery=em_particles requires --q-from-em-particles-csv" << std::endl;
        return false;
    }

    if (!validateCli(cli))
        return false;

    if (!cli.thaw_T_set)
        cli.thaw_T = (cli.mode == "baseline");
    // Doc139: radiation ambient may differ from convection ambient; default keeps historical single-Ta behaviour.
    if (!cli.t_rad_ambient_set)
        cli.t_rad_ambient = cli.t_ambient;

    // Mode defaults for the thermal tests (only when the user did not set them).
    if (isThermalTestMode(cli.mode) && !cli.end_time_set)
    {
        if (cli.mode == "T2")
            cli.end_time = 600.0;
        else if (cli.mode == "T4")
            cli.end_time = 1000.0;
        else if (cli.mode == "T3")
            cli.end_time = 0.2;
        else
            cli.end_time = 1.0;
    }
    if (cli.thermal_dt <= 0)
    {
        if (cli.mode == "T2")
            cli.thermal_dt = 1.0;
        else if (cli.mode == "T4")
            cli.thermal_dt = 5.0;
        else if (cli.mode == "T3")
            cli.thermal_dt = 0.01;
        else
            // Doc177: flowing advection steps must subcycle thermal; 0.05 left-Riemann was ~3.5% cold ledger.
            cli.thermal_dt = 0.005;
    }

    if (!cli.output_dir.empty())
    {
        if (!cli.out_csv_set)
            cli.out_csv = "doc131_" + cli.mode + ".csv";
        if (cli.summary_json.empty())
            cli.summary_json = "summary.json";
        cli.out_csv = prefixOutputDir(cli.output_dir, cli.out_csv);
        cli.summary_json = prefixOutputDir(cli.output_dir, cli.summary_json);
        std::error_code ec;
        std::filesystem::create_directories(cli.output_dir, ec);
        if (ec)
        {
            std::cerr << "[doc131] cannot create --output-dir=" << cli.output_dir << ": " << ec.message()
                      << std::endl;
            return false;
        }
    }
    return true;
}

//----------------------------------------------------------------------------------------------------
// Closed solid wall definition
//----------------------------------------------------------------------------------------------------
SolidBody &defineClosedCylinderSolid(SolidBody &body)
{
    body.defineAdaptationRatios(1.3, 1.0);
    body.defineMatterMaterial<Solid>();
    (void)body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    body.generateParticles<BaseParticles, Lattice>();
    return body;
}

/** Base hydrostatic only: p_h = rho0 g (H - zeta), zeta = z - z_bottom. */
inline Real analyticBaseHydroPressure(Real zeta, Real H, Real rho0, Real g)
{
    return rho0 * g * (H - zeta);
}

/**
 * Stratification perturbation pressure (p' only): p'(zeta) = 1/2 rho0 g beta a (zeta^2 - H^2),
 * p'(H)=0 and dp'/dz = rho0 g beta (T - Tref) balances Boussinesq buoyancy when GravityForceCK is omitted.
 */
inline Real analyticStratPerturbationPressure(Real zeta, Real H, Real rho0, Real g, Real beta, Real a)
{
    return Real(0.5) * rho0 * g * beta * a * (zeta * zeta - H * H);
}

/**
 * Doc200: p' balances force gβ(T−Tmean) with T=Tref+aζ:
 * p'(ζ)=ρ0 g β [½ a (ζ²−H²)+(Tref−Tmean)(ζ−H)], p'(H)=0.
 */
inline Real analyticStratPerturbationPressureVsTmean(Real zeta, Real H, Real rho0, Real g, Real beta, Real a,
                                                    Real t_ref, Real t_mean)
{
    return rho0 * g * beta *
           (Real(0.5) * a * (zeta * zeta - H * H) + (t_ref - t_mean) * (zeta - H));
}

/**
 * Doc186 uniform buoyancy balance (p' only):
 *   ∂p'/∂z = ρ0 g β (T−Tref)  balances Fz/m = β(T−Tref)g with g_acc=(0,0,-g), GravityForceCK omitted.
 *   p'(ζ)=ρ0 g β (T−Tref)(ζ−H), p'(H)=0. Then ρ=ρ0+p'/c0². Mass synced once.
 */
inline void hostApplyUniformBuoyancyBalanceDensity(BaseParticles &particles, Real rho0, Real c0, Real g, Real beta,
                                                   Real t_ref, Real z_bottom, Real H, bool sync_mass = true)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "Temperature");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *T = particles.getVariableDataByName<Real>("Temperature");
    const size_t n = particles.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    double pmin = 1e300, pmax = -1e300;
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real dT = T[i] - t_ref;
        const Real p = rho0 * g * beta * dT * (zeta - H);
        pmin = std::min(pmin, static_cast<double>(p));
        pmax = std::max(pmax, static_cast<double>(p));
        rho[i] = rho0 + p * inv_c2;
        if (sync_mass)
            mass[i] = rho[i] * vol[i];
    }
    syncVariableToDevice<Real>(particles, "Density");
    if (sync_mass)
        syncVariableToDevice<Real>(particles, "Mass");
    std::cout << "[doc186][uniform-buoy-p-init] p'_range=[" << pmin << "," << pmax << "] rho0=" << rho0
              << " beta=" << beta << " t_ref=" << t_ref << std::endl;
}

/** Init only: rho = rho0 + p'/c0^2 (perturbation_only=true -> Doc105/126 reduced pressure). Mass set once. */
inline void hostApplyStratifiedHydrostaticDensity(BaseParticles &particles, Real rho0, Real c0, Real g, Real beta,
                                                  Real a, Real z_bottom, Real H, bool sync_mass = true)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p = analyticStratPerturbationPressure(zeta, H, rho0, g, beta, a);
        rho[i] = rho0 + p * inv_c2;
        if (sync_mass)
            mass[i] = rho[i] * vol[i];
    }
    syncVariableToDevice<Real>(particles, "Density");
    if (sync_mass)
        syncVariableToDevice<Real>(particles, "Mass");
}

/**
 * Doc200 B1: same as stratified hydro init but p' matches force gβ(T−Tmean).
 * Tmean = mass-weighted Temperature already prescribed on particles.
 */
inline void hostApplyStratifiedHydrostaticDensityVsTmean(BaseParticles &particles, Real rho0, Real c0, Real g,
                                                        Real beta, Real a, Real t_ref, Real z_bottom, Real H,
                                                        bool sync_mass = true)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "Temperature");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *T = particles.getVariableDataByName<Real>("Temperature");
    const size_t n = particles.TotalRealParticles();
    double m_sum = 0.0, mT = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double mi = static_cast<double>(rho0 * vol[i]); // pre-hydro mass proxy (uniform ρ0 Vol)
        m_sum += mi;
        mT += mi * static_cast<double>(T[i]);
    }
    const Real t_mean = (m_sum > TinyReal) ? static_cast<Real>(mT / m_sum) : t_ref;
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p = analyticStratPerturbationPressureVsTmean(zeta, H, rho0, g, beta, a, t_ref, t_mean);
        rho[i] = rho0 + p * inv_c2;
        if (sync_mass)
            mass[i] = rho[i] * vol[i];
    }
    syncVariableToDevice<Real>(particles, "Density");
    if (sync_mass)
        syncVariableToDevice<Real>(particles, "Mass");
    std::cout << "[doc200][b1] strat_init_vs_tmean Tref=" << t_ref << " Tmean=" << t_mean
              << " (p' balances gβ(T-Tmean))" << std::endl;
}

/** Max |c0^2 (rho - rho0) - p'_expected|; p'_expected = strat p' (R1) or 0 (others). */
inline void hostEosIdentityResidual(BaseParticles &particles, Real rho0, Real c0, Real g, Real beta, Real a,
                                    Real z_bottom, Real H, bool expect_strat, double &max_res, double &p_scale,
                                    bool vs_tmean = false, Real t_ref = Real(0), Real t_mean = Real(0))
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    max_res = 0.0;
    p_scale = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const double p_an =
            !expect_strat
                ? 0.0
                : (vs_tmean ? static_cast<double>(analyticStratPerturbationPressureVsTmean(zeta, H, rho0, g, beta, a,
                                                                                          t_ref, t_mean))
                            : static_cast<double>(analyticStratPerturbationPressure(zeta, H, rho0, g, beta, a)));
        const double p_num = static_cast<double>(c0) * static_cast<double>(c0) *
                             (static_cast<double>(rho[i]) - static_cast<double>(rho0));
        max_res = std::max(max_res, std::abs(p_num - p_an));
        p_scale = std::max(p_scale, std::abs(p_an));
    }
}

/** T(z)=Tref + a (z - z_bottom). */
inline void hostPrescribeStableStratification(BaseParticles &particles, Real t_ref, Real a, Real z_bottom)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real Ti = t_ref + a * (pos[i][2] - z_bottom);
        temp[i] = Ti;
        delta_t[i] = Ti - t_ref;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/**
 * R2 radial driver: T = Tref + dT (r/R - s_mean). Hot at the wall (r->R), cold on the axis.
 * By default s_mean = particle mean of r/R (zero net buoyancy => pure circulation drive, no net pressure shift).
 */
inline Real hostPrescribeRadialTemperature(BaseParticles &particles, const Vecd &center, Real R, Real t_ref,
                                           Real dT, bool remove_mean, const double *fixed_s_mean = nullptr)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double s_sum = 0.0;
    if (fixed_s_mean == nullptr)
    {
        for (size_t i = 0; i < n; ++i)
        {
            const double dx = static_cast<double>(pos[i][0] - center[0]);
            const double dy = static_cast<double>(pos[i][1] - center[1]);
            s_sum += std::min(1.0, std::sqrt(dx * dx + dy * dy) / (static_cast<double>(R) + 1e-30));
        }
    }
    // fixed_s_mean (D2 --r2-fixed-spatial-T): the constant removed at init is reused so the spatial field T(r) is
    // time-independent (the particle-mean of r/R drifts slightly as particles move).
    const Real s_mean = fixed_s_mean != nullptr
                            ? static_cast<Real>(*fixed_s_mean)
                            : (remove_mean ? static_cast<Real>(s_sum / static_cast<double>(std::max<size_t>(n, 1)))
                                           : Real(0));
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real s = std::min(Real(1), std::sqrt(dx * dx + dy * dy) / (R + TinyReal));
        const Real Ti = t_ref + dT * (s - s_mean);
        temp[i] = Ti;
        delta_t[i] = Ti - t_ref;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    return s_mean;
}

/** T2 non-uniform start: Delta T = a (z - z_mid) + dT (r/R - <r/R>)  (delta-authoritative; T = t_ref + dT). */
inline void hostPrescribeT2Field(BaseParticles &particles, const Vecd &center, Real R, Real H, Real z_bottom,
                                 Real t_ref, Real a_z, Real dT_radial)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double s_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        s_sum += std::min(1.0, std::sqrt(dx * dx + dy * dy) / (static_cast<double>(R) + 1e-30));
    }
    const double s_mean = s_sum / static_cast<double>(std::max<size_t>(n, 1));
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double s = std::min(1.0, std::sqrt(dx * dx + dy * dy) / (static_cast<double>(R) + 1e-30));
        const double dT = static_cast<double>(a_z) * (static_cast<double>(pos[i][2] - z_bottom) - 0.5 * H) +
                          static_cast<double>(dT_radial) * (s - s_mean);
        delta_t[i] = static_cast<Real>(dT);
        temp[i] = t_ref + delta_t[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/** T4 start: Delta T = A cos(pi zeta / H), uniform in r (1D profile). */
/** Doc161: diagnostic wall-cold / core-hot axisymmetric T (init once; not a literature field). */
inline void hostPrescribeWallColdCoreHotTemperature(BaseParticles &particles, const Vecd &center, Real R, Real t_ref,
                                                    Real T_wall, Real T_core, Real delta)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    const double Rd = static_cast<double>(R);
    const double dd = std::max(static_cast<double>(delta), 1.0e-6);
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double s = std::min(1.0, std::max(0.0, r / Rd));
        const double w = 0.5 * (1.0 + std::tanh((s - (1.0 - dd / Rd)) / std::max(0.35 * dd / Rd, 1.0e-3)));
        const Real Ti = static_cast<Real>((1.0 - w) * static_cast<double>(T_core) + w * static_cast<double>(T_wall));
        temp[i] = Ti;
        delta_t[i] = Ti - t_ref;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}


/** Doc192: set Temperature/DeltaT by OriginalID from a field export (OriginalID,T_K). */
inline bool hostPrescribeTemperatureFromParticleCsv(BaseParticles &particles, Real t_ref,
                                                    const std::string &csv_path, std::string &err)
{
    std::ifstream in(csv_path);
    if (!in)
    {
        err = "cannot open init-T-particle-csv";
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        err = "empty init-T-particle-csv";
        return false;
    }
    auto strip = [](std::string s) {
        while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t'))
            s.pop_back();
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t'))
            ++i;
        return s.substr(i);
    };
    auto lower = [&](std::string s) {
        s = strip(s);
        for (char &c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };
    std::vector<std::string> hdr;
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
            hdr.push_back(lower(cell));
    }
    int i_id = -1, i_T = -1;
    for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
    {
        if (hdr[static_cast<size_t>(i)] == "originalid" || hdr[static_cast<size_t>(i)] == "id")
            i_id = i;
        else if (hdr[static_cast<size_t>(i)] == "t_k" || hdr[static_cast<size_t>(i)] == "t" ||
                 hdr[static_cast<size_t>(i)] == "temperature" || hdr[static_cast<size_t>(i)] == "temperature_k")
            i_T = i;
    }
    if (i_id < 0 || i_T < 0)
    {
        err = "init-T-particle-csv needs OriginalID,T_K";
        return false;
    }
    std::unordered_map<UnsignedInt, Real> by_id;
    while (std::getline(in, line))
    {
        line = strip(line);
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> cols;
        while (std::getline(ss, cell, ','))
            cols.push_back(strip(cell));
        if (static_cast<int>(cols.size()) <= std::max(i_id, i_T))
        {
            err = "short row in init-T-particle-csv";
            return false;
        }
        double id_d = 0.0, T_d = 0.0;
        if (!parseRealStrict(cols[static_cast<size_t>(i_id)], id_d) ||
            !parseRealStrict(cols[static_cast<size_t>(i_T)], T_d) || !(T_d > 0.0) || !std::isfinite(T_d))
        {
            err = "bad OriginalID/T in init-T-particle-csv";
            return false;
        }
        by_id[static_cast<UnsignedInt>(id_d + 0.5)] = static_cast<Real>(T_d);
    }
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    size_t hit = 0;
    for (size_t i = 0; i < n; ++i)
    {
        auto it = by_id.find(oid[i]);
        if (it == by_id.end())
        {
            err = "missing OriginalID=" + std::to_string(oid[i]) + " in init-T-particle-csv";
            return false;
        }
        temp[i] = it->second;
        delta[i] = temp[i] - t_ref;
        ++hit;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    std::cout << "[doc192][init-T-particle-csv] path=" << csv_path << " n=" << n << " matched=" << hit
              << std::endl;
    return hit == n;
}

/** Doc164 / Doc203 F1: axisymmetric T(r,z) nodes from CSV (r_m,z_m,T_K). */
struct RzTempNode
{
    double r = 0, z = 0, T = 0;
};

inline bool hostLoadRzTemperatureCsv(const std::string &csv_path, std::vector<RzTempNode> &nodes, std::string &err)
{
    nodes.clear();
    std::ifstream in(csv_path);
    if (!in)
    {
        err = "cannot open init-T-rz-csv";
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        err = "empty init-T-rz-csv";
        return false;
    }
    auto strip = [](std::string s) {
        while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t'))
            s.pop_back();
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t'))
            ++i;
        return s.substr(i);
    };
    auto lower = [&](std::string s) {
        s = strip(s);
        for (char &c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };
    std::vector<std::string> hdr;
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
            hdr.push_back(lower(cell));
    }
    int i_r = -1, i_z = -1, i_T = -1;
    for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
    {
        if (hdr[static_cast<size_t>(i)] == "r_m" || hdr[static_cast<size_t>(i)] == "r")
            i_r = i;
        else if (hdr[static_cast<size_t>(i)] == "z_m" || hdr[static_cast<size_t>(i)] == "z")
            i_z = i;
        else if (hdr[static_cast<size_t>(i)] == "t_k" || hdr[static_cast<size_t>(i)] == "t" ||
                 hdr[static_cast<size_t>(i)] == "temperature_k")
            i_T = i;
    }
    if (i_r < 0 || i_z < 0 || i_T < 0)
    {
        err = "init-T-rz-csv needs r_m,z_m,T_K";
        return false;
    }
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(ss, cell, ','))
            cells.push_back(cell);
        if (static_cast<int>(cells.size()) <= std::max(i_r, std::max(i_z, i_T)))
            continue;
        double r = 0, z = 0, Tv = 0;
        for (std::string &c : cells)
            c = strip(c);
        if (!parseRealStrict(cells[static_cast<size_t>(i_r)], r) ||
            !parseRealStrict(cells[static_cast<size_t>(i_z)], z) ||
            !parseRealStrict(cells[static_cast<size_t>(i_T)], Tv))
        {
            err = "bad numeric row in init-T-rz-csv";
            return false;
        }
        if (!(std::isfinite(r) && std::isfinite(z) && std::isfinite(Tv)) || r < 0 || !(Tv > 0))
        {
            err = "non-finite/non-positive T or bad r in init-T-rz-csv";
            return false;
        }
        nodes.push_back({r, z, Tv});
    }
    if (nodes.empty())
    {
        err = "no usable rows in init-T-rz-csv";
        return false;
    }
    return true;
}

/** Shepard / nearest sample of cached T(r,z) nodes onto particles (Doc203 F1 re-applicable). */
inline bool hostApplyRzTemperatureNodes(BaseParticles &particles, const Vecd &center, Real t_ref,
                                        const std::vector<RzTempNode> &nodes, Real support_h, std::string &err,
                                        bool quiet = false, bool demean_delta_t = false)
{
    if (nodes.empty())
    {
        err = "empty T(r,z) node set";
        return false;
    }
    const double h = std::max(static_cast<double>(support_h), 1e-6);
    const double h2 = h * h;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    size_t n_fallback = 0, n_ok = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        double num = 0, den = 0;
        for (const RzTempNode &nd : nodes)
        {
            const double d2 = (r - nd.r) * (r - nd.r) + (z - nd.z) * (z - nd.z);
            if (d2 > h2)
                continue;
            const double s = std::sqrt(d2) / (h + 1e-300);
            const double w = (1.0 - s) * (1.0 - s);
            num += w * nd.T;
            den += w;
        }
        double Ti;
        if (den > 0)
        {
            Ti = num / den;
            ++n_ok;
        }
        else
        {
            double bd = 1e300;
            Ti = nodes.front().T;
            for (const RzTempNode &nd : nodes)
            {
                const double d2 = (r - nd.r) * (r - nd.r) + (z - nd.z) * (z - nd.z);
                if (d2 < bd)
                {
                    bd = d2;
                    Ti = nd.T;
                }
            }
            ++n_fallback;
        }
        temp[i] = static_cast<Real>(Ti);
        delta_t[i] = static_cast<Real>(Ti) - t_ref;
    }
    double demean_shift = 0.0;
    if (demean_delta_t && n > 0)
    {
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i)
            sum += static_cast<double>(delta_t[i]);
        demean_shift = sum / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i)
        {
            delta_t[i] = static_cast<Real>(static_cast<double>(delta_t[i]) - demean_shift);
            temp[i] = t_ref + delta_t[i];
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    if (!quiet)
        std::cout << "[doc164][T-rz] n_nodes=" << nodes.size() << " h=" << h << " sampled_ok=" << n_ok
                  << " nearest_fallback=" << n_fallback << " demean=" << (demean_delta_t ? 1 : 0)
                  << " demean_shift_K=" << demean_shift << std::endl;
    if (n_fallback > n / 2)
    {
        err = "T-rz coverage insufficient (nearest fallback > 50%)";
        return false;
    }
    return true;
}

/** Doc164 E: load + apply once (fresh init). */
inline bool hostPrescribeTemperatureFromRzCsv(BaseParticles &particles, const Vecd &center, Real t_ref,
                                              const std::string &csv_path, Real support_h, std::string &err,
                                              std::vector<RzTempNode> *nodes_out = nullptr,
                                              bool demean_delta_t = false, bool apply_now = true)
{
    std::vector<RzTempNode> nodes;
    if (!hostLoadRzTemperatureCsv(csv_path, nodes, err))
        return false;
    if (apply_now)
    {
        if (!hostApplyRzTemperatureNodes(particles, center, t_ref, nodes, support_h, err, /*quiet=*/false,
                                         demean_delta_t))
            return false;
    }
    std::cout << "[doc164][init-T-rz-csv] path=" << csv_path << " apply_now=" << (apply_now ? 1 : 0)
              << " demean=" << (demean_delta_t ? 1 : 0) << std::endl;
    if (nodes_out)
        *nodes_out = std::move(nodes);
    return true;
}

inline void hostPrescribeT4Field(BaseParticles &particles, Real H, Real z_bottom, Real t_ref, Real amp)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const double zeta = std::min(static_cast<double>(H), std::max(0.0, static_cast<double>(pos[i][2] - z_bottom)));
        delta_t[i] = static_cast<Real>(static_cast<double>(amp) * std::cos(kPiD * zeta / static_cast<double>(H)));
        temp[i] = t_ref + delta_t[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/** BackgroundPressure = rho0 g (H - zeta) (base only; does not include stratification). */
inline void hostApplyBackgroundHydroPressure(BaseParticles &particles, Real rho0, Real g, Real z_bottom, Real H)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        p_bg[i] = analyticBaseHydroPressure(zeta, H, rho0, g);
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
}

/**
 * Doc186: p_h = ρ0 g (H−ζ) + ρ0 g β (T_mean−Tref_eos)(ζ−H).
 * Absorbs the uniform/mean Boussinesq offset like gravity; force then uses T_mean.
 */
inline void hostApplyBackgroundHydroPressureWithMeanBuoyancy(BaseParticles &particles, Real rho0, Real g, Real beta,
                                                             Real t_mean, Real t_ref_eos, Real z_bottom, Real H)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    const Real dT_mean = t_mean - t_ref_eos;
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p_base = analyticBaseHydroPressure(zeta, H, rho0, g);
        const Real p_mean_buoy = rho0 * g * beta * dT_mean * (zeta - H);
        p_bg[i] = p_base + p_mean_buoy;
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
}

/**
 * Doc189/Doc192: non-uniform T, same Boussinesq object, sequential force-Tref updates.
 * Also seeds a non-buoyancy ForcePrior component and checks ForcePrior delta / no double-count.
 * Does not advance position/velocity. Writes JSON + per-particle CSV.
 */
template <class BoussinesqDynamics>
int runDoc189DynamicForceRepro(BoussinesqDynamics &dyn, BaseParticles &particles, Real beta, Real g,
                               Real radius, const Vecd &axis_center, Real eos_tref, const std::string &out_dir)
{
    namespace fs = std::filesystem;
    fs::create_directories(out_dir);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    Vecd *fp0 = particles.getVariableDataByName<Vecd>("ForcePrior");
    const size_t n = particles.TotalRealParticles();
    const Real A = Real(25);
    const Real seed_fx = Real(1.0e-3); // non-buoyancy ForcePrior component (x)
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_center[0];
        const Real dy = pos[i][1] - axis_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real f = Real(1) - Real(2) * (r / (radius + TinyReal)) * (r / (radius + TinyReal));
        temp[i] = eos_tref + A * f;
        dT[i] = temp[i] - eos_tref;
        fp0[i] = Vecd(seed_fx, Real(0), Real(0));
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToDevice<Vecd>(particles, "ForcePrior");

    double wT = 0.0, w = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        wT += static_cast<double>(mass[i]) * static_cast<double>(temp[i]);
        w += static_cast<double>(mass[i]);
    }
    const Real Ta = static_cast<Real>(w > 0.0 ? wT / w : 0.0);
    const std::vector<Real> trefs = {Ta, Ta + Real(10), Ta + Real(10), Ta - Real(5), Ta};
    const char *tags[] = {"Ta", "Ta_plus10", "repeat_Ta_plus10", "Ta_minus5", "return_Ta"};

    std::ofstream csv(out_dir + "/per_particle_Fz.csv");
    csv << "step,tag,tref_requested,id,mass,T,Fz_field,ForcePrior_z,Fz_analytic,ForcePrior_x\n";
    std::ofstream js(out_dir + "/DYNAMIC_FORCE_REPRO.json");
    js << std::setprecision(17);
    js << "{\n  \"n\": " << n << ",\n  \"beta\": " << static_cast<double>(beta)
       << ",\n  \"g\": " << static_cast<double>(g) << ",\n  \"Ta\": " << static_cast<double>(Ta)
       << ",\n  \"eos_tref_used_only_for_T_init\": " << static_cast<double>(eos_tref)
       << ",\n  \"seed_ForcePrior_x\": " << static_cast<double>(seed_fx)
       << ",\n  \"execution_policy\": \""
#if SPHINXSYS_USE_SYCL
       << "ParallelDevicePolicy"
#else
       << "ParallelPolicy"
#endif
       << "\",\n  \"steps\": [\n";

    bool stale = false;
    bool forceprior_ok = true;
    double max_fp_seed_err = 0.0;
    double max_fp_delta_err = 0.0;
    double max_fp_repeat_err = 0.0;
    std::vector<Vecd> prev_fb(n, Vecd::Zero());
    std::vector<Vecd> prev_fp(n, Vecd::Zero());
    for (size_t s = 0; s < trefs.size(); ++s)
    {
        dyn.setForceReferenceTemperature(trefs[s]);
        dyn.exec();
        syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
        syncVariableToHost<Vecd>(particles, "ForcePrior");
        const Vecd *fb = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
        const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
        double sum_fz = 0.0, sum_m = 0.0, sum_e2 = 0.0, sum_fref2 = 0.0, max_abs = 0.0;
        double sum_fpz = 0.0;
        double step_seed_err = 0.0, step_delta_err = 0.0, step_repeat_err = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            const double fz_ref = static_cast<double>(mass[i]) * static_cast<double>(beta) *
                                  static_cast<double>(g) * (static_cast<double>(temp[i]) - static_cast<double>(trefs[s]));
            const double fz = static_cast<double>(fb[i][2]);
            const double err = fz - fz_ref;
            sum_fz += fz;
            sum_m += static_cast<double>(mass[i]);
            sum_e2 += err * err;
            sum_fref2 += fz_ref * fz_ref;
            max_abs = std::max(max_abs, std::abs(err));
            sum_fpz += static_cast<double>(fp[i][2]);
            step_seed_err = std::max(step_seed_err, std::abs(static_cast<double>(fp[i][0] - seed_fx)));
            if (s == 0)
            {
                // First exec: ForcePrior_z should become buoyancy (previous=0).
                step_delta_err =
                    std::max(step_delta_err, std::abs(static_cast<double>(fp[i][2] - fb[i][2])));
            }
            else if (tags[s] == std::string("repeat_Ta_plus10"))
            {
                step_repeat_err = std::max(
                    step_repeat_err, static_cast<double>((fp[i] - prev_fp[i]).norm()));
            }
            else
            {
                // ForcePrior change must equal buoyancy current-previous.
                const Vecd d_fp = fp[i] - prev_fp[i];
                const Vecd d_fb = fb[i] - prev_fb[i];
                step_delta_err = std::max(step_delta_err, static_cast<double>((d_fp - d_fb).norm()));
            }
            if (i < 8 || i + 8 >= n)
            {
                csv << s << "," << tags[s] << "," << static_cast<double>(trefs[s]) << "," << i << ","
                    << static_cast<double>(mass[i]) << "," << static_cast<double>(temp[i]) << "," << fz << ","
                    << static_cast<double>(fp[i][2]) << "," << fz_ref << "," << static_cast<double>(fp[i][0])
                    << "\n";
            }
            prev_fb[i] = fb[i];
            prev_fp[i] = fp[i];
        }
        max_fp_seed_err = std::max(max_fp_seed_err, step_seed_err);
        max_fp_delta_err = std::max(max_fp_delta_err, step_delta_err);
        max_fp_repeat_err = std::max(max_fp_repeat_err, step_repeat_err);
        const double inferred =
            (sum_m * static_cast<double>(beta) * static_cast<double>(g) != 0.0)
                ? (wT / w) - sum_fz / (sum_m * static_cast<double>(beta) * static_cast<double>(g))
                : 0.0;
        const bool step_stale = std::abs(inferred - static_cast<double>(trefs[s])) > 0.05;
        stale = stale || (s > 0 && step_stale);
        if (s)
            js << ",\n";
        js << "    {\"tag\": \"" << tags[s] << "\", \"tref_requested\": " << static_cast<double>(trefs[s])
           << ", \"host_member\": " << static_cast<double>(dyn.getForceReferenceTemperature())
           << ", \"tref_inferred_from_sumFz\": " << inferred << ", \"sum_Fz\": " << sum_fz
           << ", \"sum_ForcePrior_z\": " << sum_fpz << ", \"L2_rel\": " << std::sqrt(sum_e2 / std::max(sum_fref2, 1e-30))
           << ", \"max_abs_particle_err_N\": " << max_abs << ", \"step_looks_stale\": " << (step_stale ? "true" : "false")
           << ", \"max_ForcePrior_seed_x_err\": " << step_seed_err
           << ", \"max_ForcePrior_delta_err\": " << step_delta_err
           << ", \"max_ForcePrior_repeat_err\": " << step_repeat_err << "}";
        std::cout << "[doc189][repro] " << tags[s] << " requested=" << trefs[s] << " inferred=" << inferred
                  << " sumFz=" << sum_fz << " max_abs=" << max_abs << " seed_err=" << step_seed_err
                  << " delta_err=" << step_delta_err << " repeat_err=" << step_repeat_err << std::endl;
    }
    const double fp_tol = 1.0e-8;
    forceprior_ok = (max_fp_seed_err <= fp_tol) && (max_fp_delta_err <= fp_tol) && (max_fp_repeat_err <= fp_tol);
    js << "\n  ],\n  \"gpu_kept_first_reference\": " << (stale ? "true" : "false")
       << ",\n  \"forceprior_independent_passed\": " << (forceprior_ok ? "true" : "false")
       << ",\n  \"forceprior_max_seed_x_err\": " << max_fp_seed_err
       << ",\n  \"forceprior_max_delta_err\": " << max_fp_delta_err
       << ",\n  \"forceprior_max_repeat_err\": " << max_fp_repeat_err
       << ",\n  \"forceprior_tol_N\": " << fp_tol << "\n}\n";
    std::cout << "[doc189][repro] wrote " << out_dir << " stale=" << stale
              << " forceprior_ok=" << forceprior_ok << std::endl;
    if (stale)
        return 2;
    return forceprior_ok ? 0 : 3;
}

/** Mass-weighted mean Temperature (host). */
inline Real hostMassWeightedMeanTemperature(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const size_t n = particles.TotalRealParticles();
    double wT = 0.0, w = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double m = static_cast<double>(mass[i]);
        wT += m * static_cast<double>(T[i]);
        w += m;
    }
    return static_cast<Real>(w > 0.0 ? wT / w : 0.0);
}

/** Doc189: T=Tc+A[f-<f>_m], f=1-2(r/R)^2. Rewrites Temperature and ThermalDeltaT. */
inline void hostApplyDoc189MechanismTemperature(BaseParticles &particles, const Vecd &axis_center, Real radius,
                                                Real t_eos, Real amplitude)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double wf = 0.0, w = 0.0;
    std::vector<Real> f(n, Real(0));
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_center[0];
        const Real dy = pos[i][1] - axis_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real rr = r / (radius + TinyReal);
        f[i] = Real(1) - Real(2) * rr * rr;
        wf += static_cast<double>(mass[i]) * static_cast<double>(f[i]);
        w += static_cast<double>(mass[i]);
    }
    const Real fmean = static_cast<Real>(w > 0.0 ? wf / w : 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        temp[i] = t_eos + amplitude * (f[i] - fmean);
        delta[i] = temp[i] - t_eos;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/** Acoustic momentum uses p' when flag=1. */
inline void hostSetPerturbationMomentum(BaseParticles &particles, Real flag)
{
    syncVariableToHost<Real>(particles, "PerturbationMomentum");
    Real *pm = particles.getVariableDataByName<Real>("PerturbationMomentum");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        pm[i] = flag;
    syncVariableToDevice<Real>(particles, "PerturbationMomentum");
}

/** 0=legacy dens, 1=p' dens with independent wall reference. */
inline void hostSetAcousticDensMode(BaseParticles &particles, Real mode)
{
    syncVariableToHost<Real>(particles, "AcousticDensMode");
    Real *m = particles.getVariableDataByName<Real>("AcousticDensMode");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        m[i] = mode;
    syncVariableToDevice<Real>(particles, "AcousticDensMode");
}

/**
 * Doc201 B2: tag flat-lid wall particles so AcousticStep2ndHalf dens ghost uses free-slip
 * (match viscous omit). Side/bottom stay no-slip dens mirror. Requires AcousticWallFreeSlip
 * already registered (AcousticStep2ndHalf Contact<Wall> ctor) and NormalDirection computed.
 */
inline size_t hostTagAcousticWallFreeSlipLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Vecd>(wall_particles, "NormalDirection");
    syncVariableToHost<Real>(wall_particles, "AcousticWallFreeSlip");
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *fs = wall_particles.getVariableDataByName<Real>("AcousticWallFreeSlip");
    const size_t nn = wall_particles.TotalRealParticles();
    const Real z_cut = z_top - Real(0.25) * wall_thickness;
    size_t n_lid = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        fs[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_lid;
    }
    syncVariableToDevice<Real>(wall_particles, "AcousticWallFreeSlip");
    return n_lid;
}

/**
 * Doc202 H1: tag non-lid wall particles AcousticWallContactOmit=1 so 1st/2nd wall
 * interact skips those pairs. Lid (omit=0) keeps free-slip impermeability.
 * Requires AcousticWallContactOmit registered (AcousticStep Contact<Wall> ctor).
 */
inline size_t hostTagAcousticWallContactOmitSideBottom(BaseParticles &wall_particles, Real z_top,
                                                       Real wall_thickness)
{
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Vecd>(wall_particles, "NormalDirection");
    syncVariableToHost<Real>(wall_particles, "AcousticWallContactOmit");
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *omit = wall_particles.getVariableDataByName<Real>("AcousticWallContactOmit");
    const size_t nn = wall_particles.TotalRealParticles();
    const Real z_cut = z_top - Real(0.25) * wall_thickness;
    size_t n_omit = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        omit[i] = is_lid ? Real(0) : Real(1);
        if (!is_lid)
            ++n_omit;
    }
    syncVariableToDevice<Real>(wall_particles, "AcousticWallContactOmit");
    return n_omit;
}

/** Doc202 H1b: omit lid acoustic pairs only (omit=1 on lid; side/bottom =0). */
inline size_t hostTagAcousticWallContactOmitLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Vecd>(wall_particles, "NormalDirection");
    syncVariableToHost<Real>(wall_particles, "AcousticWallContactOmit");
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *omit = wall_particles.getVariableDataByName<Real>("AcousticWallContactOmit");
    const size_t nn = wall_particles.TotalRealParticles();
    const Real z_cut = z_top - Real(0.25) * wall_thickness;
    size_t n_omit = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        omit[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_omit;
    }
    syncVariableToDevice<Real>(wall_particles, "AcousticWallContactOmit");
    return n_omit;
}

/** Tag lid particles with a Real wall field =1 (helper for Doc202 H2a/H2b). */
inline size_t hostTagLidWallRealField(BaseParticles &wall_particles, Real z_top, Real wall_thickness,
                                      const char *var_name)
{
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Vecd>(wall_particles, "NormalDirection");
    syncVariableToHost<Real>(wall_particles, var_name);
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *f = wall_particles.getVariableDataByName<Real>(var_name);
    const size_t nn = wall_particles.TotalRealParticles();
    const Real z_cut = z_top - Real(0.25) * wall_thickness;
    size_t n_tag = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        f[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_tag;
    }
    syncVariableToDevice<Real>(wall_particles, var_name);
    return n_tag;
}

/** Doc202 H2a: tag lid AcousticWallContSilent=1 (skip 2nd continuity dens; keep PJump). */
inline size_t hostTagAcousticWallContSilentLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealField(wall_particles, z_top, wall_thickness, "AcousticWallContSilent");
}

/** Doc202 H2b: tag lid AcousticWallOmitPJump=1. */
inline size_t hostTagAcousticWallOmitPJumpLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealField(wall_particles, z_top, wall_thickness, "AcousticWallOmitPJump");
}

/**
 * Doc202 H1c: tag lid with half-specific omit (1st or 2nd). which="1st"|"2nd".
 * Does not set AcousticWallContactOmit (full); side/bottom stay active.
 */
inline size_t hostTagAcousticWallContactOmitLidHalf(BaseParticles &wall_particles, Real z_top,
                                                    Real wall_thickness, const char *which)
{
    const bool tag_1st = (std::strcmp(which, "1st") == 0);
    const bool tag_2nd = (std::strcmp(which, "2nd") == 0);
    const char *var = tag_1st ? "AcousticWallContactOmit1st" : "AcousticWallContactOmit2nd";
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Vecd>(wall_particles, "NormalDirection");
    syncVariableToHost<Real>(wall_particles, var);
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *omit = wall_particles.getVariableDataByName<Real>(var);
    const size_t nn = wall_particles.TotalRealParticles();
    const Real z_cut = z_top - Real(0.25) * wall_thickness;
    size_t n_omit = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        omit[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_omit;
    }
    syncVariableToDevice<Real>(wall_particles, var);
    (void)tag_2nd;
    return n_omit;
}

/** Reference body acceleration for wall dens C_ref (0: gravity is not in ForcePrior). */
inline void hostSetReferenceBodyAcceleration(BaseParticles &particles, const Vecd &a_ref)
{
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    Vecd *a = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        a[i] = a_ref;
    syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
}

inline double hostTotalMass(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "Mass");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    double m = 0.0;
    for (size_t i = 0; i < n; ++i)
        m += static_cast<double>(mass[i]);
    return m;
}

/**
 * R3: divergence-free impermeable meridional roll (axisymmetric), amplitude eps (m/s):
 *   s=r/R, xi=zeta/H
 *   u_r = -eps (pi R/H) s (1-s)^2 cos(pi xi),   u_z = 2 eps (1-s)(1-2s) sin(pi xi).
 * Cartesian: u_x = -eps (pi R/H) (x/R)(1-s)^2 cos(pi xi) (no r->0 singularity).
 */
inline void hostApplyR3VelocityPerturbation(BaseParticles &particles, const Vecd &center, Real R, Real H,
                                            Real z_bottom, Real eps)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real s = std::min(Real(1), std::sqrt(dx * dx + dy * dy) / (R + TinyReal));
        const Real xi = std::min(Real(1), std::max(Real(0), (pos[i][2] - z_bottom) / (H + TinyReal)));
        const Real one_ms = Real(1) - s;
        const Real cx = -eps * (Pi * R / H) * (one_ms * one_ms) * std::cos(Pi * xi) / R;
        vel[i][0] += cx * dx;
        vel[i][1] += cx * dy;
        vel[i][2] += Real(2) * eps * one_ms * (Real(1) - Real(2) * s) * std::sin(Pi * xi);
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
}

/** Doc200/Doc201 top free-slip shear: u_x = U0 sin(π ζ /(2H)) exp[-ν (π/(2H))² t]. */
inline Real analyticDoc201B2TopFreeSlipUx(Real zeta, Real t, Real U0, Real H, Real nu)
{
    const Real k = Pi / (Real(2) * H);
    return U0 * std::sin(k * zeta) * std::exp(-nu * k * k * t);
}

/** Doc203 F1: hard-clear fluid velocity (host+device) after packing / settle. */
inline void hostZeroFluidVelocity(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        vel[i] = Vecd::Zero();
    syncVariableToDevice<Vecd>(particles, "Velocity");
}

/**
 * Doc203 F1: reset Density→ρ0 so residual packing δρ cannot re-excite acoustic U~c0·δρ/ρ.
 * If reset_mass: also Mass=ρ0·Vol (diagnostic IC; mass_ok vs init may fail).
 * DensityChangeRate cleared when present.
 */
inline void hostResetDensityToRho0(BaseParticles &particles, Real rho0, bool reset_mass = false)
{
    syncVariableToHost<Real>(particles, "Density");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        rho[i] = rho0;
    syncVariableToDevice<Real>(particles, "Density");
    if (reset_mass)
    {
        syncVariableToHost<Real>(particles, "Mass");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        Real *mass = particles.getVariableDataByName<Real>("Mass");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        for (size_t i = 0; i < n; ++i)
            mass[i] = rho0 * vol[i];
        syncVariableToDevice<Real>(particles, "Mass");
    }
    if (particles.getVariableDataByName<Real>("DensityChangeRate") != nullptr)
    {
        syncVariableToHost<Real>(particles, "DensityChangeRate");
        Real *dr = particles.getVariableDataByName<Real>("DensityChangeRate");
        for (size_t i = 0; i < n; ++i)
            dr[i] = Real(0);
        syncVariableToDevice<Real>(particles, "DensityChangeRate");
    }
}

/** Doc201 B2 IC: replace velocity with pure top free-slip shear (no add). Positions untouched. */
inline void hostApplyDoc201B2TopShearIC(BaseParticles &particles, Real z_bottom, Real H, Real U0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        vel[i] = Vecd(analyticDoc201B2TopFreeSlipUx(zeta, Real(0), U0, H, Real(0)), Real(0), Real(0));
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
}

struct Doc201B2Metrics
{
    Real amp_num = 0;
    Real amp_ref = 0;
    Real amp_rel_err = 1;
    Real L2_rel = 1;
    Real uz_rms = 0;
    Real uz_maxabs = 0;
    Real uz_rms_over_U0 = 1;
    Real uy_rms = 0;
    Real K = 0;
    size_t n_interior = 0;
    Real disp_rms = 0; /**< |Δx| rms vs IC positions (live check) */
    Real rho_rel_drift = 0;
    Real p_prime_rms = 0;
};

/**
 * Core cylinder metric: exclude bottom/top incomplete-support bands and side wall band.
 * r_max_keep isolates core from side no-slip contamination (plan §9 wide-domain central observation).
 */
inline Doc201B2Metrics hostEvaluateDoc201B2TopShear(BaseParticles &particles, const Vecd &center, Real R,
                                                    Real H, Real z_bottom, Real t, Real U0, Real nu, Real dp,
                                                    Real r_max_keep, const Vecd *pos0_opt, Real rho0, Real c0)
{
    Doc201B2Metrics m;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Pressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *p = particles.getVariableDataByName<Real>("Pressure");
    const size_t n = particles.TotalRealParticles();
    // Keep mid-column particles: Doc200 used 2.5 dp on H=0.1/dp=0.01; at coarser dp clamp so H-2m>0.
    const Real margin_z_bot = std::min(Real(2.5) * dp, Real(0.2) * H);
    const Real margin_z_top = margin_z_bot;
    const Real k = Pi / (Real(2) * H);
    double num = 0, den = 0, err2 = 0, ref2 = 0, K = 0, uz2 = 0, uy2 = 0, vol_sum = 0;
    double disp2 = 0, rho_err2 = 0, p2 = 0;
    Real uz_max = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real zeta = pos[i][2] - z_bottom;
        if (zeta < margin_z_bot || zeta > H - margin_z_top || r > r_max_keep)
            continue;
        const Real uref = analyticDoc201B2TopFreeSlipUx(zeta, t, U0, H, nu);
        const Real unum = vel[i][0];
        const Real w = vol[i];
        const Real s = std::sin(k * zeta);
        num += static_cast<double>(unum * s * w);
        den += static_cast<double>(s * s * w);
        err2 += static_cast<double>((unum - uref) * (unum - uref) * w);
        ref2 += static_cast<double>(uref * uref * w);
        K += 0.5 * static_cast<double>(unum * unum + vel[i][1] * vel[i][1] + vel[i][2] * vel[i][2]) *
             static_cast<double>(w);
        uz2 += static_cast<double>(vel[i][2] * vel[i][2] * w);
        uy2 += static_cast<double>(vel[i][1] * vel[i][1] * w);
        vol_sum += static_cast<double>(w);
        uz_max = std::max(uz_max, std::abs(vel[i][2]));
        if (pos0_opt != nullptr)
        {
            const Real ddx = pos[i][0] - pos0_opt[i][0];
            const Real ddy = pos[i][1] - pos0_opt[i][1];
            const Real ddz = pos[i][2] - pos0_opt[i][2];
            disp2 += static_cast<double>((ddx * ddx + ddy * ddy + ddz * ddz) * w);
        }
        const Real dr = (rho[i] - rho0) / std::max(std::abs(rho0), TinyReal);
        rho_err2 += static_cast<double>(dr * dr * w);
        p2 += static_cast<double>(p[i] * p[i] * w);
        ++m.n_interior;
    }
    m.amp_num = den > TinyReal ? static_cast<Real>(num / den) : 0;
    m.amp_ref = U0 * std::exp(-nu * k * k * t);
    m.amp_rel_err = std::abs(m.amp_num - m.amp_ref) /
                    std::max(std::abs(m.amp_ref), static_cast<Real>(TinyReal));
    m.L2_rel = static_cast<Real>(std::sqrt(err2) /
                                 std::max(std::sqrt(ref2), static_cast<double>(TinyReal)));
    m.K = static_cast<Real>(K);
    m.uz_rms = vol_sum > TinyReal ? static_cast<Real>(std::sqrt(uz2 / vol_sum)) : 0;
    m.uy_rms = vol_sum > TinyReal ? static_cast<Real>(std::sqrt(uy2 / vol_sum)) : 0;
    m.uz_maxabs = uz_max;
    m.uz_rms_over_U0 = m.uz_rms / std::max(std::abs(U0), static_cast<Real>(TinyReal));
    m.disp_rms = vol_sum > TinyReal ? static_cast<Real>(std::sqrt(disp2 / vol_sum)) : 0;
    m.rho_rel_drift = vol_sum > TinyReal ? static_cast<Real>(std::sqrt(rho_err2 / vol_sum)) : 0;
    m.p_prime_rms = vol_sum > TinyReal ? static_cast<Real>(std::sqrt(p2 / vol_sum)) : 0;
    (void)R;
    (void)c0;
    return m;
}

struct FlowStats
{
    size_t n = 0;
    double U_rms = 0, Uz_rms = 0, Ur_rms = 0, Utheta_rms = 0, U_max = 0, K = 0;
    double U_mean[3] = {0, 0, 0};
    double T_mean = 0, T_min = 0, T_max = 0;
    double rho_mean = 0, rho_min = 0, rho_max = 0, p_prime_rms = 0, p_prime_mean = 0;
    double circ_ur_proxy = 0; /**< <u_r sign(z-z_mid)> (meridional roll; negative: inflow aloft) */
    double circ_uz_proxy = 0; /**< <u_z (r/R - <r/R>)> (positive: upflow at the wall) */
    bool finite = true;
};

inline FlowStats hostFlowStats(BaseParticles &particles, const Vecd &center, Real R, Real H, Real z_bottom,
                               Real rho0, Real c0)
{
    FlowStats s;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    s.n = particles.TotalRealParticles();
    if (s.n == 0)
    {
        s.finite = false;
        return s;
    }
    const double inv_n = 1.0 / static_cast<double>(s.n);
    const double z_mid = static_cast<double>(z_bottom) + 0.5 * static_cast<double>(H);
    double u2 = 0, uz2 = 0, ur2 = 0, ut2 = 0, K = 0, pp2 = 0, pp = 0;
    double sum_s = 0, sum_ur_sign = 0, sum_uz_s = 0;
    double t_sum = 0, rho_sum = 0;
    s.T_min = std::numeric_limits<double>::max();
    s.T_max = -std::numeric_limits<double>::max();
    s.rho_min = std::numeric_limits<double>::max();
    s.rho_max = -std::numeric_limits<double>::max();
    std::vector<double> s_vals(s.n, 0.0);
    for (size_t i = 0; i < s.n; ++i)
    {
        const double ux = vel[i][0], uy = vel[i][1], uz = vel[i][2];
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        double ur = 0.0, ut = 0.0;
        if (r > 1.0e-9)
        {
            ur = (dx * ux + dy * uy) / r;
            ut = (dx * uy - dy * ux) / r;
        }
        const double un2 = ux * ux + uy * uy + uz * uz;
        if (!std::isfinite(un2) || !std::isfinite(static_cast<double>(rho[i])) ||
            !std::isfinite(static_cast<double>(temp[i])))
            s.finite = false;
        u2 += un2;
        uz2 += uz * uz;
        ur2 += ur * ur;
        ut2 += ut * ut;
        s.U_max = std::max(s.U_max, std::sqrt(un2));
        K += 0.5 * static_cast<double>(mass[i]) * un2;
        s.U_mean[0] += ux;
        s.U_mean[1] += uy;
        s.U_mean[2] += uz;
        const double sn = std::min(1.0, r / (static_cast<double>(R) + 1e-30));
        s_vals[i] = sn;
        sum_s += sn;
        const double zsgn = (static_cast<double>(pos[i][2]) >= z_mid) ? 1.0 : -1.0;
        sum_ur_sign += ur * zsgn;
        const double tt = static_cast<double>(temp[i]);
        t_sum += tt;
        s.T_min = std::min(s.T_min, tt);
        s.T_max = std::max(s.T_max, tt);
        const double rr = static_cast<double>(rho[i]);
        rho_sum += rr;
        s.rho_min = std::min(s.rho_min, rr);
        s.rho_max = std::max(s.rho_max, rr);
        const double p_prime = static_cast<double>(c0) * static_cast<double>(c0) * (rr - static_cast<double>(rho0));
        pp2 += p_prime * p_prime;
        pp += p_prime;
    }
    const double s_mean = sum_s * inv_n;
    for (size_t i = 0; i < s.n; ++i)
        sum_uz_s += static_cast<double>(vel[i][2]) * (s_vals[i] - s_mean);
    s.U_rms = std::sqrt(u2 * inv_n);
    s.Uz_rms = std::sqrt(uz2 * inv_n);
    s.Ur_rms = std::sqrt(ur2 * inv_n);
    s.Utheta_rms = std::sqrt(ut2 * inv_n);
    s.K = K;
    for (int d = 0; d < 3; ++d)
        s.U_mean[d] *= inv_n;
    s.T_mean = t_sum * inv_n;
    s.rho_mean = rho_sum * inv_n;
    s.p_prime_rms = std::sqrt(pp2 * inv_n);
    s.p_prime_mean = pp * inv_n;
    s.circ_ur_proxy = sum_ur_sign * inv_n;
    s.circ_uz_proxy = sum_uz_s * inv_n;
    if (!std::isfinite(s.U_rms))
        s.finite = false;
    return s;
}

//----------------------------------------------------------------------------------------------------
// Doc133 D1: pi_mean drift probe (read-only host diagnostics)
//----------------------------------------------------------------------------------------------------
struct DriftSnapshot
{
    double rho_mean = 0, mass_sum = 0, sum_rhoV = 0, sum_rhoV_over_mass = 0, vol_sum = 0;
    double drho_dt_mean = 0, drho_dt_rms = 0, drho_dt_band_contrib = 0; /**< band_contrib = sum_band drho_dt / N */
    double speed_mean = 0;
    size_t n_band = 0;
};

/** Wall band: within `band` of the side wall, bottom or top (where the Adami wall dens term acts). */
inline bool doc133InWallBand(const Vecd &p, const Vecd &center, double R, double H, double z_bottom, double band)
{
    const double dx = static_cast<double>(p[0] - center[0]);
    const double dy = static_cast<double>(p[1] - center[1]);
    const double r = std::sqrt(dx * dx + dy * dy);
    const double z = static_cast<double>(p[2]);
    return r > R - band || z < z_bottom + band || z > z_bottom + H - band;
}

inline DriftSnapshot hostDriftSnapshot(BaseParticles &particles, const Vecd &center, Real R, Real H, Real z_bottom,
                                       Real dp)
{
    DriftSnapshot d;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
        return d;
    double sr = 0, sm = 0, srv = 0, sv = 0, sd = 0, sd2 = 0, sb = 0, ss = 0;
    for (size_t i = 0; i < n; ++i)
    {
        sr += static_cast<double>(rho[i]);
        sm += static_cast<double>(mass[i]);
        srv += static_cast<double>(rho[i]) * static_cast<double>(vol[i]);
        sv += static_cast<double>(vol[i]);
        const double dd = static_cast<double>(drho[i]);
        sd += dd;
        sd2 += dd * dd;
        ss += static_cast<double>(vel[i].norm());
        if (doc133InWallBand(pos[i], center, static_cast<double>(R), static_cast<double>(H),
                             static_cast<double>(z_bottom), 2.0 * static_cast<double>(dp)))
        {
            sb += dd;
            ++d.n_band;
        }
    }
    const double inv = 1.0 / static_cast<double>(n);
    d.rho_mean = sr * inv;
    d.mass_sum = sm;
    d.sum_rhoV = srv;
    d.sum_rhoV_over_mass = srv / (sm + 1e-300);
    d.vol_sum = sv;
    d.drho_dt_mean = sd * inv;
    d.drho_dt_rms = std::sqrt(sd2 * inv);
    d.drho_dt_band_contrib = sb * inv;
    d.speed_mean = ss * inv;
    return d;
}

/** RMS of Boussinesq acceleration magnitude |F_b|/m (proves buoyancy is active). */
inline double hostBuoyancyAccelRms(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *fb = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    double s2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double a = static_cast<double>(fb[i].norm()) / (static_cast<double>(mass[i]) + 1e-300);
        s2 += a * a;
    }
    return std::sqrt(s2 / static_cast<double>(std::max<size_t>(n, 1)));
}

//----------------------------------------------------------------------------------------------------
// Doc136 P0: density half-step increment audit (read-only; tracks OriginalID)
// Production writes (AcousticStep CK):
//   1st Initialize: rho += 0.5*dt*rate_previous ; Interact: rate = dissipative (+ wall)
//   2nd Interact:   rate += continuity          ; Update:   rho += 0.5*dt*rate_current
// Doc133 cum_drho_* used mean(rate)*full_dt after each half — NOT equal to these writes.
//----------------------------------------------------------------------------------------------------
inline double doc136UlpAround(double x)
{
    if (!std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();
    // Production Real spacing around the float representation of x.
    const float xf = static_cast<float>(x);
    const float up = std::nextafter(xf, std::numeric_limits<float>::infinity());
    const float dn = std::nextafter(xf, -std::numeric_limits<float>::infinity());
    return static_cast<double>(std::max(up - xf, xf - dn));
}

struct Doc136TrackedParticle
{
    UnsignedInt original_id = 0;
    const char *zone = "unset";
    double cum_d_act = 0, cum_d_exp = 0, cum_r_update = 0;
    double cum_lost_signed = 0, cum_lost_abs = 0;
    size_t n_zero_store = 0, n_nonzero_exp = 0, n_samples = 0;
};

struct Doc136DensityAudit
{
    bool enabled = false;
    long stride = 1;
    long max_rows = 400;
    size_t n_substep = 0;
    size_t n_rows_written = 0;
    // Global double ledgers (particle-mean of half-step weights per sampled substep).
    double cum_d_act_mean = 0, cum_d_exp_mean = 0, cum_r_update_mean = 0;
    double sum_lost_signed = 0, sum_lost_abs = 0; /**< particle-sum of discarded expected increments */
    size_t n_zero_store_events = 0, n_nonzero_exp_events = 0, n_particle_steps = 0;
    double max_abs_r_update = 0;
    UnsignedInt max_r_update_oid = 0;
    std::vector<Doc136TrackedParticle> tracked;
    std::ofstream global_csv;
    std::ofstream id_csv;
    std::string global_path, id_path, note;

    void open(const std::string &outdir)
    {
        global_path = outdir + "/density_increment_audit_global.csv";
        id_path = outdir + "/density_increment_audit_ids.csv";
        global_csv.open(global_path);
        id_csv.open(id_path);
        global_csv << "substep,phase,time_s,dt_ac,sorted_flag,"
                      "rho_mean_before,rho_mean_after1,rho_mean_after2,"
                      "rate_mean_before,rate_mean_after1,rate_mean_after2,"
                      "d_act_mean,d_exp_mean,R_update_mean,"
                      "ulp_mean,frac_exp_over_ulp_mean,n_zero_store,n_nonzero_exp,n\n";
        id_csv << "substep,phase,time_s,dt_ac,zone,OriginalID,array_index,"
                  "rho_before,rho_after1,rho_after2,rate_before,rate_after1,rate_after2,"
                  "d_act_1,d_exp_1,d_act_2,d_exp_2,d_act,d_exp,R_update,"
                  "ulp_rho,abs_exp_over_ulp,zero_store\n";
    }

    void selectTracked(BaseParticles &particles, const Vecd &center, Real R, Real H, Real z_bottom, Real dp)
    {
        tracked.clear();
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<UnsignedInt>(particles, "OriginalID");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
        const size_t n = particles.TotalRealParticles();
        auto push_best = [&](const char *zone, auto score_fn) {
            size_t best = n;
            double best_s = 1e300;
            for (size_t i = 0; i < n; ++i)
            {
                const double s = score_fn(pos[i]);
                if (s < best_s)
                {
                    best_s = s;
                    best = i;
                }
            }
            if (best < n)
            {
                Doc136TrackedParticle t;
                t.original_id = oid[best];
                t.zone = zone;
                tracked.push_back(t);
            }
        };
        const double z_mid = static_cast<double>(z_bottom + Real(0.5) * H);
        const double r_in = 0.35 * static_cast<double>(R);
        push_best("interior", [&](const Vecd &p) {
            const double dx = static_cast<double>(p[0] - center[0]);
            const double dy = static_cast<double>(p[1] - center[1]);
            const double r = std::sqrt(dx * dx + dy * dy);
            return std::abs(r - r_in) + std::abs(static_cast<double>(p[2]) - z_mid);
        });
        push_best("side", [&](const Vecd &p) {
            const double dx = static_cast<double>(p[0] - center[0]);
            const double dy = static_cast<double>(p[1] - center[1]);
            const double r = std::sqrt(dx * dx + dy * dy);
            return std::abs(r - static_cast<double>(R) + 0.5 * static_cast<double>(dp)) +
                   std::abs(static_cast<double>(p[2]) - z_mid);
        });
        push_best("bottom", [&](const Vecd &p) {
            const double dx = static_cast<double>(p[0] - center[0]);
            const double dy = static_cast<double>(p[1] - center[1]);
            const double r = std::sqrt(dx * dx + dy * dy);
            return std::abs(static_cast<double>(p[2]) - static_cast<double>(z_bottom) - 0.5 * static_cast<double>(dp)) +
                   0.2 * r;
        });
        push_best("top", [&](const Vecd &p) {
            const double dx = static_cast<double>(p[0] - center[0]);
            const double dy = static_cast<double>(p[1] - center[1]);
            const double r = std::sqrt(dx * dx + dy * dy);
            return std::abs(static_cast<double>(p[2]) - static_cast<double>(z_bottom + H) + 0.5 * static_cast<double>(dp)) +
                   0.2 * r;
        });
        push_best("corner_side_bottom", [&](const Vecd &p) {
            const double dx = static_cast<double>(p[0] - center[0]);
            const double dy = static_cast<double>(p[1] - center[1]);
            const double r = std::sqrt(dx * dx + dy * dy);
            return std::abs(r - static_cast<double>(R) + 0.5 * static_cast<double>(dp)) +
                   std::abs(static_cast<double>(p[2]) - static_cast<double>(z_bottom) - 0.5 * static_cast<double>(dp));
        });
        note = "AcousticStep1stHalf::Initialize rho+=0.5*dt*rate_prev; AcousticStep2ndHalf::Update "
               "rho+=0.5*dt*rate_curr. Expected = half-step weights only (double). Actual = stored Density.";
    }

    size_t findIndexByOid(const UnsignedInt *oid, size_t n, UnsignedInt target) const
    {
        for (size_t i = 0; i < n; ++i)
            if (oid[i] == target)
                return i;
        return n;
    }

    void sampleSubstep(BaseParticles &particles, Real dt_ac, double time_s, const char *phase, bool sorted_flag,
                       const std::vector<Real> &rho_before, const std::vector<Real> &rate_before,
                       const std::vector<Real> &rho_after1, const std::vector<Real> &rate_after1,
                       const std::vector<Real> &rho_after2, const std::vector<Real> &rate_after2)
    {
        ++n_substep;
        if (!enabled)
            return;
        if ((static_cast<long>(n_substep) - 1) % stride != 0)
            return;
        if (max_rows > 0 && static_cast<long>(n_rows_written) >= max_rows)
            return;

        const size_t n = rho_before.size();
        if (n == 0 || rho_after2.size() != n)
            return;
        syncVariableToHost<UnsignedInt>(particles, "OriginalID");
        const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");

        const double dt = static_cast<double>(dt_ac);
        double sum_rb = 0, sum_r1 = 0, sum_r2 = 0, sum_rate_b = 0, sum_rate1 = 0, sum_rate2 = 0;
        double sum_d_act = 0, sum_d_exp = 0, sum_R = 0, sum_ulp = 0, sum_frac = 0;
        size_t n_zero = 0, n_nz = 0;
        for (size_t i = 0; i < n; ++i)
        {
            const double rb = static_cast<double>(rho_before[i]);
            const double r1 = static_cast<double>(rho_after1[i]);
            const double r2 = static_cast<double>(rho_after2[i]);
            const double rate_b = static_cast<double>(rate_before[i]);
            const double rate1 = static_cast<double>(rate_after1[i]);
            const double rate2 = static_cast<double>(rate_after2[i]);
            const double d_exp_1 = 0.5 * dt * rate_b;
            const double d_exp_2 = 0.5 * dt * rate2;
            const double d_act = r2 - rb;
            const double d_exp = d_exp_1 + d_exp_2;
            const double R = d_act - d_exp;
            const double ulp = doc136UlpAround(rb);
            sum_rb += rb;
            sum_r1 += r1;
            sum_r2 += r2;
            sum_rate_b += rate_b;
            sum_rate1 += rate1;
            sum_rate2 += rate2;
            sum_d_act += d_act;
            sum_d_exp += d_exp;
            sum_R += R;
            sum_ulp += ulp;
            sum_frac += (ulp > 0 ? std::abs(d_exp) / ulp : 0.0);
            const bool nz = std::abs(d_exp) > 0.0;
            const bool zero_store = nz && (r2 == rb); // exact float equality of stored Density
            if (nz)
                ++n_nz;
            if (zero_store)
            {
                ++n_zero;
                sum_lost_signed += d_exp;
                sum_lost_abs += std::abs(d_exp);
            }
            if (std::abs(R) > max_abs_r_update)
            {
                max_abs_r_update = std::abs(R);
                max_r_update_oid = oid[i];
            }
            ++n_particle_steps;
        }
        const double inv = 1.0 / static_cast<double>(n);
        cum_d_act_mean += sum_d_act * inv;
        cum_d_exp_mean += sum_d_exp * inv;
        cum_r_update_mean += sum_R * inv;
        n_zero_store_events += n_zero;
        n_nonzero_exp_events += n_nz;

        global_csv << n_substep << "," << phase << "," << time_s << "," << dt << "," << (sorted_flag ? 1 : 0) << ","
                   << sum_rb * inv << "," << sum_r1 * inv << "," << sum_r2 * inv << "," << sum_rate_b * inv << ","
                   << sum_rate1 * inv << "," << sum_rate2 * inv << "," << sum_d_act * inv << "," << sum_d_exp * inv
                   << "," << sum_R * inv << "," << sum_ulp * inv << "," << sum_frac * inv << "," << n_zero << ","
                   << n_nz << "," << n << "\n";

        for (auto &t : tracked)
        {
            const size_t i = findIndexByOid(oid, n, t.original_id);
            if (i >= n)
                continue;
            const double rb = static_cast<double>(rho_before[i]);
            const double r1 = static_cast<double>(rho_after1[i]);
            const double r2 = static_cast<double>(rho_after2[i]);
            const double rate_b = static_cast<double>(rate_before[i]);
            const double rate1 = static_cast<double>(rate_after1[i]);
            const double rate2 = static_cast<double>(rate_after2[i]);
            const double d_exp_1 = 0.5 * dt * rate_b;
            const double d_exp_2 = 0.5 * dt * rate2;
            const double d_act_1 = r1 - rb;
            const double d_act_2 = r2 - r1;
            const double d_act = r2 - rb;
            const double d_exp = d_exp_1 + d_exp_2;
            const double R = d_act - d_exp;
            const double ulp = doc136UlpAround(rb);
            const bool nz = std::abs(d_exp) > 0.0;
            const bool zero_store = nz && (r2 == rb);
            t.cum_d_act += d_act;
            t.cum_d_exp += d_exp;
            t.cum_r_update += R;
            ++t.n_samples;
            if (nz)
                ++t.n_nonzero_exp;
            if (zero_store)
            {
                ++t.n_zero_store;
                t.cum_lost_signed += d_exp;
                t.cum_lost_abs += std::abs(d_exp);
            }
            id_csv << n_substep << "," << phase << "," << time_s << "," << dt << "," << t.zone << ","
                   << t.original_id << "," << i << "," << rb << "," << r1 << "," << r2 << "," << rate_b << ","
                   << rate1 << "," << rate2 << "," << d_act_1 << "," << d_exp_1 << "," << d_act_2 << "," << d_exp_2
                   << "," << d_act << "," << d_exp << "," << R << "," << ulp << ","
                   << (ulp > 0 ? std::abs(d_exp) / ulp : 0.0) << "," << (zero_store ? 1 : 0) << "\n";
        }
        ++n_rows_written;
    }
};

inline void hostCopyDensityRate(BaseParticles &particles, std::vector<Real> &rho, std::vector<Real> &rate)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    const Real *r = particles.getVariableDataByName<Real>("Density");
    const Real *d = particles.getVariableDataByName<Real>("DensityChangeRate");
    const size_t n = particles.TotalRealParticles();
    rho.resize(n);
    rate.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        rho[i] = r[i];
        rate[i] = d[i];
    }
}

//----------------------------------------------------------------------------------------------------
// Doc131 D4: boundary-band velocity and particle escape report (read-only)
//----------------------------------------------------------------------------------------------------
struct BoundaryBandStats
{
    size_t n_top = 0, n_side = 0, n_bottom = 0;
    double top_uz_rms = 0, top_uz_maxabs = 0, top_uz_mean = 0;
    double side_ur_rms = 0, side_ur_maxabs = 0;
    double bottom_uz_rms = 0, bottom_uz_maxabs = 0;
    size_t n_outside_strict = 0;  /**< outside the geometric cylinder z<0 | z>H | r>R */
    size_t n_outside_half_dp = 0; /**< beyond the geometric cylinder by more than dp/2 */
    double z_max = 0, z_min = 0, r_max = 0;
};

inline BoundaryBandStats hostBoundaryBandStats(BaseParticles &particles, const Vecd &center, Real R, Real H,
                                               Real z_bottom, Real dp)
{
    BoundaryBandStats b;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    const double zb = static_cast<double>(z_bottom), zt = zb + static_cast<double>(H);
    const double Rd = static_cast<double>(R), band = static_cast<double>(dp), half = 0.5 * band;
    b.z_max = -std::numeric_limits<double>::max();
    b.z_min = std::numeric_limits<double>::max();
    double t2 = 0, s2 = 0, bt2 = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        const double uz = static_cast<double>(vel[i][2]);
        const double ur = r > 1.0e-9 ? (dx * static_cast<double>(vel[i][0]) + dy * static_cast<double>(vel[i][1])) / r : 0.0;
        b.z_max = std::max(b.z_max, z);
        b.z_min = std::min(b.z_min, z);
        b.r_max = std::max(b.r_max, r);
        if (z > zt || z < zb || r > Rd)
            ++b.n_outside_strict;
        if (z > zt + half || z < zb - half || r > Rd + half)
            ++b.n_outside_half_dp;
        if (z > zt - band)
        {
            ++b.n_top;
            t2 += uz * uz;
            b.top_uz_maxabs = std::max(b.top_uz_maxabs, std::abs(uz));
            b.top_uz_mean += uz;
        }
        if (r > Rd - band)
        {
            ++b.n_side;
            s2 += ur * ur;
            b.side_ur_maxabs = std::max(b.side_ur_maxabs, std::abs(ur));
        }
        if (z < zb + band)
        {
            ++b.n_bottom;
            bt2 += uz * uz;
            b.bottom_uz_maxabs = std::max(b.bottom_uz_maxabs, std::abs(uz));
        }
    }
    if (b.n_top > 0)
    {
        b.top_uz_rms = std::sqrt(t2 / static_cast<double>(b.n_top));
        b.top_uz_mean /= static_cast<double>(b.n_top);
    }
    if (b.n_side > 0)
        b.side_ur_rms = std::sqrt(s2 / static_cast<double>(b.n_side));
    if (b.n_bottom > 0)
        b.bottom_uz_rms = std::sqrt(bt2 / static_cast<double>(b.n_bottom));
    return b;
}

//----------------------------------------------------------------------------------------------------
// Doc131 A2: EOS / Pressure field semantics
//----------------------------------------------------------------------------------------------------
/** EOS as implemented by WeaklyCompressibleFluid::getPressure (same Real arithmetic). */
inline Real eosPressureFromDensity(Real rho, Real rho0, Real c0)
{
    const Real p0 = rho0 * c0 * c0;
    return p0 * (rho / rho0 - Real(1));
}

struct PiStats
{
    bool valid = false;
    size_t n = 0;
    double mean = 0, rms_demean = 0, min = 0, max = 0;
};

/** pi = EOS(Density) (perturbation pressure reconstructed from the Density field). */
inline PiStats hostPiStatsFromDensity(BaseParticles &particles, Real rho0, Real c0)
{
    PiStats s;
    syncVariableToHost<Real>(particles, "Density");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    s.n = particles.TotalRealParticles();
    if (s.n == 0)
        return s;
    double sum = 0.0;
    s.min = std::numeric_limits<double>::max();
    s.max = -std::numeric_limits<double>::max();
    std::vector<double> pi(s.n);
    for (size_t i = 0; i < s.n; ++i)
    {
        pi[i] = static_cast<double>(eosPressureFromDensity(rho[i], rho0, c0));
        sum += pi[i];
        s.min = std::min(s.min, pi[i]);
        s.max = std::max(s.max, pi[i]);
    }
    s.mean = sum / static_cast<double>(s.n);
    double v = 0.0;
    for (size_t i = 0; i < s.n; ++i)
        v += (pi[i] - s.mean) * (pi[i] - s.mean);
    s.rms_demean = std::sqrt(v / static_cast<double>(s.n));
    s.valid = std::isfinite(s.mean) && std::isfinite(s.rms_demean);
    return s;
}

/**
 * Pressure semantics (verified against AcousticStep1stHalf::InitializeKernel):
 *   Pressure = EOS(Density after the half density update) + BackgroundPressure  (TOTAL pressure p = p' + p_h).
 * It is refreshed at the start of each acoustic 1st half; the 2nd half then updates Density again, so
 * the identity Pressure - BackgroundPressure == EOS(Density) only holds right after the 1st half. The
 * sample below is therefore taken (read-only) right after the 1st half of the last production substep.
 */
struct PressureSample
{
    bool valid = false;
    double time = 0.0;
    size_t n = 0;
    PiStats pi_from_pressure; /**< Pressure - BackgroundPressure */
    PiStats pi_from_density;  /**< EOS(Density) at the same instant */
    double resid_rms = 0, resid_max = 0, eos_rms0 = 0, tol = 0, p_total_max_abs = 0, pbg_max_abs = 0;
    bool consistent = false;
};

inline void hostSamplePressureField(BaseParticles &particles, Real rho0, Real c0, double time, PressureSample &out)
{
    syncVariableToHost<Real>(particles, "Pressure");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    syncVariableToHost<Real>(particles, "Density");
    const Real *p = particles.getVariableDataByName<Real>("Pressure");
    const Real *pbg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    out = PressureSample();
    out.time = time;
    out.n = n;
    if (n == 0)
        return;
    std::vector<double> pi_p(n), pi_e(n);
    double r2 = 0.0, e2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double eos = static_cast<double>(eosPressureFromDensity(rho[i], rho0, c0));
        pi_e[i] = eos;
        pi_p[i] = static_cast<double>(p[i]) - static_cast<double>(pbg[i]);
        const double r = pi_p[i] - eos;
        r2 += r * r;
        e2 += eos * eos;
        out.resid_max = std::max(out.resid_max, std::abs(r));
        out.p_total_max_abs = std::max(out.p_total_max_abs, std::abs(static_cast<double>(p[i])));
        out.pbg_max_abs = std::max(out.pbg_max_abs, std::abs(static_cast<double>(pbg[i])));
    }
    out.resid_rms = std::sqrt(r2 / static_cast<double>(n));
    out.eos_rms0 = std::sqrt(e2 / static_cast<double>(n));
    auto fill = [&](const std::vector<double> &v, PiStats &s) {
        s.n = n;
        double sum = 0.0;
        s.min = std::numeric_limits<double>::max();
        s.max = -std::numeric_limits<double>::max();
        for (double x : v)
        {
            sum += x;
            s.min = std::min(s.min, x);
            s.max = std::max(s.max, x);
        }
        s.mean = sum / static_cast<double>(n);
        double var = 0.0;
        for (double x : v)
            var += (x - s.mean) * (x - s.mean);
        s.rms_demean = std::sqrt(var / static_cast<double>(n));
        s.valid = std::isfinite(s.mean) && std::isfinite(s.rms_demean);
    };
    fill(pi_p, out.pi_from_pressure);
    fill(pi_e, out.pi_from_density);
    // Pressure is stored as float(EOS + p_h): the rounding of the total (|p| up to p_total_max) bounds the residual.
    const double eps = static_cast<double>(std::numeric_limits<Real>::epsilon());
    out.tol = 8.0 * eps * out.p_total_max_abs + 1.0e-3 * out.eos_rms0;
    out.consistent = std::isfinite(out.resid_rms) && out.resid_rms <= out.tol;
    out.valid = true;
}

//----------------------------------------------------------------------------------------------------
// Doc131 B: power contract
//----------------------------------------------------------------------------------------------------
struct PowerReport
{
    std::string mode;
    bool ok = false;
    double P_target = 0, P_applied = 0, V_geom = 0, sum_V_Q = 0, mass_sum = 0, ratio_sumVQ_over_Vgeom = 0;
    double N_Q = 0, norm_coeff = 0, q_min = 0, q_max = 0, gap_rel = 0;
};

/** Doc189: one EM (or pre-EM) power stage. Coefficients are not interchangeable. */
struct EmPowerLedgerRow
{
    std::string stage;
    double t_s = 0;
    double P_em_raw_W = std::numeric_limits<double>::quiet_NaN();
    double P_particles_W = std::numeric_limits<double>::quiet_NaN();
    double P_table_raw_W = std::numeric_limits<double>::quiet_NaN();
    /** P_target / P_table_raw on the Eulerian rings. This is the raw-EM to 50 kW step. */
    double scale_table_raw_to_target = std::numeric_limits<double>::quiet_NaN();
    /** N_Q = sum s_i V_i^Q after the delivered shape is sampled onto SPH volumes. */
    double P_fluid_mapped_raw_W = 0;
    /** P_target / N_Q. SPH sampling correction only; not the 6 kW to 50 kW factor. */
    double coeff_sph_to_applied = 0;
    double P_fluid_applied_W = 0;
    double P_target_W = 0;
};

/**
 * Fill Q [W/m^3] (kDoc131QField), shape s_i:
 *   uniform     : s_i = 1
 *   radial_skin : s_i = (r/R)^2  (effective_source_study; NOT EM-verified)
 *   em_rz_csv   : s_i = max(Q_interp(r,z), 0) from a prior EM r-z dump (fixed Eulerian; discrete-renorm to P)
 *   discrete  : V_i^Q = m_i/rho0; N_Q = sum s_i V_i^Q; q_i = P s_i / N_Q   => sum q_i V_i^Q = P exactly.
 *   volumetric: q = P / V_geom (diagnostic only; P_applied differs from P by V_geom / sum V^Q).
 * Refuses N_Q <= 0 / non-finite.
 */
inline double doc133QShapeValue(const std::string &shape, double r_over_R)
{
    if (shape == "radial_skin")
        return std::max(0.0, r_over_R * r_over_R); // s = (r/R)^2 >= 0 (skin-heavy; effective_source_study)
    return 1.0;
}

/** Prior EM Q(r,z) table (stirring/natural dump CSV with columns ir,iz,r_m,z_m,Q_Wpm3[,valid]). */
struct Doc134RzQGrid
{
    std::vector<double> r_axis, z_axis;
    std::vector<double> q; /**< q[iz * nr + ir] */
    size_t nr = 0, nz = 0;
    std::string path;
    double q_table_min = 0, q_table_max = 0;
    bool ok = false;

    bool load(const std::string &csv_path)
    {
        path = csv_path;
        ok = false;
        r_axis.clear();
        z_axis.clear();
        q.clear();
        std::ifstream in(csv_path);
        if (!in)
        {
            std::cerr << "[doc134] --q-from-rz-csv cannot open: " << csv_path << std::endl;
            return false;
        }
        std::string line;
        if (!std::getline(in, line))
        {
            std::cerr << "[doc134] --q-from-rz-csv empty: " << csv_path << std::endl;
            return false;
        }
        // header: require r_m,z_m,Q_Wpm3 (ir/iz optional)
        auto lower = [](std::string s) {
            for (char &c : s)
                c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            return s;
        };
        std::vector<std::string> hdr;
        {
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ','))
                hdr.push_back(lower(cell));
        }
        int i_r = -1, i_z = -1, i_q = -1, i_valid = -1;
        for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
        {
            if (hdr[i] == "r_m" || hdr[i] == "r")
                i_r = i;
            else if (hdr[i] == "z_m" || hdr[i] == "z")
                i_z = i;
            else if (hdr[i] == "q_wpm3" || hdr[i] == "q")
                i_q = i;
            else if (hdr[i] == "valid")
                i_valid = i;
        }
        if (i_r < 0 || i_z < 0 || i_q < 0)
        {
            std::cerr << "[doc134] --q-from-rz-csv needs columns r_m,z_m,Q_Wpm3: " << csv_path << std::endl;
            return false;
        }
        struct Row
        {
            double r, z, qv;
        };
        std::vector<Row> rows;
        while (std::getline(in, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);
            std::vector<std::string> cells;
            std::string cell;
            while (std::getline(ss, cell, ','))
                cells.push_back(cell);
            if (static_cast<int>(cells.size()) <= std::max(i_r, std::max(i_z, i_q)))
                continue;
            if (i_valid >= 0 && static_cast<int>(cells.size()) > i_valid)
            {
                const std::string &v = cells[static_cast<size_t>(i_valid)];
                if (v == "0" || lower(v) == "false")
                    continue;
            }
            double r = 0, z = 0, qv = 0;
            if (!parseRealStrict(cells[static_cast<size_t>(i_r)], r) ||
                !parseRealStrict(cells[static_cast<size_t>(i_z)], z) ||
                !parseRealStrict(cells[static_cast<size_t>(i_q)], qv))
            {
                std::cerr << "[doc134] --q-from-rz-csv bad numeric row: " << line << std::endl;
                return false;
            }
            if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(qv) || r < 0)
            {
                std::cerr << "[doc134] --q-from-rz-csv non-finite/negative-r row refused" << std::endl;
                return false;
            }
            rows.push_back({r, z, qv});
        }
        if (rows.empty())
        {
            std::cerr << "[doc134] --q-from-rz-csv no data rows: " << csv_path << std::endl;
            return false;
        }
        auto uniq_sorted = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end(),
                                [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                    v.end());
            return v;
        };
        std::vector<double> rs, zs;
        rs.reserve(rows.size());
        zs.reserve(rows.size());
        for (const auto &row : rows)
        {
            rs.push_back(row.r);
            zs.push_back(row.z);
        }
        r_axis = uniq_sorted(rs);
        z_axis = uniq_sorted(zs);
        nr = r_axis.size();
        nz = z_axis.size();
        if (nr < 2 || nz < 2)
        {
            std::cerr << "[doc134] --q-from-rz-csv need >=2 r and z nodes (got " << nr << "x" << nz << ")" << std::endl;
            return false;
        }
        q.assign(nr * nz, std::numeric_limits<double>::quiet_NaN());
        auto idx_nearest = [](const std::vector<double> &ax, double v) -> size_t {
            size_t best = 0;
            double bd = std::abs(ax[0] - v);
            for (size_t i = 1; i < ax.size(); ++i)
            {
                const double d = std::abs(ax[i] - v);
                if (d < bd)
                {
                    bd = d;
                    best = i;
                }
            }
            return best;
        };
        size_t filled = 0;
        q_table_min = std::numeric_limits<double>::max();
        q_table_max = -std::numeric_limits<double>::max();
        for (const auto &row : rows)
        {
            const size_t ir = idx_nearest(r_axis, row.r);
            const size_t iz = idx_nearest(z_axis, row.z);
            const size_t k = iz * nr + ir;
            if (!std::isnan(q[k]))
                continue; // first wins on duplicates
            q[k] = row.qv;
            ++filled;
            q_table_min = std::min(q_table_min, row.qv);
            q_table_max = std::max(q_table_max, row.qv);
        }
        size_t missing = 0;
        for (double v : q)
            if (std::isnan(v))
                ++missing;
        if (missing > 0)
        {
            std::cerr << "[doc134] --q-from-rz-csv incomplete grid: filled=" << filled << " missing=" << missing
                      << " expected=" << (nr * nz) << std::endl;
            return false;
        }
        ok = true;
        std::cout << "[doc134][q-rz] loaded " << csv_path << " nr=" << nr << " nz=" << nz << " r=[" << r_axis.front()
                  << "," << r_axis.back() << "] z=[" << z_axis.front() << "," << z_axis.back() << "] Q_table=["
                  << q_table_min << "," << q_table_max << "]" << std::endl;
        return true;
    }

    double sample(double r, double z) const
    {
        if (!ok || nr < 2 || nz < 2)
            return 0.0;
        r = std::max(r, r_axis.front());
        r = std::min(r, r_axis.back());
        z = std::max(z, z_axis.front());
        z = std::min(z, z_axis.back());
        auto bracket = [](const std::vector<double> &ax, double v, size_t &i0, size_t &i1, double &t) {
            if (v <= ax.front())
            {
                i0 = 0;
                i1 = 1;
                t = 0.0;
                return;
            }
            if (v >= ax.back())
            {
                i0 = ax.size() - 2;
                i1 = ax.size() - 1;
                t = 1.0;
                return;
            }
            auto it = std::upper_bound(ax.begin(), ax.end(), v);
            i1 = static_cast<size_t>(it - ax.begin());
            i0 = i1 - 1;
            t = (v - ax[i0]) / (ax[i1] - ax[i0] + 1e-300);
        };
        size_t ir0 = 0, ir1 = 0, iz0 = 0, iz1 = 0;
        double tr = 0, tz = 0;
        bracket(r_axis, r, ir0, ir1, tr);
        bracket(z_axis, z, iz0, iz1, tz);
        const double q00 = q[iz0 * nr + ir0];
        const double q10 = q[iz0 * nr + ir1];
        const double q01 = q[iz1 * nr + ir0];
        const double q11 = q[iz1 * nr + ir1];
        const double q0 = q00 * (1.0 - tr) + q10 * tr;
        const double q1 = q01 * (1.0 - tr) + q11 * tr;
        return q0 * (1.0 - tz) + q1 * tz;
    }
};

/**
 * Doc156: axisymmetric Shepard kernel sample of EM particle Q onto fluid (r,z).
 * Shared for initial Q and EM-update reloads (no intermediate coarse rz table).
 */
struct Doc156EmParticleQCloud
{
    struct Node
    {
        double r = 0, z = 0, q = 0, v = 0;
    };
    std::vector<Node> nodes;
    std::unordered_map<long, std::vector<size_t>> bins;
    double h = 0.03;
    double cell = 0.03;
    std::string path;
    double q_min = 0, q_max = 0, P_particles = 0;
    size_t n_nodes = 0;
    bool ok = false;

    static long binKey(int ir, int iz) { return (static_cast<long>(ir) << 32) ^ static_cast<unsigned>(iz); }

    void rebuildBins()
    {
        bins.clear();
        cell = std::max(h, 1e-6);
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            const int ir = static_cast<int>(std::floor(nodes[i].r / cell));
            const int iz = static_cast<int>(std::floor(nodes[i].z / cell));
            bins[binKey(ir, iz)].push_back(i);
        }
    }

    bool load(const std::string &csv_path, double support_radius)
    {
        path = csv_path;
        ok = false;
        nodes.clear();
        bins.clear();
        h = support_radius;
        std::ifstream in(csv_path);
        if (!in)
        {
            std::cerr << "[doc156] em_particles_csv cannot open: " << csv_path << std::endl;
            return false;
        }
        std::string line;
        if (!std::getline(in, line))
        {
            std::cerr << "[doc156] em_particles_csv empty: " << csv_path << std::endl;
            return false;
        }
        auto lower = [](std::string s) {
            for (char &c : s)
                c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            return s;
        };
        std::vector<std::string> hdr;
        {
            std::stringstream ss(line);
            std::string cell_s;
            while (std::getline(ss, cell_s, ','))
                hdr.push_back(lower(cell_s));
        }
        int i_r = -1, i_z = -1, i_q = -1, i_v = -1;
        for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
        {
            if (hdr[static_cast<size_t>(i)] == "r_m" || hdr[static_cast<size_t>(i)] == "r")
                i_r = i;
            else if (hdr[static_cast<size_t>(i)] == "z_m" || hdr[static_cast<size_t>(i)] == "z")
                i_z = i;
            else if (hdr[static_cast<size_t>(i)] == "q_wpm3" || hdr[static_cast<size_t>(i)] == "q")
                i_q = i;
            else if (hdr[static_cast<size_t>(i)] == "vol_m3" || hdr[static_cast<size_t>(i)] == "vol")
                i_v = i;
        }
        if (i_r < 0 || i_z < 0 || i_q < 0)
        {
            std::cerr << "[doc156] em_particles_csv needs r_m,z_m,Q_Wpm3: " << csv_path << std::endl;
            return false;
        }
        if (i_v < 0)
        {
            std::cerr << "[doc164] em_particles_csv needs vol_m3 (formal path; no default volume)" << std::endl;
            return false;
        }
        q_min = std::numeric_limits<double>::max();
        q_max = -std::numeric_limits<double>::max();
        P_particles = 0.0;
        while (std::getline(in, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);
            std::vector<std::string> cells;
            std::string cell_s;
            while (std::getline(ss, cell_s, ','))
                cells.push_back(cell_s);
            const int need = std::max(i_r, std::max(i_z, i_q));
            if (static_cast<int>(cells.size()) <= need)
                continue;
            double r = 0, z = 0, qv = 0, vv = 1.0;
            if (!parseRealStrict(cells[static_cast<size_t>(i_r)], r) ||
                !parseRealStrict(cells[static_cast<size_t>(i_z)], z) ||
                !parseRealStrict(cells[static_cast<size_t>(i_q)], qv))
            {
                std::cerr << "[doc156] em_particles_csv bad numeric row" << std::endl;
                return false;
            }
            if (i_v < 0)
            {
                std::cerr << "[doc164] em_particles_csv formal path requires vol_m3 column (no default vol=1)"
                          << std::endl;
                return false;
            }
            if (static_cast<int>(cells.size()) <= i_v ||
                !parseRealStrict(cells[static_cast<size_t>(i_v)], vv) || !(vv > 0) || !std::isfinite(vv))
            {
                std::cerr << "[doc164] em_particles_csv refused missing/non-positive vol_m3 (no silent vol=1)"
                          << std::endl;
                return false;
            }
            if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(qv) || qv < 0 || r < 0)
            {
                std::cerr << "[doc156] em_particles_csv refused non-finite/negative Q or r" << std::endl;
                return false;
            }
            nodes.push_back({r, z, qv, vv});
            q_min = std::min(q_min, qv);
            q_max = std::max(q_max, qv);
            P_particles += qv * vv;
        }
        if (nodes.empty())
        {
            std::cerr << "[doc156] em_particles_csv no rows: " << csv_path << std::endl;
            return false;
        }
        n_nodes = nodes.size();
        rebuildBins();
        ok = true;
        std::cout << "[doc156][q-em-particles] loaded " << csv_path << " n=" << n_nodes << " h=" << h
                  << " Q=[" << q_min << "," << q_max << "] P_particles=" << P_particles << std::endl;
        return true;
    }

    double sample(double r, double z) const
    {
        if (!ok || nodes.empty())
            return 0.0;
        const double h2 = h * h;
        const int ir0 = static_cast<int>(std::floor(r / cell));
        const int iz0 = static_cast<int>(std::floor(z / cell));
        double num = 0.0, den = 0.0;
        for (int di = -1; di <= 1; ++di)
        {
            for (int dj = -1; dj <= 1; ++dj)
            {
                const auto it = bins.find(binKey(ir0 + di, iz0 + dj));
                if (it == bins.end())
                    continue;
                for (size_t idx : it->second)
                {
                    const Node &n = nodes[idx];
                    const double d2 = (r - n.r) * (r - n.r) + (z - n.z) * (z - n.z);
                    if (d2 > h2)
                        continue;
                    const double s = std::sqrt(d2) / (h + 1e-300);
                    const double w = (1.0 - s) * (1.0 - s) * n.v;
                    num += w * n.q;
                    den += w;
                }
            }
        }
        if (den > 0.0)
            return num / den;
        // nearest EM particle
        double bd = 1e300, bq = 0.0;
        for (const Node &n : nodes)
        {
            const double d2 = (r - n.r) * (r - n.r) + (z - n.z) * (z - n.z);
            if (d2 < bd)
            {
                bd = d2;
                bq = n.q;
            }
        }
        return bq;
    }
};

inline bool hostApplyPowerContract(BaseParticles &particles, const std::string &power_mode, double P_target,
                                   double rho0, double V_geom, PowerReport &rep,
                                   const std::string &q_shape = "uniform", const Vecd &axis_center = Vecd::Zero(),
                                   double R = 1.0, const Doc134RzQGrid *rz_grid = nullptr,
                                   const Doc156EmParticleQCloud *em_q = nullptr)
{
    rep = PowerReport();
    rep.mode = power_mode;
    rep.P_target = P_target;
    rep.V_geom = V_geom;
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kDoc131QField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *q = particles.getVariableDataByName<Real>(kDoc131QField);
    const size_t n = particles.TotalRealParticles();
    std::vector<double> sh(n, 1.0);
    double sum_vq = 0.0, mass_sum = 0.0, n_q = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - axis_center[0]);
        const double dy = static_cast<double>(pos[i][1] - axis_center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        if (q_shape == "em_rz_csv")
        {
            if (rz_grid == nullptr || !rz_grid->ok)
            {
                std::cerr << "[doc134] em_rz_csv shape requires a loaded RzQGrid" << std::endl;
                return false;
            }
            sh[i] = std::max(0.0, rz_grid->sample(r, static_cast<double>(pos[i][2])));
        }
        else if (q_shape == "em_particles_csv")
        {
            if (em_q == nullptr || !em_q->ok)
            {
                std::cerr << "[doc156] em_particles_csv shape requires a loaded EmParticleQCloud" << std::endl;
                return false;
            }
            sh[i] = std::max(0.0, em_q->sample(r, static_cast<double>(pos[i][2])));
        }
        else
            sh[i] = doc133QShapeValue(q_shape, r / (R + 1e-300));
        const double vq = static_cast<double>(mass[i]) / rho0;
        sum_vq += vq;
        mass_sum += static_cast<double>(mass[i]);
        n_q += sh[i] * vq; // N_Q = sum s_i V_i^Q
    }
    rep.sum_V_Q = sum_vq;
    rep.mass_sum = mass_sum;
    rep.N_Q = n_q;
    rep.ratio_sumVQ_over_Vgeom = sum_vq / (V_geom + 1e-300);
    if (!(n_q > 0.0) || !std::isfinite(n_q))
    {
        std::cerr << "[doc131] power contract refused: N_Q=" << n_q << " <= 0 or non-finite" << std::endl;
        return false;
    }
    const bool discrete = power_mode == kPowerDiscrete;
    const double q_val = discrete ? P_target / n_q : P_target / V_geom;
    rep.norm_coeff = P_target / n_q;
    rep.q_min = std::numeric_limits<double>::max();
    rep.q_max = -std::numeric_limits<double>::max();
    double p_applied = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        q[i] = static_cast<Real>(q_val * sh[i]);
        const double vq = static_cast<double>(mass[i]) / rho0;
        p_applied += static_cast<double>(q[i]) * vq;
        rep.q_min = std::min(rep.q_min, static_cast<double>(q[i]));
        rep.q_max = std::max(rep.q_max, static_cast<double>(q[i]));
    }
    rep.P_applied = p_applied;
    rep.gap_rel = P_target > 0.0 ? (p_applied - P_target) / P_target : 0.0;
    rep.ok = std::isfinite(p_applied);
    syncVariableToDevice<Real>(particles, kDoc131QField);
    return rep.ok;
}

//----------------------------------------------------------------------------------------------------
// Doc131 B3: face area audit + weights for the production Robin/radiation kernel
//----------------------------------------------------------------------------------------------------
struct FaceAreaAudit
{
    size_t n = 0;
    double raw = 0, target = 0, scale = 1, after = 0, rel_err_raw = 0, rel_err = 0, a_i_min = 0, a_i_max = 0;
};

struct AreaAudit
{
    FaceAreaAudit top, side, bottom;
    size_t n_corner = 0, n_none = 0;
    double shell = 0;
    bool normalized = false;
};

/**
 * Faces: top z > z_top - shell, bottom z < z_bottom + shell, side r > R - shell. A_i = m_i/(rho0 shell).
 * Particles in several shells (corners) are assigned to ONE primary face (closest surface), so no area
 * is double counted. Writes exclusive weights into the thermal boundary weight fields consumed by
 * ApplyOphelieThermalFrenchNaturalBcCK / WallTemperatureBcCK (kernel area = VolumetricMeasure * w / shell,
 * weight w = scale_f when --area-normalize else 1; faces not enabled get weight 0).
 * Doc207 B: also writes OphelieThermalWallGap = geometric distance to the assigned primary face.
 */
inline AreaAudit hostBuildBoundaryAreaWeights(BaseParticles &particles, const Vecd &center, Real R, Real H,
                                              Real z_bottom, Real rho0, Real shell, bool normalize,
                                              bool enable_top, bool enable_side, bool enable_bottom)
{
    AreaAudit a;
    a.shell = static_cast<double>(shell);
    a.normalized = normalize;
    const Real z_top = z_bottom + H;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieThermalWallGapField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *mask = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryMaskField);
    Real *face = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryFaceField);
    Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    Real *gap = particles.getVariableDataByName<Real>(kOphelieThermalWallGapField);
    const size_t n = particles.TotalRealParticles();
    std::vector<int> face_id(n, 0); // 0 none, 1 side, 2 bottom, 3 top
    std::vector<double> a_i(n, 0.0);
    std::vector<double> d_face(n, 0.0);
    FaceAreaAudit *faces[4] = {nullptr, &a.side, &a.bottom, &a.top};
    for (int f = 1; f <= 3; ++f)
    {
        faces[f]->a_i_min = std::numeric_limits<double>::max();
        faces[f]->a_i_max = 0.0;
    }
    for (size_t i = 0; i < n; ++i)
    {
        gap[i] = Real(0);
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double d_side = static_cast<double>(R) - r;
        const double d_bottom = static_cast<double>(pos[i][2] - z_bottom);
        const double d_top = static_cast<double>(z_top - pos[i][2]);
        const double sh = static_cast<double>(shell);
        const bool in_side = d_side < sh, in_bottom = d_bottom < sh, in_top = d_top < sh;
        const int n_in = (in_side ? 1 : 0) + (in_bottom ? 1 : 0) + (in_top ? 1 : 0);
        if (n_in == 0)
        {
            ++a.n_none;
            continue;
        }
        if (n_in >= 2)
            ++a.n_corner;
        double best = std::numeric_limits<double>::max();
        int fid = 0;
        if (in_side && d_side < best)
        {
            best = d_side;
            fid = 1;
        }
        if (in_bottom && d_bottom < best)
        {
            best = d_bottom;
            fid = 2;
        }
        if (in_top && d_top < best)
        {
            best = d_top;
            fid = 3;
        }
        face_id[i] = fid;
        d_face[i] = std::max(0.0, best);
        a_i[i] = static_cast<double>(mass[i]) / (static_cast<double>(rho0) * sh);
        faces[fid]->n += 1;
        faces[fid]->raw += a_i[i];
        faces[fid]->a_i_min = std::min(faces[fid]->a_i_min, a_i[i]);
        faces[fid]->a_i_max = std::max(faces[fid]->a_i_max, a_i[i]);
    }
    const double Rd = static_cast<double>(R), Hd = static_cast<double>(H);
    a.top.target = kPiD * Rd * Rd;
    a.bottom.target = kPiD * Rd * Rd;
    a.side.target = 2.0 * kPiD * Rd * Hd;
    for (int f = 1; f <= 3; ++f)
    {
        FaceAreaAudit &fa = *faces[f];
        if (fa.n == 0)
            fa.a_i_min = 0.0;
        fa.rel_err_raw = (fa.raw - fa.target) / fa.target;
        fa.scale = (normalize && fa.raw > 1.0e-30) ? fa.target / fa.raw : 1.0;
        fa.after = fa.raw * fa.scale;
        fa.rel_err = (fa.after - fa.target) / fa.target;
    }
    for (size_t i = 0; i < n; ++i)
    {
        w_side[i] = w_bottom[i] = w_free[i] = Real(0);
        mask[i] = face_id[i] == 0 ? Real(0) : Real(1);
        face[i] = static_cast<Real>(face_id[i]);
        if (face_id[i] == 1 && enable_side)
            w_side[i] = static_cast<Real>(a.side.scale);
        else if (face_id[i] == 2 && enable_bottom)
            w_bottom[i] = static_cast<Real>(a.bottom.scale);
        else if (face_id[i] == 3 && enable_top)
            w_free[i] = static_cast<Real>(a.top.scale);
        if (face_id[i] != 0 && (w_side[i] + w_bottom[i] + w_free[i]) > TinyReal)
            gap[i] = static_cast<Real>(d_face[i]);
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryFaceField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToDevice<Real>(particles, kOphelieThermalWallGapField);
    return a;
}

//----------------------------------------------------------------------------------------------------
// Doc131 C: thermal kernels / host helpers
//----------------------------------------------------------------------------------------------------
/** Delta-authoritative heating with dt at exec time: dT += q dt/(rho cp); T = t_initial + dT. */
class ApplyDoc131SourceQCK : public LocalDynamics
{
  public:
    ApplyDoc131SourceQCK(SPHBody &sph_body, const std::string &q_field, Real rho, Real cp, Real t_initial)
        : LocalDynamics(sph_body), inv_rho_cp_(Real(1) / (rho * cp + TinyReal)), t_initial_(t_initial),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
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
            delta_t_[index_i] += q_[index_i] * dt * inv_rho_cp_;
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

struct ThermalLedger
{
    double e_source = 0.0;   /**< J, sum P_applied dt */
    double e_boundary = 0.0; /**< J, sum P_loss dt (positive = lost) */
    double t_thermal = 0.0;  /**< s, sum of thermal dt actually applied */
    size_t n_steps = 0;
    double dt_min = std::numeric_limits<double>::max();
    double dt_max = 0.0;
    double p_loss_first = 0.0;
    bool p_loss_first_set = false;
};

/** H = sum m cp DeltaT (J); DeltaT is the authoritative thermal variable (Temperature = Tref + DeltaT). */
inline double hostEnthalpyDeltaT(BaseParticles &particles, Real cp)
{
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
        e += static_cast<double>(mass[i]) * static_cast<double>(cp) * static_cast<double>(dT[i]);
    return e;
}

/** sum m cp |DeltaT| (energy scale for relative checks). */
inline double hostEnthalpyAbsScale(BaseParticles &particles, Real cp)
{
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
        e += static_cast<double>(mass[i]) * static_cast<double>(cp) * std::abs(static_cast<double>(dT[i]));
    return e;
}

struct DeltaTStats
{
    double mean_m = 0, var_m = 0, min = 0, max = 0, abs_max = 0;
    bool finite = true;
};

inline DeltaTStats hostDeltaTStats(BaseParticles &particles)
{
    DeltaTStats s;
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double msum = 0.0, sum = 0.0;
    s.min = std::numeric_limits<double>::max();
    s.max = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < n; ++i)
    {
        const double d = static_cast<double>(dT[i]);
        if (!std::isfinite(d))
            s.finite = false;
        msum += static_cast<double>(mass[i]);
        sum += static_cast<double>(mass[i]) * d;
        s.min = std::min(s.min, d);
        s.max = std::max(s.max, d);
        s.abs_max = std::max(s.abs_max, std::abs(d));
    }
    s.mean_m = sum / (msum + 1e-300);
    double v = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double d = static_cast<double>(dT[i]) - s.mean_m;
        v += static_cast<double>(mass[i]) * d * d;
    }
    s.var_m = v / (msum + 1e-300);
    return s;
}

/**
 * T4: least-squares fit DeltaT_i = c + A cos(pi zeta_i / H) (volume weights); returns c, A and the rms residual.
 */
inline void hostFitCosineProfile(BaseParticles &particles, Real H, Real z_bottom, double &c, double &A,
                                 double &rms_resid)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double s_w = 0, s_x = 0, s_xx = 0, s_y = 0, s_xy = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double zeta = std::min(static_cast<double>(H), std::max(0.0, static_cast<double>(pos[i][2] - z_bottom)));
        const double x = std::cos(kPiD * zeta / static_cast<double>(H));
        const double w = static_cast<double>(vol[i]);
        const double y = static_cast<double>(dT[i]);
        s_w += w;
        s_x += w * x;
        s_xx += w * x * x;
        s_y += w * y;
        s_xy += w * x * y;
    }
    const double det = s_w * s_xx - s_x * s_x;
    A = (s_w * s_xy - s_x * s_y) / (det + 1e-300);
    c = (s_y - A * s_x) / (s_w + 1e-300);
    double r2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double zeta = std::min(static_cast<double>(H), std::max(0.0, static_cast<double>(pos[i][2] - z_bottom)));
        const double x = std::cos(kPiD * zeta / static_cast<double>(H));
        const double r = static_cast<double>(dT[i]) - (c + A * x);
        r2 += static_cast<double>(vol[i]) * r * r;
    }
    rms_resid = std::sqrt(r2 / (s_w + 1e-300));
}

/**
 * Doc144 Marangoni on the flat free top: axisymmetric surface gradient approximation.
 * ∇_s T ≈ (dT/dr) e_r. dT/dr from radial-bin mean T(r) + central differences (not global r-sort
 * nearest-neighbor, which mixed azimuthally distant particles).
 * τ = (dγ/dT)(dT/dr); ForcePrior += τ A e_r with A = w_free * Vol / shell.
 */
struct MarangoniApplyResult
{
    std::string method = "axisymmetric_radial_bin_FD";
    size_t n_free = 0;
    size_t n_bins_used = 0;
    Real sum_area = 0;
    Real sum_Fr = 0; /**< Σ F·e_r [N] */
    Real mean_dT_dr = 0;
    Real mean_tau = 0;
    Real n_outward = 0;
    Real n_inward = 0;
    bool applied = false;
};

inline MarangoniApplyResult hostApplyMarangoniAxisymmetric(BaseParticles &particles, const Vecd &axis_center, Real R,
                                                           Real shell, Real dgamma_dT, bool add_to_force_prior,
                                                           size_t n_bins = 24)
{
    MarangoniApplyResult out;
    if (!(std::abs(dgamma_dT) > TinyReal) || !(R > TinyReal) || !(shell > TinyReal) || n_bins < 3)
        return out;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    if (add_to_force_prior)
        syncVariableToHost<Vecd>(particles, "ForcePrior");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    Vecd *fp = add_to_force_prior ? particles.getVariableDataByName<Vecd>("ForcePrior") : nullptr;
    const size_t n = particles.TotalRealParticles();
    const Real inv_shell = Real(1) / shell;

    std::vector<Real> bin_Tsum(n_bins, 0), bin_W(n_bins, 0);
    std::vector<size_t> free_ids;
    free_ids.reserve(n / 8 + 8);
    for (size_t i = 0; i < n; ++i)
    {
        if (!(w_free[i] > TinyReal))
            continue;
        free_ids.push_back(i);
        const Real dx = pos[i][0] - axis_center[0];
        const Real dy = pos[i][1] - axis_center[1];
        const Real r = std::hypot(dx, dy);
        const size_t b = std::min(n_bins - 1, static_cast<size_t>(std::max(Real(0), r) / R * static_cast<Real>(n_bins)));
        bin_Tsum[b] += temp[i];
        bin_W[b] += Real(1);
    }
    out.n_free = free_ids.size();
    if (free_ids.size() < 4)
        return out;

    std::vector<Real> bin_T(n_bins, 0), bin_dTdr(n_bins, 0);
    std::vector<char> bin_ok(n_bins, 0);
    for (size_t b = 0; b < n_bins; ++b)
    {
        if (bin_W[b] > TinyReal)
        {
            bin_T[b] = bin_Tsum[b] / bin_W[b];
            bin_ok[b] = 1;
            ++out.n_bins_used;
        }
    }
    // Fill gaps by nearest valid bin so FD is defined across the free surface.
    for (size_t b = 0; b < n_bins; ++b)
    {
        if (bin_ok[b])
            continue;
        size_t bl = b, br = b;
        while (bl > 0 && !bin_ok[bl])
            --bl;
        while (br + 1 < n_bins && !bin_ok[br])
            ++br;
        if (bin_ok[bl] && bin_ok[br] && bl != br)
            bin_T[b] = Real(0.5) * (bin_T[bl] + bin_T[br]);
        else if (bin_ok[bl])
            bin_T[b] = bin_T[bl];
        else if (bin_ok[br])
            bin_T[b] = bin_T[br];
        else
            return out;
        bin_ok[b] = 1;
    }
    const Real dr_bin = R / static_cast<Real>(n_bins);
    for (size_t b = 0; b < n_bins; ++b)
    {
        if (b == 0)
            bin_dTdr[b] = (bin_T[1] - bin_T[0]) / (dr_bin + TinyReal);
        else if (b + 1 == n_bins)
            bin_dTdr[b] = (bin_T[b] - bin_T[b - 1]) / (dr_bin + TinyReal);
        else
            bin_dTdr[b] = (bin_T[b + 1] - bin_T[b - 1]) / (Real(2) * dr_bin + TinyReal);
    }

    Real sum_dT = 0, sum_tau = 0;
    for (size_t i : free_ids)
    {
        const Real dx = pos[i][0] - axis_center[0];
        const Real dy = pos[i][1] - axis_center[1];
        const Real r = std::hypot(dx, dy);
        if (r <= TinyReal)
            continue;
        const size_t b =
            std::min(n_bins - 1, static_cast<size_t>(std::max(Real(0), r) / R * static_cast<Real>(n_bins)));
        const Real dT_dr = bin_dTdr[b];
        const Real tau = dgamma_dT * dT_dr;
        const Real area = w_free[i] * vol[i] * inv_shell;
        const Real f_mag = tau * area;
        out.sum_area += area;
        out.sum_Fr += f_mag;
        sum_dT += dT_dr;
        sum_tau += tau;
        if (f_mag > TinyReal)
            out.n_outward += Real(1);
        else if (f_mag < -TinyReal)
            out.n_inward += Real(1);
        if (fp)
        {
            fp[i][0] += f_mag * (dx / r);
            fp[i][1] += f_mag * (dy / r);
        }
    }
    const Real n_eff = static_cast<Real>(std::max<size_t>(1, free_ids.size()));
    out.mean_dT_dr = sum_dT / n_eff;
    out.mean_tau = sum_tau / n_eff;
    out.applied = true;
    if (add_to_force_prior)
        syncVariableToDevice<Vecd>(particles, "ForcePrior");
    return out;
}

/** Prescribe T=T0+dT_amp*(1-r/R) (hot center), run Marangoni once into scratch ForcePrior delta, check direction+∫. */
struct MarangoniVerifyResult
{
    bool passed = false;
    bool direction_ok = false;
    bool area_integral_ok = false;
    Real dgamma_dT = 0;
    Real dT_dr_theory = 0;
    Real sum_area = 0;
    Real sum_Fr = 0;
    Real sum_Fr_theory = 0;
    Real rel_err_integral = 0;
    Real n_outward = 0;
    Real n_inward = 0;
    std::string method;
    std::string message;
};

inline MarangoniVerifyResult hostVerifyMarangoniAxisymmetric(BaseParticles &particles, const Vecd &axis_center, Real R,
                                                             Real shell, Real dgamma_dT, Real t_ref,
                                                             Real T0 = Real(1600), Real dT_amp = Real(200))
{
    MarangoniVerifyResult v;
    v.dgamma_dT = dgamma_dT;
    v.dT_dr_theory = -dT_amp / (R + TinyReal); // T = T0 + dT_amp*(1-r/R) => dT/dr = -dT_amp/R
    v.method = "axisymmetric_radial_bin_FD";
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const size_t n = particles.TotalRealParticles();
    std::vector<Real> T_bak(n), dT_bak(n);
    std::vector<Vecd> fp_bak(n);
    for (size_t i = 0; i < n; ++i)
    {
        T_bak[i] = temp[i];
        dT_bak[i] = delta_t[i];
        fp_bak[i] = fp[i];
        const Real r = std::hypot(pos[i][0] - axis_center[0], pos[i][1] - axis_center[1]);
        const Real Ti = T0 + dT_amp * (Real(1) - std::min(r, R) / (R + TinyReal));
        temp[i] = Ti;
        delta_t[i] = Ti - t_ref;
        fp[i] = Vecd::Zero(); // isolate Marangoni contribution
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToDevice<Vecd>(particles, "ForcePrior");

    const MarangoniApplyResult app =
        hostApplyMarangoniAxisymmetric(particles, axis_center, R, shell, dgamma_dT, /*add=*/true);
    v.sum_area = app.sum_area;
    v.sum_Fr = app.sum_Fr;
    v.n_outward = app.n_outward;
    v.n_inward = app.n_inward;
    v.sum_Fr_theory = dgamma_dT * v.dT_dr_theory * app.sum_area; // constant dT/dr
    v.rel_err_integral =
        std::abs(v.sum_Fr - v.sum_Fr_theory) / (std::abs(v.sum_Fr_theory) + TinyReal);
    // dγ/dT<0, dT/dr<0 => τ>0 => outward
    v.direction_ok = (dgamma_dT < 0) && (v.sum_Fr > 0) && (v.n_outward + TinyReal >= Real(0.9) * (v.n_outward + v.n_inward));
    v.area_integral_ok = v.rel_err_integral < Real(0.15); // bin FD + discrete area tolerance
    v.passed = app.applied && v.direction_ok && v.area_integral_ok;
    v.message = v.passed ? "marangoni_axisymmetric_verify_pass"
                         : "marangoni_axisymmetric_verify_fail";

    for (size_t i = 0; i < n; ++i)
    {
        temp[i] = T_bak[i];
        delta_t[i] = dT_bak[i];
        fp[i] = fp_bak[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    return v;
}

/**
 * LEGACY (--legacy-rad-clamp, baseline only): old host top-layer radiation with the one-way
 * max(T, T_amb) clamp and layer measure rho cp * layer. Not the Doc131 formal path. Returns J lost.
 */
inline double hostLegacyTopRadiationClampStep(BaseParticles &particles, Real rho0, Real cp, Real dt, Real emissivity,
                                              Real t_ambient, Real z_top, Real layer, Real t_ref)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    const double rho_cp = static_cast<double>(rho0) * static_cast<double>(cp);
    double e_lost = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        if (!(pos[i][2] > z_top - layer))
            continue;
        const double Ti = std::max(static_cast<double>(t_ref) + static_cast<double>(delta_t[i]),
                                   static_cast<double>(t_ambient));
        const double flux = static_cast<double>(emissivity) * kStefanBoltzmann *
                            (std::pow(Ti, 4.0) - std::pow(static_cast<double>(t_ambient), 4.0));
        const double d_rad = -flux * static_cast<double>(dt) / (rho_cp * static_cast<double>(layer));
        delta_t[i] += static_cast<Real>(d_rad);
        temp[i] = t_ref + delta_t[i];
        e_lost += -static_cast<double>(mass[i]) * static_cast<double>(cp) * d_rad;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    return e_lost;
}

//----------------------------------------------------------------------------------------------------
// JSON / misc output helpers
//----------------------------------------------------------------------------------------------------
inline std::string jnum(double v)
{
    if (!std::isfinite(v))
        return "null";
    std::ostringstream os;
    os.precision(10);
    os << v;
    return os.str();
}

inline std::string jesc(const std::string &s)
{
    std::string o;
    for (char ch : s)
    {
        switch (ch)
        {
        case '"':
            o += "\\\"";
            break;
        case '\\':
            o += "\\\\";
            break;
        case '\n':
            o += "\\n";
            break;
        case '\t':
            o += "\\t";
            break;
        default:
            o += ch;
        }
    }
    return o;
}

enum class Gate
{
    Null,
    No,
    Yes
};

struct JsonWriter
{
    std::vector<std::string> items;
    void raw(const std::string &k, const std::string &v) { items.push_back("  \"" + k + "\": " + v); }
    void num(const std::string &k, double v) { raw(k, jnum(v)); }
    void integer(const std::string &k, long long v) { raw(k, std::to_string(v)); }
    void flag(const std::string &k, bool v) { raw(k, v ? "1" : "0"); }
    void str(const std::string &k, const std::string &v) { raw(k, "\"" + jesc(v) + "\""); }
    void gate(const std::string &k, Gate g, const std::string &reason)
    {
        raw(k, g == Gate::Null ? "null" : (g == Gate::Yes ? "1" : "0"));
        if (g != Gate::Yes)
            str(k + "_reason", reason);
    }
    std::string dump() const
    {
        std::string o = "{\n";
        for (size_t i = 0; i < items.size(); ++i)
            o += items[i] + (i + 1 < items.size() ? ",\n" : "\n");
        o += "}\n";
        return o;
    }
};

inline void ensureParentDir(const std::string &path)
{
    if (path.empty())
        return;
    const std::filesystem::path p(path);
    if (p.has_parent_path())
    {
        std::error_code ec;
        std::filesystem::create_directories(p.parent_path(), ec);
    }
}
} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
    {
        return 1;
    }
    if (cli.thermal_wall_model != "legacy-shell" && cli.thermal_wall_model != "wall-temperature" &&
        cli.thermal_wall_model != "wall_temperature" && cli.thermal_wall_model != "series-gap")
    {
        std::cerr << "[doc131] --thermal-wall-model must be legacy-shell|wall-temperature\n";
        return 1;
    }
    std::string argv_joined;
    for (int i = 0; i < argc; ++i)
        argv_joined += (i ? " " : "") + std::string(argv[i]);

    const bool is_R0 = cli.mode == "R0";
    const bool is_R1 = cli.mode == "R1";
    const bool is_R2 = cli.mode == "R2";
    const bool is_R3 = cli.mode == "R3";
    const bool is_baseline = cli.mode == "baseline";
    const bool is_T = isThermalTestMode(cli.mode);
    const bool is_T0 = cli.mode == "T0", is_T1 = cli.mode == "T1", is_T2 = cli.mode == "T2",
               is_T3 = cli.mode == "T3", is_T4 = cli.mode == "T4", is_T5 = cli.mode == "T5";
    const bool freeze_flow = cli.freeze_flow || is_T0 || is_T1 || is_T2 || is_T3 || is_T4;
    if (cli.doc201_b2_live_top_shear)
    {
        if (freeze_flow)
        {
            std::cerr << "[doc201][b2] --doc201-b2-live-top-shear forbids freeze_flow / T0-T4"
                      << std::endl;
            return 1;
        }
        // Doc201: allow --doc201-b1-ablate with B2 for localization (acoustic/visc/wall).
        if (!is_R0)
        {
            std::cerr << "[doc201][b2] require --mode=R0 (uniform T; no strat driver)" << std::endl;
            return 1;
        }
    }
    const bool legacy_path = cli.legacy_rad_clamp && is_baseline;
    // Doc201: --thaw-T also unlocks conduction+Robin on R1 (strat IC + evolving T short case).
    const bool thaw_thermal = cli.thaw_T && (is_baseline || is_R1);
    const bool thermal_active = is_T || thaw_thermal;
    const bool q_active = (is_T1 || is_T5 || thaw_thermal) && cli.q_power > TinyReal;
    const bool bc_active = !legacy_path && (is_T3 || is_T5 || thaw_thermal);
    const bool conduction_active = thermal_active && !legacy_path;
    // Face enable flags (Doc131 C4): T3 enables one face; T5/baseline enable side+bottom+top (top unless --no-top-loss).
    const bool face_top_on = bc_active && (is_T3 ? cli.t3_face == "top" : cli.top_loss);
    const bool face_side_on = bc_active && (is_T3 ? cli.t3_face == "side" : true);
    const bool face_bottom_on = bc_active && (is_T3 ? cli.t3_face == "bottom" : true);
    const char *conduction_path =
        !thermal_active ? "none"
                        : (legacy_path ? "none_legacy_top_radiation"
                                       : "production_ck:OpheliePairwiseLaplaceCK+ApplyLaplaceDiffusionDtCK");
    const bool use_wall_T = (cli.thermal_wall_model == "wall-temperature" ||
                             cli.thermal_wall_model == "wall_temperature" || cli.thermal_wall_model == "series-gap");
    const char *boundary_path =
        !bc_active
            ? (legacy_path ? "legacy_host_top_radiation_clamp" : "none_adiabatic")
            : (use_wall_T ? "doc207b:ApplyOphelieThermalWallTemperatureBcCK(series-gap,B3 weights)"
                          : "production_ck:ApplyOphelieThermalFrenchNaturalBcCK(B3 area weights)");

    // Geometry (Doc126): French melt cylinder, bottom at z=0.
    OphelieFrenchReducedCaseParams french;
    french.glass_radius = cli.glass_radius;
    french.glass_half_height = Real(0.5) * cli.glass_height;
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    french.dp = cli.dp;
    french.glass_mesh_resolution = cli.mesh_resolution;

    const Real dp = cli.dp;
    const Real R = cli.glass_radius;
    const Real H = cli.glass_height;
    const Real z_bottom = 0.0;
    const Real wall_thickness = cli.wall_thickness_factor * dp;
    const Vecd axis_center = french.glass_center;
    const Real pad = Real(2) * dp;
    const double v_cyl = kPiD * static_cast<double>(R) * static_cast<double>(R) * static_cast<double>(H);

    // Effective physics per mode.
    const Real beta_eff = is_T5 ? Real(0) : cli.beta; // T5: buoyancy off
    const Real a_eff = is_R1 ? cli.a_strat : Real(0);
    const Real r2_dT_eff = cli.r2_reverse ? -cli.dT_radial : cli.dT_radial; // Doc131 D2 (R2 only)
    if ((cli.r2_reverse || cli.r2_fixed_spatial_T) && !is_R2)
        std::cout << "[doc131] note: --r2-reverse / --r2-fixed-spatial-T only affect --mode=R2 (ignored here)"
                  << std::endl;
    const Real settle_eff = freeze_flow ? Real(0) : cli.settle_time;

    // Thermal frozen-flow stability (explicit pairwise Laplace + Robin), conservative estimates.
    const double alpha_th = static_cast<double>(cli.k) / (static_cast<double>(cli.rho0) * static_cast<double>(cli.cp));
    const double dt_diff_limit = static_cast<double>(dp) * static_cast<double>(dp) / (6.0 * alpha_th);
    const double h_max = std::max({static_cast<double>(cli.h_side), static_cast<double>(cli.h_bottom),
                                   static_cast<double>(cli.h_top), 4.0 * kStefanBoltzmann *
                                       std::pow(static_cast<double>(cli.t_ref), 3.0) *
                                       static_cast<double>(cli.emissivity)});
    const double dt_bc_limit =
        static_cast<double>(cli.rho0) * static_cast<double>(cli.cp) * static_cast<double>(dp) / (h_max + 1e-300);
    const double dt_thermal_limit = std::min(dt_diff_limit, dt_bc_limit);

    // ---- Effective configuration (single source of truth for stdout and JSON) ----
    std::cout << "[doc131][config] mode=" << cli.mode << " dp=" << cli.dp << " R=" << cli.glass_radius
              << " H=" << cli.glass_height << " rho0=" << cli.rho0 << " c0=" << cli.c0 << " mu=" << cli.mu
              << " g=" << cli.gravity_g << " beta=" << cli.beta << " beta_eff=" << beta_eff << " t_ref=" << cli.t_ref
              << " u_ref=" << cli.u_ref << "\n"
              << "[doc131][config] cp=" << cli.cp << " (effective_cp) k=" << cli.k
              << " (effective_k) q_power=" << cli.q_power << " power_mode=" << cli.power_mode << " q_shape=" << cli.q_shape
              << " h_side=" << cli.h_side << " h_bottom=" << cli.h_bottom << " h_top=" << cli.h_top
              << " t_cool=" << cli.t_cool << " t_ambient=" << cli.t_ambient
              << " t_rad_ambient=" << cli.t_rad_ambient << " emissivity=" << cli.emissivity
              << " variable_properties=" << (cli.variable_properties ? 1 : 0)
              << " literature_id=" << (cli.literature_id.empty() ? "-" : cli.literature_id) << "\n"
              << "[doc131][config] settle_time=" << cli.settle_time << " (effective " << settle_eff
              << ") end_time=" << cli.end_time << " thermal_dt=" << cli.thermal_dt << " (limit~" << dt_thermal_limit
              << ") freeze_flow=" << (freeze_flow ? 1 : 0) << " thermal_active=" << (thermal_active ? 1 : 0)
              << " q_active=" << (q_active ? 1 : 0) << " bc_active=" << (bc_active ? 1 : 0)
              << " faces(top/side/bottom)=" << (face_top_on ? 1 : 0) << "/" << (face_side_on ? 1 : 0) << "/"
              << (face_bottom_on ? 1 : 0) << " t3_face=" << cli.t3_face
              << " area_normalize=" << (cli.area_normalize ? 1 : 0)
              << " thermal_wall_model=" << cli.thermal_wall_model
              << " legacy_rad_clamp=" << (cli.legacy_rad_clamp ? 1 : 0)
              << " thermal_fields_follow_sort=" << (cli.legacy_thermal_sort ? 0 : 1) << "\n"
              << "[doc131][config] conduction_path=" << conduction_path << " boundary_path=" << boundary_path
              << " real_bytes=" << sizeof(Real) << " out_csv=" << cli.out_csv << " summary_json=" << cli.summary_json
              << " output_dir=" << cli.output_dir << std::endl;
    if (cli.legacy_rad_clamp && !is_baseline)
        std::cout << "[doc131] note: --legacy-rad-clamp only affects --mode=baseline (ignored here)" << std::endl;
    if (is_T && freeze_flow && cli.settle_time > 0)
        std::cout << "[doc131] note: frozen-flow thermal test: --settle-time ignored" << std::endl;
    if (freeze_flow && (is_T0 || is_T1 || is_T2 || is_T3 || is_T4) && cli.thermal_dt > dt_thermal_limit)
    {
        std::cerr << "[doc131] --thermal-dt=" << cli.thermal_dt << " exceeds the explicit stability limit "
                  << dt_thermal_limit << " s" << std::endl;
        return 1;
    }

    // Domain bounds (cylinder bbox with pad; NOT the coil-sized frenchReducedDomainBounds).
    const Real xy_ext = R + wall_thickness + pad;
    const BoundingBoxd system_bounds(
        Vecd(axis_center[0] - xy_ext, axis_center[1] - xy_ext, z_bottom - wall_thickness - pad),
        Vecd(axis_center[0] + xy_ext, axis_center[1] + xy_ext, z_bottom + H + wall_thickness + pad));
    SPHSystem sph_system(system_bounds, dp);

    FluidBody fluid(sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                                    "GlassBody", french.glass_center, french.glass_radius,
                                    french.glass_half_height, french.glass_mesh_resolution));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(cli.rho0, cli.c0);
    const bool constant_mu_implicit =
        cli.viscosity_mode == "implicit_split" && !cli.variable_properties;
    const bool use_local_particle_mu =
        (cli.variable_properties && cli.local_particle_viscosity && !cli.mean_mu_viscous) || constant_mu_implicit;
    if (use_local_particle_mu)
        fluid.addMaterialProperty<ParticleViscosity>(cli.mu);
    else
        fluid.addMaterialProperty<Viscosity>(cli.mu);
    Viscosity &fluid_viscosity = fluid.getMaterialProperty<Viscosity>();
    ParticleViscosity *fluid_particle_viscosity =
        use_local_particle_mu ? &fluid.getMaterialProperty<ParticleViscosity>() : nullptr;
    fluid.generateParticles<BaseParticles, Lattice>();

    // Acoustic pressure contact: closed crucible incl. flat lid (impermeable).
    SolidBody wall_closed(sph_system, makeShared<OphelieFrenchClosedCrucibleWallShape>(
                                          "WallClosed", french, wall_thickness, cli.mesh_resolution));
    defineClosedCylinderSolid(wall_closed);

    // Viscous no-slip contact: open cup (side + bottom), rim=0 -> lid is free-slip (no viscous contact).
    SolidBody wall_noslip(sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                                          "WallNoSlip", french, wall_thickness, cli.mesh_resolution, Real(0)));
    defineClosedCylinderSolid(wall_noslip);

    Inner<> fluid_inner(fluid);
    Contact<> fluid_closed_contact(fluid, {&wall_closed});
    Contact<> fluid_noslip_contact(fluid, {&wall_noslip});
    // Doc175: classical relations for InteractionSplit damping (CK Inner<>/Contact<> are not BaseInnerRelation).
    UniquePtr<InnerRelation> fluid_inner_classic;
    UniquePtr<ContactRelation> fluid_noslip_contact_classic;
    if (cli.viscosity_mode == "implicit_split" || cli.viscosity_operator_diag)
    {
        fluid_inner_classic = makeUnique<InnerRelation>(fluid);
        fluid_noslip_contact_classic = makeUnique<ContactRelation>(fluid, RealBodyVector{&wall_noslip});
    }

    BaseParticles &particles = fluid.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(particles, cli.t_ref);
    registerOphelieThermalDiffusionAuxFields(particles, cli.k);
    particles.registerStateVariable<Real>(kDoc131QField, Real(0));
    syncVariableToDevice<Real>(particles, kDoc131QField);
    // Thermal fields are particle-attached state: they must follow ParticleSortCK (only evolving variables are
    // permuted). Without this, Temperature/DeltaT/weights would be scrambled against positions after a sort.
    if (!cli.legacy_thermal_sort)
    {
        for (const char *name : {kOphelieTemperatureField, kOphelieThermalDeltaTField, kOphelieThermalLaplaceTField,
                                 kOphelieThermalConductivityField, kOphelieThermalBoundaryMaskField,
                                 kOphelieThermalBoundaryFaceField, kOphelieThermalBoundaryWeightSideField,
                                 kOphelieThermalBoundaryWeightBottomField, kOphelieThermalBoundaryWeightFreeField,
                                 kDoc131QField})
            particles.addEvolvingVariable<Real>(name);
    }
    if (constant_mu_implicit)
    {
        registerOphelieFrenchThermalRheologyFields(particles, cli.mu, cli.cp, cli.k);
        particles.addEvolvingVariable<Real>(kOphelieDynamicViscosityField);
        particles.addEvolvingVariable<Real>(kOphelieHeatCapacityField);
        syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
        syncVariableToDevice<Real>(particles, kOphelieHeatCapacityField);
        if (fluid_particle_viscosity)
            fluid_particle_viscosity->bindParticleViscosityField(particles, kOphelieDynamicViscosityField);
        std::cout << "[doc189][viscosity] constant particle mu=" << cli.mu
                  << " Pa s; variable_properties=0; implicit_split uses this field" << std::endl;
    }

    // ---- Thermal / density init (host, once) ----
    double r2_s_mean = 0.0;
    if (is_R1)
    {
        hostPrescribeStableStratification(particles, cli.t_ref, a_eff, z_bottom);
        // Init only: rho from p'_strat (perturbation_only). Mass synced once; never rewritten.
        if (cli.doc200_strat_init_vs_tmean)
        {
            hostApplyStratifiedHydrostaticDensityVsTmean(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff,
                                                        cli.t_ref, z_bottom, H, /*sync_mass=*/true);
        }
        else
        {
            hostApplyStratifiedHydrostaticDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom,
                                                  H, /*sync_mass=*/true);
        }
    }
    else if (is_R2)
    {
        // With --fixed-T-after-settle, keep T=Tref through packing; prescribe after settle.
        if (!(cli.r2_fixed_spatial_T && cli.fixed_T_after_settle))
        {
            r2_s_mean = static_cast<double>(hostPrescribeRadialTemperature(particles, axis_center, R, cli.t_ref,
                                                                           r2_dT_eff, !cli.r2_keep_mean));
        }
        else
        {
            std::cout << "[doc203][f1] R2 T deferred until after settle (fixed-T-after-settle)" << std::endl;
        }
    }
    else if (is_T2)
    {
        hostPrescribeT2Field(particles, axis_center, R, H, z_bottom, cli.t_ref, cli.a_strat, cli.dT_radial);
    }
    else if (is_T4)
    {
        hostPrescribeT4Field(particles, H, z_bottom, cli.t_ref, cli.t4_amp);
    }
    // R0 / R3 / baseline / T0 / T1 / T3 / T5: T = Tref (registered default).
    // Doc192: particle-wise OriginalID,T_K. Doc164 E: optional CSV T(r,z). Mutually exclusive.
    if (cli.restart_step == 0 && !cli.init_T_particle_csv.empty())
    {
        if (!cli.init_T_rz_csv.empty() || cli.diag_init_T_wall > TinyReal || cli.diag_init_T_core > TinyReal)
        {
            std::cerr << "[doc192] --init-T-particle-csv conflicts with --init-T-rz-csv / --diag-init-T-*"
                      << std::endl;
            return 1;
        }
        std::string terr;
        if (!hostPrescribeTemperatureFromParticleCsv(particles, cli.t_ref, cli.init_T_particle_csv, terr))
        {
            std::cerr << "[doc192] init-T-particle-csv failed: " << terr << std::endl;
            return 1;
        }
        std::cout << "[doc192][init-T] source_note="
                  << (cli.init_T_source_note.empty() ? "unspecified" : cli.init_T_source_note) << std::endl;
    }
    std::vector<RzTempNode> fixed_T_rz_nodes;
    if (cli.fixed_spatial_T && cli.init_T_rz_csv.empty())
    {
        std::cerr << "[doc203][f1] --fixed-spatial-T requires --init-T-rz-csv=path" << std::endl;
        return 1;
    }
    else if (cli.restart_step == 0 && !cli.init_T_rz_csv.empty())
    {
        if (cli.diag_init_T_wall > TinyReal || cli.diag_init_T_core > TinyReal)
        {
            std::cerr << "[doc164] --init-T-rz-csv conflicts with --diag-init-T-*" << std::endl;
            return 1;
        }
        std::string terr;
        const bool apply_csv_now = !(cli.fixed_spatial_T && cli.fixed_T_after_settle);
        if (!hostPrescribeTemperatureFromRzCsv(particles, axis_center, cli.t_ref, cli.init_T_rz_csv, cli.dp * Real(2),
                                               terr, cli.fixed_spatial_T ? &fixed_T_rz_nodes : nullptr,
                                               cli.fixed_spatial_T_demean, apply_csv_now))
        {
            std::cerr << "[doc164] init-T-rz-csv failed: " << terr << std::endl;
            return 1;
        }
        std::cout << "[doc164][init-T] source_note="
                  << (cli.init_T_source_note.empty() ? "unspecified" : cli.init_T_source_note)
                  << " fixed_spatial_T=" << (cli.fixed_spatial_T ? 1 : 0)
                  << " demean=" << (cli.fixed_spatial_T_demean ? 1 : 0)
                  << " after_settle=" << (cli.fixed_T_after_settle ? 1 : 0) << std::endl;
    }
    else if (cli.fixed_spatial_T && cli.restart_step > 0)
    {
        // Restart: still need the CSV nodes for re-application.
        std::string terr;
        if (!hostLoadRzTemperatureCsv(cli.init_T_rz_csv, fixed_T_rz_nodes, terr))
        {
            std::cerr << "[doc203][f1] load fixed T csv on restart failed: " << terr << std::endl;
            return 1;
        }
    }
    // Doc161 P3 diagnostic init (default off). Applied once; then normal evolution.
    // Doc164: never re-apply on --restart-step>0.
    else if (cli.restart_step == 0 && is_baseline && cli.diag_init_T_wall > TinyReal && cli.diag_init_T_core > TinyReal)
    {
        if (!(cli.diag_init_T_delta > TinyReal))
        {
            std::cerr << "[doc161] --diag-init-T-delta must be > 0 when wall/core init is set" << std::endl;
            return 1;
        }
        hostPrescribeWallColdCoreHotTemperature(particles, axis_center, R, cli.t_ref, cli.diag_init_T_wall,
                                                cli.diag_init_T_core, cli.diag_init_T_delta);
        std::cout << "[doc161][diag-init-T] wall_cold_core_hot T_wall=" << cli.diag_init_T_wall
                  << " T_core=" << cli.diag_init_T_core << " delta_m=" << cli.diag_init_T_delta
                  << " (diagnostic only; not literature field)" << std::endl;
    }
    else if (cli.restart_step > 0)
    {
        std::cout << "[doc164][restart] skip diag-init-T / fresh thermal prescription (restart_step="
                  << cli.restart_step << ")" << std::endl;
    }
    if (cli.uniform_buoyancy_balance_init && cli.absorb_mean_buoyancy_bg)
    {
        std::cerr << "[doc186] --uniform-buoyancy-balance-init conflicts with --absorb-mean-buoyancy-bg "
                     "(p' vs p_h double-count)"
                  << std::endl;
        return 1;
    }
    if (cli.uniform_buoyancy_balance_init)
    {
        if (is_R1)
        {
            std::cerr << "[doc186] --uniform-buoyancy-balance-init conflicts with R1 strat init" << std::endl;
            return 1;
        }
        hostApplyUniformBuoyancyBalanceDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, cli.t_ref,
                                               z_bottom, H, /*sync_mass=*/true);
    }
    double mass_sum_init = hostTotalMass(particles);

    double eos_res_init = 0.0, eos_pscale_init = 0.0;
    if (cli.uniform_buoyancy_balance_init)
    {
        // Residual vs particle-wise uniform-buoyancy p'
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "Density");
        syncVariableToHost<Real>(particles, "Temperature");
        const Vecd *pos_e = particles.getVariableDataByName<Vecd>("Position");
        const Real *rho_e = particles.getVariableDataByName<Real>("Density");
        const Real *T_e = particles.getVariableDataByName<Real>("Temperature");
        const size_t n_e = particles.TotalRealParticles();
        eos_res_init = 0.0;
        eos_pscale_init = 0.0;
        for (size_t i = 0; i < n_e; ++i)
        {
            const Real zeta = std::min(H, std::max(Real(0), pos_e[i][2] - z_bottom));
            const double p_an = static_cast<double>(cli.rho0) * static_cast<double>(cli.gravity_g) *
                                static_cast<double>(beta_eff) * static_cast<double>(T_e[i] - cli.t_ref) *
                                static_cast<double>(zeta - H);
            const double p_num = static_cast<double>(cli.c0) * static_cast<double>(cli.c0) *
                                 (static_cast<double>(rho_e[i]) - static_cast<double>(cli.rho0));
            eos_res_init = std::max(eos_res_init, std::abs(p_num - p_an));
            eos_pscale_init = std::max(eos_pscale_init, std::abs(p_an));
        }
    }
    else if (is_R1 && cli.doc200_strat_init_vs_tmean)
    {
        // Recompute Tmean the same way as init (ρ0·Vol proxy before hydro mass rewrite).
        syncVariableToHost<Real>(particles, "Temperature");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *T_m = particles.getVariableDataByName<Real>("Temperature");
        const Real *vol_m = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const size_t n_m = particles.TotalRealParticles();
        double m_sum = 0.0, mT = 0.0;
        for (size_t i = 0; i < n_m; ++i)
        {
            const double mi = static_cast<double>(cli.rho0 * vol_m[i]);
            m_sum += mi;
            mT += mi * static_cast<double>(T_m[i]);
        }
        const Real t_mean_eos = (m_sum > TinyReal) ? static_cast<Real>(mT / m_sum) : cli.t_ref;
        hostEosIdentityResidual(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H, true,
                                eos_res_init, eos_pscale_init, /*vs_tmean=*/true, cli.t_ref, t_mean_eos);
    }
    else
        hostEosIdentityResidual(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H, is_R1,
                                eos_res_init, eos_pscale_init);
    const double real_eps = static_cast<double>(std::numeric_limits<Real>::epsilon());
    const double eos_tol = 8.0 * real_eps * static_cast<double>(cli.rho0) * static_cast<double>(cli.c0) *
                               static_cast<double>(cli.c0) +
                           1.0e-3 * eos_pscale_init;
    const bool eos_initialization_ok = eos_res_init <= eos_tol;

    // ---- Doc131 B3 area audit + weights (host, once; shell thickness = dp) ----
    const AreaAudit area = hostBuildBoundaryAreaWeights(particles, axis_center, R, H, z_bottom, cli.rho0, dp,
                                                        cli.area_normalize, face_top_on, face_side_on,
                                                        face_bottom_on);
    std::cout << "[doc131][area] shell=" << dp << " normalize=" << (cli.area_normalize ? 1 : 0)
              << " top n/raw/target/rel_raw/rel=" << area.top.n << "/" << area.top.raw << "/" << area.top.target
              << "/" << area.top.rel_err_raw << "/" << area.top.rel_err << " side=" << area.side.n << "/"
              << area.side.raw << "/" << area.side.target << "/" << area.side.rel_err_raw << "/"
              << area.side.rel_err << " bottom=" << area.bottom.n << "/" << area.bottom.raw << "/"
              << area.bottom.target << "/" << area.bottom.rel_err_raw << "/" << area.bottom.rel_err
              << " n_corner(primary-face only)=" << area.n_corner << std::endl;

    // ---- Doc134/Doc156: prior EM Q as rz table or EM-particle cloud (discrete-renorm to P) ----
    Doc134RzQGrid rz_q_grid;
    const Doc134RzQGrid *rz_q_ptr = nullptr;
    Doc156EmParticleQCloud em_q_cloud;
    const Doc156EmParticleQCloud *em_q_ptr = nullptr;
    if (cli.q_shape == "em_rz_csv")
    {
        if (!rz_q_grid.load(cli.q_from_rz_csv))
            return 1;
        rz_q_ptr = &rz_q_grid;
    }
    else if (cli.q_shape == "em_particles_csv")
    {
        if (!em_q_cloud.load(cli.q_from_em_particles_csv, static_cast<double>(cli.em_q_support_radius)))
            return 1;
        em_q_ptr = &em_q_cloud;
        // Doc164 B1: dump production C++ samples at listed (r,z) for Python consistency checks.
        if (!cli.diag_q_sample_rz.empty())
        {
            std::ifstream pin(cli.diag_q_sample_rz);
            if (!pin)
            {
                std::cerr << "[doc164] cannot open --diag-q-sample-rz=" << cli.diag_q_sample_rz << std::endl;
                return 1;
            }
            const std::string out_samp =
                (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/diag_q_sample_cpp.csv";
            std::ofstream pout(out_samp);
            pout << "r_m,z_m,Q_Wpm3\n";
            std::string line;
            std::getline(pin, line); // header optional
            auto process = [&](const std::string &row) {
                if (row.empty())
                    return;
                std::stringstream ss(row);
                std::string a, b;
                if (!std::getline(ss, a, ',') || !std::getline(ss, b, ','))
                    return;
                double rr = 0, zz = 0;
                if (!parseRealStrict(a, rr) || !parseRealStrict(b, zz))
                    return;
                const double qq = em_q_cloud.sample(rr, zz);
                pout << std::setprecision(17) << rr << "," << zz << "," << qq << "\n";
            };
            // If first line looked numeric, process it too
            {
                std::stringstream ss(line);
                std::string a, b;
                double rr = 0, zz = 0;
                if (std::getline(ss, a, ',') && std::getline(ss, b, ',') && parseRealStrict(a, rr) &&
                    parseRealStrict(b, zz))
                    process(line);
            }
            while (std::getline(pin, line))
                process(line);
            std::cout << "[doc164][diag-q-sample] wrote " << out_samp << " h=" << cli.em_q_support_radius
                      << std::endl;
        }
    }

    // ---- Doc131 B power contract (host, once) ----
    PowerReport power;
    {
        const double p_target = q_active ? static_cast<double>(cli.q_power) : 0.0;
        if (!hostApplyPowerContract(particles, cli.power_mode, p_target, static_cast<double>(cli.rho0), v_cyl,
                                    power, cli.q_shape, axis_center, static_cast<double>(R), rz_q_ptr, em_q_ptr))
        {
            return 1;
        }
        std::cout << "[doc131][power] mode=" << power.mode << " P_target=" << power.P_target
                  << " P_applied=" << power.P_applied << " gap_rel=" << power.gap_rel << " V_geom=" << power.V_geom
                  << " sum_V_Q=" << power.sum_V_Q << " mass_sum=" << power.mass_sum
                  << " sumVQ/Vgeom=" << power.ratio_sumVQ_over_Vgeom << " q=[" << power.q_min << ","
                  << power.q_max << "] norm_coeff=" << power.norm_coeff
                  << " (norm_coeff=P_target/N_Q is the SPH sampling correction, not raw-EM to 50 kW)" << std::endl;
    }
    std::vector<EmPowerLedgerRow> em_power_ledger;
    auto note_power_stage = [&](const std::string &stage, double t_s, double p_em, double p_part, double p_table,
                                double scale_meta) {
        EmPowerLedgerRow row;
        row.stage = stage;
        row.t_s = t_s;
        row.P_em_raw_W = p_em;
        row.P_particles_W = p_part;
        row.P_table_raw_W = p_table;
        row.P_fluid_mapped_raw_W = power.N_Q;
        row.coeff_sph_to_applied = power.norm_coeff;
        row.P_fluid_applied_W = power.P_applied;
        row.P_target_W = power.P_target;
        if (std::isfinite(scale_meta))
            row.scale_table_raw_to_target = scale_meta;
        else if (std::isfinite(p_table) && std::abs(p_table) > 0.0)
            row.scale_table_raw_to_target = power.P_target / p_table;
        em_power_ledger.push_back(row);
        std::cout << "[doc189][power-stage] stage=" << stage << " t=" << t_s << " P_em_raw=" << p_em
                  << " P_particles=" << p_part << " P_table_raw=" << p_table
                  << " scale_table_raw_to_target=" << row.scale_table_raw_to_target
                  << " P_fluid_mapped_raw=" << row.P_fluid_mapped_raw_W
                  << " coeff_sph_to_applied=" << row.coeff_sph_to_applied
                  << " P_fluid_applied=" << row.P_fluid_applied_W << std::endl;
    };
    note_power_stage("init_pre_em", 0.0, std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN());

    // ---- Thermal boundary options (production struct; used by the CK kernels) ----
    OphelieThermalDiffusionOneWayOptions thermal_opts;
    thermal_opts.enable_diffusion = true;
    thermal_opts.enable_french_natural_bc = true;
    thermal_opts.boundary_shell_thickness = dp;
    thermal_opts.h_side = cli.h_side;
    thermal_opts.h_bottom = cli.h_bottom;
    thermal_opts.h_free = cli.h_top;
    thermal_opts.emissivity = cli.emissivity;
    thermal_opts.t_cool = cli.t_cool;
    thermal_opts.t_ambient = cli.t_ambient;
    thermal_opts.t_rad_ambient = cli.t_rad_ambient;
    thermal_opts.enable_free_crust = false;
    thermal_opts.thermal_wall_model = cli.thermal_wall_model;

    // ---- Dynamics ----
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cell_linked_list(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_closed_cell_linked_list(wall_closed);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_noslip_cell_linked_list(wall_noslip);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>, Contact<>> fluid_update_relation(
        fluid_inner, fluid_closed_contact, fluid_noslip_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(fluid);

    // Gravity object exists ONLY as the g-vector source for BoussinesqBuoyancyForceCK.
    // GravityForceCK is intentionally never constructed or executed (Doc126).
    Gravity gravity(Vecd(0.0, 0.0, -cli.gravity_g));
    StateDynamics<MainExecutionPolicy, BoussinesqBuoyancyForceCK> boussinesq_force(
        fluid, gravity, beta_eff, cli.t_ref, kOphelieTemperatureField);
    if (cli.doc189_dynamic_force_repro)
    {
        const std::string repro_dir =
            cli.output_dir.empty() ? std::string("doc189_dynamic_force_repro") : cli.output_dir;
        return runDoc189DynamicForceRepro(boussinesq_force, particles, beta_eff, cli.gravity_g, R, axis_center,
                                          cli.t_ref, repro_dir);
    }
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_closed_normal(wall_closed);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_noslip_normal(wall_noslip);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> advection_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_position(fluid);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex> fluid_linear_correction_matrix(
        DynamicsArgs(fluid_inner, 0.5), fluid_closed_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK> acoustic_1st(
        fluid_inner, fluid_closed_contact);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalf<
                                                   Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd_wall(fluid_closed_contact);
    acoustic_2nd.addPostContactInteraction(acoustic_2nd_wall);
    // Doc201 B1 wall ablation: Inner-only acoustic (no closed-wall Riemann contact).
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep1stHalf<
                              Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_1st_nowall(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<
                              Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd_nowall(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        acoustic_1st_corr(fluid_inner, fluid_closed_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        acoustic_2nd_corr(fluid_inner, fluid_closed_contact);

    // Optional Internal regularization path (closed domain; never FreeSurface).
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>, Contact<>>>
        density_summation(fluid_inner, fluid_closed_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, Internal>> density_reg(fluid);

    // Viscous: Inner + wall_noslip ONLY (lid omitted -> flat free-slip top).
    // Doc143: with --variable-properties (default), use ParticleViscosity pairwise μ; legacy --mean-mu-viscous kept.
    using ViscousForceLocalCK =
        fluid_dynamics::ViscousForceCK<Inner<WithUpdate, ParticleViscosity, NoKernelCorrectionCK>,
                                       Contact<Wall, ParticleViscosity, NoKernelCorrectionCK>>;
    using ViscousForceInnerLocalCK =
        fluid_dynamics::ViscousForceCK<Inner<WithUpdate, ParticleViscosity, NoKernelCorrectionCK>>;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>> viscous_force_mean;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalCK>> viscous_force_local;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceInnerLocalCK>> viscous_force_inner_diag;
    auto invalidate_viscous_kernels = [&]() {
        fluid_inner.resetComputingKernelUpdated();
        const StdVec<SPHBody *> contacts = fluid_noslip_contact.getContactBodies();
        for (size_t k = 0; k < contacts.size(); ++k)
            fluid_noslip_contact.resetComputingKernelUpdated(k);
    };
    auto exec_viscous_force = [&]() {
        if (viscous_force_local)
            viscous_force_local->exec();
        else
            viscous_force_mean->exec();
    };

    // Doc139 B1 / Doc143: optional T-dependent μ/cp/k. Local μ enters ViscousForceCK pairwise when enabled.
    OphelieFrenchGlassMaterialLaws material_laws = makeJacoutot1473KReferenceMaterialLaws();
    material_laws.rho.clear();
    material_laws.rho.setPoint(0, cli.t_ref, cli.rho0);
    material_laws.beta.clear();
    material_laws.beta.setPoint(0, cli.t_ref, cli.beta);
    material_laws.mu.clear();
    material_laws.mu.setPoint(0, cli.t_ref, cli.mu);
    material_laws.mu.logarithmic_value_ = true;
    material_laws.cp.clear();
    material_laws.cp.setPoint(0, cli.t_ref, cli.cp);
    material_laws.k.clear();
    material_laws.k.setPoint(0, cli.t_ref, cli.k);
    OphelieThermalRheologyStats rheology_stats{};
    ThesisMaterialCurveLoadReport thesis_curve_report{};
    LocalViscosityForceEvidence local_mu_evidence{};
    bool local_mu_evidence_recorded = false;
    MarangoniApplyResult last_marangoni_apply{};
    MarangoniVerifyResult marangoni_verify_result{};
    bool marangoni_verify_ran = false;
    std::string variable_props_curve_note = "constant_table_or_cli";
    if (cli.variable_properties)
    {
        const bool require_thesis_anchored = (cli.material_mode == "cep_thesis_NUM_shape_anchored");
        if (require_thesis_anchored || !cli.material_curve_dir.empty())
        {
            if (cli.material_curve_dir.empty())
            {
                std::cerr << "[doc156] material-curve-dir required for thesis material mode (CES fallback forbidden)"
                          << std::endl;
                return 1;
            }
            thesis_curve_report = installThesisIII8910ThermalRheologyFromDigitizedCsv(
                material_laws, cli.material_curve_dir, require_thesis_anchored, cli.mu, cli.cp, cli.k);
            if (!thesis_curve_report.ok)
            {
                std::cerr << "[doc139] material-curve-dir load failed: k=" << thesis_curve_report.k_msg
                          << " cp=" << thesis_curve_report.cp_msg << " mu=" << thesis_curve_report.mu_msg
                          << std::endl;
                return 1;
            }
            variable_props_curve_note =
                require_thesis_anchored
                    ? "cep_thesis_NUM_shape_anchored; III.8/9/10 digitized; CEP1473K scales; then --mu-cap"
                    : "thesis_III8_k_III9_cp_III10_mu_digitized; endpoint_hold; mu clamped to --mu-cap";
            std::cout << "[doc139][material-curves] loaded thesis NUM from " << cli.material_curve_dir
                      << " anchored=" << (thesis_curve_report.cep_anchored ? 1 : 0)
                      << " @1473K raw k/cp/mu=" << thesis_curve_report.k_1473_raw << "/"
                      << thesis_curve_report.cp_1473_raw << "/" << thesis_curve_report.mu_1473_raw
                      << " s_k/s_cp/s_mu=" << thesis_curve_report.s_k << "/" << thesis_curve_report.s_cp << "/"
                      << thesis_curve_report.s_mu << " @1473K after k/cp/mu=" << thesis_curve_report.k_1473
                      << "/" << thesis_curve_report.cp_1473 << "/" << thesis_curve_report.mu_1473
                      << " | @303K: k=" << thesis_curve_report.k_303 << " cp=" << thesis_curve_report.cp_303
                      << " mu=" << thesis_curve_report.mu_303 << " | " << thesis_curve_report.interpolation
                      << std::endl;
        }
        else if (cli.material_mode == "legacy_CES_shape_control" || cli.material_mode.empty())
        {
            installJacoutotCes2008Fig2ThermalRheology(material_laws);
            variable_props_curve_note =
                "borrowed_CES_2008_Fig2_NUM_shape; set --material-mode=cep_thesis_NUM_shape_anchored for Doc156 D2";
        }
        else
        {
            std::cerr << "[doc156] unsupported material_mode with --variable-properties: " << cli.material_mode
                      << std::endl;
            return 1;
        }
        registerOphelieFrenchThermalRheologyFields(particles, cli.mu, cli.cp, cli.k);
        if (!cli.legacy_thermal_sort)
        {
            particles.addEvolvingVariable<Real>(kOphelieDynamicViscosityField);
            particles.addEvolvingVariable<Real>(kOphelieHeatCapacityField);
        }
        if (fluid_particle_viscosity)
            fluid_particle_viscosity->bindParticleViscosityField(particles, kOphelieDynamicViscosityField);
        const Real mu_init = hostUpdateFrenchThermalRheologyFromTemperature(
            particles, material_laws, /*constant_properties=*/false, cli.t_ref, rheology_stats);
        hostClampParticleViscosityCap(particles, cli.mu_cap, rheology_stats);
        // Reference μ: for local path store mu_max (viscous CFL); for mean path store volume mean.
        const Real mu_ref_for_material =
            use_local_particle_mu
                ? (rheology_stats.mu_used_max > TinyReal ? rheology_stats.mu_used_max : mu_init)
                : (rheology_stats.mu_used_mean > TinyReal ? rheology_stats.mu_used_mean : mu_init);
        fluid_viscosity.setReferenceViscosity(mu_ref_for_material);
        std::cout << "[doc139][variable-properties] curves=" << variable_props_curve_note
                  << " mu_raw_max=" << rheology_stats.mu_raw_max
                  << " mu_used_mean=" << rheology_stats.mu_used_mean
                  << " mu_used_max=" << rheology_stats.mu_used_max
                  << " mu_cap=" << cli.mu_cap << " n_mu_capped=" << rheology_stats.n_mu_capped
                  << " cp_mean=" << rheology_stats.cp_mean << " k_mean=" << rheology_stats.k_mean
                  << " viscous_mu_mode=" << (use_local_particle_mu ? "particle_pairwise" : "volume_mean")
                  << std::endl;
    }

    if (use_local_particle_mu)
        viscous_force_local =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalCK>>(fluid_inner,
                                                                                        fluid_noslip_contact);
    else
        viscous_force_mean =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>>(
                fluid_inner, fluid_noslip_contact);
    invalidate_viscous_kernels();

    // Doc175: classical InteractionSplit damping (CPU). Constructed after ParticleViscosity field exists.
    const bool use_implicit_split = (cli.viscosity_mode == "implicit_split");
    const bool need_classic_visc =
        use_implicit_split || cli.viscosity_operator_diag;
    UniquePtr<InteractionSplit<DampingPairwiseInner<Vecd, ParticleViscosityHarmonicDampingRate>>>
        implicit_visc_inner;
    UniquePtr<InteractionSplit<DampingPairwiseContactWall<Vecd, ParticleViscosityOneSideDampingRate>>>
        implicit_visc_wall;
    size_t n_implicit_visc_exec = 0;
    double implicit_visc_wall_s = 0.0;
    Real last_dt_visc_est = 0.0;
    std::string last_active_limiter = "none";
    if (need_classic_visc)
    {
        if (!use_local_particle_mu)
        {
            std::cerr << "[doc177] implicit/diag requires bound ParticleViscosity field" << std::endl;
            return 1;
        }
        // Ensure ViscousForce / PreviousViscousForce exist (registered by ViscousForceCK ctor above).
        hostClearViscousForcePriorContribution(particles);
        if (!fluid_inner_classic || !fluid_noslip_contact_classic)
        {
            std::cerr << "[doc177] classic InnerRelation/ContactRelation missing" << std::endl;
            return 1;
        }
        implicit_visc_inner =
            makeUnique<InteractionSplit<DampingPairwiseInner<Vecd, ParticleViscosityHarmonicDampingRate>>>(
                *fluid_inner_classic, "Velocity");
        implicit_visc_wall =
            makeUnique<InteractionSplit<DampingPairwiseContactWall<Vecd, ParticleViscosityOneSideDampingRate>>>(
                *fluid_noslip_contact_classic, "Velocity");
        if (use_implicit_split)
            std::cout << "[doc177][viscosity] mode=implicit_split backend=InteractionSplit+DampingPairwise "
                         "(CPU ParallelPolicy); ViscousForceCK disabled; dt_visc limiter off; "
                         "CK-match DampingRate factors inner="
                      << kOphelieCkMatchInnerDampingFactor << " wall=" << kOphelieCkMatchWallDampingFactor
                      << "; host sync Velocity/μ each acoustic substep" << std::endl;
        else
            std::cout << "[doc177][viscosity] operator-diag: both ViscousForceCK and DampingPairwise built"
                      << std::endl;
            if (cli.viscosity_diag_inner_only)
            {
                viscous_force_inner_diag =
                    makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceInnerLocalCK>>(fluid_inner);
                std::cout << "[doc177][viscosity] operator-diag: Inner-only ViscousForceCK enabled" << std::endl;
            }
    }
    else
    {
        std::cout << "[doc175][viscosity] mode=explicit ViscousForceCK" << std::endl
                  << std::flush;
    }
    std::cout << "[doc203][f1-heartbeat] post-viscosity: building time-step/thermal/VTP objects..."
              << std::endl
              << std::flush;
    auto update_classic_visc_relations = [&]() {
        if (!need_classic_visc || !fluid_inner_classic || !fluid_noslip_contact_classic)
            return;
        // Doc177 bridge: classic InnerRelation::searchNeighborsByMesh reads cell_data_lists_,
        // while CK UpdateCellLinkedList only fills dv_particle_index_/dv_cell_offset_.
        // Without RealBody::updateCellLinkedList(), classic neighborhoods stay empty →
        // DampingPairwise is a no-op (Doc177 operator-diag classic_nb_mean=0).
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Vecd>(wall_noslip.getBaseParticles(), "Position");
        fluid.updateCellLinkedList();
        wall_noslip.updateCellLinkedList();
        fluid_inner_classic->updateConfiguration();
        fluid_noslip_contact_classic->updateConfiguration();
    };

    auto exec_implicit_viscous_split = [&](Real dt_visc_step) {
        if (!use_implicit_split || !(dt_visc_step > TinyReal))
            return;
        const auto t0 = std::chrono::steady_clock::now();
        syncVariableToHost<Vecd>(particles, "Velocity");
        syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
        syncVariableToHost<Real>(particles, "Mass");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        // Wall velocity (usually zero) for no-slip damping contact.
        syncVariableToHost<Vecd>(wall_noslip.getBaseParticles(), "Velocity");
        syncVariableToHost<Real>(wall_noslip.getBaseParticles(), "VolumetricMeasure");
        implicit_visc_inner->exec(dt_visc_step);
        implicit_visc_wall->exec(dt_visc_step);
        syncVariableToDevice<Vecd>(particles, "Velocity");
        const auto t1 = std::chrono::steady_clock::now();
        implicit_visc_wall_s += std::chrono::duration<double>(t1 - t0).count();
        ++n_implicit_visc_exec;
    };

    // Doc144: Marangoni axisymmetric verify (prescribed hot-center T; restores state). Default on if dγ/dT≠0.
    if (!cli.marangoni_verify_set)
        cli.marangoni_verify = std::abs(cli.marangoni_dgamma_dt) > TinyReal;
    if (cli.marangoni_verify && std::abs(cli.marangoni_dgamma_dt) > TinyReal)
    {
        marangoni_verify_result = hostVerifyMarangoniAxisymmetric(particles, axis_center, R, dp,
                                                                  cli.marangoni_dgamma_dt, cli.t_ref);
        marangoni_verify_ran = true;
        std::cout << "[doc144][marangoni-verify] passed=" << marangoni_verify_result.passed
                  << " direction_ok=" << marangoni_verify_result.direction_ok
                  << " area_integral_ok=" << marangoni_verify_result.area_integral_ok
                  << " rel_err=" << marangoni_verify_result.rel_err_integral
                  << " sum_Fr=" << marangoni_verify_result.sum_Fr
                  << " theory=" << marangoni_verify_result.sum_Fr_theory
                  << " sum_A=" << marangoni_verify_result.sum_area
                  << " method=" << marangoni_verify_result.method << std::endl;
    }

    std::cout << "[doc203][f1-heartbeat] constructing Advection/Acoustic/Umax ReduceDynamicsCK..."
              << std::endl
              << std::flush;
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> advection_dt(fluid, cli.u_ref,
                                                                                            cli.advection_cfl);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> acoustic_dt(fluid, cli.acoustic_cfl);
    ReduceDynamicsCK<MainExecutionPolicy, OphelieVelocityMaxReduceCK> u_max_reduce(fluid);
    std::cout << "[doc203][f1-heartbeat] ReduceDynamicsCK done; building thermal UniquePtrs..."
              << std::endl
              << std::flush;

    // Doc131 C / Doc139: production thermal kernels. Variable-cp path uses per-particle OphelieHeatCapacity.
    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace_temperature(
        fluid_inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
        thermal_opts.pair_weight_regularization);
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>> apply_diffusion;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtParticleCpCK>> apply_diffusion_var_cp;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyDoc131SourceQCK>> apply_source_q;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyJouleHeatWithParticleCpCK>> apply_source_q_var_cp;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK>> apply_natural_bc;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcParticleCpCK>> apply_natural_bc_var_cp;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcCK>> apply_wall_bc;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcParticleCpCK>> apply_wall_bc_var_cp;
    UniquePtr<ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK>> boundary_loss_legacy;
    UniquePtr<ReduceDynamicsCK<MainExecutionPolicy, OphelieWallTemperatureBoundaryLossPowerReduceCK>> boundary_loss_wall;
    if (cli.variable_properties)
    {
        apply_diffusion_var_cp =
            makeUnique<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtParticleCpCK>>(fluid, cli.rho0,
                                                                                                cli.t_ref);
        apply_source_q_var_cp = makeUnique<StateDynamics<MainExecutionPolicy, ApplyJouleHeatWithParticleCpCK>>(
            fluid, cli.rho0, cli.t_ref, kDoc131QField);
        if (use_wall_T)
            apply_wall_bc_var_cp =
                makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcParticleCpCK>>(
                    fluid, kOphelieThermalDeltaTField, kOphelieTemperatureField, cli.t_ref, cli.rho0,
                    thermal_opts.boundary_shell_thickness, thermal_opts);
        else
            apply_natural_bc_var_cp =
                makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcParticleCpCK>>(
                    fluid, kOphelieThermalDeltaTField, kOphelieTemperatureField, cli.t_ref, cli.rho0,
                    thermal_opts.boundary_shell_thickness, thermal_opts);
    }
    else
    {
        apply_diffusion =
            makeUnique<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>>(fluid, cli.rho0, cli.cp,
                                                                                      cli.t_ref);
        apply_source_q = makeUnique<StateDynamics<MainExecutionPolicy, ApplyDoc131SourceQCK>>(
            fluid, kDoc131QField, cli.rho0, cli.cp, cli.t_ref);
        if (use_wall_T)
            apply_wall_bc = makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcCK>>(
                fluid, kOphelieThermalDeltaTField, kOphelieTemperatureField, cli.t_ref, cli.rho0, cli.cp,
                thermal_opts.boundary_shell_thickness, thermal_opts);
        else
            apply_natural_bc = makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK>>(
                fluid, kOphelieThermalDeltaTField, kOphelieTemperatureField, cli.t_ref, cli.rho0, cli.cp,
                thermal_opts.boundary_shell_thickness, thermal_opts);
    }
    if (use_wall_T)
        boundary_loss_wall =
            makeUnique<ReduceDynamicsCK<MainExecutionPolicy, OphelieWallTemperatureBoundaryLossPowerReduceCK>>(
                fluid, cli.t_ref, thermal_opts.boundary_shell_thickness, thermal_opts);
    else
        boundary_loss_legacy =
            makeUnique<ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK>>(
                fluid, cli.t_ref, thermal_opts.boundary_shell_thickness, thermal_opts);

    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    std::cout << "[doc203][f1-heartbeat] thermal objects ready; constructing SPHSolver/VTP/RestartIO..."
              << std::endl
              << std::flush;

    // Doc133: VTP field output (pattern of frozen_q). Files go to ./output (IO environment) under the binary cwd.
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &fluid_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(fluid);
    fluid_state_recorder.addToWrite<Vecd>(fluid, "Velocity");
    fluid_state_recorder.addToWrite<Real>(fluid, "Pressure");
    fluid_state_recorder.addToWrite<Real>(fluid, "Density");
    fluid_state_recorder.addToWrite<Real>(fluid, "BackgroundPressure");
    fluid_state_recorder.addToWrite<Real>(fluid, kOphelieTemperatureField);
    fluid_state_recorder.addToWrite<Real>(fluid, "Mass");
    fluid_state_recorder.addToWrite<Real>(fluid, "VolumetricMeasure");
    fluid_state_recorder.addToWrite<Real>(fluid, "DensityChangeRate");
    fluid_state_recorder.addToWrite<Real>(fluid, kDoc131QField); // volumetric heat source Q [W/m^3]
    if (use_wall_T)
    {
        fluid_state_recorder.addToWrite<Real>(fluid, kOphelieThermalWallTemperatureField);
        fluid_state_recorder.addToWrite<Real>(fluid, kOphelieThermalWallGapField);
    }
    if (cli.variable_properties)
        fluid_state_recorder.addToWrite<Real>(fluid, kOphelieDynamicViscosityField);
    const bool write_vtp = cli.state_record_interval > TinyReal;
    size_t n_vtp_written = 0;
    Real next_vtp_time = 0.0;
    std::vector<std::string> vtp_times_written;
    // Doc164 D: RestartIOCK + EM sidecar (schema v1).
    if (cli.restart_step > 0)
        sph_system.setRestartStep(static_cast<size_t>(cli.restart_step));
    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);
    const bool restarting = cli.restart_step > 0;
    const std::string em_state_path =
        (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/em_coupling_state.json";
    Real next_checkpoint_time =
        (cli.checkpoint_interval > TinyReal) ? cli.checkpoint_interval : Real(1.0e300);
    size_t last_checkpoint_step = 0;

    auto write_vtp_now = [&](size_t step) {
        if (!write_vtp)
            return;
        fluid_state_recorder.writeToFile(step);
        ++n_vtp_written;
        std::ostringstream os;
        os << step << ":" << sv_physical_time->getValue();
        vtp_times_written.push_back(os.str());
    };


    std::cout << "[doc203][f1-heartbeat] running wall normal StateDynamics..." << std::endl << std::flush;
    wall_closed_normal.exec();
    wall_noslip_normal.exec();
    std::cout << "[doc203][f1-heartbeat] wall normals done; tagging free-slip lid..." << std::endl
              << std::flush;
    // Doc201 B2: acoustic dens free-slip on flat lid (viscous already omits lid via wall_noslip).
    {
        BaseParticles &wall_closed_particles = wall_closed.getBaseParticles();
        const size_t n_lid = hostTagAcousticWallFreeSlipLid(wall_closed_particles, z_bottom + H, wall_thickness);
        std::cout << "[doc201][b2] AcousticWallFreeSlip lid tags=" << n_lid
                  << " / wall_closed_n=" << wall_closed_particles.TotalRealParticles() << std::endl;
        const int n_doc202_ablate =
            (cli.doc202_b2_ablate_side_acoustic ? 1 : 0) + (cli.doc202_b2_ablate_lid_acoustic ? 1 : 0) +
            (cli.doc202_b2_ablate_lid_acoustic_1st ? 1 : 0) + (cli.doc202_b2_ablate_lid_acoustic_2nd ? 1 : 0) +
            (cli.doc202_b2_lid_continuity_silent ? 1 : 0) + (cli.doc202_b2_lid_omit_pjump ? 1 : 0);
        if (n_doc202_ablate > 1)
        {
            std::cerr << "[doc202][b2] refuse combining Doc202 ablate-*/lid-continuity-silent switches"
                      << std::endl;
            std::exit(2);
        }
        if (cli.doc202_b2_ablate_side_acoustic)
        {
            const size_t n_omit =
                hostTagAcousticWallContactOmitSideBottom(wall_closed_particles, z_bottom + H, wall_thickness);
            std::cout << "[doc202][b2] AcousticWallContactOmit side/bottom tags=" << n_omit
                      << " (lid kept; H1 ablate-side-acoustic)" << std::endl;
        }
        if (cli.doc202_b2_ablate_lid_acoustic)
        {
            const size_t n_omit =
                hostTagAcousticWallContactOmitLid(wall_closed_particles, z_bottom + H, wall_thickness);
            std::cout << "[doc202][b2] AcousticWallContactOmit lid tags=" << n_omit
                      << " (side/bottom kept; H1b ablate-lid-acoustic)" << std::endl;
        }
        if (cli.doc202_b2_ablate_lid_acoustic_1st)
        {
            const size_t n_omit = hostTagAcousticWallContactOmitLidHalf(wall_closed_particles, z_bottom + H,
                                                                       wall_thickness, "1st");
            std::cout << "[doc202][b2] AcousticWallContactOmit1st lid tags=" << n_omit
                      << " (H1c-a; 2nd-half lid kept)" << std::endl;
        }
        if (cli.doc202_b2_ablate_lid_acoustic_2nd)
        {
            const size_t n_omit = hostTagAcousticWallContactOmitLidHalf(wall_closed_particles, z_bottom + H,
                                                                       wall_thickness, "2nd");
            std::cout << "[doc202][b2] AcousticWallContactOmit2nd lid tags=" << n_omit
                      << " (H1c-b; 1st-half lid kept)" << std::endl;
        }
        if (cli.doc202_b2_lid_continuity_silent)
        {
            const size_t n_cs =
                hostTagAcousticWallContSilentLid(wall_closed_particles, z_bottom + H, wall_thickness);
            std::cout << "[doc202][b2] AcousticWallContSilent lid tags=" << n_cs
                      << " (H2a continuity-silent; PJump kept)" << std::endl;
        }
        if (cli.doc202_b2_lid_omit_pjump)
        {
            const size_t n_pj =
                hostTagAcousticWallOmitPJumpLid(wall_closed_particles, z_bottom + H, wall_thickness);
            std::cout << "[doc202][b2] AcousticWallOmitPJump lid tags=" << n_pj
                      << " (H2b omit PJump; continuity kept)" << std::endl;
        }
    }
    std::cout << "[doc203][f1-heartbeat] updating cell-linked lists + relations (first GPU kernels)..."
              << std::endl
              << std::flush;
    fluid_cell_linked_list.exec();
    wall_closed_cell_linked_list.exec();
    wall_noslip_cell_linked_list.exec();
    fluid_update_relation.exec();
    update_classic_visc_relations();
    std::cout << "[doc203][f1-heartbeat] cell/relation update done" << std::endl << std::flush;

    if (cli.viscosity_operator_diag)
    {
        using electromagnetics::ophelie::hostCompareViscosityOperators;
        // Constant-μ sinusoid probe: Ux=A*sin(kz), ∇²U=−k²U (linear shear has continuum a≡0 → wall-dominated).
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Vecd>(particles, "Velocity");
        syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
        Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
        Real *mu_f = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        const size_t nn = particles.TotalRealParticles();
        const Real k_wave = Real(2) * Pi / std::max(H, TinyReal);
        for (size_t i = 0; i < nn; ++i)
        {
            mu_f[i] = cli.mu;
            const Real phase = k_wave * (pos[i][2] - z_bottom);
            vel[i] = Vecd(cli.viscosity_diag_shear * std::sin(phase), Real(0), Real(0));
        }
        syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
        syncVariableToDevice<Vecd>(particles, "Velocity");
        fluid_viscosity.setReferenceViscosity(cli.mu);
        invalidate_viscous_kernels();
        hostClearViscousForcePriorContribution(particles);

        // Classic neighbor count smoke
        size_t nb_sum = 0, nb_max = 0;
        for (size_t i = 0; i < nn; ++i)
        {
            const size_t nnb = fluid_inner_classic->inner_configuration_[i].current_size_;
            nb_sum += nnb;
            nb_max = std::max(nb_max, nnb);
        }
        std::cout << "[doc177][viscosity-diag] classic_inner_nb mean="
                  << (nn > 0 ? static_cast<double>(nb_sum) / static_cast<double>(nn) : 0.0)
                  << " max=" << nb_max << " inner_only=" << (cli.viscosity_diag_inner_only ? 1 : 0)
                  << " wave_k=" << k_wave << std::endl;

        const std::string diag_dir = cli.output_dir.empty() ? std::string(".") : cli.output_dir;
        std::filesystem::create_directories(diag_dir);
        const std::string csv_path = diag_dir + "/viscosity_operator_diag.csv";
        Real mean_rel = 0, max_rel = 0;
        const bool ok = hostCompareViscosityOperators(
            particles,
            [&]() {
                hostClearViscousForcePriorContribution(particles);
                if (cli.viscosity_diag_inner_only)
                {
                    if (!viscous_force_inner_diag)
                        throw std::runtime_error("viscosity-diag-inner-only requires viscous_force_inner_diag");
                    viscous_force_inner_diag->exec();
                }
                else
                    exec_viscous_force();
            },
            [&](Real dt_probe) {
                syncVariableToHost<Vecd>(particles, "Velocity");
                syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
                syncVariableToHost<Real>(particles, "Mass");
                syncVariableToHost<Real>(particles, "VolumetricMeasure");
                syncVariableToHost<Vecd>(wall_noslip.getBaseParticles(), "Velocity");
                syncVariableToHost<Real>(wall_noslip.getBaseParticles(), "VolumetricMeasure");
                implicit_visc_inner->exec(dt_probe);
                if (!cli.viscosity_diag_inner_only)
                    implicit_visc_wall->exec(dt_probe);
                syncVariableToDevice<Vecd>(particles, "Velocity");
            },
            cli.viscosity_diag_dt, csv_path, /*max_rows=*/64, mean_rel, max_rel);

        const std::string verdict_path = diag_dir + "/VISCOSITY_OPERATOR_DIAG_VERDICT.json";
        std::ofstream vj(verdict_path);
        vj << "{\n"
           << "  \"ok\": " << (ok ? "true" : "false") << ",\n"
           << "  \"dt\": " << cli.viscosity_diag_dt << ",\n"
           << "  \"amplitude\": " << cli.viscosity_diag_shear << ",\n"
           << "  \"inner_only\": " << (cli.viscosity_diag_inner_only ? "true" : "false") << ",\n"
           << "  \"mu\": " << cli.mu << ",\n"
           << "  \"ck_match_inner\": " << kOphelieCkMatchInnerDampingFactor << ",\n"
           << "  \"ck_match_wall\": " << kOphelieCkMatchWallDampingFactor << ",\n"
           << "  \"mean_rel_err\": " << mean_rel << ",\n"
           << "  \"max_rel_err\": " << max_rel << ",\n"
           << "  \"classic_nb_mean\": "
           << (nn > 0 ? static_cast<double>(nb_sum) / static_cast<double>(nn) : 0.0) << ",\n"
           << "  \"classic_nb_max\": " << nb_max << ",\n"
           << "  \"csv\": \"" << csv_path << "\",\n"
           << "  \"gate_mean_rel_le_0p2\": " << (mean_rel <= Real(0.2) ? "true" : "false") << "\n"
           << "}\n";
        std::cout << "[doc177][viscosity-diag] verdict=" << verdict_path << std::endl;
        return (ok && mean_rel <= Real(0.2)) ? 0 : 2;
    }

    // Reduced-pressure contract. ReferenceBodyAcceleration=0 because g is absorbed in BackgroundPressure
    // and is NOT in ForcePrior (C_ext=0 => C_ref must be 0).
    // AcousticDensMode default 1; Doc186 mode 2 = dens as 1, wall momentum omits C_ext (use with consistent p' hydro).
    hostSetPerturbationMomentum(particles, Real(1));
    hostSetAcousticDensMode(particles, static_cast<Real>(cli.acoustic_dens_mode));
    hostSetReferenceBodyAcceleration(particles, Vecd::Zero());
    std::cout << "[doc186] acoustic_dens_mode=" << cli.acoustic_dens_mode
              << (cli.acoustic_dens_mode == 2
                      ? " (wall mom omits C_ext)"
                      : (cli.acoustic_dens_mode == 3 ? " (Doc201: omit p dens DissipativeUJump)" : ""))
              << std::endl;
    Real last_t_mean_buoy = cli.t_ref;
    Real last_force_tref = cli.t_ref;
    auto apply_background_hydro = [&]() {
        if (cli.absorb_mean_buoyancy_bg)
        {
            last_t_mean_buoy = hostMassWeightedMeanTemperature(particles);
            last_force_tref = last_t_mean_buoy;
            hostApplyBackgroundHydroPressureWithMeanBuoyancy(particles, cli.rho0, cli.gravity_g, beta_eff,
                                                             last_t_mean_buoy, cli.t_ref, z_bottom, H);
            boussinesq_force.setForceReferenceTemperature(last_force_tref);
        }
        else
        {
            last_force_tref = cli.t_ref;
            hostApplyBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z_bottom, H);
            boussinesq_force.setForceReferenceTemperature(cli.t_ref);
        }
    };
    long doc189_param_version = 0;
    double doc189_snapshot_t = 0;
    auto set_cref_to_buoyancy = [&]() {
        syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
        syncVariableToHost<Real>(particles, "Mass");
        syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
        const Vecd *fb = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
        const Real *mm = particles.getVariableDataByName<Real>("Mass");
        Vecd *aref = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
        const size_t nn = particles.TotalRealParticles();
        for (size_t i = 0; i < nn; ++i)
            aref[i] = fb[i] / mm[i];
        syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
    };
    // One snapshot: mass-mean T, force reference, background pressure, buoyancy, wall C_ref.
    // Acoustic substeps reuse this snapshot; they do not recompute T_mean.
    //
    // Doc201 B1: with --doc200-strat-init-vs-tmean and frozen T, p'/ρ balances gβ(T−Tmean).
    // ForcePrior would double-count. Once --thaw-T evolves T, that balance is obsolete → default
    // re-enable Boussinesq ForcePrior unless --doc201-keep-strat-owns-when-thaw (diagnostic).
    const bool strat_p_prime_owns_buoyancy =
        cli.doc200_strat_init_vs_tmean && !(cli.thaw_T && !cli.doc201_keep_strat_owns_when_thaw);
    if (cli.doc200_strat_init_vs_tmean)
    {
        std::cout << "[doc201][b1] strat_init=1 thaw_T=" << (cli.thaw_T ? 1 : 0)
                  << " strat_p_prime_owns_buoyancy=" << (strat_p_prime_owns_buoyancy ? 1 : 0);
        if (strat_p_prime_owns_buoyancy)
            std::cout << " (ForcePrior off; p'/ρ holds gβ(T-Tmean) while T frozen or keep-owns)";
        else if (cli.thaw_T)
            std::cout << " (ForcePrior on; recommended under thaw-T)";
        else
            std::cout << " (ForcePrior on)";
        std::cout << std::endl;
    }
    auto refresh_buoyancy_snapshot = [&](const char *reason) {
        if (cli.doc189_mechanism)
            hostApplyDoc189MechanismTemperature(particles, axis_center, R, cli.t_ref, cli.doc189_mechanism_amp);
        apply_background_hydro();
        if (strat_p_prime_owns_buoyancy)
        {
            syncVariableToHost<Vecd>(particles, "ForcePrior");
            Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
            const size_t nn = particles.TotalRealParticles();
            for (size_t i = 0; i < nn; ++i)
                fp[i] = Vecd::Zero();
            syncVariableToDevice<Vecd>(particles, "ForcePrior");
            hostSetReferenceBodyAcceleration(particles, Vecd::Zero());
        }
        else
        {
            boussinesq_force.exec();
            if (cli.cref_buoyancy)
                set_cref_to_buoyancy();
        }
        ++doc189_param_version;
        doc189_snapshot_t = static_cast<double>(sv_physical_time->getValue());
        std::cout << "[doc189][buoyancy-snapshot] reason=" << reason << " t=" << doc189_snapshot_t
                  << " T_mean=" << last_t_mean_buoy << " force_tref=" << last_force_tref
                  << " eos_tref=" << cli.t_ref << " version=" << doc189_param_version
                  << " strat_p_prime_owns_buoyancy=" << (strat_p_prime_owns_buoyancy ? 1 : 0) << std::endl;
    };
    if (cli.doc189_mechanism)
    {
        std::cout << "[doc189] mechanism amp_K=" << cli.doc189_mechanism_amp
                  << " Tc=eos_tref=" << cli.t_ref << " energy_gate=NOT_APPLICABLE" << std::endl;
    }
    refresh_buoyancy_snapshot(restarting ? "init_before_restart_load" : "init");
    if (cli.absorb_mean_buoyancy_bg)
    {
        std::cout << "[doc186] absorb-mean-buoyancy-bg T_mean=" << last_t_mean_buoy
                  << " force_tref=" << last_force_tref << " eos_tref=" << cli.t_ref << std::endl;
    }

    const OphelieBoussinesqForceAudit force_audit = hostAuditBoussinesqForcePrior(
        particles, beta_eff, last_force_tref, -cli.gravity_g, kOphelieTemperatureField);
    logOphelieBoussinesqForceAudit(force_audit, "doc126_init");
    const double buoy_a_rms_init = hostBuoyancyAccelRms(particles);

    const Real z_top = z_bottom + H;
    ThermalLedger ledger;
    auto thermal_enthalpy_J = [&]() -> double {
        if (cli.variable_properties)
            return static_cast<double>(hostThermalEnthalpyEnergy(particles, cli.rho0, material_laws.cp, cli.t_ref,
                                                                kOphelieTemperatureField));
        return hostEnthalpyDeltaT(particles, cli.cp);
    };
    double h_init = thermal_enthalpy_J();
    const double h_abs_scale_init =
        cli.variable_properties ? std::abs(h_init) + 1.0 : hostEnthalpyAbsScale(particles, cli.cp);
    const DeltaTStats dT_init = hostDeltaTStats(particles);
    double t4_c0 = 0.0, t4_a0 = 0.0, t4_res0 = 0.0;
    if (is_T4)
        hostFitCosineProfile(particles, H, z_bottom, t4_c0, t4_a0, t4_res0);

    ensureParentDir(cli.out_csv);
    ensureParentDir(cli.summary_json);
    // Doc136: half-step expected-write ledgers (only with --drift-probe).
    double cum_drho_diss = 0.0, cum_drho_total = 0.0, cum_drho_band_total = 0.0;
    double e_side_J = 0, e_bottom_J = 0, e_top_J = 0;
    double last_P_side = 0, last_P_bottom = 0, last_P_top_conv = 0, last_P_top_rad = 0, last_P_bound = 0;
    std::ofstream monitor(cli.out_csv);
    // Doc208 A: last host power formula tag (must match applied BC; monitor is read-only).
    const char *last_power_host_formula = use_wall_T ? "wall-temperature" : "legacy-shell";
    monitor << "step,time_s,U_rms,Uz_rms,Ur_rms,Utheta_rms,K,T_mean,T_min,T_max,phase,dH_J,E_source_J,E_boundary_J,"
               "P_side_W,P_bottom_W,P_top_conv_W,P_top_rad_W,P_boundary_total_W,E_side_J,E_bottom_J,E_top_J,"
               "pi_mean,pi_rms_demean,pi_min,pi_max,rho_min,rho_max,circ_uz_proxy,circ_ur_proxy,U_max,"
               "top_uz_rms,top_uz_maxabs,n_outside_strict,n_outside_half_dp,z_max,"
               "rho_mean,mass_sum,sum_rhoV,sum_rhoV_over_mass,vol_sum,drho_dt_mean,drho_dt_rms,drho_dt_band_contrib,"
               "speed_mean,cum_drho_dissipative,cum_drho_total,cum_drho_band_total,cum_drho_continuity,"
               "thermal_wall_model,power_host_formula,power_sample_note\n";
    auto write_monitor_row = [&](size_t step, const char *phase) {
        const FlowStats s = hostFlowStats(particles, axis_center, R, H, z_bottom, cli.rho0, cli.c0);
        const PiStats pis = hostPiStatsFromDensity(particles, cli.rho0, cli.c0);
        const BoundaryBandStats bb = hostBoundaryBandStats(particles, axis_center, R, H, z_bottom, dp);
        const double dh = thermal_enthalpy_J() - h_init;
        const DriftSnapshot ds = hostDriftSnapshot(particles, axis_center, R, H, z_bottom, dp);
        // Doc208 A: sample host powers with the SAME formula as the applied BC (do not always use legacy).
        // Sampling phase note: monitor is after thermal substeps for "end"/periodic rows; init is pre-loop.
        const char *power_sample_note = "host_audit_read_only_matches_applied_bc";
        if (bc_active)
        {
            Real ps = 0, pb = 0, pfc = 0, pfr = 0, ptot = 0;
            last_power_host_formula = hostOphelieThermalAppliedBcHeatLossPowers(
                particles, thermal_opts, dp, use_wall_T, ps, pb, pfc, pfr, ptot);
            last_P_side = static_cast<double>(ps);
            last_P_bottom = static_cast<double>(pb);
            last_P_top_conv = static_cast<double>(pfc);
            last_P_top_rad = static_cast<double>(pfr);
            last_P_bound = static_cast<double>(ptot);
        }
        monitor << std::setprecision(10) << step << "," << sv_physical_time->getValue() << "," << s.U_rms << ","
                << s.Uz_rms << "," << s.Ur_rms << "," << s.Utheta_rms << "," << s.K << "," << s.T_mean << ","
                << s.T_min << "," << s.T_max << "," << phase << "," << dh << "," << ledger.e_source << ","
                << ledger.e_boundary << "," << last_P_side << "," << last_P_bottom << "," << last_P_top_conv << ","
                << last_P_top_rad << "," << last_P_bound << "," << e_side_J << "," << e_bottom_J << "," << e_top_J
                << "," << pis.mean << "," << pis.rms_demean << "," << pis.min << "," << pis.max << "," << s.rho_min
                << "," << s.rho_max << "," << s.circ_uz_proxy << "," << s.circ_ur_proxy << "," << s.U_max << ","
                << bb.top_uz_rms << "," << bb.top_uz_maxabs << "," << bb.n_outside_strict << ","
                << bb.n_outside_half_dp << "," << bb.z_max << "," << ds.rho_mean << "," << ds.mass_sum << ","
                << ds.sum_rhoV << "," << ds.sum_rhoV_over_mass << "," << ds.vol_sum << "," << ds.drho_dt_mean << ","
                << ds.drho_dt_rms << "," << ds.drho_dt_band_contrib << "," << ds.speed_mean << ","
                << cum_drho_diss << "," << cum_drho_total << "," << cum_drho_band_total << ","
                << (cum_drho_total - cum_drho_diss) << "," << (use_wall_T ? "wall-temperature" : "legacy-shell")
                << "," << last_power_host_formula << "," << power_sample_note << "\n";
        monitor.flush();
        return s;
    };

    std::cout << "[doc126] start mode=" << cli.mode << " dp=" << dp << " n=" << particles.TotalRealParticles()
              << " n_wall_closed=" << wall_closed.getBaseParticles().TotalRealParticles()
              << " n_wall_noslip=" << wall_noslip.getBaseParticles().TotalRealParticles() << " R=" << R
              << " H=" << H << " beta=" << beta_eff << " a_strat=" << a_eff << " c0=" << cli.c0
              << " mu=" << cli.mu << " settle=" << settle_eff << " end_time=" << cli.end_time
              << " force_rel_err=" << force_audit.rel_err << " Fz_sum=" << force_audit.f_z_sum
              << " buoy_a_rms=" << buoy_a_rms_init << " eos_init_res=" << eos_res_init << "/tol=" << eos_tol
              << std::endl;
    if (cli.doc201_b1_ablate != "off")
    {
        std::cout << "[doc201][b1] ablate=" << cli.doc201_b1_ablate
                  << " (diagnostic skip; positions remain live; not production)" << std::endl;
    }
    std::cout << "[doc126] contract perturbation_momentum=1 acoustic_dens_mode=1 ref_body_accel=0 "
                 "gravity_force_ck=0 boussinesq=1 density_path="
              << (cli.use_density_summation ? "internal_summation" : "continuity")
              << " kernel_correction=" << (cli.kernel_correction ? 1 : 0) << " top_bc=flat_free_slip_viscous_omit"
              << std::endl;

    size_t n_bad_dt = 0;
    size_t total_iter = 0;
    PressureSample pressure_sample;

    // One thermal advance of length dt (may subcycle): Laplace -> Q -> diffusion -> BC.
    // Doc177: refresh μ/cp/k each substep; credit E_source from Σ Q_i V_i (actual deposit), not only P_applied;
    // subcycle with cap=thermal_dt so flowing advection steps do not take one huge left-Riemann cp step.
    auto thermal_step = [&](Real dt) {
        if (!thermal_active || !(dt > TinyReal))
            return;
        const Real dt_cap = (cli.thermal_dt > TinyReal) ? cli.thermal_dt : dt;
        Real dt_remain = dt;
        while (dt_remain > TinyReal)
        {
            const Real dti = std::min(dt_cap, dt_remain);
            if (cli.variable_properties)
            {
                (void)hostUpdateFrenchThermalRheologyFromTemperature(
                    particles, material_laws, /*constant_properties=*/false, cli.t_ref, rheology_stats);
                hostClampParticleViscosityCap(particles, cli.mu_cap, rheology_stats);
            }
            if (conduction_active)
                laplace_temperature.exec();
            if (q_active)
            {
                if (cli.variable_properties)
                    apply_source_q_var_cp->exec(dti);
                else
                    apply_source_q->exec(dti);
                const double p_q = static_cast<double>(hostBodyJouleHeatPower(particles, kDoc131QField));
                ledger.e_source += p_q * static_cast<double>(dti);
            }
            if (conduction_active)
            {
                if (cli.variable_properties)
                    apply_diffusion_var_cp->exec(dti);
                else
                    apply_diffusion->exec(dti);
            }
            if (bc_active)
            {
                const double p_loss = use_wall_T ? static_cast<double>(boundary_loss_wall->exec())
                                                 : static_cast<double>(boundary_loss_legacy->exec());
                if (!ledger.p_loss_first_set)
                {
                    ledger.p_loss_first = p_loss;
                    ledger.p_loss_first_set = true;
                }
                ledger.e_boundary += p_loss * static_cast<double>(dti);
                {
                    // Face integral for E_side/E_bottom/E_top — same dispatch as monitor (Doc208 A).
                    // Taken BEFORE apply_bc so it matches F02 ledger timing (flux at pre-BC state).
                    Real ps = 0, pb = 0, pfc = 0, pfr = 0, ptot = 0;
                    last_power_host_formula = hostOphelieThermalAppliedBcHeatLossPowers(
                        particles, thermal_opts, dp, use_wall_T, ps, pb, pfc, pfr, ptot);
                    const double dtd = static_cast<double>(dti);
                    e_side_J += static_cast<double>(ps) * dtd;
                    e_bottom_J += static_cast<double>(pb) * dtd;
                    e_top_J += (static_cast<double>(pfc) + static_cast<double>(pfr)) * dtd;
                    last_P_side = static_cast<double>(ps);
                    last_P_bottom = static_cast<double>(pb);
                    last_P_top_conv = static_cast<double>(pfc);
                    last_P_top_rad = static_cast<double>(pfr);
                    last_P_bound = static_cast<double>(ptot);
                }
                if (use_wall_T)
                {
                    if (cli.variable_properties)
                        apply_wall_bc_var_cp->exec(dti);
                    else
                        apply_wall_bc->exec(dti);
                }
                else if (cli.variable_properties)
                    apply_natural_bc_var_cp->exec(dti);
                else
                    apply_natural_bc->exec(dti);
            }
            else if (legacy_path)
            {
                ledger.e_boundary += hostLegacyTopRadiationClampStep(particles, cli.rho0, cli.cp, dti, cli.emissivity,
                                                                     cli.t_ambient, z_top, dp, cli.t_ref);
            }
            ledger.t_thermal += static_cast<double>(dti);
            ledger.n_steps += 1;
            ledger.dt_min = std::min(ledger.dt_min, static_cast<double>(dti));
            ledger.dt_max = std::max(ledger.dt_max, static_cast<double>(dti));
            dt_remain -= dti;
        }
    };

    // Doc133: VTP trigger by physical-time interval (first call writes step `total_iter` at t>=next_vtp_time).
    auto maybe_write_vtp = [&](size_t step) {
        if (!write_vtp)
            return;
        if (sv_physical_time->getValue() + Real(1.0e-9) >= next_vtp_time)
        {
            write_vtp_now(step);
            next_vtp_time += cli.state_record_interval;
            while (next_vtp_time <= sv_physical_time->getValue())
                next_vtp_time += cli.state_record_interval;
        }
    };

    // Doc183: lightweight particle field CSV (scheduled times + optional post-EM tags).
    const std::string field_export_root =
        cli.field_export_dir.empty()
            ? ((cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/field_exports")
            : cli.field_export_dir;
    std::vector<char> field_export_done(cli.field_export_times.size(), 0);
    size_t n_field_exports = 0;
    auto export_field_now = [&](const std::string &tag, const std::string &phase_note) {
        if (cli.field_export_times.empty())
            return;
        const OphelieFrenchGlassMaterialLaws *laws_ptr =
            cli.variable_properties ? &material_laws : nullptr;
        const Real *mu_ptr = nullptr;
        if (cli.variable_properties || constant_mu_implicit)
        {
            syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
            mu_ptr = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        }
        const auto er = hostExportOphelieParticleFieldCsv(
            particles, field_export_root, sv_physical_time->getValue(), total_iter, tag, phase_note,
            /*sorted_flag=*/!cli.legacy_thermal_sort, laws_ptr, cli.mu_cap, cli.variable_properties, mu_ptr,
            cli.mu);
        if (!er.ok)
            std::cerr << "[doc183][field-export] FAIL tag=" << tag << " err=" << er.error << std::endl;
        else
        {
            ++n_field_exports;
            std::cout << "[doc183][field-export] ok tag=" << tag << " t=" << er.time_s << " n=" << er.n
                      << " U_max=" << er.u_max << " U_rms=" << er.u_rms << " csv=" << er.path_csv << std::endl;
        }
    };
    auto maybe_export_scheduled_fields = [&]() {
        if (cli.field_export_times.empty())
            return;
        const Real t_now = sv_physical_time->getValue();
        for (size_t k = 0; k < cli.field_export_times.size(); ++k)
        {
            if (field_export_done[k])
                continue;
            const Real t_target = cli.field_export_times[k];
            if (t_now + Real(1.0e-9) >= t_target)
            {
                std::ostringstream tag;
                tag << "sched_" << k;
                export_field_now(tag.str(), "scheduled_time_ge_target");
                field_export_done[k] = 1;
            }
        }
    };

    // Doc133/134 E: re-evaluate fixed spatial Q shape from CURRENT positions (P conserved each time).
    size_t n_q_shape_reapply = 0;
    auto reapply_q_shape = [&]() {
        if (!q_active || cli.q_shape == "uniform")
            return;
        PowerReport tmp;
        (void)hostApplyPowerContract(particles, cli.power_mode, static_cast<double>(cli.q_power),
                                     static_cast<double>(cli.rho0), v_cyl, tmp, cli.q_shape, axis_center,
                                     static_cast<double>(R), rz_q_ptr, em_q_ptr);
        power.P_applied = tmp.P_applied;
        power.N_Q = tmp.N_Q;
        power.norm_coeff = tmp.norm_coeff;
        power.sum_V_Q = tmp.sum_V_Q;
        power.q_min = tmp.q_min;
        power.q_max = tmp.q_max;
        power.gap_rel = tmp.gap_rel;
        ++n_q_shape_reapply;
    };

    // Doc151/Doc156: file-driven subprocess EM refresh (keeps fluid state; reloads Eulerian Q table).
    size_t n_em_updates = 0;
    bool em_coupling_ok = true;
    const long max_em_updates_eff =
        !cli.em_update_times.empty() ? static_cast<long>(cli.em_update_times.size()) : cli.max_em_updates;
    bool em_coupling_armed = (cli.em_update_mode == "subprocess" && max_em_updates_eff > 0);
    Real next_em_due =
        cli.em_update_times.empty() ? (cli.settle_time + cli.em_update_interval)
                                    : cli.em_update_times.front();
    double last_sigma_mean = 16.0;
    double last_delta_sigma = 0.0;
    double last_shape_L2 = -1.0;
    std::string last_em_Q_table_csv = cli.q_from_rz_csv; /**< helper shape-L2 vs previous table */
    std::string last_em_Q_csv =
        !cli.q_from_em_particles_csv.empty() ? cli.q_from_em_particles_csv : cli.q_from_rz_csv;
    std::string em_failure_reason;
    std::vector<std::string> em_update_meta_paths;
    std::vector<double> em_update_scheduled_times;
    std::vector<double> em_update_actual_times;
    size_t n_fluid_steps_after_last_em = 0;
    size_t n_thermal_steps_after_last_em = 0;
    auto write_em_coupling_state = [&](size_t step_tag) {
        // Doc171: atomic write + per-step sidecar; restore face ledger + mass baseline.
        const std::string ckpt_dir =
            (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/checkpoints";
        std::filesystem::create_directories(ckpt_dir);
        const std::string step_path =
            ckpt_dir + "/em_coupling_state_step_" + std::to_string(step_tag) + ".json";
        const std::string tmp_path = em_state_path + ".tmp";
        const double mass_now = hostTotalMass(particles);
        auto relativize_under_output = [&](const std::string &p) -> std::string {
            if (p.empty() || cli.output_dir.empty())
                return p;
            std::error_code ec;
            const auto out_abs = std::filesystem::weakly_canonical(cli.output_dir, ec);
            if (ec)
                return p;
            const auto p_abs = std::filesystem::weakly_canonical(p, ec);
            if (ec)
                return p;
            const std::string outs = out_abs.string();
            const std::string ps = p_abs.string();
            if (ps.rfind(outs, 0) != 0)
                return p;
            std::string rel = ps.substr(outs.size());
            while (!rel.empty() && (rel[0] == '/' || rel[0] == '\\'))
                rel.erase(rel.begin());
            return rel.empty() ? p : rel;
        };
        const std::string q_for_state = relativize_under_output(last_em_Q_csv);
        auto dump_body = [&](std::ostream &out) {
            out << std::setprecision(17);
            out << "{\n";
            out << "  \"schema\": \"doc171_em_coupling_state_v2\",\n";
            out << "  \"restart_step\": " << step_tag << ",\n";
            out << "  \"physical_time_s\": " << sv_physical_time->getValue() << ",\n";
            out << "  \"em_schedule\": \"" << cli.em_schedule << "\",\n";
            out << "  \"em_update_origin_s\": " << cli.em_update_origin << ",\n";
            out << "  \"em_update_interval_s\": " << cli.em_update_interval << ",\n";
            out << "  \"n_em_updates\": " << n_em_updates << ",\n";
            out << "  \"max_em_updates\": " << max_em_updates_eff << ",\n";
            out << "  \"next_em_due_s\": " << next_em_due << ",\n";
            out << "  \"em_coupling_ok\": " << (em_coupling_ok ? "true" : "false") << ",\n";
            out << "  \"em_feedback_capped\": "
                << ((static_cast<long>(n_em_updates) >= max_em_updates_eff) ? "true" : "false") << ",\n";
            out << "  \"em_initial_consistent\": " << (cli.em_initial_consistent ? "true" : "false") << ",\n";
            out << "  \"last_em_Q_csv\": \"" << q_for_state << "\",\n";
            out << "  \"em_q_support_radius_m\": " << cli.em_q_support_radius << ",\n";
            out << "  \"em_particles_csv\": \"" << cli.em_particles_csv << "\",\n";
            out << "  \"material_curve_dir\": \"" << cli.material_curve_dir << "\",\n";
            out << "  \"material_mode\": \"" << cli.material_mode << "\",\n";
            out << "  \"ledger_e_source_J\": " << ledger.e_source << ",\n";
            out << "  \"ledger_e_boundary_J\": " << ledger.e_boundary << ",\n";
            out << "  \"ledger_e_side_J\": " << e_side_J << ",\n";
            out << "  \"ledger_e_bottom_J\": " << e_bottom_J << ",\n";
            out << "  \"ledger_e_top_J\": " << e_top_J << ",\n";
            out << "  \"ledger_t_thermal_s\": " << ledger.t_thermal << ",\n";
            out << "  \"ledger_n_steps\": " << ledger.n_steps << ",\n";
            out << "  \"h_init_J\": " << h_init << ",\n";
            out << "  \"mass_sum_init_kg\": " << mass_sum_init << ",\n";
            out << "  \"mass_sum_current_kg\": " << mass_now << ",\n";
            out << "  \"total_iter\": " << total_iter << "\n";
            out << "}\n";
        };
        {
            std::ofstream out(tmp_path);
            if (!out)
            {
                std::cerr << "[doc171][checkpoint] FAIL open tmp " << tmp_path << std::endl;
                return;
            }
            dump_body(out);
        }
        std::error_code ec;
        std::filesystem::rename(tmp_path, em_state_path, ec);
        if (ec)
        {
            std::filesystem::copy_file(tmp_path, em_state_path,
                                       std::filesystem::copy_options::overwrite_existing, ec);
            std::filesystem::remove(tmp_path, ec);
        }
        {
            std::ofstream out(step_path);
            dump_body(out);
        }
        std::cout << "[doc171][checkpoint] wrote " << em_state_path << " + " << step_path
                  << " step=" << step_tag << " t=" << sv_physical_time->getValue() << std::endl;
    };
    auto write_checkpoint_now = [&](size_t step_tag) {
        restart_io.writeToFile(step_tag);
        write_em_coupling_state(step_tag);
        last_checkpoint_step = step_tag;
    };
    auto load_em_coupling_state = [&]() -> bool {
        // Doc171: prefer step-bound sidecar so Restart_N cannot pair with a newer tip state.
        const std::string step_bound =
            (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/checkpoints/em_coupling_state_step_" +
            std::to_string(cli.restart_step) + ".json";
        std::string load_path = em_state_path;
        if (std::filesystem::exists(step_bound))
            load_path = step_bound;
        std::ifstream in(load_path);
        if (!in)
        {
            std::cerr << "[doc171][restart] missing em state: " << load_path << std::endl;
            return false;
        }
        std::string txt((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::cout << "[doc171][restart] loading EM sidecar " << load_path << std::endl;
        auto find_num = [&](const char *key, double &out) -> bool {
            const std::string pat = std::string("\"") + key + "\"";
            size_t p0 = txt.find(pat);
            if (p0 == std::string::npos)
                return false;
            p0 = txt.find(':', p0 + pat.size());
            if (p0 == std::string::npos)
                return false;
            ++p0;
            while (p0 < txt.size() && std::isspace(static_cast<unsigned char>(txt[p0])))
                ++p0;
            size_t end = p0;
            while (end < txt.size() &&
                   (std::isdigit(static_cast<unsigned char>(txt[end])) || txt[end] == '+' || txt[end] == '-' ||
                    txt[end] == '.' || txt[end] == 'e' || txt[end] == 'E'))
                ++end;
            return end > p0 && parseRealStrict(txt.substr(p0, end - p0), out);
        };
        auto find_str = [&](const char *key, std::string &out) -> bool {
            const std::string pat = std::string("\"") + key + "\"";
            size_t p0 = txt.find(pat);
            if (p0 == std::string::npos)
                return false;
            p0 = txt.find(':', p0 + pat.size());
            if (p0 == std::string::npos)
                return false;
            p0 = txt.find('"', p0 + 1);
            if (p0 == std::string::npos)
                return false;
            size_t p1 = txt.find('"', p0 + 1);
            if (p1 == std::string::npos)
                return false;
            out = txt.substr(p0 + 1, p1 - p0 - 1);
            return true;
        };
        double n_em = 0, max_em = 0, next_due = 0, e_src = 0, e_bnd = 0, t_th = 0, n_steps = 0, tot = 0;
        std::string sched, qpath;
        if (!find_num("n_em_updates", n_em) || !find_num("next_em_due_s", next_due) ||
            !find_str("last_em_Q_csv", qpath))
        {
            std::cerr << "[doc164][restart] em_coupling_state.json incomplete" << std::endl;
            return false;
        }
        (void)find_str("em_schedule", sched);
        (void)find_num("max_em_updates", max_em);
        (void)find_num("ledger_e_source_J", e_src);
        (void)find_num("ledger_e_boundary_J", e_bnd);
        (void)find_num("ledger_t_thermal_s", t_th);
        (void)find_num("ledger_n_steps", n_steps);
        (void)find_num("total_iter", tot);
        double h0 = 0;
        if (find_num("h_init_J", h0))
            h_init = h0;
        double e_side = 0, e_bot = 0, e_top = 0, mass0 = 0;
        if (find_num("ledger_e_side_J", e_side))
            e_side_J = e_side;
        if (find_num("ledger_e_bottom_J", e_bot))
            e_bottom_J = e_bot;
        if (find_num("ledger_e_top_J", e_top))
            e_top_J = e_top;
        if (find_num("mass_sum_init_kg", mass0) && mass0 > 0)
            mass_sum_init = mass0;
        double step_tag = -1;
        if (find_num("restart_step", step_tag) &&
            static_cast<long>(step_tag + 0.5) != cli.restart_step)
        {
            std::cerr << "[doc171][restart] sidecar restart_step=" << step_tag
                      << " != --restart-step=" << cli.restart_step
                      << " (refusing mismatched identity)" << std::endl;
            return false;
        }
        n_em_updates = static_cast<size_t>(n_em + 0.5);
        next_em_due = static_cast<Real>(next_due);
        std::string q_resolved = qpath;
        if (!qpath.empty() && !std::filesystem::path(qpath).is_absolute())
            q_resolved = (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/" + qpath;
        if (!qpath.empty() && !std::filesystem::exists(q_resolved))
        {
            std::cerr << "[doc171][restart] last_em_Q_csv missing after resolve: " << q_resolved
                      << " (from '" << qpath << "')" << std::endl;
            return false;
        }
        last_em_Q_csv = q_resolved;
        ledger.e_source = e_src;
        ledger.e_boundary = e_bnd;
        ledger.t_thermal = t_th;
        ledger.n_steps = static_cast<size_t>(n_steps + 0.5);
        total_iter = static_cast<size_t>(tot + 0.5);
        if (!q_resolved.empty() && cli.q_shape == "em_particles_csv")
        {
            if (!em_q_cloud.load(q_resolved, static_cast<double>(cli.em_q_support_radius)))
            {
                std::cerr << "[doc171][restart] failed to reload Q cloud: " << q_resolved << std::endl;
                return false;
            }
            em_q_ptr = &em_q_cloud;
        }
        std::cout << "[doc164][restart] restored EM state n_em=" << n_em_updates << " next_due=" << next_em_due
                  << " Q=" << last_em_Q_csv << std::endl;
        return true;
    };
    auto jsonFindBool = [](const std::string &txt, const char *key, bool &out) -> bool {
        const std::string pat = std::string("\"") + key + "\"";
        size_t p = txt.find(pat);
        if (p == std::string::npos)
            return false;
        p = txt.find(':', p + pat.size());
        if (p == std::string::npos)
            return false;
        while (p + 1 < txt.size() && std::isspace(static_cast<unsigned char>(txt[p + 1])))
            ++p;
        if (txt.compare(p + 1, 4, "true") == 0)
        {
            out = true;
            return true;
        }
        if (txt.compare(p + 1, 5, "false") == 0)
        {
            out = false;
            return true;
        }
        return false;
    };
    auto jsonFindNumber = [](const std::string &txt, const char *key, double &out) -> bool {
        const std::string pat = std::string("\"") + key + "\"";
        size_t p = txt.find(pat);
        if (p == std::string::npos)
            return false;
        p = txt.find(':', p + pat.size());
        if (p == std::string::npos)
            return false;
        ++p;
        while (p < txt.size() && std::isspace(static_cast<unsigned char>(txt[p])))
            ++p;
        size_t end = p;
        while (end < txt.size() &&
               (std::isdigit(static_cast<unsigned char>(txt[end])) || txt[end] == '+' || txt[end] == '-' ||
                txt[end] == '.' || txt[end] == 'e' || txt[end] == 'E'))
            ++end;
        if (end == p)
            return false;
        return parseRealStrict(txt.substr(p, end - p), out);
    };
    auto spawn_python_helper = [](const std::string &helper, const std::string &req_path) -> int {
        const std::string req_arg = "--request=" + req_path;
        char *argv[] = {const_cast<char *>("python3"), const_cast<char *>(helper.c_str()),
                        const_cast<char *>(req_arg.c_str()), nullptr};
        pid_t pid = 0;
        const int spawn_rc = posix_spawnp(&pid, "python3", nullptr, nullptr, argv, environ);
        if (spawn_rc != 0)
            return -1000 - spawn_rc;
        int status = 0;
        if (waitpid(pid, &status, 0) < 0)
            return -1;
        if (WIFSIGNALED(status))
            return -2000 - WTERMSIG(status);
        if (!WIFEXITED(status))
            return -2;
        return WEXITSTATUS(status);
    };
    auto maybe_em_subprocess_update = [&](const char *phase, bool initial_u0) {
        if (!em_coupling_armed || !em_coupling_ok)
            return;
        if (!initial_u0 && std::string(phase) != "observe")
            return;
        const Real t_now = sv_physical_time->getValue();
        if (!initial_u0 && static_cast<long>(n_em_updates) >= max_em_updates_eff)
        {
            // Doc164: cap reached → feedback stopped (do not claim continuous coupling)
            return;
        }
        if (!initial_u0)
        {
            if (cli.em_schedule == "periodic")
            {
                // next_event = origin + (n_em_updates+1)*interval  (1-based event index)
                next_em_due = cli.em_update_origin +
                              static_cast<Real>(static_cast<long>(n_em_updates) + 1) * cli.em_update_interval;
            }
            else if (!cli.em_update_times.empty())
                next_em_due = cli.em_update_times[n_em_updates];
            if (t_now + Real(1.0e-9) < next_em_due)
                return;
        }

        const std::string work_root =
            (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/em_updates";
        std::filesystem::create_directories(work_root);
        // Doc171: u0 initial consistent solve does not consume scheduled event indices.
        const size_t idx = initial_u0 ? 0 : (n_em_updates + 1);
        const std::string udir = work_root + (initial_u0 ? "/u0" : ("/u" + std::to_string(idx)));
        if (std::filesystem::exists(udir))
            std::filesystem::remove_all(udir);
        std::filesystem::create_directories(udir);

        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, kOphelieTemperatureField);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<UnsignedInt>(particles, "OriginalID");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
        const size_t nn = particles.TotalRealParticles();
        const std::string fluid_csv = udir + "/fluid_T.csv";
        {
            std::ofstream out(fluid_csv);
            out << "OriginalID,x_m,y_m,z_m,T_K,vol_m3\n";
            out.setf(std::ios::scientific);
            out.precision(10);
            for (size_t i = 0; i < nn; ++i)
            {
                out << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << T[i] << ","
                    << vol[i] << "\n";
            }
        }

        // Doc156 formal spatial_field: never apply +T. volume_mean legacy may still use last-update diagnostic.
        const bool last_update = (static_cast<long>(idx) >= max_em_updates_eff);
        const Real T_pert =
            (cli.em_sigma_reduction == "volume_mean" && last_update &&
             std::abs(cli.em_diag_T_perturb_K) > TinyReal)
                ? cli.em_diag_T_perturb_K
                : Real(0);
        const std::string pert_label =
            (std::abs(T_pert) > TinyReal) ? "diagnostic_T_perturb_last_update" : "none";
        const bool allow_T_perturb = (std::abs(T_pert) > TinyReal);
        const double scheduled = initial_u0 ? 0.0 : static_cast<double>(next_em_due);
        const std::string req_id =
            initial_u0 ? std::string("u0_initial")
                       : (std::string("u") + std::to_string(idx) + "_t" + std::to_string(scheduled));

        const std::string req_path = udir + "/request.json";
        {
            std::ofstream req(req_path);
            req.setf(std::ios::scientific);
            req.precision(10);
            req << "{\n"
                << "  \"fluid_csv\": \"" << fluid_csv << "\",\n"
                << "  \"out_dir\": \"" << udir << "\",\n"
                << "  \"update_index\": " << idx << ",\n"
                << "  \"request_id\": \"" << req_id << "\",\n"
                << "  \"picard_bin\": \"" << cli.em_solver_bin << "\",\n"
                << "  \"I_A\": " << cli.em_I_A << ",\n"
                << "  \"dp_em\": " << cli.em_dp << ",\n"
                << "  \"j_tol\": " << cli.em_j_tol << ",\n"
                << "  \"T_perturb_K\": " << T_pert << ",\n"
                << "  \"allow_T_perturb\": " << (allow_T_perturb ? "true" : "false") << ",\n"
                << "  \"T_perturb_label\": \"" << pert_label << "\",\n"
                << "  \"previous_sigma_mean\": " << last_sigma_mean << ",\n"
                << "  \"previous_q_csv\": \"" << last_em_Q_table_csv << "\",\n"
                << "  \"q_nr\": 28, \"q_nz\": 14,\n"
                << "  \"P_target\": " << cli.q_power << ",\n"
                << "  \"sigma_reduction\": \"" << cli.em_sigma_reduction << "\",\n"
                << "  \"em_q_delivery\": \"" << cli.em_q_delivery << "\",\n"
                << "  \"em_particles_csv\": \"" << cli.em_particles_csv << "\",\n"
                << "  \"clip_sigma_max\": null\n"
                << "}\n";
        }

        std::cout << "[doc156][em-update] t=" << t_now << " scheduled=" << next_em_due << " idx=" << idx
                  << " sigma_reduction=" << cli.em_sigma_reduction << " q_delivery=" << cli.em_q_delivery
                  << " helper=" << cli.em_update_helper << std::endl
                  << std::flush;
        const std::string meta_path = udir + "/meta.json";
        const std::string q_path = udir + "/Q_rz.csv";
        const std::string particles_path = udir + "/particles_rz.csv";
        const bool need_particles = (cli.em_q_delivery == "em_particles");
        const bool need_table = !need_particles;
        int rc = 0;
        if (!cli.em_replay_dir.empty())
        {
            // Doc169 verification: copy frozen uN artifacts; do not spawn Picard.
            const std::string src_dir = cli.em_replay_dir + "/u" + std::to_string(idx);
            const std::string src_part = src_dir + "/particles_rz.csv";
            const std::string src_meta = src_dir + "/meta.json";
            const std::string src_q = src_dir + "/Q_rz.csv";
            if (need_particles && !std::filesystem::exists(src_part))
            {
                em_coupling_ok = false;
                em_failure_reason = "em_replay_missing_particles:" + src_part;
                std::cerr << "[doc169][em-replay] FAIL " << em_failure_reason << std::endl;
                return;
            }
            if (need_table && !std::filesystem::exists(src_q))
            {
                em_coupling_ok = false;
                em_failure_reason = "em_replay_missing_Q_table:" + src_q;
                std::cerr << "[doc169][em-replay] FAIL " << em_failure_reason << std::endl;
                return;
            }
            std::error_code ec;
            if (std::filesystem::exists(src_part))
                std::filesystem::copy_file(src_part, particles_path,
                                           std::filesystem::copy_options::overwrite_existing, ec);
            if (std::filesystem::exists(src_q))
                std::filesystem::copy_file(src_q, q_path, std::filesystem::copy_options::overwrite_existing, ec);
            if (std::filesystem::exists(src_meta))
                std::filesystem::copy_file(src_meta, meta_path, std::filesystem::copy_options::overwrite_existing, ec);
            else
            {
                std::ofstream m(meta_path);
                m << "{\n"
                  << "  \"ok\": true,\n"
                  << "  \"em_result_ok\": true,\n"
                  << "  \"picard_converged\": true,\n"
                  << "  \"converged\": true,\n"
                  << "  \"sigma_mean\": " << last_sigma_mean << ",\n"
                  << "  \"delta_sigma_mean\": 0.0,\n"
                  << "  \"shape_L2_vs_previous_common_power\": -1.0,\n"
                  << "  \"replay\": true,\n"
                  << "  \"replay_src\": \"" << src_dir << "\"\n"
                  << "}\n";
            }
            rc = 0;
            std::cout << "[doc169][em-replay] idx=" << idx << " from " << src_dir << std::endl << std::flush;
        }
        else
        {
            if (!std::filesystem::exists(cli.em_update_helper))
            {
                em_coupling_ok = false;
                em_failure_reason = "helper_missing:" + cli.em_update_helper;
                std::cerr << "[doc156][em-update] FAIL " << em_failure_reason << std::endl;
                return;
            }
            rc = spawn_python_helper(cli.em_update_helper, req_path);
        }
        if (rc != 0 || !std::filesystem::exists(meta_path) ||
            (need_table && !std::filesystem::exists(q_path)) ||
            (need_particles && !std::filesystem::exists(particles_path)))
        {
            em_coupling_ok = false;
            em_failure_reason = cli.em_replay_dir.empty()
                                    ? ("subprocess_failed_rc=" + std::to_string(rc))
                                    : ("em_replay_incomplete_rc=" + std::to_string(rc));
            std::cerr << "[doc156][em-update] FAIL " << em_failure_reason
                      << " (refusing silent degraded_previous_q as PASS)" << std::endl;
            return;
        }
        std::ifstream mf(meta_path);
        std::stringstream ss;
        ss << mf.rdbuf();
        const std::string mtxt = ss.str();
        bool meta_ok = false, meta_conv = false, meta_picard = false, meta_em_ok = false;
        if (!jsonFindBool(mtxt, "ok", meta_ok) || !meta_ok)
        {
            em_coupling_ok = false;
            em_failure_reason = "meta_ok_false_or_missing";
            std::cerr << "[doc156][em-update] FAIL " << em_failure_reason << std::endl;
            return;
        }
        // Doc156 spatial helper writes em_result_ok; Doc151 volume-mean helper may only have picard_converged.
        const bool has_em_result_ok = jsonFindBool(mtxt, "em_result_ok", meta_em_ok);
        const bool has_picard = jsonFindBool(mtxt, "picard_converged", meta_picard);
        const bool has_conv = jsonFindBool(mtxt, "converged", meta_conv);
        if (cli.em_sigma_reduction == "spatial_field")
        {
            if (!has_em_result_ok || !meta_em_ok || !has_picard || !meta_picard)
            {
                em_coupling_ok = false;
                em_failure_reason = "strict_em_result_gate_failed";
                std::cerr << "[doc156][em-update] FAIL " << em_failure_reason << std::endl;
                return;
            }
        }
        else if (!has_picard || !meta_picard)
        {
            em_coupling_ok = false;
            em_failure_reason = "picard_converged_false_or_missing";
            std::cerr << "[doc156][em-update] FAIL " << em_failure_reason << std::endl;
            return;
        }
        (void)has_conv;
        (void)meta_conv;
        if (need_particles)
        {
            if (!em_q_cloud.load(particles_path, static_cast<double>(cli.em_q_support_radius)))
            {
                em_coupling_ok = false;
                em_failure_reason = "reload_em_particles_Q_failed";
                std::cerr << "[doc156][em-update] FAIL reload " << particles_path << std::endl;
                return;
            }
            em_q_ptr = &em_q_cloud;
            last_em_Q_csv = particles_path;
            if (std::filesystem::exists(q_path))
                last_em_Q_table_csv = q_path;
        }
        else
        {
            if (!rz_q_grid.load(q_path))
            {
                em_coupling_ok = false;
                em_failure_reason = "reload_Q_failed";
                std::cerr << "[doc156][em-update] FAIL reload " << q_path << std::endl;
                return;
            }
            rz_q_ptr = &rz_q_grid;
            last_em_Q_csv = q_path;
            last_em_Q_table_csv = q_path;
        }
        reapply_q_shape();
        {
            double p_em = std::numeric_limits<double>::quiet_NaN();
            double p_part = std::numeric_limits<double>::quiet_NaN();
            double p_tab = std::numeric_limits<double>::quiet_NaN();
            double scale = std::numeric_limits<double>::quiet_NaN();
            jsonFindNumber(mtxt, "P_joule_W", p_em);
            jsonFindNumber(mtxt, "P_particles_W", p_part);
            jsonFindNumber(mtxt, "P_table_raw_W", p_tab);
            jsonFindNumber(mtxt, "scale_table_raw_to_target", scale);
            const std::string stage = initial_u0 ? "em_u0" : ("em_u" + std::to_string(idx));
            note_power_stage(stage, static_cast<double>(t_now), p_em, p_part, p_tab, scale);
        }
        em_update_meta_paths.push_back(meta_path);
        em_update_scheduled_times.push_back(scheduled);
        em_update_actual_times.push_back(static_cast<double>(t_now));
        if (!initial_u0)
        {
            ++n_em_updates;
            if (cli.write_em_checkpoints)
                write_checkpoint_now(total_iter);
            n_fluid_steps_after_last_em = 0;
            n_thermal_steps_after_last_em = 0;
            if (cli.em_update_times.empty())
                next_em_due = t_now + cli.em_update_interval;
            else if (n_em_updates < cli.em_update_times.size())
                next_em_due = cli.em_update_times[n_em_updates];
            else
                next_em_due = std::numeric_limits<Real>::max();
        }
        else
        {
            std::cout << "[doc171][em-initial] u0 Q applied; n_em_updates remains " << n_em_updates
                      << " next_em_due=" << next_em_due << std::endl;
        }
        if (!jsonFindNumber(mtxt, "sigma_mean", last_sigma_mean))
            jsonFindNumber(mtxt, "sigma_mean_vol", last_sigma_mean);
        jsonFindNumber(mtxt, "delta_sigma_mean", last_delta_sigma);
        jsonFindNumber(mtxt, "shape_L2_vs_previous_common_power", last_shape_L2);
        std::cout << "[doc156][em-update] ok idx=" << idx << " scheduled=" << scheduled << " actual_t=" << t_now
                  << " sigma_mean=" << last_sigma_mean << " delta_sigma=" << last_delta_sigma
                  << " shape_L2_vs_prev=" << last_shape_L2 << " Q=" << q_path << std::endl
                  << std::flush;
        // Doc183: post-EM applied Q + current fluid state (recommended export phase).
        {
            std::ostringstream tag;
            tag << (initial_u0 ? "em_u0" : ("em_u" + std::to_string(idx)));
            export_field_now(tag.str(), initial_u0 ? "after_em_initial_u0" : "after_em_periodic_update");
            maybe_export_scheduled_fields();
        }
    };

    // Doc134: --cref-buoyancy (production) / Doc133 probe: match wall dens C_ref to Boussinesq C_ext.
    auto probe_zero_velocity = [&]() {
        syncVariableToHost<Vecd>(particles, "Velocity");
        Vecd *v = particles.getVariableDataByName<Vecd>("Velocity");
        const size_t nn = particles.TotalRealParticles();
        for (size_t i = 0; i < nn; ++i)
            v[i] = Vecd::Zero();
        syncVariableToDevice<Vecd>(particles, "Velocity");
    };
    // Doc136: half-step ledger matching AcousticStep writes:
    //   expected Δρ_mean = 0.5*dt*mean(rate_before) + 0.5*dt*mean(rate_after_2nd)
    // Legacy Doc133 used mean(rate)*full_dt after each half (over-counts; not the integrator).
    double probe_rate_mean_before = 0.0;
    auto probe_capture_before = [&]() {
        const DriftSnapshot ds = hostDriftSnapshot(particles, axis_center, R, H, z_bottom, dp);
        probe_rate_mean_before = ds.drho_dt_mean;
    };
    auto probe_accumulate_substep = [&](Real dt) {
        const DriftSnapshot ds = hostDriftSnapshot(particles, axis_center, R, H, z_bottom, dp);
        const double h = 0.5 * static_cast<double>(dt);
        // Dissipative rate scale after 1st interact ≈ current mean before continuity add is not stored here;
        // use (rate_after2 - continuity_proxy). Approximate dissipative half-weight with rate_before
        // (previous step's full rate) + half of post-2nd if continuity small — keep explicit split:
        cum_drho_diss += h * probe_rate_mean_before;          // 1st write uses rate_previous
        cum_drho_total += h * probe_rate_mean_before + h * ds.drho_dt_mean; // both writes
        cum_drho_band_total += h * ds.drho_dt_band_contrib;   // band of rate_after2 (2nd write weight)
    };

    // Doc136 P0 density half-step audit
    Doc136DensityAudit density_audit;
    density_audit.enabled = cli.density_increment_audit;
    density_audit.stride = cli.density_audit_stride;
    density_audit.max_rows = cli.density_audit_max_rows;
    std::vector<Real> audit_rho_b, audit_rate_b, audit_rho_1, audit_rate_1, audit_rho_2, audit_rate_2;
    if (density_audit.enabled)
    {
        const std::string adir = cli.output_dir.empty() ? std::string(".") : cli.output_dir;
        std::filesystem::create_directories(adir);
        density_audit.open(adir);
        density_audit.selectTracked(particles, axis_center, R, H, z_bottom, dp);
        std::cout << "[doc136][density-audit] enabled stride=" << density_audit.stride
                  << " max_rows=" << density_audit.max_rows << " tracked=" << density_audit.tracked.size()
                  << " global_csv=" << density_audit.global_path << std::endl;
        for (const auto &t : density_audit.tracked)
            std::cout << "[doc136][density-audit] zone=" << t.zone << " OriginalID=" << t.original_id << std::endl;
    }

    // Frozen-flow advance (T0-T4 / --freeze-flow): positions, velocity, density, pressure untouched.
    auto advance_frozen = [&](Real t_target, const char *phase) {
        size_t local_iter = 0;
        const Real eps_t = Real(1.0e-6) * std::max(Real(1), t_target);
        // Doc183/Doc206: scheduled field exports must fire during freeze thermal advance, not only at end flush.
        maybe_export_scheduled_fields();
        while (t_target - sv_physical_time->getValue() > eps_t)
        {
            const Real dt = std::min(cli.thermal_dt, t_target - sv_physical_time->getValue());
            reapply_q_shape();
            thermal_step(dt);
            sv_physical_time->incrementValue(dt);
            ++local_iter;
            ++total_iter;
            maybe_write_vtp(total_iter);
            maybe_export_scheduled_fields();
            if (local_iter % 20 == 0)
                (void)write_monitor_row(total_iter, phase);
        }
    };

    // Doc131 D2 (--r2-fixed-spatial-T) and Doc203 F1 (--fixed-spatial-T + init-T-rz-csv):
    // T is re-applied from the CURRENT particle position so the driver is a time-independent spatial field.
    // Default (off) = Lagrangian carry of the initial T. Only T/DeltaT are rewritten.
    // --fixed-T-after-settle: stay inactive (Tref) until settle ends, then arm the driver.
    size_t n_fixed_T_reapply = 0;
    bool fixed_T_live = !cli.fixed_T_after_settle;
    auto reapply_fixed_spatial_T = [&]() {
        if (!fixed_T_live)
            return;
        if (cli.fixed_spatial_T && !fixed_T_rz_nodes.empty())
        {
            std::string terr;
            if (!hostApplyRzTemperatureNodes(particles, axis_center, cli.t_ref, fixed_T_rz_nodes, cli.dp * Real(2),
                                             terr, /*quiet=*/true, cli.fixed_spatial_T_demean))
            {
                std::cerr << "[doc203][f1] fixed-spatial-T reapply failed: " << terr << std::endl;
            }
            ++n_fixed_T_reapply;
            return;
        }
        if (!(is_R2 && cli.r2_fixed_spatial_T))
            return;
        (void)hostPrescribeRadialTemperature(particles, axis_center, R, cli.t_ref, r2_dT_eff, !cli.r2_keep_mean,
                                             &r2_s_mean);
        ++n_fixed_T_reapply;
    };
    auto activate_fixed_T_after_settle = [&]() {
        if (!cli.fixed_T_after_settle || fixed_T_live)
            return;
        if (cli.fixed_spatial_T && !fixed_T_rz_nodes.empty())
        {
            std::string terr;
            if (!hostApplyRzTemperatureNodes(particles, axis_center, cli.t_ref, fixed_T_rz_nodes, cli.dp * Real(2),
                                             terr, /*quiet=*/false, cli.fixed_spatial_T_demean))
            {
                std::cerr << "[doc203][f1] fixed-T activate after settle failed: " << terr << std::endl;
            }
            else
                std::cout << "[doc203][f1] fixed-spatial-T CSV activated after settle" << std::endl;
        }
        else if (is_R2 && cli.r2_fixed_spatial_T)
        {
            r2_s_mean = static_cast<double>(hostPrescribeRadialTemperature(
                particles, axis_center, R, cli.t_ref, r2_dT_eff, !cli.r2_keep_mean));
            std::cout << "[doc203][f1] R2 fixed-spatial-T activated after settle s_mean=" << r2_s_mean
                      << std::endl;
        }
        fixed_T_live = true;
        ++n_fixed_T_reapply;
    };

    int doc189_gate_n = 0;
    double doc189_worst_L2 = 0.0;
    double doc189_worst_Enet = 0.0;
    double doc189_worst_abs = 0.0;
    bool doc189_algebraic_ok = true;
    bool doc189_dynamic_ok = true;
    OphelieBuoyancyFieldGate doc189_last_gate;
    auto advance_to = [&](Real t_target, const char *phase, bool sample_pressure) {
        size_t local_iter = 0;
        while (sv_physical_time->getValue() < t_target)
        {
            const Real t_before = sv_physical_time->getValue();
            if (cli.use_density_summation)
            {
                density_summation.exec();
                density_reg.exec(); // Internal regularization (closed domain)
            }
            advection_setup.exec();
            if (cli.kernel_correction)
            {
                fluid_linear_correction_matrix.exec();
            }
            const Real dt_adv = advection_dt.exec();
            if (!(dt_adv > TinyReal) || !std::isfinite(static_cast<double>(dt_adv)))
            {
                ++n_bad_dt;
                std::cerr << "[doc126] non-positive/non-finite dt_adv=" << dt_adv << " - abort" << std::endl;
                break;
            }
            Real dt_adv_limited = dt_adv;
            if (cli.variable_properties)
            {
                const Real mu_now = hostUpdateFrenchThermalRheologyFromTemperature(
                    particles, material_laws, /*constant_properties=*/false, cli.t_ref, rheology_stats);
                hostClampParticleViscosityCap(particles, cli.mu_cap, rheology_stats);
                const Real mu_ref_for_material =
                    use_local_particle_mu
                        ? (rheology_stats.mu_used_max > TinyReal ? rheology_stats.mu_used_max : mu_now)
                        : (rheology_stats.mu_used_mean > TinyReal ? rheology_stats.mu_used_mean : mu_now);
                fluid_viscosity.setReferenceViscosity(mu_ref_for_material);
                // Explicit viscous CFL: dt <= 0.125 h^2 rho / mu_max (local) or / mu_mean (legacy).
                // Doc175 implicit_split: estimate retained for logging only; does not limit dt_adv.
                const Real h_sm = fluid.getSPHAdaptation().ReferenceSmoothingLength();
                const Real mu_for_dt =
                    use_local_particle_mu ? rheology_stats.mu_used_max : rheology_stats.mu_used_mean;
                const Real dt_visc =
                    Real(0.125) * h_sm * h_sm * cli.rho0 / (std::max(mu_for_dt, TinyReal));
                last_dt_visc_est = dt_visc;
                if (!use_implicit_split && dt_visc < dt_adv_limited)
                {
                    dt_adv_limited = dt_visc;
                    last_active_limiter = "dt_visc_explicit";
                }
                if (!use_local_particle_mu)
                    invalidate_viscous_kernels(); // mean path caches μ by value in InteractKernel
            }
            const bool ablate_visc = (cli.doc201_b1_ablate == "visc");
            if (!ablate_visc)
            {
                if (use_implicit_split)
                {
                    hostClearViscousForcePriorContribution(particles);
                    // Doc177: optional once-per-advection to mirror explicit ForcePrior hold over acoustic subs.
                    if (cli.implicit_visc_schedule == "advection")
                        exec_implicit_viscous_split(dt_adv_limited);
                }
                else
                {
                    exec_viscous_force();
                }
            }
            else
            {
                hostClearViscousForcePriorContribution(particles);
            }
            if (use_local_particle_mu && rheology_stats.mu_mean > TinyReal)
            {
                const bool spatially =
                    (rheology_stats.mu_max / rheology_stats.mu_mean) > Real(1.05);
                // Record once early (harmonic check), then refresh when μ field becomes spatially nonuniform.
                if (!local_mu_evidence_recorded ||
                    (spatially && !local_mu_evidence.local_mu_spatially_active))
                {
                    local_mu_evidence = hostCollectLocalViscosityForceEvidence(particles);
                    local_mu_evidence_recorded = true;
                    std::cout << "[doc143][local-mu-evidence] harmonic_ok="
                              << local_mu_evidence.harmonic_mean_ok
                              << " err=" << local_mu_evidence.harmonic_max_rel_err
                              << " mu_mean/max=" << local_mu_evidence.mu_mean << "/"
                              << local_mu_evidence.mu_max
                              << " Fvisc_rms hi/lo=" << local_mu_evidence.viscous_force_rms_high_mu
                              << "/" << local_mu_evidence.viscous_force_rms_low_mu
                              << " spatially_active=" << local_mu_evidence.local_mu_spatially_active
                              << std::endl;
                }
            }
            reapply_fixed_spatial_T();
            // Doc201 B1: ablate=buoyancy must clear ForcePrior and skip the advection refresh;
            // otherwise the prior snapshot already seeds C_ext / momentum for the whole step.
            if (cli.doc201_b1_ablate == "buoyancy")
            {
                syncVariableToHost<Vecd>(particles, "ForcePrior");
                Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
                const size_t nn = particles.TotalRealParticles();
                for (size_t i = 0; i < nn; ++i)
                    fp[i] = Vecd::Zero();
                syncVariableToDevice<Vecd>(particles, "ForcePrior");
                if (cli.cref_buoyancy)
                    hostSetReferenceBodyAcceleration(particles, Vecd::Zero());
            }
            else
                refresh_buoyancy_snapshot("advection");
            if (local_iter == 0 || local_iter % 20 == 0)
            {
                const OphelieBuoyancyFieldGate gate =
                    hostAuditBuoyancyFieldGate(particles, beta_eff, last_force_tref, cli.gravity_g);
                doc189_gate_n += 1;
                doc189_worst_L2 = std::max(doc189_worst_L2, gate.E_L2);
                doc189_worst_Enet = std::max(doc189_worst_Enet, gate.E_net);
                doc189_worst_abs = std::max(doc189_worst_abs, gate.max_abs_N);
                doc189_last_gate = gate;
                if (!gate.algebraic_ok)
                    doc189_algebraic_ok = false;
                if (!gate.dynamic_reference_ok)
                    doc189_dynamic_ok = false;
                if (doc189_gate_n <= 3 || !gate.algebraic_ok || !gate.dynamic_reference_ok)
                {
                    std::cout << "[doc189][force-gate] iter=" << local_iter << " E_L2=" << gate.E_L2
                              << " E_net=" << gate.E_net << " max_abs_N=" << gate.max_abs_N
                              << " tref_req=" << gate.tref_requested_K << " tref_inf=" << gate.tref_inferred_K
                              << " ok=" << gate.algebraic_ok << "/" << gate.dynamic_reference_ok
                              << " snapshot_version=" << doc189_param_version << std::endl;
                }
            }
            // Doc144 B2: flat-top Marangoni, axisymmetric ∇_s T ≈ (dT/dr)e_r from radial-bin FD (not r-sort).
            if (std::abs(cli.marangoni_dgamma_dt) > TinyReal)
            {
                last_marangoni_apply = hostApplyMarangoniAxisymmetric(
                    particles, axis_center, R, dp, cli.marangoni_dgamma_dt, /*add_to_force_prior=*/true);
            }
            Real relaxation = 0.0;
            Real dt_ac = 0.0;
            size_t ac_sub = 0;
            const Real t_remain_phase = t_target - sv_physical_time->getValue();
            while (relaxation < dt_adv_limited && ac_sub < 5000)
            {
                dt_ac = acoustic_dt.exec();
                const Real dt_left = std::min(dt_adv_limited - relaxation, t_remain_phase - relaxation);
                if (!(dt_left > TinyReal))
                {
                    break;
                }
                if (dt_ac > dt_left)
                    dt_ac = dt_left;
                if (!(dt_ac > TinyReal) || !std::isfinite(static_cast<double>(dt_ac)))
                {
                    ++n_bad_dt;
                    std::cerr << "[doc126] bad dt_ac=" << dt_ac << " at relax=" << relaxation
                              << " dt_adv=" << dt_adv_limited << std::endl;
                    break;
                }
                const bool do_audit = density_audit.enabled;
                if (do_audit)
                    hostCopyDensityRate(particles, audit_rho_b, audit_rate_b);
                if (cli.drift_probe)
                    probe_capture_before();
                // Doc175/177: default acoustic schedule; advection schedule already applied above.
                const bool ablate_visc_ac = (cli.doc201_b1_ablate == "visc");
                const bool ablate_wall = (cli.doc201_b1_ablate == "wall");
                const bool ablate_acoustic = (cli.doc201_b1_ablate == "acoustic");
                const bool ablate_buoyancy = (cli.doc201_b1_ablate == "buoyancy");
                if (!ablate_visc_ac && use_implicit_split && cli.implicit_visc_schedule == "acoustic")
                    exec_implicit_viscous_split(dt_ac);
                if (!ablate_acoustic)
                {
                    if (cli.kernel_correction && !ablate_wall)
                        acoustic_1st_corr.exec(dt_ac);
                    else if (ablate_wall)
                        acoustic_1st_nowall.exec(dt_ac);
                    else
                        acoustic_1st.exec(dt_ac);
                }
                if (do_audit)
                    hostCopyDensityRate(particles, audit_rho_1, audit_rate_1);
                if (cli.probe_freeze_velocity)
                    probe_zero_velocity(); // continuity (2nd half) and position update then see u=0
                // A2: Pressure is current right after the 1st half (EOS(rho_half)+p_h). Read-only sample on the
                // last production substep of the phase; nothing in the production state is modified.
                if (sample_pressure && dt_ac >= (t_remain_phase - relaxation) * Real(0.999999))
                {
                    hostSamplePressureField(particles, cli.rho0, cli.c0,
                                            static_cast<double>(sv_physical_time->getValue() + dt_ac),
                                            pressure_sample);
                }
                if (!ablate_acoustic)
                {
                    if (cli.kernel_correction && !ablate_wall)
                        acoustic_2nd_corr.exec(dt_ac);
                    else if (ablate_wall)
                        acoustic_2nd_nowall.exec(dt_ac);
                    else
                        acoustic_2nd.exec(dt_ac);
                }
                if (do_audit)
                {
                    hostCopyDensityRate(particles, audit_rho_2, audit_rate_2);
                    density_audit.sampleSubstep(
                        particles, dt_ac, static_cast<double>(sv_physical_time->getValue()), phase,
                        false, audit_rho_b, audit_rate_b, audit_rho_1, audit_rate_1, audit_rho_2, audit_rate_2);
                }
                if (cli.drift_probe)
                    probe_accumulate_substep(dt_ac);
                // Same snapshot reference; do not recompute T_mean inside acoustic substeps.
                // Doc201: strat p' owns buoyancy → never re-apply ForcePrior here.
                if (!ablate_buoyancy && !strat_p_prime_owns_buoyancy)
                    boussinesq_force.exec();
                relaxation += dt_ac;
                sv_physical_time->incrementValue(dt_ac);
                ++ac_sub;
            }
            update_position.exec();
            ++local_iter;
            ++total_iter;
            const Real dt_step = sv_physical_time->getValue() - t_before;

            if (local_iter % 50 == 0)
            {
                particle_sort.exec();
            }
            fluid_cell_linked_list.exec();
            wall_closed_cell_linked_list.exec();
            wall_noslip_cell_linked_list.exec();
            fluid_update_relation.exec();
            update_classic_visc_relations();

            // Thermal update after the neighbour relation is current (Lagrangian: T is carried by the particles;
            // conduction + Q + boundary are explicit, no Eulerian convection term).
            reapply_q_shape();
            thermal_step(dt_step);
            if (n_em_updates > 0)
            {
                ++n_fluid_steps_after_last_em;
                ++n_thermal_steps_after_last_em;
            }
            maybe_em_subprocess_update(phase, false);
            maybe_write_vtp(total_iter);
            maybe_export_scheduled_fields();
            if (cli.write_em_checkpoints && cli.checkpoint_interval > TinyReal)
            {
                const Real t_now_ck = sv_physical_time->getValue();
                if (t_now_ck + Real(1.0e-12) >= next_checkpoint_time)
                {
                    write_checkpoint_now(total_iter);
                    while (next_checkpoint_time <= t_now_ck + Real(1.0e-12))
                        next_checkpoint_time += cli.checkpoint_interval;
                }
            }

            if (local_iter % 10 == 0)
            {
                const Real um = u_max_reduce.exec();
                std::cout << "[doc126][" << phase << "] N=" << local_iter << " t=" << sv_physical_time->getValue()
                          << " U_max=" << um << " dt_adv=" << dt_adv << " dt_ac=" << dt_ac << std::endl
                          << std::flush;
                if (!std::isfinite(static_cast<double>(um)))
                {
                    std::cerr << "[doc126] NON-FINITE U_max - abort" << std::endl;
                    ++n_bad_dt;
                    break;
                }
                if (um > Real(0.5) * cli.c0)
                {
                    std::cerr << "[doc126] RUNAWAY U_max=" << um << " > 0.5*c0 - abort" << std::endl;
                    ++n_bad_dt;
                    break;
                }
            }
            if (local_iter % 20 == 0)
            {
                (void)write_monitor_row(total_iter, phase);
            }
            if (n_bad_dt > 0)
                break;
        }
    };

    if (!freeze_flow && cli.kernel_correction)
    {
        fluid_linear_correction_matrix.exec();
    }

    const FlowStats s_init = write_monitor_row(0, "init");
    if (write_vtp)
    {
        write_vtp_now(0);
        next_vtp_time = cli.state_record_interval;
    }

    // Phase 1: settle (R3: quiet R0 hold before perturbation). Velocity is not zeroed afterwards.
    // Frozen-flow thermal tests skip the settle phase.
    // Doc164: restart skips settle (time already includes prior settle/observe).
    if (restarting)
    {
        const Real t_res = restart_io.readRestartFiles(static_cast<size_t>(cli.restart_step));
        sv_physical_time->setValue(t_res);
        if (!load_em_coupling_state())
            return 1;
        fluid_cell_linked_list.exec();
        wall_closed_cell_linked_list.exec();
        wall_noslip_cell_linked_list.exec();
        fluid_update_relation.exec();
        update_classic_visc_relations();
        // Doc167: do not rewrite mu/cp/k from T after restart — evolving fields already restored.
        // Only refresh rheology stats + viscous CFL reference from restored particle mu.
        if (cli.variable_properties)
        {
            // Restart: field already restored; treat current as raw then clamp → used.
            syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
            syncVariableToHost<Real>(particles, "VolumetricMeasure");
            const Real *mu_f = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
            const Real *vol_f = particles.getVariableDataByName<Real>("VolumetricMeasure");
            const size_t nn = particles.TotalRealParticles();
            hostScalarFieldStats(mu_f, vol_f, nn, rheology_stats.mu_raw_min, rheology_stats.mu_raw_max,
                                 rheology_stats.mu_raw_mean);
            rheology_stats.mu_stats_stage = "restart_pre_cap_raw";
            hostClampParticleViscosityCap(particles, cli.mu_cap, rheology_stats);
            const Real mu_ref = use_local_particle_mu
                                    ? (rheology_stats.mu_used_max > TinyReal ? rheology_stats.mu_used_max : cli.mu)
                                    : (rheology_stats.mu_used_mean > TinyReal ? rheology_stats.mu_used_mean : cli.mu);
            fluid_viscosity.setReferenceViscosity(mu_ref);
        }
        // Doc169: checkpoint cadence continues from next barrier strictly after t_res.
        if (cli.write_em_checkpoints && cli.checkpoint_interval > TinyReal)
        {
            next_checkpoint_time = cli.checkpoint_interval;
            while (next_checkpoint_time <= t_res + Real(1.0e-12))
                next_checkpoint_time += cli.checkpoint_interval;
        }
        // Doc171: keep mass_sum_init from sidecar when present; report readback separately.
        wall_closed_normal.exec();
        wall_noslip_normal.exec();
        const double mass_readback = hostTotalMass(particles);
        std::cout << "[doc171][restart] active t=" << t_res << " step=" << cli.restart_step
                  << " next_ckpt=" << next_checkpoint_time << " mass_sum_init=" << mass_sum_init
                  << " mass_readback=" << mass_readback << std::endl;
        refresh_buoyancy_snapshot("restart");
    }
    else if (!freeze_flow)
    {
        // Doc171 P3: replace prior S0 Q with T-consistent u0 before physical advance.
        if (cli.em_initial_consistent)
        {
            if (!em_coupling_armed)
            {
                std::cerr << "[doc171] --em-initial-consistent requires --em-update-mode=subprocess"
                          << std::endl;
                return 1;
            }
            if (cli.variable_properties)
            {
                (void)hostUpdateFrenchThermalRheologyFromTemperature(
                    particles, material_laws, /*constant_properties=*/false, cli.t_ref, rheology_stats);
                hostClampParticleViscosityCap(particles, cli.mu_cap, rheology_stats);
            }
            std::cout << "[doc171][em-initial] solving u0 at t=" << sv_physical_time->getValue()
                      << std::endl;
            maybe_em_subprocess_update("initial", true);
            if (!em_coupling_ok)
            {
                std::cerr << "[doc171][em-initial] FAIL " << em_failure_reason
                          << " (refusing inconsistent prior Q)" << std::endl;
                return 1;
            }
        }
        std::cout << "[doc203][f1-heartbeat] entering settle to t=" << settle_eff << std::endl
                  << std::flush;
        advance_to(settle_eff, "settle", false);
    }
    const FlowStats s_settle = write_monitor_row(total_iter, "settle");
    std::cout << "[doc203][f1-heartbeat] settle done U_rms=" << s_settle.U_rms << std::endl
              << std::flush;
    FlowStats s_quiet = s_settle;
    if (cli.zero_velocity_after_settle)
    {
        // Density first: leftover δρ after packing re-excites U~c0·δρ/ρ within ~0.2 s if only U is cleared.
        hostResetDensityToRho0(particles, cli.rho0, cli.quiet_reset_mass);
        hostZeroFluidVelocity(particles);
        s_quiet = write_monitor_row(total_iter, "quiet");
        std::cout << "[doc203][f1] quiet-after-settle (rho→rho0"
                  << (cli.quiet_reset_mass ? " + Mass=rho0*Vol" : " keepMass") << " + U→0) U_rms "
                  << s_settle.U_rms << " -> " << s_quiet.U_rms << " rho=[" << s_quiet.rho_min << ","
                  << s_quiet.rho_max << "] mass=" << hostTotalMass(particles) << std::endl
                  << std::flush;
    }
    activate_fixed_T_after_settle();

    // R3: add the velocity bump after the quiet hold.
    FlowStats s_perturb = s_quiet;
    if (is_R3)
    {
        hostApplyR3VelocityPerturbation(particles, axis_center, R, H, z_bottom, cli.perturb_u);
        s_perturb = write_monitor_row(total_iter, "perturb");
        std::cout << "[doc126] R3 perturb applied eps=" << cli.perturb_u << " U_rms_before=" << s_settle.U_rms
                  << " U_rms_after=" << s_perturb.U_rms << std::endl;
    }

    // Doc201 B2: live top-shear IC after settle (default settle=0). Positions remain live; no uz projection.
    Real b2_t0 = 0;
    std::vector<Vecd> b2_pos0;
    Doc201B2Metrics b2_core{}, b2_annulus{};
    if (cli.doc201_b2_live_top_shear)
    {
        hostApplyDoc201B2TopShearIC(particles, z_bottom, H, cli.doc201_b2_U0);
        b2_t0 = sv_physical_time->getValue();
        syncVariableToHost<Vecd>(particles, "Position");
        const Vecd *p0 = particles.getVariableDataByName<Vecd>("Position");
        b2_pos0.assign(p0, p0 + particles.TotalRealParticles());
        s_perturb = write_monitor_row(total_iter, "b2_ic");
        std::cout << "[doc201][b2] live_top_shear IC U0=" << cli.doc201_b2_U0
                  << " analytic=U0*sin(pi*z/(2H))*exp(-nu*(pi/(2H))^2*t)"
                  << " acoustic_lid=wall_closed positions=live viscous_top=omit"
                  << " U_rms_after=" << s_perturb.U_rms << std::endl;
    }

    // Doc164: on restart, --end-time is absolute physical end; fresh runs keep observe duration.
    const Real t_obs = restarting ? cli.end_time : (sv_physical_time->getValue() + cli.end_time);
    if (freeze_flow)
    {
        advance_frozen(t_obs, "observe");
    }
    else if (em_coupling_armed && cli.em_schedule == "periodic")
    {
        // Doc164/167: half-open [t_now, t_obs). Advance to the sorted union of
        //   (a) EM times origin+n*interval and
        //   (b) checkpoint barriers n*checkpoint_interval (when enabled).
        // Barriers force the same dt truncation as a segmented stop at that time (Doc167 root cause).
        // EM fires only via maybe_em_subprocess_update when t_now >= next_em_due.
        struct Hit
        {
            Real t;
            bool em;
        };
        std::vector<Hit> hits;
        for (long n = 1;; ++n)
        {
            const Real t_ev = cli.em_update_origin + static_cast<Real>(n) * cli.em_update_interval;
            if (!(t_ev + Real(1.0e-12) < t_obs))
                break;
            hits.push_back({t_ev, true});
        }
        if (cli.write_em_checkpoints && cli.checkpoint_interval > TinyReal)
        {
            for (long n = 1;; ++n)
            {
                const Real t_ck = static_cast<Real>(n) * cli.checkpoint_interval;
                if (!(t_ck + Real(1.0e-12) < t_obs))
                    break;
                hits.push_back({t_ck, false});
            }
        }
        std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
            if (a.t < b.t)
                return true;
            if (b.t < a.t)
                return false;
            return a.em && !b.em; // EM before bare barrier at same time
        });
        // Dedup near-equal times (keep EM flag if either is EM)
        std::vector<Hit> uniq;
        for (const Hit &h : hits)
        {
            if (!uniq.empty() && std::abs(static_cast<double>(uniq.back().t - h.t)) < 1.0e-12)
            {
                uniq.back().em = uniq.back().em || h.em;
                continue;
            }
            uniq.push_back(h);
        }
        const Real t_now0 = sv_physical_time->getValue();
        for (const Hit &h : uniq)
        {
            if (!em_coupling_ok)
                break;
            if (h.t + Real(1.0e-12) <= t_now0)
                continue; // already past (restart)
            if (!(h.t + Real(1.0e-12) < t_obs))
                break;
            const size_t n_em_before = n_em_updates;
            advance_to(h.t, "observe", true);
            if (h.em && static_cast<long>(n_em_updates) < max_em_updates_eff)
            {
                maybe_em_subprocess_update("observe", false);
                if (!em_coupling_ok)
                    break;
                if (n_em_updates <= n_em_before)
                {
                    em_coupling_ok = false;
                    em_failure_reason = "em_periodic_update_missed";
                    break;
                }
            }
        }
        if (em_coupling_ok)
            advance_to(t_obs, "observe", true);
        if (em_coupling_ok && static_cast<long>(n_em_updates) >= max_em_updates_eff)
            std::cout << "[doc164][em-schedule] periodic feedback stopped at max_em_updates="
                      << max_em_updates_eff << std::endl;
    }
    else if (em_coupling_armed && !cli.em_update_times.empty())
    {
        // Doc156: advance to each absolute event time, then continue to end (no end flush).
        for (size_t ie = 0; ie < cli.em_update_times.size(); ++ie)
        {
            const Real t_ev = cli.em_update_times[ie];
            if (t_ev > t_obs + Real(1.0e-12))
            {
                em_coupling_ok = false;
                em_failure_reason = "em_update_time_past_end";
                break;
            }
            advance_to(t_ev, "observe", true);
            maybe_em_subprocess_update("observe", false);
            if (!em_coupling_ok)
                break;
            if (static_cast<long>(n_em_updates) < static_cast<long>(ie + 1))
            {
                em_coupling_ok = false;
                em_failure_reason = "em_update_missed_scheduled_time";
                break;
            }
        }
        if (em_coupling_ok)
            advance_to(t_obs, "observe", true);
    }
    else
    {
        advance_to(t_obs, "observe", true);
    }

    // Doc151 legacy only: optional end flush. Doc156 default --no-em-end-flush.
    if (cli.em_end_flush && em_coupling_armed && em_coupling_ok)
    {
        while (static_cast<long>(n_em_updates) < max_em_updates_eff)
        {
            next_em_due = Real(0);
            maybe_em_subprocess_update("observe", false);
            if (!em_coupling_ok)
                break;
        }
    }

    // ---- Production end state: all production quantities are read HERE, before any diagnostic evaluation ----
    if (cli.write_em_checkpoints)
        write_checkpoint_now(total_iter);
    maybe_export_scheduled_fields();
    export_field_now("end", "production_end_state");
    const FlowStats s_end = write_monitor_row(total_iter, "end");
    write_vtp_now(total_iter + 1);
    const PiStats pi_density_end = hostPiStatsFromDensity(particles, cli.rho0, cli.c0);
    const double mass_sum_end = hostTotalMass(particles);
    const double h_end = thermal_enthalpy_J();
    const DeltaTStats dT_end = hostDeltaTStats(particles);
    const double time_production_end = static_cast<double>(sv_physical_time->getValue());

    // D4: boundary-band / escape report at production end (read-only)
    const BoundaryBandStats bb_end = hostBoundaryBandStats(particles, axis_center, R, H, z_bottom, dp);

    // Diagnostic force evaluation (after production end; does not advance time or velocity).
    reapply_fixed_spatial_T();
    refresh_buoyancy_snapshot("diagnostic_end");
    const OphelieBoussinesqForceAudit force_end = hostAuditBoussinesqForcePrior(
        particles, beta_eff, last_force_tref, -cli.gravity_g, kOphelieTemperatureField);
    logOphelieBoussinesqForceAudit(force_end, "doc126_end");
    const double buoy_a_rms_end = hostBuoyancyAccelRms(particles);

    const bool mass_conservation_pass = std::abs(mass_sum_end - mass_sum_init) <=
                                        1.0e-12 * (std::abs(mass_sum_init) + 1.0);

    // A2 semantics split.
    const bool aux_density_range_ok =
        s_end.finite && s_end.rho_min > 0.9 * static_cast<double>(cli.rho0) &&
        s_end.rho_max < 1.1 * static_cast<double>(cli.rho0);
    const bool eos_identity_ok_legacy = eos_initialization_ok && aux_density_range_ok; // NOT field consistency
    const bool pressure_field_available = pressure_sample.valid;
    const bool eos_field_consistency_ok = pressure_field_available && pressure_sample.consistent;

    const bool no_runaway = s_end.finite && s_end.U_max < 0.25 * static_cast<double>(cli.c0);
    const bool finite_all = s_end.finite && s_perturb.finite && s_settle.finite && dT_end.finite;

    // ---- Energy ledger (Doc131 §10): R_E = DeltaH - (E_source - E_boundary) ----
    const double d_h = h_end - h_init;
    const double r_e = d_h - (ledger.e_source - ledger.e_boundary);
    const double e_floor = 1.0; // J
    const double err_input = std::abs(r_e) / std::max(std::abs(ledger.e_source), e_floor);
    const double err_exchange =
        std::abs(r_e) / std::max({std::abs(ledger.e_source), std::abs(ledger.e_boundary), e_floor});
    const double e_net = ledger.e_source - ledger.e_boundary;
    const bool net_ok = std::abs(e_net) > 0.01 * std::max(std::abs(ledger.e_source), std::abs(ledger.e_boundary));
    const double err_net = net_ok ? std::abs(r_e) / std::abs(e_net) : std::numeric_limits<double>::quiet_NaN();
    // Doc177 unified residual: |R_E| / max(|E_source|, |E_boundary|, |ΔH|, floor)
    const double err_unified = std::abs(r_e) / std::max({std::abs(ledger.e_source), std::abs(ledger.e_boundary),
                                                         std::abs(d_h), e_floor});

    // ---- Pass criteria ----
    bool passed = false;
    std::string note;
    const bool reapply_flag = false; // no analytic density rewrite exists in this case
    const bool buoyancy_present = buoy_a_rms_end > 1.0e-8;
    const double circ_mag = std::max(std::abs(s_end.circ_ur_proxy), std::abs(s_end.circ_uz_proxy));
    const bool circ_nonzero = circ_mag > 1.0e-8;
    bool r2_signal = false;
    double r2_dT_sign = 0.0;
    bool r2_null_driver = false, r2_sign_uz_ok = false, r2_sign_ur_ok = false;
    std::string r2_gate_kind = "n/a";
    double u_growth = 0.0;
    double decay_ratio = 0.0;
    JsonWriter tj; // mode-specific thermal-test metrics
    Gate gate_source_power = Gate::Null;
    std::string gate_source_power_reason = "OPEN: no fixed-power source in this mode";
    Gate gate_conduction = Gate::Null;
    std::string gate_conduction_reason = "OPEN: run --mode=T2 (and T4) for the conduction check";
    Gate gate_thermal_bc = Gate::Null;
    std::string gate_thermal_bc_reason = "OPEN: run --mode=T3 (one run per face) for the boundary-flux check";

    const Real b2_t = std::max(Real(0), sv_physical_time->getValue() - b2_t0);
    const Real b2_nu = cli.mu / std::max(cli.rho0, TinyReal);
    if (cli.doc201_b2_live_top_shear)
    {
        const Real r_core = Real(0.5) * R;
        const Real r_full = std::max(Real(0), R - Real(3.5) * dp);
        const Vecd *pos0_ptr = b2_pos0.empty() ? nullptr : b2_pos0.data();
        b2_core = hostEvaluateDoc201B2TopShear(particles, axis_center, R, H, z_bottom, b2_t, cli.doc201_b2_U0,
                                               b2_nu, dp, r_core, pos0_ptr, cli.rho0, cli.c0);
        b2_annulus = hostEvaluateDoc201B2TopShear(particles, axis_center, R, H, z_bottom, b2_t, cli.doc201_b2_U0,
                                                  b2_nu, dp, r_full, pos0_ptr, cli.rho0, cli.c0);
        // Profile dump (all particles; includes uref).
        {
            const std::string pdir = cli.output_dir.empty() ? std::string(".") : cli.output_dir;
            std::filesystem::create_directories(pdir);
            syncVariableToHost<Vecd>(particles, "Position");
            syncVariableToHost<Vecd>(particles, "Velocity");
            const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
            const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
            std::ofstream pf(pdir + "/b2_profile_end.csv");
            pf << "x,y,z,r,ux,uy,uz,ux_ref\n";
            const size_t nn = particles.TotalRealParticles();
            for (size_t i = 0; i < nn; ++i)
            {
                const Real dx = pos[i][0] - axis_center[0];
                const Real dy = pos[i][1] - axis_center[1];
                const Real rr = std::sqrt(dx * dx + dy * dy);
                const Real zeta = pos[i][2] - z_bottom;
                const Real uref = analyticDoc201B2TopFreeSlipUx(zeta, b2_t, cli.doc201_b2_U0, H, b2_nu);
                pf << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << rr << "," << vel[i][0] << ","
                   << vel[i][1] << "," << vel[i][2] << "," << uref << "\n";
            }
        }
        std::cout << "[doc201][b2] t=" << b2_t << " core L2_rel=" << b2_core.L2_rel
                  << " amp_rel=" << b2_core.amp_rel_err << " uz/U0=" << b2_core.uz_rms_over_U0
                  << " n_core=" << b2_core.n_interior << " disp_rms=" << b2_core.disp_rms
                  << " | annulus L2_rel=" << b2_annulus.L2_rel << " n_ann=" << b2_annulus.n_interior
                  << std::endl;
    }

    const double t_th = ledger.t_thermal;
    const bool power_discrete_ok =
        power.mode == kPowerDiscrete && power.ok &&
        std::abs(power.P_applied - power.P_target) <= 1.0e-5 * std::max(power.P_target, 1.0);
    if (q_active)
    {
        gate_source_power = power_discrete_ok ? Gate::Yes : Gate::No;
        gate_source_power_reason = power.mode == kPowerDiscrete
                                       ? "P_applied differs from P_target by more than 1e-5 relative"
                                       : "fixed_volumetric_density is diagnostic: P_applied != P_target by V_geom/sum_V_Q";
    }

    if (cli.doc201_b2_live_top_shear)
    {
        // Plan §9: profile L2≤5%, uz_rms/U0≤1%. Positions must have moved (disp_rms>0) — not frozen.
        const bool live_ok = b2_core.disp_rms > TinyReal;
        const bool profile_ok = b2_core.L2_rel <= Real(0.05) && b2_core.amp_rel_err <= Real(0.05);
        const bool normal_ok = b2_core.uz_rms_over_U0 <= Real(0.01);
        passed = finite_all && no_runaway && n_bad_dt == 0 && live_ok && profile_ok && normal_ok &&
                 b2_core.n_interior > 0;
        note = "Doc201 B2 live acoustic top: wall_closed lid + top viscous omit; core r<0.5R vs "
               "U0*sin(pi z/(2H))*exp(-nu (pi/(2H))^2 t); gates L2/amp≤5%, uz/U0≤1%, positions live";
    }
    else if (is_R0)
    {
        passed = finite_all && no_runaway && (s_end.U_rms < 5.0e-4) && n_bad_dt == 0;
        note = "R0 quiet hold: T=Tref, U_rms<5e-4 mech budget";
    }
    else if (is_R1)
    {
        if (cli.thaw_T)
        {
            // Doc201 short thaw case: report-only finite gate (not a calibrated free-convection PASS).
            passed = finite_all && no_runaway && n_bad_dt == 0;
            note = strat_p_prime_owns_buoyancy
                       ? "R1+thaw-T with strat_owns kept (diagnostic; ForcePrior off while T evolves)"
                       : "R1+thaw-T with Boussinesq ForcePrior on (strat_init is IC only; not literature calibrated)";
        }
        else
        {
            passed = finite_all && no_runaway && (s_end.U_rms < 1.0e-3) && !reapply_flag && n_bad_dt == 0;
            note = "R1 stable strat hold: rho init from p'_strat only, no per-step analytic rewrite";
        }
    }
    else if (is_R2)
    {
        // Prefer quiet baseline when --zero-velocity-after-settle cleared packing noise.
        u_growth = s_end.U_rms / (s_quiet.U_rms + 1.0e-30);
        // Doc131 D2: sign-correct gate (replaces the old "nonzero circulation OR buoyancy present" gate).
        // Hot wall (dT>0, beta>0): uz_proxy > 0 (upflow at wall) AND ur_proxy < 0 (inflow aloft); both flip with
        // the sign of dT. beta=0 (or dT=0) has no driver: the run must stay quiet (U_rms < 0.1 x signal threshold).
        r2_dT_sign = r2_dT_eff > 0 ? 1.0 : (r2_dT_eff < 0 ? -1.0 : 0.0);
        r2_null_driver = (beta_eff <= TinyReal) || (r2_dT_sign == 0.0);
        r2_sign_uz_ok = r2_dT_sign * s_end.circ_uz_proxy > 0.0;
        r2_sign_ur_ok = r2_dT_sign * s_end.circ_ur_proxy < 0.0;
        if (r2_null_driver)
        {
            r2_signal = s_end.U_rms < 1.0e-5;
            r2_gate_kind = "null_driver_must_stay_quiet";
        }
        else
        {
            r2_signal = (s_end.U_rms > 1.0e-4) && r2_sign_uz_ok && r2_sign_ur_ok;
            r2_gate_kind = "driven_sign_correct_uz_and_ur";
        }
        passed = finite_all && no_runaway && r2_signal && n_bad_dt == 0;
        note = "R2 radial driver (hot wall for dT>0, zero-mean unless --r2-keep-mean). Gate (Doc131 D2): driven -> "
               "U_rms>1e-4 AND sign(dT)*uz_proxy>0 AND sign(dT)*ur_proxy<0; beta=0/dT=0 -> U_rms<1e-5 (quiet). "
               "Cross-run reversal / beta=0 comparison is in the D2 matrix summary.";
    }
    else if (is_R3)
    {
        decay_ratio = s_end.U_rms / (s_perturb.U_rms + 1.0e-30);
        const bool perturb_nonzero = s_perturb.U_rms > 0.0 && std::isfinite(s_perturb.U_rms);
        const bool bounded = decay_ratio < 5.0;
        passed = finite_all && no_runaway && perturb_nonzero && bounded && n_bad_dt == 0;
        note = "R3 perturbation response: nonzero, finite; decays or stays bounded (end/perturb<5)";
    }
    else if (is_baseline)
    {
        if (legacy_path)
        {
            const bool heated_ok = s_end.T_mean > s_init.T_mean - 1.0e-9 || cli.q_power <= TinyReal;
            passed = finite_all && no_runaway && heated_ok && n_bad_dt == 0;
            note = "baseline thermal_smoke_baseline (LEGACY --legacy-rad-clamp: top-only host radiation with "
                   "max(T,Tamb) clamp, no conduction, no side/bottom loss); NOT a literature baseline, NOT "
                   "a calibrated match";
        }
        else if (!cli.thaw_T)
        {
            passed = finite_all && no_runaway && n_bad_dt == 0;
            if (doc189_gate_n > 0)
                passed = passed && doc189_algebraic_ok && doc189_dynamic_ok;
            note = cli.doc189_mechanism
                       ? "Doc189 prescribed-T mechanism: flow+wall+pressure+implicit viscosity; energy_gate=NOT_APPLICABLE"
                       : "baseline thermal_smoke_baseline with T frozen (no Q, no BC); flow diagnostics only; "
                         "NOT a literature baseline, NOT a calibrated match";
        }
        else
        {
            const bool energy_ok = err_input <= 0.01;
            bool em_gate = true;
            if (em_coupling_armed)
            {
                // Doc164: periodic schedule uses max_em_updates as a safety cap, not a required count.
                // Absolute/times list still requires completing exactly the scheduled updates.
                if (cli.em_schedule == "periodic")
                {
                    em_gate = em_coupling_ok && (n_em_updates > 0) && (n_fluid_steps_after_last_em > 0);
                }
                else
                {
                    em_gate = em_coupling_ok &&
                              (static_cast<long>(n_em_updates) == max_em_updates_eff);
                    if (em_gate && !cli.em_update_times.empty())
                        em_gate = (n_fluid_steps_after_last_em > 0);
                }
                if (em_gate && cli.em_sigma_reduction == "volume_mean" &&
                    std::abs(cli.em_diag_T_perturb_K) > TinyReal)
                    em_gate = std::abs(last_delta_sigma) > 1.0e-6;
            }
            passed = finite_all && no_runaway && n_bad_dt == 0 && power_discrete_ok && energy_ok && em_gate;
            if (doc189_gate_n > 0)
            {
                passed = passed && doc189_algebraic_ok && doc189_dynamic_ok;
                note += " Doc189 in-loop buoyancy field gate required.";
            }
            note = em_coupling_armed
                       ? (cli.em_sigma_reduction == "spatial_field"
                              ? "baseline + Doc156 spatial-σ subprocess T→σ→L1 Picard Q refresh. Gate: energy<=1%, "
                                "scheduled EM (no end flush), post-update steps>0. NOT literature circulation match."
                              : "baseline + Doc151 subprocess T→σ→L1 Picard Q refresh (volume-mean σ). Gate: energy<=1%, "
                                "scheduled EM updates succeeded, nonzero Δσ when diag T-perturb set. NOT literature "
                                "circulation match.")
                       : "baseline thermal_smoke_baseline: fixed discrete Q + production conduction + Robin/radiation "
                         "BC + thawed T + Boussinesq buoyancy. Gate: energy ledger error_input<=1% and "
                         "P_applied=P_target. CODE_PROVISIONAL h/t_cool; NOT a literature baseline, NOT a calibrated "
                         "match";
        }
    }
    else if (is_T0)
    {
        const bool t_const = dT_end.abs_max <= 1.0e-6 && std::abs(d_h) <= 1.0e-3;
        passed = finite_all && n_bad_dt == 0 && t_const && ledger.n_steps > 0;
        tj.num("T0_dT_abs_max_K", dT_end.abs_max);
        tj.num("T0_dH_J", d_h);
        note = "T0 no Q, adiabatic, frozen flow: |DeltaT|max<=1e-6 K and |DeltaH|<=1e-3 J";
    }
    else if (is_T1)
    {
        const double dh_exp = static_cast<double>(cli.q_power) * t_th;
        const double dh_rel = std::abs(d_h - dh_exp) / std::max(std::abs(dh_exp), 1.0);
        const double dT_exp = dh_exp / (power.mass_sum * static_cast<double>(cli.cp));
        const double dT_rel = std::abs(dT_end.mean_m - dT_exp) / std::max(std::abs(dT_exp), 1e-30);
        const double dT_spread = dT_end.max - dT_end.min;
        const bool ok = dh_rel <= 1.0e-4 && dT_rel <= 1.0e-4 && power_discrete_ok;
        passed = finite_all && n_bad_dt == 0 && ok;
        tj.num("T1_dH_J", d_h);
        tj.num("T1_dH_expected_P_t_J", dh_exp);
        tj.num("T1_dH_rel_err", dh_rel);
        tj.num("T1_dT_mean_K", dT_end.mean_m);
        tj.num("T1_dT_expected_K", dT_exp);
        tj.num("T1_dT_rel_err", dT_rel);
        tj.num("T1_dT_spread_K", dT_spread);
        note = "T1 50 kW discrete uniform Q, no loss, frozen flow: DeltaH vs P_target*t and DeltaT vs P t/(M cp) "
               "(rel<=1e-4; independent of the code's own P_applied integral)";
    }
    else if (is_T2)
    {
        const bool conserved = std::abs(d_h) <= 1.0e-5 * std::max(h_abs_scale_init, 1.0);
        const double spread0 = dT_init.max - dT_init.min;
        const double spread1 = dT_end.max - dT_end.min;
        const bool relaxing = dT_end.var_m < dT_init.var_m && spread1 < spread0;
        passed = finite_all && n_bad_dt == 0 && conserved && relaxing;
        tj.num("T2_dH_J", d_h);
        tj.num("T2_H_abs_scale_J", h_abs_scale_init);
        tj.num("T2_dH_over_scale", d_h / std::max(h_abs_scale_init, 1.0));
        tj.flag("T2_energy_conserved", conserved);
        tj.num("T2_dT_spread_start_K", spread0);
        tj.num("T2_dT_spread_end_K", spread1);
        tj.num("T2_dT_var_start_K2", dT_init.var_m);
        tj.num("T2_dT_var_end_K2", dT_end.var_m);
        tj.num("T2_dT_mean_start_K", dT_init.mean_m);
        tj.num("T2_dT_mean_end_K", dT_end.mean_m);
        tj.flag("T2_relaxing", relaxing);
        gate_conduction = (conserved && relaxing) ? Gate::Yes : Gate::No;
        gate_conduction_reason = "T2 energy not conserved or temperature spread/variance did not decrease";
        note = "T2 no Q, adiabatic, non-uniform start, frozen flow: total energy conserved (|dH|<=1e-5 sum m cp "
               "|dT0|) and the temperature variance/spread decrease (heat flows down-gradient)";
    }
    else if (is_T3)
    {
        const double A_geom = cli.t3_face == "side" ? area.side.target
                                                    : (cli.t3_face == "bottom" ? area.bottom.target : area.top.target);
        const double T0d = static_cast<double>(cli.t_ref);
        double flux = 0.0;
        if (cli.t3_face == "side")
            flux = static_cast<double>(cli.h_side) * (T0d - static_cast<double>(cli.t_cool));
        else if (cli.t3_face == "bottom")
            flux = static_cast<double>(cli.h_bottom) * (T0d - static_cast<double>(cli.t_cool));
        else
            flux = static_cast<double>(cli.h_top) * (T0d - static_cast<double>(cli.t_ambient)) +
                   static_cast<double>(cli.emissivity) * kStefanBoltzmann *
                       (std::pow(T0d, 4.0) - std::pow(static_cast<double>(cli.t_ambient), 4.0));
        const double p_exp = flux * A_geom;
        const double e_exp = p_exp * t_th;
        const double e_loss_field = -d_h;
        const double rel_e = std::abs(e_loss_field - e_exp) / std::max(std::abs(e_exp), 1.0);
        const double rel_p0 = std::abs(ledger.p_loss_first - p_exp) / std::max(std::abs(p_exp), 1.0);
        const double tol = cli.area_normalize ? 0.02 : 0.15;
        const bool ok = rel_e <= tol && rel_p0 <= tol;
        passed = finite_all && n_bad_dt == 0 && ok;
        tj.str("T3_face", cli.t3_face);
        tj.num("T3_A_geom_m2", A_geom);
        tj.num("T3_flux_Wm2_analytic", flux);
        tj.num("T3_P_analytic_W", p_exp);
        tj.num("T3_P_loss_initial_W_ledger", ledger.p_loss_first);
        tj.num("T3_P_loss_initial_rel_err", rel_p0);
        tj.num("T3_E_expected_J", e_exp);
        tj.num("T3_E_loss_from_dH_J", e_loss_field);
        tj.num("T3_E_loss_rel_err", rel_e);
        tj.num("T3_tol", tol);
        gate_thermal_bc = ok ? Gate::Yes : Gate::No;
        gate_thermal_bc_reason = "T3 face loss differs from analytic A_geom*flux*t beyond tolerance";
        note = "T3 isothermal short run, single face (--t3-face): loss from field DeltaH vs analytic "
               "flux(T0)*A_geom*t (independent geometry/flux formula); tol 2% with --area-normalize else 15%";
    }
    else if (is_T4)
    {
        double c1 = 0.0, a1 = 0.0, res1 = 0.0;
        hostFitCosineProfile(particles, H, z_bottom, c1, a1, res1);
        const double rate = alpha_th * (kPiD / static_cast<double>(H)) * (kPiD / static_cast<double>(H));
        const double ratio_an = std::exp(-rate * t_th);
        const double ratio_num = a1 / (t4_a0 + 1e-300);
        const double decay_err = std::abs(std::log(std::max(ratio_num, 1e-300)) - std::log(ratio_an)) /
                                 std::max(std::abs(std::log(ratio_an)), 1e-12);
        const double c_drift = std::abs(c1 - t4_c0);
        const bool conserved = std::abs(d_h) <= 1.0e-4 * std::max(h_abs_scale_init, 1.0);
        const bool ok = decay_err <= 0.10 && conserved;
        passed = finite_all && n_bad_dt == 0 && ok;
        tj.num("T4_A0_fit_K", t4_a0);
        tj.num("T4_A_end_fit_K", a1);
        tj.num("T4_ratio_numeric", ratio_num);
        tj.num("T4_ratio_analytic", ratio_an);
        tj.num("T4_decay_rate_rel_err", decay_err);
        tj.num("T4_fit_rms_resid_start_K", t4_res0);
        tj.num("T4_fit_rms_resid_end_K", res1);
        tj.num("T4_mean_drift_K", c_drift);
        tj.num("T4_dH_J", d_h);
        gate_conduction = ok ? Gate::Yes : Gate::No;
        gate_conduction_reason = "T4 cosine decay-rate error > 10% or energy not conserved";
        note = "T4 1D cosine profile dT=A cos(pi z/H), all faces adiabatic, frozen flow: amplitude decay "
               "exp(-alpha (pi/H)^2 t) (LS fit, decay-rate rel err<=10%) and energy conserved; "
               "engineering gate (SPH truncation near z=0,H)";
    }
    else if (is_T5)
    {
        const bool energy_ok = err_input <= 0.01;
        passed = finite_all && no_runaway && n_bad_dt == 0 && power_discrete_ok && energy_ok;
        note = "T5 full BC + fixed discrete Q, buoyancy off (beta=0), flow runs: energy ledger error_input<=1% "
               "and P_applied=P_target; CODE_PROVISIONAL h/t_cool";
    }

    std::cout << "[doc126] end mode=" << cli.mode << " t=" << sv_physical_time->getValue()
              << " U_rms=" << s_end.U_rms << " U_settle=" << s_settle.U_rms << " U_max=" << s_end.U_max
              << " Uz_rms=" << s_end.Uz_rms << " Ur_rms=" << s_end.Ur_rms << " Utheta_rms=" << s_end.Utheta_rms
              << " K=" << s_end.K << " T_mean=" << s_end.T_mean << " circ_ur=" << s_end.circ_ur_proxy
              << " circ_uz=" << s_end.circ_uz_proxy << " buoy_a_rms=" << buoy_a_rms_end
              << " eos_init_ok=" << (eos_initialization_ok ? 1 : 0)
              << " aux_rho_range_ok=" << (aux_density_range_ok ? 1 : 0)
              << " eos_field_consistent=" << (pressure_field_available ? (eos_field_consistency_ok ? 1 : 0) : -1)
              << " mass_ok=" << (mass_conservation_pass ? 1 : 0) << " n_bad_dt=" << n_bad_dt
              << " passed=" << (passed ? 1 : 0) << std::endl;
    if (thermal_active)
        std::cout << "[doc131][ledger] t_thermal=" << ledger.t_thermal << " steps=" << ledger.n_steps
                  << " dH=" << d_h << " E_source=" << ledger.e_source << " E_boundary=" << ledger.e_boundary
                  << " R_E=" << r_e << " err_input=" << err_input << " err_exchange=" << err_exchange
                  << " err_net=" << err_net << std::endl;

    // diagnostic_consistent: no bad flags (finite fields, valid dt, mass conserved, EOS field consistent when
    // sampled). The production end state (s_end, pi, dH, mass) was read before the diagnostic force evaluation, and
    // the sample/probe reads are read-only (no state mutated before production end).
    const bool diagnostic_consistent =
        finite_all && n_bad_dt == 0 && mass_conservation_pass && (!pressure_field_available || eos_field_consistency_ok);

    if (!cli.summary_json.empty())
    {
        const std::string em_power_ledger_csv =
            (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/em_power_ledger.csv";
        {
            std::ofstream led(em_power_ledger_csv);
            led << "stage,t_s,P_em_raw_W,P_particles_W,P_table_raw_W,scale_table_raw_to_target,"
                   "P_fluid_mapped_raw_W,coeff_sph_to_applied,P_fluid_applied_W,P_target_W\n";
            led.setf(std::ios::scientific);
            led.precision(10);
            for (const EmPowerLedgerRow &row : em_power_ledger)
            {
                led << row.stage << "," << row.t_s << ",";
                auto cell = [&](double v) {
                    if (!std::isfinite(v))
                        led << "";
                    else
                        led << v;
                };
                cell(row.P_em_raw_W);
                led << ",";
                cell(row.P_particles_W);
                led << ",";
                cell(row.P_table_raw_W);
                led << ",";
                cell(row.scale_table_raw_to_target);
                led << "," << row.P_fluid_mapped_raw_W << "," << row.coeff_sph_to_applied << ","
                    << row.P_fluid_applied_W << "," << row.P_target_W << "\n";
            }
        }
        std::ostringstream power_stages_json;
        power_stages_json << "[";
        for (size_t i = 0; i < em_power_ledger.size(); ++i)
        {
            const EmPowerLedgerRow &row = em_power_ledger[i];
            if (i)
                power_stages_json << ",";
            power_stages_json << "{\"stage\":\"" << jesc(row.stage) << "\",\"t_s\":" << jnum(row.t_s)
                              << ",\"P_em_raw_W\":" << jnum(row.P_em_raw_W)
                              << ",\"P_particles_W\":" << jnum(row.P_particles_W)
                              << ",\"P_table_raw_W\":" << jnum(row.P_table_raw_W)
                              << ",\"scale_table_raw_to_target\":" << jnum(row.scale_table_raw_to_target)
                              << ",\"P_fluid_mapped_raw_W\":" << jnum(row.P_fluid_mapped_raw_W)
                              << ",\"coeff_sph_to_applied\":" << jnum(row.coeff_sph_to_applied)
                              << ",\"P_fluid_applied_W\":" << jnum(row.P_fluid_applied_W)
                              << ",\"P_target_W\":" << jnum(row.P_target_W) << "}";
        }
        power_stages_json << "]";

        JsonWriter js;
        js.str("case", "doc126_french_natural_reduced_pressure");
        js.integer("doc131", 1);
        js.integer("doc139", 1);
        js.str("literature_id", cli.literature_id);
        js.flag("variable_properties", cli.variable_properties);
        js.num("mu_cap", cli.mu_cap);
        js.num("t_rad_ambient", cli.t_rad_ambient);
        js.flag("t_rad_ambient_set", cli.t_rad_ambient_set);
        if (cli.variable_properties)
        {
            js.str("variable_properties_curves", variable_props_curve_note);
            js.num("mu_raw_min_end", rheology_stats.mu_raw_min);
            js.num("mu_raw_max_end", rheology_stats.mu_raw_max);
            js.num("mu_raw_mean_end", rheology_stats.mu_raw_mean);
            js.num("mu_used_min_end", rheology_stats.mu_used_min);
            js.num("mu_used_max_end", rheology_stats.mu_used_max);
            js.num("mu_used_mean_end", rheology_stats.mu_used_mean);
            js.num("mu_mean_end", rheology_stats.mu_used_mean); /**< Doc180: used field, not raw curve */
            js.num("mu_max_end", rheology_stats.mu_used_max);
            js.integer("n_mu_capped_end", static_cast<long long>(rheology_stats.n_mu_capped));
            js.num("vol_frac_mu_capped_end", rheology_stats.vol_frac_mu_capped);
            js.str("mu_stats_stage_end", rheology_stats.mu_stats_stage);
            js.num("cp_mean_end", rheology_stats.cp_mean);
            js.num("k_mean_end", rheology_stats.k_mean);
            js.flag("variable_cp_in_thermal_operators", true);
            js.str("variable_cp_in_thermal_operators_note",
                   "Q/diffusion/BC use per-particle OphelieHeatCapacity; energy ledger uses enthalpy integral of cp(T)");
            js.flag("variable_k_in_laplace", true);
            js.flag("variable_mu_in_viscous_mean", !use_local_particle_mu);
            js.flag("variable_mu_in_viscous_pairwise", use_local_particle_mu);
            js.str("viscous_mu_mode", use_local_particle_mu ? "particle_pairwise_harmonic" : "volume_mean");
            js.flag("enthalpy_energy_ledger", true);
            if (!cli.material_curve_dir.empty())
            {
                js.str("material_curve_dir", cli.material_curve_dir);
                js.str("material_k_path", thesis_curve_report.k_path);
                js.str("material_cp_path", thesis_curve_report.cp_path);
                js.str("material_mu_path", thesis_curve_report.mu_path);
                js.str("material_interpolation", thesis_curve_report.interpolation);
                js.str("material_endpoint_policy", thesis_curve_report.endpoint_policy);
                js.num("material_k_1473", thesis_curve_report.k_1473);
                js.num("material_cp_1473", thesis_curve_report.cp_1473);
                js.num("material_mu_1473", thesis_curve_report.mu_1473);
                js.num("material_k_303", thesis_curve_report.k_303);
                js.num("material_cp_303", thesis_curve_report.cp_303);
                js.num("material_mu_303", thesis_curve_report.mu_303);
            }
            if (local_mu_evidence_recorded)
            {
                js.flag("local_mu_harmonic_ok", local_mu_evidence.harmonic_mean_ok);
                js.num("local_mu_harmonic_max_rel_err", local_mu_evidence.harmonic_max_rel_err);
                js.flag("local_mu_spatially_active", local_mu_evidence.local_mu_spatially_active);
                js.flag("local_mu_force_contrast_ok", local_mu_evidence.force_contrast_supports_local_mu);
                js.num("local_mu_Fvisc_rms", local_mu_evidence.viscous_force_rms);
                js.num("local_mu_Fvisc_rms_high_mu", local_mu_evidence.viscous_force_rms_high_mu);
                js.num("local_mu_Fvisc_rms_low_mu", local_mu_evidence.viscous_force_rms_low_mu);
            }
        }
        js.str("viscosity_mode", cli.viscosity_mode);
        js.str("implicit_visc_schedule", cli.implicit_visc_schedule);
        js.str("viscosity_mode_note",
               use_implicit_split
                   ? "Doc175 implicit_split: InteractionSplit+DampingPairwise on CPU; ViscousForceCK off; "
                     "dt_visc estimate logged only; host sync each acoustic substep"
                   : "explicit ViscousForceCK + dt_visc CFL on advection step");
        js.flag("implicit_visc_cpu_split", use_implicit_split);
        js.integer("n_implicit_visc_exec", static_cast<long long>(n_implicit_visc_exec));
        js.num("implicit_visc_wall_s", implicit_visc_wall_s);
        js.num("dt_visc_explicit_est_last", last_dt_visc_est);
        js.str("active_limiter_last", last_active_limiter);
        js.num("marangoni_dgamma_dt", cli.marangoni_dgamma_dt);
        js.flag("marangoni_active", std::abs(cli.marangoni_dgamma_dt) > TinyReal);
        js.str("marangoni_method",
               std::abs(cli.marangoni_dgamma_dt) > TinyReal ? "axisymmetric_radial_bin_FD" : "off");
        js.str("marangoni_note",
               std::abs(cli.marangoni_dgamma_dt) > TinyReal
                   ? "flat-top axisymmetric approx: ∇_s T≈(dT/dr)e_r from free-surface radial-bin FD; "
                     "ForcePrior += (dγ/dT)(dT/dr) A e_r; not full 3D surface SPH gradient MMS"
                   : "off");
        if (marangoni_verify_ran)
        {
            js.flag("marangoni_checked", marangoni_verify_result.passed);
            js.flag("marangoni_direction_ok", marangoni_verify_result.direction_ok);
            js.flag("marangoni_area_integral_ok", marangoni_verify_result.area_integral_ok);
            js.num("marangoni_verify_rel_err_integral", marangoni_verify_result.rel_err_integral);
            js.num("marangoni_verify_sum_Fr", marangoni_verify_result.sum_Fr);
            js.num("marangoni_verify_sum_Fr_theory", marangoni_verify_result.sum_Fr_theory);
            js.num("marangoni_verify_sum_area", marangoni_verify_result.sum_area);
            js.str("marangoni_verify_message", marangoni_verify_result.message);
        }
        else
        {
            js.gate("marangoni_checked", Gate::Null,
                    std::abs(cli.marangoni_dgamma_dt) > TinyReal
                        ? "marangoni_on_but_verify_skipped"
                        : "marangoni_off");
        }
        if (last_marangoni_apply.applied)
        {
            js.num("marangoni_last_sum_Fr", last_marangoni_apply.sum_Fr);
            js.num("marangoni_last_sum_area", last_marangoni_apply.sum_area);
            js.num("marangoni_last_mean_dT_dr", last_marangoni_apply.mean_dT_dr);
            js.integer("marangoni_last_n_free", static_cast<long long>(last_marangoni_apply.n_free));
        }
        js.str("mode", cli.mode);
        js.integer("real_bytes", static_cast<long long>(sizeof(Real)));
        js.str("precision_note",
               sizeof(Real) == 4
                   ? "Real=float32 (GPU fields, kernels, Temperature, OphelieThermalDeltaT); host reductions/ledgers in "
                     "double. DeltaT is the authoritative thermal variable (T = t_ref + DeltaT). Not a float64 build."
                   : "Real=float64 throughout; host reductions in double.");
        js.str("argv", argv_joined);
        js.num("effective_cp", cli.cp);
        js.num("effective_k", cli.k);
        js.str("power_mode", cli.power_mode);
        js.flag("diagnostic_consistent", diagnostic_consistent);
        js.str("diagnostic_consistent_note",
               "1 = finite fields, valid dt, mass conserved, EOS field consistent when sampled. Production end state "
               "(U,K,T,t,pi,dH,mass) is read before the diagnostic force evaluation; the only post-end state write "
               "is the diagnostic Boussinesq force refresh (ForcePrior) (and, only with --r2-fixed-spatial-T, one T re-application "
               "at the current position), which does not advance time/velocity.");
        js.num("time_production_end_s", time_production_end);
        js.num("dp", dp);
        js.num("R", R);
        js.num("H", H);
        js.integer("perturbation_momentum", 1);
        js.integer("acoustic_dens_mode", cli.acoustic_dens_mode);
        js.str("acoustic_dens_mode_note",
               cli.acoustic_dens_mode == 3
                   ? "mode3: Doc201 omit pressure DissipativeUJump dens (inner+wall); continuity from 2nd-half u only; "
                     "for stratified p' hydro where mode1 DissipativeUJump(p'_i-p'_j)≠0 at rest"
                   : (cli.acoustic_dens_mode == 2
                          ? "mode2: dens as mode1 (C_ext/C_ref); wall momentum Adami omits C_ext (Doc186 static-balance probe)"
                          : (cli.acoustic_dens_mode == 1
                                 ? "mode1: DissipativeUJump on p'; wall dens uses C_ext-C_ref"
                                 : "mode0: legacy Adami dens on momentum pressure")));
        js.integer("reference_body_acceleration", cli.cref_buoyancy ? 1 : 0);
        js.str("reference_body_acceleration_note",
               cli.cref_buoyancy
                   ? "--cref-buoyancy: ReferenceBodyAcceleration := Boussinesq force/mass each advection step "
                     "(AcousticDensMode dens jump C_ext-C_ref cancels buoyancy; g still NOT in a_ref)"
                   : "zero (gravity is in BackgroundPressure, not in ForcePrior; buoyancy then makes C_ext≠C_ref)");
        js.flag("cref_buoyancy", cli.cref_buoyancy);
        js.flag("doc200_strat_init_vs_tmean", cli.doc200_strat_init_vs_tmean);
        js.flag("doc201_keep_strat_owns_when_thaw", cli.doc201_keep_strat_owns_when_thaw);
        js.flag("strat_p_prime_owns_buoyancy", strat_p_prime_owns_buoyancy);
        if (strat_p_prime_owns_buoyancy)
        {
            js.str("strat_p_prime_owns_buoyancy_note",
                   "Doc201: Doc200 strat p'/ρ already balances gβ(T-Tmean); Boussinesq ForcePrior kept at 0 to "
                   "avoid double-count (discrete SPH ∇p' does not cancel particle-local ForcePrior)");
        }
        js.flag("uniform_buoyancy_balance_init", cli.uniform_buoyancy_balance_init);
        js.flag("absorb_mean_buoyancy_bg", cli.absorb_mean_buoyancy_bg);
        if (cli.absorb_mean_buoyancy_bg)
        {
            js.num("absorb_t_mean_last", last_t_mean_buoy);
            js.num("absorb_force_tref_last", last_force_tref);
            js.num("absorb_eos_tref", cli.t_ref);
            js.str("absorb_mean_buoyancy_bg_note",
                   "mean buoyancy in BackgroundPressure; Boussinesq force Tref:=T_mean; EOS/material Tref unchanged");
        }
        // ---- Doc133 D1 drift probe / E Q shape / VTP ----
        js.flag("drift_probe", cli.drift_probe);
        js.flag("probe_freeze_velocity", cli.probe_freeze_velocity);
        js.flag("probe_cref_buoyancy", cli.probe_cref_buoyancy);
        // ---- Doc136 P0 density half-step audit ----
        js.flag("density_increment_audit", cli.density_increment_audit);
        if (cli.density_increment_audit)
        {
            js.integer("density_audit_stride", cli.density_audit_stride);
            js.integer("density_audit_max_rows", cli.density_audit_max_rows);
            js.integer("density_audit_n_substep", static_cast<long long>(density_audit.n_substep));
            js.integer("density_audit_n_rows", static_cast<long long>(density_audit.n_rows_written));
            js.num("density_audit_cum_d_act_mean", density_audit.cum_d_act_mean);
            js.num("density_audit_cum_d_exp_mean", density_audit.cum_d_exp_mean);
            js.num("density_audit_cum_R_update_mean", density_audit.cum_r_update_mean);
            js.num("density_audit_sum_lost_signed", density_audit.sum_lost_signed);
            js.num("density_audit_sum_lost_abs", density_audit.sum_lost_abs);
            js.integer("density_audit_n_zero_store_events",
                       static_cast<long long>(density_audit.n_zero_store_events));
            js.integer("density_audit_n_nonzero_exp_events",
                       static_cast<long long>(density_audit.n_nonzero_exp_events));
            js.integer("density_audit_n_particle_steps",
                       static_cast<long long>(density_audit.n_particle_steps));
            const double nz = static_cast<double>(std::max<size_t>(density_audit.n_nonzero_exp_events, 1));
            js.num("density_audit_zero_store_fraction_of_nonzero_exp",
                   static_cast<double>(density_audit.n_zero_store_events) / nz);
            js.num("density_audit_max_abs_R_update", density_audit.max_abs_r_update);
            js.integer("density_audit_max_R_update_OriginalID",
                       static_cast<long long>(density_audit.max_r_update_oid));
            js.str("density_audit_global_csv", density_audit.global_path);
            js.str("density_audit_id_csv", density_audit.id_path);
            js.str("density_audit_note", density_audit.note);
            js.str("density_update_budget_closed",
                   (std::abs(density_audit.cum_r_update_mean) < 1e-6 * std::max(1.0, std::abs(density_audit.cum_d_exp_mean)))
                       ? "likely_closed_within_float_noise"
                       : "open_or_large_R_update");
            js.gate("roundoff_loss_quantified",
                    density_audit.n_zero_store_events > 0 ? Gate::Yes : Gate::No,
                    density_audit.n_zero_store_events > 0
                        ? "nonzero expected half-step increments with zero stored Density change observed"
                        : "no zero-store events in sampled window");
        }
        js.num("cum_drho_dissipative_kg_m3", cum_drho_diss);
        js.num("cum_drho_total_kg_m3", cum_drho_total);
        js.num("cum_drho_continuity_kg_m3", cum_drho_total - cum_drho_diss);
        js.num("cum_drho_band_total_kg_m3", cum_drho_band_total);
        js.str("cum_drho_ledger_note",
               "Doc136: cum_drho_total = sum(0.5*dt*mean(rate_before)+0.5*dt*mean(rate_after2)) matching AcousticStep "
               "writes; cum_drho_dissipative = sum(0.5*dt*mean(rate_before)); continuity column = total-diss = 2nd-write "
               "term. Legacy Doc133 used full-dt after each half (over-count).");
        js.num("E_side_J", e_side_J);
        js.num("E_bottom_J", e_bottom_J);
        js.num("E_top_J", e_top_J);
        js.num("P_side_W_end", last_P_side);
        js.num("P_bottom_W_end", last_P_bottom);
        js.num("P_top_conv_W_end", last_P_top_conv);
        js.num("P_top_rad_W_end", last_P_top_rad);
        js.num("P_boundary_total_W_end", last_P_bound);
        js.str("power_host_formula_end", last_power_host_formula);
        js.flag("MONITOR_MATCHES_APPLIED_BC", true);
        js.str("P_boundary_avg_preferred",
               "Doc208: use (E_boundary(t1)-E_boundary(t0))/(t1-t0); do not arithmetic-mean unequally spaced P_* columns");
        // Face-sum vs total consistency at end sample (host audit).
        {
            const double face_sum =
                last_P_side + last_P_bottom + last_P_top_conv + last_P_top_rad;
            const double denom = std::abs(last_P_bound) + 1.0e-30;
            js.num("P_face_sum_W_end", face_sum);
            js.num("P_face_sum_minus_total_rel", std::abs(face_sum - last_P_bound) / denom);
        }
        {
            const double e_face = e_side_J + e_bottom_J + e_top_J;
            const double denom = std::abs(ledger.e_boundary) + 1.0e-30;
            js.num("E_face_sum_J", e_face);
            js.num("E_face_sum_minus_E_boundary_rel", std::abs(e_face - ledger.e_boundary) / denom);
        }
        js.str("q_shape", cli.q_shape);
        js.str("q_from_rz_csv", cli.q_from_rz_csv);
        js.str("q_from_em_particles_csv", cli.q_from_em_particles_csv);
        js.str("em_q_delivery", cli.em_q_delivery);
        js.num("em_q_support_radius_m", cli.em_q_support_radius);
        const bool em_rz = (cli.q_shape == "em_rz_csv" || cli.q_shape == "em_particles_csv");
        js.integer("source_study_not_em_verified", em_rz ? 0 : 1);
        js.flag("effective_source_study", q_active && cli.q_shape == "radial_skin");
        js.flag("fixed_space_em_q_shape", q_active && em_rz);
        js.integer("em_q_live_solve", em_coupling_armed ? 1 : 0);
        js.str("em_update_mode", cli.em_update_mode);
        js.str("em_schedule", cli.em_schedule);
        js.str("init_T_rz_csv", cli.init_T_rz_csv);
        js.str("init_T_particle_csv", cli.init_T_particle_csv);
        js.str("init_T_source_note", cli.init_T_source_note);
        js.str("em_replay_dir", cli.em_replay_dir);
        js.num("em_update_origin_s", cli.em_update_origin);
        js.num("em_update_interval_s", cli.em_update_interval);
        js.integer("max_em_updates", static_cast<long long>(max_em_updates_eff));
        js.flag("em_feedback_capped", static_cast<long>(n_em_updates) >= max_em_updates_eff);
        js.integer("em_updates", static_cast<long long>(n_em_updates));
        js.flag("em_coupling_ok", em_coupling_ok);
        js.str("em_failure_reason", em_failure_reason);
        js.num("em_last_sigma_mean", last_sigma_mean);
        js.num("em_last_delta_sigma_mean", last_delta_sigma);
        js.num("em_last_shape_L2_vs_prev", last_shape_L2);
        js.num("em_diag_T_perturb_K", cli.em_diag_T_perturb_K);
        js.str("em_last_Q_csv", last_em_Q_csv);
        js.str("em_sigma_reduction", em_coupling_armed ? cli.em_sigma_reduction : "n/a");
        js.str("em_particles_csv", cli.em_particles_csv);
        js.num("em_dp", cli.em_dp);
        js.flag("em_end_flush", cli.em_end_flush);
        js.flag("update_scheduler_no_end_flush", !cli.em_end_flush);
        js.integer("fluid_steps_after_last_em", static_cast<long long>(n_fluid_steps_after_last_em));
        js.integer("thermal_steps_after_last_em", static_cast<long long>(n_thermal_steps_after_last_em));
        js.flag("post_update_flow_exercised", n_fluid_steps_after_last_em > 0);
        js.str("material_mode", cli.material_mode.empty() ? "unset" : cli.material_mode);
        if (!em_update_scheduled_times.empty())
        {
            js.num("em_update_scheduled_t0", em_update_scheduled_times[0]);
            js.num("em_update_actual_t0", em_update_actual_times[0]);
        }
        if (em_update_scheduled_times.size() > 1)
        {
            js.num("em_update_scheduled_t1", em_update_scheduled_times[1]);
            js.num("em_update_actual_t1", em_update_actual_times[1]);
        }
        if (thesis_curve_report.cep_anchored)
        {
            js.flag("material_anchor_coefficients_checked", true);
            js.num("material_s_mu", thesis_curve_report.s_mu);
            js.num("material_s_cp", thesis_curve_report.s_cp);
            js.num("material_s_k", thesis_curve_report.s_k);
            js.num("material_mu_1473_raw", thesis_curve_report.mu_1473_raw);
            js.num("material_cp_1473_raw", thesis_curve_report.cp_1473_raw);
            js.num("material_k_1473_raw", thesis_curve_report.k_1473_raw);
            js.num("material_mu_1473", thesis_curve_report.mu_1473);
            js.num("material_cp_1473", thesis_curve_report.cp_1473);
            js.num("material_k_1473", thesis_curve_report.k_1473);
        }
        js.str("q_shape_note",
               cli.q_shape == "uniform"
                   ? "uniform s_i=1"
                   : (cli.q_shape == "radial_skin"
                          ? "radial_skin: s_i=(r/R)^2 >= 0 re-evaluated from CURRENT positions each advection step; "
                            "P conserved by discrete normalisation. Effective source study only: NOT EM-verified."
                          : (cli.q_shape == "em_particles_csv"
                                 ? (em_coupling_armed
                                        ? "em_particles_csv + Doc156 spatial-σ; direct EM→fluid kernel Q; discrete-renorm; "
                                          "strict em_result.json gate."
                                        : "em_particles_csv: Shepard sample from EM particle Q cloud; discrete-renorm.")
                                 : (em_coupling_armed
                                        ? (cli.em_sigma_reduction == "spatial_field"
                                               ? "em_rz_csv + Doc156 spatial-σ subprocess T→σ→L1 Picard; discrete-renorm; "
                                                 "fixed EM lattice; strict em_result.json gate."
                                               : "em_rz_csv + Doc151 subprocess T→σ→L1 Picard refresh; discrete-renorm to "
                                                 "P_target; σ reduced to volume-mean for existing Picard --sigma=.")
                                        : "em_rz_csv: s_i=max(Q_interp(r,z),0) from prior EM r-z dump; fixed Eulerian; "
                                          "discrete-renorm to P_target each reapply. NOT a live EM solve in this binary; "
                                          "provenance in q_from_rz_csv."))));
        if (cli.q_shape == "em_rz_csv" && rz_q_ptr != nullptr)
        {
            js.num("q_rz_table_min", rz_q_ptr->q_table_min);
            js.num("q_rz_table_max", rz_q_ptr->q_table_max);
            js.integer("q_rz_nr", static_cast<long long>(rz_q_ptr->nr));
            js.integer("q_rz_nz", static_cast<long long>(rz_q_ptr->nz));
        }
        if (cli.q_shape == "em_particles_csv" && em_q_ptr != nullptr)
        {
            js.integer("em_q_n_nodes", static_cast<long long>(em_q_ptr->n_nodes));
            js.num("em_q_P_particles_W", em_q_ptr->P_particles);
            js.num("em_q_min", em_q_ptr->q_min);
            js.num("em_q_max", em_q_ptr->q_max);
        }
        js.integer("q_shape_reapply_count", static_cast<long long>(n_q_shape_reapply));
        js.num("state_record_interval_s", cli.state_record_interval);
        js.integer("vtp_files_written", static_cast<long long>(n_vtp_written));
        js.str("vtp_output_dir", std::filesystem::absolute("output").string());
        js.str("vtp_fields", "Velocity,Pressure,Density,BackgroundPressure,Temperature,Mass,VolumetricMeasure,"
                             "DensityChangeRate,OphelieThermalSourceQ(volumetric Q W/m3)");
        js.integer("gravity_force_ck", 0);
        js.integer("boussinesq", 1);
        js.str("density_path", cli.use_density_summation ? "internal_summation" : "continuity");
        js.flag("reapply_hydro_rho", reapply_flag);
        js.integer("analytic_density_rewrite_per_step", 0);
        js.flag("kernel_correction", cli.kernel_correction);
        js.str("top_bc", "flat_free_slip_viscous_omit");
        js.flag("doc201_b2_live_top_shear", cli.doc201_b2_live_top_shear);
        js.flag("doc202_b2_ablate_side_acoustic", cli.doc202_b2_ablate_side_acoustic);
        js.flag("doc202_b2_ablate_lid_acoustic", cli.doc202_b2_ablate_lid_acoustic);
        js.flag("doc202_b2_ablate_lid_acoustic_1st", cli.doc202_b2_ablate_lid_acoustic_1st);
        js.flag("doc202_b2_ablate_lid_acoustic_2nd", cli.doc202_b2_ablate_lid_acoustic_2nd);
        js.flag("doc202_b2_lid_continuity_silent", cli.doc202_b2_lid_continuity_silent);
        js.flag("doc202_b2_lid_omit_pjump", cli.doc202_b2_lid_omit_pjump);
        if (cli.doc201_b2_live_top_shear)
        {
            js.num("doc201_b2_U0", cli.doc201_b2_U0);
            js.num("doc201_b2_t_s", b2_t);
            js.num("doc201_b2_nu", b2_nu);
            js.str("doc201_b2_acoustic_lid", "wall_closed_live_positions");
            js.str("doc201_b2_viscous_top", "omit_noslip_open_cup");
            js.num("doc201_b2_core_L2_rel", b2_core.L2_rel);
            js.num("doc201_b2_core_amp_rel_err", b2_core.amp_rel_err);
            js.num("doc201_b2_core_amp_num", b2_core.amp_num);
            js.num("doc201_b2_core_amp_ref", b2_core.amp_ref);
            js.num("doc201_b2_core_uz_rms_over_U0", b2_core.uz_rms_over_U0);
            js.num("doc201_b2_core_uy_rms", b2_core.uy_rms);
            js.num("doc201_b2_core_disp_rms", b2_core.disp_rms);
            js.num("doc201_b2_core_rho_rel_drift", b2_core.rho_rel_drift);
            js.num("doc201_b2_core_p_prime_rms", b2_core.p_prime_rms);
            js.integer("doc201_b2_core_n_interior", static_cast<long long>(b2_core.n_interior));
            js.num("doc201_b2_annulus_L2_rel", b2_annulus.L2_rel);
            js.num("doc201_b2_annulus_uz_rms_over_U0", b2_annulus.uz_rms_over_U0);
            js.integer("doc201_b2_annulus_n_interior", static_cast<long long>(b2_annulus.n_interior));
            js.str("doc201_b2_profile_csv",
                   (cli.output_dir.empty() ? std::string(".") : cli.output_dir) + "/b2_profile_end.csv");
            js.str("doc201_b2_note",
                   "core=r<0.5R; annulus=r<R-3.5dp excludes side band; infinite-plate analytic polluted by "
                   "cylinder side no-slip — annulus vs core shows contamination");
        }
        js.num("rho0", cli.rho0);
        js.num("g", cli.gravity_g);
        js.num("beta", beta_eff);
        js.num("t_ref", cli.t_ref);
        js.num("mu", cli.mu);
        js.num("c0", cli.c0);
        js.num("a_strat", a_eff);
        js.num("dT_radial", is_R2 ? static_cast<double>(r2_dT_eff) : (is_T2 ? cli.dT_radial : 0.0));
        js.flag("r2_mean_removed", is_R2 && !cli.r2_keep_mean);
        js.num("r2_s_mean", r2_s_mean);
        // ---- Doc131 D2 R2 driver protocol / sign gate ----
        js.flag("r2_reverse", is_R2 && cli.r2_reverse);
        js.flag("r2_fixed_spatial_T", is_R2 && cli.r2_fixed_spatial_T);
        js.flag("fixed_spatial_T", cli.fixed_spatial_T);
        js.flag("fixed_spatial_T_demean", cli.fixed_spatial_T_demean);
        js.flag("fixed_T_after_settle", cli.fixed_T_after_settle);
        js.flag("zero_velocity_after_settle", cli.zero_velocity_after_settle);
        js.flag("quiet_reset_mass", cli.quiet_reset_mass);
        js.flag("fixed_T_live_at_end", fixed_T_live);
        js.num("U_rms_after_quiet", s_quiet.U_rms);
        js.str("fixed_spatial_T_csv", cli.fixed_spatial_T ? cli.init_T_rz_csv : "");
        js.str("r2_T_driver_mode",
               !is_R2 ? "n/a"
                      : (cli.r2_fixed_spatial_T
                             ? "fixed_spatial_T: T(r) re-applied from the current position every advection step "
                               "(diagnostic; isolates flow response)"
                             : "lagrangian_carry: T(r) written once at init, then carried with the particles (no "
                               "re-application, no conduction/Q/BC in R2); NOT a fixed-space driver"));
        js.integer("r2_fixed_T_reapply_count", static_cast<long long>(n_fixed_T_reapply));
        if (is_R2)
        {
            js.num("r2_dT_sign", r2_dT_sign);
            js.flag("r2_null_driver", r2_null_driver);
            js.str("r2_gate_kind", r2_gate_kind);
            js.flag("r2_sign_uz_ok", r2_sign_uz_ok);
            js.flag("r2_sign_ur_ok", r2_sign_ur_ok);
            js.flag("r2_signal", r2_signal);
        }
        js.num("acoustic_cfl", cli.acoustic_cfl);
        js.num("advection_cfl", cli.advection_cfl);
        js.num("settle_time_s", settle_eff);
        js.num("observe_time_s", cli.end_time);
        // ---- Doc131 D4 boundary band / escape report (production end) ----
        js.num("d4_band_thickness_m", dp);
        js.integer("d4_n_top_band", static_cast<long long>(bb_end.n_top));
        js.num("d4_top_band_uz_rms", bb_end.top_uz_rms);
        js.num("d4_top_band_uz_maxabs", bb_end.top_uz_maxabs);
        js.num("d4_top_band_uz_mean", bb_end.top_uz_mean);
        js.integer("d4_n_side_band", static_cast<long long>(bb_end.n_side));
        js.num("d4_side_band_ur_rms", bb_end.side_ur_rms);
        js.num("d4_side_band_ur_maxabs", bb_end.side_ur_maxabs);
        js.integer("d4_n_bottom_band", static_cast<long long>(bb_end.n_bottom));
        js.num("d4_bottom_band_uz_rms", bb_end.bottom_uz_rms);
        js.num("d4_bottom_band_uz_maxabs", bb_end.bottom_uz_maxabs);
        js.integer("d4_n_outside_strict", static_cast<long long>(bb_end.n_outside_strict));
        js.integer("d4_n_outside_half_dp", static_cast<long long>(bb_end.n_outside_half_dp));
        js.num("d4_z_max", bb_end.z_max);
        js.num("d4_z_min", bb_end.z_min);
        js.num("d4_r_max", bb_end.r_max);
        js.num("perturb_u", is_R3 ? cli.perturb_u : 0.0);
        js.integer("n_fluid", static_cast<long long>(particles.TotalRealParticles()));
        js.num("U_rms", s_end.U_rms);
        js.num("Uz_rms", s_end.Uz_rms);
        js.num("Ur_rms", s_end.Ur_rms);
        js.num("Utheta_rms", s_end.Utheta_rms);
        js.num("U_max", s_end.U_max);
        js.num("K", s_end.K);
        js.num("U_mean_x", s_end.U_mean[0]);
        js.num("U_mean_y", s_end.U_mean[1]);
        js.num("U_mean_z", s_end.U_mean[2]);
        js.num("U_rms_settle", s_settle.U_rms);
        js.num("U_rms_perturb", s_perturb.U_rms);
        js.num("U_growth_end_over_settle", u_growth);
        js.num("decay_ratio_end_over_perturb", decay_ratio);
        js.num("circulation_ur_proxy", s_end.circ_ur_proxy);
        js.num("circulation_uz_proxy", s_end.circ_uz_proxy);
        js.num("buoyancy_accel_rms", buoy_a_rms_end);
        js.flag("buoyancy_present", buoyancy_present);
        js.num("boussinesq_force_rel_err", force_end.rel_err);
        js.str("boussinesq_force_rel_err_note",
               "legacy |sum F - sum Fref|/|sum Fref|; not the Doc189 pass criterion when net force is near zero");
        js.flag("doc189_mechanism", cli.doc189_mechanism);
        js.num("doc189_mechanism_amp_K", cli.doc189_mechanism_amp);
        js.integer("doc189_buoyancy_param_version", doc189_param_version);
        js.num("doc189_buoyancy_snapshot_t_s", doc189_snapshot_t);
        js.num("doc189_buoyancy_snapshot_T_mean_K", last_t_mean_buoy);
        js.num("doc189_buoyancy_snapshot_force_tref_K", last_force_tref);
        if (cli.doc189_mechanism)
            js.str("energy_gate", "NOT_APPLICABLE");
        if (doc189_gate_n > 0)
        {
            js.integer("doc189_force_gate_samples", doc189_gate_n);
            js.num("doc189_force_E_L2_worst", doc189_worst_L2);
            js.num("doc189_force_E_net_worst", doc189_worst_Enet);
            js.num("doc189_force_max_abs_N_worst", doc189_worst_abs);
            js.num("doc189_force_ffloor_N", doc189_last_gate.ffloor_N);
            js.num("doc189_tref_requested_last", doc189_last_gate.tref_requested_K);
            js.num("doc189_tref_inferred_last", doc189_last_gate.tref_inferred_K);
            js.flag("dynamic_force_reference_passed", doc189_dynamic_ok);
            js.flag("buoyancy_field_consistency_passed", doc189_algebraic_ok);
            // Doc192: ForcePrior double-count is not implied by the dedicated-field algebraic gate.
            js.gate("forceprior_no_double_count_passed", Gate::Null,
                    "NOT_RUN_in_production: use --doc189-dynamic-force-repro ForcePrior independent check");
            js.str("forceprior_no_double_count_note",
                   "production samples only the dedicated Boussinesq field vs formula; "
                   "ForcePrior seed/delta/repeat gate lives in DYNAMIC_FORCE_REPRO.json");
        }
        else
        {
            js.gate("dynamic_force_reference_passed", Gate::Null, "NOT_RUN: no in-loop buoyancy sample");
            js.gate("buoyancy_field_consistency_passed", Gate::Null, "NOT_RUN: no in-loop buoyancy sample");
            js.gate("forceprior_no_double_count_passed", Gate::Null, "NOT_RUN: no in-loop buoyancy sample");
        }
        js.num("T_mean", s_end.T_mean);
        js.num("T_min", s_end.T_min);
        js.num("T_max", s_end.T_max);
        js.num("T_mean_init", s_init.T_mean);
        js.num("rho_mean", s_end.rho_mean);
        js.num("rho_min", s_end.rho_min);
        js.num("rho_max", s_end.rho_max);
        js.num("p_prime_rms", s_end.p_prime_rms);
        // ---- A2 EOS semantics ----
        js.flag("eos_initialization_ok", eos_initialization_ok);
        js.flag("aux_density_range_ok", aux_density_range_ok);
        js.flag("eos_identity_ok", eos_identity_ok_legacy);
        js.str("eos_identity_ok_note",
               "LEGACY alias = eos_initialization_ok AND aux_density_range_ok. It is NOT field consistency "
               "(see eos_field_consistency_ok).");
        if (pressure_field_available)
            js.flag("eos_field_consistency_ok", eos_field_consistency_ok);
        else
            js.gate("eos_field_consistency_ok", Gate::Null,
                    "no acoustic substep sampled (frozen flow or aborted run): Pressure is not refreshed");
        js.flag("pressure_field_available", pressure_field_available);
        js.str("pressure_storage_semantics",
               "total_pressure_p=EOS(rho_half)+BackgroundPressure (p_prime=Pressure-BackgroundPressure; "
               "PerturbationMomentum=1 makes the momentum use p_prime)");
        js.str("pressure_sample_phase",
               pressure_field_available ? "after_acoustic_1st_half_of_last_production_substep (read-only)"
                                        : "not_sampled");
        js.num("eos_init_residual", eos_res_init);
        js.num("eos_init_tol", eos_tol);
        js.num("pi_mean", pi_density_end.mean);
        js.num("pi_rms_demean", pi_density_end.rms_demean);
        js.num("pi_min", pi_density_end.min);
        js.num("pi_max", pi_density_end.max);
        js.str("pi_source", "EOS(Density) at production_end (p0*(rho/rho0-1), Real arithmetic)");
        if (pressure_field_available)
        {
            js.num("pressure_sample_time_s", pressure_sample.time);
            js.num("pressure_field_pi_mean", pressure_sample.pi_from_pressure.mean);
            js.num("pressure_field_pi_rms_demean", pressure_sample.pi_from_pressure.rms_demean);
            js.num("pressure_field_pi_min", pressure_sample.pi_from_pressure.min);
            js.num("pressure_field_pi_max", pressure_sample.pi_from_pressure.max);
            js.num("density_eos_at_sample_pi_mean", pressure_sample.pi_from_density.mean);
            js.num("density_eos_at_sample_pi_rms_demean", pressure_sample.pi_from_density.rms_demean);
            js.num("density_eos_at_sample_pi_min", pressure_sample.pi_from_density.min);
            js.num("density_eos_at_sample_pi_max", pressure_sample.pi_from_density.max);
            js.num("eos_field_residual_rms", pressure_sample.resid_rms);
            js.num("eos_field_residual_max", pressure_sample.resid_max);
            js.num("eos_field_residual_tol", pressure_sample.tol);
            js.num("eos_field_scale_rms_pi", pressure_sample.eos_rms0);
            js.num("eos_field_pressure_total_max_abs", pressure_sample.p_total_max_abs);
            js.num("eos_field_background_max_abs", pressure_sample.pbg_max_abs);
        }
        js.num("mass_sum_init", mass_sum_init);
        js.num("mass_sum_end", mass_sum_end);
        js.flag("mass_conservation_pass", mass_conservation_pass);
        js.integer("n_bad_dt", static_cast<long long>(n_bad_dt));
        js.flag("no_runaway", no_runaway);
        // ---- A3/B config + power + area ----
        js.num("k", cli.k);
        js.num("cp", cli.cp);
        js.num("h_side", cli.h_side);
        js.num("h_bottom", cli.h_bottom);
        js.num("h_top", cli.h_top);
        js.num("t_cool", cli.t_cool);
        js.num("t_ambient", cli.t_ambient);
        js.num("t_rad_ambient", cli.t_rad_ambient);
        js.num("emissivity", cli.emissivity);
        js.str("thermal_param_status",
               cli.literature_id.empty()
                   ? "h_side,h_bottom,h_top,t_cool: CODE_PROVISIONAL unless set via Doc139 literature config"
                   : "Doc139 literature_id active; see config assumptions/borrowed fields");
        js.flag("freeze_flow", freeze_flow);
        js.flag("thermal_active", thermal_active);
        js.flag("legacy_rad_clamp", legacy_path);
        js.flag("thermal_fields_follow_particle_sort", !cli.legacy_thermal_sort);
        js.str("conduction_path", conduction_path);
        js.str("thermal_boundary_path", boundary_path);
        js.flag("conduction_active", conduction_active);
        js.flag("bc_active", bc_active);
        js.str("thermal_authoritative_variable", "OphelieThermalDeltaT (Temperature = t_ref + DeltaT)");
        js.num("P_target_W", power.P_target);
        js.num("P_applied_W", power.P_applied);
        js.num("P_gap_rel", power.gap_rel);
        js.num("V_geom_m3", power.V_geom);
        js.num("sum_V_Q_m3", power.sum_V_Q);
        js.num("mass_sum_kg", power.mass_sum);
        js.num("sum_V_Q_over_V_geom", power.ratio_sumVQ_over_Vgeom);
        js.num("N_Q_m3", power.N_Q);
        js.num("power_norm_coeff", power.norm_coeff);
        js.str("power_norm_coeff_meaning",
               "P_target/N_Q after the delivered shape is sampled on SPH volumes; not raw EM joule to 50 kW");
        js.str("em_power_ledger_csv", em_power_ledger_csv);
        js.raw("em_power_stages", power_stages_json.str());
        js.num("q_min_W_per_m3", power.q_min);
        js.num("q_max_W_per_m3", power.q_max);
        js.num("q_power_W", cli.q_power);
        js.num("q_vol_W_per_m3", power.q_max);
        js.flag("area_normalize", cli.area_normalize);
        js.str("thermal_wall_model", cli.thermal_wall_model);
        js.flag("RECONSTRUCTED_WALL_TEMPERATURE", use_wall_T);
        js.str("MONITOR_FIX_DOC", "Doc208_A: write_monitor_row uses hostOphelieThermalAppliedBcHeatLossPowers");
        js.flag("NOT_FRENCH_VALIDATED", true);
        js.num("area_shell_thickness_m", area.shell);
        js.integer("area_n_corner_primary_assigned", static_cast<long long>(area.n_corner));
        js.integer("n_top", static_cast<long long>(area.top.n));
        js.num("A_top_raw", area.top.raw);
        js.num("A_top_target", area.top.target);
        js.num("A_top_rel_err_raw", area.top.rel_err_raw);
        js.num("A_top_after", area.top.after);
        js.num("A_top_rel_err", area.top.rel_err);
        js.num("A_top_scale", area.top.scale);
        js.num("A_top_ai_min", area.top.a_i_min);
        js.num("A_top_ai_max", area.top.a_i_max);
        js.integer("n_side", static_cast<long long>(area.side.n));
        js.num("A_side_raw", area.side.raw);
        js.num("A_side_target", area.side.target);
        js.num("A_side_rel_err_raw", area.side.rel_err_raw);
        js.num("A_side_after", area.side.after);
        js.num("A_side_rel_err", area.side.rel_err);
        js.num("A_side_scale", area.side.scale);
        js.num("A_side_ai_min", area.side.a_i_min);
        js.num("A_side_ai_max", area.side.a_i_max);
        js.integer("n_bottom", static_cast<long long>(area.bottom.n));
        js.num("A_bottom_raw", area.bottom.raw);
        js.num("A_bottom_target", area.bottom.target);
        js.num("A_bottom_rel_err_raw", area.bottom.rel_err_raw);
        js.num("A_bottom_after", area.bottom.after);
        js.num("A_bottom_rel_err", area.bottom.rel_err);
        js.num("A_bottom_scale", area.bottom.scale);
        js.num("A_bottom_ai_min", area.bottom.a_i_min);
        js.num("A_bottom_ai_max", area.bottom.a_i_max);
        // ---- thermal ledger ----
        js.num("thermal_start_time_s", 0.0);
        js.num("thermal_duration_s", ledger.t_thermal);
        js.integer("thermal_steps", static_cast<long long>(ledger.n_steps));
        js.num("thermal_dt_min_s", ledger.n_steps ? ledger.dt_min : 0.0);
        js.num("thermal_dt_max_s", ledger.dt_max);
        js.num("thermal_dt_stability_limit_s", dt_thermal_limit);
        js.num("E_source_J", ledger.e_source);
        js.num("E_boundary_J", ledger.e_boundary);
        js.num("DeltaH_J", d_h);
        js.num("energy_residual_J", r_e);
        js.num("energy_error_input", err_input);
        js.num("energy_error_exchange", err_exchange);
        if (net_ok)
            js.num("energy_error_net", err_net);
        else
            js.gate("energy_error_net", Gate::Null, "net input |E_source-E_boundary| < 1% of the larger term");
        js.num("energy_error_unified", err_unified);
        js.flag("energy_error_unified_le_1pct", err_unified <= 0.01);
        js.num("energy_floor_J", e_floor);
        js.num("dT_mean_mass_K", dT_end.mean_m);
        js.num("dT_min_K", dT_end.min);
        js.num("dT_max_K", dT_end.max);
        js.num("P_loss_initial_W", ledger.p_loss_first);
        for (const std::string &it : tj.items)
            js.items.push_back(it);
        // ---- narrow baseline naming (Doc131 §2.2) ----
        js.flag("thermal_smoke_baseline", is_baseline);
        js.flag("literature_baseline", false);
        js.flag("calibrated_match", false);
        js.flag("thaw_T", cli.thaw_T);
        if (is_baseline)
        {
            js.flag("top_loss", cli.top_loss);
            js.num("E_in_J", ledger.e_source);
            js.num("E_rad_loss_J", ledger.e_boundary);
            js.num("E_thermal_end_J", d_h);
            js.num("energy_balance_rel", err_input);
        }
        // ---- Doc131 gate flags (per run; OPEN = null + reason, never default true) ----
        js.gate("source_power_verified", gate_source_power, gate_source_power_reason);
        const bool area_ok = std::abs(area.top.rel_err) <= 0.02 && std::abs(area.side.rel_err) <= 0.02 &&
                             std::abs(area.bottom.rel_err) <= 0.02;
        js.gate("boundary_area_verified", area_ok ? Gate::Yes : Gate::No,
                cli.area_normalize ? "face area after normalisation deviates > 2% from the geometric target"
                                   : "face area not normalised (--area-normalize off) and raw error > 2%");
        js.gate("conduction_checked", gate_conduction, gate_conduction_reason);
        js.gate("thermal_boundary_checked", gate_thermal_bc, gate_thermal_bc_reason);
        js.gate("weak_buoyancy_response_checked", Gate::Null, "OPEN: Doc131 D (R1 extended, R2 +/-/zero beta) not in this run");
        if (sizeof(Real) == 4)
            js.gate("precision_sensitivity_checked", Gate::Null, "build_is_float_only_no_double_binary");
        else
            js.gate("precision_sensitivity_checked", Gate::Null,
                    "OPEN: double build present but no paired float run recorded by this binary");
        js.gate("slip_boundary_checked", Gate::Null,
                "no_independent_plate_slip_mms_yet (top-band |uz|, side/bottom band and escape counts are reported "
                "as d4_* fields, but that is not an independent slip/shear accuracy check)");
        {
            const bool have_vtp = cli.state_record_interval > TinyReal && n_vtp_written > 0;
            const bool have_em_shape = em_rz && q_active;
            if (is_baseline && have_vtp && have_em_shape && cli.area_normalize)
                js.gate("full_thermal_flow_baseline", Gate::No,
                        "PARTIAL: VTP + prior-EM rz CSV shape + discrete power + area-normalize present; still missing "
                        "live EM solve on this lattice, longer thermal/flow window, and CEP digitised targets");
            else
                js.gate("full_thermal_flow_baseline", Gate::No,
                        "OPEN: need baseline with --q-from-rz-csv + --state-record-interval>0 + --area-normalize; "
                        "live EM solve and CEP targets still separate");
        }
        js.gate("literature_targets_available", Gate::No, "OPEN: CEP target data not digitised");
        js.str("note", note);
        js.flag("passed", passed);
        std::ofstream out(cli.summary_json);
        out << js.dump();
        std::cout << "[doc126] summary_json=" << cli.summary_json << std::endl;
    }

    return passed ? 0 : 2;
}
