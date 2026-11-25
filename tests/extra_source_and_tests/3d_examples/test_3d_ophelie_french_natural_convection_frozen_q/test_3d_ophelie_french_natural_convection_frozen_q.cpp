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
#include "electromagnetic_ophelie_client_case_params.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_doc245_material_dataset.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_fixed_space_probe.h"
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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
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

inline std::string resolveDefaultFrenchReloadFolder(bool prefer_client_nc = false)
{
    const StdVec<std::string> client_candidates = {
        "../test_3d_ophelie_client_nc_glass_relax/bin/reload",
        "../../test_3d_ophelie_client_nc_glass_relax/bin/reload",
        "./reload",
        "../reload",
    };
    const StdVec<std::string> french_candidates = {
        "./reload",
        "../reload",
        "../../../../../reload",
        "../test_3d_ophelie_client_nc_glass_relax/bin/reload",
        "../../test_3d_ophelie_client_nc_glass_relax/bin/reload",
        "../test_3d_ophelie_french_natural_glass_relax/bin/reload",
        "../../test_3d_ophelie_french_natural_glass_relax/bin/reload",
    };
    const StdVec<std::string> &candidates = prefer_client_nc ? client_candidates : french_candidates;
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return candidates.front();
}

/** Doc224 D: regular r–z Kelvin grid from CSV (bilinear sample; O(n_particles) reapply). */
struct Doc224RzTemperatureGrid
{
    std::vector<double> r_axis;
    std::vector<double> z_axis;
    std::vector<double> T; // row-major: ir + nr*iz
    double T_min = 0;
    double T_max = 0;
    size_t n_nodes = 0;
    bool ok = false;

    size_t nr() const { return r_axis.size(); }
    size_t nz() const { return z_axis.size(); }
    double at(size_t ir, size_t iz) const { return T[ir + nr() * iz]; }

    double sample(double r, double z) const
    {
        if (!ok || nr() < 2 || nz() < 2)
            return T.empty() ? 0.0 : T.front();
        r = std::max(r_axis.front(), std::min(r_axis.back(), r));
        z = std::max(z_axis.front(), std::min(z_axis.back(), z));
        auto upper = [](const std::vector<double> &ax, double v) {
            return static_cast<size_t>(
                std::max<std::ptrdiff_t>(1, std::lower_bound(ax.begin(), ax.end(), v) - ax.begin()));
        };
        const size_t ir1 = std::min(upper(r_axis, r), nr() - 1);
        const size_t iz1 = std::min(upper(z_axis, z), nz() - 1);
        const size_t ir0 = ir1 - 1;
        const size_t iz0 = iz1 - 1;
        const double r0 = r_axis[ir0], r1 = r_axis[ir1];
        const double z0 = z_axis[iz0], z1 = z_axis[iz1];
        const double tr = (r1 > r0) ? (r - r0) / (r1 - r0) : 0.0;
        const double tz = (z1 > z0) ? (z - z0) / (z1 - z0) : 0.0;
        const double t00 = at(ir0, iz0), t10 = at(ir1, iz0), t01 = at(ir0, iz1), t11 = at(ir1, iz1);
        return (1.0 - tr) * (1.0 - tz) * t00 + tr * (1.0 - tz) * t10 + (1.0 - tr) * tz * t01 + tr * tz * t11;
    }
};

inline bool hostLoadDoc224RzTemperatureGrid(const std::string &csv_path, Doc224RzTemperatureGrid &grid,
                                            std::string &err)
{
    grid = Doc224RzTemperatureGrid{};
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
    std::vector<double> rs, zs, Ts;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(ss, cell, ','))
            cells.push_back(strip(cell));
        if (static_cast<int>(cells.size()) <= std::max(i_r, std::max(i_z, i_T)))
            continue;
        const double r = std::atof(cells[static_cast<size_t>(i_r)].c_str());
        const double z = std::atof(cells[static_cast<size_t>(i_z)].c_str());
        const double Tv = std::atof(cells[static_cast<size_t>(i_T)].c_str());
        if (!(std::isfinite(r) && std::isfinite(z) && std::isfinite(Tv)) || r < 0.0 || !(Tv > 0.0))
        {
            err = "non-finite/non-positive T or bad r in init-T-rz-csv";
            return false;
        }
        rs.push_back(r);
        zs.push_back(z);
        Ts.push_back(Tv);
    }
    if (rs.empty())
    {
        err = "no usable rows in init-T-rz-csv";
        return false;
    }
    auto uniq_sorted = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(),
                            [](double a, double b) { return std::fabs(a - b) < 1e-12; }),
                v.end());
        return v;
    };
    grid.r_axis = uniq_sorted(rs);
    grid.z_axis = uniq_sorted(zs);
    if (grid.nr() < 2 || grid.nz() < 2 || grid.nr() * grid.nz() != rs.size())
    {
        err = "init-T-rz-csv is not a complete regular r-z grid";
        return false;
    }
    grid.T.assign(grid.nr() * grid.nz(), std::numeric_limits<double>::quiet_NaN());
    grid.T_min = Ts.front();
    grid.T_max = Ts.front();
    for (size_t k = 0; k < rs.size(); ++k)
    {
        const auto ir = static_cast<size_t>(
            std::lower_bound(grid.r_axis.begin(), grid.r_axis.end(), rs[k]) - grid.r_axis.begin());
        const auto iz = static_cast<size_t>(
            std::lower_bound(grid.z_axis.begin(), grid.z_axis.end(), zs[k]) - grid.z_axis.begin());
        if (ir >= grid.nr() || iz >= grid.nz() || std::fabs(grid.r_axis[ir] - rs[k]) > 1e-12 ||
            std::fabs(grid.z_axis[iz] - zs[k]) > 1e-12)
        {
            err = "init-T-rz-csv grid axis mismatch";
            return false;
        }
        grid.T[ir + grid.nr() * iz] = Ts[k];
        grid.T_min = std::min(grid.T_min, Ts[k]);
        grid.T_max = std::max(grid.T_max, Ts[k]);
    }
    for (double v : grid.T)
    {
        if (!std::isfinite(v))
        {
            err = "init-T-rz-csv missing grid node";
            return false;
        }
    }
    grid.n_nodes = rs.size();
    grid.ok = true;
    return true;
}

/**
 * Doc225 B: sample CSV absolute T; thermal_deltaT = T_abs - T_thermal_ref.
 * `legacy_demean` restores Doc224 ABSOLUTE_T_SHIFTED_LEGACY (shifts Temperature); forbidden on
 * acceptance mainline.
 */
inline void hostApplyDoc224RzTemperatureGrid(BaseParticles &particles, const Vecd &center, Real t_thermal_ref,
                                             const Doc224RzTemperatureGrid &grid, bool legacy_demean, bool quiet)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double Tmin = 1.0e300, Tmax = -1.0e300, Tsum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        const double Ti = grid.sample(r, z);
        temp[i] = static_cast<Real>(Ti);
        delta_t[i] = static_cast<Real>(Ti) - t_thermal_ref;
        Tmin = std::min(Tmin, Ti);
        Tmax = std::max(Tmax, Ti);
        Tsum += Ti;
    }
    double demean_shift = 0.0;
    if (legacy_demean && n > 0)
    {
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i)
            sum += static_cast<double>(delta_t[i]);
        demean_shift = sum / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i)
        {
            delta_t[i] = static_cast<Real>(static_cast<double>(delta_t[i]) - demean_shift);
            temp[i] = t_thermal_ref + delta_t[i];
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    if (!quiet)
    {
        std::cout << "[ophelie][doc225][T-rz] n_nodes=" << grid.n_nodes << " nr=" << grid.nr()
                  << " nz=" << grid.nz() << " T_csv=[" << grid.T_min << "," << grid.T_max
                  << "] T_particle=[" << Tmin << "," << Tmax << "] T_particle_arith_mean="
                  << (n > 0 ? Tsum / static_cast<double>(n) : 0.0)
                  << " T_thermal_ref=" << t_thermal_ref << " legacy_demean=" << (legacy_demean ? 1 : 0)
                  << " demean_shift_K=" << demean_shift
                  << (legacy_demean ? " ABSOLUTE_T_SHIFTED_LEGACY" : " ABSOLUTE_T_PRESERVED") << std::endl;
    }
}

inline Real hostMassWeightedMeanTemperature(BaseParticles &particles,
                                            const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "Mass");
    const Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    double wT = 0.0, w = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        wT += static_cast<double>(mass[i]) * static_cast<double>(temp[i]);
        w += static_cast<double>(mass[i]);
    }
    return w > 0.0 ? static_cast<Real>(wT / w) : Real(0);
}

/** Doc225 B4: bilinear sample vs analytic T=1000+100*r/R+50*z/H on particle positions. */
inline bool hostDoc225AnalyticRzSampleSelfTest(BaseParticles &particles, const Vecd &center, Real R, Real H,
                                               Real &max_abs_err_K)
{
    Doc224RzTemperatureGrid grid;
    grid.r_axis.resize(3);
    grid.z_axis.resize(3);
    grid.T.resize(9);
    grid.n_nodes = 9;
    grid.r_axis = {0.0, 0.5 * static_cast<double>(R), static_cast<double>(R)};
    grid.z_axis = {0.0, 0.5 * static_cast<double>(H), static_cast<double>(H)};
    auto setT = [&](size_t ir, size_t iz, double Ti) { grid.T[iz * 3 + ir] = Ti; };
    for (size_t ir = 0; ir < 3; ++ir)
        for (size_t iz = 0; iz < 3; ++iz)
        {
            const double r = grid.r_axis[ir], z = grid.z_axis[iz];
            setT(ir, iz, 1000.0 + 100.0 * r / static_cast<double>(R) + 50.0 * z / static_cast<double>(H));
        }
    grid.T_min = 1000.0;
    grid.T_max = 1000.0 + 100.0 + 50.0;
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    max_abs_err_K = 0.0;
    const size_t n_check = std::min(n, size_t(512));
    for (size_t i = 0; i < n_check; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        const double exact = 1000.0 + 100.0 * r / static_cast<double>(R) + 50.0 * z / static_cast<double>(H);
        const double samp = grid.sample(r, z);
        max_abs_err_K = std::max(max_abs_err_K, static_cast<Real>(std::abs(samp - exact)));
    }
    return max_abs_err_K <= Real(1.0e-3);
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
    Real hydro_settle_time = 0.15; // gravity-only settle before Q / Boussinesq (max physical time)
    /**
     * Doc245: physics-based settle exit (GUI hangs when acoustic dt collapses under pure time gate).
     * Exit early when U_max < u_quiet for quiet_steps consecutive advection steps after min_time.
     * u_quiet<=0 disables quiet exit (pure --hydro-settle-time). After exit, velocity is still zeroed.
     */
    Real hydro_settle_u_quiet = 5.0e-3;
    size_t hydro_settle_quiet_steps = 20;
    Real hydro_settle_min_time = 0.01;
    /** Hard cap on advection steps during settle (0 = disabled). Guards infinite tiny-dt loops. */
    size_t hydro_settle_max_adv_steps = 100000;
    /**
     * If acoustic_dt < floor for collapse_steps consecutive acoustic calls, end settle with WARNING.
     * floor<=0 → auto 1e-3 * nominal dt_ac after hydro setup.
     */
    Real hydro_settle_dt_floor = 0.0;
    size_t hydro_settle_dt_collapse_steps = 50;
    Real escaped_frac_max = 0.01;
    Real state_record_interval = 0.0; // 0 -> no VTP
    /** Start time-averaged meridional (r,z) sampling after this physical time (post-settle clock). */
    Real meridional_average_start = 0.5;
    /** P1-03: if >0, prescribe hot-bottom linear T: T=T0+dT*(1-ξ), ξ=(z-zb)/H (independent of EM). */
    Real prescribe_t_hot_bottom_dT = 0.0;
    /** Doc77 §3: if >0, prescribe radial hot-core T=T0+dT*(1-r/R); overrides hot-bottom if both set. */
    Real prescribe_t_radial_dT = 0.0;
    /**
     * Doc224 D: axisymmetric T(r,z) CSV (r_m,z_m,T_K). Exclusive with F2 / radial / hot-bottom.
     * With --freeze-prescribed-t, re-sampled each advection (Eulerian fixed spatial T).
     */
    std::string init_T_rz_csv;
    /**
     * Doc224 legacy: demean shifted absolute Temperature (ABSOLUTE_T_SHIFTED_LEGACY).
     * Doc225: bare `--init-T-rz-demean` errors; use `--init-T-rz-demean-legacy` to reproduce only.
     */
    bool init_T_rz_demean_legacy = false;
    /**
     * Doc225 B3: after Fig8 CSV load, set T_force_ref = mass-weighted mean(T) once (default ON for CSV).
     * F2 / radial / hot-bottom keep literature t0. Disable with `--no-t-force-ref-from-mass-mean`.
     */
    bool t_force_ref_from_mass_mean = true;
    bool t_force_ref_from_mass_mean_user_set = false;
    /** Explicit override; NaN => unset. */
    Real t_force_ref_override = std::numeric_limits<Real>::quiet_NaN();
    /** Doc225 A: refuse warm Spatial-Q runs that miss thermal diffusion / use legacy demean / etc. */
    bool doc225_warm_spatial_q_preflight = false;
    bool freeze_prescribed_temperature = false;
    bool prescribe_radial_z_taper = false;
    bool require_validation = false;
    bool thermal_diffusion_explicit = false;
    bool boundary_area_normalize = true;
    bool enable_self_induction = true;
    bool enable_target_power = true;
    /**
     * Doc231/232: development prescribed absorbed power (uniform volumetric Q).
     * Mode is armed only when --prescribed-heat-power= appears (Doc232 H), not by P>0 alone.
     * Skips EM; not generator 30 kW; not confirmed experiment Q.
     */
    bool prescribed_heat_mode = false;
    Real prescribed_heat_power_w = Real(0);
    /** Explicit CLI sightings for order-independent conflict checks (Doc232 H). */
    bool zero_q_flag_seen = false;
    bool target_power_flag_seen = false;
    bool freeze_prescribed_t_flag_seen = false;
    /** Doc219 E: host D2 for frozen Q snapshot (opt-in; matches C/D EM identity). */
    bool doc217_d2_local_poly = false;
    /**
     * Doc219 E: freeze Q as Eulerian spatial field; re-sample to particles each thermal update.
     * Default ON so particle motion does not advect Q as a material label.
     */
    bool spatial_frozen_q = true;
    /**
     * Doc95: quieter WCSPH defaults from Doc91/92 (free-surface: dens-sum normalize +
     * analytic p_h + LinearCorrectionCK). Off by default for regression; enable for lit-β SNR.
     */
    bool dens_sum_normalize = false;
    bool dens_sum_normalize_user_set = false;
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
    /**
     * Doc219 D: validated cavity hydro identity (B + Riemann correction + TVF 0.2 +
     * viscous advection dt + thermal after position/neighborhood update).
     * Default ON; --no-doc219-hydro-profile restores legacy French path.
     */
    bool doc219_hydro_profile = true;
    bool kernel_correction = false; // resolved after CLI from doc219 profile
    bool kernel_correction_user_set = false;
    bool enable_tvf = false;
    bool enable_tvf_user_set = false;
    /**
     * Doc105 P2d opt-in: AcousticDensMode for wall dens dissipation.
     * 0 = legacy (default): dens jump on Adami momentum pressure (−C_ext).
     * 1 = dens on p' with ReferenceBodyAcceleration; French keeps total-p momentum
     *     (PerturbationMomentum=0) + GravityForceCK, so a_ref must be −ĝ (g in ForcePrior).
     * Do NOT copy closed-box P2b (PM=1, GravityForceCK off, a_ref=0).
     */
    int acoustic_dens_mode = 0;
    bool acoustic_dens_mode_user_set = false;
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
     * Doc220 S0: after settle (velocity zeroed), prescribe radial hot-core if requested,
     * dump by OriginalID, run one ParticleSortCK with no physics, dump again, compare, exit.
     * Opt-in audit; never a production default.
     */
    bool doc220_s0_sort_audit = false;
    /**
     * Doc220 D/H0: same as --doc105-p2-hydro-stage-probe, plus per-particle CSV dumps
     * (OriginalID) at dens/gravity/acoustic stages for same-state operator diagnosis.
     */
    bool doc220_h0_hydro_probe = false;
    /**
     * Doc220 E2: French free-surface perturbation-pressure hydrostatic package (not closed-box
     * P2b copy). Enables BackgroundPressure p_h, PerturbationMomentum=1 (momentum on p'),
     * AcousticDensMode=1 with a_ref=0, and omits GravityForceCK (g lives in ∇p_h).
     * Free surface: depth→0 ⇒ p_h→0 ⇒ p_total≈p' ; keep Density~ρ0 so p'≈0 at rest.
     * Opt-in; never a production default until H1 gate passes.
     */
    bool doc220_pprime_hydro = false;
    /** When true, skip GravityForceCK (used by doc220_pprime_hydro). */
    bool omit_gravity_force = false;
    /**
     * Doc221 A: restore the Doc220 omit-gravity bug that wiped shared ForcePrior each
     * advection step (evidence-only). Default false = Doc221 B fixed omit (noop).
     */
    bool doc221_legacy_omit_clears_forceprior = false;
    /**
     * Doc221: dump ForcePrior ledger stages for the first N keep-loop advection steps
     * (0 = off). Writes output/doc221_force_ledger_stages.csv (+ summary log lines).
     */
    int doc221_force_ledger_probe_adv = 0;
    /** Doc221 C/D: after first keep advection force assemble, run lifecycle asserts then exit. */
    bool doc221_force_lifecycle_selftest = false;
    /**
     * Doc221 K diagnostic: before each position update, subtract mass-weighted mean Uz from all
     * fluid particles (enforce zero net vertical volume flux). Controlled closed-topology probe —
     * NOT a literature open free-surface claim; never a production default.
     */
    bool diag_zero_net_uz = false;
    /**
     * Doc221 L diagnostic: kinematic hard lid at z_free — kill upward Uz in a surface band and
     * clamp z <= z_free after position update (impermeable free-slip top). NOT literature open FS.
     * `--diag-flat-lid` enables this and also `--diag-zero-net-uz`.
     */
    bool diag_hard_lid = false;
    /**
     * Doc221 N diagnostic: soft lid on K skeleton — enables zero-net-Uz and kills upward Uz in a
     * thin surface band, but does NOT clamp z (no hard lid / no M settle reset). Goal: keep
     * axis-up/rim-down polarity while discouraging tip venting. NOT literature open FS.
     */
    bool diag_soft_lid = false;
    /**
     * Doc221 O diagnostic: physical closed crucible wall (side+bottom+flat lid via
     * OphelieFrenchClosedCrucibleWallShape). No kinematic kill-up / z-clamp. Optional pairing with
     * `--diag-zero-net-uz` for ablation. NOT literature open free-surface.
     */
    bool diag_wall_flat_lid = false;
    /**
     * Doc221 P diagnostic: quiet ClosedCrucible lid contact (Doc201 B2 pattern).
     * Implies `--diag-wall-flat-lid`. Acoustic stays on closed wall with AcousticWallFreeSlip on
     * the flat lid (mom∥n, omit dens DissipativeUJump, free-slip dens ghost); viscous no-slip
     * uses open-cup wall (rim=0) so the lid has no viscous traction. Goal: U_rms→O(10⁻³).
     * NOT literature open free-surface; never a production default.
     * P1: arm free-slip from t=0 (shown noisy). Prefer P2 post-settle arming below.
     */
    bool diag_lid_free_slip = false;
    /**
     * Doc221 P2: same geometry as `--diag-lid-free-slip`, but settle uses O-style no-slip lid
     * (no AcousticWallFreeSlip, viscous on ClosedCrucible). After settle + hostZeroVelocity,
     * tag free-slip and switch viscous to OpenCup. NOT literature open FS.
     */
    bool diag_lid_free_slip_post_settle = false;
    /**
     * Doc221 O3 ablation: ClosedCrucible acoustic (no AcousticWallFreeSlip) + viscous OpenCup
     * (rim=0, lid omitted). Isolates viscous-omit-lid vs P acoustic free-slip noise.
     * Implies `--diag-wall-flat-lid`. NOT literature open FS.
     */
    bool diag_wall_lid_viscous_omit = false;
    /**
     * Doc221 O4: quiet ClosedCrucible dens/Riemann WITHOUT AcousticWallFreeSlip.
     * Implies O3 skeleton + acoustic_dens_mode=3 (omit DissipativeUJump) + lid
     * AcousticWallContSilent + AcousticWallOmitPJump. NOT literature open FS.
     */
    bool diag_lid_dens_silence = false;
    /**
     * Doc221 O5: O3 + lid ContSilent/OmitPJump ONLY, keep dens_mode=1.
     * Split of O4 package — no dens_mode=3, no AcousticWallFreeSlip.
     * NOT literature open FS.
     */
    bool diag_lid_cont_silent = false;
    /**
     * Doc221 O6a: O3 + lid ContSilent ONLY (no OmitPJump), dens_mode=1.
     * Split of O5. NOT literature open FS.
     */
    bool diag_lid_cont_silent_only = false;
    /**
     * Doc221 O6b: O3 + lid OmitPJump ONLY (no ContSilent), dens_mode=1.
     * Split of O5. NOT literature open FS.
     */
    bool diag_lid_omit_pjump_only = false;
    /**
     * Doc221 N2: soft tip seal — K zero-net-Uz + kill upward Uz only for r < tip_r
     * in the surface band (no full-lid kill-up, no z clamp). Tests whether tip-local
     * venting seal recovers TOPO. NOT literature open FS.
     */
    bool diag_soft_lid_tip = false;
    /** Tip radius for --diag-soft-lid-tip (m). Default ~J r_eff. */
    Real diag_soft_lid_tip_r = 0.10;
    /**
     * Doc221 N3: scale applied to upward Uz inside tip seal (0=kill=N2, 1=no-op, 0.5=dampen).
     * Only used with diag_soft_lid_tip. NOT literature open FS.
     */
    Real diag_soft_lid_tip_uz_scale = 0.0;
    /**
     * Doc221 N4: tip soft physical lid — K zero-net-Uz + ForcePrior spring/dashpot in tip
     * surface band (NO host Uz kill/dampen, NO z clamp). NOT literature open FS.
     */
    bool diag_tip_soft_force = false;
    /** Tip radius for --diag-tip-soft-force (m). Default = min kill PASS ~0.05. */
    Real diag_tip_soft_force_r = 0.05;
    /**
     * Spring rate a_z -= k * max(z-z_free,0) [1/s^2].
     * Default 0: N4 ablation showed spring harms polarity; N4c dashpot-only PASS.
     */
    Real diag_tip_soft_force_k = 0.0;
    /** Upward dashpot a_z -= c * max(Uz,0) [1/s]. Default 100 = N4c PASS. */
    Real diag_tip_soft_force_c = 100.0;
    /**
     * Doc221 N5: same tip ForcePrior dashpot as N4c, but only on free-surface incomplete
     * particles (DensitySummation < rho0) inside tip cylinder — binds seal to FS path
     * rather than pure geometric surface band. k forced 0. NOT literature open FS.
     */
    bool diag_tip_fs_dashpot = false;
    /**
     * Doc222 B0: restore Doc221 tip ForcePrior direct-accumulate bug (no previous revoke).
     * Opt-in repro only; never a production/default path.
     */
    bool diag_tip_soft_force_legacy_accumulate = false;
    /**
     * Doc222 B: after particles+ForcePrior exist, run tip ForcePrior lifecycle self-test and exit.
     */
    bool doc222_tip_force_selftest = false;
    /**
     * Doc222 D2: after hydro settle, apply volume-preserving axisymmetric free-surface
     * deformation η(r)=eps*[1-2(r/R)^2] via column stretch, then re-init density/pressure
     * per active path. eps=0 disables. NOT a tip/lid intervention.
     */
    Real doc222_fs_deform_eps = 0.0;
    /**
     * Doc222 E1: BackgroundPressure p_h = ρ0 g (z_free−z) with signed depth (default true).
     * Legacy Doc220/221 used max(0,z_free−z), which zeros ∇p_h above the reference plane and
     * (with omit GravityForceCK) removes free-surface restoring — D2 FAIL. Opt-out:
     * --doc222-ph-clamp-legacy.
     */
    bool doc222_ph_signed = true;
    /**
     * Doc222 C3: after settle(+deform), compare clamp vs signed BackgroundPressure on the
     * same geometry via one acoustic-1st force probe each; write JSON and exit.
     * Implies p' path; if fs_deform_eps==0 uses 0.02.
     */
    bool doc222_c3_ph_compare = false;
    /**
     * Doc222 E2a: after free-surface deform, keep FreeSurface Density (do NOT freeze to ρ0)
     * so p' = Pressure−p_h can become nonzero. Default true. Opt-out:
     * --doc222-deform-rho-freeze (legacy E1/D2 path where p_mom≡0).
     */
    bool doc222_deform_skip_rho_freeze = true;
    /**
     * Doc222 E2b: BackgroundPressure from local free-surface height η(r) (radial bin z_max),
     * p_h=ρ0 g (η−z). Fixed z_free reference cannot restore η(r) (∇p_h has no radial part).
     * Opt-in via --doc222-ph-local-eta; --doc222-ph-fixed-zfree restores E1/E2a path.
     */
    bool doc222_ph_local_eta = false;
    int doc222_ph_local_eta_bins = 12;
    /**
     * Doc223 B: after hydro settle, dump ForcePrior ledger (post-settle → legacy
     * ForcePrior-only zero + gravity.exec → full channel reset + gravity.exec) then exit.
     * Requires GravityForceCK ON (do not combine with --doc220-pprime-hydro).
     */
    bool doc223_force_reset_probe = false;
    /**
     * Doc223 D: paper flat slip impermeable top (single entry). ClosedCrucible + P2
     * post-settle AcousticWallFreeSlip + OpenCup viscous + p' (omit GravityForceCK) +
     * Internal density (not FreeSurface). Illegal with tip/zero-net/z-clamp/ContSilent/deform.
     */
    bool doc223_paper_flat_slip_top = false;
    /**
     * Doc223 E1 ablation: keep paper-top closed domain (Internal + dens_norm + p') but
     * leave lid acoustic/viscous as ClosedCrucible no-slip (no AcousticWallFreeSlip).
     * Requires --doc223-paper-flat-slip-top. Not the paper slip model.
     */
    bool doc223_ablate_no_free_slip = false;
    /**
     * Doc224 C3: on AcousticWallFreeSlip lid, keep DissipativeUJump dens (impermeability).
     * Mom force still projected onto n; OpenCup viscous unchanged. Default ON for paper-top.
     * Disable with --doc224-no-lid-keep-dens-jump (restores Doc201 free-slip dens omit).
     */
    bool doc224_lid_keep_dens_jump = false;
    bool doc224_lid_keep_dens_jump_user_set = false;
    /**
     * Doc224 C3c ablation: arm AcousticWallFreeSlip (+keep dens) but keep wall viscous on
     * ClosedCrucible (do not switch to OpenCup). Isolates OpenCup lid viscous omit.
     */
    bool doc224_ablate_opencup = false;
    /**
     * Doc224 C3d: OpenCup viscous (no lid shear) but do NOT tag AcousticWallFreeSlip
     * (acoustic stays ClosedCrucible / impermeable). Doc221 O3-style.
     * Default ON for --doc223-paper-flat-slip-top after C3a–C3d evidence.
     * Opt out with --doc224-acoustic-freeslip (legacy P2; penetrates ~0.26dp by 10s).
     */
    bool doc224_ablate_acoustic_freeslip = false;
    bool doc224_ablate_acoustic_freeslip_user_set = false;
    /**
     * Doc224 U_max blowup gate override (m/s). <0 → auto:
     *   zero-q / prescribed-heat → 0.05·c0 (Doc230 soft);
     *   heated French NC → absolute 0.05 (legacy red-line).
     * Client long demos pass an explicit value (e.g. 50) so early acoustic
     * spikes do not abort before useful flow time; non-finite U/T still aborts.
     */
    Real doc224_u_blowup_gate = Real(-1);
    /**
     * Doc223 E2: after paper-top settle+arm, impose uniform tangential U=(U0,0,0),
     * measure wall-only ViscousForce in lid-core vs side/bot; OpenCup (slip) vs
     * ClosedCrucible (no-slip lid) contrast. Snapshot/restore; early return.
     * Requires --doc223-paper-flat-slip-top (P2 OpenCup contact available).
     */
    bool doc223_e2_lid_slip_probe = false;
    /**
     * Doc226 B: after paper-top settle+arm, fixed-state top-operator probe.
     * U=(0,0,0) vs U=(0.01,0,0); split wall 1st-half pressure, 2nd-half dissipative,
     * continuity DensityChangeRate, OpenCup vs Closed viscous; lid_core only.
     * Snapshot/restore; early return. Requires --doc223-paper-flat-slip-top.
     */
    bool doc226_top_operator_probe = false;
    /**
     * Doc226 D/E: bind ParticleViscosity → OphelieDynamicViscosity (pairwise/one-side).
     * Default OFF keeps Viscosity mean path. Does not alone change cp/k/enthalpy.
     */
    bool doc226_local_particle_viscosity = false;
    /**
     * Doc226: refresh μ(T) onto particle field while keeping cp/k at T0 diagnostic values.
     * Requires --doc226-local-particle-viscosity. Independent of constant_thermal_rheology energy path.
     */
    bool doc226_mu_T_only = false;
    /** Optional μ_cap (Pa·s); 0 = no clamp. Diagnostic only when set. */
    Real doc226_mu_cap = Real(0);
    /**
     * Doc226 E1: after settle, fixed-state local-μ probe (uniform μ=4 vs half-swap 4/40),
     * snapshot/restore; early return. Requires --doc226-local-particle-viscosity.
     */
    bool doc226_e1_local_mu_probe = false;
    /**
     * Doc227 B: matched outer advection budget — truncate each acoustic substep so
     * sum(dt_actual) ≤ min(inertial, viscous) and viscous ForcePrior age stays inside limit.
     * Isolates French keep-loop only; does not change global TriggerByInterval semantics.
     */
    bool doc227_actual_timestep = false;
    /** Optional CSV path for actual_step_audit (Doc227 B1). Empty → output/doc227_actual_step_audit.csv */
    std::string doc227_step_audit_csv;
    /** Host-only schedule truncation self-test; early return. */
    bool doc227_schedule_selftest = false;
    /**
     * Doc227 C1: after settle, same-state Viscosity vs uniform ParticleViscosity (μ=4)
     * force comparison; snapshot/restore; early return. Needs --doc226-local-particle-viscosity.
     */
    bool doc227_c_uniform_mu_probe = false;
    /**
     * Doc223 E3: after paper-top settle+arm, impose smooth near-lid normal velocity
     * (mass-mean removed, Ma≪1), then integrate; gate on no lid penetration / no
     * runaway KE / no COM lift. Requires --doc223-paper-flat-slip-top.
     */
    bool doc223_e3_normal_perturb = false;
    /**
     * Doc223 F2: qualitative annular hot ring T (plan §11 fixed formula; NOT Fig8 digitization).
     * Mutual exclusive with --prescribe-t-radial / hot-bottom.
     */
    bool doc223_f2_annular_hot = false;
    /** PerturbationMomentum flag written to particles (0=total p, 1=p'). */
    int perturbation_momentum = 0;
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
    /**
     * Doc232 B: fixed-space probe timeseries CSV (production sampler). Empty =>
     * default output/probes_timeseries.csv when --client-nc-geometry.
     */
    std::string probes_timeseries_csv;
    bool probes_timeseries_user_set = false;
    /** Doc232 B2: run uniform+linear manufacture gate on the production sampler before keep. */
    bool doc232_probe_path_verify = false;
    /** If set with path-verify, exit after gate (no warm keep). */
    bool doc232_probe_path_verify_exit = false;
    /**
     * Doc232 B2: dump host Position/Temperature(/Mass) CSV at the last probe sample for offline
     * Wendland recompute. Empty => no dump.
     */
    std::string doc232_dump_particles_csv;
    /**
     * Doc232 C: adiabatic production heat gate — prescribed Q via ApplyFrozenQToTemperatureCK,
     * no Robin/radiation/diffusion, no fluid keep (thermal-only time advance).
     */
    bool doc232_adiabatic_heat = false;
    /**
     * Desktop / client thermal wall BC (orthogonal to Doc232 C gate).
     *   french-natural | adiabatic | isothermal-wall
     * Optional --h-side=/--h-bottom=/--h-free=/--emissivity=/--t-cool= for french-natural.
     */
    std::string thermal_bc_mode = "french-natural";
    bool thermal_h_side_set = false;
    bool thermal_h_bottom_set = false;
    bool thermal_h_free_set = false;
    bool thermal_emissivity_set = false;
    bool thermal_t_cool_set = false;
    Real thermal_h_side = Real(300);
    Real thermal_h_bottom = Real(35);
    Real thermal_h_free = Real(20);
    Real thermal_emissivity = Real(0.8);
    Real thermal_t_cool = Real(300);
    /** Doc232 D: face-by-face boundary audit on production BC operator; exit after. */
    bool doc232_boundary_face_audit = false;
    bool doc232_boundary_face_audit_exit = false;
    std::string doc232_boundary_audit_dir = "output/doc232_boundary_audit";
    /**
     * Doc233 D: face-T reconstruction mode for production BC.
     * "global_dt_dz" = Doc232 mass-weighted vertical fit (experimental).
     * "band" = shell-band particle T without face extrapolation (legacy).
     */
    std::string boundary_face_t_mode = "global_dt_dz";
    /** Doc233 D: also run radial + parabolic-z manufacturing fields (scope audit). */
    bool doc233_boundary_scope_audit = false;
    /**
     * Doc233 C3: after IC (prefer non-uniform T), verify buoyancy force chain and exit.
     */
    bool doc233_buoyancy_force_chain_exit = false;
    /**
     * Doc234 B: same frozen snapshot force-chain as Doc233, plus per-particle CSV,
     * host-double reference, Float cancellation metrics, and one acoustic read probe.
     */
    bool doc234_force_precision_exit = false;
    std::string doc234_force_precision_csv;
    /**
     * Doc235 B: dump OFF-background event snapshots (settle cleanup + early keep times).
     * Directory receives event_*.csv and event_stats.jsonl.
     */
    std::string doc235_off_event_dump_dir;
    /** Extra keep dump times (seconds) when off-event dump is set. Default 0.01,0.05,0.1 */
    std::string doc235_off_event_keep_times = "0.01,0.05,0.1";
    /**
     * Doc234 D: T = T0 + dT*((r/R)^2 - 1/2); pair with --freeze-prescribed-t.
     * dT<=0 disables. Default off.
     */
    Real doc234_fixed_t_radial_dT = Real(0);
    /** Doc234 D: dump Position/Velocity/T CSV at t=0 and listed times (comma seconds). */
    std::string doc234_particle_field_csv_prefix;
    std::string doc234_field_dump_times = "0.1,0.5,1.0";
    /**
     * Doc232 H: after CLI finalize, print heat_source_mode / em_status and exit (no SPH).
     * Used for prescribed-heat conflict / invalid-value software regression.
     */
    bool doc232_h_cli_selftest_exit = false;
    Real probe_kernel_h_factor = Real(1.3);
    /**
     * Doc230: map client NC parametric melt into french geometry knobs (Reload from
     * test_3d_ophelie_client_nc_glass_relax). Prefer pairing with --client-nc-material.
     * Do not combine with bare --glass-height= (that resets center to H/2).
     */
    bool client_nc_geometry = false;
    /** Doc230: install client development materials (frozen μ @ 1000°C, σ16 borrowed). */
    bool client_nc_material = false;
    /** Optional override when --client-nc-geometry; NaN => development 0.180 m. */
    Real client_melt_height_m = std::numeric_limits<Real>::quiet_NaN();
    /**
     * Doc245: independent μ/k/cp dataset switch (does not change ρ/β/σ/power/geometry).
     * LEGACY_APPROX_TABLE | JAC_THESIS_RAW_BORROWED | JAC_THESIS_SHAPE_CEP_ANCHORED
     * Empty = leave whatever literature currently holds (usually CES Fig.2 via french default).
     */
    std::string material_dataset;
    /** Solver tables dir for JAC_THESIS_* (same CSVs as Doc239). */
    std::string material_tables_dir = "docs/ophelie/reference/doc239_runs/solver_tables";
    /** Directory for material_effective.json / material_probe.csv / curve compare (optional). */
    std::string material_evidence_dir;
    /**
     * Doc246 Q1_UPPER: after EM→fluid Q sample, apply
     *   Q_ctrl = C * Q_EM * exp[kappa*(z/H−0.5)]
     * with C from volume integral so ΣQ_ctrl V = target absorbed power.
     * kappa<=0 disables (default). Diagnostic only — not paper EM source.
     */
    Real doc246_q_upper_kappa = Real(0);
    FrenchLiteratureParameters literature = makeFrenchNaturalLiteratureParameters();
    FrenchSkullResistanceParams skull;
    FrenchFreeCrustParams free_crust;
    FrenchSkullShellParams skull_shell;
};

inline int applyLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, LocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
        }
        else if (std::strcmp(av[i], "--client-nc-geometry") == 0)
        {
            cli.client_nc_geometry = true;
        }
        else if (std::strcmp(av[i], "--client-nc-material") == 0)
        {
            cli.client_nc_material = true;
        }
        else if (std::strncmp(av[i], "--material-dataset=", 19) == 0)
        {
            cli.material_dataset = std::string(av[i] + 19);
        }
        else if (std::strncmp(av[i], "--material-tables-dir=", 22) == 0)
        {
            cli.material_tables_dir = std::string(av[i] + 22);
        }
        else if (std::strncmp(av[i], "--material-evidence-dir=", 24) == 0)
        {
            cli.material_evidence_dir = std::string(av[i] + 24);
        }
        else if (std::strncmp(av[i], "--doc246-q-upper-kappa=", 23) == 0)
        {
            cli.doc246_q_upper_kappa = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strncmp(av[i], "--melt-height=", 14) == 0)
        {
            // Consumed with --client-nc-geometry (same knob as client relax).
            cli.client_melt_height_m = static_cast<Real>(std::atof(av[i] + 14));
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
            cli.target_power_flag_seen = true;
            cli.enable_target_power = true;
            cli.target_power_w = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strncmp(av[i], "--prescribed-heat-power=", 24) == 0)
        {
            // Doc232 H: arm mode by flag presence; validate value (no silent EM fallback).
            cli.prescribed_heat_mode = true;
            const char *val = av[i] + 24;
            if (val[0] == '\0')
            {
                std::cerr << "[ophelie][doc232][H] FATAL: --prescribed-heat-power= missing value\n";
                return 2;
            }
            char *end = nullptr;
            const double parsed = std::strtod(val, &end);
            if (end == val || (end && *end != '\0') || !std::isfinite(parsed) || parsed < 0.0)
            {
                std::cerr << "[ophelie][doc232][H] FATAL: --prescribed-heat-power invalid '" << val
                          << "' (need finite >=0; use --zero-q for isothermal Q≡0)\n";
                return 2;
            }
            cli.prescribed_heat_power_w = static_cast<Real>(parsed);
            // Defer conflict/policy apply to post-parse finalize (order-independent).
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.enable_self_induction = false;
            cli.literature.enable_aind = false;
        }
        else if (std::strcmp(av[i], "--doc217-d2-local-poly") == 0)
        {
            cli.doc217_d2_local_poly = true;
        }
        else if (std::strcmp(av[i], "--spatial-frozen-q") == 0)
        {
            cli.spatial_frozen_q = true;
        }
        else if (std::strcmp(av[i], "--no-spatial-frozen-q") == 0)
        {
            cli.spatial_frozen_q = false;
        }
        else if (tryApplyFrenchLiteratureCli(av[i], cli.literature, cli.thermal_diffusion_explicit))
        {
            if (std::strcmp(av[i], "--zero-q") == 0)
            {
                cli.zero_q_flag_seen = true;
            }
            else if (std::strncmp(av[i], "--target-absorbed-power=", 24) == 0)
            {
                cli.target_power_flag_seen = true;
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
        else if (std::strncmp(av[i], "--wall-rim-height=", 18) == 0)
        {
            // Doc221 H alias of --wall-rim= (splash-guard height above melt).
            cli.wall_rim_height = static_cast<Real>(std::atof(av[i] + 18));
        }
        else if (std::strncmp(av[i], "--hydro-settle-time=", 20) == 0)
        {
            cli.hydro_settle_time = static_cast<Real>(std::atof(av[i] + 20));
        }
        else if (std::strncmp(av[i], "--hydro-settle-u-quiet=", 23) == 0)
        {
            cli.hydro_settle_u_quiet = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strncmp(av[i], "--hydro-settle-quiet-steps=", 27) == 0)
        {
            cli.hydro_settle_quiet_steps = static_cast<size_t>(std::strtoull(av[i] + 27, nullptr, 10));
        }
        else if (std::strncmp(av[i], "--hydro-settle-min-time=", 24) == 0)
        {
            cli.hydro_settle_min_time = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--hydro-settle-max-adv-steps=", 29) == 0)
        {
            cli.hydro_settle_max_adv_steps = static_cast<size_t>(std::strtoull(av[i] + 29, nullptr, 10));
        }
        else if (std::strncmp(av[i], "--hydro-settle-dt-floor=", 24) == 0)
        {
            cli.hydro_settle_dt_floor = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--hydro-settle-dt-collapse-steps=", 33) == 0)
        {
            cli.hydro_settle_dt_collapse_steps = static_cast<size_t>(std::strtoull(av[i] + 33, nullptr, 10));
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
        else if (std::strncmp(av[i], "--init-T-rz-csv=", 16) == 0)
        {
            cli.init_T_rz_csv = std::string(av[i] + 16);
        }
        else if (std::strcmp(av[i], "--init-T-rz-demean") == 0)
        {
            std::cerr
                << "[ophelie][doc225] FATAL: --init-T-rz-demean silently shifted absolute Temperature "
                   "(ABSOLUTE_T_SHIFTED_LEGACY). Mainline keeps CSV T. Re-run without demean, or use "
                   "--init-T-rz-demean-legacy only to reproduce Doc224 shifted cases (not for acceptance).\n";
            std::exit(2);
        }
        else if (std::strcmp(av[i], "--init-T-rz-demean-legacy") == 0)
        {
            cli.init_T_rz_demean_legacy = true;
        }
        else if (std::strncmp(av[i], "--t-force-ref=", 14) == 0)
        {
            cli.t_force_ref_override = static_cast<Real>(std::atof(av[i] + 14));
        }
        else if (std::strcmp(av[i], "--t-force-ref-from-mass-mean") == 0)
        {
            cli.t_force_ref_from_mass_mean = true;
            cli.t_force_ref_from_mass_mean_user_set = true;
        }
        else if (std::strcmp(av[i], "--no-t-force-ref-from-mass-mean") == 0)
        {
            cli.t_force_ref_from_mass_mean = false;
            cli.t_force_ref_from_mass_mean_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc225-warm-spatial-q-preflight") == 0)
        {
            cli.doc225_warm_spatial_q_preflight = true;
        }
        else if (std::strcmp(av[i], "--freeze-prescribed-t") == 0)
        {
            cli.freeze_prescribed_t_flag_seen = true;
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
            cli.dens_sum_normalize_user_set = true;
        }
        else if (std::strcmp(av[i], "--no-dens-sum-normalize") == 0 || std::strcmp(av[i], "--nodens") == 0 ||
                 std::strcmp(av[i], "--no-extra-normalize") == 0)
        {
            // Doc105: nodens / no-extra-normalize = disable extra dens_sum_normalize only
            // (DensitySummationCK + FreeSurface regularization still run).
            cli.dens_sum_normalize = false;
            cli.dens_sum_normalize_user_set = true;
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
        else if (std::strcmp(av[i], "--doc220-s0-sort-audit") == 0)
        {
            cli.doc220_s0_sort_audit = true;
        }
        else if (std::strcmp(av[i], "--doc220-h0-hydro-probe") == 0)
        {
            cli.doc220_h0_hydro_probe = true;
            cli.doc105_p2_hydro_stage_probe = true;
        }
        else if (std::strcmp(av[i], "--doc220-pprime-hydro") == 0)
        {
            // Coherent free-surface p' hydro: p_h in BackgroundPressure, momentum on p',
            // no GravityForceCK, wall dens on p' with a_ref=0.
            cli.doc220_pprime_hydro = true;
            cli.background_hydro_p = true;
            cli.omit_gravity_force = true;
            cli.perturbation_momentum = 1;
            cli.acoustic_dens_mode = 1;
            cli.acoustic_dens_mode_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc222-totalp-bgp-hydro") == 0)
        {
            // Doc222 E2: total-p + BackgroundPressure p_h + GravityForceCK (PM=0).
            // Flat rest: ∇p_h cancels g. FS restoring requires --doc222-ph-local-eta
            // (fixed z_free p_h has no radial ∇). Not the p' omit-g path (E2b falsified).
            cli.doc220_pprime_hydro = false;
            cli.background_hydro_p = true;
            cli.omit_gravity_force = false;
            cli.perturbation_momentum = 0;
            cli.acoustic_dens_mode = 1;
            cli.acoustic_dens_mode_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc221-legacy-omit-clears-forceprior") == 0)
        {
            cli.doc221_legacy_omit_clears_forceprior = true;
        }
        else if (std::strncmp(av[i], "--doc221-force-ledger-probe=", 28) == 0)
        {
            cli.doc221_force_ledger_probe_adv = std::atoi(av[i] + 28);
        }
        else if (std::strcmp(av[i], "--doc221-force-lifecycle-selftest") == 0)
        {
            cli.doc221_force_lifecycle_selftest = true;
        }
        else if (std::strcmp(av[i], "--diag-zero-net-uz") == 0)
        {
            // Doc221 K: volume-closure diagnostic (not literature open FS).
            cli.diag_zero_net_uz = true;
        }
        else if (std::strcmp(av[i], "--diag-hard-lid") == 0)
        {
            // Doc221 L: kinematic impermeable lid at z_free (not literature open FS).
            cli.diag_hard_lid = true;
        }
        else if (std::strcmp(av[i], "--diag-flat-lid") == 0)
        {
            // Doc221 L convenience: hard lid + zero-net-Uz volume closure.
            cli.diag_hard_lid = true;
            cli.diag_zero_net_uz = true;
        }
        else if (std::strcmp(av[i], "--diag-soft-lid") == 0)
        {
            // Doc221 N: K skeleton + thin-band kill-up, no z clamp (not literature open FS).
            cli.diag_soft_lid = true;
            cli.diag_soft_lid_tip = false;
            cli.diag_tip_soft_force = false;
            cli.diag_tip_fs_dashpot = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false; // never combine with hard clamp
        }
        else if (std::strcmp(av[i], "--diag-soft-lid-tip") == 0)
        {
            // Doc221 N2: K + tip-local kill-up only (r < tip_r), no full-band, no z clamp.
            cli.diag_soft_lid_tip = true;
            cli.diag_soft_lid = false;
            cli.diag_tip_soft_force = false;
            cli.diag_tip_fs_dashpot = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false;
            // keep existing uz_scale (default 0 = kill); N3 sets scale via dedicated flag.
        }
        else if (std::strcmp(av[i], "--diag-soft-lid-tip-dampen") == 0)
        {
            // Doc221 N3: tip-local dampen upward Uz (scale=0.5), not hard kill.
            cli.diag_soft_lid_tip = true;
            cli.diag_soft_lid = false;
            cli.diag_tip_soft_force = false;
            cli.diag_tip_fs_dashpot = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid_tip_uz_scale = Real(0.5);
        }
        else if (std::strncmp(av[i], "--diag-soft-lid-tip-r=", 22) == 0)
        {
            // Prefix length is 22 ("--diag-soft-lid-tip-r="); old 21 made atof("=0.08")==0
            // and the TinyReal fallback silently reset tip_r→0.10 (N2r scan invalid).
            cli.diag_soft_lid_tip_r = static_cast<Real>(std::atof(av[i] + 22));
        }
        else if (std::strncmp(av[i], "--diag-soft-lid-tip-uz-scale=", 29) == 0)
        {
            cli.diag_soft_lid_tip_uz_scale = static_cast<Real>(std::atof(av[i] + 29));
            // Enabling a non-trivial scale implies tip path (N3-style).
            cli.diag_soft_lid_tip = true;
            cli.diag_soft_lid = false;
            cli.diag_tip_soft_force = false;
            cli.diag_tip_fs_dashpot = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false;
        }
        else if (std::strcmp(av[i], "--diag-tip-soft-force") == 0)
        {
            // Doc221 N4: K + tip ForcePrior spring/dashpot; no host Uz edit.
            cli.diag_tip_soft_force = true;
            cli.diag_tip_fs_dashpot = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false;
            cli.diag_wall_flat_lid = false;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
        }
        else if (std::strcmp(av[i], "--diag-tip-fs-dashpot") == 0)
        {
            // Doc221 N5: K + tip dashpot only on FS-incomplete particles (DensitySummation<rho0).
            cli.diag_tip_fs_dashpot = true;
            cli.diag_tip_soft_force = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
            cli.diag_zero_net_uz = true;
            cli.diag_hard_lid = false;
            cli.diag_wall_flat_lid = false;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_tip_soft_force_k = Real(0); // spring abandoned
        }
        else if (std::strncmp(av[i], "--diag-tip-soft-force-r=", 24) == 0)
        {
            cli.diag_tip_soft_force_r = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--diag-tip-soft-force-k=", 24) == 0)
        {
            cli.diag_tip_soft_force_k = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--diag-tip-soft-force-c=", 24) == 0)
        {
            cli.diag_tip_soft_force_c = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strcmp(av[i], "--diag-tip-soft-force-legacy-accumulate") == 0)
        {
            // Doc222 B0: old ForcePrior += m*a without previous revoke (bug repro).
            cli.diag_tip_soft_force_legacy_accumulate = true;
        }
        else if (std::strcmp(av[i], "--doc222-tip-force-selftest") == 0)
        {
            cli.doc222_tip_force_selftest = true;
        }
        else if (std::strcmp(av[i], "--doc223-force-reset-probe") == 0)
        {
            cli.doc223_force_reset_probe = true;
        }
        else if (std::strcmp(av[i], "--doc223-paper-flat-slip-top") == 0)
        {
            cli.doc223_paper_flat_slip_top = true;
        }
        else if (std::strcmp(av[i], "--doc223-ablate-no-free-slip") == 0)
        {
            cli.doc223_ablate_no_free_slip = true;
        }
        else if (std::strcmp(av[i], "--doc224-lid-keep-dens-jump") == 0)
        {
            cli.doc224_lid_keep_dens_jump = true;
            cli.doc224_lid_keep_dens_jump_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc224-no-lid-keep-dens-jump") == 0)
        {
            cli.doc224_lid_keep_dens_jump = false;
            cli.doc224_lid_keep_dens_jump_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc224-ablate-opencup") == 0)
        {
            cli.doc224_ablate_opencup = true;
        }
        else if (std::strcmp(av[i], "--doc224-ablate-acoustic-freeslip") == 0)
        {
            cli.doc224_ablate_acoustic_freeslip = true;
            cli.doc224_ablate_acoustic_freeslip_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc224-acoustic-freeslip") == 0)
        {
            // Legacy Doc223 P2: tag AcousticWallFreeSlip (known lid penetration under F2).
            cli.doc224_ablate_acoustic_freeslip = false;
            cli.doc224_ablate_acoustic_freeslip_user_set = true;
        }
        else if (std::strncmp(av[i], "--doc224-u-blowup-gate=", 23) == 0)
        {
            char *end = nullptr;
            const double parsed = std::strtod(av[i] + 23, &end);
            if (end != av[i] + 23 && std::isfinite(parsed) && parsed > 0.0)
            {
                cli.doc224_u_blowup_gate = static_cast<Real>(parsed);
            }
        }
        else if (std::strcmp(av[i], "--doc223-e2-lid-slip-probe") == 0)
        {
            cli.doc223_e2_lid_slip_probe = true;
        }
        else if (std::strcmp(av[i], "--doc226-top-operator-probe") == 0)
        {
            cli.doc226_top_operator_probe = true;
        }
        else if (std::strcmp(av[i], "--doc226-local-particle-viscosity") == 0)
        {
            cli.doc226_local_particle_viscosity = true;
        }
        else if (std::strcmp(av[i], "--doc226-mu-T-only") == 0)
        {
            cli.doc226_mu_T_only = true;
        }
        else if (std::strncmp(av[i], "--doc226-mu-cap=", 16) == 0)
        {
            cli.doc226_mu_cap = std::strtod(av[i] + 16, nullptr);
        }
        else if (std::strcmp(av[i], "--doc226-e1-local-mu-probe") == 0)
        {
            cli.doc226_e1_local_mu_probe = true;
        }
        else if (std::strcmp(av[i], "--doc227-actual-timestep") == 0)
        {
            cli.doc227_actual_timestep = true;
        }
        else if (std::strncmp(av[i], "--doc227-step-audit=", 20) == 0)
        {
            cli.doc227_step_audit_csv = av[i] + 20;
        }
        else if (std::strcmp(av[i], "--doc227-schedule-selftest") == 0)
        {
            cli.doc227_schedule_selftest = true;
        }
        else if (std::strcmp(av[i], "--doc227-c-uniform-mu-probe") == 0)
        {
            cli.doc227_c_uniform_mu_probe = true;
        }
        else if (std::strcmp(av[i], "--doc223-e3-normal-perturb") == 0)
        {
            cli.doc223_e3_normal_perturb = true;
        }
        else if (std::strcmp(av[i], "--doc223-f2-annular-hot") == 0)
        {
            cli.doc223_f2_annular_hot = true;
        }
        else if (std::strncmp(av[i], "--doc222-fs-deform-eps=", 23) == 0)
        {
            cli.doc222_fs_deform_eps = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strcmp(av[i], "--doc222-ph-signed") == 0)
        {
            cli.doc222_ph_signed = true;
        }
        else if (std::strcmp(av[i], "--doc222-ph-clamp-legacy") == 0)
        {
            // Doc222 E1 A/B: restore max(0,z_free-z) truncation (pre-fix D2 FAIL path).
            cli.doc222_ph_signed = false;
        }
        else if (std::strcmp(av[i], "--doc222-c3-ph-compare") == 0)
        {
            cli.doc222_c3_ph_compare = true;
            cli.doc220_pprime_hydro = true;
            cli.background_hydro_p = true;
            cli.perturbation_momentum = 1;
            cli.omit_gravity_force = true;
            if (!cli.acoustic_dens_mode_user_set)
            {
                cli.acoustic_dens_mode = 1;
                cli.acoustic_dens_mode_user_set = true;
            }
            if (cli.doc222_fs_deform_eps <= TinyReal)
                cli.doc222_fs_deform_eps = Real(0.02);
        }
        else if (std::strcmp(av[i], "--doc222-deform-rho-freeze") == 0)
        {
            // Legacy: freeze Density=ρ0 after deform (C3 showed p_mom≡0).
            cli.doc222_deform_skip_rho_freeze = false;
        }
        else if (std::strcmp(av[i], "--doc222-deform-no-rho-freeze") == 0)
        {
            cli.doc222_deform_skip_rho_freeze = true;
        }
        else if (std::strcmp(av[i], "--doc222-ph-local-eta") == 0)
        {
            // Doc222 E2b: p_h from local free-surface η(r).
            cli.doc222_ph_local_eta = true;
        }
        else if (std::strcmp(av[i], "--doc222-ph-fixed-zfree") == 0)
        {
            cli.doc222_ph_local_eta = false;
        }
        else if (std::strncmp(av[i], "--doc222-ph-local-eta-bins=", 27) == 0)
        {
            cli.doc222_ph_local_eta_bins = std::max(4, std::atoi(av[i] + 27));
            cli.doc222_ph_local_eta = true;
        }
        else if (std::strcmp(av[i], "--diag-wall-flat-lid") == 0)
        {
            // Doc221 O: physical closed crucible (flat lid wall contact); not literature open FS.
            cli.diag_wall_flat_lid = true;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
        }
        else if (std::strcmp(av[i], "--diag-lid-free-slip") == 0)
        {
            // Doc221 P1: free-slip from t=0 (legacy; noisy — prefer post-settle).
            cli.diag_lid_free_slip = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
        }
        else if (std::strcmp(av[i], "--diag-lid-free-slip-post-settle") == 0)
        {
            // Doc221 P2: settle no-slip lid, then zero U and arm free-slip for keep.
            cli.diag_lid_free_slip = true;
            cli.diag_lid_free_slip_post_settle = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
        }
        else if (std::strcmp(av[i], "--diag-wall-lid-viscous-omit") == 0)
        {
            // Doc221 O3: ClosedCrucible acoustic + viscous OpenCup; NO AcousticWallFreeSlip.
            cli.diag_wall_lid_viscous_omit = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
        }
        else if (std::strcmp(av[i], "--diag-lid-dens-silence") == 0)
        {
            // Doc221 O4: O3 + dens_mode=3 + lid ContSilent/OmitPJump; NO AcousticWallFreeSlip.
            cli.diag_lid_dens_silence = true;
            cli.diag_lid_cont_silent = false;
            cli.diag_lid_cont_silent_only = false;
            cli.diag_lid_omit_pjump_only = false;
            cli.diag_wall_lid_viscous_omit = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
            cli.acoustic_dens_mode = 3;
            cli.acoustic_dens_mode_user_set = true;
        }
        else if (std::strcmp(av[i], "--diag-lid-cont-silent") == 0)
        {
            // Doc221 O5: O3 + lid ContSilent/OmitPJump; dens_mode stays 1 (no dens_mode=3).
            cli.diag_lid_cont_silent = true;
            cli.diag_lid_cont_silent_only = false;
            cli.diag_lid_omit_pjump_only = false;
            cli.diag_lid_dens_silence = false;
            cli.diag_wall_lid_viscous_omit = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
            // Do NOT force dens_mode=3; pprime / explicit --acoustic-dens-mode=1 keep mode 1.
        }
        else if (std::strcmp(av[i], "--diag-lid-cont-silent-only") == 0)
        {
            // Doc221 O6a: O3 + ContSilent only; dens_mode=1.
            cli.diag_lid_cont_silent_only = true;
            cli.diag_lid_omit_pjump_only = false;
            cli.diag_lid_cont_silent = false;
            cli.diag_lid_dens_silence = false;
            cli.diag_wall_lid_viscous_omit = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
        }
        else if (std::strcmp(av[i], "--diag-lid-omit-pjump-only") == 0)
        {
            // Doc221 O6b: O3 + OmitPJump only; dens_mode=1.
            cli.diag_lid_omit_pjump_only = true;
            cli.diag_lid_cont_silent_only = false;
            cli.diag_lid_cont_silent = false;
            cli.diag_lid_dens_silence = false;
            cli.diag_wall_lid_viscous_omit = true;
            cli.diag_wall_flat_lid = true;
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.diag_hard_lid = false;
            cli.diag_soft_lid = false;
            cli.diag_soft_lid_tip = false;
        }
        else if (std::strcmp(av[i], "--no-doc220-pprime-hydro") == 0)
        {
            cli.doc220_pprime_hydro = false;
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
            cli.dens_sum_normalize_user_set = true;
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
        else if (std::strcmp(av[i], "--doc219-hydro-profile") == 0)
        {
            cli.doc219_hydro_profile = true;
        }
        else if (std::strcmp(av[i], "--no-doc219-hydro-profile") == 0)
        {
            cli.doc219_hydro_profile = false;
        }
        else if (std::strcmp(av[i], "--kernel-correction") == 0)
        {
            cli.kernel_correction = true;
            cli.kernel_correction_user_set = true;
        }
        else if (std::strcmp(av[i], "--no-kernel-correction") == 0)
        {
            cli.kernel_correction = false;
            cli.kernel_correction_user_set = true;
        }
        else if (std::strcmp(av[i], "--tvf") == 0)
        {
            cli.enable_tvf = true;
            cli.enable_tvf_user_set = true;
        }
        else if (std::strcmp(av[i], "--no-tvf") == 0)
        {
            cli.enable_tvf = false;
            cli.enable_tvf_user_set = true;
        }
        else if (std::strncmp(av[i], "--acoustic-dens-mode=", 21) == 0)
        {
            cli.acoustic_dens_mode = std::atoi(av[i] + 21);
            cli.acoustic_dens_mode_user_set = true;
            // 0/1 legacy; 3 = Doc201 omit DissipativeUJump (Doc221 O4 lid dens silence).
            if (cli.acoustic_dens_mode != 0 && cli.acoustic_dens_mode != 1 && cli.acoustic_dens_mode != 3)
            {
                std::cerr << "[ophelie] --acoustic-dens-mode must be 0, 1, or 3\n";
                std::exit(1);
            }
        }
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
        {
            cli.summary_json = std::string(av[i] + 15);
        }
        else if (std::strncmp(av[i], "--probes-timeseries=", 20) == 0)
        {
            cli.probes_timeseries_csv = std::string(av[i] + 20);
            cli.probes_timeseries_user_set = true;
        }
        else if (std::strcmp(av[i], "--doc232-probe-path-verify") == 0)
        {
            cli.doc232_probe_path_verify = true;
        }
        else if (std::strcmp(av[i], "--doc232-probe-path-verify-exit") == 0)
        {
            cli.doc232_probe_path_verify = true;
            cli.doc232_probe_path_verify_exit = true;
        }
        else if (std::strncmp(av[i], "--doc232-dump-particles=", 24) == 0)
        {
            cli.doc232_dump_particles_csv = std::string(av[i] + 24);
        }
        else if (std::strcmp(av[i], "--doc232-adiabatic-heat") == 0)
        {
            cli.doc232_adiabatic_heat = true;
        }
        else if (std::strncmp(av[i], "--thermal-bc=", 13) == 0)
        {
            cli.thermal_bc_mode = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--h-side=", 9) == 0)
        {
            cli.thermal_h_side_set = true;
            cli.thermal_h_side = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--h-bottom=", 11) == 0)
        {
            cli.thermal_h_bottom_set = true;
            cli.thermal_h_bottom = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--h-free=", 9) == 0)
        {
            cli.thermal_h_free_set = true;
            cli.thermal_h_free = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--emissivity=", 13) == 0)
        {
            cli.thermal_emissivity_set = true;
            cli.thermal_emissivity = static_cast<Real>(std::atof(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--t-cool=", 9) == 0)
        {
            cli.thermal_t_cool_set = true;
            cli.thermal_t_cool = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strcmp(av[i], "--doc232-boundary-face-audit") == 0)
        {
            cli.doc232_boundary_face_audit = true;
        }
        else if (std::strcmp(av[i], "--doc232-boundary-face-audit-exit") == 0)
        {
            cli.doc232_boundary_face_audit = true;
            cli.doc232_boundary_face_audit_exit = true;
        }
        else if (std::strncmp(av[i], "--doc232-boundary-audit-dir=", 28) == 0)
        {
            cli.doc232_boundary_audit_dir = std::string(av[i] + 28);
        }
        else if (std::strncmp(av[i], "--boundary-face-t-mode=", 23) == 0)
        {
            // strlen("--boundary-face-t-mode=") == 23; offset 22 wrongly kept a leading '='.
            cli.boundary_face_t_mode = std::string(av[i] + 23);
            if (cli.boundary_face_t_mode != "global_dt_dz" && cli.boundary_face_t_mode != "band" &&
                cli.boundary_face_t_mode != "local_dt_dz")
            {
                std::cerr << "[ophelie][doc233][D] FATAL: --boundary-face-t-mode must be "
                             "global_dt_dz|local_dt_dz|band, got '"
                          << cli.boundary_face_t_mode << "'\n";
                return 2;
            }
        }
        else if (std::strcmp(av[i], "--doc233-boundary-scope-audit") == 0 ||
                 std::strcmp(av[i], "--doc233-boundary-scope-audit-exit") == 0)
        {
            cli.doc232_boundary_face_audit = true;
            cli.doc233_boundary_scope_audit = true;
            if (std::strcmp(av[i], "--doc233-boundary-scope-audit-exit") == 0)
                cli.doc232_boundary_face_audit_exit = true;
        }
        else if (std::strcmp(av[i], "--doc233-buoyancy-force-chain-exit") == 0 ||
                 std::strcmp(av[i], "--doc234-force-precision-exit") == 0)
        {
            cli.doc233_buoyancy_force_chain_exit = true;
            if (std::strcmp(av[i], "--doc234-force-precision-exit") == 0)
                cli.doc234_force_precision_exit = true;
            cli.literature.enable_boussinesq = true;
            // Ensure prescribe block runs; manufactured radial ΔT applied inside if isothermal.
            if (!(cli.prescribe_t_radial_dT > TinyReal) && !(cli.prescribe_t_hot_bottom_dT > TinyReal) &&
                cli.init_T_rz_csv.empty() && !cli.doc223_f2_annular_hot &&
                !(cli.doc234_fixed_t_radial_dT > TinyReal))
                cli.prescribe_t_radial_dT = Real(30);
        }
        else if (std::strncmp(av[i], "--doc234-force-precision-csv=", 29) == 0)
        {
            cli.doc234_force_precision_csv = std::string(av[i] + 29);
            cli.doc234_force_precision_exit = true;
            cli.doc233_buoyancy_force_chain_exit = true;
            cli.literature.enable_boussinesq = true;
        }
        else if (std::strncmp(av[i], "--doc235-off-event-dump-dir=", 28) == 0)
        {
            cli.doc235_off_event_dump_dir = std::string(av[i] + 28);
        }
        else if (std::strncmp(av[i], "--doc235-off-event-keep-times=", 30) == 0)
        {
            cli.doc235_off_event_keep_times = std::string(av[i] + 30);
        }
        else if (std::strncmp(av[i], "--doc234-fixed-t-radial-dT=", 27) == 0)
        {
            cli.doc234_fixed_t_radial_dT = static_cast<Real>(std::atof(av[i] + 27));
        }
        else if (std::strncmp(av[i], "--doc234-particle-field-csv-prefix=", 35) == 0)
        {
            cli.doc234_particle_field_csv_prefix = std::string(av[i] + 35);
        }
        else if (std::strncmp(av[i], "--doc234-field-dump-times=", 26) == 0)
        {
            cli.doc234_field_dump_times = std::string(av[i] + 26);
        }
        else if (std::strcmp(av[i], "--doc232-h-cli-selftest-exit") == 0)
        {
            cli.doc232_h_cli_selftest_exit = true;
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
    if (cli.doc223_paper_flat_slip_top)
    {
        const bool illegal =
            cli.diag_tip_soft_force || cli.diag_tip_fs_dashpot || cli.diag_soft_lid ||
            cli.diag_soft_lid_tip || cli.diag_hard_lid || cli.diag_zero_net_uz ||
            cli.diag_lid_dens_silence || cli.diag_lid_cont_silent ||
            cli.diag_lid_cont_silent_only || cli.diag_lid_omit_pjump_only ||
            cli.doc222_fs_deform_eps > TinyReal || cli.doc222_c3_ph_compare ||
            cli.doc222_ph_local_eta || cli.doc223_force_reset_probe;
        if (illegal)
        {
            std::cerr << "[ophelie][doc223][D] --doc223-paper-flat-slip-top is illegal with "
                         "tip/zero-net/z-clamp/ContSilent/fs-deform/local-eta/force-reset-probe\n";
            return 2;
        }
        cli.diag_wall_flat_lid = true;
        cli.diag_lid_free_slip = true;
        cli.diag_lid_free_slip_post_settle = true;
        cli.doc220_pprime_hydro = true;
        cli.background_hydro_p = true;
        cli.omit_gravity_force = true;
        cli.perturbation_momentum = 1;
        if (!cli.acoustic_dens_mode_user_set)
        {
            cli.acoustic_dens_mode = 1;
            cli.acoustic_dens_mode_user_set = true;
        }
        // Settle: Internal dens + packing normalize (absorb overpack into σ0/ρ).
        // Keep: do NOT re-run dens reinit every advection — Doc223 E1 showed that
        // periodic Internal+normalize snaps (Δρ∼O(10²)) are CFL-sensitive forcing
        // (halving acoustic OR advection both raised U_rms by ∼100×). Keep uses
        // acoustic continuity only after settle/flush dens.
        if (!cli.dens_sum_normalize_user_set)
        {
            cli.dens_sum_normalize = false;
            cli.settle_dens_sum_normalize = true;
        }
        // Doc223 E1: Doc219 default B-matrix amplifies ClosedCrucible side-wall noise
        // (~30×; keep_adv200 u_side 3e-4→3e-6). Mark user_set so Doc219 cannot re-enable.
        if (!cli.kernel_correction_user_set)
        {
            cli.kernel_correction = false;
            cli.kernel_correction_user_set = true;
        }
        if (cli.doc223_ablate_no_free_slip)
        {
            cli.diag_lid_free_slip = false;
            cli.diag_lid_free_slip_post_settle = false;
            cli.doc224_lid_keep_dens_jump = false;
            cli.doc224_ablate_acoustic_freeslip = false;
        }
        else
        {
            if (!cli.doc224_lid_keep_dens_jump_user_set)
            {
                // Retained for legacy --doc224-acoustic-freeslip path.
                cli.doc224_lid_keep_dens_jump = true;
            }
            if (!cli.doc224_ablate_acoustic_freeslip_user_set)
            {
                // Doc224 C3d: impermeable acoustic + OpenCup viscous (default paper-top).
                cli.doc224_ablate_acoustic_freeslip = true;
            }
        }
    }
    else if (cli.doc223_ablate_no_free_slip)
    {
        std::cerr << "[ophelie][doc223][E1] --doc223-ablate-no-free-slip requires "
                     "--doc223-paper-flat-slip-top\n";
        return 2;
    }
    if (cli.doc223_e2_lid_slip_probe)
    {
        if (!cli.doc223_paper_flat_slip_top || cli.doc223_ablate_no_free_slip)
        {
            std::cerr << "[ophelie][doc223][E2] --doc223-e2-lid-slip-probe requires "
                         "--doc223-paper-flat-slip-top (not ablate-no-free-slip)\n";
            return 2;
        }
    }
    if (cli.doc226_top_operator_probe)
    {
        if (!cli.doc223_paper_flat_slip_top || cli.doc223_ablate_no_free_slip)
        {
            std::cerr << "[ophelie][doc226][B] --doc226-top-operator-probe requires "
                         "--doc223-paper-flat-slip-top (not ablate-no-free-slip)\n";
            return 2;
        }
        if (cli.doc223_e2_lid_slip_probe || cli.doc223_e3_normal_perturb ||
            cli.doc226_e1_local_mu_probe)
        {
            std::cerr << "[ophelie][doc226][B] do not combine with E2/E3/E1 probes\n";
            return 2;
        }
    }
    if (cli.doc226_mu_T_only && !cli.doc226_local_particle_viscosity)
    {
        std::cerr << "[ophelie][doc226] --doc226-mu-T-only requires --doc226-local-particle-viscosity\n";
        return 2;
    }
    if (cli.doc227_c_uniform_mu_probe)
    {
        if (!cli.doc226_local_particle_viscosity)
        {
            std::cerr << "[ophelie][doc227][C] --doc227-c-uniform-mu-probe requires "
                         "--doc226-local-particle-viscosity\n";
            return 2;
        }
        if (cli.doc226_e1_local_mu_probe || cli.doc226_top_operator_probe)
        {
            std::cerr << "[ophelie][doc227][C] do not combine with E1/B probes\n";
            return 2;
        }
    }
    if (cli.doc226_e1_local_mu_probe)
    {
        if (!cli.doc226_local_particle_viscosity)
        {
            std::cerr << "[ophelie][doc226][E1] --doc226-e1-local-mu-probe requires "
                         "--doc226-local-particle-viscosity\n";
            return 2;
        }
        if (cli.doc226_top_operator_probe || cli.doc223_e2_lid_slip_probe ||
            cli.doc223_e3_normal_perturb || cli.doc227_c_uniform_mu_probe)
        {
            std::cerr << "[ophelie][doc226][E1] do not combine with B/E2/E3/C probes\n";
            return 2;
        }
    }
    if (cli.doc223_e3_normal_perturb)
    {
        if (!cli.doc223_paper_flat_slip_top || cli.doc223_ablate_no_free_slip)
        {
            std::cerr << "[ophelie][doc223][E3] --doc223-e3-normal-perturb requires "
                         "--doc223-paper-flat-slip-top (not ablate-no-free-slip)\n";
            return 2;
        }
        if (cli.doc223_e2_lid_slip_probe)
        {
            std::cerr << "[ophelie][doc223][E3] do not combine with --doc223-e2-lid-slip-probe\n";
            return 2;
        }
    }
    if (cli.doc223_f2_annular_hot)
    {
        if (!cli.doc223_paper_flat_slip_top)
        {
            std::cerr << "[ophelie][doc223][F2] --doc223-f2-annular-hot requires "
                         "--doc223-paper-flat-slip-top\n";
            return 2;
        }
        if (cli.prescribe_t_radial_dT > TinyReal || cli.prescribe_t_hot_bottom_dT > TinyReal ||
            !cli.init_T_rz_csv.empty())
        {
            std::cerr << "[ophelie][doc223][F2] annular-hot is exclusive with "
                         "--prescribe-t-radial / --prescribe-t-hot-bottom / --init-T-rz-csv\n";
            return 2;
        }
    }
    if (!cli.init_T_rz_csv.empty())
    {
        if (cli.prescribe_t_radial_dT > TinyReal || cli.prescribe_t_hot_bottom_dT > TinyReal)
        {
            std::cerr << "[ophelie][doc224][D] --init-T-rz-csv exclusive with "
                         "--prescribe-t-radial / --prescribe-t-hot-bottom\n";
            return 2;
        }
        if (!cli.doc223_paper_flat_slip_top)
        {
            std::cerr << "[ophelie][doc224][D] --init-T-rz-csv requires "
                         "--doc223-paper-flat-slip-top for paper-top fluid identity\n";
            return 2;
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
        cli.reload_dir = resolveDefaultFrenchReloadFolder(cli.client_nc_geometry);
    }
    // Doc232 H: prescribed-heat mode is armed by flag presence (incl. P=0); conflicts order-independent.
    if (cli.prescribed_heat_mode)
    {
        if (cli.zero_q_flag_seen || cli.literature.zero_q)
        {
            std::cerr << "[ophelie][doc232][H] FATAL: --prescribed-heat-power conflicts with --zero-q\n";
            return 2;
        }
        if (cli.freeze_prescribed_t_flag_seen || cli.freeze_prescribed_temperature)
        {
            std::cerr << "[ophelie][doc232][H] FATAL: --prescribed-heat-power conflicts with "
                         "--freeze-prescribed-t (warm path needs evolving T)\n";
            return 2;
        }
        if (cli.target_power_flag_seen)
        {
            std::cerr << "[ophelie][doc232][H] FATAL: --prescribed-heat-power conflicts with "
                         "--target-power=/--target-absorbed-power= (EM target path)\n";
            return 2;
        }
        cli.enable_target_power = false;
        cli.literature.zero_q = false;
        cli.literature.coupling = ThermoEMCouplingMode::Off; // mutual exclusion with online EM refresh
        cli.literature.enable_thermal_diffusion = true;
        cli.literature.target_glass_absorbed_power_w = cli.prescribed_heat_power_w;
        cli.target_power_w = cli.prescribed_heat_power_w;
        std::cout << "[ophelie][doc232][H] heat_source_mode=prescribed_uniform_volumetric"
                  << " P_abs=" << cli.prescribed_heat_power_w
                  << " W em_status=NOT_RUN (NOT generator "
                  << cli.literature.generator_power_w << " W; development absorbed power)\n";
    }
    if (cli.doc232_adiabatic_heat)
    {
        if (!cli.prescribed_heat_mode || !(cli.prescribed_heat_power_w > TinyReal))
        {
            std::cerr << "[ophelie][doc232][C] FATAL: --doc232-adiabatic-heat requires "
                         "--prescribed-heat-power=<P> with P>0\n";
            return 2;
        }
        // Override warm-path diffusion; adiabatic gate is production Q only.
        cli.literature.enable_thermal_diffusion = false;
        cli.literature.enable_boussinesq = false;
        cli.hydro_settle_time = Real(0);
        std::cout << "[ophelie][doc232][C] adiabatic production heat gate: BC/diffusion/buoyancy off; "
                     "hydro_settle forced 0; thermal-only keep\n";
    }
    if (cli.doc232_boundary_face_audit)
    {
        cli.literature.enable_boussinesq = false;
        cli.hydro_settle_time = Real(0);
        cli.literature.enable_thermal_diffusion = false;
        if (cli.doc232_boundary_face_audit_exit && !(cli.literature.zero_q) && !cli.prescribed_heat_mode)
        {
            // Face audit needs particles + BC only; force Q off to skip EM.
            cli.literature.zero_q = true;
            cli.enable_target_power = false;
        }
        std::cout << "[ophelie][doc232][D] boundary face audit: settle/diffusion/buoyancy off"
                  << (cli.doc232_boundary_face_audit_exit ? "; exit after gate" : "") << "\n";
    }
    if (cli.literature.zero_q)
    {
        // Isothermal / Q-off hydro: do not calibrate coil to French 50 kW.
        cli.enable_target_power = false;
    }
    if (cli.enable_target_power)
    {
        cli.literature.target_glass_absorbed_power_w = cli.target_power_w;
    }
    if (cli.prescribe_t_hot_bottom_dT > TinyReal || cli.prescribe_t_radial_dT > TinyReal ||
        cli.doc223_f2_annular_hot || !cli.init_T_rz_csv.empty() ||
        cli.doc234_fixed_t_radial_dT > TinyReal)
    {
        // Freeze-T buoyancy smoke: hold ΔT, no EM heating / diffusion (pass --zero-q explicitly too).
        // Thaw path (no --freeze-prescribed-t): prescribe is IC only; Q + diffusion + Robin evolve.
        if (cli.freeze_prescribed_temperature)
        {
            cli.literature.zero_q = true;
            cli.literature.enable_thermal_diffusion = false;
            cli.thermal_diffusion_explicit = false;
        }
    }
    finalizeFrenchPeriodicLiteratureDefaults(cli.literature, cli.thermal_diffusion_explicit);

    // Doc225 A: warm Spatial-Q preflight (auto when thaw + CSV + Q on, or explicit flag).
    const bool warm_csv_q = !cli.init_T_rz_csv.empty() && !cli.freeze_prescribed_temperature &&
                            !cli.literature.zero_q;
    if (cli.doc225_warm_spatial_q_preflight || warm_csv_q)
    {
        auto fail = [](const char *msg) {
            std::cerr << "[ophelie][doc225][warm-preflight] FATAL: " << msg << std::endl;
            return 2;
        };
        if (cli.freeze_prescribed_temperature)
            return fail("freeze_T must be 0 for warm Spatial-Q");
        if (!cli.literature.enable_thermal_diffusion)
            return fail("thermal_diffusion must be ON (--enable-thermal-diffusion)");
        // Explicit --init-T-rz-demean-legacy is allowed for F0 read-only probe only; never acceptance.
        if (cli.init_T_rz_demean_legacy)
        {
            std::cout << "[ophelie][doc225][warm-preflight] WARN ABSOLUTE_T_SHIFTED_LEGACY probe path "
                         "(explicit --init-T-rz-demean-legacy); NOT for acceptance\n";
        }
        if (cli.diag_soft_lid_tip || cli.diag_tip_soft_force || cli.diag_tip_fs_dashpot ||
            cli.diag_zero_net_uz || cli.diag_lid_cont_silent)
            return fail("diagnostic velocity interventions must be OFF");
        if (!cli.prescribed_heat_mode &&
            !(cli.target_power_w > Real(1) && std::abs(cli.target_power_w - Real(50000)) <= Real(50)))
        {
            std::cerr << "[ophelie][doc225][warm-preflight] WARN target_power_w=" << cli.target_power_w
                      << " (expected ~50000)\n";
        }
        std::cout << "[ophelie][doc225][warm-preflight] PASS freeze=0 diffusion=1 legacy_demean="
                  << (cli.init_T_rz_demean_legacy ? 1 : 0)
                  << " no tip/Uz-kill ContSilent=0 Q_target=" << cli.target_power_w << std::endl;
    }
    return 0;
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

/** Doc235 B: keep Pressure consistent with Density after dens flush (p = c0²(ρ−ρ0)). */
inline void hostSyncEosPressureFromDensity(BaseParticles &particles, Real rho0, Real c0)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Pressure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *pres = particles.getVariableDataByName<Real>("Pressure");
    const Real c2 = c0 * c0;
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        pres[i] = c2 * (rho[i] - rho0);
    syncVariableToDevice<Real>(particles, "Pressure");
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

/**
 * Doc227 B5: deterministic schedule truncation (no fluid). Verifies last acoustic
 * substep is min(candidate, remaining) for Dt ≷ / non-integer multiple of dt_ac.
 */
inline bool hostDoc227ScheduleSelftest(std::string &reason)
{
    auto run_case = [](Real Dt, Real dt_ac_fixed, Real end_cap) -> std::pair<Real, int> {
        Real integrated = 0;
        int n = 0;
        constexpr Real eps = Real(1.0e-14);
        while (integrated + eps < Dt)
        {
            const Real rem = Dt - integrated;
            const Real dt = std::min({dt_ac_fixed, rem, end_cap});
            if (dt <= TinyReal)
                break;
            if (dt > rem + Real(1.0e-12) || dt > dt_ac_fixed + Real(1.0e-12))
                return {-1, -1};
            integrated += dt;
            ++n;
        }
        return {integrated, n};
    };
    {
        const auto r = run_case(Real(1.0e-3), Real(7.0e-4), Real(1.0)); // ac > rem on last
        if (r.second < 0 || std::abs(r.first - Real(1.0e-3)) > Real(1.0e-12) || r.second != 2)
        {
            reason = "case_ac_gt_remaining_failed";
            return false;
        }
    }
    {
        const auto r = run_case(Real(2.0e-3), Real(5.0e-4), Real(1.0)); // ac < Dt
        if (r.second < 0 || std::abs(r.first - Real(2.0e-3)) > Real(1.0e-12) || r.second != 4)
        {
            reason = "case_ac_lt_Dt_failed";
            return false;
        }
    }
    {
        const auto r = run_case(Real(1.25e-3), Real(5.0e-4), Real(1.0)); // non-integer
        if (r.second < 0 || std::abs(r.first - Real(1.25e-3)) > Real(1.0e-12) || r.second != 3)
        {
            reason = "case_noninteger_Dt_failed";
            return false;
        }
    }
    reason = "ok";
    return true;
}

/** Doc227 keep-loop schedule compare: Real is float; 1e-12 is below ulp near 1e-3. */
inline Real doc227ScheduleTol(Real scale)
{
    const Real a = std::max(std::abs(scale), Real(1.0e-6));
    const Real eps = std::numeric_limits<Real>::epsilon();
    return std::max({Real(1.0e-6), Real(1.0e-4) * a, Real(64) * eps * std::max(a, Real(1))});
}

/**
 * Doc226 F: thin BaseDynamics wrapper so ViscousForce wall can sit in
 * addPostContactInteraction while still switching Closed↔OpenCup after settle.
 * Couette / Viscosity+post-contact need wall BEFORE ForcePriorCK update.
 */
class OphelieCallFnDynamicsCK : public BaseDynamics<void>
{
    std::function<void()> fn_;

  public:
    explicit OphelieCallFnDynamicsCK(std::function<void()> fn) : fn_(std::move(fn)) {}
    void exec(Real dt = 0.0) override
    {
        (void)dt;
        fn_();
    }
};

/** Doc223: ForcePriorCK channels that must be cleared as one transaction with ForcePrior. */
inline const char *const kDoc223ForcePriorCurrentNames[] = {
    "GravityForceCK", "BoussinesqBuoyancyForceCK", "ViscousForce", "TipSoftForce"};

inline void hostDoc223ZeroAllForcePriorChannels(BaseParticles &particles)
{
    hostZeroVecFieldFrench(particles, "ForcePrior");
    for (const char *cur : kDoc223ForcePriorCurrentNames)
    {
        hostZeroVecFieldFrench(particles, cur);
        hostZeroVecFieldFrench(particles, std::string("Previous") + cur);
    }
}

struct Doc223ForceFieldMoments
{
    Real mean_z = 0;
    Real rms = 0;
    Real max_abs = 0;
    Real a_z_mean = 0; // mass-weighted Fz/m
};

inline Doc223ForceFieldMoments hostDoc223ForceFieldMoments(BaseParticles &particles, const std::string &name)
{
    Doc223ForceFieldMoments m;
    syncVariableToHost<Vecd>(particles, name);
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    Real s2 = 0, sz = 0, fz_sum = 0, msum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        s2 += v[i].squaredNorm();
        sz += v[i][2];
        m.max_abs = std::max(m.max_abs, v[i].norm());
        fz_sum += v[i][2];
        msum += mass[i];
    }
    const Real inv = n > 0 ? Real(1) / Real(n) : Real(0);
    m.mean_z = sz * inv;
    m.rms = std::sqrt(s2 * inv);
    m.a_z_mean = (msum > TinyReal) ? (fz_sum / msum) : Real(0);
    return m;
}

inline void hostDoc223LogForceStage(const std::string &tag, BaseParticles &particles, Real g)
{
    const auto fp = hostDoc223ForceFieldMoments(particles, "ForcePrior");
    const auto grav = hostDoc223ForceFieldMoments(particles, "GravityForceCK");
    const auto prev = hostDoc223ForceFieldMoments(particles, "PreviousGravityForceCK");
    const auto visc = hostDoc223ForceFieldMoments(particles, "ViscousForce");
    std::cout << "[ophelie][doc223][B] " << tag
              << " a_fp_z=" << fp.a_z_mean << " (expect ~" << -g << ")"
              << " fp_z_mean=" << fp.mean_z << " fp_rms=" << fp.rms
              << " grav_z_mean=" << grav.mean_z << " prevG_z_mean=" << prev.mean_z
              << " visc_rms=" << visc.rms << "\n";
}

/** Doc223 E2: wall-viscous a=F/m band stats under imposed tangential U. */
struct Doc223E2ViscBand
{
    size_t n = 0;
    Real a_tx_rms = 0; // Fx/m
    Real a_tz_rms = 0; // Fz/m
    Real a_rms = 0;
};

inline Doc223E2ViscBand hostDoc223E2ViscBand(BaseParticles &particles, const Vecd &center, Real radius,
                                              Real z_bottom, Real z_lid, Real dp, const char *band)
{
    Doc223E2ViscBand out;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *fv = particles.getVariableDataByName<Vecd>("ViscousForce");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    const Real lid_band = Real(2.5) * dp;
    const Real wall_band = Real(2.5) * dp;
    const Real core_clear = Real(3.0) * dp;
    Real s_tx = 0, s_tz = 0, s_a = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real z = pos[i][2];
        bool take = false;
        if (std::strcmp(band, "lid_core") == 0)
            take = (z > z_lid - lid_band) && (r < radius - core_clear);
        else if (std::strcmp(band, "side") == 0)
            take = (r > radius - wall_band) && (z > z_bottom + wall_band) && (z < z_lid - lid_band);
        else if (std::strcmp(band, "bot") == 0)
            take = (z < z_bottom + wall_band) && (r < radius - wall_band);
        if (!take)
            continue;
        const Real m = std::max(mass[i], TinyReal);
        const Real ax = fv[i][0] / m;
        const Real az = fv[i][2] / m;
        s_tx += ax * ax;
        s_tz += az * az;
        s_a += (fv[i] / m).squaredNorm();
        ++out.n;
    }
    if (out.n > 0)
    {
        const Real inv = Real(1) / Real(out.n);
        out.a_tx_rms = std::sqrt(s_tx * inv);
        out.a_tz_rms = std::sqrt(s_tz * inv);
        out.a_rms = std::sqrt(s_a * inv);
    }
    return out;
}

/** Doc226 B: lid_core band stats for Force / DensityChangeRate / ViscousForce. */
struct Doc226TopOpBand
{
    size_t n = 0;
    Real a_tx_rms = 0; // Fx/m
    Real a_tz_rms = 0; // Fz/m
    Real a_t_rms = 0;  // sqrt(ax^2+ay^2)/m tangential RMS
    Real a_rms = 0;
    Real drho_dt_rms = 0;
    Real work_tx = 0; // Σ F·U (tangential U imposed)
};

inline Doc226TopOpBand hostDoc226LidCoreForceBand(BaseParticles &particles, const Vecd &center,
                                                   Real radius, Real z_bottom, Real z_lid, Real dp,
                                                   const std::string &force_name, Real u_tx_for_work)
{
    Doc226TopOpBand out;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, force_name);
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *fv = particles.getVariableDataByName<Vecd>(force_name);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const size_t n = particles.TotalRealParticles();
    const Real lid_band = Real(2.5) * dp;
    const Real core_clear = Real(3.0) * dp;
    Real s_tx = 0, s_tz = 0, s_t = 0, s_a = 0, s_dr = 0, work = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real z = pos[i][2];
        if (!(z > z_lid - lid_band) || !(r < radius - core_clear))
            continue;
        const Real m = std::max(mass[i], TinyReal);
        const Real ax = fv[i][0] / m;
        const Real ay = fv[i][1] / m;
        const Real az = fv[i][2] / m;
        s_tx += ax * ax;
        s_tz += az * az;
        s_t += ax * ax + ay * ay;
        s_a += (fv[i] / m).squaredNorm();
        s_dr += drho[i] * drho[i];
        work += fv[i][0] * u_tx_for_work;
        ++out.n;
    }
    if (out.n > 0)
    {
        const Real inv = Real(1) / Real(out.n);
        out.a_tx_rms = std::sqrt(s_tx * inv);
        out.a_tz_rms = std::sqrt(s_tz * inv);
        out.a_t_rms = std::sqrt(s_t * inv);
        out.a_rms = std::sqrt(s_a * inv);
        out.drho_dt_rms = std::sqrt(s_dr * inv);
        out.work_tx = work;
    }
    return out;
}

inline void hostDoc223SnapshotVecField(BaseParticles &particles, const std::string &name,
                                       std::vector<Vecd> &buf)
{
    syncVariableToHost<Vecd>(particles, name);
    const Vecd *v = particles.getVariableDataByName<Vecd>(name);
    buf.assign(v, v + particles.TotalRealParticles());
}

inline void hostDoc223RestoreVecField(BaseParticles &particles, const std::string &name,
                                      const std::vector<Vecd> &buf)
{
    syncVariableToHost<Vecd>(particles, name);
    Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = buf[i];
    syncVariableToDevice<Vecd>(particles, name);
}

inline void hostDoc223SetUniformTangentialVelocity(BaseParticles &particles, Real u_tx)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        vel[i] = Vecd(u_tx, Real(0), Real(0));
    syncVariableToDevice<Vecd>(particles, "Velocity");
}

/** Doc223 E3 IC + end metrics. */
struct Doc223E3State
{
    bool active = false;
    Real amp = 0;
    Real ma = 0;
    Real u_rms0 = 0;
    Real uz_mean0 = 0;
    Real ke0 = 0;
    Real z_com0 = 0;
    Real z_max0 = 0;
    size_t n_above0 = 0;
};

inline Real hostDoc223MassComZ(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    Real mz = 0, msum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        mz += mass[i] * pos[i][2];
        msum += mass[i];
    }
    return (msum > TinyReal) ? (mz / msum) : Real(0);
}

inline Real hostDoc223KineticEnergy(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    Real ke = 0;
    for (size_t i = 0; i < n; ++i)
        ke += Real(0.5) * mass[i] * vel[i].squaredNorm();
    return ke;
}

/**
 * Near-lid normal bump Uz=A·exp(-((z_lid-z)/(3dp))²)·(1-(r/R)²), then remove
 * mass-weighted mean velocity (all components). Amp default 1e-3 ⇒ Ma≪1 for c0~13.
 */
inline Doc223E3State hostDoc223ApplyE3NormalPerturb(BaseParticles &particles, const Vecd &center,
                                                    Real radius, Real z_lid, Real dp, Real c0,
                                                    Real amp)
{
    Doc223E3State st;
    st.active = true;
    st.amp = amp;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    const Real inv_r2 = Real(1) / std::max(radius * radius, TinyReal);
    const Real sig = Real(3) * dp;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r2 = dx * dx + dy * dy;
        const Real radial = std::max(Real(0), Real(1) - r2 * inv_r2);
        const Real s = (z_lid - pos[i][2]) / sig;
        const Real uz = amp * std::exp(-s * s) * radial;
        vel[i] = Vecd(Real(0), Real(0), uz);
    }
    Real msum = 0, mx = 0, my = 0, mz = 0;
    for (size_t i = 0; i < n; ++i)
    {
        msum += mass[i];
        mx += mass[i] * vel[i][0];
        my += mass[i] * vel[i][1];
        mz += mass[i] * vel[i][2];
    }
    const Vecd umean = (msum > TinyReal) ? Vecd(mx / msum, my / msum, mz / msum) : Vecd::Zero();
    Real s2 = 0, umax = 0, zmax = -Real(1e30);
    size_t n_above = 0;
    for (size_t i = 0; i < n; ++i)
    {
        vel[i] -= umean;
        s2 += vel[i].squaredNorm();
        umax = std::max(umax, vel[i].norm());
        zmax = std::max(zmax, pos[i][2]);
        if (pos[i][2] > z_lid)
            ++n_above;
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    st.u_rms0 = (n > 0) ? std::sqrt(s2 / Real(n)) : Real(0);
    st.uz_mean0 = Real(0); // removed
    st.ma = (c0 > TinyReal) ? (umax / c0) : Real(0);
    st.ke0 = hostDoc223KineticEnergy(particles);
    st.z_com0 = hostDoc223MassComZ(particles);
    st.z_max0 = zmax;
    st.n_above0 = n_above;
    return st;
}

inline void hostDoc223LogE3LidAudit(const char *tag, BaseParticles &particles, Real z_lid, Real dp,
                                    Real ke0, Real z_com0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    const Real band = Real(2.5) * dp;
    size_t n_lid = 0, n_above = 0;
    Real un2 = 0, un_max = 0, z_max = -Real(1e30), overshoot = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        z_max = std::max(z_max, z);
        if (z > z_lid)
        {
            ++n_above;
            overshoot = std::max(overshoot, z - z_lid);
        }
        if (z > z_lid - band)
        {
            ++n_lid;
            const Real un = vel[i][2];
            un2 += un * un;
            un_max = std::max(un_max, std::abs(un));
        }
    }
    const Real un_rms = (n_lid > 0) ? std::sqrt(un2 / Real(n_lid)) : Real(0);
    const Real ke = hostDoc223KineticEnergy(particles);
    const Real z_com = hostDoc223MassComZ(particles);
    std::cout << "[ophelie][doc223][E3][" << tag << "]"
              << " n_lid_band=" << n_lid << " un_rms=" << un_rms << " un_abs_max=" << un_max
              << " n_above_lid=" << n_above << " overshoot_max=" << overshoot
              << " z_max=" << z_max << " z_lid=" << z_lid
              << " ke=" << ke << " ke/ke0=" << ((ke0 > TinyReal) ? (ke / ke0) : Real(0))
              << " dz_com=" << (z_com - z_com0) << " support~" << (Real(2) * dp) << "\n";
}

/** Doc221: ForcePrior ledger vs registered non-pressure channels (Gravity+Boussinesq+Viscous). */
struct Doc221ForceLedgerStats
{
    size_t n = 0;
    Real fp_rms = 0;
    Real grav_rms = 0;
    Real buoy_rms = 0;
    Real visc_rms = 0;
    Real sum_rms = 0;
    Real residual_rms = 0;
    Real residual_max = 0;
    Real residual_rel_rms = 0;
    Real fp_z_mean = 0;
    Real buoy_z_mean = 0;
    Real visc_z_mean = 0;
};

inline Doc221ForceLedgerStats hostDoc221ForceLedgerStats(BaseParticles &particles)
{
    Doc221ForceLedgerStats st;
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "GravityForceCK");
    syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *fg = particles.getVariableDataByName<Vecd>("GravityForceCK");
    const Vecd *fb = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
    const Vecd *fv = particles.getVariableDataByName<Vecd>("ViscousForce");
    st.n = particles.TotalRealParticles();
    Real s_fp = 0, s_g = 0, s_b = 0, s_v = 0, s_sum = 0, s_r = 0, s_scale = 0;
    Real s_fpz = 0, s_bz = 0, s_vz = 0;
    for (size_t i = 0; i < st.n; ++i)
    {
        const Vecd sum = fg[i] + fb[i] + fv[i];
        const Vecd r = fp[i] - sum;
        const Real nfp = fp[i].squaredNorm();
        const Real ng = fg[i].squaredNorm();
        const Real nb = fb[i].squaredNorm();
        const Real nv = fv[i].squaredNorm();
        const Real ns = sum.squaredNorm();
        const Real nr = r.squaredNorm();
        s_fp += nfp;
        s_g += ng;
        s_b += nb;
        s_v += nv;
        s_sum += ns;
        s_r += nr;
        st.residual_max = std::max(st.residual_max, std::sqrt(nr));
        const Real scale = std::sqrt(ng + nb + nv) + Real(1.0e-12);
        s_scale += (nr / (scale * scale));
        s_fpz += fp[i][2];
        s_bz += fb[i][2];
        s_vz += fv[i][2];
    }
    const Real inv = st.n > 0 ? Real(1) / Real(st.n) : Real(0);
    st.fp_rms = std::sqrt(s_fp * inv);
    st.grav_rms = std::sqrt(s_g * inv);
    st.buoy_rms = std::sqrt(s_b * inv);
    st.visc_rms = std::sqrt(s_v * inv);
    st.sum_rms = std::sqrt(s_sum * inv);
    st.residual_rms = std::sqrt(s_r * inv);
    st.residual_rel_rms = std::sqrt(s_scale * inv);
    st.fp_z_mean = s_fpz * inv;
    st.buoy_z_mean = s_bz * inv;
    st.visc_z_mean = s_vz * inv;
    return st;
}

inline void hostDoc221AppendForceLedgerCsv(const std::string &path, const std::string &tag, int adv_step,
                                           Real physical_time, const Doc221ForceLedgerStats &st, bool write_header)
{
    std::ofstream out(path, write_header ? std::ios::out : std::ios::app);
    if (!out.good())
        return;
    if (write_header)
    {
        out << "tag,adv_step,physical_time_s,n,fp_rms,grav_rms,buoy_rms,visc_rms,sum_rms,"
               "residual_rms,residual_max,residual_rel_rms,fp_z_mean,buoy_z_mean,visc_z_mean\n";
    }
    out << tag << "," << adv_step << "," << physical_time << "," << st.n << "," << st.fp_rms << ","
        << st.grav_rms << "," << st.buoy_rms << "," << st.visc_rms << "," << st.sum_rms << ","
        << st.residual_rms << "," << st.residual_max << "," << st.residual_rel_rms << ","
        << st.fp_z_mean << "," << st.buoy_z_mean << "," << st.visc_z_mean << "\n";
}

inline void logDoc221ForceLedger(const std::string &tag, const Doc221ForceLedgerStats &st)
{
    std::cout << "[ophelie][doc221][ledger] " << tag << " n=" << st.n << " fp_rms=" << st.fp_rms
              << " grav_rms=" << st.grav_rms << " buoy_rms=" << st.buoy_rms << " visc_rms=" << st.visc_rms
              << " sum_rms=" << st.sum_rms << " residual_rms=" << st.residual_rms
              << " residual_max=" << st.residual_max << " residual_rel_rms=" << st.residual_rel_rms
              << " fp_z_mean=" << st.fp_z_mean << " buoy_z_mean=" << st.buoy_z_mean
              << " visc_z_mean=" << st.visc_z_mean << std::endl;
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

/**
 * Doc201 B2 / Doc221 P: tag flat-lid wall particles AcousticWallFreeSlip=1 so acoustic dens ghost
 * is free-slip and 1st-half mom force projects onto n. Side/bottom stay no-slip. Requires
 * AcousticWallFreeSlip already registered (AcousticStep Contact<Wall> ctor) and normals set.
 */
inline size_t hostTagLidWallRealFieldFrench(BaseParticles &wall_particles, Real z_top, Real wall_thickness,
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
    size_t n_lid = 0;
    for (size_t i = 0; i < nn; ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        f[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_lid;
    }
    syncVariableToDevice<Real>(wall_particles, var_name);
    return n_lid;
}

inline size_t hostTagAcousticWallFreeSlipLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealFieldFrench(wall_particles, z_top, wall_thickness, "AcousticWallFreeSlip");
}

/** Doc224 C3: keep dens jump on free-slip lid faces (mom still ∥ n). */
inline size_t hostTagAcousticWallKeepDensJumpLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealFieldFrench(wall_particles, z_top, wall_thickness, "AcousticWallKeepDensJump");
}

inline size_t hostArmDoc224FreeSlipLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness,
                                       bool keep_dens_jump)
{
    const size_t n_lid = hostTagAcousticWallFreeSlipLid(wall_particles, z_top, wall_thickness);
    if (keep_dens_jump)
    {
        const size_t n_keep = hostTagAcousticWallKeepDensJumpLid(wall_particles, z_top, wall_thickness);
        std::cout << "[ophelie][doc224][C3] AcousticWallKeepDensJump lid tags=" << n_keep
                  << " (free-slip mom∥n + dens jump ON)\n";
    }
    return n_lid;
}

/** Doc202 H2a / Doc221 O4/O5: lid skip 2nd-half continuity dens (keep PJump). */
inline size_t hostTagAcousticWallContSilentLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealFieldFrench(wall_particles, z_top, wall_thickness, "AcousticWallContSilent");
}

/** Doc202 H2b / Doc221 O4/O5: lid omit DissipativePJump (keep continuity). */
inline size_t hostTagAcousticWallOmitPJumpLid(BaseParticles &wall_particles, Real z_top, Real wall_thickness)
{
    return hostTagLidWallRealFieldFrench(wall_particles, z_top, wall_thickness, "AcousticWallOmitPJump");
}

/** Doc220 E2: PerturbationMomentum 1 ⇒ acoustic momentum uses p' = p − BackgroundPressure. */
inline void hostSetPerturbationMomentumFrench(BaseParticles &particles, Real flag)
{
    syncVariableToHost<Real>(particles, "PerturbationMomentum");
    Real *pm = particles.getVariableDataByName<Real>("PerturbationMomentum");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        pm[i] = flag;
    syncVariableToDevice<Real>(particles, "PerturbationMomentum");
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

inline void hostLogDoc223PaperTopHydroAudit(const char *tag, BaseParticles &particles, Real rho0, Real c0,
                                            Real gravity_g, const Vecd &center, Real radius, Real z_bottom,
                                            Real z_lid, Real dp)
{
    const Real band = Real(2) * dp;
    const FrenchHydroStageSnapshot st =
        hostFrenchHydroStageSnapshot(particles, rho0, c0, gravity_g, center, radius, z_bottom, z_lid, band,
                                     Real(4) * dp);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    size_t n_above = 0;
    size_t n_lid_band = 0;
    Real z_max = -Real(1e30);
    Real gap_mean = Real(0);
    Real un_rms = Real(0);
    Real un_abs_max = Real(0);
    Real overshoot_max = Real(0);
    Real mass_tot = Real(0);
    Real mass_above = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        const Real m = mass[i];
        mass_tot += m;
        z_max = SMAX(z_max, z);
        if (z > z_lid)
        {
            ++n_above;
            mass_above += m;
            overshoot_max = SMAX(overshoot_max, z - z_lid);
        }
        if (z > z_lid - band)
        {
            ++n_lid_band;
            const Real un = vel[i][2];
            un_rms += un * un;
            un_abs_max = SMAX(un_abs_max, std::abs(un));
            gap_mean += z_lid - z;
        }
    }
    if (n_lid_band > 0)
    {
        un_rms = std::sqrt(un_rms / static_cast<Real>(n_lid_band));
        gap_mean /= static_cast<Real>(n_lid_band);
    }
    const Real f_over = (mass_tot > TinyReal) ? (mass_above / mass_tot) : Real(0);
    const Real overshoot_dp = (dp > TinyReal) ? (overshoot_max / dp) : Real(0);
    std::cout << "[ophelie][doc223][E1][" << tag << "]"
              << " n_int=" << st.interior.n << " n_surf=" << st.surface.n
              << " n_side=" << st.wall_side.n << " n_bot=" << st.wall_bottom.n
              << " u_int=" << st.interior.u_rms << " u_surf=" << st.surface.u_rms
              << " u_side=" << st.wall_side.u_rms << " u_bot=" << st.wall_bottom.u_rms
              << " u_side_inc=" << st.wall_side_incomplete.u_rms
              << " u_side_comp=" << st.wall_side_complete.u_rms
              << " a_f_side=" << st.wall_side.a_force_z_rms
              << " a_f_r_side=" << st.wall_side.a_force_r_rms
              << " pack_int=" << st.interior.packing_mean << " pack_surf=" << st.surface.packing_mean
              << " pack_side=" << st.wall_side.packing_mean << " pack_bot=" << st.wall_bottom.packing_mean
              << " p_eos_int=" << st.interior.p_eos_rms << " p_eos_surf=" << st.surface.p_eos_rms
              << " a_net_int=" << st.interior.a_net_z_rms << " a_net_surf=" << st.surface.a_net_z_rms
              << " a_f_int=" << st.interior.a_force_z_rms << " a_f_surf=" << st.surface.a_force_z_rms
              << " n_inc_surf=" << st.surface.n_incomplete << " n_inc_side=" << st.wall_side.n_incomplete
              << " z_max=" << z_max << " z_lid=" << z_lid << " n_above_lid=" << n_above
              << " f_over=" << f_over << " overshoot_max=" << overshoot_max
              << " overshoot_dp=" << overshoot_dp
              << " n_lid_band=" << n_lid_band << " gap_mean=" << gap_mean << " un_rms=" << un_rms
              << " un_abs_max=" << un_abs_max << "\n";
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

/**
 * Doc221 K: subtract mass-weighted mean Uz so ∫ Uz dV ≈ 0 before position update.
 * Returns the subtracted mean (diagnostic log). Opt-in only.
 */
inline Real hostSubtractMassWeightedMeanUz(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    Real m_sum = 0;
    Real muz = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real m = mass[i];
        m_sum += m;
        muz += m * vel[i][2];
    }
    const Real uz_mean = m_sum > TinyReal ? muz / m_sum : Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        vel[i][2] -= uz_mean;
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return uz_mean;
}

/**
 * Doc221 L: kinematic hard lid — in band z > z_lid - band, set Uz = min(Uz, 0).
 * Returns number of particles with upward velocity killed. Opt-in only.
 */
inline size_t hostKillUpwardVelocityNearLid(BaseParticles &particles, Real z_lid, Real band)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_lo = z_lid - band;
    size_t n_kill = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] > z_lo && vel[i][2] > Real(0))
        {
            vel[i][2] = Real(0);
            ++n_kill;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return n_kill;
}

/**
 * Doc221 N4/N5 / Doc222 B: tip soft physical lid via ForcePrior (no host Uz edit, no z clamp).
 * In tip cylinder + surface band: a_z = -k*max(z-z_lid,0) - c*max(Uz,0).
 * Doc222 default: ForcePriorCK-style channel
 *   TipSoftForce = m*a_tip; ForcePrior += TipSoftForce - PreviousTipSoftForce; Previous=TipSoftForce.
 * All particles updated (outside tip / Uz<=0 write zero and revoke). No extra *dt.
 * legacy_accumulate=true restores Doc221 direct ForcePrior += m*a (B0 repro only).
 * If fs_incomplete_only: also require DensitySummation < rho0 (N5 FS-path bind).
 * Returns number of particles with non-zero tip acceleration this call.
 */
inline size_t hostApplyTipSoftForcePrior(BaseParticles &particles, Real z_lid, Real band, Real tip_r,
                                          Real k_spring, Real c_damp, bool fs_incomplete_only = false,
                                          Real rho0 = Real(0), bool legacy_accumulate = false)
{
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, "TipSoftForce");
    syncVariableToHost<Vecd>(particles, "PreviousTipSoftForce");
    const Real *rho_sum = nullptr;
    if (fs_incomplete_only)
    {
        syncVariableToHost<Real>(particles, "DensitySummation");
        rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    }
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Vecd *tip_cur = particles.getVariableDataByName<Vecd>("TipSoftForce");
    Vecd *tip_prev = particles.getVariableDataByName<Vecd>("PreviousTipSoftForce");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    const Real z_lo = z_lid - band;
    const Real tip_r2 = tip_r * tip_r;
    size_t n_touch = 0;
    for (size_t i = 0; i < n; ++i)
    {
        Real a_z = Real(0);
        const bool in_band = pos[i][2] > z_lo;
        const Real r2 = pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1];
        const bool in_tip = in_band && r2 < tip_r2;
        const bool fs_ok =
            !fs_incomplete_only || rho_sum == nullptr || (rho_sum[i] + TinyReal < rho0);
        if (in_tip && fs_ok)
        {
            const Real pen = pos[i][2] - z_lid;
            const Real uz_up = vel[i][2] > Real(0) ? vel[i][2] : Real(0);
            a_z = -k_spring * (pen > Real(0) ? pen : Real(0)) - c_damp * uz_up;
        }
        Vecd f_new = Vecd::Zero();
        f_new[2] = mass[i] * a_z;
        if (std::abs(a_z) > TinyReal)
            ++n_touch;
        if (legacy_accumulate)
        {
            // Doc221 bug path: only add when nonzero; never revoke previous contribution.
            if (std::abs(a_z) > TinyReal)
                fp[i][2] += f_new[2];
            tip_cur[i] = f_new;
            // leave tip_prev untouched in legacy (not a real channel)
        }
        else
        {
            fp[i] += f_new - tip_prev[i];
            tip_prev[i] = f_new;
            tip_cur[i] = f_new;
        }
    }
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    syncVariableToDevice<Vecd>(particles, "TipSoftForce");
    syncVariableToDevice<Vecd>(particles, "PreviousTipSoftForce");
    return n_touch;
}

/** Doc222 B: read-only tip touch count — never mutates ForcePrior / tip channels. */
inline size_t hostCountTipSoftForceTouches(BaseParticles &particles, Real z_lid, Real band, Real tip_r,
                                            Real k_spring, Real c_damp, bool fs_incomplete_only = false,
                                            Real rho0 = Real(0))
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *rho_sum = nullptr;
    if (fs_incomplete_only)
    {
        syncVariableToHost<Real>(particles, "DensitySummation");
        rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    }
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_lo = z_lid - band;
    const Real tip_r2 = tip_r * tip_r;
    size_t n_touch = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] <= z_lo)
            continue;
        const Real r2 = pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1];
        if (r2 >= tip_r2)
            continue;
        if (fs_incomplete_only && rho_sum != nullptr && rho_sum[i] + TinyReal >= rho0)
            continue;
        const Real pen = pos[i][2] - z_lid;
        const Real uz_up = vel[i][2] > Real(0) ? vel[i][2] : Real(0);
        const Real a_z = -k_spring * (pen > Real(0) ? pen : Real(0)) - c_damp * uz_up;
        if (std::abs(a_z) > TinyReal)
            ++n_touch;
    }
    return n_touch;
}

/** Doc222 B: revoke tip channel from ForcePrior and zero tip fields (disable / selftest reset). */
inline void hostClearTipSoftForcePriorContribution(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "TipSoftForce");
    syncVariableToHost<Vecd>(particles, "PreviousTipSoftForce");
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Vecd *tip_cur = particles.getVariableDataByName<Vecd>("TipSoftForce");
    Vecd *tip_prev = particles.getVariableDataByName<Vecd>("PreviousTipSoftForce");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        fp[i] -= tip_prev[i];
        tip_prev[i] = Vecd::Zero();
        tip_cur[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    syncVariableToDevice<Vecd>(particles, "TipSoftForce");
    syncVariableToDevice<Vecd>(particles, "PreviousTipSoftForce");
}

struct Doc222TipForceSelftestResult
{
    bool pass = false;
    bool repeat10_ok = false;
    bool exit_zone_ok = false;
    bool uz_down_ok = false;
    bool clear_ok = false;
    bool count_no_side_effect = false;
    bool legacy_grows = false;
    Real fp_z_after_1 = 0;
    Real fp_z_after_10 = 0;
    Real fp_z_after_exit = 0;
    Real fp_z_after_clear = 0;
    Real fp_z_legacy_10 = 0;
    size_t n_touch_1 = 0;
    std::string message;
};

inline Real hostSumForcePriorZ(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const size_t n = particles.TotalRealParticles();
    Real s = 0;
    for (size_t i = 0; i < n; ++i)
        s += fp[i][2];
    return s;
}

/**
 * Doc222 B: fixed-state tip ForcePrior lifecycle checks (host). Mutates a few tip particles then restores.
 */
inline Doc222TipForceSelftestResult hostDoc222TipForceLifecycleSelftest(BaseParticles &particles, Real z_lid,
                                                                        Real band, Real tip_r, Real c_damp)
{
    Doc222TipForceSelftestResult r;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        r.message = "no_particles";
        return r;
    }

    // Pick / synthesize one tip-surface particle with Uz>0.
    size_t i_tip = n;
    const Real tip_r2 = tip_r * tip_r;
    for (size_t i = 0; i < n; ++i)
    {
        const Real r2 = pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1];
        if (r2 < tip_r2 && pos[i][2] > z_lid - band)
        {
            i_tip = i;
            break;
        }
    }
    if (i_tip >= n)
        i_tip = 0;

    const Vecd pos0 = pos[i_tip];
    const Vecd vel0 = vel[i_tip];
    std::vector<Vecd> fp0(n);
    for (size_t i = 0; i < n; ++i)
        fp0[i] = fp[i];

    pos[i_tip] = Vecd(Real(0), Real(0), z_lid + Real(0.25) * band);
    vel[i_tip] = Vecd(Real(0), Real(0), Real(0.02));
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");

    hostClearTipSoftForcePriorContribution(particles);
    // Baseline ForcePrior after clear (tip channel zero).
    const Real fp_base = hostSumForcePriorZ(particles);

    r.n_touch_1 = hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0),
                                             false);
    r.fp_z_after_1 = hostSumForcePriorZ(particles) - fp_base;
    for (int k = 0; k < 9; ++k)
        hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0), false);
    r.fp_z_after_10 = hostSumForcePriorZ(particles) - fp_base;
    const Real tol = Real(1.0e-9) * (std::abs(r.fp_z_after_1) + Real(1));
    r.repeat10_ok = std::abs(r.fp_z_after_10 - r.fp_z_after_1) <= tol && std::abs(r.fp_z_after_1) > TinyReal;

    // Count must not change ForcePrior.
    const Real before_count = hostSumForcePriorZ(particles);
    const size_t n_cnt =
        hostCountTipSoftForceTouches(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0));
    const Real after_count = hostSumForcePriorZ(particles);
    r.count_no_side_effect =
        (n_cnt == r.n_touch_1) && (std::abs(after_count - before_count) <= Real(1.0e-12) * (std::abs(before_count) + Real(1)));

    // Exit tip cylinder → revoke.
    pos[i_tip][0] = tip_r + band + Real(0.05);
    syncVariableToDevice<Vecd>(particles, "Position");
    hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0), false);
    r.fp_z_after_exit = hostSumForcePriorZ(particles) - fp_base;
    r.exit_zone_ok = std::abs(r.fp_z_after_exit) <= tol;

    // Restore tip, set Uz downward → contribution zero.
    pos[i_tip] = Vecd(Real(0), Real(0), z_lid + Real(0.25) * band);
    vel[i_tip] = Vecd(Real(0), Real(0), Real(-0.02));
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0), false);
    const Real fp_down = hostSumForcePriorZ(particles) - fp_base;
    r.uz_down_ok = std::abs(fp_down) <= tol;

    hostClearTipSoftForcePriorContribution(particles);
    r.fp_z_after_clear = hostSumForcePriorZ(particles) - fp_base;
    r.clear_ok = std::abs(r.fp_z_after_clear) <= tol;

    // Legacy accumulate grows ~10x on fixed state (bug repro).
    hostClearTipSoftForcePriorContribution(particles);
    for (size_t i = 0; i < n; ++i)
        fp[i] = fp0[i];
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    pos[i_tip] = Vecd(Real(0), Real(0), z_lid + Real(0.25) * band);
    vel[i_tip] = Vecd(Real(0), Real(0), Real(0.02));
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    const Real fp_leg_base = hostSumForcePriorZ(particles);
    hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0), true);
    const Real fp_leg_1 = hostSumForcePriorZ(particles) - fp_leg_base;
    for (int k = 0; k < 9; ++k)
        hostApplyTipSoftForcePrior(particles, z_lid, band, tip_r, Real(0), c_damp, false, Real(0), true);
    r.fp_z_legacy_10 = hostSumForcePriorZ(particles) - fp_leg_base;
    r.legacy_grows = std::abs(r.fp_z_legacy_10) > Real(5) * std::abs(fp_leg_1) && std::abs(fp_leg_1) > TinyReal;

    // Restore particle state and zero tip channel.
    pos[i_tip] = pos0;
    vel[i_tip] = vel0;
    for (size_t i = 0; i < n; ++i)
        fp[i] = fp0[i];
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    hostClearTipSoftForcePriorContribution(particles);

    r.pass = r.repeat10_ok && r.exit_zone_ok && r.uz_down_ok && r.clear_ok && r.count_no_side_effect &&
             r.legacy_grows;
    std::ostringstream oss;
    oss << "repeat10=" << (r.repeat10_ok ? 1 : 0) << " exit=" << (r.exit_zone_ok ? 1 : 0)
        << " uz_down=" << (r.uz_down_ok ? 1 : 0) << " clear=" << (r.clear_ok ? 1 : 0)
        << " count_safe=" << (r.count_no_side_effect ? 1 : 0) << " legacy_grows=" << (r.legacy_grows ? 1 : 0)
        << " fp1=" << r.fp_z_after_1 << " fp10=" << r.fp_z_after_10 << " leg10=" << r.fp_z_legacy_10
        << " n_touch=" << r.n_touch_1;
    r.message = oss.str();
    return r;
}

/** Doc222 D2: η(r)=eps*[1-2(r/R)^2]; ∫_0^R 2π r η dr = 0. */
inline Real doc222FsEtaOfR(Real r, Real R, Real eps)
{
    if (R <= TinyReal)
        return Real(0);
    const Real rr = r / R;
    return eps * (Real(1) - Real(2) * rr * rr);
}

struct Doc222FsDeformReport
{
    size_t n = 0;
    Real eps = 0;
    Real eps_over_dp = 0;
    Real eps_over_H = 0;
    Real vol_sum_before = 0;
    Real vol_sum_after = 0;
    Real mass_over_rho_before = 0;
    Real mass_over_rho_after = 0;
    Real z_max_before = 0;
    Real z_max_after = 0;
    Real eta_center = 0;
    Real eta_rim = 0;
};

inline Real hostSumParticleVolumeFrench(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real s = 0;
    for (size_t i = 0; i < n; ++i)
        s += Vol[i];
    return s;
}

inline Real hostSumMassOverDensityFrench(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    Real s = 0;
    for (size_t i = 0; i < n; ++i)
        s += mass[i] / (rho[i] + TinyReal);
    return s;
}

/**
 * Doc222 D2: column-stretch free-surface deformation about fixed (R,H) cylinder.
 * z' = z_bot + (z-z_bot)/H * (H+η(r)), η=eps*[1-2(r/R)^2]. Keeps particle mass; caller
 * must re-init density/pressure and rebuild neighbor lists.
 */
inline Doc222FsDeformReport hostApplyDoc222AxisymmetricFsDeform(BaseParticles &particles, const Vecd &center,
                                                                Real R, Real z_bot, Real z_free, Real eps)
{
    Doc222FsDeformReport rep;
    rep.eps = eps;
    rep.eta_center = doc222FsEtaOfR(Real(0), R, eps);
    rep.eta_rim = doc222FsEtaOfR(R, R, eps);
    const Real H = z_free - z_bot;
    rep.eps_over_H = (H > TinyReal) ? eps / H : Real(0);
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    rep.n = particles.TotalRealParticles();
    rep.vol_sum_before = hostSumParticleVolumeFrench(particles);
    rep.mass_over_rho_before = hostSumMassOverDensityFrench(particles);
    rep.z_max_before = -Real(1.0e30);
    for (size_t i = 0; i < rep.n; ++i)
        rep.z_max_before = std::max(rep.z_max_before, pos[i][2]);
    if (H <= TinyReal || std::abs(eps) <= TinyReal)
        return rep;
    for (size_t i = 0; i < rep.n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real eta = doc222FsEtaOfR(std::min(r, R), R, eps);
        Real Hloc = H + eta;
        if (Hloc < Real(0.25) * H)
            Hloc = Real(0.25) * H;
        Real z_rel = (pos[i][2] - z_bot) / H;
        if (z_rel < Real(0))
            z_rel = Real(0);
        pos[i][2] = z_bot + z_rel * Hloc;
    }
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Position");
    pos = particles.getVariableDataByName<Vecd>("Position");
    rep.z_max_after = -Real(1.0e30);
    for (size_t i = 0; i < rep.n; ++i)
        rep.z_max_after = std::max(rep.z_max_after, pos[i][2]);
    rep.vol_sum_after = hostSumParticleVolumeFrench(particles);
    rep.mass_over_rho_after = hostSumMassOverDensityFrench(particles);
    return rep;
}

struct Doc222FsSurfaceStats
{
    Real mode_amp = 0;     // projection onto φ=1-2(r/R)^2
    Real z_max = 0;
    Real z_p10 = 0;
    Real z_p50 = 0;
    Real z_p90 = 0;
    Real z_mean_bins = 0;
    Real ke = 0;
    Real u_rms = 0;
    Real mass_over_rho = 0;
    int n_bins = 0;
    int n_bins_valid = 0;
};

inline Doc222FsSurfaceStats hostMeasureDoc222FsSurfaceStats(BaseParticles &particles, const Vecd &center, Real R,
                                                            Real z_ref, int n_bins = 12)
{
    Doc222FsSurfaceStats st;
    st.n_bins = n_bins;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    std::vector<Real> bin_zmax(static_cast<size_t>(n_bins), -Real(1.0e30));
    std::vector<int> bin_count(static_cast<size_t>(n_bins), 0);
    Real ke = 0, msum = 0;
    Real u2 = 0;
    st.z_max = -Real(1.0e30);
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        st.z_max = std::max(st.z_max, pos[i][2]);
        ke += Real(0.5) * mass[i] * vel[i].squaredNorm();
        msum += mass[i];
        u2 += mass[i] * vel[i].squaredNorm();
        if (r >= R || n_bins <= 0)
            continue;
        int b = static_cast<int>(std::floor(r / R * Real(n_bins)));
        if (b < 0)
            b = 0;
        if (b >= n_bins)
            b = n_bins - 1;
        bin_zmax[static_cast<size_t>(b)] = std::max(bin_zmax[static_cast<size_t>(b)], pos[i][2]);
        bin_count[static_cast<size_t>(b)]++;
    }
    st.ke = ke;
    st.u_rms = (msum > TinyReal) ? std::sqrt(u2 / msum) : Real(0);
    st.mass_over_rho = hostSumMassOverDensityFrench(particles);

    Real num = 0, den = 0, zsum = 0;
    std::vector<Real> valid_z;
    valid_z.reserve(static_cast<size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b)
    {
        if (bin_count[static_cast<size_t>(b)] <= 0 || bin_zmax[static_cast<size_t>(b)] < -Real(1.0e20))
            continue;
        const Real r0 = (Real(b) + Real(0.5)) / Real(n_bins) * R;
        const Real phi = Real(1) - Real(2) * (r0 / R) * (r0 / R);
        const Real w = r0; // ∝ annulus weight (dr const)
        const Real eta = bin_zmax[static_cast<size_t>(b)] - z_ref;
        num += w * eta * phi;
        den += w * phi * phi;
        zsum += bin_zmax[static_cast<size_t>(b)];
        valid_z.push_back(bin_zmax[static_cast<size_t>(b)]);
        st.n_bins_valid++;
    }
    st.mode_amp = (den > TinyReal) ? (num / den) : Real(0);
    st.z_mean_bins = (st.n_bins_valid > 0) ? (zsum / Real(st.n_bins_valid)) : Real(0);
    if (!valid_z.empty())
    {
        std::sort(valid_z.begin(), valid_z.end());
        auto pct = [&](Real p) {
            const Real idx = p * Real(valid_z.size() - 1);
            const size_t i0 = static_cast<size_t>(std::floor(idx));
            const size_t i1 = std::min(i0 + 1, valid_z.size() - 1);
            const Real a = idx - Real(i0);
            return (Real(1) - a) * valid_z[i0] + a * valid_z[i1];
        };
        st.z_p10 = pct(Real(0.10));
        st.z_p50 = pct(Real(0.50));
        st.z_p90 = pct(Real(0.90));
    }
    return st;
}

inline void hostDoc222AppendFsSurfaceCsv(const std::string &path, Real t, int adv, const Doc222FsSurfaceStats &st,
                                         bool write_header)
{
    std::ofstream out(path, write_header ? std::ios::out : std::ios::app);
    if (!out.good())
        return;
    if (write_header)
    {
        out << "time_s,adv_step,mode_amp,z_max,z_p10,z_p50,z_p90,z_mean_bins,U_rms,KE,mass_over_rho,"
               "n_bins_valid\n";
    }
    out << t << "," << adv << "," << st.mode_amp << "," << st.z_max << "," << st.z_p10 << "," << st.z_p50 << ","
        << st.z_p90 << "," << st.z_mean_bins << "," << st.u_rms << "," << st.ke << "," << st.mass_over_rho << ","
        << st.n_bins_valid << "\n";
}

/** Doc222 C3: pressure/acceleration budget on fixed deformed geometry (one acoustic-1st snapshot). */
struct Doc222C3PhProbeStats
{
    bool signed_depth = false;
    size_t n_all = 0;
    size_t n_above = 0; // z > z_free
    size_t n_surface = 0;
    Real p_bg_mean_all = 0;
    Real p_bg_mean_above = 0;
    Real p_bg_mean_below = 0;
    Real p_bg_min_above = 0;
    Real p_bg_max_above = 0;
    Real p_eos_mean_above = 0;
    Real p_mom_mean_above = 0; // Pressure - BackgroundPressure (α=1)
    Real a_force_z_mean_all = 0;
    Real a_force_z_mean_above = 0;
    Real a_force_z_rms_all = 0;
    Real a_force_z_rms_above = 0;
    Real a_prior_z_mean_all = 0;
    Real a_net_z_mean_all = 0;
    Real a_net_z_mean_above = 0;
    Real a_net_z_rms_above = 0;
    Real z_max = 0;
};

inline Doc222C3PhProbeStats hostMeasureDoc222C3PhProbe(BaseParticles &particles, Real rho0, Real c0, Real z_free,
                                                       const Vecd &center, Real radius, Real z_bottom, Real band,
                                                       bool signed_depth_flag)
{
    Doc222C3PhProbeStats st;
    st.signed_depth = signed_depth_flag;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Pressure");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *p = particles.getVariableDataByName<Real>("Pressure");
    const Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    st.n_all = n;
    st.p_bg_min_above = Real(1.0e30);
    st.p_bg_max_above = -Real(1.0e30);
    st.z_max = -Real(1.0e30);
    Real s_pbg = 0, s_pbg_a = 0, s_pbg_b = 0, s_peos_a = 0, s_pmom_a = 0;
    Real s_afz = 0, s_afz_a = 0, s_afz2 = 0, s_afz2_a = 0;
    Real s_ap = 0, s_an = 0, s_an_a = 0, s_an2_a = 0;
    size_t n_below = 0;
    const Real c2 = c0 * c0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real m = mass[i] + TinyReal;
        const Real afz = f[i][2] / m;
        const Real apz = fp[i][2] / m;
        const Real anz = apz + afz;
        const Real p_eos = c2 * (rho[i] - rho0);
        // Prefer stored Pressure if present; else EOS proxy.
        const Real p_use = p[i];
        const Real p_mom = p_use - p_bg[i];
        st.z_max = std::max(st.z_max, pos[i][2]);
        s_pbg += p_bg[i];
        s_afz += afz;
        s_afz2 += afz * afz;
        s_ap += apz;
        s_an += anz;
        const bool above = pos[i][2] > z_free + TinyReal;
        const FrenchHydroZone zone =
            hostClassifyFrenchHydroZone(pos[i], center, radius, z_bottom, z_free, band);
        if (zone == FrenchHydroZone::Surface || zone == FrenchHydroZone::WallSurface)
            ++st.n_surface;
        if (above)
        {
            ++st.n_above;
            s_pbg_a += p_bg[i];
            s_peos_a += p_eos;
            s_pmom_a += p_mom;
            s_afz_a += afz;
            s_afz2_a += afz * afz;
            s_an_a += anz;
            s_an2_a += anz * anz;
            st.p_bg_min_above = std::min(st.p_bg_min_above, p_bg[i]);
            st.p_bg_max_above = std::max(st.p_bg_max_above, p_bg[i]);
        }
        else
        {
            ++n_below;
            s_pbg_b += p_bg[i];
        }
    }
    const Real inv = n > 0 ? Real(1) / Real(n) : Real(0);
    const Real inv_a = st.n_above > 0 ? Real(1) / Real(st.n_above) : Real(0);
    const Real inv_b = n_below > 0 ? Real(1) / Real(n_below) : Real(0);
    st.p_bg_mean_all = s_pbg * inv;
    st.p_bg_mean_above = s_pbg_a * inv_a;
    st.p_bg_mean_below = s_pbg_b * inv_b;
    st.p_eos_mean_above = s_peos_a * inv_a;
    st.p_mom_mean_above = s_pmom_a * inv_a;
    st.a_force_z_mean_all = s_afz * inv;
    st.a_force_z_mean_above = s_afz_a * inv_a;
    st.a_force_z_rms_all = std::sqrt(s_afz2 * inv);
    st.a_force_z_rms_above = std::sqrt(s_afz2_a * inv_a);
    st.a_prior_z_mean_all = s_ap * inv;
    st.a_net_z_mean_all = s_an * inv;
    st.a_net_z_mean_above = s_an_a * inv_a;
    st.a_net_z_rms_above = std::sqrt(s_an2_a * inv_a);
    if (st.n_above == 0)
    {
        st.p_bg_min_above = 0;
        st.p_bg_max_above = 0;
    }
    return st;
}

inline void hostDoc222WriteC3PhCompareJson(const std::string &path, const Doc222C3PhProbeStats &signed_st,
                                           const Doc222C3PhProbeStats &clamp_st, Real eps, Real z_free)
{
    auto dump = [&](std::ostream &js, const char *key, const Doc222C3PhProbeStats &s) {
        js << "  \"" << key << "\": {\n"
           << "    \"signed_depth\": " << (s.signed_depth ? 1 : 0) << ",\n"
           << "    \"n_all\": " << s.n_all << ",\n"
           << "    \"n_above\": " << s.n_above << ",\n"
           << "    \"n_surface\": " << s.n_surface << ",\n"
           << "    \"z_max\": " << s.z_max << ",\n"
           << "    \"p_bg_mean_all\": " << s.p_bg_mean_all << ",\n"
           << "    \"p_bg_mean_above\": " << s.p_bg_mean_above << ",\n"
           << "    \"p_bg_mean_below\": " << s.p_bg_mean_below << ",\n"
           << "    \"p_bg_min_above\": " << s.p_bg_min_above << ",\n"
           << "    \"p_bg_max_above\": " << s.p_bg_max_above << ",\n"
           << "    \"p_eos_mean_above\": " << s.p_eos_mean_above << ",\n"
           << "    \"p_mom_mean_above\": " << s.p_mom_mean_above << ",\n"
           << "    \"a_force_z_mean_all\": " << s.a_force_z_mean_all << ",\n"
           << "    \"a_force_z_mean_above\": " << s.a_force_z_mean_above << ",\n"
           << "    \"a_force_z_rms_all\": " << s.a_force_z_rms_all << ",\n"
           << "    \"a_force_z_rms_above\": " << s.a_force_z_rms_above << ",\n"
           << "    \"a_prior_z_mean_all\": " << s.a_prior_z_mean_all << ",\n"
           << "    \"a_net_z_mean_all\": " << s.a_net_z_mean_all << ",\n"
           << "    \"a_net_z_mean_above\": " << s.a_net_z_mean_above << ",\n"
           << "    \"a_net_z_rms_above\": " << s.a_net_z_rms_above << "\n"
           << "  }";
    };
    std::ofstream js(path);
    js << "{\n"
       << "  \"case\": \"doc222_c3_ph_compare\",\n"
       << "  \"eps\": " << eps << ",\n"
       << "  \"z_free\": " << z_free << ",\n";
    dump(js, "signed", signed_st);
    js << ",\n";
    dump(js, "clamp", clamp_st);
    js << ",\n"
       << "  \"delta_signed_minus_clamp\": {\n"
       << "    \"p_bg_mean_above\": " << (signed_st.p_bg_mean_above - clamp_st.p_bg_mean_above) << ",\n"
       << "    \"p_mom_mean_above\": " << (signed_st.p_mom_mean_above - clamp_st.p_mom_mean_above) << ",\n"
       << "    \"a_force_z_mean_above\": " << (signed_st.a_force_z_mean_above - clamp_st.a_force_z_mean_above)
       << ",\n"
       << "    \"a_net_z_mean_above\": " << (signed_st.a_net_z_mean_above - clamp_st.a_net_z_mean_above) << ",\n"
       << "    \"a_force_z_rms_above\": " << (signed_st.a_force_z_rms_above - clamp_st.a_force_z_rms_above) << "\n"
       << "  }\n"
       << "}\n";
}

/**
 * Doc221 N2/N3: tip-local soft seal — for r < tip_r in surface band, scale upward Uz
 * (uz_scale=0 → kill/N2; 0<scale<1 → dampen/N3; 1 → no-op). No z clamp. Opt-in only.
 */
inline size_t hostScaleUpwardVelocityNearTip(BaseParticles &particles, Real z_lid, Real band,
                                              Real tip_r, Real uz_scale)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_lo = z_lid - band;
    const Real tip_r2 = tip_r * tip_r;
    const Real scale = uz_scale < Real(0) ? Real(0) : (uz_scale > Real(1) ? Real(1) : uz_scale);
    size_t n_touch = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] <= z_lo || vel[i][2] <= Real(0))
            continue;
        const Real r2 = pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1];
        if (r2 < tip_r2)
        {
            vel[i][2] *= scale;
            ++n_touch;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
    return n_touch;
}

/**
 * Doc221 L: clamp fluid z <= z_lid after position update (hard lid).
 * Returns number of particles clamped.
 */
inline size_t hostClampZToLid(BaseParticles &particles, Real z_lid)
{
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    size_t n_clamp = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (pos[i][2] > z_lid)
        {
            pos[i][2] = z_lid;
            ++n_clamp;
        }
    }
    syncVariableToDevice<Vecd>(particles, "Position");
    return n_clamp;
}

/**
 * Doc221 M: one-shot post-settle hard-lid reset — clamp z<=z_lid then zero all velocity.
 * Avoids keep-start polarity inversion from slamming a raised free surface into the lid.
 */
inline size_t hostResetHardLidGeometry(BaseParticles &particles, Real z_lid)
{
    const size_t n_clamp = hostClampZToLid(particles, z_lid);
    hostZeroVelocity(particles);
    return n_clamp;
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

/** Doc220 H0: per-particle hydrostatic fields keyed by OriginalID (same-state stages). */
inline void hostDoc220DumpHydroParticleCsv(BaseParticles &particles, const std::string &path, Real rho0, Real c0,
                                          Real gravity_g, Real z_free, const Vecd &center, Real radius,
                                          Real z_bottom, Real band, bool has_B)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const Matd *B = nullptr;
    if (has_B)
    {
        syncVariableToHost<Matd>(particles, "LinearCorrectionMatrix");
        B = particles.getVariableDataByName<Matd>("LinearCorrectionMatrix");
    }
    const size_t n = particles.TotalRealParticles();
    const Real c2 = c0 * c0;
    std::ofstream out(path);
    out << "OriginalID,x,y,z,depth,r,zone,rho,rho_sum,packing,p_eos,p_bg,p_total,p_h_ref,dp_h,"
           "mass,Fx,Fy,Fz,FPx,FPy,FPz,ax_net,ay_net,az_net,ux,uy,uz,drho_dt,B00,B11,B22,incomplete\n";
    out.setf(std::ios::scientific);
    out.precision(17);
    for (size_t i = 0; i < n; ++i)
    {
        const Real m = mass[i] + TinyReal;
        const Real depth = SMAX(z_free - pos[i][2], Real(0));
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real p_eos = c2 * (rho[i] - rho0);
        const Real p_tot = p_eos + p_bg[i];
        const Real p_h = rho0 * gravity_g * depth;
        const Real packing = rho_sum[i] / (rho0 + TinyReal);
        const bool incomplete = rho_sum[i] + TinyReal < rho0;
        const FrenchHydroZone zone =
            hostClassifyFrenchHydroZone(pos[i], center, radius, z_bottom, z_free, band);
        const char *zone_name = "Interior";
        if (zone == FrenchHydroZone::Surface)
            zone_name = "Surface";
        else if (zone == FrenchHydroZone::WallSide)
            zone_name = "WallSide";
        else if (zone == FrenchHydroZone::WallBottom)
            zone_name = "WallBottom";
        else if (zone == FrenchHydroZone::WallCorner)
            zone_name = "WallCorner";
        else if (zone == FrenchHydroZone::WallSurface)
            zone_name = "WallSurface";
        const Real B00 = B ? B[i](0, 0) : Real(1);
        const Real B11 = B ? B[i](1, 1) : Real(1);
        const Real B22 = B ? B[i](2, 2) : Real(1);
        out << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << depth << "," << r
            << "," << zone_name << "," << rho[i] << "," << rho_sum[i] << "," << packing << "," << p_eos << ","
            << p_bg[i] << "," << p_tot << "," << p_h << "," << (p_tot - p_h) << "," << mass[i] << "," << f[i][0]
            << "," << f[i][1] << "," << f[i][2] << "," << fp[i][0] << "," << fp[i][1] << "," << fp[i][2] << ","
            << (f[i][0] + fp[i][0]) / m << "," << (f[i][1] + fp[i][1]) / m << "," << (f[i][2] + fp[i][2]) / m
            << "," << vel[i][0] << "," << vel[i][1] << "," << vel[i][2] << "," << drho[i] << "," << B00 << ","
            << B11 << "," << B22 << "," << (incomplete ? 1 : 0) << "\n";
    }
}

/** Doc220 S0: dump particle state keyed for OriginalID join. */
inline void hostDoc220DumpSortAuditCsv(BaseParticles &particles, const std::string &path, Real t0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const size_t n = particles.TotalRealParticles();
    std::ofstream out(path);
    out << "index,OriginalID,x,y,z,T,dT,T_minus_T0_dT,Mass,Density\n";
    out.setf(std::ios::scientific);
    out.precision(17);
    for (size_t i = 0; i < n; ++i)
    {
        const Real id_err = T[i] - (t0 + dT[i]);
        out << i << "," << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << T[i] << ","
            << dT[i] << "," << id_err << "," << mass[i] << "," << rho[i] << "\n";
    }
}

struct Doc220SortAuditCompare
{
    size_t n = 0;
    size_t n_matched = 0;
    size_t n_pos_mismatch = 0;
    size_t n_T_mismatch = 0;
    size_t n_dT_mismatch = 0;
    size_t n_TdT_identity_fail = 0;
    size_t n_mass_mismatch = 0;
    Real max_abs_dT = 0;
    Real max_abs_dPos = 0;
    Real T_mean_before = 0;
    Real T_mean_after = 0;
    Real T_min_before = 0;
    Real T_max_before = 0;
    Real T_min_after = 0;
    Real T_max_after = 0;
    bool temperature_in_evolving = false;
    bool delta_t_in_evolving = false;
    bool passed = false;
};

inline bool hostDiscreteVarInEvolving(BaseParticles &particles, const std::string &name)
{
    auto &evolving = particles.EvolvingVariables();
    constexpr int type_index = DataTypeIndex<Real>::value;
    auto &keepers = std::get<type_index>(evolving);
    for (auto *dv : keepers)
    {
        if (dv != nullptr && dv->Name() == name)
            return true;
    }
    return false;
}

inline Doc220SortAuditCompare hostDoc220CompareSortAuditCsv(const std::string &before_path,
                                                            const std::string &after_path, BaseParticles &particles,
                                                            Real t0)
{
    Doc220SortAuditCompare c;
    c.temperature_in_evolving = hostDiscreteVarInEvolving(particles, kOphelieTemperatureField);
    c.delta_t_in_evolving = hostDiscreteVarInEvolving(particles, kOphelieThermalDeltaTField);

    struct Row
    {
        Vecd pos{};
        Real T = 0;
        Real dT = 0;
        Real mass = 0;
        Real rho = 0;
    };
    auto load = [](const std::string &path, std::unordered_map<UnsignedInt, Row> &by_id, Real &tmin, Real &tmax,
                   Real &tmean) {
        std::ifstream in(path);
        std::string line;
        std::getline(in, line); // header
        size_t n = 0;
        double tsum = 0;
        tmin = std::numeric_limits<Real>::max();
        tmax = -std::numeric_limits<Real>::max();
        while (std::getline(in, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);
            std::string cell;
            UnsignedInt oid = 0;
            Row r;
            std::getline(ss, cell, ','); // index
            std::getline(ss, cell, ',');
            oid = static_cast<UnsignedInt>(std::stoull(cell));
            std::getline(ss, cell, ',');
            r.pos[0] = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ',');
            r.pos[1] = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ',');
            r.pos[2] = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ',');
            r.T = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ',');
            r.dT = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ','); // identity err
            std::getline(ss, cell, ',');
            r.mass = static_cast<Real>(std::stod(cell));
            std::getline(ss, cell, ',');
            r.rho = static_cast<Real>(std::stod(cell));
            by_id[oid] = r;
            tsum += r.T;
            tmin = std::min(tmin, r.T);
            tmax = std::max(tmax, r.T);
            ++n;
        }
        tmean = n ? static_cast<Real>(tsum / static_cast<double>(n)) : Real(0);
        return n;
    };

    std::unordered_map<UnsignedInt, Row> before, after;
    c.n = load(before_path, before, c.T_min_before, c.T_max_before, c.T_mean_before);
    const size_t n_after = load(after_path, after, c.T_min_after, c.T_max_after, c.T_mean_after);
    (void)n_after;
    const Real pos_tol = Real(0); // pure permute: expect bit-identical after ID join
    const Real t_tol = Real(0);
    for (const auto &kv : before)
    {
        auto it = after.find(kv.first);
        if (it == after.end())
            continue;
        ++c.n_matched;
        const Row &a = kv.second;
        const Row &b = it->second;
        const Real dpos = (a.pos - b.pos).norm();
        c.max_abs_dPos = std::max(c.max_abs_dPos, dpos);
        if (dpos > pos_tol)
            ++c.n_pos_mismatch;
        const Real dT = std::abs(a.T - b.T);
        c.max_abs_dT = std::max(c.max_abs_dT, dT);
        if (dT > t_tol)
            ++c.n_T_mismatch;
        if (std::abs(a.dT - b.dT) > t_tol)
            ++c.n_dT_mismatch;
        if (std::abs(b.T - (t0 + b.dT)) > Real(1.0e-6) * (std::abs(t0) + Real(1)))
            ++c.n_TdT_identity_fail;
        if (std::abs(a.mass - b.mass) > t_tol)
            ++c.n_mass_mismatch;
    }
    c.passed = c.temperature_in_evolving && c.delta_t_in_evolving && c.n_matched == c.n && c.n_pos_mismatch == 0 &&
               c.n_T_mismatch == 0 && c.n_dT_mismatch == 0 && c.n_TdT_identity_fail == 0 && c.n_mass_mismatch == 0;
    return c;
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

/**
 * Doc246 Q1_UPPER diagnostic: keep raw EM Q in JouleHeatEM; write controlled Q to JouleHeat.
 * z measured from glass bottom; H = melt height. Does not modify EM solver.
 */
inline void applyDoc246QUpperControlled(BaseParticles &particles, Real kappa, Real glass_bottom_z, Real melt_h,
                                        Real target_power_w, bool verbose = true)
{
    if (!(kappa > TinyReal) || !(melt_h > TinyReal) || !(target_power_w > TinyReal))
        return;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "JouleHeat");
    syncVariableToHost<Real>(particles, "JouleHeatEM");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *q = particles.getVariableDataByName<Real>("JouleHeat");
    Real *q_em = particles.getVariableDataByName<Real>("JouleHeatEM");
    const size_t n = particles.TotalRealParticles();
    Real p_em = Real(0), mom_z_em = Real(0), p_w = Real(0), p_upper_em = Real(0);
    const Real z_mid = glass_bottom_z + Real(0.5) * melt_h;
    for (size_t i = 0; i < n; ++i)
    {
        q_em[i] = q[i];
        const Real z = pos[i][2];
        const Real zi = (z - glass_bottom_z) / melt_h; // 0 at bottom, 1 at top
        const Real w = std::exp(kappa * (zi - Real(0.5)));
        p_em += q_em[i] * vol[i];
        mom_z_em += q_em[i] * vol[i] * z;
        p_w += q_em[i] * w * vol[i];
        if (z >= z_mid)
            p_upper_em += q_em[i] * vol[i];
    }
    const Real C = (p_w > TinyReal) ? (target_power_w / p_w) : Real(0);
    Real p_ctrl = Real(0), mom_z_ctrl = Real(0), p_upper_ctrl = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        const Real zi = (z - glass_bottom_z) / melt_h;
        const Real w = std::exp(kappa * (zi - Real(0.5)));
        q[i] = C * q_em[i] * w;
        p_ctrl += q[i] * vol[i];
        mom_z_ctrl += q[i] * vol[i] * z;
        if (z >= z_mid)
            p_upper_ctrl += q[i] * vol[i];
    }
    syncVariableToDevice<Real>(particles, "JouleHeat");
    syncVariableToDevice<Real>(particles, "JouleHeatEM");
    if (verbose)
    {
        const Real zq_em = (p_em > TinyReal) ? mom_z_em / p_em : Real(0);
        const Real zq_c = (p_ctrl > TinyReal) ? mom_z_ctrl / p_ctrl : Real(0);
        std::cout << "[ophelie][doc246][Q1_UPPER] kappa=" << kappa << " C=" << C
                  << " P_EM=" << p_em << " P_ctrl=" << p_ctrl << " target=" << target_power_w
                  << " zQ_EM=" << zq_em << " zQ_ctrl=" << zq_c
                  << " upper_frac_EM=" << ((p_em > TinyReal) ? p_upper_em / p_em : Real(0))
                  << " upper_frac_ctrl=" << ((p_ctrl > TinyReal) ? p_upper_ctrl / p_ctrl : Real(0))
                  << " (controlled_EM_shape; NOT paper Maxwell source)\n";
    }
}

/** Doc219 E: re-sample spatial Q; optionally renormalize ΣViQ to target (diagnostic protocol). */
inline Real hostSampleSpatialFrozenQ(BaseParticles &particles, const rh200::Rh200ScalarEulerianGrid &grid,
                                     Real target_power_w, bool renormalize, size_t &n_uncovered,
                                     bool verbose = true)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "JouleHeat");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *q = particles.getVariableDataByName<Real>("JouleHeat");
    const size_t n = particles.TotalRealParticles();
    n_uncovered = 0;
    Real p = Real(0);
    const Vecd lo = grid.spec_.lower_bound_;
    const Vecd hi = grid.spec_.lower_bound_ +
                    Vecd(static_cast<Real>(grid.spec_.nx_ - 1), static_cast<Real>(grid.spec_.ny_ - 1),
                         static_cast<Real>(grid.spec_.nz_ - 1)) *
                        grid.spec_.spacing_;
    for (size_t i = 0; i < n; ++i)
    {
        Vecd x = pos[i];
        if (!grid.inside(x))
        {
            ++n_uncovered;
            // Clamp to grid box and sample edge value — never silent Q=0 for free-surface overshoot.
            for (int d = 0; d < 3; ++d)
                x[d] = std::min(hi[d], std::max(lo[d], x[d]));
        }
        q[i] = grid.sampleHost(x);
        p += q[i] * Vol[i];
    }
    Real scale = Real(1);
    if (renormalize && p > TinyReal && target_power_w > TinyReal)
    {
        scale = target_power_w / p;
        for (size_t i = 0; i < n; ++i)
            q[i] *= scale;
        p = target_power_w;
    }
    syncVariableToDevice<Real>(particles, "JouleHeat");
    if (verbose)
    {
        std::cout << "[ophelie][doc219][spatial-q] P_sample=" << p << " remap_scale=" << scale
                  << " n_uncovered=" << n_uncovered << "/" << n
                  << " renormalize=" << (renormalize ? 1 : 0) << std::endl;
    }
    return p;
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

/**
 * Free-surface melt: p_h = ρ0 g (z_free − z).
 * Doc222 E1: signed_depth=true (default) keeps ∇p_h above the reference plane for restoring.
 * signed_depth=false restores legacy max(0,·) clamp (D2 FAIL path).
 */
inline void hostApplyFrenchBackgroundHydroPressure(BaseParticles &particles, Real rho0, Real g, Real z_free,
                                                   bool signed_depth = true)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real raw = z_free - pos[i][2];
        const Real depth = signed_depth ? raw : std::max(Real(0), raw);
        p_bg[i] = rho0 * g * depth;
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
}

/** Doc222 E2b diagnostics for local-η BackgroundPressure write. */
struct Doc222LocalEtaPhReport
{
    int n_bins = 0;
    int n_bins_valid = 0;
    Real eta_min = 0;
    Real eta_max = 0;
    Real eta_mean = 0;
    Real eta_spread = 0; // eta_max - eta_min
};

/**
 * Doc222 E2b: p_h = ρ0 g (η(r)−z) with η = radial-bin free-surface z_max.
 * Empty bins inherit nearest valid η (else z_free_fallback). Gives ∂_r p_h ∝ ∂_r η
 * so free-surface slopes can restore under omit-GravityForceCK.
 */
inline Doc222LocalEtaPhReport hostApplyFrenchBackgroundHydroPressureLocalEta(
    BaseParticles &particles, Real rho0, Real g, Real z_free_fallback, const Vecd &center, Real R,
    int n_bins, bool signed_depth)
{
    Doc222LocalEtaPhReport rep;
    if (n_bins < 4)
        n_bins = 4;
    rep.n_bins = n_bins;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    std::vector<Real> bin_eta(static_cast<size_t>(n_bins), -Real(1.0e30));
    std::vector<int> bin_count(static_cast<size_t>(n_bins), 0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        if (r >= R || R <= TinyReal)
            continue;
        int b = static_cast<int>(std::floor(r / R * Real(n_bins)));
        if (b < 0)
            b = 0;
        if (b >= n_bins)
            b = n_bins - 1;
        bin_eta[static_cast<size_t>(b)] = std::max(bin_eta[static_cast<size_t>(b)], pos[i][2]);
        bin_count[static_cast<size_t>(b)]++;
    }
    // Fill empty bins from nearest valid neighbor; else fall back to fixed z_free.
    for (int b = 0; b < n_bins; ++b)
    {
        if (bin_count[static_cast<size_t>(b)] > 0 && bin_eta[static_cast<size_t>(b)] > -Real(1.0e20))
        {
            rep.n_bins_valid++;
            continue;
        }
        Real best = z_free_fallback;
        int best_dist = n_bins + 1;
        for (int j = 0; j < n_bins; ++j)
        {
            if (bin_count[static_cast<size_t>(j)] <= 0 || bin_eta[static_cast<size_t>(j)] < -Real(1.0e20))
                continue;
            const int d = std::abs(j - b);
            if (d < best_dist)
            {
                best_dist = d;
                best = bin_eta[static_cast<size_t>(j)];
            }
        }
        bin_eta[static_cast<size_t>(b)] = best;
    }
    Real eta_sum = 0;
    rep.eta_min = Real(1.0e30);
    rep.eta_max = -Real(1.0e30);
    for (int b = 0; b < n_bins; ++b)
    {
        const Real e = bin_eta[static_cast<size_t>(b)];
        eta_sum += e;
        rep.eta_min = std::min(rep.eta_min, e);
        rep.eta_max = std::max(rep.eta_max, e);
    }
    rep.eta_mean = eta_sum / Real(n_bins);
    rep.eta_spread = rep.eta_max - rep.eta_min;

    auto eta_at_r = [&](Real r) -> Real {
        if (R <= TinyReal || n_bins <= 0)
            return z_free_fallback;
        Real rr = std::min(std::max(r, Real(0)), R * (Real(1) - Real(1.0e-6)));
        const Real fb = rr / R * Real(n_bins) - Real(0.5);
        int b0 = static_cast<int>(std::floor(fb));
        if (b0 < 0)
            b0 = 0;
        if (b0 >= n_bins - 1)
            return bin_eta[static_cast<size_t>(n_bins - 1)];
        const Real a = fb - Real(b0);
        return (Real(1) - a) * bin_eta[static_cast<size_t>(b0)] + a * bin_eta[static_cast<size_t>(b0 + 1)];
    };

    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real eta = eta_at_r(r);
        const Real raw = eta - pos[i][2];
        const Real depth = signed_depth ? raw : std::max(Real(0), raw);
        p_bg[i] = rho0 * g * depth;
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
    return rep;
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

/** Doc232 B + Doc248: fixed-space NC probes (legacy tip + midradius depth-12cm + interiors). */
inline StdVec<OphelieFixedSpaceProbeDef> makeDoc232ClientNcProbes(const OphelieFrenchReducedCaseParams &french)
{
    StdVec<OphelieFixedSpaceProbeDef> probes;
    const Real z_bot = french.glass_center[2] - french.glass_half_height;
    const Real H = Real(2) * french.glass_half_height;
    const Real z_mid = french.glass_center[2];
    const Real R = french.glass_radius;
    const OphelieClientProbeDefinition thermo =
        makeClientNominalThermocoupleProbe(makeClientNaturalConvectionAssemblyCatalog());
    probes.push_back({"legacy_centerline_probes_tip_proxy",
                      Vecd(thermo.solver_m.x, thermo.solver_m.y, thermo.solver_m.z), thermo.status, true,
                      "CAD AABB tip; isolate from Doc248 midradius (legacy naming)"});
    const OphelieClientProbeDefinition mid =
        makeClientMidradiusDepth12cmProbe(R, z_bot, H, Real(1273.15));
    probes.push_back({"client_midradius_probe", Vecd(mid.solver_m.x, mid.solver_m.y, mid.solver_m.z), mid.status,
                      true, "Doc248 r=0.5R depth=0.12m below free surface; +X azimuth ASSUMED"});
    probes.push_back({"dev_axis_mid", Vecd(0, 0, z_mid), "development", true, "axis mid-height"});
    probes.push_back(
        {"dev_r_half_mid", Vecd(Real(0.5) * R, 0, z_mid), "development", true, "mid-radius mid-height"});
    probes.push_back({"dev_near_wall_mid", Vecd(std::max(R - Real(2) * french.dp, Real(0.5) * R), 0, z_mid),
                      "development", true, "near-wall diagnostic (not surface T)"});
    probes.push_back({"dev_axis_lower_quartile", Vecd(0, 0, z_bot + Real(0.25) * H), "development", true,
                      "lower-quartile axis"});
    return probes;
}

inline void writeProbesTimeseriesHeader(std::ostream &os)
{
    os << "time_s,probe_id,x,y,z,T_K,T_C,n_neighbors,support_h_m,weight_sum,nearest_dist_m,"
          "coverage_valid,location_status,config_note\n";
}

inline void appendProbesTimeseries(std::ostream &os, Real t, BaseParticles &particles,
                                   const StdVec<OphelieFixedSpaceProbeDef> &probes, Real kernel_h)
{
    syncFixedSpaceProbeFieldsToHost(particles, kOphelieTemperatureField);
    for (const auto &p : probes)
    {
        const auto s = hostSampleFixedSpaceTemperature(particles, p.pos, kernel_h, kOphelieTemperatureField);
        const Real T_C = s.valid ? (s.T - Real(273.15)) : std::numeric_limits<Real>::quiet_NaN();
        os << std::setprecision(10) << t << "," << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2]
           << "," << s.T << "," << T_C << "," << s.n_neighbors << "," << kernel_h << "," << s.weight_sum << ","
           << s.nearest_dist << "," << (s.valid ? 1 : 0) << "," << p.location_status << "," << p.note << "\n";
    }
    os.flush();
}

/**
 * Doc232 B2 manufacture gate on the production sampler. Restores Temperature/ΔT afterward.
 * Returns true if uniform+linear interior gates pass.
 */
inline bool runDoc232ProbePathVerify(BaseParticles &particles, const StdVec<OphelieFixedSpaceProbeDef> &probes,
                                     Real kernel_h, Real z_bot, Real melt_h, Real t0_restore,
                                     const std::string &csv_path, bool &uniform_ok, bool &linear_ok)
{
    uniform_ok = true;
    linear_ok = true;
    syncFixedSpaceProbeFieldsToHost(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Vecd>(particles, "Position");
    Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> T_backup(n), dT_backup(n);
    for (size_t i = 0; i < n; ++i)
    {
        T_backup[i] = T[i];
        dT_backup[i] = dT[i];
    }

    const Real T_uniform = Real(1273.15);
    for (size_t i = 0; i < n; ++i)
    {
        T[i] = T_uniform;
        dT[i] = T_uniform - t0_restore;
    }

    fs::create_directories(fs::path(csv_path).parent_path());
    std::ofstream csv(csv_path);
    csv << std::setprecision(10);
    csv << "probe_id,x,y,z,field,T_exact,T_sample,abs_err,rel_err_vs_amp,n_neighbors,weight_sum,"
           "nearest_dist_m,valid,pass_row,location_status\n";

    const Real T_bot = Real(1100);
    const Real T_top = Real(1300);
    const Real amp = T_top - T_bot;
    auto T_linear = [&](Real z) {
        const Real a = (z - z_bot) / std::max(melt_h, TinyReal);
        return T_bot + amp * std::min(Real(1), std::max(Real(0), a));
    };

    for (const auto &p : probes)
    {
        const auto s = hostSampleFixedSpaceTemperature(particles, p.pos, kernel_h, kOphelieTemperatureField);
        const Real abs_err = std::abs(s.T - T_uniform);
        const bool row_pass = s.valid && (abs_err <= Real(1.0e-3)) && (s.n_neighbors >= 4);
        if (!row_pass)
            uniform_ok = false;
        csv << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ",uniform," << T_uniform << ","
            << s.T << "," << abs_err << "," << (abs_err / std::max(amp, Real(1))) << "," << s.n_neighbors << ","
            << s.weight_sum << "," << s.nearest_dist << "," << (s.valid ? 1 : 0) << "," << (row_pass ? 1 : 0)
            << "," << p.location_status << "\n";
    }

    for (size_t i = 0; i < n; ++i)
    {
        T[i] = T_linear(pos[i][2]);
        dT[i] = T[i] - t0_restore;
    }
    for (const auto &p : probes)
    {
        if (!p.expect_interior)
            continue;
        const auto s = hostSampleFixedSpaceTemperature(particles, p.pos, kernel_h, kOphelieTemperatureField);
        const Real T_ex = T_linear(p.pos[2]);
        const Real abs_err = std::abs(s.T - T_ex);
        const Real rel = abs_err / std::max(amp, TinyReal);
        const bool row_pass = s.valid && (rel <= Real(0.01)) && (s.n_neighbors >= 4);
        if (!row_pass)
            linear_ok = false;
        csv << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ",linear_z," << T_ex << ","
            << s.T << "," << abs_err << "," << rel << "," << s.n_neighbors << "," << s.weight_sum << ","
            << s.nearest_dist << "," << (s.valid ? 1 : 0) << "," << (row_pass ? 1 : 0) << ","
            << p.location_status << "\n";
    }

    for (size_t i = 0; i < n; ++i)
    {
        T[i] = T_backup[i];
        dT[i] = dT_backup[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    std::cout << "[ophelie][doc232][B2] probe_path_verify uniform_ok=" << (uniform_ok ? 1 : 0)
              << " linear_ok=" << (linear_ok ? 1 : 0) << " csv=" << csv_path << "\n";
    return uniform_ok && linear_ok;
}

/** Doc232 B2: dump host Position/Temperature/Mass for offline probe recompute. */
inline void dumpDoc232ParticleSnapshotCsv(BaseParticles &particles, const std::string &csv_path, Real t_s,
                                          Real kernel_h)
{
    syncFixedSpaceProbeFieldsToHost(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "Mass");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    fs::create_directories(fs::path(csv_path).parent_path());
    std::ofstream os(csv_path);
    os << std::setprecision(10);
    os << "# t_s=" << t_s << " kernel_h_m=" << kernel_h << " n=" << n << "\n";
    os << "x,y,z,T_K,Mass\n";
    for (size_t i = 0; i < n; ++i)
    {
        os << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << temp[i] << "," << mass[i] << "\n";
    }
    std::cout << "[ophelie][doc232][B2] particle_snapshot csv=" << csv_path << " t=" << t_s << " n=" << n
              << "\n";
}

/** Doc234 D: particle field with velocity for ON/OFF Eulerian postprocess. */
inline void dumpDoc234ParticleFieldCsv(BaseParticles &particles, const std::string &csv_path, Real t_s)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    fs::create_directories(fs::path(csv_path).parent_path());
    std::ofstream os(csv_path);
    os << std::setprecision(10);
    os << "# t_s=" << t_s << " n=" << n << "\n";
    os << "particle_i,x,y,z,vx,vy,vz,T_K,Mass,Vol\n";
    for (size_t i = 0; i < n; ++i)
    {
        os << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << vel[i][0] << ","
           << vel[i][1] << "," << vel[i][2] << "," << temp[i] << "," << mass[i] << "," << vol[i] << "\n";
    }
    std::cout << "[ophelie][doc234][D] particle_field csv=" << csv_path << " t=" << t_s << " n=" << n << "\n";
}

/**
 * Doc235 B: event snapshot for OFF-background attribution (read-only dump).
 * Core mechanics fields used in settle cleanup / first production steps.
 */
inline void dumpDoc235OffEventCsv(BaseParticles &particles, const std::string &csv_path, const std::string &event,
                                  Real t_s, const Vecd &axis_origin, Real /*radius_m*/, Real z_bot, Real z_top)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Pressure");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *dens = particles.getVariableDataByName<Real>("Density");
    const Real *pres = particles.getVariableDataByName<Real>("Pressure");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *ff = particles.getVariableDataByName<Vecd>("Force");
    const Real *dcr = particles.getVariableDataByName<Real>("DensityChangeRate");
    const size_t n = particles.TotalRealParticles();
    fs::create_directories(fs::path(csv_path).parent_path());
    std::ofstream os(csv_path);
    os << std::setprecision(10);
    os << "# event=" << event << " t_s=" << t_s << " n=" << n << "\n";
    os << "particle_i,x,y,z,r,vx,vy,vz,Mass,Vol,Density,Pressure,FPx,FPy,FPz,Fx,Fy,Fz,DensityChangeRate\n";
    double ke = 0.0, sum_m = 0.0, sum_v = 0.0, umax2 = 0.0;
    double rho_min = 0.0, rho_max = 0.0;
    Vecd mom = Vecd::Zero();
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - axis_origin[0];
        const Real dy = pos[i][1] - axis_origin[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        os << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << r << "," << vel[i][0]
           << "," << vel[i][1] << "," << vel[i][2] << "," << mass[i] << "," << vol[i] << "," << dens[i] << ","
           << pres[i] << "," << fp[i][0] << "," << fp[i][1] << "," << fp[i][2] << "," << ff[i][0] << ","
           << ff[i][1] << "," << ff[i][2] << "," << dcr[i] << "\n";
        const double m = static_cast<double>(mass[i]);
        const double u2 = static_cast<double>(vel[i].squaredNorm());
        ke += 0.5 * m * u2;
        sum_m += m;
        sum_v += static_cast<double>(vol[i]);
        umax2 = std::max(umax2, u2);
        mom += mass[i] * vel[i];
        if (i == 0)
        {
            rho_min = dens[i];
            rho_max = dens[i];
        }
        else
        {
            rho_min = std::min(rho_min, static_cast<double>(dens[i]));
            rho_max = std::max(rho_max, static_cast<double>(dens[i]));
        }
    }
    const double U_rms = (sum_m > 0.0) ? std::sqrt(2.0 * ke / sum_m) : 0.0;
    const double U_max = std::sqrt(umax2);
    std::cout << "[ophelie][doc235][B] event=" << event << " t=" << t_s << " csv=" << csv_path
              << " U_rms=" << U_rms << " U_max=" << U_max << " KE=" << ke << " sum_m=" << sum_m
              << " sum_V=" << sum_v << " rho=[" << rho_min << "," << rho_max << "] mom=(" << mom[0] << ","
              << mom[1] << "," << mom[2] << ") z_span_geom=[" << z_bot << "," << z_top << "]\n";
    const fs::path stats = fs::path(csv_path).parent_path() / "event_stats.jsonl";
    std::ofstream js(stats, std::ios::app);
    js << std::setprecision(10);
    js << "{\"event\":\"" << event << "\",\"t_s\":" << t_s << ",\"U_rms\":" << U_rms << ",\"U_max\":" << U_max
       << ",\"KE\":" << ke << ",\"sum_m\":" << sum_m << ",\"sum_V\":" << sum_v << ",\"rho_min\":" << rho_min
       << ",\"rho_max\":" << rho_max << ",\"csv\":\"" << csv_path << "\"}\n";
}

inline Real hostMassWeightedMeanDeltaT(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "Mass");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    double sw = 0.0, s = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double m = static_cast<double>(mass[i]);
        sw += m;
        s += m * static_cast<double>(dT[i]);
    }
    return (sw > 0.0) ? static_cast<Real>(s / sw) : Real(0);
}

/** Independent energy from Mass * cp * ΔT (Doc232 C reference path). */
inline Real hostThermalEnergyFromDeltaTMass(BaseParticles &particles, Real cp)
{
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "Mass");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    const double cp_d = static_cast<double>(cp);
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
        e += static_cast<double>(mass[i]) * cp_d * static_cast<double>(dT[i]);
    return static_cast<Real>(e);
}

inline Real hostSumParticleMass(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "Mass");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    double s = 0.0;
    for (size_t i = 0; i < n; ++i)
        s += static_cast<double>(mass[i]);
    return static_cast<Real>(s);
}

/**
 * Doc232 D: production BC face isolation + single-step energy vs host power ledger.
 * Independent geometric reference is written for offline Python; this gate also reports
 * production vs closed-form analytic for immediate triage (≤2% power, ≤1% energy).
 */
inline bool runDoc232BoundaryFaceAudit(SPHBody &glass_body, BaseParticles &particles,
                                       OphelieThermalDiffusionOneWayOptions thermal_bc_base,
                                       const OphelieFrenchReducedCaseParams &french, Real rho0, Real cp, Real t0,
                                       const std::string &out_dir,
                                       const std::string &face_t_mode = "global_dt_dz",
                                       bool doc233_scope_fields = false)
{
    fs::create_directories(out_dir);
    const Real z_bot = french.glass_center[2] - french.glass_half_height;
    const Real H = Real(2) * french.glass_half_height;
    const Real R = french.glass_radius;
    const Real A_side = Real(2) * Real(M_PI) * R * H;
    const Real A_disk = Real(M_PI) * R * R;
    const Real shell = thermal_bc_base.boundary_shell_thickness;
    const Real dt = Real(1.0e-3);
    const Real T_uniform = Real(1273.15);
    const Real T_lo = Real(1100);
    const Real T_hi = Real(1300);
    const bool extrapolate_face = (face_t_mode != "band");
    std::cout << "[ophelie][doc233][D] boundary_face_t_mode=" << face_t_mode
              << " extrapolate_face_t=" << (extrapolate_face ? 1 : 0)
              << " scope_fields=" << (doc233_scope_fields ? 1 : 0) << "\n";

    const OphelieFrenchBoundaryAreaAudit area =
        hostOphelieThermalFrenchNaturalBoundaryAreaAudit(particles, thermal_bc_base, french, shell);
    logOphelieFrenchBoundaryAreaAudit(area, "doc232_D_preface");
    {
        std::ofstream as(fs::path(out_dir) / "boundary_area_audit.csv");
        as << std::setprecision(10);
        as << "quantity,value\n";
        as << "R_m," << R << "\nH_m," << H << "\nshell_m," << shell << "\n";
        as << "A_geom_side," << area.geometric_area_side << "\nA_geom_bottom," << area.geometric_area_bottom
           << "\nA_geom_free," << area.geometric_area_free << "\n";
        as << "A_disc_side," << area.discrete_area_side << "\nA_disc_bottom," << area.discrete_area_bottom
           << "\nA_disc_free," << area.discrete_area_free << "\n";
        as << "A_ratio_side," << area.area_side_ratio << "\nA_ratio_bottom," << area.area_bottom_ratio
           << "\nA_ratio_free," << area.area_free_ratio << "\n";
        as << "n_side," << area.n_side << "\nn_bottom," << area.n_bottom << "\nn_free," << area.n_free
           << "\nn_corner," << area.n_corner << "\nn_multi_face," << area.n_multi_face << "\n";
        as << "note,shell_band_T_is_not_exact_wall_T; flux uses particle T in V/delta shell\n";
    }
    {
        std::ofstream gj(fs::path(out_dir) / "geometry_ref.json");
        gj << std::setprecision(10);
        gj << "{\n"
           << "  \"R_m\": " << R << ",\n"
           << "  \"H_m\": " << H << ",\n"
           << "  \"z_bot_m\": " << z_bot << ",\n"
           << "  \"h_side\": " << thermal_bc_base.h_side << ",\n"
           << "  \"h_bottom\": " << thermal_bc_base.h_bottom << ",\n"
           << "  \"h_free\": " << thermal_bc_base.h_free << ",\n"
           << "  \"emissivity\": " << thermal_bc_base.emissivity << ",\n"
           << "  \"t_cool\": " << thermal_bc_base.t_cool << ",\n"
           << "  \"t_ambient\": " << thermal_bc_base.t_ambient << ",\n"
           << "  \"t_rad_ambient\": " << thermal_bc_base.t_rad_ambient << ",\n"
           << "  \"stefan_boltzmann\": " << thermal_bc_base.stefan_boltzmann << ",\n"
           << "  \"T_uniform\": " << T_uniform << ",\n"
           << "  \"T_lo\": " << T_lo << ",\n"
           << "  \"T_hi\": " << T_hi << ",\n"
           << "  \"dt_s\": " << dt << ",\n"
           << "  \"rho0\": " << rho0 << ",\n"
           << "  \"cp\": " << cp << "\n"
           << "}\n";
    }

    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Vecd>(particles, "Position");
    Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> T_backup(n), dT_backup(n);
    for (size_t i = 0; i < n; ++i)
    {
        T_backup[i] = T[i];
        dT_backup[i] = dT[i];
    }

    auto set_uniform = [&]() {
        for (size_t i = 0; i < n; ++i)
        {
            T[i] = T_uniform;
            dT[i] = T_uniform - t0;
        }
        syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
        syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    };
    auto set_linear_z = [&]() {
        for (size_t i = 0; i < n; ++i)
        {
            const Real a = (pos[i][2] - z_bot) / std::max(H, TinyReal);
            const Real Ti = T_lo + (T_hi - T_lo) * std::min(Real(1), std::max(Real(0), a));
            T[i] = Ti;
            dT[i] = Ti - t0;
        }
        syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
        syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    };
    // Doc233 D: T=1200+100*(r/R)^2
    auto set_radial_quad = [&]() {
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - french.glass_center[0];
            const Real dy = pos[i][1] - french.glass_center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real rr = std::min(Real(1), r / std::max(R, TinyReal));
            const Real Ti = Real(1200) + Real(100) * rr * rr;
            T[i] = Ti;
            dT[i] = Ti - t0;
        }
        syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
        syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    };
    // Doc233 D: T=1200+100*(2*(z-zb)/H-1)^2 — symmetric; global dT/dz ~ 0 but faces have slope.
    auto set_parabolic_z = [&]() {
        for (size_t i = 0; i < n; ++i)
        {
            const Real a = (pos[i][2] - z_bot) / std::max(H, TinyReal);
            const Real s = Real(2) * std::min(Real(1), std::max(Real(0), a)) - Real(1);
            const Real Ti = Real(1200) + Real(100) * s * s;
            T[i] = Ti;
            dT[i] = Ti - t0;
        }
        syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
        syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    };
    // Doc235 F: T=1200+100*(1+0.5*(r/R)^2)*(2*(z-zb)/H-1)^2 — vertical slope varies with r.
    auto set_rz_quad = [&]() {
        for (size_t i = 0; i < n; ++i)
        {
            const Real dx = pos[i][0] - french.glass_center[0];
            const Real dy = pos[i][1] - french.glass_center[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real rr = std::min(Real(1), r / std::max(R, TinyReal));
            const Real a = (pos[i][2] - z_bot) / std::max(H, TinyReal);
            const Real s = Real(2) * std::min(Real(1), std::max(Real(0), a)) - Real(1);
            const Real amp = Real(1) + Real(0.5) * rr * rr;
            const Real Ti = Real(1200) + Real(100) * amp * s * s;
            T[i] = Ti;
            dT[i] = Ti - t0;
        }
        syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
        syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    };

    struct FaceCase
    {
        const char *face;
        Real h_side, h_bottom, h_free, emissivity;
    };
    const FaceCase faces[] = {
        {"side", thermal_bc_base.h_side, Real(0), Real(0), Real(0)},
        {"bottom", Real(0), thermal_bc_base.h_bottom, Real(0), Real(0)},
        {"free_conv", Real(0), Real(0), thermal_bc_base.h_free, Real(0)},
        {"free_rad", Real(0), Real(0), Real(0), thermal_bc_base.emissivity},
    };

    auto analytic = [&](const std::string &field, const FaceCase &fc) -> Real {
        const Real tc = thermal_bc_base.t_cool;
        const Real ta = thermal_bc_base.t_ambient;
        const Real tr = thermal_bc_base.t_rad_ambient;
        const Real sb = thermal_bc_base.stefan_boltzmann;
        if (field == "uniform")
        {
            if (std::strcmp(fc.face, "side") == 0)
                return fc.h_side * A_side * (T_uniform - tc);
            if (std::strcmp(fc.face, "bottom") == 0)
                return fc.h_bottom * A_disk * (T_uniform - tc);
            if (std::strcmp(fc.face, "free_conv") == 0)
                return fc.h_free * A_disk * (T_uniform - ta);
            const Real Tu2 = T_uniform * T_uniform;
            const Real Tr2 = tr * tr;
            return fc.emissivity * sb * A_disk * (Tu2 * Tu2 - Tr2 * Tr2);
        }
        if (field == "radial_quad")
        {
            // Geometric face T: side at r=R → 1300; disks mean T=1250.
            const Real T_side = Real(1300);
            const Real T_disk = Real(1250);
            if (std::strcmp(fc.face, "side") == 0)
                return fc.h_side * A_side * (T_side - tc);
            if (std::strcmp(fc.face, "bottom") == 0)
                return fc.h_bottom * A_disk * (T_disk - tc);
            if (std::strcmp(fc.face, "free_conv") == 0)
                return fc.h_free * A_disk * (T_disk - ta);
            // Radiation uses area integral of T^4; for T=1200+100(r/R)^2 closed form:
            // ∫ T^4 dA / A is not (mean T)^4 — use numerical ring sum for reference.
            const int nbin = 400;
            double acc = 0.0;
            for (int k = 0; k < nbin; ++k)
            {
                const double r0 = R * (static_cast<double>(k) / nbin);
                const double r1 = R * (static_cast<double>(k + 1) / nbin);
                const double rm = 0.5 * (r0 + r1);
                const double Tk = 1200.0 + 100.0 * (rm / R) * (rm / R);
                const double dA = M_PI * (r1 * r1 - r0 * r0);
                acc += dA * Tk * Tk * Tk * Tk;
            }
            const double Tr2 = static_cast<double>(tr) * static_cast<double>(tr);
            return static_cast<Real>(fc.emissivity * sb * (acc - A_disk * Tr2 * Tr2));
        }
        if (field == "parabolic_z")
        {
            // Side mean: 1200+100/3; faces at zbot/ztop: 1300.
            const Real T_side_mean = Real(1200) + Real(100) / Real(3);
            const Real T_face = Real(1300);
            if (std::strcmp(fc.face, "side") == 0)
                return fc.h_side * A_side * (T_side_mean - tc);
            if (std::strcmp(fc.face, "bottom") == 0)
                return fc.h_bottom * A_disk * (T_face - tc);
            if (std::strcmp(fc.face, "free_conv") == 0)
                return fc.h_free * A_disk * (T_face - ta);
            const Real Tf2 = T_face * T_face;
            const Real Tr2 = tr * tr;
            return fc.emissivity * sb * A_disk * (Tf2 * Tf2 - Tr2 * Tr2);
        }
        if (field == "rz_quad")
        {
            // Face T(r)=1300+50(r/R)^2; disk mean=1325. Side: 1200+150*s^2 → mean 1250.
            const Real T_side_mean = Real(1250);
            const Real T_disk_mean = Real(1325);
            if (std::strcmp(fc.face, "side") == 0)
                return fc.h_side * A_side * (T_side_mean - tc);
            if (std::strcmp(fc.face, "bottom") == 0)
                return fc.h_bottom * A_disk * (T_disk_mean - tc);
            if (std::strcmp(fc.face, "free_conv") == 0)
                return fc.h_free * A_disk * (T_disk_mean - ta);
            const int nbin = 400;
            double acc = 0.0;
            for (int k = 0; k < nbin; ++k)
            {
                const double r0 = R * (static_cast<double>(k) / nbin);
                const double r1 = R * (static_cast<double>(k + 1) / nbin);
                const double rm = 0.5 * (r0 + r1);
                const double Tk = 1300.0 + 50.0 * (rm / R) * (rm / R);
                const double dA = M_PI * (r1 * r1 - r0 * r0);
                acc += dA * Tk * Tk * Tk * Tk;
            }
            const double Tr2 = static_cast<double>(tr) * static_cast<double>(tr);
            return static_cast<Real>(fc.emissivity * sb * (acc - A_disk * Tr2 * Tr2));
        }
        // linear_z: surface integral on geometric cylinder (not particle weights).
        if (std::strcmp(fc.face, "side") == 0)
        {
            const Real T_mean = Real(0.5) * (T_lo + T_hi);
            return fc.h_side * A_side * (T_mean - tc);
        }
        if (std::strcmp(fc.face, "bottom") == 0)
            return fc.h_bottom * A_disk * (T_lo - tc);
        if (std::strcmp(fc.face, "free_conv") == 0)
            return fc.h_free * A_disk * (T_hi - ta);
        const Real Th2 = T_hi * T_hi;
        const Real Tr2 = tr * tr;
        return fc.emissivity * sb * A_disk * (Th2 * Th2 - Tr2 * Tr2);
    };

    std::ofstream pow_csv(fs::path(out_dir) / "boundary_face_powers.csv");
    pow_csv << std::setprecision(10);
    pow_csv << "field,face,P_production_W,P_analytic_geom_W,power_rel,DeltaH_J,E_ledger_J,energy_rel,"
               "pass_power,pass_energy,h_side,h_bottom,h_free,emissivity,dt_s\n";

    bool all_ok = true;
    const Real power_tol = Real(0.02);
    const Real energy_tol = Real(0.01);
    StdVec<std::string> fields = {"uniform", "linear_z"};
    if (doc233_scope_fields)
    {
        fields.push_back("radial_quad");
        fields.push_back("parabolic_z");
        fields.push_back("rz_quad"); // Doc235 F cross r–z manufactured field
    }

    for (const auto &field : fields)
    {
        for (const auto &fc : faces)
        {
            if (field == "uniform")
                set_uniform();
            else if (field == "linear_z")
                set_linear_z();
            else if (field == "radial_quad")
                set_radial_quad();
            else if (field == "rz_quad")
                set_rz_quad();
            else
                set_parabolic_z();

            OphelieThermalDiffusionOneWayOptions face_bc = thermal_bc_base;
            face_bc.h_side = fc.h_side;
            face_bc.h_bottom = fc.h_bottom;
            face_bc.h_free = fc.h_free;
            face_bc.emissivity = fc.emissivity;
            face_bc.enable_free_crust = false;

            Real dTdz_bot = Real(0), dTdz_top = Real(0);
            if (extrapolate_face)
            {
                if (face_t_mode == "local_dt_dz")
                {
                    const Real band = std::max(Real(2) * french.dp, shell);
                    const OphelieLocalFaceDTdz loc =
                        hostFitLocalFaceBandDTdz(particles, z_bot, z_bot + H, band, kOphelieTemperatureField);
                    dTdz_bot = loc.bottom_ok ? loc.dTdz_bottom : Real(0);
                    dTdz_top = loc.top_ok ? loc.dTdz_top : Real(0);
                    std::cout << "[ophelie][doc234][E] local_dt_dz field=" << field << " face=" << fc.face
                              << " band=" << band << " n_bot=" << loc.n_bottom << " ok_bot=" << (loc.bottom_ok ? 1 : 0)
                              << " dTdz_bot=" << dTdz_bot << " z_var_bot=" << loc.bottom_detail.z_var
                              << " det_bot=" << loc.bottom_detail.det
                              << " fail_bot=" << loc.bottom_detail.fail_reason << " n_top=" << loc.n_top
                              << " ok_top=" << (loc.top_ok ? 1 : 0) << " dTdz_top=" << dTdz_top
                              << " z_var_top=" << loc.top_detail.z_var << " det_top=" << loc.top_detail.det
                              << " fail_top=" << loc.top_detail.fail_reason << "\n";
                    if (!loc.bottom_ok)
                        std::cout << "[ophelie][doc235][F] fit_fallback face=bottom field=" << field
                                  << " reason=" << loc.bottom_detail.fail_reason << " slope_forced=0\n";
                    if (!loc.top_ok)
                        std::cout << "[ophelie][doc235][F] fit_fallback face=top field=" << field
                                  << " reason=" << loc.top_detail.fail_reason << " slope_forced=0\n";
                }
                else
                {
                    dTdz_bot = hostFitMassWeightedTemperatureDTdz(particles, kOphelieTemperatureField);
                    dTdz_top = dTdz_bot;
                }
            }
            Real side_w = 0, bottom_w = 0, free_c = 0, free_r = 0, total_w = 0;
            hostOphelieThermalFrenchNaturalHeatLossPowers(particles, face_bc, shell, side_w, bottom_w, free_c, free_r,
                                                          total_w, z_bot, z_bot + H, dTdz_bot, extrapolate_face,
                                                          dTdz_top);
            Real P_prod = total_w;
            if (std::strcmp(fc.face, "side") == 0)
                P_prod = side_w;
            else if (std::strcmp(fc.face, "bottom") == 0)
                P_prod = bottom_w;
            else if (std::strcmp(fc.face, "free_conv") == 0)
                P_prod = free_c;
            else
                P_prod = free_r;

            const Real P_ref = analytic(field, fc);
            const Real power_scale = std::max({std::abs(P_ref), std::abs(P_prod), Real(1)});
            const Real power_rel = std::abs(P_prod - P_ref) / power_scale;

            // Incremental Σ ρ0 V cp Δ(dT): matches Apply BC (inv_rho_cp_delta) and host power (V/δ).
            // Avoids float32 cancellation of two large absolute H sums; Mass≠ρ0V must not fail the gate.
            syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
            syncVariableToHost<Real>(particles, "VolumetricMeasure");
            const Real *dT_before = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
            const Real *vol_e = particles.getVariableDataByName<Real>("VolumetricMeasure");
            StdVec<Real> dT_pre(n);
            for (size_t i = 0; i < n; ++i)
                dT_pre[i] = dT_before[i];

            Real dTdz_bot_live = dTdz_bot;
            Real dTdz_top_live = dTdz_top;
            Real *dTdz_bot_ptr = extrapolate_face ? &dTdz_bot_live : nullptr;
            Real *dTdz_top_ptr = extrapolate_face ? &dTdz_top_live : nullptr;
            StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_face(
                glass_body, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, shell, face_bc, z_bot,
                z_bot + H, dTdz_bot_ptr, dTdz_top_ptr);
            apply_face.exec(dt);

            syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
            syncVariableToHost<Real>(particles, "VolumetricMeasure");
            const Real *dT_after = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
            vol_e = particles.getVariableDataByName<Real>("VolumetricMeasure");
            double dE = 0.0;
            const double rho_cp = static_cast<double>(rho0) * static_cast<double>(cp);
            for (size_t i = 0; i < n; ++i)
                dE += rho_cp * static_cast<double>(vol_e[i]) *
                      (static_cast<double>(dT_after[i]) - static_cast<double>(dT_pre[i]));
            const Real DeltaH = static_cast<Real>(dE);
            const Real E_ledger = -P_prod * dt; // outward loss cools melt
            const Real energy_scale = std::max({std::abs(E_ledger), std::abs(DeltaH), Real(1.0e-6)});
            const Real energy_rel = std::abs(DeltaH - E_ledger) / energy_scale;

            const bool pass_power = power_rel <= power_tol;
            const bool pass_energy = energy_rel <= energy_tol;
            all_ok = all_ok && pass_power && pass_energy;

            pow_csv << field << "," << fc.face << "," << P_prod << "," << P_ref << "," << power_rel << "," << DeltaH
                    << "," << E_ledger << "," << energy_rel << "," << (pass_power ? 1 : 0) << ","
                    << (pass_energy ? 1 : 0) << "," << fc.h_side << "," << fc.h_bottom << "," << fc.h_free << ","
                    << fc.emissivity << "," << dt << "\n";
            std::cout << "[ophelie][doc232][D] " << field << "/" << fc.face << " P_prod=" << P_prod
                      << " P_ref=" << P_ref << " power_rel=" << power_rel << " energy_rel=" << energy_rel
                      << " pass_p=" << (pass_power ? 1 : 0) << " pass_e=" << (pass_energy ? 1 : 0) << "\n";
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        T[i] = T_backup[i];
        dT[i] = dT_backup[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);

    std::cout << "[ophelie][doc232][D] PRODUCTION_BOUNDARY_FACE_AUDIT=" << (all_ok ? 1 : 0)
              << " boundary_face_t_mode=" << face_t_mode << " out_dir=" << out_dir << "\n";
    if (doc233_scope_fields)
    {
        std::cout << "[ophelie][doc233][D] BOUNDARY_RECONSTRUCTION_SCOPE="
                  << (all_ok ? "PASS_WITHIN_TOL" : "LIMITED_OR_FAIL")
                  << " note=band|global_dt_dz|local_dt_dz; parabolic_z+rz_quad stress face-slope limits; see CSV\n";
        std::cout << "[ophelie][doc235][F] rz_quad=T=1200+100*(1+0.5*(r/R)^2)*(2*(z-zb)/H-1)^2 "
                     "local_dt_dz=per-face-single-slope (not per-point 2D)\n";
    }
    return all_ok;
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

    if (const int cli_rc = applyLocalCli(ac, av, french, local_cli))
        return cli_rc;
    if (local_cli.client_nc_geometry)
    {
        OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
        if (std::isfinite(static_cast<double>(local_cli.client_melt_height_m)) &&
            local_cli.client_melt_height_m > TinyReal)
        {
            melt.melt_height_m = local_cli.client_melt_height_m;
            melt.glass_half_height = Real(0.5) * melt.melt_height_m;
            melt.glass_center_z_m = melt.inner_bottom_z_m + melt.glass_half_height;
            melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);
            melt.melt_height_status = "user_cli_override";
        }
        // Client development default 10 mm (matches relax). Later filterFrenchReduced applies --dp= if set.
        applyClientNcMeltToFrenchReduced(french, melt, Real(0.010));
        std::cout << "[ophelie][client-nc] geometry: R=" << french.glass_radius
                  << " H=" << (Real(2) * french.glass_half_height) << " center_z=" << french.glass_center[2]
                  << " bottom_z=" << (french.glass_center[2] - french.glass_half_height)
                  << " dp=" << french.dp
                  << " coil_loops=" << french.coil.num_loops << " f_Hz=" << french.frequency_hz
                  << " (avoid bare --glass-height= after this; use --melt-height=)\n";
    }
    if (local_cli.client_nc_material)
    {
        const OphelieClientCaseParams preset = makeClientNaturalConvectionDevelopmentParams();
        applyClientNcDevelopmentToLiterature(local_cli.literature, preset);
        local_cli.enable_target_power = false;
        local_cli.beta = preset.beta;
        french.sigma_glass = local_cli.literature.material.sigma.evaluate(preset.t_initial_k);
        french.frequency_hz = preset.frequency_hz;
        // Isothermal hydro short-run default: no EM Q unless user later enables heating.
        if (!local_cli.literature.zero_q && !local_cli.prescribed_heat_mode)
        {
            std::cout << "[ophelie][client-nc] tip: pass --zero-q for isothermal hydro (no EM Q / no 50 kW)\n";
        }
        std::cout << "[ophelie][client-nc] material: id=" << preset.material_id
                  << " T0_K=" << preset.t_initial_k
                  << " mu=" << local_cli.literature.material.mu.evaluate(preset.t_initial_k)
                  << " sigma=" << french.sigma_glass
                  << " rho=" << local_cli.literature.material.rho.evaluate(preset.t_initial_k)
                  << " boussinesq=" << (local_cli.literature.enable_boussinesq ? 1 : 0)
                  << " generator_input_W=" << local_cli.literature.generator_power_w << "\n";
    }
    // Doc245: material dataset AFTER case identity so μ/k/cp can be swapped without
    // rewriting 50 kW / β / ρ. Missing JAC thesis CSV → hard fail (no CES fallback).
    if (!local_cli.material_dataset.empty())
    {
        if (local_cli.client_nc_material)
        {
            std::cerr << "[ophelie][doc245] FATAL: --material-dataset= cannot combine with "
                         "--client-nc-material\n";
            return 2;
        }
        const Real power_for_check =
            local_cli.enable_target_power
                ? local_cli.target_power_w
                : local_cli.literature.target_glass_absorbed_power_w;
        const auto mat = applyDoc245MaterialDataset(local_cli.material_dataset, local_cli.literature.material,
                                                    power_for_check, local_cli.material_tables_dir,
                                                    /*apply_cep1473_anchors=*/false);
        if (!mat.ok)
            return 2;
        // Force particle-local thermal rheology when thesis/legacy tables are explicitly selected.
        local_cli.literature.constant_thermal_rheology = false;
        std::string ev_dir = local_cli.material_evidence_dir;
        if (ev_dir.empty() && !local_cli.summary_json.empty())
        {
            const auto p = fs::path(local_cli.summary_json).parent_path();
            if (!p.empty())
                ev_dir = p.string();
        }
        if (!ev_dir.empty())
        {
            fs::create_directories(ev_dir);
            OphelieFrenchGlassMaterialLaws legacy_snapshot;
            installJacoutotCes2008Fig2ThermalRheology(legacy_snapshot);
            writeDoc245MaterialEffectiveJson(ev_dir + "/material_effective.json", mat,
                                             local_cli.literature.material, local_cli.literature.t_initial_k);
            writeDoc245MaterialProbeCsv(ev_dir + "/material_probe.csv", mat.dataset_id,
                                        local_cli.literature.material, &legacy_snapshot);
            writeDoc245MaterialCurveCompareCsv(ev_dir + "/material_curve_compare.csv",
                                               local_cli.literature.material, legacy_snapshot);
            std::cout << "[ophelie][doc245][material] evidence_dir=" << ev_dir
                      << " wrote material_effective.json material_probe.csv material_curve_compare.csv\n";
        }
    }
    if (local_cli.prescribed_heat_mode)
    {
        // client-nc material zeroes target_glass_absorbed_power_w; restore prescribed identity.
        local_cli.enable_target_power = false;
        local_cli.literature.coupling = ThermoEMCouplingMode::Off;
        local_cli.literature.target_glass_absorbed_power_w = local_cli.prescribed_heat_power_w;
        local_cli.target_power_w = local_cli.prescribed_heat_power_w;
    }
    if (local_cli.doc232_h_cli_selftest_exit)
    {
        // After client-nc geometry/material so generator_power_w reflects the production 30 kW identity.
        const char *heat_mode = local_cli.prescribed_heat_mode
                                    ? "prescribed_uniform_volumetric"
                                    : (local_cli.literature.zero_q ? "zero_q"
                                                                  : (local_cli.enable_target_power
                                                                         ? "em_target_power"
                                                                         : "em_or_frozen_q"));
        const char *em_status = local_cli.prescribed_heat_mode
                                    ? "NOT_RUN"
                                    : (local_cli.literature.zero_q ? "DISABLED" : "RUN");
        std::cout << "[ophelie][doc232][H] cli_selftest heat_source_mode=" << heat_mode
                  << " em_status=" << em_status
                  << " prescribed_heat_power_w=" << local_cli.prescribed_heat_power_w
                  << " generator_power_w=" << local_cli.literature.generator_power_w
                  << " target_power_w=" << local_cli.target_power_w
                  << " enable_target_power=" << (local_cli.enable_target_power ? 1 : 0) << "\n";
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << std::setprecision(10);
            js << "{\n"
               << "  \"case\": \"doc232_h_cli_selftest\",\n"
               << "  \"heat_source_mode\": \"" << heat_mode << "\",\n"
               << "  \"em_status\": \"" << em_status << "\",\n"
               << "  \"prescribed_heat_power_w\": " << local_cli.prescribed_heat_power_w << ",\n"
               << "  \"generator_power_w\": " << local_cli.literature.generator_power_w << ",\n"
               << "  \"target_power_w\": " << local_cli.target_power_w << ",\n"
               << "  \"enable_target_power\": " << (local_cli.enable_target_power ? 1 : 0) << ",\n"
               << "  \"em_solver_ok\": null,\n"
               << "  \"SOFTWARE_DELTA_READY_cli\": 1,\n"
               << "  \"note\": \"CLI-only; EM not executed; em_solver_ok must not be read as convergence\"\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return 0;
    }
    if (local_cli.doc227_schedule_selftest)
    {
        std::string reason;
        const bool ok = hostDoc227ScheduleSelftest(reason);
        std::cout << "[ophelie][doc227][B5] schedule_selftest PASS=" << (ok ? 1 : 0)
                  << " reason=" << reason << std::endl;
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n  \"case\": \"doc227_schedule_selftest\",\n  \"passed\": " << (ok ? 1 : 0)
               << ",\n  \"reason\": \"" << reason << "\"\n}\n";
        }
        return ok ? 0 : 1;
    }
    std::cout << "[ophelie][cli] doc223_paper_flat_slip_top=" << (local_cli.doc223_paper_flat_slip_top ? 1 : 0)
              << " diag_wall_flat_lid=" << (local_cli.diag_wall_flat_lid ? 1 : 0)
              << " pprime=" << (local_cli.doc220_pprime_hydro ? 1 : 0)
              << " omit_g=" << (local_cli.omit_gravity_force ? 1 : 0)
              << " doc227_actual_timestep=" << (local_cli.doc227_actual_timestep ? 1 : 0) << "\n";
    // Doc219 D/F: cavity-validated B + thermal-after-update; TVF stays opt-in on French FS.
    // E N0/N1@1s with default TVF gave U_rms≈0.15 ≫ U_buoy≈0.012; Doc105 quiet floor was ~5e-3 without TVF.
    if (local_cli.doc219_hydro_profile)
    {
        if (!local_cli.kernel_correction_user_set)
            local_cli.kernel_correction = true;
        if (!local_cli.enable_tvf_user_set)
            local_cli.enable_tvf = false;
    }
    // Doc220 E: with p_h in BackgroundPressure and GravityForceCK still on (total-p momentum),
    // wall dens DissipativeUJump must use p' + a_ref=−g (AcousticDensMode=1). Mode 0 keeps
    // injecting −C_ext dens at walls even when hydrostatic pressure is already consistent.
    // Doc220 E2 p' package already sets dens_mode=1 with a_ref=0 (g not in ForcePrior).
    if (local_cli.background_hydro_p && !local_cli.acoustic_dens_mode_user_set &&
        !local_cli.omit_gravity_force)
    {
        local_cli.acoustic_dens_mode = 1;
        std::cout << "[ophelie][doc220][e] background_hydro_p => acoustic_dens_mode=1 "
                     "(wall dens on p', ReferenceBodyAcceleration=-g)\n";
    }
    if (local_cli.doc220_pprime_hydro)
    {
        std::cout << "[ophelie][doc220][e2] pprime_hydro=1 BackgroundPressure=p_h "
                     "PerturbationMomentum=1 omit_GravityForceCK=1 acoustic_dens_mode=1 a_ref=0 "
                     "(French free-surface; not closed-box P2b copy)\n";
    }
    if (local_cli.doc223_paper_flat_slip_top)
    {
        std::cout << "[ophelie][doc223][D] hydro_model=paper_flat_slip_top "
                     "pressure_model=perturbation_fixed_top dens=Internal "
                     "dens_sum_normalize=" << (local_cli.dens_sum_normalize ? 1 : 0)
                  << " settle_dens_sum_normalize="
                  << (local_cli.settle_dens_sum_normalize ? 1 : 0)
                  << " keep_dens_reinit=0"
                  << " kernel_correction=" << (local_cli.kernel_correction ? 1 : 0)
                  << (local_cli.doc223_ablate_no_free_slip
                          ? " lid=ClosedCrucible_noslip_ABLATION"
                          : (local_cli.doc224_ablate_acoustic_freeslip
                                 ? " ClosedCrucible_acoustic+OpenCup_viscous_O3"
                                 : " ClosedCrucible+P2_AcousticWallFreeSlip+OpenCup_viscous"))
                  << " acoustic_fs=" << (local_cli.doc224_ablate_acoustic_freeslip ? 0 : 1)
                  << " keep_dens_jump=" << (local_cli.doc224_lid_keep_dens_jump ? 1 : 0)
                  << " NO tip/zero-net/z-clamp/ContSilent (paper 2.1 slip top engineering model)\n";
        if (local_cli.doc223_ablate_no_free_slip)
        {
            std::cout << "[ophelie][doc223][E1] ablate_no_free_slip=1 (AcousticWallFreeSlip OFF; "
                         "not paper slip; isolates free-slip noise)\n";
        }
        else if (local_cli.doc224_ablate_acoustic_freeslip)
        {
            std::cout << "[ophelie][doc224][C3d] paper-top default: acoustic no-slip + OpenCup "
                         "viscous (impermeable slip; AcousticWallFreeSlip OFF)\n";
        }
        else if (local_cli.doc224_lid_keep_dens_jump)
        {
            std::cout << "[ophelie][doc224][C3] lid_keep_dens_jump=1 "
                         "(legacy acoustic free-slip mom∥n; DissipativeUJump dens ON)\n";
        }
    }
    if (local_cli.background_hydro_p && !local_cli.doc220_pprime_hydro && !local_cli.omit_gravity_force &&
        local_cli.perturbation_momentum == 0)
    {
        std::cout << "[ophelie][doc222][E2] totalp_bgp_hydro=1 BackgroundPressure=p_h "
                     "PerturbationMomentum=0 GravityForceCK=1 acoustic_dens_mode=1 a_ref=-g "
                     "(flat: ∇p_h↔g; deform restoring needs local_eta)\n";
    }
    if (local_cli.background_hydro_p || local_cli.doc220_pprime_hydro)
    {
        std::cout << "[ophelie][doc222][E1] BackgroundPressure p_h signed_depth="
                  << (local_cli.doc222_ph_signed ? 1 : 0)
                  << (local_cli.doc222_ph_signed
                          ? " (ρ0*g*(z_ref-z); restoring above ref plane)"
                          : " (LEGACY max(0,z_ref-z); D2 FAIL path)")
                  << "\n";
        if (local_cli.doc222_ph_local_eta)
        {
            std::cout << "[ophelie][doc222][E2b] BackgroundPressure local_eta=1 bins="
                      << local_cli.doc222_ph_local_eta_bins
                      << " (p_h=ρ0*g*(η(r)-z); η=radial bin z_max; fixed-z_free path OFF)\n";
        }
    }
    if (local_cli.diag_zero_net_uz)
    {
        std::cout << "[ophelie][doc221][K] diag_zero_net_uz=1 (mass-weighted mean Uz removed each "
                     "advection; volume-closure DIAGNOSTIC only — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_hard_lid)
    {
        std::cout << "[ophelie][doc221][L] diag_hard_lid=1 (kinematic impermeable lid at z_free; "
                     "DIAGNOSTIC only — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_soft_lid)
    {
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_tip_soft_force = false;
        std::cout << "[ophelie][doc221][N] diag_soft_lid=1 (zero-net-Uz + thin-band kill-up, NO z "
                     "clamp; DIAGNOSTIC only — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_soft_lid_tip)
    {
        local_cli.diag_soft_lid = false;
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_zero_net_uz = true;
        local_cli.diag_hard_lid = false;
        if (local_cli.diag_soft_lid_tip_r <= TinyReal)
            local_cli.diag_soft_lid_tip_r = Real(0.10);
        if (local_cli.diag_soft_lid_tip_uz_scale < Real(0))
            local_cli.diag_soft_lid_tip_uz_scale = Real(0);
        if (local_cli.diag_soft_lid_tip_uz_scale > Real(1))
            local_cli.diag_soft_lid_tip_uz_scale = Real(1);
        const bool dampen = local_cli.diag_soft_lid_tip_uz_scale > TinyReal;
        std::cout << "[ophelie][doc221][" << (dampen ? "N3" : "N2")
                  << "] diag_soft_lid_tip=1 (zero-net-Uz + tip-local "
                  << (dampen ? "dampen" : "kill-up") << " r<" << local_cli.diag_soft_lid_tip_r
                  << " m uz_scale=" << local_cli.diag_soft_lid_tip_uz_scale
                  << "; NO full-lid kill-up; NO z clamp; DIAGNOSTIC — NOT literature open "
                     "free-surface)\n";
    }
    if (local_cli.diag_tip_soft_force)
    {
        local_cli.diag_tip_fs_dashpot = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_zero_net_uz = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_wall_flat_lid = false;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        if (local_cli.diag_tip_soft_force_r <= TinyReal)
            local_cli.diag_tip_soft_force_r = Real(0.05);
        if (local_cli.diag_tip_soft_force_k < Real(0))
            local_cli.diag_tip_soft_force_k = Real(0);
        if (local_cli.diag_tip_soft_force_c < Real(0))
            local_cli.diag_tip_soft_force_c = Real(0);
        std::cout << "[ophelie][doc221][N4] diag_tip_soft_force=1 (zero-net-Uz + tip ForcePrior "
                     "dashpot"
                  << (local_cli.diag_tip_soft_force_k > TinyReal ? "+spring" : "") << " r<"
                  << local_cli.diag_tip_soft_force_r << " m k=" << local_cli.diag_tip_soft_force_k
                  << " c=" << local_cli.diag_tip_soft_force_c
                  << "; NO host Uz kill; NO z clamp; DIAGNOSTIC — NOT literature open "
                     "free-surface)\n";
    }
    if (local_cli.diag_tip_fs_dashpot)
    {
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_zero_net_uz = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_wall_flat_lid = false;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_tip_soft_force_k = Real(0);
        if (local_cli.diag_tip_soft_force_r <= TinyReal)
            local_cli.diag_tip_soft_force_r = Real(0.05);
        if (local_cli.diag_tip_soft_force_c < Real(0))
            local_cli.diag_tip_soft_force_c = Real(0);
        std::cout << "[ophelie][doc221][N5] diag_tip_fs_dashpot=1 (zero-net-Uz + tip ForcePrior "
                     "dashpot ONLY on DensitySummation<rho0 in tip r<"
                  << local_cli.diag_tip_soft_force_r << " m c=" << local_cli.diag_tip_soft_force_c
                  << "; NO spring; NO host Uz kill; NO z clamp; DIAGNOSTIC — NOT literature open "
                     "free-surface)\n";
    }
    if (local_cli.diag_lid_dens_silence)
    {
        local_cli.diag_lid_cont_silent = false;
        local_cli.diag_lid_cont_silent_only = false;
        local_cli.diag_lid_omit_pjump_only = false;
        local_cli.diag_wall_lid_viscous_omit = true;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.acoustic_dens_mode = 3;
        local_cli.acoustic_dens_mode_user_set = true;
        std::cout << "[ophelie][doc221][O4] diag_lid_dens_silence=1 (O3 skeleton + dens_mode=3 + "
                     "lid ContSilent/OmitPJump; NO AcousticWallFreeSlip; DIAGNOSTIC — NOT "
                     "literature open free-surface)\n";
    }
    if (local_cli.diag_lid_cont_silent)
    {
        local_cli.diag_lid_dens_silence = false;
        local_cli.diag_lid_cont_silent_only = false;
        local_cli.diag_lid_omit_pjump_only = false;
        local_cli.diag_wall_lid_viscous_omit = true;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        // Keep dens_mode=1 (do not override to 3).
        std::cout << "[ophelie][doc221][O5] diag_lid_cont_silent=1 (O3 skeleton + lid "
                     "ContSilent/OmitPJump ONLY; dens_mode stays 1; NO AcousticWallFreeSlip; "
                     "DIAGNOSTIC — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_lid_cont_silent_only)
    {
        local_cli.diag_lid_dens_silence = false;
        local_cli.diag_lid_cont_silent = false;
        local_cli.diag_lid_omit_pjump_only = false;
        local_cli.diag_wall_lid_viscous_omit = true;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        std::cout << "[ophelie][doc221][O6a] diag_lid_cont_silent_only=1 (O3 + ContSilent ONLY; "
                     "OmitPJump=0 dens_mode=1; NO AcousticWallFreeSlip; DIAGNOSTIC — NOT "
                     "literature open free-surface)\n";
    }
    if (local_cli.diag_lid_omit_pjump_only)
    {
        local_cli.diag_lid_dens_silence = false;
        local_cli.diag_lid_cont_silent = false;
        local_cli.diag_lid_cont_silent_only = false;
        local_cli.diag_wall_lid_viscous_omit = true;
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        std::cout << "[ophelie][doc221][O6b] diag_lid_omit_pjump_only=1 (O3 + OmitPJump ONLY; "
                     "ContSilent=0 dens_mode=1; NO AcousticWallFreeSlip; DIAGNOSTIC — NOT "
                     "literature open free-surface)\n";
    }
    if (local_cli.diag_wall_lid_viscous_omit)
    {
        // O3 ablation must not combine with P acoustic free-slip.
        local_cli.diag_lid_free_slip = false;
        local_cli.diag_lid_free_slip_post_settle = false;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_tip_fs_dashpot = false;
        std::cout << "[ophelie][doc221][O3] diag_wall_lid_viscous_omit=1 (ClosedCrucible acoustic "
                     "NO AcousticWallFreeSlip; viscous OpenCup omit lid; DIAGNOSTIC — NOT "
                     "literature open free-surface)\n";
    }
    if (local_cli.diag_lid_free_slip_post_settle)
    {
        local_cli.diag_lid_free_slip = true;
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_tip_fs_dashpot = false;
        std::cout << "[ophelie][doc221][P2] diag_lid_free_slip_post_settle=1 (settle=O no-slip lid; "
                     "after settle zero-U + arm AcousticWallFreeSlip + viscous OpenCup; "
                     "DIAGNOSTIC — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_lid_free_slip)
    {
        local_cli.diag_wall_flat_lid = true;
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_tip_fs_dashpot = false;
        std::cout << "[ophelie][doc221][P] diag_lid_free_slip=1 (ClosedCrucible acoustic + "
                     "AcousticWallFreeSlip lid; viscous open-cup omit lid"
                  << (local_cli.diag_lid_free_slip_post_settle ? "; arm AFTER settle" : "; arm from t=0")
                  << "; DIAGNOSTIC — NOT literature open free-surface)\n";
    }
    if (local_cli.diag_wall_flat_lid)
    {
        // Physical lid supersedes kinematic soft/hard lids.
        local_cli.diag_hard_lid = false;
        local_cli.diag_soft_lid = false;
        local_cli.diag_soft_lid_tip = false;
        local_cli.diag_tip_soft_force = false;
        local_cli.diag_tip_fs_dashpot = false;
        std::cout << "[ophelie][doc221][O] diag_wall_flat_lid=1 (ClosedCrucible wall incl. flat lid; "
                     "no kinematic clamp/kill-up; DIAGNOSTIC — NOT literature open free-surface)\n";
    }
    std::cout << "[ophelie][doc219][hydro] profile=" << (local_cli.doc219_hydro_profile ? 1 : 0)
              << " kernel_correction=" << (local_cli.kernel_correction ? 1 : 0)
              << " tvf=" << (local_cli.enable_tvf ? 1 : 0)
              << " (French: Riemann+B+thermal-after-pos; TVF opt-in only — cavity closed-box TVF not default on FS)"
              << std::endl;
    std::cout << "[ophelie][doc219][spatial-q] enabled=" << (local_cli.spatial_frozen_q ? 1 : 0)
              << " doc217_d2=" << (local_cli.doc217_d2_local_poly ? 1 : 0)
              << " self_induction=" << (local_cli.enable_self_induction ? 1 : 0)
              << " (Q_ref(x) Euler resample each thermal; not particle-ID advection)" << std::endl;
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
    if (local_cli.doc217_d2_local_poly)
    {
        params.doc217_d2_local_poly_ = true;
        params.doc216_e5_p1_projection_ = false;
        params.doc217_d2_radius_ = french.glass_radius;
        params.doc217_d2_half_height_ = french.glass_half_height;
        params.doc217_d2_center_ = french.glass_center;
        params.doc217_d2_dp_ = french.dp;
        std::cout << "[ophelie][doc217][d2] frozen_q enable local quadratic Neumann "
                     "(host SparseLU; opt-in)"
                  << std::endl;
    }
    if (local_cli.enable_target_power)
    {
        // Must update french.target_joule_power too: syncFrenchReducedToParameters copies
        // french → params and would otherwise wipe CLI P35/P65 back to the 50 kW default.
        french.target_joule_power = local_cli.target_power_w;
        params.target_joule_power_ = local_cli.target_power_w;
        params.enable_power_scaling_ = false; // French NC uses coil-I calibration, not field scaling
    }
    if (params.self_induction_phi_eq_res_tolerance_ > Real(1.0e-3))
    {
        params.self_induction_phi_eq_res_tolerance_ = Real(2.0e-4);
    }
    syncFrenchReducedToParameters(french, params);
    if (local_cli.enable_target_power)
    {
        // Re-assert after sync (and before calibrate) so FixedAbsorbedPower hits CLI watts.
        french.target_joule_power = local_cli.target_power_w;
        params.target_joule_power_ = local_cli.target_power_w;
    }
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

        BaseParticles &particles = glass_em->getBaseParticles();
        const size_t n = particles.TotalRealParticles();
        frozen_q.assign(n, Real(0));

        if (local_cli.literature.zero_q)
        {
            // Isothermal hydro: keep EM SPHSystem alive for GPU lifetime, but do not solve EM.
            p_joule = Real(0);
            phi_eq = Real(0);
            coupling_state.calibrated_current_per_loop = french.coil.current_per_loop;
            coupling_state.coil_current_per_loop = french.coil.current_per_loop;
            coupling_state.current_frozen = true;
            coupling_state.generator_power_w = local_cli.literature.generator_power_w;
            coupling_state.target_glass_absorbed_power_w = Real(0);
            coupling_state.reconstructed_glass_power_w = Real(0);
            coupling_state.last_phi_eq_res = Real(0);
            std::cout << "[ophelie][stage4.1] zero_q: skip EM; frozen_q=0 n=" << n
                      << " (isothermal hydro path)\n";
        }
        else if (local_cli.prescribed_heat_mode)
        {
            // Doc232 H / Doc231 C1: uniform volumetric absorbed heat; EM NOT_RUN (current unknown).
            // Armed by --prescribed-heat-power= presence; P=0 is valid (Q≡0) and must not fall back to EM.
            syncVariableToHost<Real>(particles, "VolumetricMeasure");
            const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
            Real V = Real(0);
            for (size_t i = 0; i < n; ++i)
                V += vol[i];
            const Real q_uni = local_cli.prescribed_heat_power_w / std::max(V, TinyReal);
            const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_em_names, params);
            syncVariableToHost<Real>(particles, q_field);
            Real *q = particles.getVariableDataByName<Real>(q_field);
            for (size_t i = 0; i < n; ++i)
            {
                q[i] = q_uni;
                frozen_q[i] = q_uni;
            }
            syncVariableToDevice<Real>(particles, q_field);
            p_joule = local_cli.prescribed_heat_power_w;
            phi_eq = Real(0);
            coupling_state.calibrated_current_per_loop = french.coil.current_per_loop;
            coupling_state.coil_current_per_loop = french.coil.current_per_loop;
            coupling_state.current_frozen = true;
            coupling_state.generator_power_w = local_cli.literature.generator_power_w;
            coupling_state.target_glass_absorbed_power_w = local_cli.prescribed_heat_power_w;
            coupling_state.reconstructed_glass_power_w = p_joule;
            coupling_state.last_phi_eq_res = Real(0);
            std::cout << "[ophelie][doc232][H] prescribed_uniform_q: em_status=NOT_RUN; n=" << n << " V=" << V
                      << " Q=" << q_uni << " P_abs=" << p_joule
                      << " W generator_power_w=" << coupling_state.generator_power_w << "\n";

            const bool build_spatial_q = local_cli.spatial_frozen_q || periodic_coupling;
            if (build_spatial_q)
            {
                const Vecd &gc_em = french.glass_center;
                const FrenchCylindricalFrame frame = makeFrenchCylinderFrameFromGlass(
                    gc_em, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
                const BoundingBoxd melt_bounds(
                    Vecd(gc_em[0] - french.glass_radius, gc_em[1] - french.glass_radius,
                         gc_em[2] - french.glass_half_height),
                    Vecd(gc_em[0] + french.glass_radius, gc_em[1] + french.glass_radius,
                         gc_em[2] + french.glass_half_height));
                q_grid.spec_ = rh200::makeRh200JouleHeatGridSpecFromBounds(melt_bounds, french.dp, Real(8));
                coupling_state.sigma_rz.frame = frame;
                coupling_state.q_rz.frame = frame;
                FrenchLiteratureParameters lit_q = local_cli.literature;
                lit_q.q_conservative_remap = true;
                depositEmQAndRecord(particles, glass_em_names, params, frame, lit_q, q_grid, coupling_state,
                                    p_joule);
                std::cout << "[ophelie][doc231] prescribed_uniform_q spatial grid built"
                          << " spatial_frozen_q=" << (local_cli.spatial_frozen_q ? 1 : 0) << "\n";
            }
        }
        else
        {
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

            const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_em_names, params);
            syncVariableToHost<Real>(particles, q_field);
            const Real *q = particles.getVariableDataByName<Real>(q_field);
            for (size_t i = 0; i < n; ++i)
            {
                frozen_q[i] = q[i];
            }
            std::cout << "[ophelie][stage4.1] EM freeze: n=" << n << " P_joule_W=" << p_joule
                      << " phi_eq_res_vol=" << phi_eq
                      << " self_induction=" << (probe.used_self_induction ? 1 : 0)
                      << " generator_power=" << coupling_state.generator_power_w
                      << " target_glass_absorbed_power=" << coupling_state.target_glass_absorbed_power_w
                      << " reconstructed_glass_power=" << p_joule
                      << " coil_current=" << coupling_state.coil_current_per_loop << std::endl;

            // Doc219 E: always deposit frozen Q to Euler grid when spatial_frozen_q (default).
            const bool build_spatial_q = local_cli.spatial_frozen_q || periodic_coupling;
            if (build_spatial_q)
            {
                const Vecd &gc_em = french.glass_center;
                const FrenchCylindricalFrame frame = makeFrenchCylinderFrameFromGlass(
                    gc_em, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
                const BoundingBoxd melt_bounds(
                    Vecd(gc_em[0] - french.glass_radius, gc_em[1] - french.glass_radius,
                         gc_em[2] - french.glass_half_height),
                    Vecd(gc_em[0] + french.glass_radius, gc_em[1] + french.glass_radius,
                         gc_em[2] + french.glass_half_height));
                // Free-surface particles can leave the geometric melt AABB; keep a wide Euler pad
                // and clamp-sample (below) so Doc219 E never silently writes Q=0 on coverage miss.
                q_grid.spec_ = rh200::makeRh200JouleHeatGridSpecFromBounds(melt_bounds, french.dp, Real(8));
                coupling_state.sigma_rz.frame = frame;
                coupling_state.q_rz.frame = frame;
                // Fixed-shape / fixed-total-power diagnostic remap (not live EM feedback).
                FrenchLiteratureParameters lit_q = local_cli.literature;
                lit_q.q_conservative_remap = true;
                depositEmQAndRecord(particles, glass_em_names, params, frame, lit_q, q_grid, coupling_state,
                                    p_joule);
                if (periodic_coupling)
                {
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
                }
                writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
                writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
                appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                          local_cli.literature.em_control);
                writeFrenchEnergyBudgetCsvHeader("output/french_energy_budget.csv");
                std::cout << "[ophelie][doc219][spatial-q] euler_q_grid built; periodic="
                          << (periodic_coupling ? 1 : 0) << " spatial_frozen_q="
                          << (local_cli.spatial_frozen_q ? 1 : 0)
                          << " mode=" << emControlModeName(local_cli.literature.em_control) << std::endl;
            }
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
    // Doc221 O: closed flat-lid wall seals at melt top — tall splash rim unused for geometry.
    const Real rim_height = local_cli.diag_wall_flat_lid ? Real(0) : local_cli.wall_rim_height;
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
    const bool use_particle_mu = local_cli.doc226_local_particle_viscosity;
    if (use_particle_mu)
        glass.addMaterialProperty<ParticleViscosity>(mu);
    else
        glass.addMaterialProperty<Viscosity>(mu);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());
    Viscosity &glass_viscosity = glass.getMaterialProperty<Viscosity>();
    ParticleViscosity *glass_particle_viscosity =
        use_particle_mu ? &glass.getMaterialProperty<ParticleViscosity>() : nullptr;

    // SPHinXsys cylinder wall path: TriangleMeshShapeCylinder + LevelSet + NormalFromBodyShapeCK
    // (same pattern as taylor_bar_sycl / french_natural_glass_relax).
    // Doc221 O: optional ClosedCrucible (side+bottom+flat lid) instead of open cup + rim.
    // Doc221 P: + open-cup viscous partner (rim=0) so lid has free-slip (no viscous contact).
    const int wall_mesh_resolution =
        french.glass_mesh_resolution > 0 ? french.glass_mesh_resolution : 20;
    SharedPtr<Shape> wall_shape;
    if (local_cli.diag_wall_flat_lid)
    {
        wall_shape = makeShared<OphelieFrenchClosedCrucibleWallShape>("WallBoundary", french, wall_thickness,
                                                                     wall_mesh_resolution);
    }
    else
    {
        wall_shape = makeShared<OphelieFrenchNaturalCrucibleWallShape>(
            "WallBoundary", french, wall_thickness, wall_mesh_resolution, rim_height);
    }
    SolidBody wall(sph_system, wall_shape);
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    (void)wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall.generateParticles<BaseParticles, Lattice>();
    if (local_cli.diag_wall_flat_lid)
    {
        std::cout << "[ophelie][doc221][O] wall=ClosedCrucible n_wall="
                  << wall.getBaseParticles().TotalRealParticles() << " wall_thickness=" << wall_thickness
                  << " (flat lid at melt top)\n";
    }

    // Doc221 P/O3: viscous open-cup (side+bottom only). Lifetime must cover the whole run.
    const bool need_wall_noslip =
        local_cli.diag_lid_free_slip || local_cli.diag_wall_lid_viscous_omit;
    UniquePtr<SolidBody> wall_noslip_body;
    if (need_wall_noslip)
    {
        wall_noslip_body = makeUnique<SolidBody>(
            sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                            "WallNoSlip", french, wall_thickness, wall_mesh_resolution, Real(0)));
        wall_noslip_body->defineAdaptationRatios(1.3, 1.0);
        wall_noslip_body->defineMatterMaterial<Solid>();
        (void)wall_noslip_body->defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
        wall_noslip_body->generateParticles<BaseParticles, Lattice>();
        std::cout << "[ophelie][doc221][" << (local_cli.diag_wall_lid_viscous_omit ? "O3" : "P")
                  << "] wall_noslip=OpenCup n_wall="
                  << wall_noslip_body->getBaseParticles().TotalRealParticles()
                  << " (viscous; lid omitted → free-slip top)\n";
    }

    Inner<> glass_inner(glass);
    Contact<> glass_wall_contact(glass, {&wall});
    UniquePtr<Contact<>> glass_wall_noslip_contact;
    if (need_wall_noslip)
    {
        // Brace-init {ptr} does not forward through makeUnique; use an explicit StdVec.
        StdVec<RealBody *> noslip_targets{wall_noslip_body.get()};
        glass_wall_noslip_contact = makeUnique<Contact<>>(glass, noslip_targets);
    }
    // Acoustic/dens always on closed (or open) wall.
    // P1 / O3: viscous OpenCup from t=0.
    // P2 (post-settle): viscous ClosedCrucible during settle; OpenCup after arming.
    const bool viscous_open_cup_from_t0 =
        local_cli.diag_wall_lid_viscous_omit ||
        (local_cli.diag_lid_free_slip && !local_cli.diag_lid_free_slip_post_settle);
    Contact<> &viscous_wall_contact =
        viscous_open_cup_from_t0 ? *glass_wall_noslip_contact : glass_wall_contact;
    Inner<> wall_inner(wall);
    Contact<> wall_glass_contact(wall, {&glass});

    BaseParticles &glass_particles = glass.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(glass_particles, t0);
    registerOphelieThermalDiffusionAuxFields(glass_particles, k_th);
    registerOphelieFrenchThermalRheologyFields(glass_particles, mu, cp, k_th);
    // Doc226: μ/cp must follow ParticleSortCK with positions when local viscosity is on.
    glass_particles.addEvolvingVariable<Real>(kOphelieDynamicViscosityField);
    glass_particles.addEvolvingVariable<Real>(kOphelieHeatCapacityField);
    if (glass_particle_viscosity)
    {
        glass_particle_viscosity->bindParticleViscosityField(glass_particles, kOphelieDynamicViscosityField);
        std::cout << "[ophelie][doc226] ParticleViscosity bound to " << kOphelieDynamicViscosityField
                  << " (pairwise/one-side from particle field)\n";
    }
    glass_particles.registerStateVariable<Real>("DensSumSigma0");
    glass_particles.addEvolvingVariable<Real>("DensSumSigma0");
    glass_particles.registerStateVariable<Real>("DensSumZ0");
    glass_particles.addEvolvingVariable<Real>("DensSumZ0");
    // Doc222 B: tip ForcePrior channel (ForcePriorCK-style current/previous); sorted with particles.
    glass_particles.registerStateVariable<Vecd>("TipSoftForce");
    glass_particles.addEvolvingVariable<Vecd>("TipSoftForce");
    glass_particles.registerStateVariable<Vecd>("PreviousTipSoftForce");
    glass_particles.addEvolvingVariable<Vecd>("PreviousTipSoftForce");
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
    // Doc246 Q1: raw EM deposit kept separately when upper-kappa weighting is enabled.
    if (local_cli.doc246_q_upper_kappa > TinyReal)
        glass_particles.registerStateVariable<Real>("JouleHeatEM", Real(0));
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
        // Doc230 isothermal hydro (--zero-q): Q≡0 and no Euler Q grid — skip Doc219 spatial-q gates.
        if (local_cli.literature.zero_q)
        {
            std::cout << "[ophelie][doc219][spatial-q] skip self_check (zero_q isothermal; Q=0 by design)\n";
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
        }
        else if (local_cli.spatial_frozen_q || periodic_coupling)
        {
            // Doc219 E: initial map is Eulerian; self-check vs particle-ID snapshot at same layout.
            size_t n_uncovered = 0;
            const Real p_raw =
                hostSampleSpatialFrozenQ(glass_particles, q_grid, p_joule, /*renormalize=*/false, n_uncovered);
            syncVariableToHost<Real>(glass_particles, "JouleHeat");
            syncVariableToHost<Real>(glass_particles, "VolumetricMeasure");
            const Real *q_s = glass_particles.getVariableDataByName<Real>("JouleHeat");
            const Real *vol_s = glass_particles.getVariableDataByName<Real>("VolumetricMeasure");
            Real num = Real(0), den = Real(0), p_id = Real(0);
            for (size_t i = 0; i < n; ++i)
            {
                const Real dq = q_s[i] - frozen_q[i];
                num += dq * dq * vol_s[i];
                den += frozen_q[i] * frozen_q[i] * vol_s[i];
                p_id += frozen_q[i] * vol_s[i];
            }
            const Real q_rel = (den > TinyReal) ? std::sqrt(num / den) : Real(0);
            const Real p_raw_rel =
                (std::abs(p_joule) > TinyReal) ? std::abs(p_raw - p_joule) / std::abs(p_joule) : Real(0);
            const Real p_id_rel =
                (std::abs(p_joule) > TinyReal) ? std::abs(p_id - p_joule) / std::abs(p_joule) : Real(0);
            // Hard gates: coverage + power after optional diagnostic renormalize.
            const Real p_final =
                hostSampleSpatialFrozenQ(glass_particles, q_grid, p_joule, /*renormalize=*/true, n_uncovered);
            const Real p_final_rel =
                (std::abs(p_joule) > TinyReal) ? std::abs(p_final - p_joule) / std::abs(p_joule) : Real(0);
            const bool q_ok = q_rel <= Real(0.01);
            const bool p_ok = p_final_rel <= Real(0.01) && n_uncovered == 0;
            std::cout << "[ophelie][doc219][spatial-q] self_check q_L2_rel=" << q_rel
                      << " P_raw_rel=" << p_raw_rel << " P_final_rel=" << p_final_rel
                      << " P_particleID_rel=" << p_id_rel << " n_uncovered=" << n_uncovered << "/" << n
                      << " gate_q<=0.01 gate_P<=0.01"
                      << " passed=" << ((q_ok && p_ok) ? 1 : 0) << std::endl;
            if (n_uncovered > 0)
            {
                std::cerr << "[ophelie][doc219][spatial-q] abort: uncovered particles get silent Q=0\n";
                return 1;
            }
            if (!p_ok)
            {
                std::cerr << "[ophelie][doc219][spatial-q] abort: power after remap exceeds 1% gate\n";
                return 1;
            }
            if (!q_ok)
            {
                // Deposit/sample smoothing can exceed 1% vs particle-ID; keep going but flag.
                std::cout << "[ophelie][doc219][spatial-q] WARNING q_L2_rel>" << Real(0.01)
                          << " (Euler deposit/sample vs particle-ID; power gate still held)\n";
            }
            // Doc246 Q1_UPPER: reshape fluid Q after EM/Euler sample (keeps JouleHeatEM = raw).
            if (local_cli.doc246_q_upper_kappa > TinyReal)
            {
                const Real melt_h = Real(2) * french.glass_half_height;
                const Real z_bot = french.glass_center[2] - french.glass_half_height;
                applyDoc246QUpperControlled(glass_particles, local_cli.doc246_q_upper_kappa, z_bot, melt_h,
                                            local_cli.target_power_w);
            }
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
            if (periodic_coupling)
            {
                writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
                appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                          local_cli.literature.em_control);
            }
        }
        else
        {
            if (local_cli.doc246_q_upper_kappa > TinyReal)
            {
                const Real melt_h = Real(2) * french.glass_half_height;
                const Real z_bot = french.glass_center[2] - french.glass_half_height;
                applyDoc246QUpperControlled(glass_particles, local_cli.doc246_q_upper_kappa, z_bot, melt_h,
                                            local_cli.target_power_w);
            }
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
        }
    }

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = !local_cli.doc232_adiabatic_heat;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = local_cli.literature.enable_thermal_diffusion;
    thermal_bc.boundary_width_factor = params.phi_boundary_distance_factor_;
    thermal_bc.renormalize_boundary_area_to_geometry = local_cli.boundary_area_normalize;
    thermal_bc.h_side = local_cli.doc232_adiabatic_heat ? Real(0) : Real(300);
    thermal_bc.h_bottom = local_cli.doc232_adiabatic_heat ? Real(0) : Real(35);
    thermal_bc.h_free = local_cli.doc232_adiabatic_heat ? Real(0) : Real(20);
    thermal_bc.emissivity = local_cli.doc232_adiabatic_heat ? Real(0) : Real(0.8);
    thermal_bc.t_cool = Real(300);
    thermal_bc.t_ambient = Real(300);
    thermal_bc.t_rad_ambient = Real(300);
    // Desktop --thermal-bc= (keeps EM/flow). Doc232 C gate still wins when armed.
    if (!local_cli.doc232_adiabatic_heat)
    {
        auto norm_bc = [](std::string s) {
            for (char &c : s)
                if (c == '_')
                    c = '-';
            return s;
        };
        const std::string bc_mode = norm_bc(local_cli.thermal_bc_mode);
        if (bc_mode == "adiabatic" || bc_mode == "adiabatic-walls")
        {
            thermal_bc.enable_french_natural_bc = false;
            thermal_bc.enable_cold_wall_dirichlet = false;
            thermal_bc.h_side = Real(0);
            thermal_bc.h_bottom = Real(0);
            thermal_bc.h_free = Real(0);
            thermal_bc.emissivity = Real(0);
            std::cout << "[ophelie][thermal-bc] mode=adiabatic (no Robin/radiation; EM/diffusion kept)\n";
        }
        else if (bc_mode == "isothermal-wall" || bc_mode == "isothermal" || bc_mode == "dirichlet-wall")
        {
            thermal_bc.enable_french_natural_bc = false;
            thermal_bc.enable_cold_wall_dirichlet = true;
            thermal_bc.h_side = Real(0);
            thermal_bc.h_bottom = Real(0);
            thermal_bc.h_free = Real(0);
            thermal_bc.emissivity = Real(0);
            std::cout << "[ophelie][thermal-bc] mode=isothermal-wall (Dirichlet boundary T→T0)\n";
        }
        else
        {
            thermal_bc.enable_french_natural_bc = true;
            thermal_bc.enable_cold_wall_dirichlet = false;
            if (local_cli.thermal_h_side_set)
                thermal_bc.h_side = local_cli.thermal_h_side;
            if (local_cli.thermal_h_bottom_set)
                thermal_bc.h_bottom = local_cli.thermal_h_bottom;
            if (local_cli.thermal_h_free_set)
                thermal_bc.h_free = local_cli.thermal_h_free;
            if (local_cli.thermal_emissivity_set)
                thermal_bc.emissivity = local_cli.thermal_emissivity;
            if (local_cli.thermal_t_cool_set)
            {
                thermal_bc.t_cool = local_cli.thermal_t_cool;
                thermal_bc.t_ambient = local_cli.thermal_t_cool;
                thermal_bc.t_rad_ambient = local_cli.thermal_t_cool;
            }
            std::cout << "[ophelie][thermal-bc] mode=french-natural"
                      << " h_side/bottom/free=" << thermal_bc.h_side << "/" << thermal_bc.h_bottom << "/"
                      << thermal_bc.h_free << " eps=" << thermal_bc.emissivity
                      << " Tcool=" << thermal_bc.t_cool << "\n";
        }
    }
    if (local_cli.doc232_adiabatic_heat)
    {
        std::cout << "[ophelie][doc232][C] french natural BC disabled (adiabatic production heat)\n";
    }
    if (!local_cli.doc232_adiabatic_heat)
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
    if (need_wall_noslip)
        host_methods.addStateDynamics<NormalFromBodyShapeCK>(*wall_noslip_body).exec();

    auto &update_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(glass);
    auto &update_wall_cell_linked_list = main_methods.addParticleDynamicsGroup();
    update_wall_cell_linked_list.add(&main_methods.addCellLinkedListDynamics(wall));
    if (need_wall_noslip)
        update_wall_cell_linked_list.add(&main_methods.addCellLinkedListDynamics(*wall_noslip_body));
    auto &update_glass_relations = main_methods.addParticleDynamicsGroup();
    if (need_wall_noslip)
    {
        update_glass_relations.add(
            &main_methods.addRelationDynamics(glass_inner, glass_wall_contact, *glass_wall_noslip_contact));
    }
    else
    {
        update_glass_relations.add(&main_methods.addRelationDynamics(glass_inner, glass_wall_contact));
    }
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
    InteractionDynamicsCK<MainExecutionPolicy, KernelGradientIntegralCorrectedComplex>
        kernel_grad_integral(glass_inner, glass_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionCK<SPHBody, TruncatedLinear>>
        transport_correction(glass); // Doc210/219 coefficient 0.2
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
    struct Doc223DensExec
    {
        std::function<void()> fn;
        void exec() { fn(); }
    } density_regularization;
    if (local_cli.doc223_paper_flat_slip_top)
    {
        auto &dr =
            main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
                .addPostContactInteraction(glass_wall_contact)
                .addPostStateDynamics<fluid_dynamics::DensityRegularization, Internal>(glass);
        density_regularization.fn = [&dr]() { dr.exec(); };
    }
    else
    {
        auto &dr =
            main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
                .addPostContactInteraction(glass_wall_contact)
                .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(glass);
        density_regularization.fn = [&dr]() { dr.exec(); };
    }

    // Doc221 P: viscous wall = open-cup (omit lid) when free-slip; else same as acoustic wall.
    // Doc221 P2: closed during settle, OpenCup after arming (no post-contact chain).
    // Doc226: optional ParticleViscosity path (pairwise field); default Viscosity unchanged.
    using ViscousForceLocalInnerCK =
        fluid_dynamics::ViscousForceCK<Inner<WithUpdate, ParticleViscosity, NoKernelCorrectionCK>>;
    using ViscousForceLocalWallCK =
        fluid_dynamics::ViscousForceCK<Contact<Wall, ParticleViscosity, NoKernelCorrectionCK>>;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalInnerCK>> fluid_viscous_force_local;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalWallCK>>
        fluid_viscous_force_wall_closed_local;
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalWallCK>>
        fluid_viscous_force_wall_opencup_local;

    bool lid_free_slip_armed =
        local_cli.diag_lid_free_slip && !local_cli.diag_lid_free_slip_post_settle;
    std::function<void()> exec_viscous;
    std::function<void()> exec_viscous_wall_closed;
    std::function<void()> exec_viscous_wall_opencup;

    UniquePtr<OphelieCallFnDynamicsCK> viscous_wall_post;
    if (use_particle_mu)
    {
        // Couette isomorphism: wall must be post-contact so ForcePriorCK update sees wall shear.
        // Prior F1 path ran wall AFTER Inner WithUpdate → wall never entered ForcePrior.
        fluid_viscous_force_local =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalInnerCK>>(glass_inner);
        fluid_viscous_force_wall_closed_local =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalWallCK>>(
                viscous_wall_contact);
        exec_viscous_wall_closed = [&]() { fluid_viscous_force_wall_closed_local->exec(); };
        if (local_cli.diag_lid_free_slip_post_settle)
        {
            fluid_viscous_force_wall_opencup_local =
                makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceLocalWallCK>>(
                    *glass_wall_noslip_contact);
            exec_viscous_wall_opencup = [&]() { fluid_viscous_force_wall_opencup_local->exec(); };
            viscous_wall_post = makeUnique<OphelieCallFnDynamicsCK>([&]() {
                if (lid_free_slip_armed)
                    fluid_viscous_force_wall_opencup_local->exec();
                else
                    fluid_viscous_force_wall_closed_local->exec();
            });
            fluid_viscous_force_local->addPostContactInteraction(*viscous_wall_post);
            exec_viscous = [&]() { fluid_viscous_force_local->exec(); };
        }
        else
        {
            fluid_viscous_force_local->addPostContactInteraction(*fluid_viscous_force_wall_closed_local);
            exec_viscous = [&]() { fluid_viscous_force_local->exec(); };
            if (local_cli.diag_lid_free_slip || local_cli.diag_wall_lid_viscous_omit)
                exec_viscous_wall_opencup = exec_viscous_wall_closed;
        }
        std::cout << "[ophelie][doc226] ViscousForceCK templated on ParticleViscosity "
                     "(post-contact wall → ForcePrior)\n";
    }
    else
    {
        // Capture raw pointers into main_methods-owned dynamics (local refs would dangle).
        auto *fluid_viscous_force = &main_methods.addInteractionDynamicsWithUpdate<
            fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(glass_inner);
        auto *fluid_viscous_force_from_wall = &main_methods.addInteractionDynamics<
            fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(viscous_wall_contact);
        exec_viscous_wall_closed = [fluid_viscous_force_from_wall]() {
            fluid_viscous_force_from_wall->exec();
        };
        if (local_cli.diag_lid_free_slip_post_settle)
        {
            // Doc227 C2: same post-contact wall → ForcePrior as ParticleViscosity path.
            // Prior paper-top Viscosity path ran wall AFTER Inner WithUpdate → wall shear
            // never entered ForcePrior; C2a/C2b then diverged ~25% despite C1 force match.
            auto *noslip_visc = &main_methods.addInteractionDynamics<
                fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(
                *glass_wall_noslip_contact);
            exec_viscous_wall_opencup = [noslip_visc]() { noslip_visc->exec(); };
            viscous_wall_post = makeUnique<OphelieCallFnDynamicsCK>(
                [noslip_visc, fluid_viscous_force_from_wall, &lid_free_slip_armed]() {
                    if (lid_free_slip_armed)
                        noslip_visc->exec();
                    else
                        fluid_viscous_force_from_wall->exec();
                });
            fluid_viscous_force->addPostContactInteraction(*viscous_wall_post);
            exec_viscous = [fluid_viscous_force]() { fluid_viscous_force->exec(); };
            std::cout << "[ophelie][doc227][C2] Viscosity wall post-contact → ForcePrior "
                         "(aligned with ParticleViscosity)\n";
        }
        else
        {
            fluid_viscous_force->addPostContactInteraction(*fluid_viscous_force_from_wall);
            exec_viscous = [fluid_viscous_force]() { fluid_viscous_force->exec(); };
            if (local_cli.diag_lid_free_slip || local_cli.diag_wall_lid_viscous_omit)
                exec_viscous_wall_opencup = exec_viscous_wall_closed;
        }
    }

    // Doc221 P1: tag lid AFTER AcousticStep Contact registers AcousticWallFreeSlip (from t=0).
    // Doc221 P2: defer tagging until after settle + zero-U.
    if (local_cli.diag_lid_free_slip && !local_cli.diag_lid_free_slip_post_settle)
    {
        const Real z_lid = french.glass_center[2] + french.glass_half_height;
        const size_t n_lid = hostArmDoc224FreeSlipLid(wall.getBaseParticles(), z_lid, wall_thickness,
                                                      local_cli.doc224_lid_keep_dens_jump);
        std::cout << "[ophelie][doc221][P] AcousticWallFreeSlip lid tags=" << n_lid
                  << " / wall_closed_n=" << wall.getBaseParticles().TotalRealParticles()
                  << " z_lid=" << z_lid << " (armed from t=0)"
                  << " keep_dens_jump=" << (local_cli.doc224_lid_keep_dens_jump ? 1 : 0) << "\n";
    }
    else if (local_cli.diag_lid_free_slip_post_settle)
    {
        std::cout << "[ophelie][doc221][P2] AcousticWallFreeSlip deferred until post-settle arming\n";
    }
    else if (local_cli.diag_wall_lid_viscous_omit)
    {
        std::cout << "[ophelie][doc221][O3] AcousticWallFreeSlip=0 (viscous OpenCup only; acoustic "
                     "ClosedCrucible stays no-slip dens/mom)\n";
    }
    if (local_cli.diag_lid_dens_silence || local_cli.diag_lid_cont_silent ||
        local_cli.diag_lid_cont_silent_only || local_cli.diag_lid_omit_pjump_only)
    {
        const Real z_lid = french.glass_center[2] + french.glass_half_height;
        size_t n_silent = 0;
        size_t n_omit_pj = 0;
        if (local_cli.diag_lid_dens_silence || local_cli.diag_lid_cont_silent ||
            local_cli.diag_lid_cont_silent_only)
        {
            n_silent =
                hostTagAcousticWallContSilentLid(wall.getBaseParticles(), z_lid, wall_thickness);
        }
        if (local_cli.diag_lid_dens_silence || local_cli.diag_lid_cont_silent ||
            local_cli.diag_lid_omit_pjump_only)
        {
            n_omit_pj =
                hostTagAcousticWallOmitPJumpLid(wall.getBaseParticles(), z_lid, wall_thickness);
        }
        if (local_cli.diag_lid_dens_silence)
        {
            std::cout << "[ophelie][doc221][O4] lid dens silence tags: ContSilent=" << n_silent
                      << " OmitPJump=" << n_omit_pj
                      << " dens_mode=3 AcousticWallFreeSlip=0 z_lid=" << z_lid << "\n";
        }
        else if (local_cli.diag_lid_cont_silent)
        {
            std::cout << "[ophelie][doc221][O5] lid ContSilent/OmitPJump tags: ContSilent="
                      << n_silent << " OmitPJump=" << n_omit_pj
                      << " dens_mode=" << local_cli.acoustic_dens_mode
                      << " AcousticWallFreeSlip=0 z_lid=" << z_lid << "\n";
        }
        else if (local_cli.diag_lid_cont_silent_only)
        {
            std::cout << "[ophelie][doc221][O6a] ContSilent-only tags: ContSilent=" << n_silent
                      << " OmitPJump=0 dens_mode=" << local_cli.acoustic_dens_mode
                      << " AcousticWallFreeSlip=0 z_lid=" << z_lid << "\n";
        }
        else
        {
            std::cout << "[ophelie][doc221][O6b] OmitPJump-only tags: ContSilent=0 OmitPJump="
                      << n_omit_pj << " dens_mode=" << local_cli.acoustic_dens_mode
                      << " AcousticWallFreeSlip=0 z_lid=" << z_lid << "\n";
        }
    }

    const Real acoustic_cfl_eff = local_cli.acoustic_cfl * local_cli.dt_scale;
    const Real advection_cfl_eff = local_cli.advection_cfl * local_cli.dt_scale;
    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(
        glass, local_cli.u_ref, advection_cfl_eff);
    auto &advection_viscous_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionViscousTimeStepCK>(glass, local_cli.u_ref,
                                                                                    advection_cfl_eff);
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
    // Doc232/233/234: live dT/dz for geometric-face T extrapolation (shell-band T ≠ surface T).
    // band: no extrapolate; global_dt_dz: single slope; local_dt_dz: near-face bands (Doc234 E).
    Real face_dTdz_bottom_live = Real(0);
    Real face_dTdz_top_live = Real(0);
    const bool face_t_extrapolate = (local_cli.boundary_face_t_mode != "band");
    const bool face_t_local = (local_cli.boundary_face_t_mode == "local_dt_dz");
    Real *face_dTdz_bot_ptr = face_t_extrapolate ? &face_dTdz_bottom_live : nullptr;
    Real *face_dTdz_top_ptr = face_t_extrapolate ? &face_dTdz_top_live : nullptr;
    const Real z_bot_face = gc[2] - french.glass_half_height;
    const Real z_top_face = gc[2] + french.glass_half_height;
    const Real face_local_band = std::max(Real(2) * french.dp, thermal_bc.boundary_shell_thickness);
    std::cout << "[ophelie][doc233][D] live boundary_face_t_mode=" << local_cli.boundary_face_t_mode
              << " extrapolate=" << (face_t_extrapolate ? 1 : 0) << " local_band_m=" << face_local_band << "\n";
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc, z_bot_face, z_top_face, face_dTdz_bot_ptr, face_dTdz_top_ptr);
    ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK> boundary_loss_power(
        glass, t0, thermal_bc.boundary_shell_thickness, thermal_bc, z_bot_face, z_top_face, face_dTdz_bot_ptr,
        face_dTdz_top_ptr);
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
    const bool refresh_mu_from_T =
        local_cli.doc226_mu_T_only || !local_cli.literature.constant_thermal_rheology;
    // Doc226 F: keep μ(T) (+ optional cp/k lock) + ReferenceViscosity CFL in one place.
    // Must re-run after Fig8 prescribe — init below still sees near-T0 isothermal field.
    auto sync_mu_from_temperature = [&](const char *tag, bool log_particle) -> Real {
        const Real mu_now = hostUpdateFrenchThermalRheologyFromTemperature(
            glass_particles, laws,
            /*constant_properties=*/!refresh_mu_from_T, t0, rheology_stats);
        if (local_cli.doc226_mu_T_only)
        {
            // Isolate viscosity: restore cp/k particle fields to T0 diagnostic constants.
            syncVariableToHost<Real>(glass_particles, kOphelieHeatCapacityField);
            syncVariableToHost<Real>(glass_particles, "OphelieThermalConductivity");
            Real *cp_f = glass_particles.getVariableDataByName<Real>(kOphelieHeatCapacityField);
            Real *k_f = glass_particles.getVariableDataByName<Real>("OphelieThermalConductivity");
            const size_t n_th = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n_th; ++i)
            {
                cp_f[i] = cp;
                k_f[i] = k_th;
            }
            syncVariableToDevice<Real>(glass_particles, kOphelieHeatCapacityField);
            syncVariableToDevice<Real>(glass_particles, "OphelieThermalConductivity");
            rheology_stats.cp_min = rheology_stats.cp_max = rheology_stats.cp_mean = cp;
            rheology_stats.k_min = rheology_stats.k_max = rheology_stats.k_mean = k_th;
        }
        if (local_cli.doc226_mu_cap > TinyReal)
            hostClampParticleViscosityCap(glass_particles, local_cli.doc226_mu_cap, rheology_stats);
        else
        {
            rheology_stats.mu_used_min = rheology_stats.mu_raw_min;
            rheology_stats.mu_used_max = rheology_stats.mu_raw_max;
            rheology_stats.mu_used_mean = rheology_stats.mu_raw_mean;
        }
        const Real mu_ref = use_particle_mu
                                ? (rheology_stats.mu_used_max > TinyReal ? rheology_stats.mu_used_max
                                                                       : (rheology_stats.mu_raw_max > TinyReal
                                                                              ? rheology_stats.mu_raw_max
                                                                              : mu_now))
                                : (rheology_stats.mu_used_mean > TinyReal ? rheology_stats.mu_used_mean : mu_now);
        glass_viscosity.setReferenceViscosity(mu_ref);
        invalidateViscosityDependentCkKernels(glass_inner, viscous_wall_contact);
        // AdvectionViscousTimeStepCK baked μ into speed_ref_ at construction; refresh live.
        advection_viscous_time_step.refreshViscousSpeedReference();
        if (log_particle && use_particle_mu)
        {
            std::cout << "[ophelie][doc226] viscosity_mode=particle_pairwise tag=" << tag
                      << " mu_T_only=" << (local_cli.doc226_mu_T_only ? 1 : 0)
                      << " mu_used_mean=" << rheology_stats.mu_used_mean
                      << " mu_used_max=" << rheology_stats.mu_used_max
                      << " mu_cap=" << local_cli.doc226_mu_cap
                      << " n_capped=" << rheology_stats.n_mu_capped
                      << " mu_ref_cfl=" << mu_ref << std::endl;
        }
        else if (log_particle && !local_cli.literature.constant_thermal_rheology)
        {
            std::cout << "[ophelie][stage4.1] viscosity_mode=mean_mu_T_volume_average tag=" << tag
                      << " mu_ref=" << mu_ref << std::endl;
        }
        return mu_ref;
    };
    const Real mu_init = sync_mu_from_temperature("pre_settle_T0", /*log_particle=*/true);
    (void)mu_init;
    // Doc220 C: provisional H0 at current T (usually isothermal). If prescribe-T runs later,
    // H0/prev_energy/ledger accumulators are rebuilt AFTER the final IC (see flow-start block).
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
    Doc222LocalEtaPhReport doc222_local_eta_ph_last{};
    auto apply_french_background_ph = [&]() {
        if (local_cli.doc222_ph_local_eta)
        {
            doc222_local_eta_ph_last = hostApplyFrenchBackgroundHydroPressureLocalEta(
                glass_particles, rho0, local_cli.gravity_g, z_free, gc, french.glass_radius,
                local_cli.doc222_ph_local_eta_bins, local_cli.doc222_ph_signed);
        }
        else
        {
            hostApplyFrenchBackgroundHydroPressure(glass_particles, rho0, local_cli.gravity_g, z_free,
                                                   local_cli.doc222_ph_signed);
        }
    };
    if (local_cli.doc220_pprime_hydro || local_cli.background_hydro_p)
    {
        // p_h path (p' or total-p+g): Density=ρ0 so EOS≈0; hydro lives in BackgroundPressure.
        hostFreezeDensityToRho0(glass_particles, rho0);
        apply_french_background_ph();
    }
    else
    {
        hostApplyHydrostaticDensity(glass_particles, rho0, local_cli.c0, local_cli.gravity_g, z_free);
    }
    fluid_advection_step_setup.exec();
    if (!local_cli.omit_gravity_force)
        constant_gravity.exec();
    else
    {
        // Doc221 B: one-time init — zero gravity channel only. Do not wipe shared ForcePrior.
        hostZeroVecFieldFrench(glass_particles, "GravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousGravityForceCK");
    }

    if (local_cli.doc222_tip_force_selftest)
    {
        const Real tip_r = local_cli.diag_tip_soft_force_r > TinyReal ? local_cli.diag_tip_soft_force_r : Real(0.05);
        const Real c_damp = local_cli.diag_tip_soft_force_c > TinyReal ? local_cli.diag_tip_soft_force_c : Real(100);
        const auto tip_st =
            hostDoc222TipForceLifecycleSelftest(glass_particles, z_free, Real(0.5) * french.dp, tip_r, c_damp);
        std::cout << "[ophelie][doc222][B] tip_force_selftest pass=" << (tip_st.pass ? 1 : 0) << " " << tip_st.message
                  << "\n";
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc222_tip_force_selftest\",\n"
               << "  \"pass\": " << (tip_st.pass ? 1 : 0) << ",\n"
               << "  \"repeat10_ok\": " << (tip_st.repeat10_ok ? 1 : 0) << ",\n"
               << "  \"exit_zone_ok\": " << (tip_st.exit_zone_ok ? 1 : 0) << ",\n"
               << "  \"uz_down_ok\": " << (tip_st.uz_down_ok ? 1 : 0) << ",\n"
               << "  \"clear_ok\": " << (tip_st.clear_ok ? 1 : 0) << ",\n"
               << "  \"count_no_side_effect\": " << (tip_st.count_no_side_effect ? 1 : 0) << ",\n"
               << "  \"legacy_grows\": " << (tip_st.legacy_grows ? 1 : 0) << ",\n"
               << "  \"fp_z_after_1\": " << tip_st.fp_z_after_1 << ",\n"
               << "  \"fp_z_after_10\": " << tip_st.fp_z_after_10 << ",\n"
               << "  \"fp_z_after_exit\": " << tip_st.fp_z_after_exit << ",\n"
               << "  \"fp_z_legacy_10\": " << tip_st.fp_z_legacy_10 << ",\n"
               << "  \"n_touch_1\": " << tip_st.n_touch_1 << ",\n"
               << "  \"message\": \"" << tip_st.message << "\"\n"
               << "}\n";
        }
        return tip_st.pass ? 0 : 2;
    }

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
    bool doc223_logged_keep_dens_skip = false;
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
        // Never during hydro settle — fighting packing collapses acoustic dt / hangs settle.
        dens_fs_cap_n_last = 0;
        if (local_cli.dens_fs_cap_overpack && !do_dens_norm && !in_hydro_settle)
        {
            dens_fs_cap_n_last = hostCapOverpackDensityToRho0(glass_particles, rho0);
        }
        // Doc222 E2 total-p+bgp: keep ρ=ρ0 after dens so p_mom≈p_h (∇p_h↔g). Not during settle.
        if (local_cli.background_hydro_p && !local_cli.doc220_pprime_hydro && !in_hydro_settle)
            hostFreezeDensityToRho0(glass_particles, rho0);
        if (local_cli.background_hydro_p)
            apply_french_background_ph();
        // Doc105 P2b: after p_h write, restore atm condition on incomplete-support particles.
        if (local_cli.dens_fs_atm_consistent)
            fs_atm_last = hostApplyFrenchFsAtmConsistentDensity(glass_particles, rho0, local_cli.c0, z_free,
                                                                dens_audit_surface_depth);
        // Note: hydro_rho_freeze_rho0 is applied only in stage/keep probes (not settle),
        // so layout matches the unfrozen baseline.
    };
    auto run_density_update = [&]() {
        // Doc223 E1: paper-top keep skips dens reinit (continuity only). Settle/flush still dens.
        if (local_cli.doc223_paper_flat_slip_top && !in_hydro_settle)
        {
            if (!doc223_logged_keep_dens_skip)
            {
                std::cout << "[ophelie][doc223][E1] keep dens reinit OFF "
                             "(acoustic continuity only; settle/flush dens absorbed packing)\n";
                doc223_logged_keep_dens_skip = true;
            }
            if (local_cli.background_hydro_p)
                apply_french_background_ph();
            return;
        }
        density_regularization.exec();
        apply_quiet_hydro_after_dens();
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
    if (local_cli.kernel_correction || local_cli.enable_tvf)
        fluid_linear_correction_matrix.exec();
    if (local_cli.background_hydro_p)
        apply_french_background_ph();

    // Doc105 P2d / Doc220 E2: wall dens + momentum pressure convention.
    // total-p + GravityForceCK: dens_mode 1 ⇒ a_ref=−g (g in ForcePrior).
    // p' + omit GravityForceCK: dens_mode 1 ⇒ a_ref=0 (g absorbed in BackgroundPressure).
    {
        hostSetAcousticDensModeFrench(glass_particles, static_cast<Real>(local_cli.acoustic_dens_mode));
        hostSetPerturbationMomentumFrench(glass_particles, static_cast<Real>(local_cli.perturbation_momentum));
        const bool dens1_needs_g_ref =
            local_cli.acoustic_dens_mode == 1 && !local_cli.omit_gravity_force;
        const Vecd a_ref = dens1_needs_g_ref ? Vecd(0.0, 0.0, -local_cli.gravity_g) : Vecd::Zero();
        hostSetReferenceBodyAccelerationFrench(glass_particles, a_ref);
        const char *path =
            local_cli.doc220_pprime_hydro
                ? "p_prime_plus_BackgroundPressure_no_GravityForceCK"
                : (local_cli.background_hydro_p ? "total_p_plus_BackgroundPressure_plus_GravityForceCK"
                                                : "total_p_plus_GravityForceCK");
        std::cout << "[ophelie][doc105-p2d] acoustic_dens_mode=" << local_cli.acoustic_dens_mode
                  << " reference_body_acc_z=" << a_ref[2]
                  << " perturbation_momentum=" << local_cli.perturbation_momentum
                  << " omit_gravity_force=" << (local_cli.omit_gravity_force ? 1 : 0)
                  << " french_pressure_path=" << path << " closed_box_p2b_copied=0\n";
    }
    auto apply_gravity_force_prior = [&]() {
        if (local_cli.omit_gravity_force)
        {
            // Doc221 B: g lives in BackgroundPressure; GravityForceCK was never advanced.
            // Do not touch ForcePrior / other force histories each advection step.
            // Doc221 A evidence-only: optional legacy wipe reproduces the Doc220 bug.
            if (local_cli.doc221_legacy_omit_clears_forceprior)
            {
                hostZeroVecFieldFrench(glass_particles, "GravityForceCK");
                hostZeroVecFieldFrench(glass_particles, "PreviousGravityForceCK");
                hostZeroVecFieldFrench(glass_particles, "ForcePrior");
            }
        }
        else
        {
            constant_gravity.exec();
        }
    };
    auto doc223_reinit_force_priors = [&]() {
        hostDoc223ZeroAllForcePriorChannels(glass_particles);
        apply_gravity_force_prior();
        if (local_cli.literature.enable_boussinesq)
            boussinesq_force.exec();
        exec_viscous();
    };

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    // Doc219 D / Doc226 F: min(inertial, viscous) when hydro profile or local-μ path is on.
    const bool use_viscous_adv_cfl =
        local_cli.doc219_hydro_profile || use_particle_mu || refresh_mu_from_T;
    const Real dt_adv_inertial0 = advection_time_step.exec();
    const Real dt_adv_viscous0 = advection_viscous_time_step.exec();
    const Real dt_adv0 =
        use_viscous_adv_cfl ? std::min(dt_adv_inertial0, dt_adv_viscous0) : dt_adv_inertial0;
    auto &advection_step = time_stepper.addTriggerByInterval(dt_adv0);
    auto &state_recording =
        time_stepper.addTriggerByInterval(std::max(local_cli.state_record_interval, Real(1.0e-6)));
    std::cout << "[ophelie][doc219][hydro] dt_adv0=" << dt_adv0 << " inertial=" << dt_adv_inertial0
              << " viscous=" << dt_adv_viscous0
              << " viscous_active="
              << (use_viscous_adv_cfl && dt_adv_viscous0 + TinyReal < dt_adv_inertial0 ? 1 : 0)
              << std::endl;

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
              << " doc219_viscous_dt=" << (local_cli.doc219_hydro_profile ? 1 : 0) << "\n";

    std::ofstream monitor("output/french_natural_convection_monitor.csv");
    monitor << "advection_step,physical_time_s,U_max,T_min,T_mean,T_max,dT,U_buoy,escaped,escaped_frac,"
               "P_joule,dE_dt,energy_residual,energy_residual_rel,total_heat_loss,em_updates,"
               "U_rms,U_theta_rms,U_z_rms,T_std,T_outer,T_center,T_bottom,T_top\n";
    writeFrenchSpatialStatsCsvHeader("output/french_spatial_stats.csv");

    if (local_cli.doc223_paper_flat_slip_top)
    {
        run_density_update();
        hostLogDoc223PaperTopHydroAudit("pre_settle", glass_particles, rho0, local_cli.c0, local_cli.gravity_g,
                                        gc, french.glass_radius, z_bottom, z_free, french.dp);
    }

    // Client GUI: settle was silent for minutes → looked "stuck after EM".
    // Doc245: exit on U_quiet / dt collapse / max adv steps, not only max physical time
    // (tiny acoustic dt can make a fixed settle_time take wall-clock forever).
    size_t settle_heartbeat_last = 0;
    size_t settle_quiet_streak = 0;
    size_t settle_dt_collapse_streak = 0;
    std::string settle_exit_reason = "max_time";
    Real settle_dt_floor = local_cli.hydro_settle_dt_floor;
    if (settle_dt_floor <= TinyReal && dt_ac_est > TinyReal)
        settle_dt_floor = Real(1.0e-3) * dt_ac_est;
    std::cout << "[ophelie][doc245][settle] criteria: max_t=" << local_cli.hydro_settle_time
              << " min_t=" << local_cli.hydro_settle_min_time
              << " U_quiet=" << local_cli.hydro_settle_u_quiet
              << " quiet_steps=" << local_cli.hydro_settle_quiet_steps
              << " max_adv_steps=" << local_cli.hydro_settle_max_adv_steps
              << " dt_floor=" << settle_dt_floor
              << " dt_collapse_steps=" << local_cli.hydro_settle_dt_collapse_steps << std::endl;
    while (time_stepper.getPhysicalTime() + TinyReal < local_cli.hydro_settle_time)
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        if (settle_dt_floor > TinyReal && acoustic_dt + TinyReal < settle_dt_floor)
            ++settle_dt_collapse_streak;
        else
            settle_dt_collapse_streak = 0;
        if (local_cli.hydro_settle_dt_collapse_steps > 0 &&
            settle_dt_collapse_streak >= local_cli.hydro_settle_dt_collapse_steps)
        {
            settle_exit_reason = "dt_ac_collapse";
            std::cout << "[ophelie][doc245][settle] WARNING acoustic_dt collapsed: dt_ac=" << acoustic_dt
                      << " floor=" << settle_dt_floor << " streak=" << settle_dt_collapse_streak
                      << " → end settle early (will zero-U and continue)\n";
            break;
        }
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
            run_density_update();
            maybe_dens_sum_audit(time_stepper.getPhysicalTime(), advection_steps);
            fluid_advection_step_setup.exec();
            exec_viscous();
            apply_gravity_force_prior();
            // Doc219: B (+ TVF) after ForcePrior, before next acoustic accumulation.
            if (local_cli.kernel_correction || local_cli.enable_tvf)
                fluid_linear_correction_matrix.exec();
            if (local_cli.enable_tvf)
            {
                kernel_grad_integral.exec();
                transport_correction.exec();
            }

            const Real t_settle = time_stepper.getPhysicalTime();
            const Real u_settle_now = u_max_reduce.exec();
            if (local_cli.hydro_settle_u_quiet > TinyReal &&
                t_settle + TinyReal >= local_cli.hydro_settle_min_time)
            {
                if (u_settle_now < local_cli.hydro_settle_u_quiet)
                    ++settle_quiet_streak;
                else
                    settle_quiet_streak = 0;
                if (local_cli.hydro_settle_quiet_steps > 0 &&
                    settle_quiet_streak >= local_cli.hydro_settle_quiet_steps)
                {
                    settle_exit_reason = "u_quiet";
                    std::cout << "[ophelie][doc245][settle] U_quiet reached: U_max=" << u_settle_now
                              << " < " << local_cli.hydro_settle_u_quiet << " for "
                              << settle_quiet_streak << " adv steps at t=" << t_settle << std::endl;
                    break;
                }
            }
            if (local_cli.hydro_settle_max_adv_steps > 0 &&
                advection_steps >= local_cli.hydro_settle_max_adv_steps)
            {
                settle_exit_reason = "max_adv_steps";
                std::cout << "[ophelie][doc245][settle] WARNING hit max_adv_steps="
                          << local_cli.hydro_settle_max_adv_steps << " at t=" << t_settle
                          << " U_max=" << u_settle_now << " dt_ac~" << acoustic_dt
                          << " → end settle early\n";
                break;
            }

            if (advection_steps == 1 || advection_steps - settle_heartbeat_last >= 100)
            {
                settle_heartbeat_last = advection_steps;
                const Real t_target = local_cli.hydro_settle_time;
                const Real pct = t_target > TinyReal ? Real(100) * t_settle / t_target : Real(100);
                std::cout << "[ophelie][stage4.1] hydro settle: N=" << advection_steps
                          << " t=" << t_settle << " / " << t_target << " (" << pct << "%)"
                          << " U_max=" << u_settle_now << " dt_ac=" << acoustic_dt
                          << " quiet_streak=" << settle_quiet_streak << std::endl;
            }
        }
    }
    if (local_cli.hydro_settle_time > TinyReal)
    {
        const Real u_settle = u_max_reduce.exec();
        const GlassEscapeCensus settle_escape = hostCountEscapedGlassParticles(
            glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
        std::cout << "[ophelie][stage4.1] hydro settle done: reason=" << settle_exit_reason
                  << " t=" << time_stepper.getPhysicalTime()
                  << " U_max=" << u_settle << " escaped=" << settle_escape.escaped << "/" << settle_escape.total
                  << " r_max=" << settle_escape.r_max << " z=[" << settle_escape.z_min << "," << settle_escape.z_max
                  << "]" << std::endl;
        if (local_cli.doc223_paper_flat_slip_top)
        {
            hostLogDoc223PaperTopHydroAudit("post_settle", glass_particles, rho0, local_cli.c0,
                                            local_cli.gravity_g, gc, french.glass_radius, z_bottom, z_free,
                                            french.dp);
        }
        if (!local_cli.doc235_off_event_dump_dir.empty())
        {
            dumpDoc235OffEventCsv(glass_particles,
                                  (fs::path(local_cli.doc235_off_event_dump_dir) / "event_01_settle_done_pre_cleanup.csv")
                                      .string(),
                                  "settle_done_pre_cleanup", time_stepper.getPhysicalTime(), gc, french.glass_radius,
                                  z_bottom, z_free);
        }
        hostZeroVelocity(glass_particles);
        if (!local_cli.doc235_off_event_dump_dir.empty())
        {
            dumpDoc235OffEventCsv(glass_particles,
                                  (fs::path(local_cli.doc235_off_event_dump_dir) / "event_02_after_zero_U.csv").string(),
                                  "after_zero_U", time_stepper.getPhysicalTime(), gc, french.glass_radius, z_bottom,
                                  z_free);
        }
        // Doc223 E1 / Doc230 isothermal: zero-U alone leaves DensityChangeRate/Force/Displacement
        // from the last acoustic cycle; keep then re-excites U (client NC tripped Doc224 0.05 gate).
        // Doc235 B fix1: dens flush + EOS Pressure sync + ForcePrior(gravity) — fixed stale (ρ,p,FP)
        // but did NOT reduce first-acoustic / t=0.01 U_rms (bit-match baseline).
        // Doc235 B fix2 (zero_q): skip dens summation flush; keep settle Density/Pressure, only clear
        // DCR/Force/Displacement and rebuild ForcePrior(gravity). dens flush itself injects an
        // EOS ∇p that is not hydrostatically balanced → acoustic kick.
        if (local_cli.doc223_paper_flat_slip_top || local_cli.literature.zero_q)
        {
            hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
            hostZeroVecFieldFrench(glass_particles, "Force");
            hostZeroVecFieldFrench(glass_particles, "Displacement");
            const bool doc235_skip_dens_flush = local_cli.literature.zero_q && !local_cli.doc223_paper_flat_slip_top;
            if (!doc235_skip_dens_flush)
            {
                run_density_update();
                fluid_advection_step_setup.exec();
                if (local_cli.kernel_correction || local_cli.enable_tvf)
                    fluid_linear_correction_matrix.exec();
            }
            // Rebuild mechanics (keep U=0). Skip EOS rewrite when dens was not refreshed.
            if (local_cli.doc220_pprime_hydro || local_cli.background_hydro_p)
            {
                hostFreezeDensityToRho0(glass_particles, rho0);
                apply_french_background_ph();
            }
            else if (!doc235_skip_dens_flush)
            {
                hostSyncEosPressureFromDensity(glass_particles, rho0, local_cli.c0);
            }
            hostDoc223ZeroAllForcePriorChannels(glass_particles);
            apply_gravity_force_prior();
            hostZeroVelocity(glass_particles);
            hostZeroVecFieldFrench(glass_particles, "Force");
            if (local_cli.doc223_paper_flat_slip_top)
            {
                hostLogDoc223PaperTopHydroAudit("post_settle_flush", glass_particles, rho0, local_cli.c0,
                                                local_cli.gravity_g, gc, french.glass_radius, z_bottom, z_free,
                                                french.dp);
                std::cout << "[ophelie][doc223][E1] post-settle flush: dens/adv + EOS-Pressure sync + "
                             "ForcePrior(gravity) rebuild (Doc235)\n";
            }
            else if (doc235_skip_dens_flush)
            {
                std::cout << "[ophelie][doc235][B] post-settle flush: skip dens; keep settle Density/Pressure; "
                             "clear DCR/Force/Displacement + ForcePrior(gravity); U/Force cleared\n";
            }
            else
            {
                std::cout << "[ophelie][doc235][B] post-settle flush: dens/adv + EOS-Pressure sync + "
                             "ForcePrior(gravity) rebuild; U/Force cleared\n";
            }
            if (!local_cli.doc235_off_event_dump_dir.empty())
            {
                dumpDoc235OffEventCsv(
                    glass_particles,
                    (fs::path(local_cli.doc235_off_event_dump_dir) / "event_03_after_post_settle_flush.csv").string(),
                    "after_post_settle_flush", time_stepper.getPhysicalTime(), gc, french.glass_radius, z_bottom,
                    z_free);
            }
        }
        if (local_cli.diag_lid_free_slip_post_settle)
        {
            // Doc221 P2: arm free-slip only after settle packing + zero-U (O-quiet settle).
            const Real z_lid = french.glass_center[2] + french.glass_half_height;
            size_t n_lid = 0;
            if (!local_cli.doc224_ablate_acoustic_freeslip)
            {
                n_lid = hostArmDoc224FreeSlipLid(wall.getBaseParticles(), z_lid, wall_thickness,
                                                 local_cli.doc224_lid_keep_dens_jump);
            }
            // Doc224 C3c: acoustic free-slip can arm without OpenCup.
            // Doc224 C3d: OpenCup without AcousticWallFreeSlip tags.
            if (!local_cli.doc224_ablate_opencup)
            {
                lid_free_slip_armed = true;
                invalidateViscosityDependentCkKernels(glass_inner, *glass_wall_noslip_contact);
            }
            std::cout << "[ophelie][doc221][P2] armed free-slip after settle: AcousticWallFreeSlip lid tags="
                      << n_lid << " / wall_closed_n=" << wall.getBaseParticles().TotalRealParticles()
                      << " z_lid=" << z_lid << " U_cleared=1"
                      << " acoustic_fs=" << (local_cli.doc224_ablate_acoustic_freeslip ? 0 : 1)
                      << " viscous=" << (local_cli.doc224_ablate_opencup ? "ClosedCrucible_ABLATION" : "OpenCup")
                      << " keep_dens_jump=" << (local_cli.doc224_lid_keep_dens_jump ? 1 : 0) << "\n";
            if (local_cli.doc224_ablate_opencup)
            {
                std::cout << "[ophelie][doc224][C3c] ablate_opencup=1 "
                             "(acoustic free-slip ON; wall viscous stays ClosedCrucible)\n";
            }
            if (local_cli.doc224_ablate_acoustic_freeslip)
            {
                std::cout << "[ophelie][doc224][C3d] ablate_acoustic_freeslip=1 "
                             "(AcousticWallFreeSlip OFF; viscous OpenCup — O3-style impermeable slip)\n";
            }
        }
        if (local_cli.diag_hard_lid)
        {
            // Doc221 M: settle often leaves tip above z_free under p'; hard-lid keep must
            // start from a flushed geometry or the first clamps invert meridional polarity.
            const size_t n_reset = hostResetHardLidGeometry(glass_particles, z_free);
            update_glass_cell_linked_list.exec();
            update_wall_cell_linked_list.exec();
            update_glass_relations.exec();
            update_wall_relations.exec();
            density_regularization.exec();
            if (local_cli.doc220_pprime_hydro || local_cli.background_hydro_p)
            {
                hostFreezeDensityToRho0(glass_particles, rho0);
                apply_french_background_ph();
            }
            const GlassEscapeCensus after_reset = hostCountEscapedGlassParticles(
                glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
            std::cout << "[ophelie][doc221][M] hard-lid post-settle reset n_clamp=" << n_reset
                      << " z_lid=" << z_free << " z_max_after=" << after_reset.z_max
                      << " (vel zeroed; NOT literature open FS)\n";
        }
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
        if (!local_cli.doc235_off_event_dump_dir.empty())
        {
            dumpDoc235OffEventCsv(glass_particles,
                                  (fs::path(local_cli.doc235_off_event_dump_dir) / "event_04_t0_ready.csv").string(),
                                  "t0_ready", time_stepper.getPhysicalTime(), gc, french.glass_radius, z_bottom,
                                  z_free);
        }
        in_hydro_settle = false;
    }
    else
    {
        in_hydro_settle = false;
        if (local_cli.diag_lid_free_slip_post_settle && !lid_free_slip_armed)
        {
            hostZeroVelocity(glass_particles);
            const Real z_lid = french.glass_center[2] + french.glass_half_height;
            size_t n_lid = 0;
            if (!local_cli.doc224_ablate_acoustic_freeslip)
            {
                n_lid = hostArmDoc224FreeSlipLid(wall.getBaseParticles(), z_lid, wall_thickness,
                                                 local_cli.doc224_lid_keep_dens_jump);
            }
            if (!local_cli.doc224_ablate_opencup)
            {
                lid_free_slip_armed = true;
                invalidateViscosityDependentCkKernels(glass_inner, *glass_wall_noslip_contact);
            }
            std::cout << "[ophelie][doc221][P2] armed free-slip (no settle): tags=" << n_lid
                      << " z_lid=" << z_lid
                      << " acoustic_fs=" << (local_cli.doc224_ablate_acoustic_freeslip ? 0 : 1)
                      << " viscous=" << (local_cli.doc224_ablate_opencup ? "ClosedCrucible_ABLATION" : "OpenCup")
                      << " keep_dens_jump=" << (local_cli.doc224_lid_keep_dens_jump ? 1 : 0) << "\n";
        }
    }

    // Doc223 B: ForcePrior reset probe (settle → legacy zero-FP+g.exec → full channel reinit).
    if (local_cli.doc223_force_reset_probe)
    {
        if (local_cli.omit_gravity_force)
        {
            std::cerr << "[ophelie][doc223][B] --doc223-force-reset-probe needs GravityForceCK "
                         "(do not use --doc220-pprime-hydro)\n";
            return 2;
        }
        auto dump_stage = [&](const char *tag) {
            hostDoc223LogForceStage(tag, glass_particles, local_cli.gravity_g);
            return hostDoc223ForceFieldMoments(glass_particles, "ForcePrior");
        };
        const auto st_settle = dump_stage("01_post_settle");
        const auto grav_settle = hostDoc223ForceFieldMoments(glass_particles, "GravityForceCK");
        const auto prev_settle = hostDoc223ForceFieldMoments(glass_particles, "PreviousGravityForceCK");
        hostZeroVecFieldFrench(glass_particles, "ForcePrior");
        const auto st_zero_fp = dump_stage("02_legacy_zero_ForcePrior_only");
        constant_gravity.exec();
        const auto st_legacy_g = dump_stage("03_legacy_gravity_exec");
        const auto grav_legacy = hostDoc223ForceFieldMoments(glass_particles, "GravityForceCK");
        const auto prev_legacy = hostDoc223ForceFieldMoments(glass_particles, "PreviousGravityForceCK");
        doc223_reinit_force_priors();
        const auto st_reinit = dump_stage("04_full_channel_reinit");
        const auto grav_reinit = hostDoc223ForceFieldMoments(glass_particles, "GravityForceCK");
        const auto prev_reinit = hostDoc223ForceFieldMoments(glass_particles, "PreviousGravityForceCK");
        const Real g = local_cli.gravity_g;
        const bool lost = std::abs(st_legacy_g.a_z_mean) < Real(0.1) * g &&
                          std::abs(grav_legacy.a_z_mean + g) < Real(0.1) * g;
        const bool restored = std::abs(st_reinit.a_z_mean + g) < Real(0.05) * g;
        std::cout << "[ophelie][doc223][B] verdict gravity_lost_legacy=" << (lost ? 1 : 0)
                  << " gravity_restored_reinit=" << (restored ? 1 : 0)
                  << " a_fp_z settle=" << st_settle.a_z_mean << " legacy=" << st_legacy_g.a_z_mean
                  << " reinit=" << st_reinit.a_z_mean << "\n";
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc223_force_reset_probe\",\n"
               << "  \"gravity_g\": " << g << ",\n"
               << "  \"post_settle_a_fp_z\": " << st_settle.a_z_mean << ",\n"
               << "  \"post_settle_grav_a_z\": " << grav_settle.a_z_mean << ",\n"
               << "  \"post_settle_prevG_a_z\": " << prev_settle.a_z_mean << ",\n"
               << "  \"legacy_zero_fp_a_fp_z\": " << st_zero_fp.a_z_mean << ",\n"
               << "  \"legacy_after_g_a_fp_z\": " << st_legacy_g.a_z_mean << ",\n"
               << "  \"legacy_after_g_grav_a_z\": " << grav_legacy.a_z_mean << ",\n"
               << "  \"legacy_after_g_prevG_a_z\": " << prev_legacy.a_z_mean << ",\n"
               << "  \"reinit_a_fp_z\": " << st_reinit.a_z_mean << ",\n"
               << "  \"reinit_grav_a_z\": " << grav_reinit.a_z_mean << ",\n"
               << "  \"reinit_prevG_a_z\": " << prev_reinit.a_z_mean << ",\n"
               << "  \"gravity_lost_legacy\": " << (lost ? 1 : 0) << ",\n"
               << "  \"gravity_restored_reinit\": " << (restored ? 1 : 0) << ",\n"
               << "  \"pass\": " << (lost && restored ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return (lost && restored) ? 0 : 1;
    }

    // Doc223 E2: lid tangential slip operator (wall viscous only; snapshot/restore).
    if (local_cli.doc223_e2_lid_slip_probe)
    {
        if (!lid_free_slip_armed || !exec_viscous_wall_opencup)
        {
            std::cerr << "[ophelie][doc223][E2] need post-settle free-slip armed + OpenCup viscous\n";
            return 2;
        }
        const Real u_tx = Real(0.01);
        const Real z_lid = french.glass_center[2] + french.glass_half_height;
        std::vector<Vecd> snap_vel, snap_visc, snap_prev_visc, snap_fp;
        hostDoc223SnapshotVecField(glass_particles, "Velocity", snap_vel);
        hostDoc223SnapshotVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223SnapshotVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223SnapshotVecField(glass_particles, "ForcePrior", snap_fp);

        hostDoc223SetUniformTangentialVelocity(glass_particles, u_tx);

        hostZeroVecFieldFrench(glass_particles, "ViscousForce");
        exec_viscous_wall_opencup();
        const auto slip_lid =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp,
                                 "lid_core");
        const auto slip_side =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp, "side");
        const auto slip_bot =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp, "bot");

        hostZeroVecFieldFrench(glass_particles, "ViscousForce");
        exec_viscous_wall_closed();
        const auto closed_lid =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp,
                                 "lid_core");
        const auto closed_side =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp, "side");
        const auto closed_bot =
            hostDoc223E2ViscBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid, french.dp, "bot");

        hostDoc223RestoreVecField(glass_particles, "Velocity", snap_vel);
        hostDoc223RestoreVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223RestoreVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223RestoreVecField(glass_particles, "ForcePrior", snap_fp);

        const Real lid_ratio =
            (closed_lid.a_tx_rms > TinyReal) ? (slip_lid.a_tx_rms / closed_lid.a_tx_rms) : Real(1e30);
        const Real slip_lid_vs_side =
            (slip_side.a_tx_rms > TinyReal) ? (slip_lid.a_tx_rms / slip_side.a_tx_rms) : Real(0);
        // OpenCup lid traction ~0 vs Closed no-slip lid; side still no-slip under both.
        const bool pass_lid = (slip_lid.n >= 10) && (closed_lid.n >= 10) && (lid_ratio < Real(0.05)) &&
                              (closed_lid.a_tx_rms > Real(1e-3));
        const bool pass_side_alive = (slip_side.n >= 10) && (slip_side.a_tx_rms > Real(1e-3));
        const bool pass_lid_vs_side = slip_lid_vs_side < Real(0.1);
        const bool pass = pass_lid && pass_side_alive && pass_lid_vs_side;
        std::cout << "[ophelie][doc223][E2] u_tx=" << u_tx
                  << " slip_lid n=" << slip_lid.n << " a_tx=" << slip_lid.a_tx_rms
                  << " a_tz=" << slip_lid.a_tz_rms
                  << " closed_lid n=" << closed_lid.n << " a_tx=" << closed_lid.a_tx_rms
                  << " lid_ratio=" << lid_ratio
                  << " slip_side a_tx=" << slip_side.a_tx_rms
                  << " closed_side a_tx=" << closed_side.a_tx_rms
                  << " slip_bot a_tx=" << slip_bot.a_tx_rms
                  << " closed_bot a_tx=" << closed_bot.a_tx_rms
                  << " slip_lid/side=" << slip_lid_vs_side
                  << " pass_lid=" << (pass_lid ? 1 : 0)
                  << " pass_side=" << (pass_side_alive ? 1 : 0)
                  << " pass_lid_vs_side=" << (pass_lid_vs_side ? 1 : 0)
                  << " pass=" << (pass ? 1 : 0) << "\n";
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc223_e2_lid_slip_probe\",\n"
               << "  \"u_tx\": " << u_tx << ",\n"
               << "  \"slip_lid_n\": " << slip_lid.n << ",\n"
               << "  \"slip_lid_a_tx_rms\": " << slip_lid.a_tx_rms << ",\n"
               << "  \"slip_lid_a_tz_rms\": " << slip_lid.a_tz_rms << ",\n"
               << "  \"closed_lid_n\": " << closed_lid.n << ",\n"
               << "  \"closed_lid_a_tx_rms\": " << closed_lid.a_tx_rms << ",\n"
               << "  \"lid_ratio_slip_over_closed\": " << lid_ratio << ",\n"
               << "  \"slip_side_a_tx_rms\": " << slip_side.a_tx_rms << ",\n"
               << "  \"closed_side_a_tx_rms\": " << closed_side.a_tx_rms << ",\n"
               << "  \"slip_bot_a_tx_rms\": " << slip_bot.a_tx_rms << ",\n"
               << "  \"closed_bot_a_tx_rms\": " << closed_bot.a_tx_rms << ",\n"
               << "  \"slip_lid_over_side\": " << slip_lid_vs_side << ",\n"
               << "  \"pass_lid\": " << (pass_lid ? 1 : 0) << ",\n"
               << "  \"pass_side_alive\": " << (pass_side_alive ? 1 : 0) << ",\n"
               << "  \"pass_lid_vs_side\": " << (pass_lid_vs_side ? 1 : 0) << ",\n"
               << "  \"pass\": " << (pass ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return pass ? 0 : 1;
    }

    // Doc226 B: top total operator fixed-state probe (acoustic + continuity + viscous).
    if (local_cli.doc226_top_operator_probe)
    {
        if (!lid_free_slip_armed || !exec_viscous_wall_opencup)
        {
            std::cerr << "[ophelie][doc226][B] need post-settle OpenCup viscous armed\n";
            return 2;
        }
        const Real u_tx = Real(0.01);
        const Real z_lid = french.glass_center[2] + french.glass_half_height;
        const Real h_sm = glass.getSPHAdaptation().ReferenceSmoothingLength();
        const Real rho0_b = rho0;
        const Real dt_probe = Real(1.0e-6); // Interact only; Update not required for force dump
        FrenchHydroFieldSnapshot snap_hydro;
        hostCaptureFrenchHydroFields(glass_particles, snap_hydro);
        std::vector<Vecd> snap_visc, snap_prev_visc;
        hostDoc223SnapshotVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223SnapshotVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);

        auto run_state = [&](Real u_set, const char *tag, Doc226TopOpBand &ac1, Doc226TopOpBand &ac2,
                             Doc226TopOpBand &visc_open, Doc226TopOpBand &visc_closed) {
            hostRestoreFrenchHydroFields(glass_particles, snap_hydro);
            hostDoc223SetUniformTangentialVelocity(glass_particles, u_set);

            hostZeroVecFieldFrench(glass_particles, "Force");
            hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
            if (local_cli.kernel_correction)
                acoustic_1st_wall_corr_only.exec(dt_probe);
            else
                acoustic_1st_wall_only.exec(dt_probe);
            ac1 = hostDoc226LidCoreForceBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid,
                                             french.dp, "Force", u_set);
            std::cout << "[ophelie][doc226][B] " << tag << " ac1st_wall lid n=" << ac1.n
                      << " a_t=" << ac1.a_t_rms << " a_tx=" << ac1.a_tx_rms << " a_tz=" << ac1.a_tz_rms
                      << " drho=" << ac1.drho_dt_rms << " work_tx=" << ac1.work_tx << "\n";

            hostZeroVecFieldFrench(glass_particles, "Force");
            hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
            acoustic_step_2nd_half_wall.exec(dt_probe);
            ac2 = hostDoc226LidCoreForceBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid,
                                             french.dp, "Force", u_set);
            std::cout << "[ophelie][doc226][B] " << tag << " ac2nd_wall lid n=" << ac2.n
                      << " a_t=" << ac2.a_t_rms << " a_tx=" << ac2.a_tx_rms << " a_tz=" << ac2.a_tz_rms
                      << " drho=" << ac2.drho_dt_rms << " work_tx=" << ac2.work_tx << "\n";

            hostZeroVecFieldFrench(glass_particles, "ViscousForce");
            exec_viscous_wall_opencup();
            visc_open = hostDoc226LidCoreForceBand(glass_particles, gc, french.glass_radius, z_bottom,
                                                   z_lid, french.dp, "ViscousForce", u_set);
            std::cout << "[ophelie][doc226][B] " << tag << " visc_OpenCup lid n=" << visc_open.n
                      << " a_t=" << visc_open.a_t_rms << " a_tx=" << visc_open.a_tx_rms
                      << " work_tx=" << visc_open.work_tx << "\n";

            hostZeroVecFieldFrench(glass_particles, "ViscousForce");
            exec_viscous_wall_closed();
            visc_closed = hostDoc226LidCoreForceBand(glass_particles, gc, french.glass_radius, z_bottom,
                                                     z_lid, french.dp, "ViscousForce", u_set);
            std::cout << "[ophelie][doc226][B] " << tag << " visc_Closed lid n=" << visc_closed.n
                      << " a_t=" << visc_closed.a_t_rms << " a_tx=" << visc_closed.a_tx_rms
                      << " work_tx=" << visc_closed.work_tx << "\n";
        };

        Doc226TopOpBand u0_ac1{}, u0_ac2{}, u0_vo{}, u0_vc{};
        Doc226TopOpBand ut_ac1{}, ut_ac2{}, ut_vo{}, ut_vc{};
        run_state(Real(0), "U0", u0_ac1, u0_ac2, u0_vo, u0_vc);
        run_state(u_tx, "Utx", ut_ac1, ut_ac2, ut_vo, ut_vc);

        // Aggregate tangential accel under slip path (ac1+ac2+OpenCup) vs noslip viscous ref.
        const Real a_t_slip =
            std::sqrt(ut_ac1.a_t_rms * ut_ac1.a_t_rms + ut_ac2.a_t_rms * ut_ac2.a_t_rms +
                      ut_vo.a_t_rms * ut_vo.a_t_rms);
        const Real a_t_slip_u0 =
            std::sqrt(u0_ac1.a_t_rms * u0_ac1.a_t_rms + u0_ac2.a_t_rms * u0_ac2.a_t_rms +
                      u0_vo.a_t_rms * u0_vo.a_t_rms);
        const Real a_t_delta = std::sqrt(std::max(Real(0), a_t_slip * a_t_slip - a_t_slip_u0 * a_t_slip_u0));
        const Real a_noslip_ref = std::max(ut_vc.a_t_rms, TinyReal);
        const Real a_thresh = std::max(Real(1.0e-6), Real(0.01) * a_noslip_ref);
        const Real a_ratio = a_t_delta / a_noslip_ref;
        std::string a_gate = "FAIL";
        if (a_t_delta <= a_thresh)
            a_gate = "PASS";
        else if (a_ratio <= Real(0.05))
            a_gate = "WARN";

        const Real cont_scale = rho0_b * u_tx / std::max(h_sm, TinyReal);
        const Real cont_delta =
            std::sqrt(std::max(Real(0), ut_ac2.drho_dt_rms * ut_ac2.drho_dt_rms -
                                            u0_ac2.drho_dt_rms * u0_ac2.drho_dt_rms));
        const Real cont_norm = cont_delta / std::max(cont_scale, TinyReal);
        std::string cont_gate = "FAIL";
        if (cont_norm <= Real(0.01))
            cont_gate = "PASS";
        else if (cont_norm <= Real(0.05))
            cont_gate = "WARN";

        const Real visc_ratio =
            (ut_vc.a_t_rms > TinyReal) ? (ut_vo.a_t_rms / ut_vc.a_t_rms) : Real(1e30);
        const bool visc_pass = (ut_vo.n >= 10) && (ut_vc.n >= 10) && (visc_ratio < Real(0.05)) &&
                               (ut_vc.a_t_rms > Real(1e-3));

        // B3: small normal-perturb fixed-state — wall ac2 dens from Uz only; restore after.
        hostRestoreFrenchHydroFields(glass_particles, snap_hydro);
        {
            syncVariableToHost<Vecd>(glass_particles, "Position");
            syncVariableToHost<Vecd>(glass_particles, "Velocity");
            Vecd *vel = glass_particles.getVariableDataByName<Vecd>("Velocity");
            const Vecd *pos = glass_particles.getVariableDataByName<Vecd>("Position");
            const size_t n = glass_particles.TotalRealParticles();
            const Real amp = Real(1.0e-3);
            for (size_t i = 0; i < n; ++i)
            {
                const Real dz = z_lid - pos[i][2];
                const Real w = std::exp(-dz * dz / (Real(4) * french.dp * french.dp));
                vel[i] = Vecd(Real(0), Real(0), amp * w);
            }
            syncVariableToDevice<Vecd>(glass_particles, "Velocity");
        }
        hostZeroVecFieldFrench(glass_particles, "Force");
        hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
        acoustic_step_2nd_half_wall.exec(dt_probe);
        const auto npert =
            hostDoc226LidCoreForceBand(glass_particles, gc, french.glass_radius, z_bottom, z_lid,
                                       french.dp, "Force", Real(0));
        std::cout << "[ophelie][doc226][B] Npert ac2nd lid a_tz=" << npert.a_tz_rms
                  << " a_t=" << npert.a_t_rms << " drho=" << npert.drho_dt_rms << "\n";

        hostRestoreFrenchHydroFields(glass_particles, snap_hydro);
        hostDoc223RestoreVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223RestoreVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);

        const bool pass = (a_gate == "PASS" || a_gate == "WARN") &&
                          (cont_gate == "PASS" || cont_gate == "WARN") && visc_pass;
        // Strict TOP_TOTAL_OPERATOR_CHECKED needs a_gate+cont_gate PASS; WARN still delivers data.
        const bool checked_pass = (a_gate == "PASS") && (cont_gate == "PASS") && visc_pass;

        std::cout << "[ophelie][doc226][B] SUMMARY acoustic_fs="
                  << (local_cli.doc224_ablate_acoustic_freeslip ? 0 : 1)
                  << " h=" << h_sm << " rho0=" << rho0_b
                  << " a_t_delta=" << a_t_delta << " a_noslip_ref=" << a_noslip_ref
                  << " a_ratio=" << a_ratio << " a_gate=" << a_gate
                  << " cont_delta=" << cont_delta << " cont_scale=" << cont_scale
                  << " cont_norm=" << cont_norm << " cont_gate=" << cont_gate
                  << " visc_ratio=" << visc_ratio << " visc_pass=" << (visc_pass ? 1 : 0)
                  << " checked_pass=" << (checked_pass ? 1 : 0) << " pass=" << (pass ? 1 : 0)
                  << "\n";

        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << std::scientific << std::setprecision(8);
            js << "{\n"
               << "  \"case\": \"doc226_top_operator_probe\",\n"
               << "  \"u_tx\": " << u_tx << ",\n"
               << "  \"h\": " << h_sm << ",\n"
               << "  \"rho0\": " << rho0_b << ",\n"
               << "  \"acoustic_freeslip_tags\": "
               << (local_cli.doc224_ablate_acoustic_freeslip ? 0 : 1) << ",\n"
               << "  \"U0_ac1_a_t_rms\": " << u0_ac1.a_t_rms << ",\n"
               << "  \"U0_ac2_a_t_rms\": " << u0_ac2.a_t_rms << ",\n"
               << "  \"U0_ac2_drho_dt_rms\": " << u0_ac2.drho_dt_rms << ",\n"
               << "  \"U0_visc_open_a_t_rms\": " << u0_vo.a_t_rms << ",\n"
               << "  \"Utx_ac1_a_t_rms\": " << ut_ac1.a_t_rms << ",\n"
               << "  \"Utx_ac1_work_tx\": " << ut_ac1.work_tx << ",\n"
               << "  \"Utx_ac2_a_t_rms\": " << ut_ac2.a_t_rms << ",\n"
               << "  \"Utx_ac2_drho_dt_rms\": " << ut_ac2.drho_dt_rms << ",\n"
               << "  \"Utx_ac2_work_tx\": " << ut_ac2.work_tx << ",\n"
               << "  \"Utx_visc_open_a_t_rms\": " << ut_vo.a_t_rms << ",\n"
               << "  \"Utx_visc_closed_a_t_rms\": " << ut_vc.a_t_rms << ",\n"
               << "  \"Utx_visc_open_work_tx\": " << ut_vo.work_tx << ",\n"
               << "  \"a_t_delta\": " << a_t_delta << ",\n"
               << "  \"a_noslip_ref\": " << a_noslip_ref << ",\n"
               << "  \"a_ratio\": " << a_ratio << ",\n"
               << "  \"a_gate\": \"" << a_gate << "\",\n"
               << "  \"cont_delta\": " << cont_delta << ",\n"
               << "  \"cont_scale_rho0_U_h\": " << cont_scale << ",\n"
               << "  \"cont_norm\": " << cont_norm << ",\n"
               << "  \"cont_gate\": \"" << cont_gate << "\",\n"
               << "  \"visc_ratio_open_over_closed\": " << visc_ratio << ",\n"
               << "  \"visc_pass\": " << (visc_pass ? 1 : 0) << ",\n"
               << "  \"npert_a_tz_rms\": " << npert.a_tz_rms << ",\n"
               << "  \"npert_drho_dt_rms\": " << npert.drho_dt_rms << ",\n"
               << "  \"TOP_TOTAL_OPERATOR_CHECKED\": \""
               << (checked_pass ? "PASS" : (pass ? "WARN" : "FAIL")) << "\",\n"
               << "  \"pass\": " << (pass ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return pass ? 0 : 1;
    }

    // Doc226 E1: local-μ read-path probe (uniform + half-domain swap; snapshot/restore).
    if (local_cli.doc226_e1_local_mu_probe)
    {
        if (!use_particle_mu || !exec_viscous)
        {
            std::cerr << "[ophelie][doc226][E1] need ParticleViscosity + exec_viscous\n";
            return 2;
        }
        const Real u_tx = Real(0.01);
        const Real mu_lo = Real(4);
        const Real mu_hi = Real(40);
        const Real z_mid = gc[2];
        FrenchHydroFieldSnapshot snap_hydro;
        hostCaptureFrenchHydroFields(glass_particles, snap_hydro);
        std::vector<Vecd> snap_visc, snap_prev_visc, snap_fp;
        hostDoc223SnapshotVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223SnapshotVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223SnapshotVecField(glass_particles, "ForcePrior", snap_fp);
        syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
        const Real *mu0p = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        std::vector<Real> snap_mu(mu0p, mu0p + glass_particles.TotalRealParticles());

        auto set_mu_all = [&](Real m) {
            syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
            Real *mu_f = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
            const size_t n = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n; ++i)
                mu_f[i] = m;
            syncVariableToDevice<Real>(glass_particles, kOphelieDynamicViscosityField);
        };
        auto set_mu_half_z = [&](Real below, Real above) {
            syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
            syncVariableToHost<Vecd>(glass_particles, "Position");
            Real *mu_f = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
            const Vecd *pos = glass_particles.getVariableDataByName<Vecd>("Position");
            const size_t n = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n; ++i)
                mu_f[i] = (pos[i][2] < z_mid) ? below : above;
            syncVariableToDevice<Real>(glass_particles, kOphelieDynamicViscosityField);
        };
        auto visc_rms = [&]() {
            syncVariableToHost<Vecd>(glass_particles, "ViscousForce");
            const Vecd *fv = glass_particles.getVariableDataByName<Vecd>("ViscousForce");
            const size_t n = glass_particles.TotalRealParticles();
            Real s2 = 0;
            for (size_t i = 0; i < n; ++i)
                s2 += fv[i].squaredNorm();
            return std::sqrt(s2 / Real(std::max<size_t>(n, 1)));
        };

        hostDoc223SetUniformTangentialVelocity(glass_particles, u_tx);

        set_mu_all(mu_lo);
        hostZeroVecFieldFrench(glass_particles, "ViscousForce");
        exec_viscous();
        const Real rms_uniform = visc_rms();
        const auto ev_uniform = hostCollectLocalViscosityForceEvidence(glass_particles);

        set_mu_half_z(mu_lo, mu_hi);
        hostZeroVecFieldFrench(glass_particles, "ViscousForce");
        exec_viscous();
        const Real rms_swap_a = visc_rms();
        const auto ev_a = hostCollectLocalViscosityForceEvidence(glass_particles);

        set_mu_half_z(mu_hi, mu_lo); // swap blocks; volume-mean μ unchanged if equal volumes
        hostZeroVecFieldFrench(glass_particles, "ViscousForce");
        exec_viscous();
        const Real rms_swap_b = visc_rms();
        const auto ev_b = hostCollectLocalViscosityForceEvidence(glass_particles);

        // Restore
        hostRestoreFrenchHydroFields(glass_particles, snap_hydro);
        hostDoc223RestoreVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223RestoreVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223RestoreVecField(glass_particles, "ForcePrior", snap_fp);
        syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
        Real *mu_f = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        for (size_t i = 0; i < snap_mu.size(); ++i)
            mu_f[i] = snap_mu[i];
        syncVariableToDevice<Real>(glass_particles, kOphelieDynamicViscosityField);

        const Real rel_ab =
            std::abs(rms_swap_a - rms_swap_b) / (std::max(rms_swap_a, rms_swap_b) + TinyReal);
        const bool force_changes_vs_uniform =
            std::abs(rms_swap_a - rms_uniform) > Real(1.0e-3) * (rms_uniform + Real(1.0e-12));
        const bool swap_changes = rel_ab > Real(1.0e-4) ||
                                  std::abs(ev_a.viscous_force_rms_high_mu - ev_b.viscous_force_rms_high_mu) >
                                      Real(1.0e-6);
        const bool harmonic_ok = ev_uniform.harmonic_mean_ok;
        const bool pass = harmonic_ok && force_changes_vs_uniform && (rms_uniform > TinyReal);

        std::cout << "[ophelie][doc226][E1] u_tx=" << u_tx
                  << " rms_uniform=" << rms_uniform << " rms_swap_lo_hi=" << rms_swap_a
                  << " rms_swap_hi_lo=" << rms_swap_b << " rel_ab=" << rel_ab
                  << " harmonic_ok=" << (harmonic_ok ? 1 : 0)
                  << " force_vs_uniform=" << (force_changes_vs_uniform ? 1 : 0)
                  << " swap_changes=" << (swap_changes ? 1 : 0)
                  << " pass=" << (pass ? 1 : 0) << "\n";
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << std::scientific << std::setprecision(8);
            js << "{\n"
               << "  \"case\": \"doc226_e1_local_mu_probe\",\n"
               << "  \"rms_uniform_mu4\": " << rms_uniform << ",\n"
               << "  \"rms_swap_mu4_below_mu40_above\": " << rms_swap_a << ",\n"
               << "  \"rms_swap_mu40_below_mu4_above\": " << rms_swap_b << ",\n"
               << "  \"rel_swap_ab\": " << rel_ab << ",\n"
               << "  \"harmonic_mean_ok\": " << (harmonic_ok ? 1 : 0) << ",\n"
               << "  \"force_changes_vs_uniform\": " << (force_changes_vs_uniform ? 1 : 0) << ",\n"
               << "  \"swap_spatially_changes\": " << (swap_changes ? 1 : 0) << ",\n"
               << "  \"LOCAL_MU_READ_PATH_VERIFIED\": \"" << (pass ? "PASS" : "FAIL") << "\",\n"
               << "  \"pass\": " << (pass ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return pass ? 0 : 1;
    }

    // Doc227 C1: Viscosity (Fixed μ_ref) vs ParticleViscosity (uniform field μ=4) same-state forces.
    if (local_cli.doc227_c_uniform_mu_probe)
    {
        if (!use_particle_mu || !exec_viscous)
        {
            std::cerr << "[ophelie][doc227][C] need ParticleViscosity + exec_viscous\n";
            return 2;
        }
        const Real mu_eq = Real(4);
        const Real u_tx = Real(0.01);
        FrenchHydroFieldSnapshot snap_hydro;
        hostCaptureFrenchHydroFields(glass_particles, snap_hydro);
        std::vector<Vecd> snap_visc, snap_prev_visc, snap_fp;
        hostDoc223SnapshotVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223SnapshotVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223SnapshotVecField(glass_particles, "ForcePrior", snap_fp);
        syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
        const Real *mu0p = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        std::vector<Real> snap_mu(mu0p, mu0p + glass_particles.TotalRealParticles());

        // Fixed-μ ViscousForceCK (base Viscosity APIs) with same wall post-contact order.
        using ViscousForceFixedInnerCK =
            fluid_dynamics::ViscousForceCK<Inner<WithUpdate, Viscosity, NoKernelCorrectionCK>>;
        using ViscousForceFixedWallCK =
            fluid_dynamics::ViscousForceCK<Contact<Wall, Viscosity, NoKernelCorrectionCK>>;
        auto fixed_inner =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceFixedInnerCK>>(glass_inner);
        Contact<> &wall_for_c =
            (lid_free_slip_armed && glass_wall_noslip_contact) ? *glass_wall_noslip_contact
                                                              : viscous_wall_contact;
        auto fixed_wall =
            makeUnique<InteractionDynamicsCK<MainExecutionPolicy, ViscousForceFixedWallCK>>(wall_for_c);
        fixed_inner->addPostContactInteraction(*fixed_wall);

        auto set_mu_all = [&](Real m) {
            syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
            Real *mu_f = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
            const size_t n = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n; ++i)
                mu_f[i] = m;
            syncVariableToDevice<Real>(glass_particles, kOphelieDynamicViscosityField);
            glass_viscosity.setReferenceViscosity(m);
            invalidateViscosityDependentCkKernels(glass_inner, wall_for_c);
        };
        auto field_rms = [&](const char *name) {
            syncVariableToHost<Vecd>(glass_particles, name);
            const Vecd *v = glass_particles.getVariableDataByName<Vecd>(name);
            Real s2 = 0;
            const size_t n = glass_particles.TotalRealParticles();
            for (size_t i = 0; i < n; ++i)
                s2 += v[i].squaredNorm();
            return std::sqrt(s2 / Real(std::max<size_t>(n, 1)));
        };

        hostDoc223SetUniformTangentialVelocity(glass_particles, u_tx);
        set_mu_all(mu_eq);

        auto run_and_snap = [&](bool particle_path, std::vector<Vecd> &out_visc, std::vector<Vecd> &out_fp) {
            // Same baseline each path: restore ForcePrior so WithUpdate does not
            // accumulate particle-path viscous onto the fixed-path second exec.
            hostDoc223RestoreVecField(glass_particles, "ForcePrior", snap_fp);
            hostZeroVecFieldFrench(glass_particles, "ViscousForce");
            hostZeroVecFieldFrench(glass_particles, "PreviousViscousForce");
            if (particle_path)
                exec_viscous();
            else
                fixed_inner->exec();
            hostDoc223SnapshotVecField(glass_particles, "ViscousForce", out_visc);
            hostDoc223SnapshotVecField(glass_particles, "ForcePrior", out_fp);
        };

        std::vector<Vecd> visc_p, fp_p, visc_f, fp_f;
        run_and_snap(true, visc_p, fp_p);
        const Real rms_visc_p = field_rms("ViscousForce");
        const Real rms_fp_p = field_rms("ForcePrior");
        run_and_snap(false, visc_f, fp_f);
        const Real rms_visc_f = field_rms("ViscousForce");
        const Real rms_fp_f = field_rms("ForcePrior");

        // Restore particle-path result into host arrays already in visc_p; compare Fixed vs Particle snaps.
        Real max_visc = 0, rms_visc = 0, max_fp = 0, rms_fp = 0;
        {
            const size_t n = visc_p.size();
            Real s_v = 0, s_f = 0;
            for (size_t i = 0; i < n; ++i)
            {
                const Real dv = (visc_p[i] - visc_f[i]).norm();
                const Real df = (fp_p[i] - fp_f[i]).norm();
                max_visc = std::max(max_visc, dv);
                max_fp = std::max(max_fp, df);
                s_v += dv * dv;
                s_f += df * df;
            }
            rms_visc = std::sqrt(s_v / Real(std::max<size_t>(n, 1)));
            rms_fp = std::sqrt(s_f / Real(std::max<size_t>(n, 1)));
        }
        const Real scale_visc = std::max({rms_visc_p, rms_visc_f, Real(1.0e-12)});
        const Real scale_fp = std::max({rms_fp_p, rms_fp_f, Real(1.0e-12)});
        const Real rel_visc = max_visc / scale_visc;
        const Real rel_fp = max_fp / scale_fp;
        // Plan: single-precision non-zero force relative ≤1e-5; near-zero use abs+scale.
        const Real rel_tol = Real(1.0e-5);
        const Real abs_tol = Real(1.0e-8);
        const bool visc_ok = (max_visc <= abs_tol) || (rel_visc <= rel_tol);
        const bool fp_ok = (max_fp <= abs_tol) || (rel_fp <= rel_tol);
        const bool pass = visc_ok && fp_ok && (rms_visc_p > TinyReal) && (rms_visc_f > TinyReal);

        std::cout << "[ophelie][doc227][C1] uniform_mu=" << mu_eq << " u_tx=" << u_tx
                  << " rms_visc_particle=" << rms_visc_p << " rms_visc_fixed=" << rms_visc_f
                  << " max_abs_visc=" << max_visc << " rel_visc=" << rel_visc
                  << " max_abs_fp=" << max_fp << " rel_fp=" << rel_fp
                  << " wall=" << (lid_free_slip_armed ? "OpenCup" : "Closed")
                  << " pass=" << (pass ? 1 : 0) << "\n";

        hostRestoreFrenchHydroFields(glass_particles, snap_hydro);
        hostDoc223RestoreVecField(glass_particles, "ViscousForce", snap_visc);
        hostDoc223RestoreVecField(glass_particles, "PreviousViscousForce", snap_prev_visc);
        hostDoc223RestoreVecField(glass_particles, "ForcePrior", snap_fp);
        syncVariableToHost<Real>(glass_particles, kOphelieDynamicViscosityField);
        Real *mu_f = glass_particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
        for (size_t i = 0; i < snap_mu.size(); ++i)
            mu_f[i] = snap_mu[i];
        syncVariableToDevice<Real>(glass_particles, kOphelieDynamicViscosityField);

        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << std::scientific << std::setprecision(8);
            js << "{\n"
               << "  \"case\": \"doc227_c_uniform_mu_probe\",\n"
               << "  \"mu\": " << mu_eq << ",\n"
               << "  \"rms_visc_particle\": " << rms_visc_p << ",\n"
               << "  \"rms_visc_fixed\": " << rms_visc_f << ",\n"
               << "  \"max_abs_visc_diff\": " << max_visc << ",\n"
               << "  \"rel_visc_diff\": " << rel_visc << ",\n"
               << "  \"max_abs_forceprior_diff\": " << max_fp << ",\n"
               << "  \"rel_forceprior_diff\": " << rel_fp << ",\n"
               << "  \"CONSTANT_MU_PATH_EQUIVALENT\": \"" << (pass ? "PASS" : "FAIL") << "\",\n"
               << "  \"pass\": " << (pass ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return pass ? 0 : 1;
    }

    // Doc223 E3: near-lid normal velocity IC (mean-zero); integrate in main loop.
    Doc223E3State doc223_e3_state{};
    if (local_cli.doc223_e3_normal_perturb)
    {
        if (!lid_free_slip_armed)
        {
            std::cerr << "[ophelie][doc223][E3] need post-settle free-slip armed\n";
            return 2;
        }
        const Real z_lid = french.glass_center[2] + french.glass_half_height;
        const Real amp = Real(1.0e-3);
        doc223_e3_state = hostDoc223ApplyE3NormalPerturb(glass_particles, gc, french.glass_radius, z_lid,
                                                         french.dp, local_cli.c0, amp);
        hostDoc223LogE3LidAudit("ic", glass_particles, z_lid, french.dp, doc223_e3_state.ke0,
                                doc223_e3_state.z_com0);
        std::cout << "[ophelie][doc223][E3] applied near-lid Uz bump amp=" << amp
                  << " Ma=" << doc223_e3_state.ma << " U_rms0=" << doc223_e3_state.u_rms0
                  << " ke0=" << doc223_e3_state.ke0 << " z_com0=" << doc223_e3_state.z_com0
                  << " (mass-mean removed; no stepwise Uz kill)\n";
    }

    // Doc222 D2: post-settle free-surface deformation (R0; no tip). Re-init dens/p per path.
    Doc222FsDeformReport doc222_fs_deform_rep{};
    Doc222FsSurfaceStats doc222_fs_surf0{};
    Doc222FsSurfaceStats doc222_fs_surf_last{};
    bool doc222_fs_deform_on = local_cli.doc222_fs_deform_eps > TinyReal;
    bool doc222_fs_csv_header = true;
    const std::string doc222_fs_csv_path = "output/doc222_fs_deform_surface.csv";
    if (doc222_fs_deform_on)
    {
        if (local_cli.diag_tip_soft_force || local_cli.diag_tip_fs_dashpot || local_cli.diag_soft_lid ||
            local_cli.diag_soft_lid_tip || local_cli.diag_hard_lid || local_cli.diag_zero_net_uz)
        {
            std::cout << "[ophelie][doc222][D2] WARNING: tip/lid/zero-net still on — D2 requires R0; "
                         "continuing but results are not a clean restoring test\n";
        }
        doc222_fs_deform_rep = hostApplyDoc222AxisymmetricFsDeform(
            glass_particles, gc, french.glass_radius, z_bottom, z_free, local_cli.doc222_fs_deform_eps);
        doc222_fs_deform_rep.eps_over_dp =
            (french.dp > TinyReal) ? (local_cli.doc222_fs_deform_eps / french.dp) : Real(0);
        hostZeroVelocity(glass_particles);
        hostZeroVecFieldFrench(glass_particles, "Force");
        update_glass_cell_linked_list.exec();
        update_wall_cell_linked_list.exec();
        update_glass_relations.exec();
        update_wall_relations.exec();
        density_regularization.exec();
        if (local_cli.doc220_pprime_hydro)
        {
            // Doc222 E2a: default keep FreeSurface Density after deform so p' can be nonzero.
            // Legacy --doc222-deform-rho-freeze restores ρ→ρ0 (p_mom≡0, E1 ineffective).
            if (!local_cli.doc222_deform_skip_rho_freeze)
                hostFreezeDensityToRho0(glass_particles, rho0);
            apply_french_background_ph();
        }
        else if (local_cli.background_hydro_p)
        {
            // Doc222 E2 total-p+bgp: ρ=ρ0 ⇒ p_mom≈p_h; with local_eta, −∇p_h/ρ+g ⇒ −g∇η.
            hostFreezeDensityToRho0(glass_particles, rho0);
            apply_french_background_ph();
        }
        else
        {
            hostApplyHydrostaticDensity(glass_particles, rho0, local_cli.c0, local_cli.gravity_g, z_free);
        }
        fluid_advection_step_setup.exec();
        doc223_reinit_force_priors();
        fs::create_directories("output");
        doc222_fs_surf0 =
            hostMeasureDoc222FsSurfaceStats(glass_particles, gc, french.glass_radius, z_free, 12);
        doc222_fs_surf_last = doc222_fs_surf0;
        hostDoc222AppendFsSurfaceCsv(doc222_fs_csv_path, Real(0), 0, doc222_fs_surf0, true);
        doc222_fs_csv_header = false;
        std::cout << "[ophelie][doc222][D2] fs_deform eps=" << local_cli.doc222_fs_deform_eps
                  << " eps/dp=" << doc222_fs_deform_rep.eps_over_dp
                  << " eps/H=" << doc222_fs_deform_rep.eps_over_H << " eta_c=" << doc222_fs_deform_rep.eta_center
                  << " eta_rim=" << doc222_fs_deform_rep.eta_rim
                  << " z_max " << doc222_fs_deform_rep.z_max_before << "->" << doc222_fs_deform_rep.z_max_after
                  << " mode_amp0=" << doc222_fs_surf0.mode_amp
                  << " n_bins_valid=" << doc222_fs_surf0.n_bins_valid
                  << " path=" << (local_cli.doc220_pprime_hydro ? "pprime" : "total_p")
                  << " skip_rho_freeze=" << (local_cli.doc222_deform_skip_rho_freeze ? 1 : 0)
                  << " (column stretch; dens/p reinit; NO tip)\n";
        if (local_cli.doc220_pprime_hydro)
        {
            std::cout << "[ophelie][doc222][E2a] after_deform "
                      << (local_cli.doc222_deform_skip_rho_freeze
                              ? "Density=FreeSurface (allow p'); "
                              : "Density frozen to rho0 (legacy p_mom~0); ")
                      << "p_h signed=" << (local_cli.doc222_ph_signed ? 1 : 0)
                      << " local_eta=" << (local_cli.doc222_ph_local_eta ? 1 : 0);
            if (local_cli.doc222_ph_local_eta)
            {
                std::cout << " eta_spread=" << doc222_local_eta_ph_last.eta_spread
                          << " eta_mean=" << doc222_local_eta_ph_last.eta_mean
                          << " bins_valid=" << doc222_local_eta_ph_last.n_bins_valid;
            }
            std::cout << "\n";
        }
    }

    // Doc222 C3: same deformed geometry — clamp vs signed p_h acoustic-1st force compare; then exit.
    if (local_cli.doc222_c3_ph_compare)
    {
        if (!local_cli.doc220_pprime_hydro)
        {
            std::cerr << "[ophelie][doc222][C3] requires --doc220-pprime-hydro (or --doc222-c3-ph-compare)\n";
            return 2;
        }
        if (!doc222_fs_deform_on)
        {
            std::cerr << "[ophelie][doc222][C3] fs deform did not run (eps<=0)\n";
            return 2;
        }
        const Real zone_band_c3 = Real(2) * french.dp;
        FrenchHydroFieldSnapshot c3_base;
        hostCaptureFrenchHydroFields(glass_particles, c3_base);
        auto run_c3_probe = [&](bool signed_ph) -> Doc222C3PhProbeStats {
            hostRestoreFrenchHydroFields(glass_particles, c3_base);
            hostZeroVelocity(glass_particles);
            hostZeroVecFieldFrench(glass_particles, "Force");
            hostZeroRealFieldFrench(glass_particles, "DensityChangeRate");
            density_regularization.exec();
            // Match E2a deform path: skip ρ→ρ0 unless legacy --doc222-deform-rho-freeze.
            if (!local_cli.doc222_deform_skip_rho_freeze)
                hostFreezeDensityToRho0(glass_particles, rho0);
            hostApplyFrenchBackgroundHydroPressure(glass_particles, rho0, local_cli.gravity_g, z_free, signed_ph);
            fluid_advection_step_setup.exec();
            if (local_cli.kernel_correction)
                fluid_linear_correction_matrix.exec();
            hostDoc223ZeroAllForcePriorChannels(glass_particles);
            apply_gravity_force_prior();
            hostZeroVecFieldFrench(glass_particles, "Force");
            const Real dt_ac = acoustic_time_step.exec();
            if (local_cli.kernel_correction)
                acoustic_1st_corr.exec(dt_ac);
            else
                acoustic_step_1st_half.exec(dt_ac);
            return hostMeasureDoc222C3PhProbe(glass_particles, rho0, local_cli.c0, z_free, gc, french.glass_radius,
                                              z_bottom, zone_band_c3, signed_ph);
        };
        const Doc222C3PhProbeStats st_signed = run_c3_probe(true);
        const Doc222C3PhProbeStats st_clamp = run_c3_probe(false);
        hostRestoreFrenchHydroFields(glass_particles, c3_base);
        fs::create_directories("output");
        const std::string c3_json =
            local_cli.summary_json.empty() ? "output/doc222_c3_ph_compare.json" : local_cli.summary_json;
        hostDoc222WriteC3PhCompareJson(c3_json, st_signed, st_clamp, local_cli.doc222_fs_deform_eps, z_free);
        const Real d_pbg = st_signed.p_bg_mean_above - st_clamp.p_bg_mean_above;
        const Real d_af = st_signed.a_force_z_mean_above - st_clamp.a_force_z_mean_above;
        const Real d_an = st_signed.a_net_z_mean_above - st_clamp.a_net_z_mean_above;
        std::cout << "[ophelie][doc222][C3] ph_compare n_above=" << st_signed.n_above
                  << " z_max=" << st_signed.z_max << " z_free=" << z_free << "\n"
                  << "[ophelie][doc222][C3] signed: p_bg_above=" << st_signed.p_bg_mean_above
                  << " p_mom_above=" << st_signed.p_mom_mean_above
                  << " a_force_z_above=" << st_signed.a_force_z_mean_above
                  << " a_net_z_above=" << st_signed.a_net_z_mean_above << "\n"
                  << "[ophelie][doc222][C3] clamp:  p_bg_above=" << st_clamp.p_bg_mean_above
                  << " p_mom_above=" << st_clamp.p_mom_mean_above
                  << " a_force_z_above=" << st_clamp.a_force_z_mean_above
                  << " a_net_z_above=" << st_clamp.a_net_z_mean_above << "\n"
                  << "[ophelie][doc222][C3] delta(signed-clamp): p_bg_above=" << d_pbg
                  << " a_force_z_above=" << d_af << " a_net_z_above=" << d_an << "\n"
                  << "[ophelie][doc222][C3] summary_json=" << c3_json << "\n";
        // p_bg must differ above z_free if signed works; |d_pbg| ~ O(ρ0 g * η)
        const bool pbg_differs = std::abs(d_pbg) > Real(1.0); // Pa-scale threshold
        const bool a_differs = std::abs(d_af) > Real(1.0e-4) || std::abs(d_an) > Real(1.0e-4);
        std::cout << "[ophelie][doc222][C3] verdict p_bg_differs=" << (pbg_differs ? 1 : 0)
                  << " a_force_differs=" << (a_differs ? 1 : 0)
                  << (pbg_differs && !a_differs
                          ? " (p_h reaches field but acoustic Force insensitive — look beyond E1)"
                          : (pbg_differs && a_differs
                                 ? " (p_h changes Force; D2 FAIL driven by time-integrated/other terms)"
                                 : (!pbg_differs ? " (p_bg ABOVE identical — apply/path bug)"
                                                 : "")))
                  << "\n";
        return (st_signed.n_above > 0) ? 0 : 2;
    }

    // Doc220 S0: zero-physics sort identity audit (exit after compare).
    if (local_cli.doc220_s0_sort_audit)
    {
        fs::create_directories("output");
        if (local_cli.prescribe_t_radial_dT <= TinyReal)
            local_cli.prescribe_t_radial_dT = Real(200);
        hostPrescribeRadialHotCoreTemperature(glass_particles, french.glass_center, french.glass_radius,
                                              french.glass_half_height, t0, local_cli.prescribe_t_radial_dT,
                                              local_cli.prescribe_radial_z_taper, kOphelieTemperatureField,
                                              kOphelieThermalDeltaTField);
        const std::string before_csv = "output/doc220_s0_before_sort.csv";
        const std::string after_csv = "output/doc220_s0_after_sort.csv";
        hostDoc220DumpSortAuditCsv(glass_particles, before_csv, t0);
        std::cout << "[ophelie][doc220][s0] before_sort written " << before_csv
                  << " temperature_in_evolving="
                  << (hostDiscreteVarInEvolving(glass_particles, kOphelieTemperatureField) ? 1 : 0)
                  << " delta_t_in_evolving="
                  << (hostDiscreteVarInEvolving(glass_particles, kOphelieThermalDeltaTField) ? 1 : 0) << std::endl;
        fluid_particle_sorting.exec();
        update_glass_cell_linked_list.exec();
        update_glass_relations.exec();
        hostDoc220DumpSortAuditCsv(glass_particles, after_csv, t0);
        const Doc220SortAuditCompare cmp =
            hostDoc220CompareSortAuditCsv(before_csv, after_csv, glass_particles, t0);
        std::cout << "[ophelie][doc220][s0] after_sort n=" << cmp.n << " matched=" << cmp.n_matched
                  << " pos_mismatch=" << cmp.n_pos_mismatch << " T_mismatch=" << cmp.n_T_mismatch
                  << " dT_mismatch=" << cmp.n_dT_mismatch << " TdT_id_fail=" << cmp.n_TdT_identity_fail
                  << " mass_mismatch=" << cmp.n_mass_mismatch << " max_abs_dT=" << cmp.max_abs_dT
                  << " max_abs_dPos=" << cmp.max_abs_dPos << " T_mean " << cmp.T_mean_before << "->"
                  << cmp.T_mean_after << " T_range_before=[" << cmp.T_min_before << "," << cmp.T_max_before
                  << "] after=[" << cmp.T_min_after << "," << cmp.T_max_after << "] evolving_T="
                  << (cmp.temperature_in_evolving ? 1 : 0) << " evolving_dT=" << (cmp.delta_t_in_evolving ? 1 : 0)
                  << " SORT_IDENTITY_PASS=" << (cmp.passed ? 1 : 0) << std::endl;
        {
            std::ofstream js("output/doc220_s0_compare.json");
            js << "{\n"
               << "  \"case\": \"doc220_s0_sort_audit\",\n"
               << "  \"n\": " << cmp.n << ",\n"
               << "  \"n_matched\": " << cmp.n_matched << ",\n"
               << "  \"n_pos_mismatch\": " << cmp.n_pos_mismatch << ",\n"
               << "  \"n_T_mismatch\": " << cmp.n_T_mismatch << ",\n"
               << "  \"n_dT_mismatch\": " << cmp.n_dT_mismatch << ",\n"
               << "  \"n_TdT_identity_fail\": " << cmp.n_TdT_identity_fail << ",\n"
               << "  \"n_mass_mismatch\": " << cmp.n_mass_mismatch << ",\n"
               << "  \"max_abs_dT\": " << cmp.max_abs_dT << ",\n"
               << "  \"max_abs_dPos\": " << cmp.max_abs_dPos << ",\n"
               << "  \"temperature_in_evolving\": " << (cmp.temperature_in_evolving ? 1 : 0) << ",\n"
               << "  \"delta_t_in_evolving\": " << (cmp.delta_t_in_evolving ? 1 : 0) << ",\n"
               << "  \"SORT_IDENTITY_PASS\": " << (cmp.passed ? 1 : 0) << "\n"
               << "}\n";
        }
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc220_s0_sort_audit\",\n"
               << "  \"SORT_IDENTITY_PASS\": " << (cmp.passed ? 1 : 0) << ",\n"
               << "  \"temperature_in_evolving\": " << (cmp.temperature_in_evolving ? 1 : 0) << ",\n"
               << "  \"delta_t_in_evolving\": " << (cmp.delta_t_in_evolving ? 1 : 0) << ",\n"
               << "  \"n_T_mismatch\": " << cmp.n_T_mismatch << ",\n"
               << "  \"max_abs_dT\": " << cmp.max_abs_dT << ",\n"
               << "  \"passed\": " << (cmp.passed ? 1 : 0) << "\n"
               << "}\n";
        }
        return cmp.passed ? 0 : 2;
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
        hostDoc223ZeroAllForcePriorChannels(glass_particles);
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
        if (local_cli.doc220_h0_hydro_probe)
        {
            hostDoc220DumpHydroParticleCsv(glass_particles, "output/doc220_h0_01_after_dens.csv", rho0,
                                           local_cli.c0, local_cli.gravity_g, z_free, glass_c, french.glass_radius,
                                           z_bottom, zone_band, local_cli.kernel_correction);
        }

        hostZeroVecFieldFrench(glass_particles, "BoussinesqBuoyancyForceCK");
        hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
        apply_gravity_force_prior();
        // Intentionally no boussinesq here — probe isolates gravity/pressure path.
        // Doc220 E2 p' path: ForcePrior stays 0 (g in BackgroundPressure).
        FrenchHydroStageSnapshot s_prior = snap("02_after_gravity");
        if (local_cli.doc220_h0_hydro_probe)
        {
            hostDoc220DumpHydroParticleCsv(glass_particles, "output/doc220_h0_02_after_gravity.csv", rho0,
                                           local_cli.c0, local_cli.gravity_g, z_free, glass_c, french.glass_radius,
                                           z_bottom, zone_band, local_cli.kernel_correction);
        }

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
        if (local_cli.doc220_h0_hydro_probe)
        {
            hostDoc220DumpHydroParticleCsv(glass_particles, "output/doc220_h0_03a_ff_only.csv", rho0, local_cli.c0,
                                           local_cli.gravity_g, z_free, glass_c, french.glass_radius, z_bottom,
                                           zone_band, local_cli.kernel_correction);
        }

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
        if (local_cli.doc220_h0_hydro_probe)
        {
            hostDoc220DumpHydroParticleCsv(glass_particles, "output/doc220_h0_03b_wall_only.csv", rho0, local_cli.c0,
                                           local_cli.gravity_g, z_free, glass_c, french.glass_radius, z_bottom,
                                           zone_band, local_cli.kernel_correction);
        }

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
        if (local_cli.doc220_h0_hydro_probe)
        {
            hostDoc220DumpHydroParticleCsv(glass_particles, "output/doc220_h0_03_acoustic_full.csv", rho0,
                                           local_cli.c0, local_cli.gravity_g, z_free, glass_c, french.glass_radius,
                                           z_bottom, zone_band, local_cli.kernel_correction);
        }
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
               << "  \"case\": \""
               << (local_cli.doc220_h0_hydro_probe ? "doc220_h0_hydro_stage" : "french_doc105_p2_hydro_stage")
               << "\",\n"
               << "  \"doc220_h0\": " << (local_cli.doc220_h0_hydro_probe ? 1 : 0) << ",\n"
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
               << ((local_cli.acoustic_dens_mode == 1 && !local_cli.omit_gravity_force)
                       ? -local_cli.gravity_g
                       : Real(0))
               << ",\n"
               << "  \"french_pressure_path\": \""
               << (local_cli.doc220_pprime_hydro
                       ? "p_prime_plus_BackgroundPressure_no_GravityForceCK"
                       : (local_cli.background_hydro_p
                              ? "total_p_plus_BackgroundPressure_plus_GravityForceCK"
                              : "total_p_plus_GravityForceCK"))
               << "\",\n"
               << "  \"perturbation_momentum\": " << local_cli.perturbation_momentum << ",\n"
               << "  \"doc220_pprime_hydro\": " << (local_cli.doc220_pprime_hydro ? 1 : 0) << ",\n"
               << "  \"omit_gravity_force\": " << (local_cli.omit_gravity_force ? 1 : 0) << ",\n"
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
                if (local_cli.kernel_correction || local_cli.enable_tvf)
                    fluid_linear_correction_matrix.exec();
                // Doc219 F: quiet-keep must exercise TVF when enabled (F0/F1 A/B was a false
                // negative before — only settle/main flow ran transport correction).
                if (local_cli.enable_tvf)
                {
                    kernel_grad_integral.exec();
                    transport_correction.exec();
                }
                exec_viscous();
                apply_gravity_force_prior();
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
        hostDoc223ZeroAllForcePriorChannels(glass_particles);
        apply_gravity_force_prior();
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
               << "  \"enable_tvf\": " << (local_cli.enable_tvf ? 1 : 0) << ",\n"
               << "  \"doc219_hydro_profile\": " << (local_cli.doc219_hydro_profile ? 1 : 0) << ",\n"
               << "  \"acoustic_dens_mode\": " << local_cli.acoustic_dens_mode << ",\n"
               << "  \"reference_body_acc_z\": "
               << ((local_cli.acoustic_dens_mode == 1 && !local_cli.omit_gravity_force)
                       ? -local_cli.gravity_g
                       : Real(0))
               << ",\n"
               << "  \"doc220_pprime_hydro\": " << (local_cli.doc220_pprime_hydro ? 1 : 0) << ",\n"
               << "  \"perturbation_momentum\": " << local_cli.perturbation_momentum << ",\n"
               << "  \"omit_gravity_force\": " << (local_cli.omit_gravity_force ? 1 : 0) << ",\n"
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
    Doc224RzTemperatureGrid doc224_T_rz_grid;
    if (!local_cli.init_T_rz_csv.empty())
    {
        std::string terr;
        if (!hostLoadDoc224RzTemperatureGrid(local_cli.init_T_rz_csv, doc224_T_rz_grid, terr))
        {
            std::cerr << "[ophelie][doc224][D] " << terr << " path=" << local_cli.init_T_rz_csv << std::endl;
            return 2;
        }
        std::cout << "[ophelie][doc225][B] loaded init-T-rz-csv=" << local_cli.init_T_rz_csv
                  << " ASSUMED_absolute_Tmax_not_literature legacy_demean="
                  << (local_cli.init_T_rz_demean_legacy ? 1 : 0) << std::endl;
        if (local_cli.init_T_rz_demean_legacy)
        {
            std::cout << "[ophelie][doc225][B] WARNING ABSOLUTE_T_SHIFTED_LEGACY — not for acceptance mainline\n";
        }
    }
    auto apply_prescribed_temperature = [&](bool quiet = false) {
        if (!local_cli.init_T_rz_csv.empty() && doc224_T_rz_grid.ok)
        {
            hostApplyDoc224RzTemperatureGrid(glass_particles, french.glass_center, t0, doc224_T_rz_grid,
                                             local_cli.init_T_rz_demean_legacy, quiet);
        }
        else if (local_cli.doc234_fixed_t_radial_dT > TinyReal)
        {
            hostPrescribeDoc234FixedRadialMeanZeroTemperature(
                glass_particles, french.glass_center, french.glass_radius, t0,
                local_cli.doc234_fixed_t_radial_dT, kOphelieTemperatureField, kOphelieThermalDeltaTField);
            if (!quiet)
                std::cout << "[ophelie][doc234][D] fixed T=T0+dT*((r/R)^2-1/2) dT="
                          << local_cli.doc234_fixed_t_radial_dT << " T0=" << t0 << "\n";
        }
        else if (local_cli.doc223_f2_annular_hot)
        {
            hostPrescribeDoc223AnnularHotRingTemperature(glass_particles, french.glass_center,
                                                         french.glass_half_height, t0, Real(200),
                                                         kOphelieTemperatureField, kOphelieThermalDeltaTField);
        }
        else if (local_cli.prescribe_t_radial_dT > TinyReal)
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
    const bool has_prescribed_t = !local_cli.init_T_rz_csv.empty() || local_cli.doc223_f2_annular_hot ||
                                  local_cli.prescribe_t_radial_dT > TinyReal ||
                                  local_cli.prescribe_t_hot_bottom_dT > TinyReal ||
                                  local_cli.doc234_fixed_t_radial_dT > TinyReal;
    if (has_prescribed_t)
    {
        apply_prescribed_temperature(/*quiet=*/false);
        const char *init_t_id = "prescribed_other";
        if (!local_cli.init_T_rz_csv.empty())
        {
            init_t_id = local_cli.freeze_prescribed_temperature
                            ? "prescribed_doc224_fig8_rz_csv_frozen_ASSUMED"
                            : "prescribed_doc224_fig8_rz_csv_then_thaw_ASSUMED";
        }
        else if (local_cli.doc223_f2_annular_hot)
        {
            init_t_id = local_cli.freeze_prescribed_temperature
                            ? "prescribed_doc223_f2_annular_hot_frozen"
                            : "prescribed_doc223_f2_annular_hot_then_thaw";
        }
        else if (local_cli.prescribe_t_radial_dT > TinyReal)
        {
            init_t_id = local_cli.freeze_prescribed_temperature ? "prescribed_radial_hot_core_frozen"
                                                                : "prescribed_radial_hot_core_then_thaw";
        }
        else if (local_cli.prescribe_t_hot_bottom_dT > TinyReal)
        {
            init_t_id = local_cli.freeze_prescribed_temperature ? "prescribed_hot_bottom_frozen"
                                                               : "prescribed_hot_bottom_then_thaw";
        }
        const Real dT_log = !local_cli.init_T_rz_csv.empty()
                                ? static_cast<Real>(doc224_T_rz_grid.T_max - doc224_T_rz_grid.T_min)
                                : (local_cli.doc223_f2_annular_hot
                                       ? Real(200)
                                       : (local_cli.prescribe_t_radial_dT > TinyReal
                                              ? local_cli.prescribe_t_radial_dT
                                              : local_cli.prescribe_t_hot_bottom_dT));
        std::cout << "[ophelie][init_T] " << init_t_id << " dT=" << dT_log
                  << " freeze=" << (local_cli.freeze_prescribed_temperature ? 1 : 0)
                  << " zero_q=" << (local_cli.literature.zero_q ? 1 : 0)
                  << " thermal_diffusion=" << (local_cli.literature.enable_thermal_diffusion ? 1 : 0)
                  << (!local_cli.init_T_rz_csv.empty()
                          ? " (Doc224 D Fig8 ASSUMED Tmax; absolute T NOT_EVALUATED)"
                          : (local_cli.doc223_f2_annular_hot
                                 ? " (Doc223 F2 qualitative ring; NOT Fig8 digitization)"
                                 : " (not literature steady IC; thaw keeps Q+diffusion+Robin)"))
                  << std::endl;
        // Doc225 B3: Fig8 CSV → T_force_ref = mass-weighted mean once (absolute T unchanged).
        // F2 / radial / hot-bottom keep literature t0 as T_force_ref unless overridden.
        Real t_force_ref = t0;
        if (!std::isnan(local_cli.t_force_ref_override))
        {
            t_force_ref = local_cli.t_force_ref_override;
        }
        else if (!local_cli.init_T_rz_csv.empty() && local_cli.t_force_ref_from_mass_mean &&
                 !local_cli.init_T_rz_demean_legacy)
        {
            t_force_ref = hostMassWeightedMeanTemperature(glass_particles, kOphelieTemperatureField);
        }
        if (local_cli.literature.enable_boussinesq)
        {
            boussinesq_force.setForceReferenceTemperature(t_force_ref);
            boussinesq_force.exec();
            const Real t_force_used = boussinesq_force.getForceReferenceTemperature();
            const OphelieBoussinesqForceAudit audit = hostAuditBoussinesqForcePrior(
                glass_particles, beta, t_force_used, -local_cli.gravity_g, kOphelieTemperatureField);
            const OphelieBuoyancyFieldGate gate =
                hostAuditBuoyancyFieldGate(glass_particles, beta, t_force_used, local_cli.gravity_g,
                                           kOphelieTemperatureField);
            const char *tag = !local_cli.init_T_rz_csv.empty()
                                  ? "prescribed_doc225_fig8_rz"
                                  : (local_cli.doc223_f2_annular_hot
                                         ? "prescribed_doc223_f2_annular"
                                         : (local_cli.prescribe_t_radial_dT > TinyReal ? "prescribed_radial"
                                                                                      : "prescribed_hot_bottom"));
            logOphelieBoussinesqForceAudit(audit, tag);
            // Doc225 C: particle L2 + absolute net; do not FAIL near-zero ΣF via rel_err.
            const bool force_ok = gate.algebraic_ok && gate.dynamic_reference_ok;
            std::cout << "[ophelie][doc225][refs] T_thermal_ref=" << t0 << " T_property_ref=" << t0
                      << " T_force_ref=" << t_force_used << " T_mean_mass=" << gate.tmean_K
                      << " legacy_demean=" << (local_cli.init_T_rz_demean_legacy ? 1 : 0) << std::endl;
            std::cout << "[ophelie][doc225][C] E_particle=" << gate.E_L2 << " E_net_abs=" << gate.net_abs_N
                      << " max_abs_N=" << gate.max_abs_N << " sum_|Fref|=" << gate.sum_abs_fref_N
                      << " near_zero=" << (gate.near_zero_force ? 1 : 0)
                      << " algebraic_ok=" << (gate.algebraic_ok ? 1 : 0)
                      << " tref_ok=" << (gate.dynamic_reference_ok ? 1 : 0) << std::endl;
            std::cout << "[ophelie][p1-03] prescribe_mode=" << tag << " dT=" << dT_log
                      << " freeze=" << (local_cli.freeze_prescribed_temperature ? 1 : 0)
                      << " force_sign_up=" << (audit.f_z_sum > 0 ? 1 : 0) << " force_ok=" << (force_ok ? 1 : 0)
                      << " Fz_sum=" << audit.f_z_sum << " legacy_net_rel_err=" << audit.rel_err
                      << " (informational; gate uses E_particle)" << std::endl;
        }
        else
        {
            std::cout << "[ophelie][doc225][refs] T_thermal_ref=" << t0 << " T_property_ref=" << t0
                      << " T_force_ref=" << t_force_ref << " (boussinesq off)\n";
        }
        if (local_cli.doc233_buoyancy_force_chain_exit)
        {
            // Doc233/234: frozen-state force-chain (+ Doc234 precision CSV / Float analysis).
            if (!local_cli.literature.enable_boussinesq)
            {
                std::cerr << "[ophelie][doc233][C] FATAL: buoyancy force-chain requires Boussinesq ON\n";
                return 2;
            }
            syncVariableToHost<Real>(glass_particles, kOphelieTemperatureField);
            syncVariableToHost<Vecd>(glass_particles, "Position");
            Real *T_fc = glass_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
            Real *dT_fc = glass_particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
            const Vecd *pos_fc = glass_particles.getVariableDataByName<Vecd>("Position");
            const size_t n_fc = glass_particles.TotalRealParticles();
            Real tmin_fc = T_fc[0], tmax_fc = T_fc[0];
            for (size_t i = 0; i < n_fc; ++i)
            {
                tmin_fc = std::min(tmin_fc, T_fc[i]);
                tmax_fc = std::max(tmax_fc, T_fc[i]);
            }
            if (tmax_fc - tmin_fc < Real(1))
            {
                const Real R_fc = french.glass_radius;
                for (size_t i = 0; i < n_fc; ++i)
                {
                    const Real dx = pos_fc[i][0] - gc[0];
                    const Real dy = pos_fc[i][1] - gc[1];
                    const Real r = std::sqrt(dx * dx + dy * dy);
                    const Real Ti = t0 + Real(30) * (r / std::max(R_fc, TinyReal));
                    T_fc[i] = Ti;
                    dT_fc[i] = Ti - t0;
                }
                syncVariableToDevice<Real>(glass_particles, kOphelieTemperatureField);
                syncVariableToDevice<Real>(glass_particles, kOphelieThermalDeltaTField);
                std::cout << "[ophelie][doc233][C] manufactured radial ΔT=30K for force-chain (was isothermal)\n";
            }
            const Real t_ref_fc = hostMassWeightedMeanTemperature(glass_particles, kOphelieTemperatureField);
            boussinesq_force.setForceReferenceTemperature(t_ref_fc);
            // Pre-buoyancy ForcePrior (viscous + gravity only).
            hostZeroVecFieldFrench(glass_particles, "BoussinesqBuoyancyForceCK");
            hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
            hostDoc223ZeroAllForcePriorChannels(glass_particles);
            fluid_advection_step_setup.exec();
            exec_viscous();
            apply_gravity_force_prior();
            syncVariableToHost<Vecd>(glass_particles, "ForcePrior");
            syncVariableToHost<Real>(glass_particles, "Mass");
            syncVariableToHost<Real>(glass_particles, kOphelieTemperatureField);
            const Vecd *fp_pre_ptr = glass_particles.getVariableDataByName<Vecd>("ForcePrior");
            const Real *mass_fc = glass_particles.getVariableDataByName<Real>("Mass");
            const Real *T_now = glass_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
            StdVec<Vecd> fp_pre(n_fc);
            for (size_t i = 0; i < n_fc; ++i)
                fp_pre[i] = fp_pre_ptr[i];
            boussinesq_force.exec();
            const OphelieBuoyancyFieldGate gate_on =
                hostAuditBuoyancyFieldGate(glass_particles, beta, t_ref_fc, local_cli.gravity_g,
                                           kOphelieTemperatureField);
            syncVariableToHost<Vecd>(glass_particles, "ForcePrior");
            syncVariableToHost<Vecd>(glass_particles, "BoussinesqBuoyancyForceCK");
            const Vecd *fp_on = glass_particles.getVariableDataByName<Vecd>("ForcePrior");
            const Vecd *fb = glass_particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
            StdVec<Vecd> fp_on_copy(n_fc), fb_copy(n_fc);
            for (size_t i = 0; i < n_fc; ++i)
            {
                fp_on_copy[i] = fp_on[i];
                fb_copy[i] = fb[i];
            }
            // OFF: clear then rebuild viscous+gravity only.
            hostZeroVecFieldFrench(glass_particles, "BoussinesqBuoyancyForceCK");
            hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
            hostDoc223ZeroAllForcePriorChannels(glass_particles);
            fluid_advection_step_setup.exec();
            exec_viscous();
            apply_gravity_force_prior();
            syncVariableToHost<Vecd>(glass_particles, "ForcePrior");
            const Vecd *fp_off = glass_particles.getVariableDataByName<Vecd>("ForcePrior");
            StdVec<Vecd> fp_off_copy(n_fc);
            for (size_t i = 0; i < n_fc; ++i)
                fp_off_copy[i] = fp_off[i];

            const double g_pos = static_cast<double>(local_cli.gravity_g);
            const double beta_d = static_cast<double>(beta);
            const double tref_d = static_cast<double>(t_ref_fc);
            const double eps_f = static_cast<double>(std::numeric_limits<float>::epsilon());
            // Doc235 D: separate Float-native squaredNorm accumulate vs component→double promote.
            double sum_e2_float_sq = 0.0, sum_b2_float_sq = 0.0, sum_df2 = 0.0;
            double sum_e2_dprom = 0.0, sum_b2_dprom = 0.0;
            double sum_e2_vs_dref = 0.0, sum_dref2 = 0.0;
            double sum_abs_fb = 0.0, sum_abs_fp = 0.0, sum_abs_df_err = 0.0;
            double sum_ulp_scale = 0.0;
            double max_abs_ez = 0.0;
            size_t n_nonzero_fb = 0;
            for (size_t i = 0; i < n_fc; ++i)
            {
                const Vecd df = fp_on_copy[i] - fp_off_copy[i];
                const Vecd eb = df - fb_copy[i];
                sum_df2 += static_cast<double>(df.squaredNorm());
                sum_b2_float_sq += static_cast<double>(fb_copy[i].squaredNorm());
                sum_e2_float_sq += static_cast<double>(eb.squaredNorm());
                const double dfx = static_cast<double>(fp_on_copy[i][0]) - static_cast<double>(fp_off_copy[i][0]);
                const double dfy = static_cast<double>(fp_on_copy[i][1]) - static_cast<double>(fp_off_copy[i][1]);
                const double dfz = static_cast<double>(fp_on_copy[i][2]) - static_cast<double>(fp_off_copy[i][2]);
                const double fbx = static_cast<double>(fb_copy[i][0]);
                const double fby = static_cast<double>(fb_copy[i][1]);
                const double fbz = static_cast<double>(fb_copy[i][2]);
                const double ex = dfx - fbx;
                const double ey = dfy - fby;
                const double ez = dfz - fbz;
                sum_e2_dprom += ex * ex + ey * ey + ez * ez;
                sum_b2_dprom += fbx * fbx + fby * fby + fbz * fbz;
                max_abs_ez = std::max(max_abs_ez, std::abs(ez));
                const double m = static_cast<double>(mass_fc[i]);
                const double dT = static_cast<double>(T_now[i]) - tref_d;
                // F_b_ref = -m*beta*(T-Tref)*g_acc with g_acc=(0,0,-g) ⇒ Fz = +m*beta*g*dT
                const double fb_ref_z = m * beta_d * g_pos * dT;
                const double ez_d = fbz - fb_ref_z;
                sum_e2_vs_dref += ez_d * ez_d + fbx * fbx + fby * fby;
                sum_dref2 += fb_ref_z * fb_ref_z;
                const double abs_fb = std::sqrt(fbx * fbx + fby * fby + fbz * fbz);
                const double abs_fp = std::sqrt(static_cast<double>(fp_on_copy[i][0]) * static_cast<double>(fp_on_copy[i][0]) +
                                               static_cast<double>(fp_on_copy[i][1]) * static_cast<double>(fp_on_copy[i][1]) +
                                               static_cast<double>(fp_on_copy[i][2]) * static_cast<double>(fp_on_copy[i][2]));
                sum_abs_fb += abs_fb;
                sum_abs_fp += abs_fp;
                sum_abs_df_err += std::sqrt(ex * ex + ey * ey + ez * ez);
                const double F_g_scale = m * g_pos;
                sum_ulp_scale += eps_f * std::max(abs_fp, F_g_scale);
                if (abs_fb > 1.0e-12)
                    ++n_nonzero_fb;
            }
            // Legacy cast path (can collapse tiny sum_e2); keep for diagnosis only.
            const Real E_chain_float_cast =
                std::sqrt(static_cast<Real>(sum_e2_float_sq) / (static_cast<Real>(sum_b2_float_sq) + TinyReal));
            const double E_chain_float_native =
                std::sqrt(sum_e2_float_sq / (sum_b2_float_sq + 1.0e-300));
            const double E_chain_double_promoted =
                std::sqrt(sum_e2_dprom / (sum_b2_dprom + 1.0e-300));
            // Doc235 primary chain metric: component→double promote (not 50×floor gate).
            const Real E_chain = static_cast<Real>(E_chain_double_promoted);
            const Real E_chain_abs = static_cast<Real>(std::sqrt(sum_e2_dprom));
            const Real E_fb_vs_double =
                std::sqrt(static_cast<Real>(sum_e2_vs_dref) / (static_cast<Real>(sum_dref2) + TinyReal));
            const double mean_abs_fb = sum_abs_fb / std::max<size_t>(1, n_fc);
            const double mean_abs_fp = sum_abs_fp / std::max<size_t>(1, n_fc);
            const double mean_ulp = sum_ulp_scale / std::max<size_t>(1, n_fc);
            const double cancel_ratio = mean_abs_fb > 0.0 ? mean_abs_fp / mean_abs_fb : 0.0;
            const double expected_rel_floor =
                mean_abs_fb > 0.0 ? (eps_f * mean_abs_fp) / mean_abs_fb : 0.0;

            // CSV (Doc234 B)
            std::string csv_path = local_cli.doc234_force_precision_csv;
            if (csv_path.empty() && local_cli.doc234_force_precision_exit)
                csv_path = "output/doc234_force_precision_particles.csv";
            if (!csv_path.empty())
            {
                fs::create_directories(fs::path(csv_path).parent_path());
                std::ofstream csv(csv_path);
                csv << std::setprecision(10);
                csv << "particle_i,x,y,z,T_K,T_force_ref,mass,g_pos,Fb_ref_z_double,Fb_z,Fp_pre_z,Fp_on_z,Fp_off_z,"
                       "df_z,err_df_minus_Fb_z,Fp_on_norm,Fb_norm,ulp_scale\n";
                for (size_t i = 0; i < n_fc; ++i)
                {
                    const double m = static_cast<double>(mass_fc[i]);
                    const double dT = static_cast<double>(T_now[i]) - tref_d;
                    const double fb_ref_z = m * beta_d * g_pos * dT;
                    const Vecd df = fp_on_copy[i] - fp_off_copy[i];
                    const double ulp = eps_f * std::max(std::sqrt(static_cast<double>(fp_on_copy[i].squaredNorm())),
                                                         m * g_pos);
                    csv << i << "," << pos_fc[i][0] << "," << pos_fc[i][1] << "," << pos_fc[i][2] << ","
                        << T_now[i] << "," << t_ref_fc << "," << mass_fc[i] << "," << g_pos << "," << fb_ref_z
                        << "," << fb_copy[i][2] << "," << fp_pre[i][2] << "," << fp_on_copy[i][2] << ","
                        << fp_off_copy[i][2] << "," << df[2] << "," << (df[2] - fb_copy[i][2]) << ","
                        << std::sqrt(static_cast<double>(fp_on_copy[i].squaredNorm())) << ","
                        << std::sqrt(static_cast<double>(fb_copy[i].squaredNorm())) << "," << ulp << "\n";
                }
                std::cout << "[ophelie][doc234][B] force_precision_csv=" << csv_path << " n=" << n_fc << "\n";
            }

            // Repeat buoyancy — no double-add.
            for (size_t i = 0; i < n_fc; ++i)
            {
                glass_particles.getVariableDataByName<Vecd>("ForcePrior")[i] = fp_off_copy[i];
                glass_particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK")[i] = Vecd::Zero();
            }
            syncVariableToDevice<Vecd>(glass_particles, "ForcePrior");
            syncVariableToDevice<Vecd>(glass_particles, "BoussinesqBuoyancyForceCK");
            hostZeroVecFieldFrench(glass_particles, "PreviousBoussinesqBuoyancyForceCK");
            boussinesq_force.exec();
            syncVariableToHost<Vecd>(glass_particles, "ForcePrior");
            const Vecd *fp1 = glass_particles.getVariableDataByName<Vecd>("ForcePrior");
            StdVec<Vecd> fp_after1(n_fc);
            for (size_t i = 0; i < n_fc; ++i)
                fp_after1[i] = fp1[i];
            boussinesq_force.exec();
            syncVariableToHost<Vecd>(glass_particles, "ForcePrior");
            const Vecd *fp2 = glass_particles.getVariableDataByName<Vecd>("ForcePrior");
            double sum_d2 = 0.0, sum_f2 = 0.0;
            for (size_t i = 0; i < n_fc; ++i)
            {
                const Vecd d = fp2[i] - fp_after1[i];
                sum_d2 += d.squaredNorm();
                sum_f2 += fp_after1[i].squaredNorm();
            }
            const Real repeat_rel =
                std::sqrt(static_cast<Real>(sum_d2) / (static_cast<Real>(sum_f2) + TinyReal));

            // Doc234 B5: one acoustic 1st/2nd half; compare Δvz to ForcePrior_z/m*dt (pressure also acts).
            Real momentum_read_rel = Real(-1);
            Real momentum_read_abs = Real(-1);
            size_t momentum_n_interior = 0;
            {
                for (size_t i = 0; i < n_fc; ++i)
                    glass_particles.getVariableDataByName<Vecd>("ForcePrior")[i] = fp_on_copy[i];
                syncVariableToDevice<Vecd>(glass_particles, "ForcePrior");
                syncVariableToHost<Vecd>(glass_particles, "Velocity");
                StdVec<Vecd> vel0(n_fc);
                const Vecd *vel_p = glass_particles.getVariableDataByName<Vecd>("Velocity");
                for (size_t i = 0; i < n_fc; ++i)
                    vel0[i] = vel_p[i];
                const Real dt_probe = Real(1.0e-5);
                run_acoustic(dt_probe);
                syncVariableToHost<Vecd>(glass_particles, "Velocity");
                vel_p = glass_particles.getVariableDataByName<Vecd>("Velocity");
                double se2 = 0.0, sp2 = 0.0;
                for (size_t i = 0; i < n_fc; ++i)
                {
                    const Real dx = pos_fc[i][0] - gc[0];
                    const Real dy = pos_fc[i][1] - gc[1];
                    const Real r = std::sqrt(dx * dx + dy * dy);
                    if (r > Real(0.8) * french.glass_radius)
                        continue; // skip near-wall where pressure/BC dominate
                    ++momentum_n_interior;
                    const double dvz = static_cast<double>(vel_p[i][2] - vel0[i][2]);
                    const double pred =
                        static_cast<double>(fp_on_copy[i][2] / (mass_fc[i] + TinyReal) * dt_probe);
                    // Residual after removing ForcePrior prediction — pressure etc. remain.
                    // We only check that |dv| is not ≪ |pred| when |Fb| is non-trivial (read happened).
                    const double fb_n = std::sqrt(static_cast<double>(fb_copy[i].squaredNorm()));
                    if (fb_n < 1.0e-10)
                        continue;
                    se2 += (dvz - pred) * (dvz - pred);
                    sp2 += pred * pred;
                }
                momentum_read_abs = std::sqrt(static_cast<Real>(se2));
                momentum_read_rel =
                    std::sqrt(static_cast<Real>(se2) / (static_cast<Real>(sp2) + TinyReal));
            }

            const bool algebraic_ok = gate_on.algebraic_ok && gate_on.dynamic_reference_ok;
            // Absolute residual gate for near-zero Fb particles is covered by E_chain_abs + max_abs_ez.
            const bool chain_ok =
                E_chain_double_promoted <= 1.0e-3 || (sum_b2_dprom < 1.0e-16 && sum_e2_dprom < 1.0e-16);
            const bool no_double = repeat_rel <= Real(1.0e-6);
            // Doc235: expected_rel_floor is DIAGNOSTIC only — never an automatic PASS multiplier.
            const bool float_cancel_diag_only = true;
            const bool channel_verified = algebraic_ok && no_double && chain_ok;
            // Momentum read remains PARTIAL until ForcePrior contribution is recorded inside update.
            const bool momentum_read_verified = false;
            const bool pass_strict = channel_verified;
            std::cout << "[ophelie][doc233][C] BUOYANCY_FORCE_CHAIN_VERIFIED=" << (pass_strict ? 1 : 0)
                      << " T_force_ref=" << t_ref_fc << " E_particle=" << gate_on.E_L2
                      << " E_chain_ON_minus_OFF_vs_Fb=" << E_chain << " repeat_rel=" << repeat_rel
                      << " algebraic_ok=" << (algebraic_ok ? 1 : 0) << " chain_ok=" << (chain_ok ? 1 : 0)
                      << " no_double=" << (no_double ? 1 : 0) << "\n";
            std::cout << "[ophelie][doc235][D] E_chain_double_promoted=" << E_chain_double_promoted
                      << " E_chain_float_native=" << E_chain_float_native
                      << " E_chain_float_cast_legacy=" << E_chain_float_cast
                      << " E_chain_abs=" << E_chain_abs << " max_abs_ez=" << max_abs_ez
                      << " E_Fb_vs_host_double=" << E_fb_vs_double << " mean_|Fp|/mean_|Fb|=" << cancel_ratio
                      << " expected_rel_floor_DIAG_ONLY=" << expected_rel_floor
                      << " mean_ulp_scale=" << mean_ulp << " n_nonzero_Fb=" << n_nonzero_fb
                      << " CHANNEL_VERIFIED=" << (channel_verified ? 1 : 0)
                      << " MOMENTUM_READ_VERIFIED=" << (momentum_read_verified ? 1 : 0)
                      << " FLOAT_FLOOR_AUTO_PASS_DISABLED=" << (float_cancel_diag_only ? 1 : 0)
                      << " momentum_probe_rel=" << momentum_read_rel
                      << " momentum_probe_abs=" << momentum_read_abs
                      << " momentum_n_interior=" << momentum_n_interior
                      << " note=probe presence-not-equality; 50x floor not used for PASS\n";
            if (!local_cli.summary_json.empty())
            {
                std::ofstream js(local_cli.summary_json);
                js << std::setprecision(10);
                js << "{\n"
                   << "  \"case\": \"doc235_d_force_gate\",\n"
                   << "  \"FORCE_FORMULA_VERIFIED\": " << (algebraic_ok ? 1 : 0) << ",\n"
                   << "  \"CHANNEL_VERIFIED\": " << (channel_verified ? 1 : 0) << ",\n"
                   << "  \"MOMENTUM_READ_VERIFIED\": " << (momentum_read_verified ? 1 : 0) << ",\n"
                   << "  \"MOMENTUM_READ_STATUS\": \"PARTIAL\",\n"
                   << "  \"FLOAT_FLOOR_AUTO_PASS_DISABLED\": 1,\n"
                   << "  \"BUOYANCY_FORCE_CHAIN_VERIFIED\": " << (pass_strict ? 1 : 0) << ",\n"
                   << "  \"T_force_ref\": " << t_ref_fc << ",\n"
                   << "  \"beta\": " << beta << ",\n"
                   << "  \"E_particle\": " << gate_on.E_L2 << ",\n"
                   << "  \"E_chain\": " << E_chain << ",\n"
                   << "  \"E_chain_double_promoted\": " << E_chain_double_promoted << ",\n"
                   << "  \"E_chain_float_native\": " << E_chain_float_native << ",\n"
                   << "  \"E_chain_float_cast_legacy\": " << E_chain_float_cast << ",\n"
                   << "  \"E_chain_abs\": " << E_chain_abs << ",\n"
                   << "  \"max_abs_ez\": " << max_abs_ez << ",\n"
                   << "  \"E_Fb_vs_host_double\": " << E_fb_vs_double << ",\n"
                   << "  \"cancel_ratio_Fp_over_Fb\": " << cancel_ratio << ",\n"
                   << "  \"expected_rel_floor_DIAG_ONLY\": " << expected_rel_floor << ",\n"
                   << "  \"repeat_rel\": " << repeat_rel << ",\n"
                   << "  \"momentum_probe_rel\": " << momentum_read_rel << ",\n"
                   << "  \"csv\": \"" << csv_path << "\",\n"
                   << "  \"note\": \"Primary E_chain is component-promoted double; 50x floor disabled\"\n"
                   << "}\n";
                std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
            }
            return pass_strict ? 0 : 1;
        }
    }
    // Doc220 C: H0 must use the final IC temperature (after prescribe). Clear evolution ledgers.
    {
        const Real h_pre_ic = prev_energy;
        prev_energy = local_cli.literature.constant_thermal_rheology
                          ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                          : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
        prev_energy_time = 0.0;
        const Real e_prescribe = prev_energy - h_pre_ic;
        coupling_state.ledger.e_q_applied = 0.0;
        coupling_state.ledger.e_boundary_out = 0.0;
        coupling_state.ledger.setH0(prev_energy);
        if (local_cli.freeze_prescribed_temperature)
        {
            coupling_state.ledger.noteLiveLoss(Real(0), Real(0), Real(0));
        }
        else
        {
            Real side0 = 0, bottom0 = 0, free_conv0 = 0, free_rad0 = 0, total0 = 0;
            hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc,
                                                          thermal_bc.boundary_shell_thickness, side0, bottom0,
                                                          free_conv0, free_rad0, total0);
            coupling_state.ledger.noteLiveLoss(side0, bottom0, free_conv0 + free_rad0);
        }
        std::cout << "[ophelie][doc220][energy] THERMAL_INITIALIZATION H0=" << prev_energy
                  << " E_prescribe=" << e_prescribe << " (thaw: hot-core in H0; freeze: constraint energy)"
                  << " residual_t0=" << coupling_state.ledger.residualH(prev_energy) << std::endl;
    }
    if (write_vtp)
    {
        glass_state_recorder.writeToFile(0);
    }

    if (local_cli.diag_hard_lid)
    {
        // Doc221 M: re-flush after T prescribe / dens path so keep starts with z<=z_free, U=0.
        const size_t n_reset = hostResetHardLidGeometry(glass_particles, z_free);
        update_glass_cell_linked_list.exec();
        update_wall_cell_linked_list.exec();
        update_glass_relations.exec();
        update_wall_relations.exec();
        density_regularization.exec();
        if (local_cli.doc220_pprime_hydro || local_cli.background_hydro_p)
        {
            hostFreezeDensityToRho0(glass_particles, rho0);
            apply_french_background_ph();
        }
        std::cout << "[ophelie][doc221][M] hard-lid pre-keep reset n_clamp=" << n_reset
                  << " z_lid=" << z_free << "\n";
    }
    // Doc226 F: Fig8/prescribe changes T after settle — refresh μ before first keep step,
    // and shrink the already-created advection trigger (otherwise first interval stays T0/μ=4).
    if (has_prescribed_t && refresh_mu_from_T)
    {
        const Real mu_ref_post = sync_mu_from_temperature("post_prescribe_T", /*log_particle=*/true);
        const Real dt_in_post = advection_time_step.exec();
        const Real dt_vis_post = advection_viscous_time_step.exec();
        const Real dt_adv_post =
            use_viscous_adv_cfl ? std::min(dt_in_post, dt_vis_post) : dt_in_post;
        if (dt_adv_post + TinyReal < advection_step.getInterval())
            advection_step.setInterval(dt_adv_post);
        std::cout << "[ophelie][doc226] post_prescribe_adv_cfl dt_adv=" << advection_step.getInterval()
                  << " inertial=" << dt_in_post << " viscous=" << dt_vis_post
                  << " mu_ref=" << mu_ref_post << std::endl;
    }
    std::cout << "[ophelie][stage4.1] flow start: end_time=" << local_cli.end_time
              << " n_glass=" << glass_particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " beta=" << beta << " mu=" << mu
              << " k=" << k_th << " c0=" << local_cli.c0 << " T0=" << t0
              << " boussinesq=" << (local_cli.literature.enable_boussinesq ? 1 : 0)
              << " meridional_avg_start=" << local_cli.meridional_average_start
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0) << std::endl;

    // Doc232 B: production fixed-space probes (client NC).
    StdVec<OphelieFixedSpaceProbeDef> doc232_probes;
    const Real probe_kernel_h = local_cli.probe_kernel_h_factor * french.dp;
    std::unique_ptr<std::ofstream> probes_ts_stream;
    int probes_ts_samples = 0;
    bool doc232_probe_path_ok = true;
    if (local_cli.client_nc_geometry || local_cli.probes_timeseries_user_set ||
        local_cli.doc232_probe_path_verify)
    {
        doc232_probes = makeDoc232ClientNcProbes(french);
        if (!local_cli.probes_timeseries_user_set && local_cli.probes_timeseries_csv.empty() &&
            local_cli.client_nc_geometry)
            local_cli.probes_timeseries_csv = "output/probes_timeseries.csv";
        if (local_cli.doc232_probe_path_verify)
        {
            const Real z_bot_pv = french.glass_center[2] - french.glass_half_height;
            const Real melt_h_pv = Real(2) * french.glass_half_height;
            bool u_ok = false, l_ok = false;
            doc232_probe_path_ok =
                runDoc232ProbePathVerify(glass_particles, doc232_probes, probe_kernel_h, z_bot_pv, melt_h_pv, t0,
                                         "output/doc232_probe_path_verify.csv", u_ok, l_ok);
            std::cout << "[ophelie][doc232] PRODUCTION_PROBE_PATH_VERIFY="
                      << (doc232_probe_path_ok ? 1 : 0) << "\n";
            if (local_cli.doc232_probe_path_verify_exit)
            {
                if (!local_cli.summary_json.empty())
                {
                    std::ofstream js(local_cli.summary_json);
                    js << "{\n  \"case\": \"doc232_probe_path_verify\",\n"
                       << "  \"PRODUCTION_PROBE_PATH_VERIFY\": " << (doc232_probe_path_ok ? 1 : 0) << ",\n"
                       << "  \"uniform_ok\": " << (u_ok ? 1 : 0) << ",\n"
                       << "  \"linear_ok\": " << (l_ok ? 1 : 0) << ",\n"
                       << "  \"probe_kernel_h_m\": " << probe_kernel_h << ",\n"
                       << "  \"n_probes\": " << doc232_probes.size() << ",\n"
                       << "  \"note\": \"manufacture gate on production sampler; sensor location unconfirmed\"\n}\n";
                    std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
                }
                return doc232_probe_path_ok ? 0 : 1;
            }
        }
        if (!local_cli.probes_timeseries_csv.empty())
        {
            fs::create_directories(fs::path(local_cli.probes_timeseries_csv).parent_path());
            probes_ts_stream = std::make_unique<std::ofstream>(local_cli.probes_timeseries_csv);
            writeProbesTimeseriesHeader(*probes_ts_stream);
            appendProbesTimeseries(*probes_ts_stream, time_stepper.getPhysicalTime(), glass_particles,
                                   doc232_probes, probe_kernel_h);
            ++probes_ts_samples;
            std::cout << "[ophelie][doc232] probes_timeseries=" << local_cli.probes_timeseries_csv
                      << " n_probes=" << doc232_probes.size() << " kernel_h=" << probe_kernel_h << "\n";
        }
    }

    if (local_cli.doc232_boundary_face_audit)
    {
        if (!thermal_bc.enable_french_natural_bc)
        {
            std::cerr << "[ophelie][doc232][D] FATAL: french natural BC disabled; cannot audit faces\n";
            return 2;
        }
        const bool d_ok =
            runDoc232BoundaryFaceAudit(glass, glass_particles, thermal_bc, french, rho0, cp, t0,
                                       local_cli.doc232_boundary_audit_dir, local_cli.boundary_face_t_mode,
                                       local_cli.doc233_boundary_scope_audit);
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc233_d_boundary_face_audit\",\n"
               << "  \"boundary_face_t_mode\": \"" << local_cli.boundary_face_t_mode << "\",\n"
               << "  \"doc233_scope_fields\": " << (local_cli.doc233_boundary_scope_audit ? 1 : 0) << ",\n"
               << "  \"PRODUCTION_BOUNDARY_FACE_AUDIT\": " << (d_ok ? 1 : 0) << ",\n"
               << "  \"passed\": " << (d_ok ? 1 : 0) << ",\n"
               << "  \"out_dir\": \"" << local_cli.doc232_boundary_audit_dir << "\",\n"
               << "  \"power_tol\": 0.02,\n"
               << "  \"energy_tol\": 0.01,\n"
               << "  \"note\": \"production host power + ApplyOphelieThermalFrenchNaturalBcCK; "
                  "analytic ref is geometric cylinder (also offline Python)\"\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        if (local_cli.doc232_boundary_face_audit_exit)
            return d_ok ? 0 : 1;
    }

    FrenchMeridionalFlowAccumulator meridional;
    meridional.reset(cyl_frame);
    bool meridional_ts_header = true;
    const std::string meridional_ts_path = "output/french_meridional_timeseries.csv";
    const Real r_axis_hi_ts = Real(0.08);
    const Real r_rim_lo_ts = Real(0.18);

    Real last_adv_phys_time = time_stepper.getPhysicalTime();
    size_t spatial_q_resample_count = 0;
    /** Doc224 E4: abort instead of spinning after NaN / catastrophic U_max. */
    bool numerical_abort = false;
    std::string numerical_abort_reason;
    bool doc221_ledger_csv_header = true;
    const std::string doc221_ledger_csv = "output/doc221_force_ledger_stages.csv";
    if (local_cli.doc221_legacy_omit_clears_forceprior)
    {
        std::cout << "[ophelie][doc221] WARNING legacy omit clears ForcePrior each adv (Doc220 bug repro)\n";
    }
    else if (local_cli.omit_gravity_force)
    {
        std::cout << "[ophelie][doc221] omit_gravity=fixed_noop (ForcePrior not wiped)\n";
    }
    // Doc227 B: matched outer viscous/inertial budget for French keep (opt-in).
    Real doc227_outer_budget = Real(0);
    Real doc227_outer_elapsed = Real(0);
    Real doc227_outer_t0 = time_stepper.getPhysicalTime();
    Real doc227_viscous_force_last_t = time_stepper.getPhysicalTime();
    Real doc227_dt_visc_cand = dt_adv_viscous0;
    Real doc227_dt_ac_cand_last = Real(0);
    size_t doc227_position_updates = 0;
    size_t doc227_neighbor_updates = 0;
    size_t doc227_viscous_updates = 0;
    bool doc227_outer_need_init = local_cli.doc227_actual_timestep;
    bool doc227_audit_header = true;
    // Doc233 B: audit CSV on the *effective* path — default/legacy when --doc227-step-audit= is set
    // without --doc227-actual-timestep (do not grant PASS to default path from matched-path-only logs).
    const bool doc233_timestep_audit =
        local_cli.doc227_actual_timestep || !local_cli.doc227_step_audit_csv.empty();
    Real doc233_thermal_dt_last = Real(0);
    Real doc233_thermal_time_sum = Real(0);
    Real doc233_temp_refresh_last_t = time_stepper.getPhysicalTime();
    Real doc233_max_viscous_age = Real(0);
    Real doc233_max_dt_over_ac = Real(0);
    Real doc233_max_age_over_visc = Real(0);
    size_t doc233_n_dt_ac_violate = 0;
    size_t doc233_n_age_violate = 0;
    size_t doc233_n_end_time_mismatch = 0;
    // Doc234 C: count at most one age exceedance per viscous-hold interval (not every acoustic sample).
    bool doc234_age_flagged_this_hold = false;
    Real doc234_prev_viscous_last_t = -Real(1);
    const std::string doc227_audit_path =
        local_cli.doc227_step_audit_csv.empty() ? "output/doc227_actual_step_audit.csv"
                                               : local_cli.doc227_step_audit_csv;
    auto doc227_append_step_audit = [&](Real t_start, Real t_end, Real dt_actual, Real dt_ac_cand,
                                        Real outer_start, Real outer_budget, Real outer_elapsed,
                                        Real outer_remaining, Real dt_visc_cand, Real mu_raw_max,
                                        Real mu_used_max, Real thermal_dt, Real thermal_time_sum,
                                        Real temp_refresh_last_t, const char *schedule_path) {
        if (!doc233_timestep_audit)
            return;
        std::ofstream csv(doc227_audit_path, std::ios::app);
        if (doc227_audit_header)
        {
            csv << "t_start,t_end,dt_actual,dt_acoustic_candidate,"
                   "outer_start,outer_budget,outer_elapsed,outer_remaining,"
                   "dt_viscous_candidate,mu_raw_max,mu_used_max,"
                   "viscous_force_last_update_time,viscous_force_age,"
                   "age_minus_dt_visc_cand,age_exceed_hold,"
                   "position_update_count,neighbor_update_count,viscous_update_count,"
                   "thermal_dt,thermal_time_sum,temp_refresh_last_t,schedule_path\n";
            doc227_audit_header = false;
            std::cout << "[ophelie][doc234][C] audit defs: viscous_force_age=t_end-viscous_force_last_update_time; "
                         "dt_viscous_candidate=advection_viscous_time_step.exec(); "
                         "age_minus_dt_visc_cand=age-dt_visc (NOT age/dt); "
                         "legacy path does not truncate acoustic to remaining viscous budget; "
                         "n_age_violate counts once per hold interval when age>dt_visc+tol\n";
        }
        if (doc227_viscous_force_last_t != doc234_prev_viscous_last_t)
        {
            doc234_age_flagged_this_hold = false;
            doc234_prev_viscous_last_t = doc227_viscous_force_last_t;
        }
        const Real age = t_end - doc227_viscous_force_last_t;
        const Real age_minus = (dt_visc_cand > TinyReal) ? (age - dt_visc_cand) : Real(0);
        const bool age_exceed =
            dt_visc_cand > TinyReal && age > dt_visc_cand + doc227ScheduleTol(dt_visc_cand);
        csv << t_start << "," << t_end << "," << dt_actual << "," << dt_ac_cand << ","
            << outer_start << "," << outer_budget << "," << outer_elapsed << "," << outer_remaining
            << "," << dt_visc_cand << "," << mu_raw_max << "," << mu_used_max << ","
            << doc227_viscous_force_last_t << "," << age << "," << age_minus << ","
            << (age_exceed ? 1 : 0) << "," << doc227_position_updates << ","
            << doc227_neighbor_updates << "," << doc227_viscous_updates << "," << thermal_dt << ","
            << thermal_time_sum << "," << temp_refresh_last_t << "," << schedule_path << "\n";
        doc233_max_viscous_age = std::max(doc233_max_viscous_age, age);
        if (dt_ac_cand > TinyReal && dt_actual > dt_ac_cand + doc227ScheduleTol(dt_ac_cand))
        {
            ++doc233_n_dt_ac_violate;
            doc233_max_dt_over_ac =
                std::max(doc233_max_dt_over_ac, dt_actual - dt_ac_cand);
        }
        if (age_exceed)
        {
            doc233_max_age_over_visc = std::max(doc233_max_age_over_visc, age_minus);
            if (!doc234_age_flagged_this_hold)
            {
                ++doc233_n_age_violate;
                doc234_age_flagged_this_hold = true;
            }
        }
    };
    auto doc227_begin_outer_if_needed = [&]() {
        if (!local_cli.doc227_actual_timestep)
            return;
        if (!doc227_outer_need_init && doc227_outer_elapsed + TinyReal < doc227_outer_budget)
            return;
        advection_viscous_time_step.refreshViscousSpeedReference();
        const Real dt_in = advection_time_step.exec();
        const Real dt_vis = advection_viscous_time_step.exec();
        doc227_dt_visc_cand = dt_vis;
        const Real t_now = time_stepper.getPhysicalTime();
        doc227_outer_budget = std::min({dt_in, dt_vis, local_cli.end_time - t_now});
        if (doc227_outer_budget < TinyReal)
            doc227_outer_budget = TinyReal;
        doc227_outer_elapsed = Real(0);
        doc227_outer_t0 = t_now;
        doc227_outer_need_init = false;
        // Lock viscous ForcePrior for this outer interval (age starts here).
        fluid_advection_step_setup.exec();
        exec_viscous();
        apply_gravity_force_prior();
        if (local_cli.literature.enable_boussinesq)
            boussinesq_force.exec();
        doc227_viscous_force_last_t = t_now;
        ++doc227_viscous_updates;
    };

    auto apply_thermal_interval = [&](Real dt_th) {
        if (dt_th <= TinyReal)
            return;
        // Doc219 E: re-sample fixed spatial Q each thermal update (not particle-ID advection).
        if (local_cli.spatial_frozen_q && !local_cli.literature.zero_q)
        {
            size_t n_uncovered = 0;
            const bool verbose = (spatial_q_resample_count % 50 == 0);
            const Real p_sample = hostSampleSpatialFrozenQ(glass_particles, q_grid, p_joule,
                                                           /*renormalize=*/true, n_uncovered, verbose);
            ++spatial_q_resample_count;
            recordFrenchFluidHeatPower(coupling_state, glass_particles);
            if (!std::isfinite(p_sample) || !std::isfinite(coupling_state.last_power.p_heat_fluid))
            {
                numerical_abort = true;
                numerical_abort_reason = "non_finite_P_heat_or_P_sample_in_spatial_q";
                std::cerr << "[ophelie][doc224][numerical-abort] " << numerical_abort_reason
                          << " t=" << time_stepper.getPhysicalTime() << " P_sample=" << p_sample
                          << " P_heat_fluid=" << coupling_state.last_power.p_heat_fluid << std::endl;
                return;
            }
            if (n_uncovered > 0 && verbose)
            {
                std::cerr << "[ophelie][doc219][spatial-q] note: n_outside_grid=" << n_uncovered
                          << " clamped_to_edge (not Q=0) P_sample=" << p_sample << std::endl;
            }
        }
        if (!local_cli.literature.zero_q)
        {
            if (local_cli.literature.constant_thermal_rheology)
                apply_q_const.exec(dt_th);
            else
                apply_q_variable.exec(dt_th);
            coupling_state.ledger.accumulateQ(coupling_state.last_power.p_heat_fluid, dt_th);
        }
        if (refresh_mu_from_T)
            sync_mu_from_temperature("thermal_interval", /*log_particle=*/false);
        if (!local_cli.freeze_prescribed_temperature)
        {
            if (use_skull_shell_thermal)
            {
                glass_shell_thermal->exec(dt_th);
                wall_shell_thermal->exec(dt_th);
            }
            else if (laplace_temperature && apply_diffusion)
            {
                laplace_temperature->exec();
                apply_diffusion->exec(dt_th);
            }
            if (thermal_bc.enable_french_natural_bc && !local_cli.doc232_adiabatic_heat)
            {
                if (face_t_extrapolate)
                {
                    if (face_t_local)
                    {
                        const OphelieLocalFaceDTdz loc = hostFitLocalFaceBandDTdz(
                            glass_particles, z_bot_face, z_top_face, face_local_band, kOphelieTemperatureField);
                        face_dTdz_bottom_live = loc.bottom_ok ? loc.dTdz_bottom : Real(0);
                        face_dTdz_top_live = loc.top_ok ? loc.dTdz_top : Real(0);
                    }
                    else
                    {
                        const Real g =
                            hostFitMassWeightedTemperatureDTdz(glass_particles, kOphelieTemperatureField);
                        face_dTdz_bottom_live = g;
                        face_dTdz_top_live = g;
                    }
                }
                const Real p_boundary_now = boundary_loss_power.exec();
                coupling_state.ledger.p_loss_live = p_boundary_now;
                coupling_state.ledger.accumulateBoundary(p_boundary_now, dt_th);
                apply_natural_bc.exec(dt_th);
            }
        }
    };

    if (doc233_timestep_audit)
    {
        std::cout << "[ophelie][doc233][B] timestep_audit enabled schedule_path="
                  << (local_cli.doc227_actual_timestep ? "matched_outer" : "legacy_default")
                  << " audit=" << doc227_audit_path << std::endl;
        std::ofstream wipe(doc227_audit_path, std::ios::trunc);
        wipe.close();
        doc227_audit_header = true;
    }

    // Doc232 C: production ApplyFrozenQ path, thermal-only (no fluid keep / no BC).
    if (local_cli.doc232_adiabatic_heat)
    {
        const Real t0_gate = time_stepper.getPhysicalTime();
        const Real M = hostSumParticleMass(glass_particles);
        const Real E0 = hostThermalEnergyFromDeltaTMass(glass_particles, cp);
        const Real dT0 = hostMassWeightedMeanDeltaT(glass_particles);
        const Real P = local_cli.prescribed_heat_power_w;
        const Real dt_th = Real(1.0e-3);
        int heat_steps = 0;
        std::cout << "[ophelie][doc232][C] thermal-only keep start t=" << t0_gate << " end=" << local_cli.end_time
                  << " P=" << P << " M=" << M << " cp=" << cp << " dt=" << dt_th
                  << " operator=ApplyFrozenQToTemperatureCK"
#if SPHINXSYS_USE_SYCL
                  << " exec=SYCL_MainExecutionPolicy"
#else
                  << " exec=Host_MainExecutionPolicy"
#endif
                  << "\n";
        while (!time_stepper.isEndTime(local_cli.end_time))
        {
            const Real t_now = time_stepper.getPhysicalTime();
            Real dt_use = dt_th;
            if (t_now + dt_use > local_cli.end_time)
                dt_use = local_cli.end_time - t_now;
            if (dt_use <= TinyReal)
                break;
            time_stepper.incrementPhysicalTime(dt_use);
            apply_thermal_interval(dt_use);
            ++heat_steps;
        }
        const Real t_end = time_stepper.getPhysicalTime();
        const Real t_heat = t_end - t0_gate;
        const Real E_end = hostThermalEnergyFromDeltaTMass(glass_particles, cp);
        const Real DeltaH = E_end - E0;
        const Real E_Q = coupling_state.ledger.e_q_applied;
        const Real E_b = coupling_state.ledger.e_boundary_out;
        const Real dT_num = hostMassWeightedMeanDeltaT(glass_particles) - dT0;
        const Real dT_exact = (M * cp > TinyReal) ? (P * t_heat / (M * cp)) : Real(0);
        const Real dT_rel =
            std::abs(dT_exact) > TinyReal ? std::abs(dT_num - dT_exact) / std::abs(dT_exact) : Real(1);
        const Real energy_scale = std::max({std::abs(E_Q), std::abs(DeltaH), Real(1)});
        const Real energy_rel = std::abs(DeltaH - E_Q + E_b) / energy_scale;
        const Real abs_dT_err = std::abs(dT_num - dT_exact);
        const bool finite_ok =
            std::isfinite(static_cast<double>(dT_num)) && std::isfinite(static_cast<double>(energy_rel));
        const bool dT_ok = (dT_rel <= Real(0.01)) && (dT_num > TinyReal);
        const bool energy_ok = (energy_rel <= Real(0.01)) && (std::abs(E_b) <= Real(1.0e-6) * energy_scale);
        const bool time_ok = std::abs(t_heat - local_cli.end_time) <= Real(1.0e-9) + Real(1.0e-6) * local_cli.end_time;
        const bool pass = finite_ok && dT_ok && energy_ok && time_ok && (heat_steps > 0);
        std::cout << "[ophelie][doc232][C] dT_num=" << dT_num << " dT_exact=" << dT_exact
                  << " dT_rel=" << dT_rel << " abs_dT_err=" << abs_dT_err << " DeltaH=" << DeltaH
                  << " E_Q=" << E_Q << " E_boundary=" << E_b << " energy_rel=" << energy_rel
                  << " t_heat=" << t_heat << " steps=" << heat_steps << "\n";
        std::cout << "[ophelie][doc232][C] PRODUCTION_HEAT_SOURCE_VERIFIED=" << (pass ? 1 : 0)
                  << " (production ApplyFrozenQ; not hostApplyUniformHeat)\n";
        if (!local_cli.doc232_dump_particles_csv.empty())
        {
            dumpDoc232ParticleSnapshotCsv(glass_particles, local_cli.doc232_dump_particles_csv, t_end,
                                          probe_kernel_h);
        }
        if (!local_cli.summary_json.empty())
        {
            std::ofstream js(local_cli.summary_json);
            js << std::setprecision(10);
            js << "{\n"
               << "  \"case\": \"doc232_c_adiabatic_production_heat\",\n"
               << "  \"PRODUCTION_HEAT_SOURCE_VERIFIED\": " << (pass ? 1 : 0) << ",\n"
               << "  \"passed\": " << (pass ? 1 : 0) << ",\n"
               << "  \"heat_source_mode\": \"prescribed_uniform_volumetric\",\n"
               << "  \"em_status\": \"NOT_RUN\",\n"
               << "  \"em_solver_ok\": null,\n"
               << "  \"generator_power_w\": " << local_cli.literature.generator_power_w << ",\n"
               << "  \"operator\": \"ApplyFrozenQToTemperatureCK\",\n"
               << "  \"not_hostApplyUniformHeat\": 1,\n"
#if SPHINXSYS_USE_SYCL
               << "  \"execution\": \"SYCL_MainExecutionPolicy\",\n"
#else
               << "  \"execution\": \"Host_MainExecutionPolicy\",\n"
#endif
               << "  \"prescribed_heat_power_w\": " << P << ",\n"
               << "  \"end_time\": " << local_cli.end_time << ",\n"
               << "  \"t_heat_s\": " << t_heat << ",\n"
               << "  \"heat_steps\": " << heat_steps << ",\n"
               << "  \"dt_thermal_s\": " << dt_th << ",\n"
               << "  \"M_kg\": " << M << ",\n"
               << "  \"cp\": " << cp << ",\n"
               << "  \"rho0\": " << rho0 << ",\n"
               << "  \"dT_num_K\": " << dT_num << ",\n"
               << "  \"dT_exact_K\": " << dT_exact << ",\n"
               << "  \"dT_rel\": " << dT_rel << ",\n"
               << "  \"abs_dT_err_K\": " << abs_dT_err << ",\n"
               << "  \"DeltaH_J\": " << DeltaH << ",\n"
               << "  \"E_Q_J\": " << E_Q << ",\n"
               << "  \"E_boundary_J\": " << E_b << ",\n"
               << "  \"energy_rel\": " << energy_rel << ",\n"
               << "  \"bc_disabled\": 1,\n"
               << "  \"diffusion_disabled\": 1,\n"
               << "  \"fluid_keep_disabled\": 1,\n"
               << "  \"note\": \"independent ref DeltaT=P*t/(M*cp); DeltaH=sum(m*cp*dT); not EM\"\n"
               << "}\n";
            std::cout << "[ophelie] summary_json=" << local_cli.summary_json << "\n";
        }
        return pass ? 0 : 1;
    }

    StdVec<Real> doc234_dump_times;
    StdVec<char> doc234_dump_done;
    if (!local_cli.doc234_particle_field_csv_prefix.empty())
    {
        std::stringstream ss(local_cli.doc234_field_dump_times);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            if (tok.empty())
                continue;
            doc234_dump_times.push_back(static_cast<Real>(std::atof(tok.c_str())));
        }
        doc234_dump_done.assign(doc234_dump_times.size(), 0);
        dumpDoc234ParticleFieldCsv(glass_particles,
                                   local_cli.doc234_particle_field_csv_prefix + "_t0.csv",
                                   time_stepper.getPhysicalTime());
    }
    StdVec<Real> doc235_keep_dump_times;
    StdVec<char> doc235_keep_dump_done;
    bool doc235_first_acoustic_dumped = false;
    if (!local_cli.doc235_off_event_dump_dir.empty())
    {
        std::stringstream ss(local_cli.doc235_off_event_keep_times);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            if (tok.empty())
                continue;
            doc235_keep_dump_times.push_back(static_cast<Real>(std::atof(tok.c_str())));
        }
        doc235_keep_dump_done.assign(doc235_keep_dump_times.size(), 0);
    }

    while (!time_stepper.isEndTime(local_cli.end_time) && !numerical_abort)
    {
        doc227_begin_outer_if_needed();

        const Real t_step_start = time_stepper.getPhysicalTime();
        const Real dt_ac_cand = acoustic_time_step.exec();
        doc227_dt_ac_cand_last = dt_ac_cand;
        Real acoustic_dt = dt_ac_cand;
        if (local_cli.doc227_actual_timestep)
        {
            const Real rem_outer = std::max(Real(0), doc227_outer_budget - doc227_outer_elapsed);
            const Real rem_end = std::max(Real(0), local_cli.end_time - t_step_start);
            acoustic_dt = std::min({dt_ac_cand, rem_outer, rem_end});
            if (acoustic_dt <= TinyReal)
            {
                doc227_outer_need_init = true;
                continue;
            }
            // Hard gate: never exceed outer viscous/inertial budget (float tol only).
            const Real rem_tol = doc227ScheduleTol(doc227_outer_budget);
            if (acoustic_dt > rem_outer + rem_tol)
            {
                numerical_abort = true;
                numerical_abort_reason = "doc227_dt_exceeds_outer_remaining";
                break;
            }
            time_stepper.incrementPhysicalTime(acoustic_dt);
        }
        else
        {
            acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        }
        run_acoustic(acoustic_dt);

        if (!doc234_dump_times.empty())
        {
            const Real t_now_dump = time_stepper.getPhysicalTime();
            for (size_t id = 0; id < doc234_dump_times.size(); ++id)
            {
                if (doc234_dump_done[id])
                    continue;
                if (t_now_dump + TinyReal >= doc234_dump_times[id])
                {
                    std::ostringstream name;
                    name << local_cli.doc234_particle_field_csv_prefix << "_t"
                         << doc234_dump_times[id] << ".csv";
                    dumpDoc234ParticleFieldCsv(glass_particles, name.str(), t_now_dump);
                    doc234_dump_done[id] = 1;
                }
            }
        }
        if (!local_cli.doc235_off_event_dump_dir.empty())
        {
            const Real t_now_dump = time_stepper.getPhysicalTime();
            if (!doc235_first_acoustic_dumped)
            {
                dumpDoc235OffEventCsv(
                    glass_particles,
                    (fs::path(local_cli.doc235_off_event_dump_dir) / "event_05_after_first_acoustic.csv").string(),
                    "after_first_acoustic", t_now_dump, gc, french.glass_radius, z_bottom, z_free);
                doc235_first_acoustic_dumped = true;
            }
            for (size_t id = 0; id < doc235_keep_dump_times.size(); ++id)
            {
                if (doc235_keep_dump_done[id])
                    continue;
                if (t_now_dump + TinyReal >= doc235_keep_dump_times[id])
                {
                    std::ostringstream name;
                    name << "event_keep_t" << doc235_keep_dump_times[id] << ".csv";
                    dumpDoc235OffEventCsv(glass_particles,
                                          (fs::path(local_cli.doc235_off_event_dump_dir) / name.str()).string(),
                                          "keep_t", t_now_dump, gc, french.glass_radius, z_bottom, z_free);
                    doc235_keep_dump_done[id] = 1;
                }
            }
        }

        // Legacy path: thermal on every acoustic substep (pre–Doc219).
        // Doc227 matched path: thermal once per outer (after position), not per acoustic.
        if (!local_cli.doc219_hydro_profile && !local_cli.doc227_actual_timestep)
        {
            apply_thermal_interval(acoustic_dt);
            doc233_thermal_dt_last = acoustic_dt;
            doc233_thermal_time_sum += acoustic_dt;
            doc233_temp_refresh_last_t = time_stepper.getPhysicalTime();
        }

        // Refresh buoyancy each acoustic substep (cavity-validated).
        if (local_cli.literature.enable_boussinesq)
            boussinesq_force.exec();
        if (local_cli.diag_tip_soft_force)
        {
            // Doc221 N4 / Doc222 B: tip soft force channel into ForcePrior (no host Uz edit).
            hostApplyTipSoftForcePrior(glass_particles, z_free, Real(0.5) * french.dp,
                                       local_cli.diag_tip_soft_force_r, local_cli.diag_tip_soft_force_k,
                                       local_cli.diag_tip_soft_force_c, false, Real(0),
                                       local_cli.diag_tip_soft_force_legacy_accumulate);
        }
        else if (local_cli.diag_tip_fs_dashpot)
        {
            // Doc221 N5 / Doc222 B: tip dashpot only on FS-incomplete particles.
            hostApplyTipSoftForcePrior(glass_particles, z_free, Real(0.5) * french.dp,
                                       local_cli.diag_tip_soft_force_r, Real(0),
                                       local_cli.diag_tip_soft_force_c, /*fs_incomplete_only=*/true,
                                       rho0, local_cli.diag_tip_soft_force_legacy_accumulate);
        }

        if (local_cli.doc227_actual_timestep)
        {
            doc227_outer_elapsed += acoustic_dt;
            const Real t_step_end = time_stepper.getPhysicalTime();
            doc227_append_step_audit(t_step_start, t_step_end, acoustic_dt, dt_ac_cand, doc227_outer_t0,
                                     doc227_outer_budget, doc227_outer_elapsed,
                                     doc227_outer_budget - doc227_outer_elapsed, doc227_dt_visc_cand,
                                     rheology_stats.mu_raw_max, rheology_stats.mu_used_max,
                                     doc233_thermal_dt_last, doc233_thermal_time_sum,
                                     doc233_temp_refresh_last_t, "matched_outer");
            // Viscous force age ≤ outer budget; equality allowed (one truncated step = full Dt).
            // Age = t_end - t0 can sit ~1 ulp above budget after multi-acoustic fill (C2 abort).
            const Real age = t_step_end - doc227_viscous_force_last_t;
            const Real age_tol = doc227ScheduleTol(doc227_outer_budget);
            if (age > doc227_outer_budget + age_tol)
            {
                numerical_abort = true;
                numerical_abort_reason = "doc227_viscous_force_age_exceeds_budget";
                break;
            }
        }
        else if (doc233_timestep_audit)
        {
            // Doc233 B: record the currently effective legacy schedule (not matched-outer).
            advection_viscous_time_step.refreshViscousSpeedReference();
            const Real dt_in_leg = advection_time_step.exec();
            const Real dt_vis_leg = advection_viscous_time_step.exec();
            doc227_dt_visc_cand = dt_vis_leg;
            const Real t_step_end = time_stepper.getPhysicalTime();
            doc227_append_step_audit(t_step_start, t_step_end, acoustic_dt, dt_ac_cand, t_step_start,
                                     std::min(dt_in_leg, dt_vis_leg), acoustic_dt,
                                     std::min(dt_in_leg, dt_vis_leg) - acoustic_dt, dt_vis_leg,
                                     rheology_stats.mu_raw_max, rheology_stats.mu_used_max,
                                     doc233_thermal_dt_last, doc233_thermal_time_sum,
                                     doc233_temp_refresh_last_t, "legacy_default");
        }

        const bool do_outer_adv =
            local_cli.doc227_actual_timestep
                ? (doc227_outer_elapsed + doc227ScheduleTol(doc227_outer_budget) >=
                       doc227_outer_budget ||
                   time_stepper.isEndTime(local_cli.end_time))
                : advection_step(use_viscous_adv_cfl
                                     ? static_cast<BaseDynamics<Real> &>(advection_viscous_time_step)
                                     : static_cast<BaseDynamics<Real> &>(advection_time_step));
        if (do_outer_adv)
        {
            ++advection_steps;
            if (local_cli.doc227_actual_timestep)
            {
                ++doc227_position_updates;
                ++doc227_neighbor_updates;
            }
            if (local_cli.diag_zero_net_uz)
            {
                const Real uz_mean = hostSubtractMassWeightedMeanUz(glass_particles);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-zero-net-uz] adv=" << advection_steps
                              << " subtracted_Uz_mean=" << uz_mean
                              << " (volume-closure diagnostic; NOT literature open FS)\n";
                }
            }
            if (local_cli.diag_hard_lid)
            {
                // Doc221 M: thin band (0.5 dp) — only stop lid penetration; do not wipe near-surface upwelling.
                const size_t n_kill =
                    hostKillUpwardVelocityNearLid(glass_particles, z_free, Real(0.5) * french.dp);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-hard-lid] adv=" << advection_steps
                              << " n_kill_up=" << n_kill
                              << " z_lid=" << z_free
                              << " band=0.5*dp (kinematic hard lid; NOT literature open FS)\n";
                }
            }
            else if (local_cli.diag_soft_lid)
            {
                // Doc221 N: same thin-band kill-up as soft top, but never clamp z.
                const size_t n_kill =
                    hostKillUpwardVelocityNearLid(glass_particles, z_free, Real(0.5) * french.dp);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-soft-lid] adv=" << advection_steps
                              << " n_kill_up=" << n_kill << " z_lid=" << z_free
                              << " band=0.5*dp no_clamp (NOT literature open FS)\n";
                }
            }
            else if (local_cli.diag_soft_lid_tip)
            {
                // Doc221 N2/N3: tip-local scale upward Uz; never clamp z; never full-lid wipe.
                const size_t n_touch = hostScaleUpwardVelocityNearTip(
                    glass_particles, z_free, Real(0.5) * french.dp, local_cli.diag_soft_lid_tip_r,
                    local_cli.diag_soft_lid_tip_uz_scale);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-soft-lid-tip] adv=" << advection_steps
                              << " n_touch_up=" << n_touch << " z_lid=" << z_free
                              << " tip_r=" << local_cli.diag_soft_lid_tip_r
                              << " uz_scale=" << local_cli.diag_soft_lid_tip_uz_scale
                              << " band=0.5*dp no_clamp (NOT literature open FS)\n";
                }
            }
            else if (local_cli.diag_tip_soft_force)
            {
                // Doc222 B: log-only count (never re-apply ForcePrior here).
                const size_t n_touch = hostCountTipSoftForceTouches(
                    glass_particles, z_free, Real(0.5) * french.dp, local_cli.diag_tip_soft_force_r,
                    local_cli.diag_tip_soft_force_k, local_cli.diag_tip_soft_force_c);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-tip-soft-force] adv=" << advection_steps
                              << " n_touch=" << n_touch << " z_lid=" << z_free
                              << " tip_r=" << local_cli.diag_tip_soft_force_r
                              << " k=" << local_cli.diag_tip_soft_force_k
                              << " c=" << local_cli.diag_tip_soft_force_c
                              << " band=0.5*dp no_Uz_kill no_clamp count_only"
                              << " legacy_acc=" << (local_cli.diag_tip_soft_force_legacy_accumulate ? 1 : 0)
                              << " (NOT literature open FS)\n";
                }
            }
            else if (local_cli.diag_tip_fs_dashpot)
            {
                const size_t n_touch = hostCountTipSoftForceTouches(
                    glass_particles, z_free, Real(0.5) * french.dp, local_cli.diag_tip_soft_force_r,
                    Real(0), local_cli.diag_tip_soft_force_c, /*fs_incomplete_only=*/true, rho0);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-tip-fs-dashpot] adv=" << advection_steps
                              << " n_touch_fs=" << n_touch << " z_lid=" << z_free
                              << " tip_r=" << local_cli.diag_tip_soft_force_r
                              << " c=" << local_cli.diag_tip_soft_force_c
                              << " fs_incomplete_only=1 no_Uz_kill no_clamp count_only"
                              << " legacy_acc=" << (local_cli.diag_tip_soft_force_legacy_accumulate ? 1 : 0)
                              << " (NOT literature open FS)\n";
                }
            }
            fluid_update_particle_position.exec();
            if (local_cli.diag_hard_lid)
            {
                const size_t n_clamp = hostClampZToLid(glass_particles, z_free);
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc221][diag-hard-lid] adv=" << advection_steps
                              << " n_clamp_z=" << n_clamp << "\n";
                }
            }
            // Keep geometric freeze-T after motion (Doc77-C style); else residual advection scrambles drive.
            if (local_cli.freeze_prescribed_temperature && has_prescribed_t)
            {
                apply_prescribed_temperature(/*quiet=*/true);
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
            run_density_update();
            maybe_dens_sum_audit(time_stepper.getPhysicalTime(), advection_steps);
            if (local_cli.doc223_paper_flat_slip_top &&
                (advection_steps == 1 || advection_steps % 100 == 0))
            {
                hostLogDoc223PaperTopHydroAudit(
                    ("keep_adv" + std::to_string(advection_steps)).c_str(), glass_particles, rho0,
                    local_cli.c0, local_cli.gravity_g, gc, french.glass_radius, z_bottom, z_free, french.dp);
                if (doc223_e3_state.active)
                {
                    hostDoc223LogE3LidAudit(("keep_adv" + std::to_string(advection_steps)).c_str(),
                                           glass_particles, z_free, french.dp, doc223_e3_state.ke0,
                                           doc223_e3_state.z_com0);
                }
            }
            if (doc222_fs_deform_on && (advection_steps == 1 || advection_steps % 10 == 0))
            {
                doc222_fs_surf_last = hostMeasureDoc222FsSurfaceStats(glass_particles, gc, french.glass_radius,
                                                                      z_free, 12);
                hostDoc222AppendFsSurfaceCsv(doc222_fs_csv_path, time_stepper.getPhysicalTime(),
                                             static_cast<int>(advection_steps), doc222_fs_surf_last,
                                             doc222_fs_csv_header);
                doc222_fs_csv_header = false;
                if (advection_steps == 1 || advection_steps % 100 == 0)
                {
                    std::cout << "[ophelie][doc222][D2] adv=" << advection_steps
                              << " t=" << time_stepper.getPhysicalTime()
                              << " mode_amp=" << doc222_fs_surf_last.mode_amp
                              << " mode_amp0=" << doc222_fs_surf0.mode_amp
                              << " z_max=" << doc222_fs_surf_last.z_max
                              << " U_rms=" << doc222_fs_surf_last.u_rms << " KE=" << doc222_fs_surf_last.ke
                              << "\n";
                }
            }
            fluid_advection_step_setup.exec();
            // Doc227 matched path: viscous ForcePrior locked at outer START (not here).
            // Legacy trigger path keeps end-of-adv viscous refresh.
            if (!local_cli.doc227_actual_timestep)
            {
                exec_viscous();
                doc227_viscous_force_last_t = time_stepper.getPhysicalTime();
                ++doc227_viscous_updates;
                ++doc227_position_updates;
                ++doc227_neighbor_updates;
            }
            else
            {
                // Thermal once per outer with actual integrated span; then reopen budget.
                if (!local_cli.doc219_hydro_profile)
                {
                    apply_thermal_interval(doc227_outer_elapsed);
                    doc233_thermal_dt_last = doc227_outer_elapsed;
                    doc233_thermal_time_sum += doc227_outer_elapsed;
                    doc233_temp_refresh_last_t = time_stepper.getPhysicalTime();
                }
                else
                {
                    const Real t_now = time_stepper.getPhysicalTime();
                    const Real dt_adv = std::max(Real(0), t_now - last_adv_phys_time);
                    last_adv_phys_time = t_now;
                    apply_thermal_interval(dt_adv);
                    doc233_thermal_dt_last = dt_adv;
                    doc233_thermal_time_sum += dt_adv;
                    doc233_temp_refresh_last_t = t_now;
                }
                if (refresh_mu_from_T)
                    sync_mu_from_temperature("doc227_outer_end", /*log_particle=*/false);
                doc227_outer_need_init = true;
            }
            const bool doc221_probe_this =
                local_cli.doc221_force_ledger_probe_adv > 0 &&
                static_cast<int>(advection_steps) <= local_cli.doc221_force_ledger_probe_adv;
            if (doc221_probe_this)
            {
                const auto st_v = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("after_viscous_adv" + std::to_string(advection_steps), st_v);
                hostDoc221AppendForceLedgerCsv(doc221_ledger_csv, "after_viscous",
                                               static_cast<int>(advection_steps),
                                               time_stepper.getPhysicalTime(), st_v, doc221_ledger_csv_header);
                doc221_ledger_csv_header = false;
            }
            apply_gravity_force_prior();
            if (doc221_probe_this)
            {
                const auto st_g = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("after_omit_or_gravity_adv" + std::to_string(advection_steps), st_g);
                hostDoc221AppendForceLedgerCsv(doc221_ledger_csv, "after_omit_or_gravity",
                                               static_cast<int>(advection_steps),
                                               time_stepper.getPhysicalTime(), st_g, doc221_ledger_csv_header);
                doc221_ledger_csv_header = false;
            }
            if (local_cli.literature.enable_boussinesq)
            {
                boussinesq_force.exec();
            }
            if (local_cli.diag_tip_soft_force)
            {
                hostApplyTipSoftForcePrior(glass_particles, z_free, Real(0.5) * french.dp,
                                           local_cli.diag_tip_soft_force_r,
                                           local_cli.diag_tip_soft_force_k,
                                           local_cli.diag_tip_soft_force_c, false, Real(0),
                                           local_cli.diag_tip_soft_force_legacy_accumulate);
            }
            else if (local_cli.diag_tip_fs_dashpot)
            {
                hostApplyTipSoftForcePrior(glass_particles, z_free, Real(0.5) * french.dp,
                                           local_cli.diag_tip_soft_force_r, Real(0),
                                           local_cli.diag_tip_soft_force_c,
                                           /*fs_incomplete_only=*/true, rho0,
                                           local_cli.diag_tip_soft_force_legacy_accumulate);
            }
            if (doc221_probe_this)
            {
                const auto st_b = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("after_boussinesq_adv" + std::to_string(advection_steps), st_b);
                hostDoc221AppendForceLedgerCsv(doc221_ledger_csv, "after_boussinesq",
                                               static_cast<int>(advection_steps),
                                               time_stepper.getPhysicalTime(), st_b, doc221_ledger_csv_header);
                doc221_ledger_csv_header = false;
            }
            if (local_cli.doc221_force_lifecycle_selftest && advection_steps == 1)
            {
                // Doc221 C/D: ledger + repeat buoyancy + omit helper must not destroy other channels.
                const auto st0 = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("selftest_assembled", st0);
                const Real fp_before_omit = st0.fp_rms;
                const Real visc_before = st0.visc_rms;
                const Real buoy_before = st0.buoy_rms;
                apply_gravity_force_prior();
                const auto st_omit = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("selftest_after_omit_helper", st_omit);
                if (local_cli.literature.enable_boussinesq)
                    boussinesq_force.exec();
                const auto st_rep = hostDoc221ForceLedgerStats(glass_particles);
                logDoc221ForceLedger("selftest_after_repeat_boussinesq", st_rep);
                const bool ledger_ok = st0.residual_rel_rms <= Real(1.0e-5) && st0.residual_max <= Real(1.0e-4);
                const bool omit_preserves =
                    local_cli.doc221_legacy_omit_clears_forceprior
                        ? (st_omit.fp_rms < Real(1.0e-8))
                        : (std::abs(st_omit.fp_rms - fp_before_omit) <= Real(1.0e-6) * (fp_before_omit + Real(1)) &&
                           std::abs(st_omit.visc_rms - visc_before) <= Real(1.0e-6) * (visc_before + Real(1)));
                const bool repeat_ok =
                    !local_cli.literature.enable_boussinesq ||
                    (std::abs(st_rep.fp_rms - st_omit.fp_rms) <= Real(1.0e-5) * (st_omit.fp_rms + Real(1)) &&
                     std::abs(st_rep.buoy_rms - buoy_before) <= Real(1.0e-5) * (buoy_before + Real(1)));
                const bool pass = ledger_ok && omit_preserves && repeat_ok;
                std::cout << "[ophelie][doc221][selftest] ledger_ok=" << (ledger_ok ? 1 : 0)
                          << " omit_preserves=" << (omit_preserves ? 1 : 0) << " repeat_ok=" << (repeat_ok ? 1 : 0)
                          << " legacy_wipe=" << (local_cli.doc221_legacy_omit_clears_forceprior ? 1 : 0)
                          << " PASS=" << (pass ? 1 : 0) << std::endl;
                if (!local_cli.summary_json.empty())
                {
                    std::ofstream js(local_cli.summary_json);
                    js << "{\n"
                       << "  \"doc221_force_lifecycle_selftest\": 1,\n"
                       << "  \"legacy_omit_clears_forceprior\": "
                       << (local_cli.doc221_legacy_omit_clears_forceprior ? 1 : 0) << ",\n"
                       << "  \"ledger_ok\": " << (ledger_ok ? 1 : 0) << ",\n"
                       << "  \"omit_preserves\": " << (omit_preserves ? 1 : 0) << ",\n"
                       << "  \"repeat_ok\": " << (repeat_ok ? 1 : 0) << ",\n"
                       << "  \"residual_rel_rms\": " << st0.residual_rel_rms << ",\n"
                       << "  \"residual_max\": " << st0.residual_max << ",\n"
                       << "  \"fp_rms\": " << st0.fp_rms << ",\n"
                       << "  \"visc_rms\": " << st0.visc_rms << ",\n"
                       << "  \"buoy_rms\": " << st0.buoy_rms << ",\n"
                       << "  \"fp_rms_after_omit\": " << st_omit.fp_rms << ",\n"
                       << "  \"fp_rms_after_repeat_bouss\": " << st_rep.fp_rms << ",\n"
                       << "  \"passed\": " << (pass ? 1 : 0) << "\n"
                       << "}\n";
                }
                return pass ? 0 : 2;
            }
            // Doc219: B (+ TVF) after ForcePrior; thermal after neighborhood rebuild.
            if (local_cli.kernel_correction || local_cli.enable_tvf)
                fluid_linear_correction_matrix.exec();
            if (local_cli.enable_tvf)
            {
                kernel_grad_integral.exec();
                transport_correction.exec();
            }
            if (local_cli.doc219_hydro_profile && !local_cli.doc227_actual_timestep)
            {
                const Real t_now = time_stepper.getPhysicalTime();
                const Real dt_adv = std::max(Real(0), t_now - last_adv_phys_time);
                last_adv_phys_time = t_now;
                apply_thermal_interval(dt_adv);
                // Doc234 C3: doc219 path previously applied thermal but never accumulated audit sum.
                doc233_thermal_dt_last = dt_adv;
                doc233_thermal_time_sum += dt_adv;
                doc233_temp_refresh_last_t = t_now;
            }
            if (periodic_coupling && !local_cli.prescribed_heat_mode)
            {
                hostSampleEulerQOntoParticles(glass_particles, q_grid);
                if (local_cli.doc246_q_upper_kappa > TinyReal)
                {
                    const Real melt_h = Real(2) * french.glass_half_height;
                    const Real z_bot = french.glass_center[2] - french.glass_half_height;
                    applyDoc246QUpperControlled(glass_particles, local_cli.doc246_q_upper_kappa, z_bot, melt_h,
                                                local_cli.target_power_w, /*verbose=*/false);
                }
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
                    if (local_cli.doc246_q_upper_kappa > TinyReal)
                    {
                        const Real melt_h = Real(2) * french.glass_half_height;
                        const Real z_bot = french.glass_center[2] - french.glass_half_height;
                        applyDoc246QUpperControlled(glass_particles, local_cli.doc246_q_upper_kappa, z_bot, melt_h,
                                                    local_cli.target_power_w);
                    }
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
                {
                    char snap[128];
                    std::snprintf(snap, sizeof(snap), "output/french_T_rz_t%07.2f.csv",
                                  static_cast<double>(t_now));
                    writeFluidTemperatureRzCsv(snap, glass_particles, cyl_frame, kOphelieTemperatureField);
                }
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
                if (probes_ts_stream)
                {
                    appendProbesTimeseries(*probes_ts_stream, t_now, glass_particles, doc232_probes, probe_kernel_h);
                    ++probes_ts_samples;
                }
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
                // Doc224 heated French NC: absolute U_max≫0.05 m/s was a red-line for quiet smokes.
                // Doc230 isothermal (--zero-q) and Doc231 prescribed-heat client warm NC: use Ma-style
                // gate vs c0 (settle/acoustic residual ~0.05–0.2 m/s is not heated French runaway).
                // Explicit --doc224-u-blowup-gate= overrides both (client long demos).
                const bool soft_u_blowup =
                    local_cli.literature.zero_q || local_cli.prescribed_heat_mode;
                const bool gate_overridden = local_cli.doc224_u_blowup_gate > Real(0);
                const Real u_blowup_gate =
                    gate_overridden
                        ? local_cli.doc224_u_blowup_gate
                        : (soft_u_blowup ? Real(0.05) * local_cli.c0 : Real(0.05));
                if (!std::isfinite(u_now) || !std::isfinite(spatial.u_rms) || !std::isfinite(t_mean_now) ||
                    u_now > u_blowup_gate)
                {
                    numerical_abort = true;
                    numerical_abort_reason = !std::isfinite(u_now) || !std::isfinite(spatial.u_rms) ||
                                                    !std::isfinite(t_mean_now)
                                                ? "non_finite_U_or_T_at_monitor"
                                                : (gate_overridden ? "U_max_above_doc224_u_blowup_gate"
                                                                   : (soft_u_blowup
                                                                          ? "U_max_above_0p05_c0_blowup_gate"
                                                                          : "U_max_above_0p05_blowup_gate"));
                    std::cerr << "[ophelie][doc224][numerical-abort] " << numerical_abort_reason
                              << " t=" << t_now << " U_max=" << u_now << " U_rms=" << spatial.u_rms
                              << " gate=" << u_blowup_gate << std::endl;
                }
            }
        }
    }

    if (numerical_abort)
    {
        std::cerr << "[ophelie][doc224][numerical-abort] exiting early: " << numerical_abort_reason
                  << " t=" << time_stepper.getPhysicalTime() << std::endl;
        return 3;
    }

    if (!local_cli.doc232_dump_particles_csv.empty())
    {
        const Real t_dump = time_stepper.getPhysicalTime();
        dumpDoc232ParticleSnapshotCsv(glass_particles, local_cli.doc232_dump_particles_csv, t_dump,
                                      probe_kernel_h);
        const fs::path online_path =
            fs::path(local_cli.doc232_dump_particles_csv).parent_path() / "online_probe_samples.csv";
        std::ofstream online(online_path);
        online << std::setprecision(10);
        online << "t_s,probe_id,x,y,z,T_K,n_neighbors,support_h_m,weight_sum,nearest_dist_m,coverage_valid\n";
        syncFixedSpaceProbeFieldsToHost(glass_particles, kOphelieTemperatureField);
        if (doc232_probes.empty() && local_cli.client_nc_geometry)
            doc232_probes = makeDoc232ClientNcProbes(french);
        for (const auto &p : doc232_probes)
        {
            const auto s =
                hostSampleFixedSpaceTemperature(glass_particles, p.pos, probe_kernel_h, kOphelieTemperatureField);
            online << t_dump << "," << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << "," << s.T
                   << "," << s.n_neighbors << "," << probe_kernel_h << "," << s.weight_sum << "," << s.nearest_dist
                   << "," << (s.valid ? 1 : 0) << "\n";
        }
        std::cout << "[ophelie][doc232][B2] online_probe_samples=" << online_path << "\n";
    }

    const Real u_max = u_max_reduce.exec();
    const Real t_mean = t_mean_reduce.exec();
    const Real t_max = t_max_reduce.exec();
    const Real t_min = t_min_reduce.exec();
    const Real dT = t_max - t_min;
    const Real u_buoy = frenchBuoyancySpeedScale(local_cli.gravity_g, beta, dT, glass_h);
    const FrenchMeltSpatialStats spatial_end =
        hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
    if (doc222_fs_deform_on)
    {
        doc222_fs_surf_last =
            hostMeasureDoc222FsSurfaceStats(glass_particles, gc, french.glass_radius, z_free, 12);
        hostDoc222AppendFsSurfaceCsv(doc222_fs_csv_path, time_stepper.getPhysicalTime(),
                                     static_cast<int>(advection_steps), doc222_fs_surf_last, doc222_fs_csv_header);
        doc222_fs_csv_header = false;
        std::cout << "[ophelie][doc222][D2] end mode_amp=" << doc222_fs_surf_last.mode_amp
                  << " mode_amp0=" << doc222_fs_surf0.mode_amp
                  << " ratio="
                  << (std::abs(doc222_fs_surf0.mode_amp) > TinyReal
                          ? doc222_fs_surf_last.mode_amp / doc222_fs_surf0.mode_amp
                          : Real(0))
                  << " z_max=" << doc222_fs_surf_last.z_max << " U_rms=" << doc222_fs_surf_last.u_rms
                  << " KE=" << doc222_fs_surf_last.ke << "\n";
    }
    writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame, kOphelieTemperatureField);
    const GlassEscapeCensus escape = hostCountEscapedGlassParticles(glass_particles, gc, french.glass_radius, z_bottom,
                                                                   z_rim, french.dp);
    const Real escaped_frac =
        escape.total > 0 ? static_cast<Real>(escape.escaped) / static_cast<Real>(escape.total) : Real(0);
    bool doc223_e3_pass = false;
    Real doc223_e3_ke_end = 0, doc223_e3_dz_com = 0, doc223_e3_overshoot = 0;
    size_t doc223_e3_n_above = 0;
    if (doc223_e3_state.active)
    {
        hostDoc223LogE3LidAudit("end", glass_particles, z_free, french.dp, doc223_e3_state.ke0,
                                doc223_e3_state.z_com0);
        syncVariableToHost<Vecd>(glass_particles, "Position");
        const Vecd *pos_e3 = glass_particles.getVariableDataByName<Vecd>("Position");
        const size_t n_e3 = glass_particles.TotalRealParticles();
        for (size_t i = 0; i < n_e3; ++i)
        {
            if (pos_e3[i][2] > z_free)
            {
                ++doc223_e3_n_above;
                doc223_e3_overshoot = std::max(doc223_e3_overshoot, pos_e3[i][2] - z_free);
            }
        }
        doc223_e3_ke_end = hostDoc223KineticEnergy(glass_particles);
        doc223_e3_dz_com = hostDoc223MassComZ(glass_particles) - doc223_e3_state.z_com0;
        const Real ke_ratio =
            (doc223_e3_state.ke0 > TinyReal) ? (doc223_e3_ke_end / doc223_e3_state.ke0) : Real(0);
        const bool pass_escape = escaped_frac <= TinyReal && doc223_e3_n_above == 0;
        const bool pass_overshoot = doc223_e3_overshoot < Real(0.25) * french.dp;
        const bool pass_ke = ke_ratio < Real(2); // no strong self-excitation
        const bool pass_com = std::abs(doc223_e3_dz_com) < french.dp;
        const bool pass_ma = doc223_e3_state.ma < Real(0.05);
        doc223_e3_pass = pass_escape && pass_overshoot && pass_ke && pass_com && pass_ma;
        std::cout << "[ophelie][doc223][E3] verdict Ma=" << doc223_e3_state.ma
                  << " ke_ratio=" << ke_ratio << " dz_com=" << doc223_e3_dz_com
                  << " n_above=" << doc223_e3_n_above << " overshoot=" << doc223_e3_overshoot
                  << " escaped_frac=" << escaped_frac
                  << " pass_escape=" << (pass_escape ? 1 : 0)
                  << " pass_ke=" << (pass_ke ? 1 : 0) << " pass_com=" << (pass_com ? 1 : 0)
                  << " pass=" << (doc223_e3_pass ? 1 : 0) << "\n";
    }
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
    if (doc233_timestep_audit)
    {
        const Real end_tol = kFrenchValidationEndTimeTol;
        if (std::abs(physical_time - local_cli.end_time) > end_tol &&
            physical_time + end_tol < local_cli.end_time)
            ++doc233_n_end_time_mismatch;
        const Real thermal_gap =
            std::abs(doc233_thermal_time_sum - std::min(physical_time, local_cli.end_time));
        const bool thermal_accounted =
            thermal_gap <= end_tol + Real(1.0e-3) * std::max(local_cli.end_time, Real(1));
        // Doc234: n_age_violate is once-per-hold; legacy may still mildly exceed dt_visc_cand
        // because acoustic is not truncated to remaining viscous budget (see matched_outer).
        const bool audited_ok =
            (doc233_n_dt_ac_violate == 0) && (doc233_n_age_violate == 0) &&
            (doc233_n_end_time_mismatch == 0) && thermal_accounted;
        std::cout << "[ophelie][doc233][B] ACTUAL_TIMESTEP_AUDITED=" << (audited_ok ? 1 : 0)
                  << " schedule_path="
                  << (local_cli.doc227_actual_timestep ? "matched_outer" : "legacy_default")
                  << " audit_csv=" << doc227_audit_path
                  << " max_viscous_age=" << doc233_max_viscous_age
                  << " max_dt_over_ac=" << doc233_max_dt_over_ac
                  << " max_age_minus_dt_visc_cand=" << doc233_max_age_over_visc
                  << " n_dt_ac_violate=" << doc233_n_dt_ac_violate
                  << " n_age_violate_per_hold=" << doc233_n_age_violate
                  << " n_end_time_mismatch=" << doc233_n_end_time_mismatch
                  << " thermal_time_sum=" << doc233_thermal_time_sum
                  << " thermal_gap=" << thermal_gap
                  << " THERMAL_TIME_ACCOUNTED=" << (thermal_accounted ? 1 : 0) << "\n";
    }
    const Real energy_end = local_cli.literature.constant_thermal_rheology
                                ? hostThermalEnergy(glass_particles, rho0, cp, t0)
                                : hostThermalEnthalpyEnergy(glass_particles, rho0, laws.cp, t0);
    // Doc233 E: mass enthalpy vs rho0*V operator ledger side-by-side.
    const Real mass_sum_e = hostSumParticleMass(glass_particles);
    const Real H_mass_end = hostThermalEnergyFromDeltaTMass(glass_particles, cp);
    Real sum_V = Real(0);
    Real rho_min = Real(0), rho_max = Real(0);
    {
        syncVariableToHost<Real>(glass_particles, "VolumetricMeasure");
        syncVariableToHost<Real>(glass_particles, "Density");
        const Real *vol_e = glass_particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Real *dens = glass_particles.getVariableDataByName<Real>("Density");
        const size_t n_e = glass_particles.TotalRealParticles();
        rho_min = dens[0];
        rho_max = dens[0];
        for (size_t i = 0; i < n_e; ++i)
        {
            sum_V += vol_e[i];
            rho_min = std::min(rho_min, dens[i]);
            rho_max = std::max(rho_max, dens[i]);
        }
    }
    const Real H_rhoV_end = energy_end; // hostThermalEnergy uses rho0*V*cp*dT
    const Real mass_vs_rhoV_rel =
        std::abs(H_mass_end - H_rhoV_end) /
        std::max({std::abs(H_mass_end), std::abs(H_rhoV_end), Real(1)});
    // Doc234 E1: explicit dual-ledger residuals (do not confuse rho0V residual with mass enthalpy).
    const Real E_source_J = coupling_state.ledger.e_q_applied;
    const Real E_boundary_out_J = coupling_state.ledger.e_boundary_out;
    const Real residual_mass_J = H_mass_end - (E_source_J - E_boundary_out_J);
    const Real residual_rho0V_J = H_rhoV_end - (E_source_J - E_boundary_out_J);
    const Real energy_norm =
        std::max({std::abs(E_source_J), std::abs(E_boundary_out_J), Real(1)});
    const Real residual_mass_rel = residual_mass_J / energy_norm;
    const Real residual_rho0V_rel = residual_rho0V_J / energy_norm;
    std::cout << "[ophelie][doc234][E] MASS_ENTHALPY_ACCOUNTING"
              << " sum_m=" << mass_sum_e << " sum_V=" << sum_V << " rho0*sum_V=" << (rho0 * sum_V)
              << " rho_range=[" << rho_min << "," << rho_max << "]"
              << " deltaH_mass_J=" << H_mass_end << " deltaH_rho0V_J=" << H_rhoV_end
              << " mass_vs_rhoV_rel=" << mass_vs_rhoV_rel
              << " E_source_J=" << E_source_J << " E_boundary_out_J=" << E_boundary_out_J
              << " residual_mass_J=" << residual_mass_J << " residual_mass_rel=" << residual_mass_rel
              << " residual_rho0V_J=" << residual_rho0V_J << " residual_rho0V_rel=" << residual_rho0V_rel
              << " residual_norm=max(|E_source|,|E_boundary|,1)"
              << " legacy_ledger_residualH_rel=" << coupling_state.ledger.residualHRel(energy_end) << "\n";
    // Doc246 §2.2: residual_mass_J above ignores H0 (IC enthalpy after prescribe). Correct gate uses H0.
    {
        const Real H0_ic = static_cast<Real>(coupling_state.ledger.h0);
        const Real dH_from_H0 = energy_end - H0_ic;
        const Real residual_vs_H0 = coupling_state.ledger.residualH(energy_end);
        std::cout << "[ophelie][doc246][energy] H0_set=" << (coupling_state.ledger.h0_set ? 1 : 0)
                  << " H0=" << H0_ic << " H_end=" << energy_end << " dH_from_H0=" << dH_from_H0
                  << " E_Q=" << E_source_J << " E_boundary_out=" << E_boundary_out_J
                  << " residual_vs_H0=" << residual_vs_H0
                  << " residual_vs_H0_rel=" << coupling_state.ledger.residualHRel(energy_end)
                  << " (doc234 residual_mass_J ignores H0 — do not treat as operator failure)\n";
    }    const bool periodic_coupling_end = periodic_coupling;
    FrenchCouplingGateInputs gate_in;
    gate_in.periodic_coupling = periodic_coupling_end && !local_cli.prescribed_heat_mode;
    gate_in.zero_q = local_cli.literature.zero_q;
    gate_in.prescribed_heat = local_cli.prescribed_heat_mode;
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
    const char *heat_source_mode = local_cli.prescribed_heat_mode
                                       ? "prescribed_uniform_volumetric"
                                       : (local_cli.literature.zero_q
                                              ? "zero_q"
                                              : (local_cli.enable_target_power ? "em_target_power"
                                                                              : "em_or_frozen_q"));
    const char *em_status = local_cli.prescribed_heat_mode
                                ? "NOT_RUN"
                                : (local_cli.literature.zero_q ? "DISABLED" : "RUN");

    std::cout << "test_3d_ophelie_french_natural_convection_frozen_q"
              << " heat_source_mode=" << heat_source_mode << " em_status=" << em_status
              << " prescribed_heat_power_w=" << local_cli.prescribed_heat_power_w
              << " generator_power_w=" << local_cli.literature.generator_power_w
              << " P_heat_fluid_W=" << coupling_state.last_power.p_heat_fluid
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
              << " em_solver_ok="
              << ((local_cli.prescribed_heat_mode || local_cli.literature.zero_q)
                      ? "null"
                      : (gates.em_solver_ok ? "1" : "0"))
              << " q_transfer_ok=" << (gates.q_transfer_ok ? 1 : 0)
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
    // Doc224 B1: fixed 5 s windows; extend list when end_time grows (20–25 / 25–30).
    std::vector<std::pair<Real, Real>> default_windows = {
        {Real(0), Real(5)}, {Real(5), Real(10)}, {Real(10), Real(15)}, {Real(15), Real(20)}};
    if (local_cli.end_time + Real(1e-9) >= Real(25))
        default_windows.push_back({Real(20), Real(25)});
    if (local_cli.end_time + Real(1e-9) >= Real(30))
        default_windows.push_back({Real(25), Real(30)});
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
           << "  \"heat_source_mode\": \"" << heat_source_mode << "\",\n"
           << "  \"em_status\": \"" << em_status << "\",\n"
           << "  \"prescribed_heat_power_w\": " << local_cli.prescribed_heat_power_w << ",\n"
           << "  \"generator_power_w\": " << local_cli.literature.generator_power_w << ",\n"
           << "  \"P_heat_fluid_W\": " << coupling_state.last_power.p_heat_fluid << ",\n"
           << "  \"P_joule_W\": " << p_joule << ",\n"
           << "  \"target_power_w\": " << local_cli.target_power_w << ",\n"
           << "  \"em_solver_ok\": "
           << ((local_cli.prescribed_heat_mode || local_cli.literature.zero_q)
                   ? "null"
                   : (gates.em_solver_ok ? "1" : "0"))
           << ",\n"
           << "  \"doc233_timestep_audit\": " << (doc233_timestep_audit ? 1 : 0) << ",\n"
           << "  \"ACTUAL_TIMESTEP_AUDITED\": "
           << (doc233_timestep_audit
                   ? ((doc233_n_dt_ac_violate == 0 && doc233_n_age_violate == 0 &&
                       doc233_n_end_time_mismatch == 0)
                          ? "1"
                          : "0")
                   : "null")
           << ",\n"
           << "  \"doc233_schedule_path\": \""
           << (local_cli.doc227_actual_timestep ? "matched_outer" : "legacy_default") << "\",\n"
           << "  \"doc233_max_viscous_age\": " << doc233_max_viscous_age << ",\n"
           << "  \"doc233_n_dt_ac_violate\": " << doc233_n_dt_ac_violate << ",\n"
           << "  \"doc233_n_age_violate\": " << doc233_n_age_violate << ",\n"
           << "  \"doc233_thermal_time_sum\": " << doc233_thermal_time_sum << ",\n"
           << "  \"dp\": " << french.dp << ",\n"
           << "  \"end_time\": " << local_cli.end_time << ",\n"
           << "  \"beta\": " << beta << ",\n"
           << "  \"requested_beta\": " << requested_beta << ",\n"
           << "  \"effective_beta\": " << effective_beta << ",\n"
           << "  \"boussinesq\": " << (local_cli.literature.enable_boussinesq ? 1 : 0) << ",\n"
           << "  \"freeze_prescribed_t\": " << (local_cli.freeze_prescribed_temperature ? 1 : 0) << ",\n"
           << "  \"doc232_probes_timeseries\": \"" << local_cli.probes_timeseries_csv << "\",\n"
           << "  \"doc232_probe_samples\": " << probes_ts_samples << ",\n"
           << "  \"doc232_probe_path_verify\": " << (local_cli.doc232_probe_path_verify ? 1 : 0) << ",\n"
           << "  \"PRODUCTION_PROBE_PATH_VERIFY\": " << (doc232_probe_path_ok ? 1 : 0) << ",\n"
           << "  \"diag_zero_net_uz\": " << (local_cli.diag_zero_net_uz ? 1 : 0) << ",\n"
           << "  \"diag_hard_lid\": " << (local_cli.diag_hard_lid ? 1 : 0) << ",\n"
           << "  \"diag_soft_lid\": " << (local_cli.diag_soft_lid ? 1 : 0) << ",\n"
           << "  \"diag_soft_lid_tip\": " << (local_cli.diag_soft_lid_tip ? 1 : 0) << ",\n"
           << "  \"diag_soft_lid_tip_r\": " << local_cli.diag_soft_lid_tip_r << ",\n"
           << "  \"diag_soft_lid_tip_uz_scale\": " << local_cli.diag_soft_lid_tip_uz_scale << ",\n"
           << "  \"diag_tip_soft_force\": " << (local_cli.diag_tip_soft_force ? 1 : 0) << ",\n"
           << "  \"diag_tip_fs_dashpot\": " << (local_cli.diag_tip_fs_dashpot ? 1 : 0) << ",\n"
           << "  \"diag_tip_soft_force_r\": " << local_cli.diag_tip_soft_force_r << ",\n"
           << "  \"diag_tip_soft_force_k\": " << local_cli.diag_tip_soft_force_k << ",\n"
           << "  \"diag_tip_soft_force_c\": " << local_cli.diag_tip_soft_force_c << ",\n"
           << "  \"diag_wall_flat_lid\": " << (local_cli.diag_wall_flat_lid ? 1 : 0) << ",\n"
           << "  \"diag_lid_free_slip\": " << (local_cli.diag_lid_free_slip ? 1 : 0) << ",\n"
           << "  \"diag_lid_free_slip_post_settle\": "
           << (local_cli.diag_lid_free_slip_post_settle ? 1 : 0) << ",\n"
           << "  \"diag_wall_lid_viscous_omit\": " << (local_cli.diag_wall_lid_viscous_omit ? 1 : 0)
           << ",\n"
           << "  \"diag_lid_dens_silence\": " << (local_cli.diag_lid_dens_silence ? 1 : 0) << ",\n"
           << "  \"diag_lid_cont_silent\": " << (local_cli.diag_lid_cont_silent ? 1 : 0) << ",\n"
           << "  \"diag_lid_cont_silent_only\": " << (local_cli.diag_lid_cont_silent_only ? 1 : 0)
           << ",\n"
           << "  \"diag_lid_omit_pjump_only\": " << (local_cli.diag_lid_omit_pjump_only ? 1 : 0)
           << ",\n"
           << "  \"prescribe_t_radial_dT\": " << local_cli.prescribe_t_radial_dT << ",\n"
           << "  \"init_T_rz_csv\": \"" << local_cli.init_T_rz_csv << "\",\n"
           << "  \"init_T_rz_demean_legacy\": " << (local_cli.init_T_rz_demean_legacy ? 1 : 0) << ",\n"
           << "  \"t_force_ref_from_mass_mean\": " << (local_cli.t_force_ref_from_mass_mean ? 1 : 0) << ",\n"
           << "  \"doc223_f2_annular_hot\": " << (local_cli.doc223_f2_annular_hot ? 1 : 0) << ",\n"
           << "  \"init_T\": \""
           << (!local_cli.init_T_rz_csv.empty()
                   ? (local_cli.freeze_prescribed_temperature ? "prescribed_doc224_fig8_rz_csv_frozen_ASSUMED"
                                                             : "prescribed_doc224_fig8_rz_csv_then_thaw_ASSUMED")
                   : (local_cli.doc223_f2_annular_hot
                          ? (local_cli.freeze_prescribed_temperature ? "prescribed_doc223_f2_annular_hot_frozen"
                                                                    : "prescribed_doc223_f2_annular_hot_then_thaw")
                          : (local_cli.prescribe_t_radial_dT > TinyReal
                                 ? (local_cli.freeze_prescribed_temperature
                                        ? "prescribed_radial_hot_core_frozen"
                                        : "prescribed_radial_hot_core_then_thaw")
                                 : (local_cli.prescribe_t_hot_bottom_dT > TinyReal
                                        ? (local_cli.freeze_prescribed_temperature
                                               ? "prescribed_hot_bottom_frozen"
                                               : "prescribed_hot_bottom_then_thaw")
                                        : "isothermal_T0"))))
           << "\",\n"
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
           << "  \"french_pressure_path\": \""
           << (local_cli.doc220_pprime_hydro
                   ? "p_prime_plus_BackgroundPressure_no_GravityForceCK"
                   : (local_cli.background_hydro_p ? "total_p_plus_BackgroundPressure_plus_GravityForceCK"
                                                   : "total_p_plus_GravityForceCK"))
           << "\",\n"
           << "  \"perturbation_momentum\": " << local_cli.perturbation_momentum << ",\n"
           << "  \"doc220_pprime_hydro\": " << (local_cli.doc220_pprime_hydro ? 1 : 0) << ",\n"
           << "  \"omit_gravity_force\": " << (local_cli.omit_gravity_force ? 1 : 0) << ",\n"
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
           << "  \"doc222_ph_signed\": " << (local_cli.doc222_ph_signed ? 1 : 0) << ",\n"
           << "  \"doc222_deform_skip_rho_freeze\": " << (local_cli.doc222_deform_skip_rho_freeze ? 1 : 0)
           << ",\n"
           << "  \"doc222_ph_local_eta\": " << (local_cli.doc222_ph_local_eta ? 1 : 0) << ",\n"
           << "  \"doc222_ph_local_eta_bins\": " << local_cli.doc222_ph_local_eta_bins << ",\n"
           << "  \"doc222_ph_local_eta_spread\": " << doc222_local_eta_ph_last.eta_spread << ",\n"
           << "  \"doc222_ph_local_eta_mean\": " << doc222_local_eta_ph_last.eta_mean << ",\n"
           << "  \"doc222_fs_deform_eps\": " << local_cli.doc222_fs_deform_eps << ",\n"
           << "  \"doc222_fs_deform_on\": " << (doc222_fs_deform_on ? 1 : 0) << ",\n"
           << "  \"doc222_fs_mode_amp0\": " << doc222_fs_surf0.mode_amp << ",\n"
           << "  \"doc222_fs_mode_amp_end\": " << doc222_fs_surf_last.mode_amp << ",\n"
           << "  \"doc222_fs_z_max0\": " << doc222_fs_surf0.z_max << ",\n"
           << "  \"doc222_fs_z_max_end\": " << doc222_fs_surf_last.z_max << ",\n"
           << "  \"doc222_fs_ke0\": " << doc222_fs_surf0.ke << ",\n"
           << "  \"doc222_fs_ke_end\": " << doc222_fs_surf_last.ke << ",\n"
           << "  \"doc222_fs_eps_over_dp\": " << doc222_fs_deform_rep.eps_over_dp << ",\n"
           << "  \"doc222_fs_eps_over_H\": " << doc222_fs_deform_rep.eps_over_H << ",\n"
           << "  \"hydro_model\": \""
           << (local_cli.doc223_paper_flat_slip_top
                   ? "paper_flat_slip_top"
                   : (local_cli.diag_wall_flat_lid ? "closed_crucible_diag" : "open_free_surface"))
           << "\",\n"
           << "  \"pressure_model\": \""
           << (local_cli.doc223_paper_flat_slip_top ? "perturbation_fixed_top"
                                                    : (local_cli.doc220_pprime_hydro ? "pprime_open_fs"
                                                                                     : "total_p"))
           << "\",\n"
           << "  \"doc223_paper_flat_slip_top\": " << (local_cli.doc223_paper_flat_slip_top ? 1 : 0)
           << ",\n"
           << "  \"doc223_ablate_no_free_slip\": " << (local_cli.doc223_ablate_no_free_slip ? 1 : 0)
           << ",\n"
           << "  \"doc224_lid_keep_dens_jump\": " << (local_cli.doc224_lid_keep_dens_jump ? 1 : 0)
           << ",\n"
           << "  \"doc224_ablate_opencup\": " << (local_cli.doc224_ablate_opencup ? 1 : 0) << ",\n"
           << "  \"doc224_ablate_acoustic_freeslip\": "
           << (local_cli.doc224_ablate_acoustic_freeslip ? 1 : 0) << ",\n"
           << "  \"doc223_e3_normal_perturb\": " << (doc223_e3_state.active ? 1 : 0) << ",\n"
           << "  \"doc223_e3_Ma\": " << doc223_e3_state.ma << ",\n"
           << "  \"doc223_e3_ke0\": " << doc223_e3_state.ke0 << ",\n"
           << "  \"doc223_e3_ke_end\": " << doc223_e3_ke_end << ",\n"
           << "  \"doc223_e3_dz_com\": " << doc223_e3_dz_com << ",\n"
           << "  \"doc223_e3_n_above_lid\": " << doc223_e3_n_above << ",\n"
           << "  \"doc223_e3_overshoot_max\": " << doc223_e3_overshoot << ",\n"
           << "  \"doc223_e3_pass\": "
           << (doc223_e3_state.active ? (doc223_e3_pass ? "1" : "0") : "null") << ",\n"
           << "  \"smoke_passed\": " << (gates.smoke_passed ? 1 : 0) << ",\n"
           << "  \"passed\": " << (passed ? 1 : 0) << ",\n"
           << "  \"hydrostatic_pass\": "
           << ((local_cli.doc223_paper_flat_slip_top && std::abs(local_cli.beta) <= TinyReal &&
                !doc223_e3_state.active)
                   ? ((spatial_end.u_rms < Real(1e-4) && escaped_frac <= TinyReal) ? "1" : "0")
                   : "null")
           << ",\n"
           << "  \"meridional_timeseries\": \"" << meridional_ts_path << "\",\n";
        writeFrenchMeridionalWindowsJsonArray(js, meridional_windows);
        js << ",\n"
           << "  \"note\": \"Doc97 P3: lit-β SNR + non-overlap windows; not u_max gate; signal_verified OPEN until repeats\"\n"
           << "}\n";
        std::cout << "[ophelie] summary_json=" << local_cli.summary_json << std::endl;
    }
    if (doc223_e3_state.active)
        return (passed && doc223_e3_pass) ? 0 : 1;
    return passed ? 0 : 1;
}
