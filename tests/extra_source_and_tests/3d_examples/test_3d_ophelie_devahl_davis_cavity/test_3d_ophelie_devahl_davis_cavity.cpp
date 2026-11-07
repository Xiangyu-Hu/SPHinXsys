/**
 * @file test_3d_ophelie_devahl_davis_cavity.cpp
 * @brief Doc77-C / Doc84-P0: differentially heated closed cavity (de Vahl Davis path).
 *
 * Geometry: square (x,z), extruded in y (3D prism). Left hot / right cold; top/bottom
 * adiabatic fluid; all solid walls no-slip. Gravity −z.
 *
 * Doc84 P0 changes:
 *  - Default bw=0: wall-particle Dirichlet via Contact Laplace (fluid shell optional).
 *  - Primary Nu: constraint heat (fluid shell) or Contact wall heat (bw=0) / (A_wall·Δt).
 *  - Circulation bands fixed, non-overlapping, independent of Nu metric.
 *  - Thermal substeps fill actual advance time; freeze-velocity skips fluid density/acoustic.
 *
 * Label: 3D prism / mid-y diagnostic — NOT pure 2D benchmark until documented.
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
struct Cli
{
    Real L = 0.2;
    Real depth_factor = 0.5; /**< Ly = depth_factor * L */
    Real dp = 0.025;
    Real Ra = 1000.0;
    Real Pr = 0.71;
    Real rho0 = 1.0;
    Real cp = 1000.0;
    Real gravity_g = 9.81;
    Real beta = 0.01;
    Real dT = 10.0;
    Real t_c = 300.0;
    Real settle_time = 2.0;
    Real end_time = 12.0;
    Real c0_factor = 50.0;
    Real mu_floor = 1.0e-4;
    /**
     * Fluid-side Dirichlet shell in units of dp. Doc84-P0+: default 0 → wall-particle BC
     * via Contact Laplace (left/right Th/Tc; other faces k=0 adiabatic).
     */
    Real boundary_width_factor = 0.0;
    /** Diagnostic grad-Nu layer = factor × shell (capped to avoid overlap). */
    Real nu_probe_layer_factor = 1.0;
    /** Circulation band width as fraction of L (hot/cold, non-overlapping). */
    Real circ_band_factor = 0.1;
    Real report_every = 1.0;
    bool enable_diffusion = false;
    bool enable_boussinesq = true;
    bool freeze_velocity = false;
    /**
     * Doc87 R1-A: single-hot-wall energy ledger (freeze-vel, no buoyancy).
     * Left wall Dirichlet Th; all other faces adiabatic (k=0). Uniform cold init.
     */
    bool energy_audit = false;
    /** Scaled Contact heat for diagnostic Nu only; formal Nu / energy pass use raw. */
    Real contact_pair_factor = 2.0;
    /** Doc93: Adami Dirichlet dummy φ_j=2φ_bc−φ_i on Contact Laplace. Default on. */
    bool dirichlet_ghost = true;
    /** Wall k / k_fluid; harmonicMean(k_f, k_w→∞)→2 k_f. Default 1e6 (legacy Dirichlet drive). */
    Real k_wall_factor = Real(1.0e6);
    /** Dummy wall thickness in units of dp (ClosedWall BW). Default 3. */
    Real wall_bw_factor = Real(3.0);
    /**
     * Doc197 corner-rim: on Dirichlet side walls, set k=0 within this many dp of
     * z top/bottom (DVD adiabatic). 0 = off. Optional y via --corner-rim-y.
     */
    Real corner_rim_dp = Real(0.0);
    bool corner_rim_y = false;
    /** k multiplier on rim Dirichlet dummies (1=full; 0=off). Default 0 when rim on. */
    Real corner_rim_k_factor = Real(0.0);
    /** Laplace pair weight: 2 r dW V / (r² + ε h²). Match class/1D default 0.01 (Doc93). */
    Real pair_weight_regularization = Real(0.01);
    /**
     * Doc91/92/105 production hydro path (Doc87 OPEN hydrostatic): dens-sum packing
     * normalize + analytic p_h in BackgroundPressure + LinearCorrectionCK.
     * When background_hydro_p: PerturbationMomentum + AcousticDensMode=1, omit GravityForceCK
     * (match stratified hydrostatic / French closed-box). Legacy = total-p + GravityForceCK.
     * Freeze-velocity conduction skips fluid operators regardless of these flags.
     */
    bool dens_sum_normalize = true;
    bool background_hydro_p = true;
    bool kernel_correction = true;
    /** Particle sort period in fluid steps (Doc197 spike A/B). Default 50. */
    int sort_every = 50;
    std::string out_csv = "output/doc77c_devahl_davis.csv";
    std::string summary_json;
};

bool parseCli(int argc, char *argv[], Cli &cli)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto take = [&](const std::string &key, Real &dst) {
            const std::string p = key + "=";
            if (arg.rfind(p, 0) == 0)
            {
                dst = static_cast<Real>(std::atof(arg.c_str() + p.size()));
                return true;
            }
            return false;
        };
        auto take_s = [&](const std::string &key, std::string &dst) {
            const std::string p = key + "=";
            if (arg.rfind(p, 0) == 0)
            {
                dst = arg.substr(p.size());
                return true;
            }
            return false;
        };
        if (take("--L", cli.L) || take("--depth-factor", cli.depth_factor) || take("--dp", cli.dp) ||
            take("--Ra", cli.Ra) || take("--Pr", cli.Pr) || take("--rho0", cli.rho0) || take("--cp", cli.cp) ||
            take("--g", cli.gravity_g) || take("--beta", cli.beta) || take("--dT", cli.dT) ||
            take("--Tc", cli.t_c) || take("--settle-time", cli.settle_time) || take("--end-time", cli.end_time) ||
            take("--c0-factor", cli.c0_factor) || take("--mu-floor", cli.mu_floor) ||
            take("--boundary-width-factor", cli.boundary_width_factor) ||
            take("--nu-probe-layer-factor", cli.nu_probe_layer_factor) ||
            take("--circ-band-factor", cli.circ_band_factor) || take("--report-every", cli.report_every) ||
            take("--contact-pair-factor", cli.contact_pair_factor) ||
            take("--k-wall-factor", cli.k_wall_factor) ||
            take("--wall-bw-factor", cli.wall_bw_factor) ||
            take("--corner-rim-dp", cli.corner_rim_dp) ||
            take("--corner-rim-k-factor", cli.corner_rim_k_factor) ||
            take("--pair-weight-regularization", cli.pair_weight_regularization) || take_s("--out-csv", cli.out_csv) ||
            take_s("--summary-json", cli.summary_json))
        {
            continue;
        }
        {
            const std::string p = "--sort-every=";
            if (arg.rfind(p, 0) == 0)
            {
                cli.sort_every = std::max(1, std::atoi(arg.c_str() + p.size()));
                continue;
            }
        }
        if (arg == "--enable-diffusion")
        {
            cli.enable_diffusion = true;
            continue;
        }
        if (arg == "--freeze-t")
        {
            cli.enable_diffusion = false;
            continue;
        }
        if (arg == "--zero-buoyancy")
        {
            cli.enable_boussinesq = false;
            continue;
        }
        if (arg == "--freeze-velocity")
        {
            cli.freeze_velocity = true;
            continue;
        }
        if (arg == "--energy-audit")
        {
            cli.energy_audit = true;
            cli.enable_diffusion = true;
            cli.enable_boussinesq = false;
            cli.freeze_velocity = true;
            cli.boundary_width_factor = 0.0;
            continue;
        }
        if (arg == "--dens-sum-normalize")
        {
            cli.dens_sum_normalize = true;
            continue;
        }
        if (arg == "--no-dens-sum-normalize")
        {
            cli.dens_sum_normalize = false;
            continue;
        }
        if (arg == "--background-hydro-p")
        {
            cli.background_hydro_p = true;
            continue;
        }
        if (arg == "--no-background-hydro-p")
        {
            cli.background_hydro_p = false;
            continue;
        }
        if (arg == "--kernel-correction")
        {
            cli.kernel_correction = true;
            continue;
        }
        if (arg == "--corner-rim-y")
        {
            cli.corner_rim_y = true;
            continue;
        }
        if (arg == "--no-kernel-correction")
        {
            cli.kernel_correction = false;
            continue;
        }
        if (arg == "--dirichlet-ghost")
        {
            cli.dirichlet_ghost = true;
            continue;
        }
        if (arg == "--no-dirichlet-ghost")
        {
            cli.dirichlet_ghost = false;
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            std::cout << "Doc77-C / Doc84-P0 / Doc87-R1-A / Doc92 de Vahl Davis cavity\n"
                      << "  --settle-time= --end-time= --freeze-t|--enable-diffusion --zero-buoyancy\n"
                      << "  --freeze-velocity  --L= --dp= --Ra= --beta= --mu-floor=\n"
                      << "  --boundary-width-factor= (default 0) --circ-band-factor= (default 0.1)\n"
                      << "  --nu-probe-layer-factor= (diagnostic grad-Nu only)\n"
                      << "  --energy-audit  (single-hot-wall thermal energy ledger)\n"
                      << "  --contact-pair-factor= (diagnostic scaled Nu only; formal=raw)\n"
                      << "  --dirichlet-ghost / --no-dirichlet-ghost (default on; Adami 2φ_bc−φ_i)\n"
                      << "  --k-wall-factor= (default 1e6; harmonic → 2 k_f)\n"
                      << "  --wall-bw-factor= (dummy wall thickness / dp; default 3)\n"
                      << "  --corner-rim-dp= (Dirichlet side rim near z edges; 0=off)\n"
                      << "  --corner-rim-y  (also rim y front/back)\n"
                      << "  --corner-rim-k-factor= (rim k / k_side; default 0)\n"
                      << "  --pair-weight-regularization= (Laplace ε; default 0.01)\n"
                      << "  --dens-sum-normalize / --no-dens-sum-normalize (default on)\n"
                      << "  --background-hydro-p / --no-background-hydro-p (default on;\n"
                      << "      on ⇒ p' momentum + omit GravityForceCK, Doc105)\n"
                      << "  --kernel-correction / --no-kernel-correction (default on)\n"
                      << "  --sort-every= (particle sort period; default 50)\n"
                      << "  --out-csv= --summary-json=\n";
            return false;
        }
        std::cerr << "[doc77c] unknown arg: " << arg << std::endl;
        return false;
    }
    return true;
}

class FluidBox : public ComplexShape
{
  public:
    FluidBox(const std::string &name, const Vecd &c, const Vecd &h) : ComplexShape(name)
    {
        add<GeometricShapeBox>(Transform(c), h);
    }
};

class ClosedWall : public ComplexShape
{
  public:
    ClosedWall(const std::string &name, const Vecd &c, const Vecd &fluid_h, Real bw) : ComplexShape(name)
    {
        add<GeometricShapeBox>(Transform(c), fluid_h + Vecd(bw, bw, bw));
        subtract<GeometricShapeBox>(Transform(c), fluid_h);
    }
};

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

/** Capture σ0 = DensitySummation into evolving DensSumSigma0 (survives particle sort). */
inline void hostCaptureDensSumSigma0Field(BaseParticles &particles, double &sigma0_mean)
{
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const size_t n = particles.TotalRealParticles();
    sigma0_mean = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sigma0[i] = std::max(rho_sum[i], TinyReal);
        sigma0_mean += static_cast<double>(sigma0[i]);
    }
    sigma0_mean /= static_cast<double>(std::max<size_t>(n, 1));
    syncVariableToDevice<Real>(particles, "DensSumSigma0");
}

inline void hostApplyDensSumNormalizeField(BaseParticles &particles, Real rho0)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    syncVariableToHost<Real>(particles, "DensSumSigma0");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *sigma0 = particles.getVariableDataByName<Real>("DensSumSigma0");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        rho[i] = rho0 * (rho_sum[i] / (sigma0[i] + TinyReal));
    syncVariableToDevice<Real>(particles, "Density");
}

/** Closed-cavity: p_h = ρ0 g (z_top − z) → BackgroundPressure; buoyancy stays in ForcePrior. */
inline void hostApplyCavityBackgroundHydroPressure(BaseParticles &particles, Real rho0, Real g, Real z_bottom,
                                                   Real H)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        p_bg[i] = rho0 * g * (H - zeta);
    }
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
}

inline void hostZeroVecFieldCavity(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Vecd>(particles, name);
    Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Vecd::Zero();
    syncVariableToDevice<Vecd>(particles, name);
}

/** Doc98 P2b: Acoustic momentum uses p' when flag=1 (closed-box + BackgroundPressure). */
inline void hostSetPerturbationMomentumCavity(BaseParticles &particles, Real flag)
{
    syncVariableToHost<Real>(particles, "PerturbationMomentum");
    Real *pm = particles.getVariableDataByName<Real>("PerturbationMomentum");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        pm[i] = flag;
    syncVariableToDevice<Real>(particles, "PerturbationMomentum");
}

/** Doc105: 0=legacy dens, 1=p' dens with independent wall reference. */
inline void hostSetAcousticDensModeCavity(BaseParticles &particles, Real mode)
{
    syncVariableToHost<Real>(particles, "AcousticDensMode");
    Real *m = particles.getVariableDataByName<Real>("AcousticDensMode");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        m[i] = mode;
    syncVariableToDevice<Real>(particles, "AcousticDensMode");
}

inline void hostSetReferenceBodyAccelerationCavity(BaseParticles &particles, const Vecd &a_ref)
{
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    Vecd *a = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        a[i] = a_ref;
    syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
}

inline void hostSetUniformT(BaseParticles &particles, Real t_val, Real t0)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        temp[i] = t_val;
        delta[i] = t_val - t0;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/** Energy injected into the fluid by clamping side Dirichlet bands [J]. */
struct ConstraintHeat
{
    double Q_hot = 0.0;
    double Q_cold = 0.0;
};

inline ConstraintHeat hostApplySideDirichlet(BaseParticles &particles, Real x0, Real L, Real shell, Real th, Real tc,
                                             Real t0, Real rho0, Real cp)
{
    ConstraintHeat out;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real rho_cp = rho0 * cp;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - x0;
        if (x <= shell)
        {
            out.Q_hot += static_cast<double>(rho_cp * (th - temp[i]) * vol[i]);
            temp[i] = th;
            delta[i] = th - t0;
        }
        else if (x >= L - shell)
        {
            out.Q_cold += static_cast<double>(rho_cp * (tc - temp[i]) * vol[i]);
            temp[i] = tc;
            delta[i] = tc - t0;
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    return out;
}

inline void hostInitLinearT(BaseParticles &particles, Real x0, Real L, Real th, Real tc, Real t0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real xi = std::min(Real(1), std::max(Real(0), (pos[i][0] - x0) / (L + TinyReal)));
        const Real Ti = th + (tc - th) * xi;
        temp[i] = Ti;
        delta[i] = Ti - t0;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

/**
 * Contact-only Laplace heat into fluid over dt [J].
 * Raw pair energy ΔE_i = −laplace_contact_i·Vol_i·dt (matches ApplyLaplaceDiffusionDtCK).
 * Left/right attribution by proximity to wall faces (not midplane; Doc87).
 *
 * pair_factor scales diagnostic Nu only; formal Nu / thermal_energy_pass use raw (Doc87).
 */
struct ContactHeatDiag
{
    ConstraintHeat heat; /**< scaled by pair_factor (diagnostic) */
    ConstraintHeat heat_raw;
    double Q_inner_raw = 0.0;
    double Q_contact_raw = 0.0;
    double sum_vol = 0.0;
    double sum_vol_lap_c = 0.0;
    double t_left_mean = 0.0;
    double t_right_mean = 0.0;
    size_t n_left = 0;
    size_t n_right = 0;
    /** Particle with largest |dE_raw| (Doc197 spike localization). */
    double max_abs_dE = 0.0;
    double max_dE = 0.0;
    double max_lap_c = 0.0;
    double max_T = 0.0;
    double max_x = 0.0, max_y = 0.0, max_z = 0.0;
    size_t max_i = 0;
};

inline ContactHeatDiag hostContactWallHeat(BaseParticles &particles, const std::vector<Real> &laplace_inner,
                                           Real x0, Real L, Real dt, Real dp, Real pair_factor,
                                           bool single_hot_wall = false)
{
    ContactHeatDiag out;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalLaplaceTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *lap = particles.getVariableDataByName<Real>(kOphelieThermalLaplaceTField);
    const size_t n = particles.TotalRealParticles();
    const Real band = std::max(Real(2) * dp, Real(0.1) * L);
    double tl = 0.0, tr = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real lap_in = (i < laplace_inner.size() ? laplace_inner[i] : Real(0));
        const Real lap_c = lap[i] - lap_in;
        const double dE_inner = -static_cast<double>(lap_in * vol[i] * dt);
        const double dE_raw = -static_cast<double>(lap_c * vol[i] * dt);
        const double dE = dE_raw * static_cast<double>(pair_factor);
        out.Q_inner_raw += dE_inner;
        out.Q_contact_raw += dE_raw;
        out.sum_vol += static_cast<double>(vol[i]);
        out.sum_vol_lap_c += static_cast<double>(lap_c * vol[i]);
        const double abs_dE = std::abs(dE_raw);
        if (abs_dE > out.max_abs_dE)
        {
            out.max_abs_dE = abs_dE;
            out.max_dE = dE_raw;
            out.max_lap_c = static_cast<double>(lap_c);
            out.max_T = static_cast<double>(temp[i]);
            out.max_x = static_cast<double>(pos[i][0]);
            out.max_y = static_cast<double>(pos[i][1]);
            out.max_z = static_cast<double>(pos[i][2]);
            out.max_i = i;
        }
        const Real x = pos[i][0] - x0;
        const Real dist_left = x;
        const Real dist_right = L - x;
        const bool near_left = dist_left <= band;
        const bool near_right = dist_right <= band;
        if (single_hot_wall || (near_left && !near_right) || (near_left && near_right && dist_left <= dist_right))
        {
            out.heat_raw.Q_hot += dE_raw;
            out.heat.Q_hot += dE;
        }
        else if (near_right)
        {
            out.heat_raw.Q_cold += dE_raw;
            out.heat.Q_cold += dE;
        }
        else if (dist_left <= dist_right)
        {
            out.heat_raw.Q_hot += dE_raw;
            out.heat.Q_hot += dE;
        }
        else
        {
            out.heat_raw.Q_cold += dE_raw;
            out.heat.Q_cold += dE;
        }
        if (x <= dp)
        {
            tl += static_cast<double>(temp[i]);
            ++out.n_left;
        }
        if (x >= L - dp)
        {
            tr += static_cast<double>(temp[i]);
            ++out.n_right;
        }
    }
    out.t_left_mean = out.n_left ? tl / static_cast<double>(out.n_left) : 0.0;
    out.t_right_mean = out.n_right ? tr / static_cast<double>(out.n_right) : 0.0;
    return out;
}

/** Doc87: H = Σ ρ0 cp V ΔT (matches ApplyLaplaceDiffusionDtCK accumulation). */
inline double hostThermalEnthalpyDelta(BaseParticles &particles, Real rho0, Real cp)
{
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double H = 0.0;
    const double rho_cp = static_cast<double>(rho0 * cp);
    for (size_t i = 0; i < n; ++i)
        H += rho_cp * static_cast<double>(vol[i]) * static_cast<double>(delta[i]);
    return H;
}

inline void hostCopyLaplaceField(BaseParticles &particles, std::vector<Real> &out)
{
    syncVariableToHost<Real>(particles, kOphelieThermalLaplaceTField);
    const Real *lap = particles.getVariableDataByName<Real>(kOphelieThermalLaplaceTField);
    const size_t n = particles.TotalRealParticles();
    out.resize(n);
    for (size_t i = 0; i < n; ++i)
        out[i] = lap[i];
}

/**
 * Left/right walls: Dirichlet T. Other faces: k=0 (adiabatic via harmonic mean).
 *
 * Side-wall conductivity uses k_wall >> k_fluid so harmonicMean(k_f, k_w) → 2 k_f
 * unless --k-wall-factor=1. Dummy T stored here is the face BC φ_bc; Contact Laplace
 * may convert to Adami ghost 2φ_bc−φ_i (--dirichlet-ghost).
 *
 * Doc197: optional corner rim — Dirichlet side dummies within corner_rim of an
 * adiabatic edge (z top/bottom, y front/back) get k=0 to suppress |laplace| spikes.
 */
inline void hostPrepareWallThermalBC(BaseParticles &wall_particles, Real x0, Real L, Real th, Real tc, Real t0,
                                     Real k_fluid, Real k_wall_factor = Real(1.0e6), bool single_hot_wall = false,
                                     Real z0 = Real(0), Real y_lo = Real(0), Real y_hi = Real(0),
                                     Real corner_rim = Real(0), Real corner_rim_k_factor = Real(0))
{
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(wall_particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(wall_particles, kOphelieThermalConductivityField);
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    Real *temp = wall_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = wall_particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    Real *cond = wall_particles.getVariableDataByName<Real>(kOphelieThermalConductivityField);
    const Real k_side = k_fluid * k_wall_factor;
    const size_t n = wall_particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0];
        const Real y = pos[i][1];
        const Real z = pos[i][2];
        const bool on_left = x < x0 + TinyReal;
        const bool on_right = x > x0 + L - TinyReal;
        // y_lo/y_hi unused unless caller sets corner_rim and enables y via y_hi>y_lo span check:
        // pass y_hi<=y_lo to disable y rim; otherwise include front/back.
        const bool use_y_rim = (y_hi > y_lo + TinyReal);
        const bool near_adi_edge =
            (corner_rim > TinyReal) &&
            ((z < z0 + corner_rim) || (z > z0 + L - corner_rim) ||
             (use_y_rim && ((y < y_lo + corner_rim) || (y > y_hi - corner_rim))));
        const Real k_rim = k_side * std::max(corner_rim_k_factor, Real(0));
        if (on_left)
        {
            temp[i] = th;
            delta[i] = th - t0;
            cond[i] = near_adi_edge ? k_rim : k_side;
        }
        else if (on_right && !single_hot_wall)
        {
            temp[i] = tc;
            delta[i] = tc - t0;
            cond[i] = near_adi_edge ? k_rim : k_side;
        }
        else
        {
            // Adiabatic faces (and right wall in single-hot energy audit).
            temp[i] = single_hot_wall ? tc : t0;
            delta[i] = (single_hot_wall ? tc : t0) - t0;
            cond[i] = Real(0);
        }
    }
    syncVariableToDevice<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(wall_particles, kOphelieThermalDeltaTField);
    syncVariableToDevice<Real>(wall_particles, kOphelieThermalConductivityField);
}

struct CavityMetrics
{
    Real u_max_star = 0.0;
    Real v_max_star = 0.0;
    Real y_u_max = 0.0;
    Real x_v_max = 0.0;
    /** Secondary |u| peak on vertical midplane (often opposite wall-jet). */
    Real u_max_star_2 = 0.0;
    Real y_u_max_2 = 0.0;
    Real Nu_grad_hot = 0.0;  /**< diagnostic (1-θ)/(x/L) in first non-overlap layer */
    Real Nu_grad_cold = 0.0;
    Real Nu_bulk = 0.0; /**< |dT/dx|_fit * L / ΔT over 0.2–0.8 L (conduction primary) */
    Real U_max = 0.0;
    Real U_rms = 0.0;
    Real Ma = 0.0;
    Real circ_hot_up = 0.0;
    Real circ_cold_dn = 0.0;
    Real t_star = 0.0;
    Real circ_band = 0.0;
    Real grad_layer = 0.0;
    size_t n_grad_hot = 0;
    size_t n_grad_cold = 0;
    size_t n_circ_hot = 0;
    size_t n_circ_cold = 0;
};

/**
 * Mid-y plane metrics. Circ bands and grad-Nu layer are non-overlapping by construction:
 *   hot circ:  (shell, shell+band]
 *   cold circ: [L-shell-band, L-shell)
 *   grad-Nu:   first layer after shell of thickness ≤ band/2 (diagnostic only).
 * DVD midlines (Doc197): |ux|_max on x=L/2 band; |uz|_max on z=L/2 band; thin mid-y slab ~dp.
 */
inline CavityMetrics hostMetricsFull(BaseParticles &particles, Real x0, Real z0, Real L, Real Ly_half, Real y_center,
                                     Real shell, Real circ_band, Real grad_layer, Real th, Real tc, Real alpha,
                                     Real c0, Real t_phys, Real dp = Real(0))
{
    CavityMetrics m;
    const Real dT = th - tc;
    m.t_star = t_phys * alpha / (L * L + TinyReal);
    m.circ_band = circ_band;
    m.grad_layer = grad_layer;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real U_scale = alpha / (L + TinyReal);
    const Real dp_eff = (dp > TinyReal) ? dp : std::max(Real(1.0e-3) * L, Real(1.5) * shell);
    // Thin mid-y (was ≤0.5*Ly_half ≈ half-depth — too thick for 2D-equivalent midline).
    const Real y_mid_tol = std::max(Real(1.5) * dp_eff, Real(1.5) * shell);
    const Real mid_frac = std::max(Real(1.5) * dp_eff / (L + TinyReal), Real(0.02));

    Real u_max = 0.0, v_max = 0.0, U_max = 0.0;
    Real y_um = 0.0, x_vm = 0.0;
    Real u_max_lo = 0.0, u_max_hi = 0.0, y_um_lo = 0.0, y_um_hi = 0.0;
    Real v_max_lo = 0.0, v_max_hi = 0.0, x_vm_lo = 0.0, x_vm_hi = 0.0;
    double Nu_h = 0.0, Nu_c = 0.0, wh = 0.0, wc = 0.0;
    double uz_hot = 0.0, uz_cold = 0.0, mh = 0.0, mc = 0.0;
    double sx = 0.0, sT = 0.0, sx2 = 0.0, sxT = 0.0, sw = 0.0;
    double sum_u2 = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - x0;
        const Real z = pos[i][2] - z0;
        const Real y = pos[i][1];
        const Real U_i = vel[i].norm();
        U_max = std::max(U_max, U_i);
        sum_u2 += static_cast<double>(U_i * U_i);
        const Real xs = x / (L + TinyReal);
        const Real zs = z / (L + TinyReal);
        const Real u_star = vel[i][0] / (U_scale + TinyReal);
        const Real v_star = vel[i][2] / (U_scale + TinyReal);
        const bool mid_y = std::abs(y - y_center) <= y_mid_tol;
        if (mid_y)
        {
            // Vertical midline x≈L/2: DVD U (=ux) extremum vs z.
            // Track lower/upper half separately — single-cell |ux| has dual wall jets.
            if (std::abs(xs - Real(0.5)) < mid_frac)
            {
                const Real au = std::abs(u_star);
                if (au > u_max)
                {
                    u_max = au;
                    y_um = zs;
                }
                if (zs < Real(0.5))
                {
                    if (au > u_max_lo)
                    {
                        u_max_lo = au;
                        y_um_lo = zs;
                    }
                }
                else
                {
                    if (au > u_max_hi)
                    {
                        u_max_hi = au;
                        y_um_hi = zs;
                    }
                }
            }
            // Horizontal midline z≈L/2: DVD V (=uz) extremum vs x (dual jets left/right).
            if (std::abs(zs - Real(0.5)) < mid_frac)
            {
                const Real av = std::abs(v_star);
                if (av > v_max)
                {
                    v_max = av;
                    x_vm = xs;
                }
                if (xs < Real(0.5))
                {
                    if (av > v_max_lo)
                    {
                        v_max_lo = av;
                        x_vm_lo = xs;
                    }
                }
                else
                {
                    if (av > v_max_hi)
                    {
                        v_max_hi = av;
                        x_vm_hi = xs;
                    }
                }
            }
        }
        // Doc84: full wall height for Nu/circ (not mid-y only).
        const Real theta = (temp[i] - tc) / (dT + TinyReal);
        if (x > shell && x <= shell + circ_band)
        {
            uz_hot += vel[i][2] * vol[i];
            mh += vol[i];
            ++m.n_circ_hot;
        }
        if (x < L - shell && x >= L - shell - circ_band)
        {
            uz_cold += vel[i][2] * vol[i];
            mc += vol[i];
            ++m.n_circ_cold;
        }
        // Grad diagnostic: first fluid layer after shell; open on far side to catch lattice.
        if (grad_layer > TinyReal && x > shell && x < shell + grad_layer + TinyReal)
        {
            const Real dxs = std::max(x, Real(0.5) * grad_layer) / (L + TinyReal);
            Nu_h += ((Real(1) - theta) / (dxs + TinyReal)) * vol[i];
            wh += vol[i];
            ++m.n_grad_hot;
        }
        if (grad_layer > TinyReal && x < L - shell && x > L - shell - grad_layer - TinyReal)
        {
            const Real dxs = std::max(L - x, Real(0.5) * grad_layer) / (L + TinyReal);
            Nu_c += (theta / (dxs + TinyReal)) * vol[i];
            wc += vol[i];
            ++m.n_grad_cold;
        }
        // Bulk slope fit (exclude near-wall SPH layer).
        if (xs > Real(0.2) && xs < Real(0.8))
        {
            const double w = static_cast<double>(vol[i]);
            const double xx = static_cast<double>(x);
            const double TT = static_cast<double>(temp[i]);
            sw += w;
            sx += w * xx;
            sT += w * TT;
            sx2 += w * xx * xx;
            sxT += w * xx * TT;
        }
    }
    if (sw > TinyReal)
    {
        const double denom = sw * sx2 - sx * sx;
        const double slope = (std::abs(denom) > TinyReal) ? (sw * sxT - sx * sT) / denom : 0.0;
        m.Nu_bulk = static_cast<Real>(std::abs(slope) * static_cast<double>(L) / (static_cast<double>(dT) + TinyReal));
    }
    // Prefer DVD-style upper-half |U| max when dual jets are comparable (ref Y≈0.813).
    if (u_max_hi > TinyReal && u_max_lo > TinyReal && u_max_hi >= Real(0.9) * u_max_lo)
    {
        m.u_max_star = u_max_hi;
        m.y_u_max = y_um_hi;
        m.u_max_star_2 = u_max_lo;
        m.y_u_max_2 = y_um_lo;
    }
    else
    {
        m.u_max_star = u_max;
        m.y_u_max = y_um;
        m.u_max_star_2 = (y_um >= Real(0.5)) ? u_max_lo : u_max_hi;
        m.y_u_max_2 = (y_um >= Real(0.5)) ? y_um_lo : y_um_hi;
    }
    // Prefer DVD-style left-half |V| max when dual jets comparable (ref X≈0.178).
    if (v_max_lo > TinyReal && v_max_hi > TinyReal && v_max_lo >= Real(0.9) * v_max_hi)
    {
        m.v_max_star = v_max_lo;
        m.x_v_max = x_vm_lo;
    }
    else
    {
        m.v_max_star = v_max;
        m.x_v_max = x_vm;
    }
    m.U_max = U_max;
    m.U_rms = (n > 0) ? static_cast<Real>(std::sqrt(sum_u2 / static_cast<double>(n))) : Real(0);
    m.Ma = U_max / (c0 + TinyReal);
    m.Nu_grad_hot = wh > TinyReal ? static_cast<Real>(Nu_h / wh) : Real(0);
    m.Nu_grad_cold = wc > TinyReal ? static_cast<Real>(Nu_c / wc) : Real(0);
    m.circ_hot_up = mh > TinyReal ? static_cast<Real>(uz_hot / mh) : Real(0);
    m.circ_cold_dn = mc > TinyReal ? static_cast<Real>(uz_cold / mc) : Real(0);
    return m;
}

/**
 * Doc197 scheme-1 diagnostic: wall Nu from SPH ∇T (optional B-matrix) including wall dummies
 * with Adami T. Not Contact-Laplace flux; parallel to formal Contact Nu (no ×2 / no β).
 * Nu = |-(L/ΔT) ∇T·n_into_fluid| averaged in the first fluid layer.
 */
struct CorrectedWallNu
{
    Real Nu_hot = 0.0;
    Real Nu_cold = 0.0;
    size_t n_hot = 0;
    size_t n_cold = 0;
};

inline Real cubicSplineDw3d(Real r, Real h)
{
    // Match SPHXsys KernelCubicBSpline 3D: W = (1/(π h^3)) W1D(q), dW/dr = (1/(π h^4)) dW1D/dq.
    if (!(h > TinyReal) || !(r > TinyReal) || r >= Real(2) * h)
        return Real(0);
    const Real q = r / h;
    const Real inv_h = Real(1) / h;
    const Real factor_W_3D = inv_h * inv_h * inv_h / Pi;
    Real dWdq = Real(0);
    if (q < Real(1))
        dWdq = Real(9) * q * q / Real(4) - Real(3) * q;
    else
        dWdq = -Real(3) * (Real(2) - q) * (Real(2) - q) / Real(4);
    return factor_W_3D * inv_h * dWdq;
}

inline CorrectedWallNu hostCorrectedWallNuFromGradT(BaseParticles &fluid_particles, BaseParticles &wall_particles,
                                                    Real x0, Real L, Real th, Real tc, Real layer, Real h,
                                                    bool use_B, bool dirichlet_ghost)
{
    CorrectedWallNu out;
    syncVariableToHost<Vecd>(fluid_particles, "Position");
    syncVariableToHost<Real>(fluid_particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(fluid_particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(wall_particles, "VolumetricMeasure");
    const Vecd *xpos = fluid_particles.getVariableDataByName<Vecd>("Position");
    const Real *T = fluid_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *Vol = fluid_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *wpos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Real *Tw = wall_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *wVol = wall_particles.getVariableDataByName<Real>("VolumetricMeasure");
    Matd *B = nullptr;
    if (use_B)
    {
        syncVariableToHost<Matd>(fluid_particles, "LinearCorrectionMatrix");
        B = fluid_particles.getVariableDataByName<Matd>("LinearCorrectionMatrix");
    }
    const size_t nf = fluid_particles.TotalRealParticles();
    const size_t nw = wall_particles.TotalRealParticles();
    const Real dT = th - tc;
    const Real cut = Real(2) * h;
    const Real cut2 = cut * cut;
    double sum_h = 0.0, sum_c = 0.0, wh = 0.0, wc = 0.0;

    for (size_t i = 0; i < nf; ++i)
    {
        const Real x = xpos[i][0] - x0;
        const bool near_hot = (x > TinyReal && x <= layer);
        const bool near_cold = (x < L - TinyReal && x >= L - layer);
        if (!near_hot && !near_cold)
            continue;

        Vecd grad = Vecd::Zero();
        const Real Ti = T[i];
        for (size_t j = 0; j < nf; ++j)
        {
            if (j == i)
                continue;
            const Vecd rij = xpos[i] - xpos[j];
            const Real r2 = rij.squaredNorm();
            if (r2 <= TinyReal || r2 > cut2)
                continue;
            const Real r = std::sqrt(r2);
            const Real dWdr = cubicSplineDw3d(r, h);
            // SPH ∇_i W_ij ~ (dW/dr) * (r_i - r_j)/r ; grad T ~ Σ (T_j-T_i) ∇_i W V_j
            grad += (T[j] - Ti) * (dWdr / r) * rij * Vol[j];
        }
        for (size_t j = 0; j < nw; ++j)
        {
            const Vecd rij = xpos[i] - wpos[j];
            const Real r2 = rij.squaredNorm();
            if (r2 <= TinyReal || r2 > cut2)
                continue;
            const Real r = std::sqrt(r2);
            const Real dWdr = cubicSplineDw3d(r, h);
            const Real Tbc = Tw[j];
            const Real Tj = dirichlet_ghost ? (Real(2) * Tbc - Ti) : Tbc;
            grad += (Tj - Ti) * (dWdr / r) * rij * wVol[j];
        }
        if (B != nullptr)
            grad = B[i] * grad;

        // n_into_fluid: hot +x, cold -x. Signed heat into fluid / q_ref; flip cold to match Contact |Q|.
        const Vecd n_into = near_hot ? Vecd(1.0, 0.0, 0.0) : Vecd(-1.0, 0.0, 0.0);
        const Real Nu_signed = -(L / (dT + TinyReal)) * grad.dot(n_into);
        const Real Nu_i = near_hot ? Nu_signed : -Nu_signed;
        if (near_hot)
        {
            sum_h += static_cast<double>(Nu_i) * static_cast<double>(Vol[i]);
            wh += static_cast<double>(Vol[i]);
            ++out.n_hot;
        }
        else
        {
            sum_c += static_cast<double>(Nu_i) * static_cast<double>(Vol[i]);
            wc += static_cast<double>(Vol[i]);
            ++out.n_cold;
        }
    }
    out.Nu_hot = wh > TinyReal ? static_cast<Real>(sum_h / wh) : Real(0);
    out.Nu_cold = wc > TinyReal ? static_cast<Real>(sum_c / wc) : Real(0);
    return out;
}

/** Export DVD midlines: vertical (x≈L/2) ux(z) and horizontal (z≈L/2) uz(x), mid-y slab. */
inline void hostExportMidlineProfiles(BaseParticles &particles, const std::string &path, Real x0, Real z0, Real L,
                                      Real y_center, Real alpha, Real dp)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const size_t n = particles.TotalRealParticles();
    const Real U_scale = alpha / (L + TinyReal);
    const Real y_tol = std::max(Real(1.5) * dp, Real(1.0e-9));
    const Real mid_tol = std::max(Real(1.5) * dp, Real(0.02) * L);
    std::ofstream ofs(path);
    ofs << "line,s_norm,u_star,v_star,T,x,y,z\n";
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - x0;
        const Real z = pos[i][2] - z0;
        const Real y = pos[i][1];
        if (std::abs(y - y_center) > y_tol)
            continue;
        const Real xs = x / (L + TinyReal);
        const Real zs = z / (L + TinyReal);
        const Real u_star = vel[i][0] / (U_scale + TinyReal);
        const Real v_star = vel[i][2] / (U_scale + TinyReal);
        if (std::abs(x - Real(0.5) * L) <= mid_tol)
            ofs << "vertical_xmid," << zs << "," << u_star << "," << v_star << "," << temp[i] << ","
                << pos[i][0] << "," << y << "," << pos[i][2] << "\n";
        if (std::abs(z - Real(0.5) * L) <= mid_tol)
            ofs << "horizontal_zmid," << xs << "," << u_star << "," << v_star << "," << temp[i] << ","
                << pos[i][0] << "," << y << "," << pos[i][2] << "\n";
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

    const Real L = cli.L;
    const Real Ly = cli.depth_factor * L;
    const Real dp = cli.dp;
    const Real BW = std::max(cli.wall_bw_factor, Real(1)) * dp;
    const Real shell = cli.boundary_width_factor * dp;
    const Real th = cli.t_c + cli.dT;
    const Real t0 = Real(0.5) * (th + cli.t_c);
    const Real beta_mat = std::max(cli.beta, Real(1.0e-3));
    const Real beta_force = cli.enable_boussinesq ? cli.beta : Real(0);
    const Real beta_dT = beta_force * cli.dT;

    const Real nu_alpha = cli.gravity_g * beta_mat * cli.dT * L * L * L / (cli.Ra + TinyReal);
    const Real alpha = std::sqrt(nu_alpha / (cli.Pr + TinyReal));
    const Real nu_ra = cli.Pr * alpha;
    const Real mu_ra = cli.rho0 * nu_ra;
    const Real mu = std::max(mu_ra, cli.mu_floor);
    const Real k = alpha * cli.rho0 * cli.cp;
    const Real U_char = std::max(Real(4) * alpha / (L + TinyReal), Real(0.05));
    const Real c0 = std::max(cli.c0_factor * U_char, Real(10.0));
    const Real q_scale = (k * cli.dT) / (L + TinyReal);

    // Non-overlapping bands: leave a clear mid-cavity gap.
    const Real band_cap = std::max(Real(0), Real(0.45) * L - shell);
    const Real circ_band = std::min(std::max(cli.circ_band_factor * L, dp), band_cap);
    // Grad layer = O(dp), capped below circ_band so Nu diagnostic ≠ circ band overlap.
    const Real grad_layer =
        std::min(std::max(cli.nu_probe_layer_factor * dp, dp), std::max(Real(0.5) * circ_band, dp));

    const Vecd fluid_half(Real(0.5) * L, Real(0.5) * Ly, Real(0.5) * L);
    const Vecd fluid_center(Real(0.5) * L, Real(0.5) * Ly, Real(0.5) * L);
    const Real x0 = 0.0;
    const Real z0 = 0.0;
    const Real y_lo = fluid_center[1] - Real(0.5) * Ly;
    const Real y_hi = fluid_center[1] + Real(0.5) * Ly;
    const Real corner_rim = std::max(cli.corner_rim_dp, Real(0)) * dp;
    const Real y_lo_bc = cli.corner_rim_y ? y_lo : Real(1);
    const Real y_hi_bc = cli.corner_rim_y ? y_hi : Real(0); // y_hi<=y_lo disables y rim
    const Real A_wall_geom = L * Ly;
    const Real L_eff = std::max(L - Real(2) * corner_rim, Real(0.5) * L);
    const Real Ly_eff =
        cli.corner_rim_y ? std::max(Ly - Real(2) * corner_rim, Real(0.5) * Ly) : Ly;
    // Formal Nu always uses geometric side area (DVD wall); A_eff is diagnostic only.
    const Real A_wall = A_wall_geom;
    const Real A_wall_eff = (corner_rim > TinyReal) ? (L_eff * Ly_eff) : A_wall_geom;

    const BoundingBoxd bounds(fluid_center - fluid_half - Vecd(BW, BW, BW),
                              fluid_center + fluid_half + Vecd(BW, BW, BW));
    SPHSystem sph_system(bounds, dp);

    FluidBody fluid(sph_system, makeShared<FluidBox>("Fluid", fluid_center, fluid_half));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(cli.rho0, c0);
    fluid.addMaterialProperty<Viscosity>(mu);
    fluid.generateParticles<BaseParticles, Lattice>();

    SolidBody wall(sph_system, makeShared<ClosedWall>("Wall", fluid_center, fluid_half, BW));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall(fluid, {&wall});
    BaseParticles &particles = fluid.getBaseParticles();
    BaseParticles &wall_particles = wall.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(particles, t0);
    registerOphelieThermalDiffusionAuxFields(particles, k);
    registerOphelieJouleHeatTemperatureField(wall_particles, t0);
    registerOphelieThermalDiffusionAuxFields(wall_particles, k);
    particles.registerStateVariable<Real>("DensSumSigma0");
    particles.addEvolvingVariable<Real>("DensSumSigma0");
    hostSetUniformT(particles, t0, t0);
    auto prepare_wall_bc = [&](bool single_hot = false) {
        hostPrepareWallThermalBC(wall_particles, x0, L, th, cli.t_c, t0, k, cli.k_wall_factor, single_hot, z0,
                                 y_lo_bc, y_hi_bc, corner_rim, cli.corner_rim_k_factor);
    };
    prepare_wall_bc();
    const bool use_fluid_dirichlet = shell > Real(0.5) * dp;

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cells(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cells(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> fluid_rel(fluid_inner, fluid_wall);
    ParticleSortCK<MainExecutionPolicy> particle_sort(fluid);

    Gravity gravity(Vecd(0.0, 0.0, -cli.gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(fluid, gravity);
    StateDynamics<MainExecutionPolicy, BoussinesqBuoyancyForceCK> boussinesq(fluid, gravity, beta_force, t0,
                                                                             kOphelieTemperatureField);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_pos(fluid);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(fluid_inner, 0.5), fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK> acoustic_1st(
        fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalf<
                                                   Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd_wall(fluid_wall);
    acoustic_2nd.addPostContactInteraction(acoustic_2nd_wall);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        acoustic_1st_corr(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        acoustic_2nd_corr(fluid_inner, fluid_wall);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>, Contact<>>> density_sum(
        fluid_inner, fluid_wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, Internal>> density_reg(fluid);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK> viscous(fluid_inner, fluid_wall);

    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace_inner(
        fluid_inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
        cli.pair_weight_regularization);
    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Contact<>>> laplace_wall(
        fluid_wall, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
        cli.pair_weight_regularization, cli.dirichlet_ghost);
    StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK> apply_diff(fluid, cli.rho0, cli.cp, t0);
    UpdateRelation<MainExecutionPolicy, Inner<>> fluid_inner_only(fluid_inner);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> adv_dt(fluid, U_char);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> ac_dt(fluid);
    ReduceDynamicsCK<MainExecutionPolicy, OphelieVelocityMaxReduceCK> u_max_reduce(fluid);

    SingleVariable<Real> *sv_t = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    const bool thermal_only = cli.freeze_velocity && cli.enable_diffusion;

    wall_normal.exec();
    fluid_cells.exec();
    wall_cells.exec();
    fluid_rel.exec();
    // Doc105 closed-box: BackgroundPressure absorbs ρ0 g → omit GravityForceCK; dens/mom on p'.
    if (cli.background_hydro_p && !thermal_only)
    {
        hostSetPerturbationMomentumCavity(particles, Real(1));
        hostSetAcousticDensModeCavity(particles, Real(1));
        hostSetReferenceBodyAccelerationCavity(particles, Vecd::Zero());
        hostZeroVecFieldCavity(particles, "GravityForceCK");
        hostZeroVecFieldCavity(particles, "PreviousGravityForceCK");
    }
    else if (!thermal_only)
    {
        constant_gravity.exec();
    }

    std::cout << "[doc77c] start n=" << particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " L=" << L << " Ly=" << Ly
              << " dp=" << dp << " shell=" << shell << " circ_band=" << circ_band << " grad_layer=" << grad_layer
              << " A_wall=" << A_wall << " Ra=" << cli.Ra << " Pr=" << cli.Pr << " beta_force=" << beta_force
              << " beta_dT=" << beta_dT << " alpha=" << alpha << " nu_ra=" << nu_ra << " mu=" << mu
              << " k=" << k << " c0=" << c0 << " U_char=" << U_char << " settle=" << cli.settle_time
              << " end_time=" << cli.end_time << " boussinesq=" << (cli.enable_boussinesq ? 1 : 0)
              << " diffusion=" << (cli.enable_diffusion ? 1 : 0)
              << " fluid_dirichlet=" << (use_fluid_dirichlet ? 1 : 0)
              << " wall_thermal_bc=1"
              << " dirichlet_ghost=" << (cli.dirichlet_ghost ? 1 : 0)
              << " k_wall_factor=" << cli.k_wall_factor
              << " pair_weight_reg=" << cli.pair_weight_regularization
              << " dens_sum_normalize=" << (cli.dens_sum_normalize ? 1 : 0)
              << " background_hydro_p=" << (cli.background_hydro_p ? 1 : 0)
              << " kernel_correction=" << (cli.kernel_correction ? 1 : 0)
              << " label=3D_prism_doc92_wall_contact_laplace_hydro_defaults" << std::endl
              << std::flush;

    const Real Nu_ref = Real(1.118);
    const Real u_ref = Real(3.649);
    const Real v_ref = Real(3.697);

    std::ofstream ofs(cli.out_csv);
    ofs << "phase,t,t_star,U_max,Ma,u_max_star,v_max_star,y_umax,x_vmax,"
           "Nu_wall_hot,Nu_wall_cold,Nu_grad_hot,Nu_grad_cold,Nu_gradcorr_hot,Nu_gradcorr_cold,"
           "circ_hot_uz,circ_cold_uz,dt_fluid,dt_thermal,n_th_sub,n_th_capped\n";

    double acc_hot = 0.0, acc_cold = 0.0, acc_w = 0.0;
    double late_Nu_grad_h = 0.0, late_Nu_grad_c = 0.0, late_Nu_bulk = 0.0, late_w = 0.0;
    double late_Nu_gradcorr_h = 0.0, late_Nu_gradcorr_c = 0.0;
    const Real h_sph = fluid.getSPHAdaptation().ReferenceSmoothingLength();
    double Q_hot_late = 0.0, Q_cold_late = 0.0, t_thermal_late = 0.0;
    double Q_hot_early = 0.0, Q_cold_early = 0.0, t_thermal_early = 0.0;
    double Q_hot_obs = 0.0, Q_cold_obs = 0.0, t_thermal_obs = 0.0;
    double t_fluid_total = 0.0, t_thermal_total = 0.0;
    size_t n_th_sub_total = 0, n_th_capped_total = 0, n_adv_steps = 0;
    size_t acc_n = 0;
    Real t_late_gate = Real(1.0e300);

    auto metrics_now = [&]() -> CavityMetrics {
        return hostMetricsFull(particles, x0, z0, L, Real(0.5) * Ly, fluid_center[1], shell, circ_band, grad_layer, th,
                               cli.t_c, alpha, c0, sv_t->getValue(), dp);
    };

    std::vector<Real> laplace_inner_buf;
    double dens_sigma0_mean = 0.0;
    bool dens_sigma0_ready = false;

    // Need B also for freeze-velocity conduction (gradcorr Nu diagnostic).
    if (cli.kernel_correction)
    {
        fluid_linear_correction_matrix.exec();
    }
    if (cli.background_hydro_p && !thermal_only)
    {
        hostApplyCavityBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z0, L);
    }

    auto advance_to = [&](Real t_target, const char *phase, bool use_boussinesq, bool use_side_t,
                          bool accumulate_obs) {
        size_t local_iter = 0;
        Real next_report = sv_t->getValue();
        while (sv_t->getValue() < t_target)
        {
            Real dt_fluid = 0.0;
            Real dt_adv = 0.0;
            if (!thermal_only)
            {
                density_sum.exec();
                if (cli.dens_sum_normalize && !dens_sigma0_ready)
                {
                    hostCaptureDensSumSigma0Field(particles, dens_sigma0_mean);
                    dens_sigma0_ready = true;
                    std::cout << "[doc77c] dens_sum_normalize sigma0_mean=" << dens_sigma0_mean
                              << " packing=" << (dens_sigma0_mean / (cli.rho0 + TinyReal)) << std::endl
                              << std::flush;
                }
                density_reg.exec();
                if (cli.dens_sum_normalize)
                    hostApplyDensSumNormalizeField(particles, cli.rho0);
                if (cli.background_hydro_p)
                    hostApplyCavityBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z0, L);
                adv_setup.exec();
                if (cli.kernel_correction)
                    fluid_linear_correction_matrix.exec();
                dt_adv = adv_dt.exec();
                viscous.exec();
                if (!cli.background_hydro_p)
                    constant_gravity.exec();
                if (use_boussinesq && beta_force > TinyReal)
                {
                    boussinesq.exec();
                }
                Real relax = 0.0;
                Real dt_ac = 0.0;
                size_t ac_sub = 0;
                const Real t_remain_phase = t_target - sv_t->getValue();
                while (relax < dt_adv && ac_sub < 2000)
                {
                    dt_ac = ac_dt.exec();
                    const Real dt_left = std::min(dt_adv - relax, t_remain_phase - relax);
                    if (!(dt_left > TinyReal))
                        break;
                    if (dt_ac > dt_left)
                        dt_ac = dt_left;
                    if (!(dt_ac > TinyReal) || !std::isfinite(static_cast<double>(dt_ac)))
                        break;
                    if (cli.kernel_correction)
                        acoustic_1st_corr.exec(dt_ac);
                    else
                        acoustic_1st.exec(dt_ac);
                    if (cli.kernel_correction)
                        acoustic_2nd_corr.exec(dt_ac);
                    else
                        acoustic_2nd.exec(dt_ac);
                    if (use_boussinesq && beta_force > TinyReal)
                    {
                        boussinesq.exec();
                    }
                    relax += dt_ac;
                    sv_t->incrementValue(dt_ac);
                    ++ac_sub;
                }
                dt_fluid = relax;
            }
            else
            {
                // Pure conduction: advance time with thermal-stable chunks (no density/acoustic).
                const Real dt_th_max = Real(0.1) * dp * dp / (alpha + TinyReal);
                dt_fluid = Real(8) * dt_th_max;
                dt_adv = dt_fluid;
                sv_t->incrementValue(dt_fluid);
            }
            t_fluid_total += static_cast<double>(dt_fluid);
            ++n_adv_steps;

            // Thermal subcycle fills the same physical time (Doc84).
            Real dt_thermal_step = 0.0;
            size_t n_th = 0, n_th_cap = 0;
            double Qh = 0.0, Qc = 0.0;
            if (cli.enable_diffusion && use_side_t && dt_fluid > TinyReal)
            {
                const Real dt_th_max = Real(0.1) * dp * dp / (alpha + TinyReal);
                Real remain = dt_fluid;
                while (remain > TinyReal && n_th < 10000)
                {
                    Real dt_th = std::min(remain, dt_th_max);
                    if (dt_th + TinyReal < remain)
                    {
                        ++n_th_cap;
                    }
                    fluid_cells.exec();
                    wall_cells.exec();
                    fluid_rel.exec();
                    fluid_inner_only.exec();
                    prepare_wall_bc();
                    laplace_inner.exec();
                    if (!use_fluid_dirichlet)
                    {
                        hostCopyLaplaceField(particles, laplace_inner_buf);
                        laplace_wall.exec(); // wall Contact only when fluid shell Dirichlet is off
                    }
                    apply_diff.exec(dt_th);
                    if (use_fluid_dirichlet)
                    {
                        const ConstraintHeat ch =
                            hostApplySideDirichlet(particles, x0, L, shell, th, cli.t_c, t0, cli.rho0, cli.cp);
                        Qh += ch.Q_hot;
                        Qc += ch.Q_cold;
                    }
                    else
                    {
                        const ContactHeatDiag ch =
                            hostContactWallHeat(particles, laplace_inner_buf, x0, L, dt_th, dp,
                                               cli.contact_pair_factor, cli.energy_audit);
                        // Doc87: accumulate raw Contact heat for formal Nu / energy bookkeeping.
                        Qh += ch.heat_raw.Q_hot;
                        Qc += ch.heat_raw.Q_cold;
                        // Periodic sample + always log when |Qh| spikes (Doc197 localization).
                        const bool periodic_ch = (n_th == 0 && local_iter % 40 == 0);
                        // Baseline |Qh|~0.016; mild~0.06; mega-spikes ≳0.3 (corner particles).
                        const bool spike_ch =
                            (n_th == 0) && (std::abs(ch.heat_raw.Q_hot) > Real(0.3) ||
                                            std::abs(ch.heat_raw.Q_cold) > Real(0.3));
                        if (periodic_ch || spike_ch)
                        {
                            std::cout << "[doc77c][contact-heat] pair_factor_diag=" << cli.contact_pair_factor
                                      << " sum_vol=" << ch.sum_vol << " Qh_raw=" << ch.heat_raw.Q_hot
                                      << " Qh_scaled=" << ch.heat.Q_hot << " Qc_raw=" << ch.heat_raw.Q_cold
                                      << " T_left=" << ch.t_left_mean << " T_right=" << ch.t_right_mean
                                      << " nL=" << ch.n_left << " nR=" << ch.n_right
                                      << " spike=" << (spike_ch ? 1 : 0) << " max_abs_dE=" << ch.max_abs_dE
                                      << " max_dE=" << ch.max_dE << " max_lap_c=" << ch.max_lap_c
                                      << " max_T=" << ch.max_T << " max_xyz=" << ch.max_x << "," << ch.max_y << ","
                                      << ch.max_z << " max_i=" << ch.max_i << std::endl
                                      << std::flush;
                        }
                    }
                    remain -= dt_th;
                    dt_thermal_step += dt_th;
                    ++n_th;
                }
            }
            else if (use_side_t && use_fluid_dirichlet)
            {
                const ConstraintHeat ch =
                    hostApplySideDirichlet(particles, x0, L, shell, th, cli.t_c, t0, cli.rho0, cli.cp);
                Qh += ch.Q_hot;
                Qc += ch.Q_cold;
            }

            t_thermal_total += static_cast<double>(dt_thermal_step);
            n_th_sub_total += n_th;
            n_th_capped_total += n_th_cap;

            if (accumulate_obs && cli.enable_diffusion)
            {
                Q_hot_obs += Qh;
                Q_cold_obs += Qc;
                t_thermal_obs += static_cast<double>(dt_thermal_step);
                if (sv_t->getValue() >= t_late_gate - TinyReal)
                {
                    Q_hot_late += Qh;
                    Q_cold_late += Qc;
                    t_thermal_late += static_cast<double>(dt_thermal_step);
                }
                else
                {
                    Q_hot_early += Qh;
                    Q_cold_early += Qc;
                    t_thermal_early += static_cast<double>(dt_thermal_step);
                }
            }

            if (cli.freeze_velocity)
            {
                hostZeroVelocity(particles);
            }
            else
            {
                update_pos.exec();
            }
            ++local_iter;
            if (local_iter % 20 == 0)
            {
                const Real um = thermal_only ? Real(0) : u_max_reduce.exec();
                std::cout << "[doc77c][" << phase << "] N=" << local_iter << " t=" << sv_t->getValue()
                          << " U_max=" << um << " dt_adv=" << dt_adv << " dt_fluid=" << dt_fluid
                          << " dt_th=" << dt_thermal_step << " n_th=" << n_th << " n_th_cap=" << n_th_cap
                          << " thermal_only=" << (thermal_only ? 1 : 0)
                          << " kc=" << (cli.kernel_correction ? 1 : 0)
                          << " bg_p=" << (cli.background_hydro_p ? 1 : 0) << std::endl
                          << std::flush;
            }
            if (!thermal_only && (local_iter % cli.sort_every == 0))
            {
                particle_sort.exec();
            }
            if (!thermal_only)
            {
                fluid_cells.exec();
                wall_cells.exec();
                fluid_rel.exec();
            }

            if (sv_t->getValue() >= next_report - TinyReal)
            {
                const CavityMetrics mm = metrics_now();
                const CorrectedWallNu ngc = hostCorrectedWallNuFromGradT(
                    particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, cli.kernel_correction,
                    cli.dirichlet_ghost);
                const Real Nu_wh =
                    (dt_thermal_step > TinyReal)
                        ? static_cast<Real>((Qh / static_cast<double>(dt_thermal_step)) / (A_wall + TinyReal) /
                                            (q_scale + TinyReal))
                        : Real(0);
                const Real Nu_wc =
                    (dt_thermal_step > TinyReal)
                        ? static_cast<Real>((-Qc / static_cast<double>(dt_thermal_step)) / (A_wall + TinyReal) /
                                            (q_scale + TinyReal))
                        : Real(0);
                ofs << phase << "," << sv_t->getValue() << "," << mm.t_star << "," << mm.U_max << "," << mm.Ma << ","
                    << mm.u_max_star << "," << mm.v_max_star << "," << mm.y_u_max << "," << mm.x_v_max << "," << Nu_wh
                    << "," << Nu_wc << "," << mm.Nu_grad_hot << "," << mm.Nu_grad_cold << "," << ngc.Nu_hot << ","
                    << ngc.Nu_cold << "," << mm.circ_hot_up << "," << mm.circ_cold_dn << "," << dt_fluid << ","
                    << dt_thermal_step << "," << n_th << "," << n_th_cap << "\n";
                ofs.flush();
                std::cout << "[doc77c][" << phase << "] t=" << sv_t->getValue() << " circ_hot_uz=" << mm.circ_hot_up
                          << " circ_cold_uz=" << mm.circ_cold_dn << " Umax=" << mm.U_max << " Nu_wall_h=" << Nu_wh
                          << " Nu_grad_h=" << mm.Nu_grad_hot << " Nu_gradcorr_h=" << ngc.Nu_hot
                          << " Nu_bulk=" << mm.Nu_bulk << std::endl
                          << std::flush;
                if (accumulate_obs)
                {
                    acc_hot += mm.circ_hot_up;
                    acc_cold += mm.circ_cold_dn;
                    acc_w += 1.0;
                    ++acc_n;
                    if (sv_t->getValue() >= t_late_gate - TinyReal)
                    {
                        late_Nu_grad_h += mm.Nu_grad_hot;
                        late_Nu_grad_c += mm.Nu_grad_cold;
                        late_Nu_bulk += mm.Nu_bulk;
                        late_Nu_gradcorr_h += ngc.Nu_hot;
                        late_Nu_gradcorr_c += ngc.Nu_cold;
                        late_w += 1.0;
                    }
                }
                next_report += cli.report_every;
            }
        }
    };

    // -------------------------------------------------------------------------
    // Doc87 R1-A: single-hot-wall thermal energy ledger (raw Contact ↔ ΔH).
    // -------------------------------------------------------------------------
    if (cli.energy_audit)
    {
        hostSetUniformT(particles, cli.t_c, t0); // cold uniform init → large resolvable ΔH
        hostZeroVelocity(particles);
        prepare_wall_bc(true);
        fluid_cells.exec();
        wall_cells.exec();
        fluid_rel.exec();
        fluid_inner_only.exec();

        const std::string energy_csv =
            cli.out_csv.find("output/") == 0
                ? (cli.out_csv.substr(0, cli.out_csv.find_last_of('.')) + "_energy.csv")
                : "output/doc77c_energy_audit.csv";
        std::ofstream eofs(energy_csv);
        eofs << "step,time,dt,H_before,H_after,dH,Q_inner_raw,Q_contact_raw,Q_contact_scaled,"
                "energy_residual,cum_abs_Q,cum_abs_resid,n_th\n";

        const Real dt_th_max = Real(0.1) * dp * dp / (alpha + TinyReal);
        const Real t_end = cli.end_time > TinyReal ? cli.end_time : Real(2.0);
        double cum_abs_Q = 0.0, cum_abs_resid = 0.0, cum_dH = 0.0, cum_Q = 0.0;
        size_t step = 0;
        size_t n_bad_dt = 0;
        Real t_acc = 0.0;
        std::vector<Real> lap_buf;
        const size_t max_steps = 200000;
        while (t_acc < t_end - TinyReal && step < max_steps)
        {
            Real remain = std::min(Real(8) * dt_th_max, t_end - t_acc);
            size_t n_th = 0;
            while (remain > TinyReal && n_th < 10000)
            {
                Real dt_th = std::min(remain, dt_th_max);
                if (!(dt_th > TinyReal) || !std::isfinite(static_cast<double>(dt_th)))
                {
                    ++n_bad_dt;
                    break;
                }
                const double H_before = hostThermalEnthalpyDelta(particles, cli.rho0, cli.cp);
                fluid_cells.exec();
                wall_cells.exec();
                fluid_rel.exec();
                fluid_inner_only.exec();
                prepare_wall_bc(true);
                laplace_inner.exec();
                hostCopyLaplaceField(particles, lap_buf);
                laplace_wall.exec();
                apply_diff.exec(dt_th);
                const ContactHeatDiag ch =
                    hostContactWallHeat(particles, lap_buf, x0, L, dt_th, dp, cli.contact_pair_factor, true);
                const double H_after = hostThermalEnthalpyDelta(particles, cli.rho0, cli.cp);
                const double dH = H_after - H_before;
                const double Q_sum = ch.Q_inner_raw + ch.Q_contact_raw;
                const double resid = dH - Q_sum;
                cum_abs_Q += std::abs(ch.Q_inner_raw) + std::abs(ch.Q_contact_raw);
                cum_abs_resid += std::abs(resid);
                cum_dH += dH;
                cum_Q += Q_sum;
                eofs << step << "," << t_acc << "," << dt_th << "," << H_before << "," << H_after << "," << dH << ","
                     << ch.Q_inner_raw << "," << ch.Q_contact_raw << ","
                     << (ch.Q_contact_raw * static_cast<double>(cli.contact_pair_factor)) << "," << resid << ","
                     << cum_abs_Q << "," << cum_abs_resid << "," << n_th << "\n";
                remain -= dt_th;
                t_acc += dt_th;
                sv_t->incrementValue(dt_th);
                ++n_th;
                ++step;
                if (step % 200 == 0)
                {
                    std::cout << "[doc77c][energy-audit] step=" << step << " t=" << t_acc << " dH=" << dH
                              << " Q_in=" << ch.Q_inner_raw << " Q_c=" << ch.Q_contact_raw << " resid=" << resid
                              << " cum_rel=" << (cum_abs_resid / (cum_abs_Q + 1e-30)) << std::endl
                              << std::flush;
                }
            }
            if (n_bad_dt > 0)
                break;
            hostZeroVelocity(particles);
        }
        eofs.flush();
        const double scale = std::max({cum_abs_Q, std::abs(cum_dH), 1e-30});
        const double resid_rel = cum_abs_resid / scale;
        const double abs_tol = 1e-8 * std::max(std::abs(cum_dH), 1.0);
        const bool energy_ok = (n_bad_dt == 0) && (resid_rel <= 0.001 || cum_abs_resid <= abs_tol);
        std::cout << "[doc77c][energy-audit] done t=" << t_acc << " steps=" << step
                  << " cum_dH=" << cum_dH << " cum_Q=" << cum_Q << " cum_abs_Q=" << cum_abs_Q
                  << " cum_abs_resid=" << cum_abs_resid << " resid_rel=" << resid_rel
                  << " thermal_energy_pass=" << (energy_ok ? 1 : 0) << " csv=" << energy_csv << std::endl;
        if (!cli.summary_json.empty())
        {
            std::ofstream js(cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc77c_energy_audit\",\n"
               << "  \"doc87_r1a\": 1,\n"
               << "  \"dp\": " << dp << ",\n"
               << "  \"end_time\": " << t_end << ",\n"
               << "  \"t_acc\": " << t_acc << ",\n"
               << "  \"steps\": " << step << ",\n"
               << "  \"cum_dH\": " << cum_dH << ",\n"
               << "  \"cum_Q\": " << cum_Q << ",\n"
               << "  \"cum_abs_Q\": " << cum_abs_Q << ",\n"
               << "  \"cum_abs_resid\": " << cum_abs_resid << ",\n"
               << "  \"resid_rel\": " << resid_rel << ",\n"
               << "  \"resid_threshold\": 0.001,\n"
               << "  \"n_bad_dt\": " << n_bad_dt << ",\n"
               << "  \"contact_pair_factor_diag\": " << cli.contact_pair_factor << ",\n"
               << "  \"dirichlet_ghost\": " << (cli.dirichlet_ghost ? 1 : 0) << ",\n"
               << "  \"k_wall_factor\": " << cli.k_wall_factor << ",\n"
               << "  \"single_hot_wall\": 1,\n"
               << "  \"energy_csv\": \"" << energy_csv << "\",\n"
               << "  \"smoke_pass\": null,\n"
               << "  \"thermal_energy_pass\": " << (energy_ok ? 1 : 0) << ",\n"
               << "  \"force_update_consistency_pass\": null,\n"
               << "  \"mass_conservation_pass\": null,\n"
               << "  \"thermal_time_ok\": null,\n"
               << "  \"viscous_dt_checked\": null,\n"
               << "  \"zero_g_rest_pass\": null,\n"
               << "  \"hydrostatic_pass\": null,\n"
               << "  \"conduction_pass\": null,\n"
               << "  \"benchmark_pass\": null,\n"
               << "  \"passed\": " << (energy_ok ? 1 : 0) << "\n"
               << "}\n";
            std::cout << "[doc77c] summary_json=" << cli.summary_json << std::endl;
        }
        return energy_ok ? 0 : 2;
    }

    advance_to(cli.settle_time, "settle", false, false, false);
    hostZeroVelocity(particles);
    std::cout << "[doc77c] settle done; velocity zeroed; applying linear T + thin side Dirichlet" << std::endl
              << std::flush;

    hostInitLinearT(particles, x0, L, th, cli.t_c, t0);
    prepare_wall_bc();
    if (use_fluid_dirichlet)
    {
        hostApplySideDirichlet(particles, x0, L, shell, th, cli.t_c, t0, cli.rho0, cli.cp);
    }
    {
        const CavityMetrics m0 = metrics_now();
        std::cout << "[doc77c] after_linear_T_init Nu_grad_h=" << m0.Nu_grad_hot << " Nu_grad_c=" << m0.Nu_grad_cold
                  << " n_grad_h=" << m0.n_grad_hot << " n_grad_c=" << m0.n_grad_cold << " n_circ_h=" << m0.n_circ_hot
                  << " n_circ_c=" << m0.n_circ_cold << " shell=" << shell << " circ_band=" << circ_band
                  << " grad_layer=" << grad_layer << " fluid_dirichlet=" << (use_fluid_dirichlet ? 1 : 0)
                  << " (grad-Nu; expect ~1 for linear T)" << std::endl
                  << std::flush;
    }
    if (cli.enable_boussinesq && beta_force > TinyReal)
    {
        boussinesq.exec();
        const OphelieBoussinesqForceAudit fa =
            hostAuditBoussinesqForcePrior(particles, beta_force, t0, -cli.gravity_g, kOphelieTemperatureField);
        logOphelieBoussinesqForceAudit(fa, "doc77c_after_linear_T");
    }

    const Real t_obs_end = sv_t->getValue() + cli.end_time;
    t_late_gate = sv_t->getValue() + Real(0.5) * cli.end_time;
    const double t_fluid_before_obs = t_fluid_total;
    const double t_thermal_before_obs = t_thermal_total;
    advance_to(t_obs_end, "observe", cli.enable_boussinesq, true, true);

    CavityMetrics mend = metrics_now();
    const Real mean_hot = acc_n > 0 ? static_cast<Real>(acc_hot / acc_w) : mend.circ_hot_up;
    const Real mean_cold = acc_n > 0 ? static_cast<Real>(acc_cold / acc_w) : mend.circ_cold_dn;
    const Real mean_delta = mean_hot - mean_cold;
    const Real late_mean_Nu_grad_h = late_w > 0.0 ? static_cast<Real>(late_Nu_grad_h / late_w) : mend.Nu_grad_hot;
    const Real late_mean_Nu_grad_c = late_w > 0.0 ? static_cast<Real>(late_Nu_grad_c / late_w) : mend.Nu_grad_cold;
    const Real late_mean_Nu_bulk = late_w > 0.0 ? static_cast<Real>(late_Nu_bulk / late_w) : mend.Nu_bulk;
    if (cli.kernel_correction)
        fluid_linear_correction_matrix.exec();
    // Scheme-1 primary: match evolution flags (default Adami + B).
    const CorrectedWallNu ngc_end = hostCorrectedWallNuFromGradT(
        particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, cli.kernel_correction, cli.dirichlet_ghost);
    // End-snapshot ablation (diagnostic only; evolution unchanged): ±Adami × ±B.
    const CorrectedWallNu ngc_ab = hostCorrectedWallNuFromGradT(
        particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, /*use_B=*/true, /*ghost=*/true);
    const CorrectedWallNu ngc_a0 = hostCorrectedWallNuFromGradT(
        particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, /*use_B=*/false, /*ghost=*/true);
    const CorrectedWallNu ngc_0b = hostCorrectedWallNuFromGradT(
        particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, /*use_B=*/true, /*ghost=*/false);
    const CorrectedWallNu ngc_00 = hostCorrectedWallNuFromGradT(
        particles, wall_particles, x0, L, th, cli.t_c, grad_layer, h_sph, /*use_B=*/false, /*ghost=*/false);
    const Real late_mean_Nu_gradcorr_h =
        late_w > 0.0 ? static_cast<Real>(late_Nu_gradcorr_h / late_w) : ngc_end.Nu_hot;
    const Real late_mean_Nu_gradcorr_c =
        late_w > 0.0 ? static_cast<Real>(late_Nu_gradcorr_c / late_w) : ngc_end.Nu_cold;

    const Real Nu_wall_h_early =
        (t_thermal_early > TinyReal)
            ? static_cast<Real>((Q_hot_early / t_thermal_early) / (A_wall + TinyReal) / (q_scale + TinyReal))
            : Real(0);
    const Real Nu_wall_c_early =
        (t_thermal_early > TinyReal)
            ? static_cast<Real>((-Q_cold_early / t_thermal_early) / (A_wall + TinyReal) / (q_scale + TinyReal))
            : Real(0);
    const Real Nu_wall_h =
        (t_thermal_late > TinyReal)
            ? static_cast<Real>((Q_hot_late / t_thermal_late) / (A_wall + TinyReal) / (q_scale + TinyReal))
            : ((t_thermal_obs > TinyReal) ? static_cast<Real>((Q_hot_obs / t_thermal_obs) / (A_wall + TinyReal) /
                                                             (q_scale + TinyReal))
                                          : Real(0));
    const Real Nu_wall_c =
        (t_thermal_late > TinyReal)
            ? static_cast<Real>((-Q_cold_late / t_thermal_late) / (A_wall + TinyReal) / (q_scale + TinyReal))
            : ((t_thermal_obs > TinyReal) ? static_cast<Real>((-Q_cold_obs / t_thermal_obs) / (A_wall + TinyReal) /
                                                             (q_scale + TinyReal))
                                          : Real(0));

    const bool circ_ok = mend.circ_hot_up > TinyReal && mend.circ_cold_dn < -TinyReal;
    const Real circ_eps = Real(1.0e-4);
    const bool circ_rel_ok = mean_delta > circ_eps;
    const Real Nu_target = (cli.enable_boussinesq && beta_force > TinyReal) ? Nu_ref : Real(1.0);
    const Real Nu_gradcorr_err =
        std::abs(late_mean_Nu_gradcorr_h - Nu_target) / (Nu_target + TinyReal);
    /**
     * Doc87: formal Nu = raw Contact heat / (A k ΔT/L); pair_factor only diagnostic scaled Nu.
     * Conduction_pass requires raw Nu≈1 + balance + thermal time (not scaled factor).
     * Nu_gradcorr_* is Doc197 scheme-1 diagnostic (corrected ∇T); does not replace formal Contact.
     */
    const bool conduction_mode =
        cli.enable_diffusion && (!cli.enable_boussinesq) && cli.freeze_velocity;
    const Real Nu_wall_h_scaled = Nu_wall_h * cli.contact_pair_factor;
    const Real Nu_wall_c_scaled = Nu_wall_c * cli.contact_pair_factor;
    const Real Nu_use = cli.enable_diffusion ? Nu_wall_h : late_mean_Nu_grad_h;
    const Real Nu_err = std::abs(Nu_use - Nu_target) / (Nu_target + TinyReal);
    const Real Nu_bal_wall = std::abs(Nu_wall_h - Nu_wall_c) / (Nu_target + TinyReal);
    const Real Nu_bal = Nu_bal_wall;
    const Real U_scale = alpha / (L + TinyReal);
    const Real mean_delta_star = mean_delta / (U_scale + TinyReal);
    const Real circ_hot_star = mean_hot / (U_scale + TinyReal);
    const Real circ_cold_star = mean_cold / (U_scale + TinyReal);
    const Real U_rms_star = mend.U_rms / (U_scale + TinyReal);
    const Real u_err = std::abs(mend.u_max_star - u_ref) / (u_ref + TinyReal);
    const Real v_err = std::abs(mend.v_max_star - v_ref) / (v_ref + TinyReal);
    const bool no_runaway = mend.Ma < Real(0.25);
    const bool Nu_soft = (!cli.enable_diffusion) || (Nu_err < Real(0.35));
    const bool Nu_tight = cli.enable_diffusion && (Nu_err < Real(0.03));
    const double t_fluid_obs = t_fluid_total - t_fluid_before_obs;
    const double t_thermal_obs_phase = t_thermal_total - t_thermal_before_obs;
    const double t_time_rel =
        std::abs(t_fluid_obs - t_thermal_obs_phase) / (std::max(t_fluid_obs, t_thermal_obs_phase) + TinyReal);
    const bool thermal_time_ok = (!cli.enable_diffusion) || (t_time_rel < Real(1.0e-6));

    const bool smoke_pass =
        no_runaway && Nu_soft && (cli.enable_boussinesq ? circ_rel_ok : true);
    // Refresh buoyancy from final T before audit (Doc87 R1-D).
    if (cli.enable_boussinesq && beta_force > TinyReal)
        boussinesq.exec();
    const bool force_audit_pass =
        (!cli.enable_boussinesq) ||
        (hostAuditBoussinesqForcePrior(particles, beta_force, t0, -cli.gravity_g, kOphelieTemperatureField).rel_err <
         Real(1.0e-4));
    const bool conduction_pass = conduction_mode && thermal_time_ok && (Nu_err < Real(0.05)) && (Nu_bal < Real(0.05));
    const bool benchmark_pass = false;
    const bool thermal_energy_pass = false; // set only by --energy-audit path
    const bool passed = smoke_pass;

    std::cout << "[doc77c] end t=" << sv_t->getValue() << " t*=" << mend.t_star << " Nu_wall_h=" << Nu_wall_h
              << " Nu_wall_c=" << Nu_wall_c << " Nu_grad_h=" << late_mean_Nu_grad_h
              << " Nu_grad_c=" << late_mean_Nu_grad_c               << " Nu_gradcorr_h=" << late_mean_Nu_gradcorr_h
              << " Nu_gradcorr_c=" << late_mean_Nu_gradcorr_c << " Nu_gradcorr_err=" << Nu_gradcorr_err
              << " Nu_gc_end_AB=" << ngc_ab.Nu_hot << " Nu_gc_end_A0=" << ngc_a0.Nu_hot
              << " Nu_gc_end_0B=" << ngc_0b.Nu_hot << " Nu_gc_end_00=" << ngc_00.Nu_hot
              << " Nu_bulk=" << late_mean_Nu_bulk << " Nu_use=" << Nu_use << " Nu_target=" << Nu_target
              << " Nu_ref_dvd=" << Nu_ref << " Nu_err=" << Nu_err << " Nu_bal=" << Nu_bal
              << " Nu_tight_ok=" << (Nu_tight ? 1 : 0)
              << " freeze_vel=" << (cli.freeze_velocity ? 1 : 0) << " shell_bw=" << cli.boundary_width_factor
              << " circ_band_f=" << cli.circ_band_factor << " u*=" << mend.u_max_star << " u*_ref=" << u_ref
              << " u_err=" << u_err << " v*=" << mend.v_max_star << " v*_ref=" << v_ref << " v_err=" << v_err
              << " U_rms=" << mend.U_rms << " U_rms*=" << U_rms_star << " U_scale=" << U_scale
              << " circ_hot_uz=" << mend.circ_hot_up << " circ_cold_uz=" << mend.circ_cold_dn
              << " circ_ok=" << (circ_ok ? 1 : 0) << " mean_hot_uz=" << mean_hot << " mean_cold_uz=" << mean_cold
              << " mean_delta_uz=" << mean_delta << " mean_delta*=" << mean_delta_star
              << " circ_rel_ok=" << (circ_rel_ok ? 1 : 0) << " acc_n=" << acc_n
              << " no_runaway=" << (no_runaway ? 1 : 0) << " Nu_soft_ok=" << (Nu_soft ? 1 : 0)
              << " boussinesq=" << (cli.enable_boussinesq ? 1 : 0) << " smoke_pass=" << (smoke_pass ? 1 : 0)
              << " conduction_pass=" << (conduction_pass ? 1 : 0) << " thermal_time_ok=" << (thermal_time_ok ? 1 : 0)
              << " t_fluid_obs=" << t_fluid_obs << " t_thermal_obs=" << t_thermal_obs_phase
              << " t_time_rel=" << t_time_rel << " n_th_sub=" << n_th_sub_total
              << " n_th_capped=" << n_th_capped_total << " passed=" << (passed ? 1 : 0)
              << " dens_sum_normalize=" << (cli.dens_sum_normalize ? 1 : 0)
              << " background_hydro_p=" << (cli.background_hydro_p ? 1 : 0)
              << " kernel_correction=" << (cli.kernel_correction ? 1 : 0)
              << " Nu_wall_h_early=" << Nu_wall_h_early << " Nu_wall_h_late=" << Nu_wall_h
              << " hydrostatic_pass=null" << std::endl;
    std::cout << "[doc77c] reuse_path="
              << (cli.background_hydro_p ? "BackgroundHydroP+omitGravityForceCK+" : "GravityForceCK+")
              << "BoussinesqBuoyancyForceCK+Viscous+AcousticRiemann+"
                 "DensitySummation+OpheliePairwiseLaplaceCK+ApplyLaplaceDiffusionDtCK"
              << std::endl;
    std::cout << "[doc77c] circ_delta_uz=" << (mend.circ_hot_up - mend.circ_cold_dn) << std::endl;
    const Real q_hot = (t_thermal_late > TinyReal)
                           ? static_cast<Real>((Q_hot_late / t_thermal_late) / (A_wall + TinyReal))
                           : Real(0);
    const Real q_an = q_scale;
    std::cout << "[doc77c] wall_heat q_hot_raw=" << q_hot << " q_an_kdt_L=" << q_an
              << " q_ratio_raw=" << (q_hot / (q_an + TinyReal))
              << " Nu_wall_h_scaled=" << Nu_wall_h_scaled << " Nu_wall_c_scaled=" << Nu_wall_c_scaled
              << " L_eff=" << (L - Real(2) * shell) << " contact_pair_factor_diag=" << cli.contact_pair_factor
              << " Nu_bulk=" << late_mean_Nu_bulk << std::endl;
    std::cout << "[doc77c] note=3D_prism_doc84_p0_constraint_Nu_nonoverlap_circ_thermal_timebook" << std::endl;
    if (!cli.summary_json.empty())
    {
        std::ofstream js(cli.summary_json);
        js << "{\n"
           << "  \"case\": \"doc77c\",\n"
           << "  \"doc84_p0\": 1,\n"
           << "  \"dp\": " << dp << ",\n"
           << "  \"depth_factor\": " << cli.depth_factor << ",\n"
           << "  \"Ly\": " << Ly << ",\n"
           << "  \"alpha\": " << alpha << ",\n"
           << "  \"U_scale\": " << U_scale << ",\n"
           << "  \"Ra\": " << cli.Ra << ",\n"
           << "  \"boussinesq\": " << (cli.enable_boussinesq ? 1 : 0) << ",\n"
           << "  \"diffusion\": " << (cli.enable_diffusion ? 1 : 0) << ",\n"
           << "  \"freeze_velocity\": " << (cli.freeze_velocity ? 1 : 0) << ",\n"
           << "  \"boundary_width_factor\": " << cli.boundary_width_factor << ",\n"
           << "  \"shell\": " << shell << ",\n"
           << "  \"circ_band\": " << circ_band << ",\n"
           << "  \"grad_layer\": " << grad_layer << ",\n"
           << "  \"A_wall\": " << A_wall << ",\n"
           << "  \"mean_hot_uz\": " << mean_hot << ",\n"
           << "  \"mean_cold_uz\": " << mean_cold << ",\n"
           << "  \"mean_delta_uz\": " << mean_delta << ",\n"
           << "  \"mean_delta_star\": " << mean_delta_star << ",\n"
           << "  \"circ_hot_star\": " << circ_hot_star << ",\n"
           << "  \"circ_cold_star\": " << circ_cold_star << ",\n"
           << "  \"U_max\": " << mend.U_max << ",\n"
           << "  \"U_rms\": " << mend.U_rms << ",\n"
           << "  \"U_rms_star\": " << U_rms_star << ",\n"
           << "  \"circ_ok\": " << (circ_ok ? 1 : 0) << ",\n"
           << "  \"circ_rel_ok\": " << (circ_rel_ok ? 1 : 0) << ",\n"
           << "  \"Nu_wall_h\": " << Nu_wall_h << ",\n"
           << "  \"Nu_wall_c\": " << Nu_wall_c << ",\n"
           << "  \"Nu_wall_h_early\": " << Nu_wall_h_early << ",\n"
           << "  \"Nu_wall_c_early\": " << Nu_wall_c_early << ",\n"
           << "  \"Nu_wall_h_late\": " << Nu_wall_h << ",\n"
           << "  \"Nu_window_drift\": " << std::abs(Nu_wall_h - Nu_wall_h_early) << ",\n"
           << "  \"Nu_wall_h_scaled\": " << Nu_wall_h_scaled << ",\n"
           << "  \"Nu_wall_c_scaled\": " << Nu_wall_c_scaled << ",\n"
           << "  \"contact_pair_factor_diag\": " << cli.contact_pair_factor << ",\n"
           << "  \"dirichlet_ghost\": " << (cli.dirichlet_ghost ? 1 : 0) << ",\n"
           << "  \"k_wall_factor\": " << cli.k_wall_factor << ",\n"
           << "  \"wall_bw_factor\": " << cli.wall_bw_factor << ",\n"
           << "  \"corner_rim_dp\": " << cli.corner_rim_dp << ",\n"
           << "  \"corner_rim\": " << corner_rim << ",\n"
           << "  \"corner_rim_y\": " << (cli.corner_rim_y ? 1 : 0) << ",\n"
           << "  \"corner_rim_k_factor\": " << cli.corner_rim_k_factor << ",\n"
           << "  \"A_wall_geom\": " << A_wall_geom << ",\n"
           << "  \"A_wall_eff_diag\": " << A_wall_eff << ",\n"
           << "  \"Nu_wall_h_Aeff_diag\": " << (Nu_wall_h * A_wall / (A_wall_eff + TinyReal)) << ",\n"
           << "  \"A_wall_note\": \"formal Nu uses geometric A_wall; A_eff/Nu_Aeff diagnostic only\",\n"
           << "  \"pair_weight_regularization\": " << cli.pair_weight_regularization << ",\n"
           << "  \"thermal_energy_pass\": " << (thermal_energy_pass ? 1 : 0) << ",\n"
           << "  \"late_Nu_grad_h\": " << late_mean_Nu_grad_h << ",\n"
           << "  \"late_Nu_grad_c\": " << late_mean_Nu_grad_c << ",\n"
           << "  \"Nu_gradcorr_h\": " << late_mean_Nu_gradcorr_h << ",\n"
           << "  \"Nu_gradcorr_c\": " << late_mean_Nu_gradcorr_c << ",\n"
           << "  \"Nu_gradcorr_err\": " << Nu_gradcorr_err << ",\n"
           << "  \"Nu_gradcorr_note\": \"scheme1 SPHXsys cubic B-spline dW + Adami wall + optional B; diagnostic only\",\n"
           << "  \"Nu_gc_end_Adami_B\": " << ngc_ab.Nu_hot << ",\n"
           << "  \"Nu_gc_end_Adami_noB\": " << ngc_a0.Nu_hot << ",\n"
           << "  \"Nu_gc_end_noAdami_B\": " << ngc_0b.Nu_hot << ",\n"
           << "  \"Nu_gc_end_raw\": " << ngc_00.Nu_hot << ",\n"
           << "  \"late_Nu_bulk\": " << late_mean_Nu_bulk << ",\n"
           << "  \"Nu_use\": " << Nu_use << ",\n"
           << "  \"Nu_ref\": " << Nu_ref << ",\n"
           << "  \"Nu_target\": " << Nu_target << ",\n"
           << "  \"Nu_err\": " << Nu_err << ",\n"
           << "  \"Nu_bal\": " << Nu_bal << ",\n"
           << "  \"Nu_tight_ok\": " << (Nu_tight ? 1 : 0) << ",\n"
           << "  \"u_max_star\": " << mend.u_max_star << ",\n"
           << "  \"v_max_star\": " << mend.v_max_star << ",\n"
           << "  \"y_u_max_star\": " << mend.y_u_max << ",\n"
           << "  \"x_v_max_star\": " << mend.x_v_max << ",\n"
           << "  \"u_max_star_2\": " << mend.u_max_star_2 << ",\n"
           << "  \"y_u_max_star_2\": " << mend.y_u_max_2 << ",\n"
           << "  \"sort_every\": " << cli.sort_every << ",\n"
           << "  \"midline_note\": \"DVD-style: prefer upper-half |ux|max when dual jets; |uz|max on z~L/2 mid-y\",\n"
           << "  \"t_fluid_obs\": " << t_fluid_obs << ",\n"
           << "  \"t_thermal_obs\": " << t_thermal_obs_phase << ",\n"
           << "  \"t_time_rel\": " << t_time_rel << ",\n"
           << "  \"n_th_sub_total\": " << n_th_sub_total << ",\n"
           << "  \"n_th_capped_total\": " << n_th_capped_total << ",\n"
           << "  \"n_adv_steps\": " << n_adv_steps << ",\n"
           << "  \"smoke_pass\": " << (smoke_pass ? 1 : 0) << ",\n"
           << "  \"force_audit_pass\": " << (force_audit_pass ? 1 : 0) << ",\n"
           << "  \"dens_sum_normalize\": " << (cli.dens_sum_normalize ? 1 : 0) << ",\n"
           << "  \"background_hydro_p\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
           << "  \"kernel_correction\": " << (cli.kernel_correction ? 1 : 0) << ",\n"
           << "  \"sigma0_mean\": " << dens_sigma0_mean << ",\n"
           << "  \"hydrostatic_pass\": null,\n"
           << "  \"hydrostatic_pass_reason\": \"Doc91 OPEN; stage diagnostic not French SNR gate\",\n"
           << "  \"conduction_pass\": " << (conduction_pass ? 1 : 0) << ",\n"
           << "  \"thermal_time_ok\": " << (thermal_time_ok ? 1 : 0) << ",\n"
           << "  \"benchmark_pass\": " << (benchmark_pass ? 1 : 0) << ",\n"
           << "  \"no_runaway\": " << (no_runaway ? 1 : 0) << ",\n"
           << "  \"passed\": " << (passed ? 1 : 0) << ",\n"
           << "  \"compare\": {\"Nu_ref\": 1.118, \"Nu_err_soft\": 0.35, \"Nu_err_tight\": 0.03, "
              "\"conduction_Nu_err\": 0.05, \"u_ref\": 3.649, \"v_ref\": 3.697, "
              "\"primary_Nu\": \"contact_wall_heat_raw\",\n"
              "\"contact_pair_factor_diag\": " << cli.contact_pair_factor << ",\n"
              "\"wall_k_factor\": \"1e6_harmonic_to_2kf\",\n"
              "\"doc87\": \"formal_Nu_raw_scaled_diagnostic_only\"}\n"
           << "}\n";
        std::cout << "[doc77c] summary_json=" << cli.summary_json << std::endl;
    }
    {
        std::string mid_path = cli.out_csv;
        const auto slash = mid_path.find_last_of("/\\");
        if (slash != std::string::npos)
            mid_path = mid_path.substr(0, slash + 1) + "midline_profiles.csv";
        else
            mid_path = "midline_profiles.csv";
        hostExportMidlineProfiles(particles, mid_path, x0, z0, L, fluid_center[1], alpha, dp);
        std::cout << "[doc77c] midline_profiles=" << mid_path << " y_u_max=" << mend.y_u_max
                  << " x_v_max=" << mend.x_v_max << " u*=" << mend.u_max_star << " v*=" << mend.v_max_star
                  << std::endl;
    }
    return passed ? 0 : 2;
}
