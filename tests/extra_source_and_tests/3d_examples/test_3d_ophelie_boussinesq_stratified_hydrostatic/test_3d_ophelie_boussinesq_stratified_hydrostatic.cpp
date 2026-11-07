/**
 * @file test_3d_ophelie_boussinesq_stratified_hydrostatic.cpp
 * @brief Doc77-A: closed-box stratified hydrostatic with production BoussinesqBuoyancyForceCK.
 *
 * Freeze T(z)=Tref+a z (a>0, hot top / cold bottom). No Q, diffusion, or thermal BC.
 * Companion --mode=uniform: T≡Tref (zero buoyancy).
 *
 * Shared production path: GravityForceCK + BoussinesqBuoyancyForceCK + WCSPH acoustic/viscous.
 * Closed domain: DensityRegularization<Internal> (not FreeSurface).
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
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
    Real dp = 0.025;
    Real end_time = 0.5; /**< observation window after settle */
    Real settle_time = 0.5;
    Real box_l = 0.2;
    Real box_h = 0.2;
    Real rho0 = 1000.0;
    Real gravity_g = 9.81;
    Real beta = 1.0e-3;
    Real t_ref = 300.0;
    Real a_strat = 50.0; /**< K/m; ΔT = a*H = 10 K at H=0.2 */
    Real mu = 5.0;
    Real c0 = 20.0;
    Real u_ref = 0.5;
    bool enable_boussinesq = true;
    bool use_density_summation = true; /**< production-like; EOS hydrostatic profile is diagnostic only */
    /**
     * After DensitySummation, rewrite ρ(z) from analytic hydrostatic EOS (diagnostic only).
     * Doc87: never sync Mass in the reapply path; init Mass is fixed after setup.
     */
    bool reapply_hydro_rho = false;
    bool use_viscous_dt = false; /**< Doc87 R1-C: AdvectionViscousTimeStepCK */
    bool force_consistency_probe = false; /**< Doc87 R1-B: short Δu/dt vs Force bookkeeping */
    bool density_stage_probe = false; /**< Doc87 R2: freeze pos, stage dens_sum dump */
    /**
     * Doc97 P1: one real advection/acoustic update with stage momentum/density dumps.
     * Distinguishes dens_sum / ForcePrior / acoustic_1st / acoustic_2nd contributions.
     */
    bool hydrostatic_step_probe = false;
    /**
     * Doc87 R3-B: after dens_sum+reg, ρ_i ← ρ0 · ρ_sum_i / σ0_i
     * with σ0_i captured once on the initial lattice (reference-state kernel packing).
     * Diagnostic prototype — not a permanent force-ρ=ρ0.
     */
    bool dens_sum_normalize = false;
    /**
     * Doc87 R3-A / Doc98 P2b / Doc105 P1 (closed box): write base p_h=ρ0 g (H−z) into
     * BackgroundPressure (does NOT absorb stratification). Acoustic dens mode 1 +
     * PerturbationMomentum=1; GravityForceCK omitted. Buoyancy (Boussinesq) still applied;
     * stratified balance is carried by initial/restored p' via EOS density.
     * Free-surface French keeps legacy total-p + gravity (PerturbationMomentum default 0).
     */
    bool background_hydro_p = false;
    /**
     * Doc91: first-order kernel correction on WCSPH acoustic steps (LinearCorrectionCK),
     * with LinearCorrectionMatrixComplex refreshed each advection step. Not δ-SPH.
     */
    bool kernel_correction = false;
    /** Doc105 P0: host algebra selftest only (no SPH). */
    bool wall_dens_algebra_selftest = false;
    size_t force_probe_steps = 5;
    std::string mode = "stratified"; /**< stratified | uniform */
    std::string out_csv = "output/doc77a_stratified_hydrostatic.csv";
    std::string summary_json; /**< optional machine-readable end summary */
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
        if (take("--dp", cli.dp) || take("--end-time", cli.end_time) || take("--settle-time", cli.settle_time) ||
            take("--box-l", cli.box_l) || take("--box-h", cli.box_h) || take("--rho0", cli.rho0) ||
            take("--g", cli.gravity_g) || take("--beta", cli.beta) || take("--t-ref", cli.t_ref) ||
            take("--a", cli.a_strat) || take("--mu", cli.mu) || take("--c0", cli.c0) || take("--u-ref", cli.u_ref) ||
            take_s("--mode", cli.mode) || take_s("--out-csv", cli.out_csv) ||
            take_s("--summary-json", cli.summary_json))
        {
            continue;
        }
        if (arg.rfind("--force-probe-steps=", 0) == 0)
        {
            cli.force_probe_steps = static_cast<size_t>(std::atoll(arg.c_str() + 20));
            continue;
        }
        if (arg == "--no-boussinesq")
        {
            cli.enable_boussinesq = false;
            continue;
        }
        if (arg == "--no-density-summation")
        {
            cli.use_density_summation = false;
            continue;
        }
        if (arg == "--density-summation")
        {
            cli.use_density_summation = true;
            continue;
        }
        if (arg == "--reapply-hydro-rho")
        {
            cli.reapply_hydro_rho = true;
            continue;
        }
        if (arg == "--viscous-dt")
        {
            cli.use_viscous_dt = true;
            continue;
        }
        if (arg == "--force-consistency-probe")
        {
            cli.force_consistency_probe = true;
            continue;
        }
        if (arg == "--density-stage-probe")
        {
            cli.density_stage_probe = true;
            cli.use_density_summation = true;
            continue;
        }
        if (arg == "--hydrostatic-step-probe")
        {
            cli.hydrostatic_step_probe = true;
            continue;
        }
        if (arg == "--dens-sum-normalize")
        {
            cli.dens_sum_normalize = true;
            cli.use_density_summation = true;
            continue;
        }
        if (arg == "--background-hydro-p")
        {
            cli.background_hydro_p = true;
            cli.dens_sum_normalize = true;
            cli.use_density_summation = true;
            continue;
        }
        if (arg == "--doc105-p0-wall-dens-algebra")
        {
            cli.wall_dens_algebra_selftest = true;
            continue;
        }
        if (arg == "--kernel-correction")
        {
            cli.kernel_correction = true;
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            std::cout << "Doc77-A / Doc87 R1–R3 stratified hydrostatic\n"
                      << "  --dp= --settle-time= --end-time= --mode=stratified|uniform --no-boussinesq\n"
                      << "  --no-density-summation (continuity ρ; often unstable — diagnostic)\n"
                      << "  --reapply-hydro-rho (ρ only; Mass fixed after init — Doc87)\n"
                      << "  --viscous-dt  --force-consistency-probe  --force-probe-steps=N\n"
                      << "  --density-stage-probe  (freeze pos; dump dens_sum stages; use --g=0)\n"
                      << "  --hydrostatic-step-probe  (Doc97 P1: one update; stage Δu/drho dumps)\n"
                      << "  --dens-sum-normalize  (R3-B: ρ=ρ0·ρ_sum/σ0, σ0 from initial lattice)\n"
                      << "  --background-hydro-p  (closed-box P2b/P1: base p_h only + p'-momentum; omits GravityForceCK)\n"
                      << "  --doc105-p0-wall-dens-algebra  (Doc105 P0: host algebra selftest; no SPH)\n"
                      << "  --kernel-correction  (Doc91: LinearCorrectionCK on acoustic + B matrix)\n"
                      << "  --a= --beta= --c0= --mu= --out-csv= --summary-json=\n";
            return false;
        }
        std::cerr << "[doc77a] unknown arg: " << arg << std::endl;
        return false;
    }
    return true;
}

class FluidBoxShape : public ComplexShape
{
  public:
    FluidBoxShape(const std::string &name, const Vecd &center, const Vecd &halfsize) : ComplexShape(name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class ClosedWallShape : public ComplexShape
{
  public:
    ClosedWallShape(const std::string &name, const Vecd &center, const Vecd &fluid_half, Real bw)
        : ComplexShape(name)
    {
        const Vecd outer_half = fluid_half + Vecd(bw, bw, bw);
        add<GeometricShapeBox>(Transform(center), outer_half);
        subtract<GeometricShapeBox>(Transform(center), fluid_half);
    }
};

/** Base hydrostatic only (Doc105 P1): p_h = ρ0 g (H−ζ), ζ=z−z_bottom. Does NOT include stratification. */
inline Real analyticBaseHydroPressure(Real zeta, Real H, Real rho0, Real g)
{
    return rho0 * g * (H - zeta);
}

/**
 * Stratification perturbation pressure for closed-box P1 (absorb only base g):
 * p'(ζ) = ½ ρ0 g β a (ζ² − H²), so p'(H)=0 and dp'/dz = ρ0 g β (T−Tref).
 * Balances Boussinesq buoyancy when momentum uses p' and GravityForceCK is omitted.
 */
inline Real analyticStratPerturbationPressure(Real zeta, Real H, Real rho0, Real g, Real beta, Real a)
{
    return Real(0.5) * rho0 * g * beta * a * (zeta * zeta - H * H);
}

/** Legacy total hydrostatic: p = p_h_base + p'_strat (used when GravityForceCK cancels ∇p_total). */
inline Real analyticPressureFromBottom(Real zeta, Real H, Real rho0, Real g, Real beta, Real a)
{
    return analyticBaseHydroPressure(zeta, H, rho0, g) +
           analyticStratPerturbationPressure(zeta, H, rho0, g, beta, a);
}

/**
 * Set Density (+ optional Mass) from a hydrostatic pressure profile via EOS p = c0²(ρ−ρ0).
 * perturbation_only: use p'_strat only (Doc105 P1 + background_hydro_p).
 * else: use full analyticPressureFromBottom (legacy total-p path).
 */
inline void hostApplyStratifiedHydrostaticDensity(BaseParticles &particles, Real rho0, Real c0, Real g, Real beta,
                                                  Real a, Real z_bottom, Real H, bool sync_mass = true,
                                                  bool perturbation_only = false)
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
        const Real p = perturbation_only ? analyticStratPerturbationPressure(zeta, H, rho0, g, beta, a)
                                         : analyticPressureFromBottom(zeta, H, rho0, g, beta, a);
        rho[i] = rho0 + p * inv_c2;
        if (sync_mass)
            mass[i] = rho[i] * vol[i];
    }
    syncVariableToDevice<Real>(particles, "Density");
    if (sync_mass)
        syncVariableToDevice<Real>(particles, "Mass");
}

inline void hostPrescribeStableStratification(BaseParticles &particles, Real t_ref, Real a, Real z_bottom,
                                              bool uniform)
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
        const Real Ti = uniform ? t_ref : (t_ref + a * (pos[i][2] - z_bottom));
        temp[i] = Ti;
        delta_t[i] = Ti - t_ref;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

struct PressureFit
{
    Real c0 = 0.0; /**< constant */
    Real c1 = 0.0; /**< · ζ */
    Real c2 = 0.0; /**< · ζ² */
    Real c1_an = 0.0;
    Real c2_an = 0.0;
    Real c1_rel = 0.0;
    Real c2_rel = 0.0;
    Real rms_resid = 0.0;
};

/**
 * Least-squares p ≈ c0 + c1 ζ + c2 ζ² (ζ from bottom).
 * use_background_field: fit BackgroundPressure (Doc105 P1: base hydro only → c2_an=0);
 * else fit EOS p = c²(ρ−ρ0). When eos_is_perturbation (bgp path), analytic target is p'_strat.
 */
inline PressureFit hostFitPressureQuadratic(BaseParticles &particles, Real rho0, Real c_sound, Real g, Real beta,
                                            Real a, Real z_bottom, Real H, bool use_background_field = false,
                                            bool eos_is_perturbation = false)
{
    PressureFit f;
    if (use_background_field)
    {
        f.c1_an = -rho0 * g;
        f.c2_an = Real(0); // Doc105 P1: p_h absorbs base g only
    }
    else if (eos_is_perturbation)
    {
        f.c1_an = Real(0);
        f.c2_an = Real(0.5) * rho0 * g * beta * a;
    }
    else
    {
        f.c1_an = -rho0 * g;
        f.c2_an = Real(0.5) * rho0 * g * beta * a;
    }
    syncVariableToHost<Vecd>(particles, "Position");
    if (use_background_field)
        syncVariableToHost<Real>(particles, "BackgroundPressure");
    else
        syncVariableToHost<Real>(particles, "Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *rho = use_background_field ? nullptr : particles.getVariableDataByName<Real>("Density");
    const Real *p_bg =
        use_background_field ? particles.getVariableDataByName<Real>("BackgroundPressure") : nullptr;
    const size_t n = particles.TotalRealParticles();
    // Normal equations for [1, ζ, ζ²]
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double b0 = 0, b1 = 0, b2 = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double z = static_cast<double>(std::min(H, std::max(Real(0), pos[i][2] - z_bottom)));
        const double p = use_background_field ? static_cast<double>(p_bg[i])
                                              : static_cast<double>(c_sound * c_sound * (rho[i] - rho0));
        const double z2 = z * z;
        s0 += 1.0;
        s1 += z;
        s2 += z2;
        s3 += z2 * z;
        s4 += z2 * z2;
        b0 += p;
        b1 += p * z;
        b2 += p * z2;
    }
    // Solve 3x3 with Cramer's / Gaussian elimination
    const double A[3][4] = {{s0, s1, s2, b0}, {s1, s2, s3, b1}, {s2, s3, s4, b2}};
    double M[3][4];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            M[i][j] = A[i][j];
    for (int col = 0; col < 3; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < 3; ++r)
            if (std::abs(M[r][col]) > std::abs(M[piv][col]))
                piv = r;
        for (int j = 0; j < 4; ++j)
            std::swap(M[col][j], M[piv][j]);
        const double div = M[col][col] + 1e-30;
        for (int j = col; j < 4; ++j)
            M[col][j] /= div;
        for (int r = 0; r < 3; ++r)
        {
            if (r == col)
                continue;
            const double fct = M[r][col];
            for (int j = col; j < 4; ++j)
                M[r][j] -= fct * M[col][j];
        }
    }
    f.c0 = static_cast<Real>(M[0][3]);
    f.c1 = static_cast<Real>(M[1][3]);
    f.c2 = static_cast<Real>(M[2][3]);
    f.c1_rel = std::abs(f.c1 - f.c1_an) / (std::abs(f.c1_an) + TinyReal);
    f.c2_rel = std::abs(f.c2 - f.c2_an) / (std::abs(f.c2_an) + TinyReal);
    double rs = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double z = static_cast<double>(std::min(H, std::max(Real(0), pos[i][2] - z_bottom)));
        const double p = use_background_field ? static_cast<double>(p_bg[i])
                                              : static_cast<double>(c_sound * c_sound * (rho[i] - rho0));
        const double pred = f.c0 + f.c1 * z + f.c2 * z * z;
        rs += (p - pred) * (p - pred);
    }
    f.rms_resid = static_cast<Real>(std::sqrt(rs / static_cast<double>(n)));
    return f;
}

struct HydroStats
{
    Real u_max = 0.0;
    Real u_rms = 0.0;
    Real p_l2_rel = 0.0;
    Real p_linf = 0.0;
    Real rho_mean = 0.0;
    size_t n = 0;
};

inline HydroStats hostHydroStats(BaseParticles &particles, Real rho0, Real c0, Real g, Real beta, Real a,
                                 Real z_bottom, Real H)
{
    HydroStats s;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    s.n = particles.TotalRealParticles();
    double u2 = 0.0, p_err2 = 0.0, p_ref2 = 0.0, rho_sum = 0.0;
    Real u_max = 0.0, p_linf = 0.0;
    // Subtract mean pressure so EOS offset does not dominate.
    double p_num_mean = 0.0, p_an_mean = 0.0;
    for (size_t i = 0; i < s.n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p_an = analyticPressureFromBottom(zeta, H, rho0, g, beta, a);
        const Real p_num = c0 * c0 * (rho[i] - rho0);
        p_an_mean += p_an;
        p_num_mean += p_num;
    }
    p_an_mean /= static_cast<double>(s.n);
    p_num_mean /= static_cast<double>(s.n);
    for (size_t i = 0; i < s.n; ++i)
    {
        const Real un = vel[i].norm();
        u_max = std::max(u_max, un);
        u2 += static_cast<double>(un) * un;
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p_an = analyticPressureFromBottom(zeta, H, rho0, g, beta, a) - static_cast<Real>(p_an_mean);
        const Real p_num = c0 * c0 * (rho[i] - rho0) - static_cast<Real>(p_num_mean);
        const Real pe = p_num - p_an;
        p_err2 += static_cast<double>(pe) * pe;
        p_ref2 += static_cast<double>(p_an) * p_an;
        p_linf = std::max(p_linf, std::abs(pe));
        rho_sum += rho[i];
    }
    s.u_max = u_max;
    s.u_rms = static_cast<Real>(std::sqrt(u2 / static_cast<double>(s.n)));
    s.p_l2_rel = static_cast<Real>(std::sqrt(p_err2) / (std::sqrt(p_ref2) + TinyReal));
    s.p_linf = p_linf;
    s.rho_mean = static_cast<Real>(rho_sum / static_cast<double>(s.n));
    return s;
}

inline void hostWritePressureProfileCsv(const std::string &path, BaseParticles &particles, Real rho0, Real c0,
                                        Real g, Real beta, Real a, Real z_bottom, Real H)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    std::ofstream ofs(path);
    ofs << "z,zeta,rho,p_eos,p_analytic,p_err\n";
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        const Real p_an = analyticPressureFromBottom(zeta, H, rho0, g, beta, a);
        const Real p_eos = c0 * c0 * (rho[i] - rho0);
        ofs << pos[i][2] << "," << zeta << "," << rho[i] << "," << p_eos << "," << p_an << "," << (p_eos - p_an)
            << "\n";
    }
}

/**
 * Doc87 R1-B: z-binned bookkeeping.
 * a_external_gb_z = grav+buoy only (NOT proof of pressure imbalance).
 * When Force/ForcePrior/ViscousForce present, also report a_force_prior / a_force / a_total_from_fields.
 */
inline void hostWriteAccelProfileCsv(const std::string &path, BaseParticles &particles, Real g, Real beta, Real a,
                                     Real t_ref, Real z_bottom, Real H, size_t n_bins = 20)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Vecd>(particles, "GravityForceCK");
    syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Vecd *f_g = particles.getVariableDataByName<Vecd>("GravityForceCK");
    const Vecd *f_b = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
    const Vecd *f_prior = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *f_ac = particles.getVariableDataByName<Vecd>("Force");
    const Vecd *f_visc = particles.getVariableDataByName<Vecd>("ViscousForce");
    const size_t n = particles.TotalRealParticles();
    n_bins = std::max<size_t>(n_bins, 4);
    std::vector<double> w(n_bins, 0.0), uz(n_bins, 0.0), ag(n_bins, 0.0), ab(n_bins, 0.0), dT(n_bins, 0.0);
    std::vector<double> afp(n_bins, 0.0), af(n_bins, 0.0), av(n_bins, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][2] - z_bottom));
        size_t b = static_cast<size_t>(zeta / (H + TinyReal) * static_cast<Real>(n_bins));
        if (b >= n_bins)
            b = n_bins - 1;
        const double m = static_cast<double>(mass[i]);
        const double inv_m = 1.0 / (m + TinyReal);
        w[b] += m;
        uz[b] += m * static_cast<double>(vel[i][2]);
        dT[b] += m * static_cast<double>(temp[i] - t_ref);
        ag[b] += m * static_cast<double>(f_g[i][2]) * inv_m;
        ab[b] += m * static_cast<double>(f_b[i][2]) * inv_m;
        afp[b] += m * static_cast<double>(f_prior[i][2]) * inv_m;
        af[b] += m * static_cast<double>(f_ac[i][2]) * inv_m;
        av[b] += m * static_cast<double>(f_visc[i][2]) * inv_m;
    }
    std::ofstream ofs(path);
    ofs << "bin,zeta_c,uz_mean,a_grav_z,a_buoy_z,a_external_gb_z,a_force_prior_z,a_force_z,"
           "a_viscous_z,a_total_fields_z,dT_mean,a_buoy_analytic,g\n";
    for (size_t b = 0; b < n_bins; ++b)
    {
        if (w[b] <= TinyReal)
            continue;
        const Real zeta_c = (static_cast<Real>(b) + Real(0.5)) * H / static_cast<Real>(n_bins);
        const double uzm = uz[b] / w[b];
        const double agm = ag[b] / w[b];
        const double abm = ab[b] / w[b];
        const double afpm = afp[b] / w[b];
        const double afm = af[b] / w[b];
        const double avm = av[b] / w[b];
        const double dTm = dT[b] / w[b];
        const double ab_an = static_cast<double>(beta * a * zeta_c * g);
        ofs << b << "," << zeta_c << "," << uzm << "," << agm << "," << abm << "," << (agm + abm) << "," << afpm << ","
            << afm << "," << avm << "," << (afpm + afm) << "," << dTm << "," << ab_an << "," << g << "\n";
    }
}

/** Doc87 R2: per-particle dens_sum / EOS snapshot (positions frozen). */
inline void hostDumpDensityStageCsv(const std::string &path, BaseParticles &particles, Real rho0, Real c0,
                                    Real z_bottom, Real H, Real L, const std::string &stage)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const bool has_sum = true;
    syncVariableToHost<Real>(particles, "DensitySummation");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    std::ofstream ofs(path);
    ofs << "stage,i,x,y,z,zeta,dist_wall,rho,rho_sum,mass,vol,delta_rho,delta_p_eos\n";
    double sum_dr = 0.0, sum_dr2 = 0.0, max_adr = 0.0;
    size_t n_near = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        const Real zeta = z - z_bottom;
        const Real dist_wall = std::min({x + Real(0.5) * L, Real(0.5) * L - x, y + Real(0.5) * L, Real(0.5) * L - y,
                                         zeta, H - zeta});
        const Real dr = rho[i] - rho0;
        const Real dp = c0 * c0 * dr;
        sum_dr += static_cast<double>(dr);
        sum_dr2 += static_cast<double>(dr * dr);
        max_adr = std::max(max_adr, std::abs(static_cast<double>(dr)));
        if (dist_wall < Real(2.5) * (H / Real(8))) // rough near-wall flag; dp sets finer later
            ++n_near;
        ofs << stage << "," << i << "," << x << "," << y << "," << z << "," << zeta << "," << dist_wall << ","
            << rho[i] << "," << (has_sum ? rho_sum[i] : rho[i]) << "," << mass[i] << "," << vol[i] << "," << dr << ","
            << dp << "\n";
    }
    const double mean = sum_dr / static_cast<double>(n);
    const double rms = std::sqrt(sum_dr2 / static_cast<double>(n));
    std::cout << "[doc77a][density-stage] " << stage << " n=" << n << " mean_dr=" << mean << " rms_dr=" << rms
              << " max_|dr|=" << max_adr << " n_near~=" << n_near << " csv=" << path << std::endl;
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

/** Capture per-particle σ0 = DensitySummation on the reference (initial) configuration. */
inline void hostCaptureDensSumSigma0(BaseParticles &particles, std::vector<Real> &sigma0, double &sigma0_mean)
{
    syncVariableToHost<Real>(particles, "DensitySummation");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    sigma0.resize(n);
    sigma0_mean = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sigma0[i] = std::max(rho_sum[i], TinyReal);
        sigma0_mean += static_cast<double>(sigma0[i]);
    }
    sigma0_mean /= static_cast<double>(n);
}

/**
 * R3-B: ρ_i = ρ0 · ρ_sum_i / σ0_i  (Mass unchanged; Vol updated later by AdvectionStepSetup).
 * Uses Density after Internal regularization (=ρ_sum) or DensitySummation if preferred.
 */
inline void hostApplyDensSumNormalize(BaseParticles &particles, Real rho0, const std::vector<Real> &sigma0)
{
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensitySummation");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *rho_sum = particles.getVariableDataByName<Real>("DensitySummation");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real s0 = (i < sigma0.size()) ? sigma0[i] : TinyReal;
        rho[i] = rho0 * (rho_sum[i] / (s0 + TinyReal));
    }
    syncVariableToDevice<Real>(particles, "Density");
}

/** Doc97 P1: stage snapshot of hydrostatic update residuals (inner vs wall band). */
struct HydroStepStageStats
{
    size_t n = 0, n_inner = 0, n_wall = 0;
    double uz_rms = 0, uz_max = 0;
    double a_prior_z_rms = 0, a_force_z_rms = 0, a_net_z_rms = 0;
    double a_prior_z_mean = 0, a_force_z_mean = 0, a_net_z_mean = 0;
    double drho_dt_rms = 0, drho_dt_max = 0;
    double p_bg_rms = 0, p_eos_rms = 0, p_total_rms = 0;
    double p_bg_mean = 0, p_eos_mean = 0;
    double rho_mean = 0, delta_rho_rms = 0;
    double a_net_inner_rms = 0, a_net_wall_rms = 0;
    double drho_dt_inner_rms = 0, drho_dt_wall_rms = 0;
};

inline HydroStepStageStats hostHydroStepStageStats(BaseParticles &particles, Real rho0, Real c0, Real L,
                                                   Real z_bottom, Real H, Real wall_band)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Force");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "DensityChangeRate");
    syncVariableToHost<Real>(particles, "Pressure");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
    const Real *p = particles.getVariableDataByName<Real>("Pressure");
    const Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    HydroStepStageStats s;
    s.n = particles.TotalRealParticles();
    double s_uz2 = 0, s_ap2 = 0, s_af2 = 0, s_an2 = 0, s_ap = 0, s_af = 0, s_an = 0;
    double s_dr2 = 0, s_pbg2 = 0, s_peos2 = 0, s_pt2 = 0, s_pbg = 0, s_peos = 0, s_rho = 0, s_drho2 = 0;
    double s_an_i2 = 0, s_an_w2 = 0, s_dr_i2 = 0, s_dr_w2 = 0;
    for (size_t i = 0; i < s.n; ++i)
    {
        const double m = static_cast<double>(mass[i]) + 1e-300;
        const double ap = static_cast<double>(fp[i][2]) / m;
        const double af = static_cast<double>(f[i][2]) / m;
        const double an = ap + af;
        const double uz = static_cast<double>(vel[i][2]);
        const double pbg = static_cast<double>(p_bg[i]);
        const double peos = static_cast<double>(c0 * c0 * (rho[i] - rho0));
        const double pt = static_cast<double>(p[i]);
        const double dr = static_cast<double>(drho[i]);
        const double drho_err = static_cast<double>(rho[i] - rho0);
        s_uz2 += uz * uz;
        s.uz_max = std::max(s.uz_max, std::abs(uz));
        s_ap2 += ap * ap;
        s_af2 += af * af;
        s_an2 += an * an;
        s_ap += ap;
        s_af += af;
        s_an += an;
        s_dr2 += dr * dr;
        s.drho_dt_max = std::max(s.drho_dt_max, std::abs(dr));
        s_pbg2 += pbg * pbg;
        s_peos2 += peos * peos;
        s_pt2 += pt * pt;
        s_pbg += pbg;
        s_peos += peos;
        s_rho += static_cast<double>(rho[i]);
        s_drho2 += drho_err * drho_err;
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        const Real dist = std::min({x + Real(0.5) * L, Real(0.5) * L - x, y + Real(0.5) * L, Real(0.5) * L - y,
                                    z - z_bottom, H - (z - z_bottom)});
        if (dist < wall_band)
        {
            ++s.n_wall;
            s_an_w2 += an * an;
            s_dr_w2 += dr * dr;
        }
        else
        {
            ++s.n_inner;
            s_an_i2 += an * an;
            s_dr_i2 += dr * dr;
        }
    }
    const double inv = 1.0 / static_cast<double>(std::max<size_t>(s.n, 1));
    s.uz_rms = std::sqrt(s_uz2 * inv);
    s.a_prior_z_rms = std::sqrt(s_ap2 * inv);
    s.a_force_z_rms = std::sqrt(s_af2 * inv);
    s.a_net_z_rms = std::sqrt(s_an2 * inv);
    s.a_prior_z_mean = s_ap * inv;
    s.a_force_z_mean = s_af * inv;
    s.a_net_z_mean = s_an * inv;
    s.drho_dt_rms = std::sqrt(s_dr2 * inv);
    s.p_bg_rms = std::sqrt(s_pbg2 * inv);
    s.p_eos_rms = std::sqrt(s_peos2 * inv);
    s.p_total_rms = std::sqrt(s_pt2 * inv);
    s.p_bg_mean = s_pbg * inv;
    s.p_eos_mean = s_peos * inv;
    s.rho_mean = s_rho * inv;
    s.delta_rho_rms = std::sqrt(s_drho2 * inv);
    s.a_net_inner_rms = s.n_inner ? std::sqrt(s_an_i2 / static_cast<double>(s.n_inner)) : 0.0;
    s.a_net_wall_rms = s.n_wall ? std::sqrt(s_an_w2 / static_cast<double>(s.n_wall)) : 0.0;
    s.drho_dt_inner_rms = s.n_inner ? std::sqrt(s_dr_i2 / static_cast<double>(s.n_inner)) : 0.0;
    s.drho_dt_wall_rms = s.n_wall ? std::sqrt(s_dr_w2 / static_cast<double>(s.n_wall)) : 0.0;
    return s;
}

inline void hostAppendHydroStepStageCsv(std::ofstream &ofs, const std::string &stage, const HydroStepStageStats &s)
{
    ofs << stage << "," << s.n << "," << s.n_inner << "," << s.n_wall << "," << s.uz_rms << "," << s.uz_max << ","
        << s.a_prior_z_rms << "," << s.a_force_z_rms << "," << s.a_net_z_rms << "," << s.a_prior_z_mean << ","
        << s.a_force_z_mean << "," << s.a_net_z_mean << "," << s.drho_dt_rms << "," << s.drho_dt_max << ","
        << s.p_bg_rms << "," << s.p_eos_rms << "," << s.p_total_rms << "," << s.p_bg_mean << "," << s.p_eos_mean
        << "," << s.rho_mean << "," << s.delta_rho_rms << "," << s.a_net_inner_rms << "," << s.a_net_wall_rms << ","
        << s.drho_dt_inner_rms << "," << s.drho_dt_wall_rms << "\n";
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

inline void hostZeroRealField(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Real>(particles, name);
    Real *v = particles.getVariableDataByName<Real>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Real(0);
    syncVariableToDevice<Real>(particles, name);
}

/** Analytic base p_h(z) → BackgroundPressure (Doc105 P1: absorb only ρ0 g, not stratification). */
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

/** Doc98 P2b: Acoustic momentum uses p' when flag=1 (closed-box with BackgroundPressure). */
inline void hostSetPerturbationMomentum(BaseParticles &particles, Real flag)
{
    syncVariableToHost<Real>(particles, "PerturbationMomentum");
    Real *pm = particles.getVariableDataByName<Real>("PerturbationMomentum");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        pm[i] = flag;
    syncVariableToDevice<Real>(particles, "PerturbationMomentum");
}

/** Doc105: 0=legacy dens (Adami total/mom), 1=p' dens with independent wall reference. */
inline void hostSetAcousticDensMode(BaseParticles &particles, Real mode)
{
    syncVariableToHost<Real>(particles, "AcousticDensMode");
    Real *m = particles.getVariableDataByName<Real>("AcousticDensMode");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        m[i] = mode;
    syncVariableToDevice<Real>(particles, "AcousticDensMode");
}

/** Doc105: reference body acceleration for wall dens C_ref (e.g. (0,0,-g)). */
inline void hostSetReferenceBodyAcceleration(BaseParticles &particles, const Vecd &a_ref)
{
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    Vecd *a = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        a[i] = a_ref;
    syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
}

/**
 * Doc105 P0 algebra probe (host): wall dens pressure jump for legacy vs broken vs mode-1.
 * Does not execute SPH; verifies C_ext / C_ref formulas used in AcousticStep1stHalf Wall.
 *
 * Sign convention: DissipativeUJump(p_i − p_j); Adami p_j = p_i + C ⇒ jump = −C.
 * Mode-1 closed-box P2b (no g in ForcePrior): C_ref must be 0, else spurious +C_ref dens.
 * Mode-1 with g in ForcePrior: C_ref = C_g, jump = −C_extra for dynamic extras.
 */
inline void hostDoc105P0WallDensAlgebraSelftest(std::ostream &out, bool &pass)
{
    const Real rho = Real(1000);
    const Real r = Real(0.01);
    const Real C_g = rho * r * Real(9.81); // hydrostatic contribution alone
    const Real C_extra = rho * r * Real(1.0);
    // Keep pressures O(10^2) so float Real does not lose C_g when added to p_i.
    const Real p_i = Real(100);
    const Real p_h_i = Real(20);
    const Real p_i_prime = p_i - p_h_i;
    const Real tol = Real(1.0e-5);

    // Broken Doc98 wall dens (same C in total and p_h ghost):
    const Real jump_broken = p_i_prime - ((p_i + C_g) - (p_h_i + C_g));

    // Legacy Adami dens on total: jump = −C
    const Real jump_legacy = p_i - (p_i + C_g);

    // Mode-1 + g in ForcePrior (C_ext=C_g+C_extra, C_ref=C_g): jump = −C_extra
    const Real jump_m1_dyn = p_i_prime - ((p_i + (C_g + C_extra)) - (p_h_i + C_g));

    // Mode-1 closed-box P2b (C_ext=0, C_ref=0): jump = 0
    const Real jump_m1_p2b = p_i_prime - ((p_i + Real(0)) - (p_h_i + Real(0)));

    // Mode-1 P2b mis-wired (C_ext=0, C_ref=C_g): spurious dens jump = +C_g (must NOT use)
    const Real jump_m1_p2b_bad = p_i_prime - ((p_i + Real(0)) - (p_h_i + C_g));

    const bool ok_broken_zero = std::abs(jump_broken) <= tol * (p_i + TinyReal);
    const bool ok_legacy = std::abs(jump_legacy + C_g) <= tol * (C_g + TinyReal);
    const bool ok_m1_dyn = std::abs(jump_m1_dyn + C_extra) <= tol * (C_extra + TinyReal);
    const bool ok_m1_p2b = std::abs(jump_m1_p2b) <= tol * (p_i + TinyReal);
    const bool ok_m1_p2b_bad = std::abs(jump_m1_p2b_bad - C_g) <= tol * (C_g + TinyReal);
    pass = ok_broken_zero && ok_legacy && ok_m1_dyn && ok_m1_p2b && ok_m1_p2b_bad;

    out << "{\n"
        << "  \"case\": \"doc105_p0_wall_dens_algebra\",\n"
        << "  \"jump_broken\": " << jump_broken << ",\n"
        << "  \"jump_legacy\": " << jump_legacy << ",\n"
        << "  \"jump_mode1_dynamic\": " << jump_m1_dyn << ",\n"
        << "  \"jump_mode1_p2b\": " << jump_m1_p2b << ",\n"
        << "  \"jump_mode1_p2b_bad_Cref\": " << jump_m1_p2b_bad << ",\n"
        << "  \"C_g\": " << C_g << ",\n"
        << "  \"C_extra\": " << C_extra << ",\n"
        << "  \"expect_broken_zero\": " << (ok_broken_zero ? 1 : 0) << ",\n"
        << "  \"expect_legacy_neg_C\": " << (ok_legacy ? 1 : 0) << ",\n"
        << "  \"expect_mode1_dyn_neg_Cextra\": " << (ok_m1_dyn ? 1 : 0) << ",\n"
        << "  \"expect_mode1_p2b_zero\": " << (ok_m1_p2b ? 1 : 0) << ",\n"
        << "  \"expect_mode1_p2b_bad_Cg\": " << (ok_m1_p2b_bad ? 1 : 0) << ",\n"
        << "  \"passed\": " << (pass ? 1 : 0) << "\n"
        << "}\n";
}

/** Snapshot particle uz before acoustic_1st; after, compare Δuz/dt vs (Fprior+F)/m. */
inline void hostCompareDuDtVsForce(BaseParticles &particles, const std::vector<Real> &uz_before, Real dt_ac,
                                   double &du_dt_rms, double &a_force_rms, double &rel_err)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "Force");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    const Vecd *f = particles.getVariableDataByName<Vecd>("Force");
    const size_t n = particles.TotalRealParticles();
    double s_du2 = 0.0, s_a2 = 0.0, s_e2 = 0.0, s_ref2 = 0.0;
    const double inv_dt = 1.0 / (static_cast<double>(dt_ac) + 1e-300);
    for (size_t i = 0; i < n; ++i)
    {
        const double du_dt = (static_cast<double>(vel[i][2]) - static_cast<double>(uz_before[i])) * inv_dt;
        // Note: after acoustic_1st Update, Force still holds pressure contribution used in Update;
        // acoustic_2nd has not yet overwritten Force with dissipation.
        const double a_f = static_cast<double>((fp[i][2] + f[i][2]) / (mass[i] + TinyReal));
        s_du2 += du_dt * du_dt;
        s_a2 += a_f * a_f;
        const double e = du_dt - a_f;
        s_e2 += e * e;
        s_ref2 += a_f * a_f;
    }
    du_dt_rms = std::sqrt(s_du2 / static_cast<double>(n));
    a_force_rms = std::sqrt(s_a2 / static_cast<double>(n));
    rel_err = std::sqrt(s_e2) / (std::sqrt(s_ref2) + 1e-30);
}
} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
    {
        return 1;
    }
    if (cli.wall_dens_algebra_selftest)
    {
        bool pass = false;
        if (!cli.summary_json.empty())
        {
            std::ofstream js(cli.summary_json);
            hostDoc105P0WallDensAlgebraSelftest(js, pass);
            std::cout << "[doc105_p0] wall_dens_algebra summary_json=" << cli.summary_json
                      << " passed=" << (pass ? 1 : 0) << std::endl;
        }
        else
        {
            hostDoc105P0WallDensAlgebraSelftest(std::cout, pass);
        }
        return pass ? 0 : 2;
    }
    const bool uniform = (cli.mode == "uniform");
    const Real a_eff = (uniform || !cli.enable_boussinesq) ? Real(0) : cli.a_strat;
    const Real beta_eff = cli.enable_boussinesq ? cli.beta : Real(0);

    const Real L = cli.box_l;
    const Real H = cli.box_h;
    const Real dp = cli.dp;
    const Real BW = Real(3) * dp;
    const Vecd fluid_center(0.0, 0.0, Real(0.5) * H);
    const Vecd fluid_half(Real(0.5) * L, Real(0.5) * L, Real(0.5) * H);
    const Real z_bottom = 0.0;

    const BoundingBoxd system_bounds(fluid_center - fluid_half - Vecd(BW, BW, BW),
                                     fluid_center + fluid_half + Vecd(BW, BW, BW));
    SPHSystem sph_system(system_bounds, dp);

    FluidBody fluid(sph_system, makeShared<FluidBoxShape>("FluidBody", fluid_center, fluid_half));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(cli.rho0, cli.c0);
    fluid.addMaterialProperty<Viscosity>(cli.mu);
    fluid.generateParticles<BaseParticles, Lattice>();

    SolidBody wall(sph_system, makeShared<ClosedWallShape>("WallBoundary", fluid_center, fluid_half, BW));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall_contact(fluid, {&wall});

    BaseParticles &particles = fluid.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(particles, cli.t_ref);
    hostPrescribeStableStratification(particles, cli.t_ref, cli.a_strat, z_bottom, uniform);
    // Init: set ρ(z) and Mass once. Doc87: Mass stays fixed for the rest of the run.
    // Doc105 P1 + bgp: ρ from p'_strat only (base g lives in BackgroundPressure).
    hostApplyStratifiedHydrostaticDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H,
                                          /*sync_mass=*/true, /*perturbation_only=*/cli.background_hydro_p);
    const double mass_sum_init = hostTotalMass(particles);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cell_linked_list(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> fluid_update_relation(fluid_inner, fluid_wall_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(fluid);

    Gravity gravity(Vecd(0.0, 0.0, -cli.gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(fluid, gravity);
    StateDynamics<MainExecutionPolicy, BoussinesqBuoyancyForceCK> boussinesq_force(
        fluid, gravity, beta_eff, cli.t_ref, kOphelieTemperatureField);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> advection_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_position(fluid);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(fluid_inner, 0.5), fluid_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
        acoustic_1st(fluid_inner, fluid_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalf<
                                                   Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd_wall(fluid_wall_contact);
    acoustic_2nd.addPostContactInteraction(acoustic_2nd_wall);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        acoustic_1st_corr(fluid_inner, fluid_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        acoustic_2nd_corr(fluid_inner, fluid_wall_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>, Contact<>>>
        density_summation(fluid_inner, fluid_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, Internal>> density_reg(fluid);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK> viscous_force(fluid_inner,
                                                                                                    fluid_wall_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> advection_dt(fluid, cli.u_ref);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionViscousTimeStepCK> advection_visc_dt(
        fluid, cli.u_ref);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> acoustic_dt(fluid);
    ReduceDynamicsCK<MainExecutionPolicy, OphelieVelocityMaxReduceCK> u_max_reduce(fluid);

    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");

    wall_normal.exec();
    fluid_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    fluid_update_relation.exec();
    // Closed-box P2b: p'-momentum + omit GravityForceCK when BackgroundPressure is active.
    // Doc105 dens mode 1: dens on p'. ReferenceBodyAcceleration = 0 here because gravity
    // is absorbed into BackgroundPressure and is NOT in ForcePrior — setting C_ref from g
    // would inject spurious wall dens (C_ext=0, C_ref=C_g ⇒ jump=C_g). Use C_ref=g only
    // when GravityForceCK remains in ForcePrior (total-p momentum + p' dens path).
    // French / default path keeps AcousticDensMode=0 (legacy wall dens on Adami mom pressure).
    if (cli.background_hydro_p)
    {
        hostSetPerturbationMomentum(particles, Real(1));
        hostSetAcousticDensMode(particles, Real(1));
        hostSetReferenceBodyAcceleration(particles, Vecd::Zero());
    }
    if (!cli.background_hydro_p)
        constant_gravity.exec();
    if (cli.enable_boussinesq)
    {
        boussinesq_force.exec();
    }

    const OphelieBoussinesqForceAudit force_audit =
        hostAuditBoussinesqForcePrior(particles, beta_eff, cli.t_ref, -cli.gravity_g, kOphelieTemperatureField);
    logOphelieBoussinesqForceAudit(force_audit, "doc77a_init");

    // Stable strat. a>0 → T>Tref above bottom → Fz_bouss = +m β ΔT g > 0 (upward).
    const bool force_sign_ok =
        uniform || !cli.enable_boussinesq ||
        (force_audit.f_z_sum > TinyReal && force_audit.rel_err < Real(1.0e-6));
    const bool force_uniform_ok =
        !uniform || !cli.enable_boussinesq ||
        (std::abs(force_audit.f_z_sum) < Real(1.0e-8) * (std::abs(force_audit.mass_sum) + TinyReal));

    std::ofstream monitor(cli.out_csv);
    monitor << "step,time_s,U_max,U_rms,p_l2_rel,p_linf,rho_mean,phase\n";

    const HydroStats s0 =
        hostHydroStats(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H);
    monitor << "0,0," << s0.u_max << "," << s0.u_rms << "," << s0.p_l2_rel << "," << s0.p_linf << ","
            << s0.rho_mean << ",init\n";

    std::cout << "[doc77a] start mode=" << cli.mode << " dp=" << dp << " n=" << particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " beta=" << beta_eff
              << " a=" << a_eff << " c0=" << cli.c0 << " settle=" << cli.settle_time
              << " end_time=" << cli.end_time << " force_rel_err=" << force_audit.rel_err
              << " Fz_sum=" << force_audit.f_z_sum << std::endl;

    double force_cons_rel = 1.0;
    bool force_update_consistency_pass = false;
    bool viscous_dt_checked = cli.use_viscous_dt;
    size_t n_bad_dt = 0;
    std::vector<Real> dens_sigma0;
    double dens_sigma0_mean = 0.0;
    bool dens_sigma0_ready = false;
    std::ofstream force_ofs;
    if (cli.force_consistency_probe)
    {
        force_ofs.open("output/doc77a_force_consistency.csv");
        force_ofs << "step,sub,dt_ac,du_dt_z_rms,a_force_z_rms,rel_err\n";
    }

    auto capture_sigma0_if_needed = [&]() {
        if (!cli.dens_sum_normalize || dens_sigma0_ready)
            return;
        // DensitySummation must already be current.
        hostCaptureDensSumSigma0(particles, dens_sigma0, dens_sigma0_mean);
        dens_sigma0_ready = true;
        std::cout << "[doc77a] dens_sum_normalize sigma0_mean=" << dens_sigma0_mean
                  << " rho0=" << cli.rho0 << " packing=" << (dens_sigma0_mean / (cli.rho0 + TinyReal))
                  << std::endl;
    };

    auto advance_to = [&](Real t_target, const char *phase) {
        size_t local_iter = 0;
        while (sv_physical_time->getValue() < t_target)
        {
            if (cli.use_density_summation)
            {
                density_summation.exec();
                capture_sigma0_if_needed();
                density_reg.exec();
                if (cli.dens_sum_normalize)
                {
                    // ρ = ρ0 · ρ_sum / σ0  (overwrites Internal reg copy of ρ_sum).
                    hostApplyDensSumNormalize(particles, cli.rho0, dens_sigma0);
                }
                if (cli.reapply_hydro_rho && !cli.background_hydro_p)
                {
                    // Doc87: ρ only — do not rewrite Mass. Skip when R3-A uses p_h field.
                    hostApplyStratifiedHydrostaticDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff,
                                                          a_eff, z_bottom, H, /*sync_mass=*/false,
                                                          /*perturbation_only=*/false);
                }
                // Doc105 P1: dens_sum_normalize → ρ≈ρ0 wipes p'_strat; restore for stratified bgp.
                if (cli.background_hydro_p && cli.enable_boussinesq && std::abs(a_eff) > TinyReal)
                {
                    hostApplyStratifiedHydrostaticDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff,
                                                          a_eff, z_bottom, H, /*sync_mass=*/false,
                                                          /*perturbation_only=*/true);
                }
            }
            if (cli.background_hydro_p)
            {
                hostApplyBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z_bottom, H);
            }
            advection_setup.exec();
            if (cli.kernel_correction)
            {
                fluid_linear_correction_matrix.exec();
            }
            const Real dt_adv = cli.use_viscous_dt ? advection_visc_dt.exec() : advection_dt.exec();
            if (!(dt_adv > TinyReal) || !std::isfinite(static_cast<double>(dt_adv)))
            {
                ++n_bad_dt;
                std::cerr << "[doc77a] non-positive/non-finite dt_adv=" << dt_adv << " — abort" << std::endl;
                break;
            }
            const Real h_min = dp; // lattice spacing proxy for C_nu diagnostic
            const Real nu = cli.mu / (cli.rho0 + TinyReal);
            const Real C_nu = nu * dt_adv / (h_min * h_min + TinyReal);
            viscous_force.exec();
            if (!cli.background_hydro_p)
                constant_gravity.exec();
            if (cli.enable_boussinesq)
            {
                boussinesq_force.exec();
            }

            Real relaxation = 0.0;
            Real dt_ac = 0.0;
            size_t ac_sub = 0;
            const Real t_remain_phase = t_target - sv_physical_time->getValue();
            while (relaxation < dt_adv && ac_sub < 5000)
            {
                dt_ac = acoustic_dt.exec();
                const Real dt_left = std::min(dt_adv - relaxation, t_remain_phase - relaxation);
                if (!(dt_left > TinyReal))
                {
                    break; // finished this advection chunk / phase remainder
                }
                if (dt_ac > dt_left)
                    dt_ac = dt_left;
                if (!(dt_ac > TinyReal) || !std::isfinite(static_cast<double>(dt_ac)))
                {
                    ++n_bad_dt;
                    std::cerr << "[doc77a] bad dt_ac=" << dt_ac << " at relax=" << relaxation
                              << " dt_adv=" << dt_adv << std::endl;
                    break;
                }
                const bool probe_here =
                    cli.force_consistency_probe && local_iter < cli.force_probe_steps && ac_sub == 0;
                std::vector<Real> uz_before;
                if (probe_here)
                {
                    syncVariableToHost<Vecd>(particles, "Velocity");
                    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
                    const size_t n = particles.TotalRealParticles();
                    uz_before.resize(n);
                    for (size_t i = 0; i < n; ++i)
                        uz_before[i] = vel[i][2];
                }
                if (cli.kernel_correction)
                {
                    acoustic_1st_corr.exec(dt_ac);
                }
                else
                {
                    acoustic_1st.exec(dt_ac);
                }
                if (probe_here)
                {
                    double du_rms = 0.0, a_rms = 0.0, rerr = 1.0;
                    hostCompareDuDtVsForce(particles, uz_before, dt_ac, du_rms, a_rms, rerr);
                    force_cons_rel = rerr;
                    force_update_consistency_pass = (rerr < 1.0e-6);
                    force_ofs << local_iter << "," << ac_sub << "," << dt_ac << "," << du_rms << "," << a_rms << ","
                              << rerr << "\n";
                    force_ofs.flush();
                    std::cout << "[doc77a][force-consistency] step=" << local_iter << " du_dt_rms=" << du_rms
                              << " a_force_rms=" << a_rms << " rel_err=" << rerr << std::endl
                              << std::flush;
                }
                if (cli.kernel_correction)
                {
                    acoustic_2nd_corr.exec(dt_ac);
                }
                else
                {
                    acoustic_2nd.exec(dt_ac);
                }
                if (cli.enable_boussinesq)
                {
                    boussinesq_force.exec();
                }
                relaxation += dt_ac;
                sv_physical_time->incrementValue(dt_ac);
                ++ac_sub;
            }
            update_position.exec();
            ++local_iter;
            if (local_iter % 50 == 0)
            {
                particle_sort.exec();
            }
            fluid_cell_linked_list.exec();
            wall_cell_linked_list.exec();
            fluid_update_relation.exec();

            if (local_iter % 10 == 0)
            {
                const Real um = u_max_reduce.exec();
                std::cout << "[doc77a][" << phase << "] N=" << local_iter
                          << " t=" << sv_physical_time->getValue() << " U_max=" << um << " dt_adv=" << dt_adv
                          << " dt_ac=" << dt_ac << " C_nu=" << C_nu
                          << " viscous_dt=" << (cli.use_viscous_dt ? 1 : 0)
                          << " bg_p=" << (cli.background_hydro_p ? 1 : 0)
                          << " kc=" << (cli.kernel_correction ? 1 : 0) << std::endl
                          << std::flush;
                // Abort runaway early (avoid dt-collapse hang).
                if (um > Real(0.5) * cli.c0)
                {
                    std::cerr << "[doc77a] RUNAWAY U_max=" << um << " > 0.5*c0 — abort" << std::endl;
                    ++n_bad_dt;
                    break;
                }
            }
            if (cli.force_consistency_probe && local_iter >= cli.force_probe_steps &&
                phase == std::string("force_probe"))
            {
                break;
            }
            if (n_bad_dt > 0)
                break;
        }
    };

    // Ensure BackgroundPressure exists before first write (registered by AcousticStep ctor).
    if (cli.background_hydro_p)
    {
        hostApplyBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z_bottom, H);
    }
    if (cli.kernel_correction)
    {
        fluid_linear_correction_matrix.exec();
    }

    std::cout << "[doc77a] density_summation=" << (cli.use_density_summation ? 1 : 0)
              << " reapply_hydro_rho=" << (cli.reapply_hydro_rho ? 1 : 0)
              << " dens_sum_normalize=" << (cli.dens_sum_normalize ? 1 : 0)
              << " background_hydro_p=" << (cli.background_hydro_p ? 1 : 0)
              << " acoustic_dens_mode=" << (cli.background_hydro_p ? 1 : 0)
              << " kernel_correction=" << (cli.kernel_correction ? 1 : 0)
              << " viscous_dt=" << (cli.use_viscous_dt ? 1 : 0)
              << " mass_sum_init=" << mass_sum_init
              << " mu=" << cli.mu << " c0=" << cli.c0
              << " (0=continuity/EOS hydrostatic path for Doc77-A)" << std::endl;

    if (cli.density_stage_probe)
    {
        // Doc87 R2/R3-B: freeze particles; dens_sum → dens_reg [→ normalize] → EOS diagnostics.
        fluid_cell_linked_list.exec();
        wall_cell_linked_list.exec();
        fluid_update_relation.exec();
        advection_setup.exec();
        hostDumpDensityStageCsv("output/doc87_r2_dens_init.csv", particles, cli.rho0, cli.c0, z_bottom, H, L, "init");
        density_summation.exec();
        capture_sigma0_if_needed();
        hostDumpDensityStageCsv("output/doc87_r2_dens_after_sum.csv", particles, cli.rho0, cli.c0, z_bottom, H, L,
                                "after_sum");
        density_reg.exec();
        if (cli.dens_sum_normalize)
            hostApplyDensSumNormalize(particles, cli.rho0, dens_sigma0);
        advection_setup.exec(); // Vol = m/rho
        hostDumpDensityStageCsv("output/doc87_r2_dens_after_reg.csv", particles, cli.rho0, cli.c0, z_bottom, H, L,
                                cli.dens_sum_normalize ? "after_norm" : "after_reg");
        // Repeat once after "sort identity" check without moving particles.
        particle_sort.exec();
        fluid_cell_linked_list.exec();
        wall_cell_linked_list.exec();
        fluid_update_relation.exec();
        density_summation.exec();
        density_reg.exec();
        if (cli.dens_sum_normalize)
            hostApplyDensSumNormalize(particles, cli.rho0, dens_sigma0);
        hostDumpDensityStageCsv("output/doc87_r2_dens_after_sort.csv", particles, cli.rho0, cli.c0, z_bottom, H, L,
                                "after_sort_sum_reg");
        // Summarize from after_reg
        syncVariableToHost<Real>(particles, "Density");
        syncVariableToHost<Vecd>(particles, "Position");
        const Real *rho = particles.getVariableDataByName<Real>("Density");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const size_t n = particles.TotalRealParticles();
        double max_adr = 0.0, rms = 0.0, mean = 0.0, max_dp = 0.0;
        double max_adr_inner = 0.0, max_adr_wall = 0.0;
        size_t n_wall = 0, n_inner = 0;
        const Real wall_band = Real(2.0) * dp;
        for (size_t i = 0; i < n; ++i)
        {
            const Real dr = rho[i] - cli.rho0;
            const Real dp_eos = cli.c0 * cli.c0 * dr;
            mean += static_cast<double>(dr);
            rms += static_cast<double>(dr * dr);
            max_adr = std::max(max_adr, std::abs(static_cast<double>(dr)));
            max_dp = std::max(max_dp, std::abs(static_cast<double>(dp_eos)));
            const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
            const Real dist = std::min({x + Real(0.5) * L, Real(0.5) * L - x, y + Real(0.5) * L, Real(0.5) * L - y,
                                        z - z_bottom, H - (z - z_bottom)});
            if (dist < wall_band)
            {
                max_adr_wall = std::max(max_adr_wall, std::abs(static_cast<double>(dr)));
                ++n_wall;
            }
            else
            {
                max_adr_inner = std::max(max_adr_inner, std::abs(static_cast<double>(dr)));
                ++n_inner;
            }
        }
        mean /= static_cast<double>(n);
        rms = std::sqrt(rms / static_cast<double>(n));
        // g=0 rest pass: density error so small that δp cannot drive measurable motion; diagnostic only here.
        const bool zero_g_rest_pass = (cli.gravity_g <= TinyReal) && (max_adr / (cli.rho0 + TinyReal) < Real(1.0e-4));
        if (!cli.summary_json.empty())
        {
            std::ofstream js(cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc77a_density_stage\",\n"
               << "  \"doc87_r2\": 1,\n"
               << "  \"g\": " << cli.gravity_g << ",\n"
               << "  \"rho0\": " << cli.rho0 << ",\n"
               << "  \"c0\": " << cli.c0 << ",\n"
               << "  \"mean_delta_rho\": " << mean << ",\n"
               << "  \"rms_delta_rho\": " << rms << ",\n"
               << "  \"max_abs_delta_rho\": " << max_adr << ",\n"
               << "  \"max_abs_delta_rho_inner\": " << max_adr_inner << ",\n"
               << "  \"max_abs_delta_rho_wall\": " << max_adr_wall << ",\n"
               << "  \"max_abs_delta_p_eos\": " << max_dp << ",\n"
               << "  \"n_inner\": " << n_inner << ",\n"
               << "  \"n_wall_band\": " << n_wall << ",\n"
               << "  \"wall_rho0_solid_default\": 1.0,\n"
               << "  \"dens_sum_normalize\": " << (cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"sigma0_mean\": " << dens_sigma0_mean << ",\n"
               << "  \"note\": \"Contact dens_sum scales by wall rho0; formula cancels if Mass=rho0*Vol\",\n"
               << "  \"zero_g_rest_pass\": " << (zero_g_rest_pass ? 1 : 0) << ",\n"
               << "  \"mass_conservation_pass\": 1,\n"
               << "  \"passed\": " << (zero_g_rest_pass ? 1 : 0) << "\n"
               << "}\n";
        }
        std::cout << "[doc77a] density_stage_probe done mean_dr=" << mean << " rms=" << rms << " max=" << max_adr
                  << " max_inner=" << max_adr_inner << " max_wall=" << max_adr_wall << " max_dp=" << max_dp
                  << " normalize=" << (cli.dens_sum_normalize ? 1 : 0)
                  << " zero_g_rest_pass=" << (zero_g_rest_pass ? 1 : 0) << std::endl;
        return zero_g_rest_pass ? 0 : 2;
    }

    if (cli.hydrostatic_step_probe)
    {
        // Doc97 P1: one real advection + one acoustic sub-step; do not clear state between stages.
        const Real wall_band = Real(2.0) * dp;
        std::ofstream stage_ofs("output/doc97_hydro_step_stages.csv");
        stage_ofs << "stage,n,n_inner,n_wall,uz_rms,uz_max,a_prior_z_rms,a_force_z_rms,a_net_z_rms,"
                     "a_prior_z_mean,a_force_z_mean,a_net_z_mean,drho_dt_rms,drho_dt_max,p_bg_rms,p_eos_rms,"
                     "p_total_rms,p_bg_mean,p_eos_mean,rho_mean,delta_rho_rms,a_net_inner_rms,a_net_wall_rms,"
                     "drho_dt_inner_rms,drho_dt_wall_rms\n";

        fluid_cell_linked_list.exec();
        wall_cell_linked_list.exec();
        fluid_update_relation.exec();
        hostZeroVecField(particles, "Velocity");
        hostZeroVecField(particles, "Force");
        hostZeroVecField(particles, "ForcePrior");
        hostZeroRealField(particles, "DensityChangeRate");

        auto snap = [&](const char *stage) {
            const HydroStepStageStats st =
                hostHydroStepStageStats(particles, cli.rho0, cli.c0, L, z_bottom, H, wall_band);
            hostAppendHydroStepStageCsv(stage_ofs, stage, st);
            std::cout << "[doc77a][hydro-step] " << stage << " uz_rms=" << st.uz_rms
                      << " a_net_rms=" << st.a_net_z_rms << " a_prior_rms=" << st.a_prior_z_rms
                      << " a_force_rms=" << st.a_force_z_rms << " drho_dt_rms=" << st.drho_dt_rms
                      << " p_eos_rms=" << st.p_eos_rms << " p_bg_rms=" << st.p_bg_rms
                      << " a_net_inner=" << st.a_net_inner_rms << " a_net_wall=" << st.a_net_wall_rms
                      << std::endl
                      << std::flush;
            return st;
        };

        HydroStepStageStats s_init = snap("00_init_zero");

        if (cli.use_density_summation)
        {
            density_summation.exec();
            capture_sigma0_if_needed();
            density_reg.exec();
            if (cli.dens_sum_normalize)
                hostApplyDensSumNormalize(particles, cli.rho0, dens_sigma0);
            // Doc105 P1: restore p'_strat after dens_sum_normalize for stratified bgp.
            if (cli.background_hydro_p && cli.enable_boussinesq && std::abs(a_eff) > TinyReal)
                hostApplyStratifiedHydrostaticDensity(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff,
                                                      z_bottom, H, /*sync_mass=*/false,
                                                      /*perturbation_only=*/true);
        }
        if (cli.background_hydro_p)
            hostApplyBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z_bottom, H);
        advection_setup.exec();
        if (cli.kernel_correction)
            fluid_linear_correction_matrix.exec();
        HydroStepStageStats s_dens = snap("01_after_dens_bgp");

        // Skip viscous for hydrostatic reference-state probe (isolates pressure/gravity).
        // ForcePriorCK accumulates (current - previous); zero BOTH ForcePrior and the named
        // current/previous buckets or a second gravity exec cancels to ~0.
        hostZeroVecField(particles, "ForcePrior");
        hostZeroVecField(particles, "GravityForceCK");
        hostZeroVecField(particles, "PreviousGravityForceCK");
        hostZeroVecField(particles, "BoussinesqBuoyancyForceCK");
        hostZeroVecField(particles, "PreviousBoussinesqBuoyancyForceCK");
        if (!cli.background_hydro_p)
            constant_gravity.exec();
        if (cli.enable_boussinesq)
            boussinesq_force.exec();
        HydroStepStageStats s_prior = snap("02_after_force_prior");

        hostZeroVecField(particles, "Force");
        hostZeroRealField(particles, "DensityChangeRate");
        const Real dt_ac = acoustic_dt.exec();
        if (!(dt_ac > TinyReal) || !std::isfinite(static_cast<double>(dt_ac)))
        {
            std::cerr << "[doc77a] hydro-step bad dt_ac=" << dt_ac << std::endl;
            return 2;
        }
        if (cli.kernel_correction)
            acoustic_1st_corr.exec(dt_ac);
        else
            acoustic_1st.exec(dt_ac);
        HydroStepStageStats s_ac1 = snap("03_after_acoustic_1st");

        if (cli.kernel_correction)
            acoustic_2nd_corr.exec(dt_ac);
        else
            acoustic_2nd.exec(dt_ac);
        HydroStepStageStats s_ac2 = snap("04_after_acoustic_2nd");

        // Doc97: background_write_error from p_h field fit (not hydrostatic_update_pass).
        const PressureFit pfit_bg =
            hostFitPressureQuadratic(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H,
                                     /*use_background_field=*/true);
        const bool background_write_pass =
            cli.background_hydro_p && (pfit_bg.c1_rel < Real(0.05)) && (pfit_bg.c1 < Real(0));
        // Doc105: dens path mode 1 = p' with independent C_ref; momentum p' + no GravityForceCK.
        const Real a_ref = std::abs(cli.gravity_g) + TinyReal;
        const Real drho_ref = cli.rho0 * cli.c0 / (H + TinyReal);
        const Real buoyancy_scale =
            std::abs(cli.gravity_g) * std::abs(beta_eff) * std::abs(a_eff) * H; // ~|g β ΔT|
        const bool dens_update_pass = (s_ac1.drho_dt_rms < Real(1.0e-3) * drho_ref);
        const bool inner_momentum_pass =
            (cli.gravity_g <= TinyReal) ? (s_ac1.a_net_inner_rms < Real(1.0e-3))
                                        : (s_ac1.a_net_inner_rms < Real(0.05) * a_ref);
        const bool wall_momentum_pass =
            (cli.gravity_g <= TinyReal) ? (s_ac1.a_net_wall_rms < Real(1.0e-3))
                                        : (s_ac1.a_net_wall_rms < Real(0.05) * a_ref);
        // Doc105 P1: stratified residual must beat buoyancy scale, not only 0.05|g|.
        const bool p1_strat_balance_pass =
            (!cli.background_hydro_p || !cli.enable_boussinesq || std::abs(a_eff) <= TinyReal)
                ? true
                : (s_ac1.a_net_inner_rms < Real(0.25) * (buoyancy_scale + TinyReal));
        const bool hydrostatic_update_pass =
            dens_update_pass && inner_momentum_pass && wall_momentum_pass && p1_strat_balance_pass;
        const bool hydrostatic_update_inner_pass = dens_update_pass && inner_momentum_pass && p1_strat_balance_pass;

        if (!cli.summary_json.empty())
        {
            std::ofstream js(cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc77a_hydrostatic_step\",\n"
               << "  \"doc97_p1\": 1,\n"
               << "  \"mode\": \"" << cli.mode << "\",\n"
               << "  \"g\": " << cli.gravity_g << ",\n"
               << "  \"dp\": " << dp << ",\n"
               << "  \"mu\": " << cli.mu << ",\n"
               << "  \"c0\": " << cli.c0 << ",\n"
               << "  \"dt_ac\": " << dt_ac << ",\n"
               << "  \"dens_sum_normalize\": " << (cli.dens_sum_normalize ? 1 : 0) << ",\n"
               << "  \"background_hydro_p\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
               << "  \"acoustic_dens_mode\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
               << "  \"reference_body_acceleration_z\": 0,\n"
               << "  \"reference_body_acceleration_note\": \"P2b: C_ref=0 because g not in ForcePrior\",\n"
               << "  \"pressure_convention\": \""
               << (cli.background_hydro_p ? "base_ph_plus_pprime_strat" : "legacy_total_p") << "\",\n"
               << "  \"buoyancy_scale\": " << buoyancy_scale << ",\n"
               << "  \"p1_strat_balance_pass\": " << (p1_strat_balance_pass ? 1 : 0) << ",\n"
               << "  \"kernel_correction\": " << (cli.kernel_correction ? 1 : 0) << ",\n"
               << "  \"enable_boussinesq\": " << (cli.enable_boussinesq ? 1 : 0) << ",\n"
               << "  \"a_prior_z_rms\": " << s_prior.a_prior_z_rms << ",\n"
               << "  \"a_force_z_rms_ac1\": " << s_ac1.a_force_z_rms << ",\n"
               << "  \"a_net_z_rms_ac1\": " << s_ac1.a_net_z_rms << ",\n"
               << "  \"a_net_inner_rms_ac1\": " << s_ac1.a_net_inner_rms << ",\n"
               << "  \"a_net_wall_rms_ac1\": " << s_ac1.a_net_wall_rms << ",\n"
               << "  \"uz_rms_ac1\": " << s_ac1.uz_rms << ",\n"
               << "  \"drho_dt_rms_ac1\": " << s_ac1.drho_dt_rms << ",\n"
               << "  \"drho_dt_inner_rms_ac1\": " << s_ac1.drho_dt_inner_rms << ",\n"
               << "  \"drho_dt_wall_rms_ac1\": " << s_ac1.drho_dt_wall_rms << ",\n"
               << "  \"p_bg_rms\": " << s_dens.p_bg_rms << ",\n"
               << "  \"p_eos_rms_ac1\": " << s_ac1.p_eos_rms << ",\n"
               << "  \"p_total_rms_ac1\": " << s_ac1.p_total_rms << ",\n"
               << "  \"delta_rho_rms_ac2\": " << s_ac2.delta_rho_rms << ",\n"
               << "  \"uz_rms_ac2\": " << s_ac2.uz_rms << ",\n"
               << "  \"background_write_error\": " << pfit_bg.c1_rel << ",\n"
               << "  \"background_write_pass\": " << (background_write_pass ? 1 : 0) << ",\n"
               << "  \"pfit_c1_rel_note\": \"legacy alias of background field fit; not update accuracy\",\n"
               << "  \"dens_update_pass\": " << (dens_update_pass ? 1 : 0) << ",\n"
               << "  \"inner_momentum_pass\": " << (inner_momentum_pass ? 1 : 0) << ",\n"
               << "  \"wall_momentum_pass\": " << (wall_momentum_pass ? 1 : 0) << ",\n"
               << "  \"hydrostatic_update_inner_pass\": " << (hydrostatic_update_inner_pass ? 1 : 0) << ",\n"
               << "  \"hydrostatic_update_pass\": " << (hydrostatic_update_pass ? 1 : 0) << ",\n"
               << "  \"hydrostatic_update_pass_note\": \"full pass needs dens+inner+wall; P2 dens uses p'; wall mom OPEN\",\n"
               << "  \"p2_prime_dissipation\": 1,\n"
               << "  \"perturbation_momentum\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
               << "  \"gravity_force_ck\": " << (cli.background_hydro_p ? 0 : 1) << ",\n"
               << "  \"stage_csv\": \"output/doc97_hydro_step_stages.csv\",\n"
               << "  \"passed\": " << (hydrostatic_update_pass ? 1 : 0) << "\n"
               << "}\n";
        }
        std::cout << "[doc77a] hydrostatic_step_probe done dt_ac=" << dt_ac
                  << " a_net_ac1=" << s_ac1.a_net_z_rms << " drho_dt_ac1=" << s_ac1.drho_dt_rms
                  << " dens_pass=" << (dens_update_pass ? 1 : 0)
                  << " inner_mom_pass=" << (inner_momentum_pass ? 1 : 0)
                  << " wall_mom_pass=" << (wall_momentum_pass ? 1 : 0)
                  << " background_write_pass=" << (background_write_pass ? 1 : 0)
                  << " hydrostatic_update_pass=" << (hydrostatic_update_pass ? 1 : 0)
                  << " pprime_mom=" << (cli.background_hydro_p ? 1 : 0) << std::endl;
        (void)s_init;
        (void)hydrostatic_update_inner_pass;
        return hydrostatic_update_pass ? 0 : 2;
    }

    if (cli.force_consistency_probe)
    {
        advance_to(cli.settle_time + cli.end_time, "force_probe");
        const double mass_end_probe = hostTotalMass(particles);
        const bool mass_ok =
            std::abs(mass_end_probe - mass_sum_init) <= 1e-12 * (std::abs(mass_sum_init) + 1.0);
        if (!cli.summary_json.empty())
        {
            std::ofstream js(cli.summary_json);
            js << "{\n"
               << "  \"case\": \"doc77a_force_consistency\",\n"
               << "  \"doc87_r1b\": 1,\n"
               << "  \"force_cons_rel\": " << force_cons_rel << ",\n"
               << "  \"force_update_consistency_pass\": " << (force_update_consistency_pass ? 1 : 0) << ",\n"
               << "  \"mass_sum_init\": " << mass_sum_init << ",\n"
               << "  \"mass_sum_end\": " << mass_end_probe << ",\n"
               << "  \"mass_conservation_pass\": " << (mass_ok ? 1 : 0) << ",\n"
               << "  \"viscous_dt_checked\": " << (viscous_dt_checked ? 1 : 0) << ",\n"
               << "  \"n_bad_dt\": " << n_bad_dt << ",\n"
               << "  \"passed\": " << (force_update_consistency_pass && mass_ok && n_bad_dt == 0 ? 1 : 0) << "\n"
               << "}\n";
        }
        std::cout << "[doc77a] force_probe done rel=" << force_cons_rel
                  << " pass=" << (force_update_consistency_pass ? 1 : 0) << " mass_ok=" << (mass_ok ? 1 : 0)
                  << std::endl;
        return (force_update_consistency_pass && mass_ok && n_bad_dt == 0) ? 0 : 2;
    }

    // Phase 1: settle under gravity (+ buoyancy). Do NOT zero velocity — rebuild noise after
    // artificial zero falsely fails "no growth". Observation continues from settle state.
    advance_to(cli.settle_time, "settle");
    HydroStats s_settle =
        hostHydroStats(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H);
    monitor << "settle," << sv_physical_time->getValue() << "," << s_settle.u_max << "," << s_settle.u_rms << ","
            << s_settle.p_l2_rel << "," << s_settle.p_linf << "," << s_settle.rho_mean << ",settle\n";

    const Real t_obs1 = sv_physical_time->getValue() + cli.end_time;
    advance_to(t_obs1, "observe");
    HydroStats s_end =
        hostHydroStats(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H);
    monitor << "end," << sv_physical_time->getValue() << "," << s_end.u_max << "," << s_end.u_rms << ","
            << s_end.p_l2_rel << "," << s_end.p_linf << "," << s_end.rho_mean << ",observe\n";

    hostWritePressureProfileCsv("output/doc77a_pressure_profile.csv", particles, cli.rho0, cli.c0, cli.gravity_g,
                                beta_eff, a_eff, z_bottom, H);
    if (cli.background_hydro_p)
    {
        hostApplyBackgroundHydroPressure(particles, cli.rho0, cli.gravity_g, z_bottom, H);
    }
    const PressureFit pfit =
        hostFitPressureQuadratic(particles, cli.rho0, cli.c0, cli.gravity_g, beta_eff, a_eff, z_bottom, H,
                                 cli.background_hydro_p);

    if (!cli.background_hydro_p)
        constant_gravity.exec();
    if (cli.enable_boussinesq)
    {
        boussinesq_force.exec();
    }
    hostWriteAccelProfileCsv("output/doc77a_accel_profile.csv", particles, cli.gravity_g, beta_eff, a_eff, cli.t_ref,
                             z_bottom, H);
    // Doc87 R1-D: refresh buoyancy from final T before end audit.
    if (cli.enable_boussinesq)
    {
        boussinesq_force.exec();
    }
    const OphelieBoussinesqForceAudit force_end =
        hostAuditBoussinesqForcePrior(particles, beta_eff, cli.t_ref, -cli.gravity_g, kOphelieTemperatureField);
    logOphelieBoussinesqForceAudit(force_end, "doc77a_end");
    const double mass_sum_end = hostTotalMass(particles);
    const bool mass_conservation_pass =
        std::abs(mass_sum_end - mass_sum_init) <= 1e-12 * (std::abs(mass_sum_init) + 1.0);

    // Pass gates: force + residual floor + no growth + no runaway.
    // With --reapply-hydro-rho, EOS profile is restored after dens_sum → pfit is meaningful.
    const bool eos_profile_live =
        cli.background_hydro_p || (!cli.use_density_summation) || cli.reapply_hydro_rho;
    const Real u_hydro = std::sqrt(cli.gravity_g * H + TinyReal);
    const Real u_rms_tol = (cli.gravity_g > TinyReal) ? Real(0.25) * u_hydro : Real(1.0e-3);
    const Real u_max_tol = (cli.gravity_g > TinyReal) ? Real(0.45) * u_hydro : Real(1.0e-3);
    const bool no_runaway = s_end.u_max < Real(0.25) * cli.c0;
    const bool residual_ok = s_end.u_rms < u_rms_tol && s_end.u_max < u_max_tol;
    // g=0: absolute floor (relative growth on ~1e-7 noise is meaningless).
    const bool no_growth =
        (cli.gravity_g <= TinyReal)
            ? (s_end.u_rms < Real(1.0e-4) && s_end.u_max < Real(1.0e-3))
            : (s_end.u_rms <= Real(1.25) * s_settle.u_rms + TinyReal);
    const bool c1_diag = eos_profile_live && (pfit.c1_rel < Real(0.35) && pfit.c1 < Real(0));
    // Doc105 P1: with background_hydro_p, p_h is base-only (c2_an=0); skip strat c2 on p_h fit.
    const bool c2_diag =
        eos_profile_live && (uniform || !cli.enable_boussinesq || cli.background_hydro_p ||
                             (pfit.c2 > Real(0) && pfit.c2_rel < Real(0.75)));
    const bool force_end_ok =
        force_end.rel_err < Real(1.0e-6) &&
        (uniform || !cli.enable_boussinesq || force_end.f_z_sum > TinyReal);
    // Doc84 hydrostatic_pass: tight residual vs weak French signal scale + EOS slope when live.
    const bool hydrostatic_pass =
        force_end_ok && no_runaway && no_growth && (s_end.u_rms < Real(0.05) * u_hydro + Real(1.0e-4)) &&
        (!eos_profile_live || (c1_diag && (uniform || !cli.enable_boussinesq || c2_diag)));
    const bool passed =
        force_sign_ok && force_uniform_ok && force_end_ok && no_runaway && residual_ok && no_growth &&
        mass_conservation_pass && n_bad_dt == 0;

    std::cout << "[doc77a] end t=" << sv_physical_time->getValue() << " U_max=" << s_end.u_max
              << " U_rms=" << s_end.u_rms << " U_settle=" << s_settle.u_rms << " p_l2_rel=" << s_end.p_l2_rel
              << " pfit_c1=" << pfit.c1 << " c1_an=" << pfit.c1_an << " c1_rel=" << pfit.c1_rel
              << " pfit_c2=" << pfit.c2 << " c2_an=" << pfit.c2_an << " c2_rel=" << pfit.c2_rel
              << " pfit_diag_ok=" << (c1_diag && c2_diag ? 1 : 0) << " u_rms_tol=" << u_rms_tol
              << " u_max_tol=" << u_max_tol
              << " force_ok=" << (force_sign_ok && force_uniform_ok && force_end_ok ? 1 : 0)
              << " residual_ok=" << (residual_ok ? 1 : 0) << " no_growth=" << (no_growth ? 1 : 0)
              << " no_runaway=" << (no_runaway ? 1 : 0)
              << " mass_ok=" << (mass_conservation_pass ? 1 : 0)
              << " hydrostatic_pass=" << (hydrostatic_pass ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    std::cout << "[doc77a] reuse_path=GravityForceCK+BoussinesqBuoyancyForceCK+AcousticRiemann+ViscousForceWithWall"
                 "+(optional DensitySummation)"
              << " density_summation=" << (cli.use_density_summation ? 1 : 0)
              << " reapply_hydro_rho=" << (cli.reapply_hydro_rho ? 1 : 0) << std::endl;

    if (!cli.summary_json.empty())
    {
        std::ofstream js(cli.summary_json);
        js << "{\n"
           << "  \"case\": \"doc77a\",\n"
           << "  \"mode\": \"" << cli.mode << "\",\n"
           << "  \"dp\": " << dp << ",\n"
           << "  \"density_summation\": " << (cli.use_density_summation ? 1 : 0) << ",\n"
           << "  \"reapply_hydro_rho\": " << (cli.reapply_hydro_rho ? 1 : 0) << ",\n"
           << "  \"dens_sum_normalize\": " << (cli.dens_sum_normalize ? 1 : 0) << ",\n"
           << "  \"background_hydro_p\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
           << "  \"acoustic_dens_mode\": " << (cli.background_hydro_p ? 1 : 0) << ",\n"
           << "  \"kernel_correction\": " << (cli.kernel_correction ? 1 : 0) << ",\n"
           << "  \"sigma0_mean\": " << dens_sigma0_mean << ",\n"
           << "  \"g\": " << cli.gravity_g << ",\n"
           << "  \"mu\": " << cli.mu << ",\n"
           << "  \"c0\": " << cli.c0 << ",\n"
           << "  \"zero_g_rest_pass\": "
           << ((cli.gravity_g <= TinyReal)
                   ? ((s_end.u_rms < Real(1.0e-4) && s_end.u_max < Real(1.0e-3)) ? 1 : 0)
                   : -1)
           << ",\n"
           << "  \"force_ok\": " << (force_sign_ok && force_uniform_ok && force_end_ok ? 1 : 0) << ",\n"
           << "  \"force_rel_err\": " << force_end.rel_err << ",\n"
           << "  \"Fz_sum\": " << force_end.f_z_sum << ",\n"
           << "  \"U_max\": " << s_end.u_max << ",\n"
           << "  \"U_rms\": " << s_end.u_rms << ",\n"
           << "  \"U_rms_settle\": " << s_settle.u_rms << ",\n"
           << "  \"residual_ok\": " << (residual_ok ? 1 : 0) << ",\n"
           << "  \"no_growth\": " << (no_growth ? 1 : 0) << ",\n"
           << "  \"no_runaway\": " << (no_runaway ? 1 : 0) << ",\n"
           << "  \"p_l2_rel\": " << s_end.p_l2_rel << ",\n"
           << "  \"pfit_c1_rel\": " << pfit.c1_rel << ",\n"
           << "  \"pfit_c2_rel\": " << pfit.c2_rel << ",\n"
           << "  \"background_write_error\": " << (cli.background_hydro_p ? pfit.c1_rel : -1) << ",\n"
           << "  \"background_write_pass\": "
           << (cli.background_hydro_p ? ((pfit.c1_rel < Real(0.05) && pfit.c1 < Real(0)) ? 1 : 0) : -1) << ",\n"
           << "  \"hydrostatic_pass\": " << (hydrostatic_pass ? 1 : 0) << ",\n"
           << "  \"hydrostatic_pass_note\": \"not hydrostatic_update_pass; see Doc97 P1 step probe\",\n"
           << "  \"mass_sum_init\": " << mass_sum_init << ",\n"
           << "  \"mass_sum_end\": " << mass_sum_end << ",\n"
           << "  \"mass_conservation_pass\": " << (mass_conservation_pass ? 1 : 0) << ",\n"
           << "  \"viscous_dt_checked\": " << (viscous_dt_checked ? 1 : 0) << ",\n"
           << "  \"force_update_consistency_pass\": "
           << (cli.force_consistency_probe ? (force_update_consistency_pass ? 1 : 0) : -1) << ",\n"
           << "  \"n_bad_dt\": " << n_bad_dt << ",\n"
           << "  \"smoke_pass\": " << (passed ? 1 : 0) << ",\n"
           << "  \"passed\": " << (passed ? 1 : 0) << "\n"
           << "}\n";
        std::cout << "[doc77a] summary_json=" << cli.summary_json << std::endl;
    }

    return passed ? 0 : 2;
}
