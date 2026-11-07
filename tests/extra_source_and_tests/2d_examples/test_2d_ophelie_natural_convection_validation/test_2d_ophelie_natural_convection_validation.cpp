/**
 * Doc203 C1 smoke: 2D differentially heated cavity, central WCSPH (NoRiemann) + Boussinesq.
 * Coarse L/dp≈40; --tstar-end for nondimensional end time (default 0.05 smoke / 0.5 formal).
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_laplace.h"
#include "ophelie_central_riemann_ck.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
constexpr const char *kT = "Temperature";
constexpr const char *kDT = "OphelieThermalDeltaT";
constexpr const char *kLap = "OphelieThermalLaplaceT";
constexpr const char *kCond = "OphelieThermalConductivity";

template <typename DataType>
inline void syncH(BaseParticles &p, const std::string &n)
{
#if SPHINXSYS_USE_SYCL
    p.template getVariableByName<DataType>(n)->prepareForOutput(execution::par_device);
#else
    (void)p;
    (void)n;
#endif
}
template <typename DataType>
inline void syncD(BaseParticles &p, const std::string &n)
{
#if SPHINXSYS_USE_SYCL
    p.template getVariableByName<DataType>(n)->finalizeLoadIn(execution::par_device);
#else
    (void)p;
    (void)n;
#endif
}

/** T += dt * (1/(ρ cp)) * (−L) with L stored as pairwise Laplace of T (includes k). */
class ApplyThermalDtCK : public LocalDynamics
{
  public:
    ApplyThermalDtCK(SPHBody &body, Real rho, Real cp, Real t0)
        : LocalDynamics(body), factor_(-Real(1) / (rho * cp + TinyReal)), t0_(t0),
          dv_lap_(particles_->getVariableByName<Real>(kLap)),
          dv_T_(particles_->getVariableByName<Real>(kT)),
          dv_dT_(particles_->getVariableByName<Real>(kDT))
    {
    }
    class UpdateKernel
    {
      public:
        template <class EP, class Enc>
        UpdateKernel(const EP &ep, Enc &e)
            : f_(e.factor_), t0_(e.t0_), lap_(e.dv_lap_->DelegatedData(ep)), T_(e.dv_T_->DelegatedData(ep)),
              dT_(e.dv_dT_->DelegatedData(ep))
        {
        }
        void update(size_t i, Real dt)
        {
            dT_[i] += dt * f_ * lap_[i];
            T_[i] = t0_ + dT_[i];
        }
        Real f_, t0_;
        Real *lap_, *T_, *dT_;
    };
    Real factor_, t0_;
    DiscreteVariable<Real> *dv_lap_, *dv_T_, *dv_dT_;
};

void zeroVec(BaseParticles &p, const std::string &name)
{
    syncH<Vecd>(p, name);
    Vecd *v = p.getVariableDataByName<Vecd>(name);
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
        v[i] = Vecd::Zero();
    syncD<Vecd>(p, name);
}

void setWallThermal(BaseParticles &wall, Real L, Real Th, Real Tc, Real T0, Real k_fluid, Real k_wall_factor)
{
    syncH<Vecd>(wall, "Position");
    syncH<Real>(wall, kT);
    syncH<Real>(wall, kDT);
    syncH<Real>(wall, kCond);
    Vecd *pos = wall.getVariableDataByName<Vecd>("Position");
    Real *T = wall.getVariableDataByName<Real>(kT);
    Real *dT = wall.getVariableDataByName<Real>(kDT);
    Real *cond = wall.getVariableDataByName<Real>(kCond);
    const Real eps = Real(1.0e-6) * L;
    for (size_t i = 0; i < wall.TotalRealParticles(); ++i)
    {
        const Real x = pos[i][0], y = pos[i][1];
        Real Tbc = T0;
        Real k = Real(0); // default adiabatic
        if (x < eps)
        {
            Tbc = Th;
            k = k_wall_factor * k_fluid;
        }
        else if (x > L - eps)
        {
            Tbc = Tc;
            k = k_wall_factor * k_fluid;
        }
        // top/bottom remain k=0
        (void)y;
        T[i] = Tbc;
        dT[i] = Tbc - T0;
        cond[i] = k;
    }
    syncD<Real>(wall, kT);
    syncD<Real>(wall, kDT);
    syncD<Real>(wall, kCond);
}

/** Doc209 §10.1: Eulerian fixed linear T(x)=Th+(Tc-Th)x/L; remap from Position (no thermal PDE). */
void applyFixedLinearTemperature(BaseParticles &fluid, Real L, Real Th, Real Tc, Real T0)
{
    syncH<Vecd>(fluid, "Position");
    syncH<Real>(fluid, kT);
    syncH<Real>(fluid, kDT);
    Vecd *pos = fluid.getVariableDataByName<Vecd>("Position");
    Real *T = fluid.getVariableDataByName<Real>(kT);
    Real *dTv = fluid.getVariableDataByName<Real>(kDT);
    for (size_t i = 0; i < fluid.TotalRealParticles(); ++i)
    {
        Real xi = pos[i][0] / (L + TinyReal);
        if (xi < Real(0))
            xi = Real(0);
        if (xi > Real(1))
            xi = Real(1);
        const Real Tbc = Th + (Tc - Th) * xi;
        T[i] = Tbc;
        dTv[i] = Tbc - T0;
    }
    syncD<Real>(fluid, kT);
    syncD<Real>(fluid, kDT);
}

/** Doc209 §10.3: host-side kinetic energy Σ ½ m |v|² (sync each call — use sparsely). */
inline double computeParticleKE(BaseParticles &p)
{
    syncH<Vecd>(p, "Velocity");
    syncH<Real>(p, "Mass");
    Vecd *vel = p.getVariableDataByName<Vecd>("Velocity");
    Real *mass = p.getVariableDataByName<Real>("Mass");
    double ke = 0.0;
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
        ke += 0.5 * double(mass[i]) * double(vel[i].squaredNorm());
    return ke;
}

/** Doc209 §5.2: IDW sample of signed field components on a geometric centerline. */
struct CenterlineSample
{
    Real s_over_L = 0;
    Real u_star = 0;
    Real v_star = 0;
    Real theta = 0;
    Real nearest_dist_over_L = 0;
    int n_support = 0;
};

inline CenterlineSample sampleCenterlineIDW(const Vecd *pos, const Vecd *vel, const Real *T, size_t n,
                                             Real L, Real Tc, Real dT, Real U_char, Real qx, Real qy,
                                             Real support_radius)
{
    CenterlineSample out;
    out.s_over_L = (std::abs(qx - Real(0.5) * L) < std::abs(qy - Real(0.5) * L)) ? (qy / L) : (qx / L);
    double wu = 0, wv = 0, wt = 0, wsum = 0;
    Real dmin = L;
    int ns = 0;
    const Real R2 = support_radius * support_radius;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - qx;
        const Real dy = pos[i][1] - qy;
        const Real d2 = dx * dx + dy * dy;
        if (d2 > R2)
            continue;
        const Real d = std::sqrt(d2 + TinyReal);
        dmin = std::min(dmin, d);
        const double w = 1.0 / (double(d2) + 1e-18);
        wu += w * double(vel[i][0]);
        wv += w * double(vel[i][1]);
        wt += w * double(T[i]);
        wsum += w;
        ++ns;
    }
    out.n_support = ns;
    out.nearest_dist_over_L = dmin / L;
    if (wsum > 0.0)
    {
        out.u_star = Real(wu / wsum) / U_char;
        out.v_star = Real(wv / wsum) / U_char;
        out.theta = (Real(wt / wsum) - Tc) / (dT + TinyReal);
    }
    return out;
}

} // namespace

int main(int argc, char **argv)
{
    // Doc203 §12.1 defaults (C1 air); --c2 retargets Pr/Ra via nu,alpha (§13).
    const Real L = Real(0.0386);
    const Real Th = Real(303), Tc = Real(283), T0 = Real(293), Tref = Real(293);
    const Real rho0 = Real(1.204);
    const Real cp = Real(1006);
    const Real beta = Real(3.43e-3);
    const Real g = Real(9.81);
    const Real dT = Th - Tc;
    Real nu = Real(1.506e-5);
    Real k = Real(0.02587);

    Real n_div = Real(40); // coarse smoke
    Real tstar_end = Real(0.05);
    // riemann | noriemann | central (AverageP+dens UJump) | dissipative (NoLimiter AV)
    std::string scheme = "riemann";
    std::string case_id = "C1";
    // Doc77c default: DensitySummation + Internal reg + packing normalize.
    bool enable_dens_sum = true;
    bool enable_dens_reg = true;
    bool enable_dens_norm = true;
    bool dens_reg_cli = false;
    bool enable_kernel_corr = false;
    int dens_mode_cli = -1; // -1 → scheme default
    bool thermal_1d = false; // Doc204 §4.1: freeze flow, pure conduction Nu audit
    bool fixed_T_linear = false; // Doc209 §10.1: freeze Eulerian linear T; hydro+Boussinesq only
    bool ke_budget = false;      // Doc209 §10.3: sparse acoustic half-step ΔKE
    std::string out_json = "c1_smoke_report.json";
    std::string centerline_csv = "c1_centerline_profiles.csv";
    std::string ke_budget_csv = "c1_ke_budget.csv";
    Real centerline_support_factor = Real(2.0); // ×dp; Doc209 sensitivity also uses 1.5
    std::vector<char *> sph_argv;
    sph_argv.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--n" && i + 1 < argc)
            n_div = Real(std::atof(argv[++i]));
        else if (a == "--tstar-end" && i + 1 < argc)
            tstar_end = Real(std::atof(argv[++i]));
        else if (a == "--scheme" && i + 1 < argc)
            scheme = argv[++i];
        else if (a == "--out" && i + 1 < argc)
            out_json = argv[++i];
        else if (a == "--centerline-csv" && i + 1 < argc)
            centerline_csv = argv[++i];
        else if (a == "--centerline-support-dp" && i + 1 < argc)
            centerline_support_factor = Real(std::atof(argv[++i]));
        else if (a == "--thermal-1d")
        {
            // Same cavity walls + OpheliePairwiseLaplace contact as C1; flow frozen.
            // Analytic steady: T=Th+(Tc-Th)x/L, Nu=1 both sides.
            thermal_1d = true;
            case_id = "C1_thermal_1d";
        }
        else if (a == "--fixed-T-linear")
        {
            // Doc209 §10.1: no thermal PDE; remap T from Eulerian linear field each advection.
            // Diagnosis only — not a DVD Table V claim (linear T ≠ cavity convection T).
            fixed_T_linear = true;
        }
        else if (a == "--ke-budget")
            ke_budget = true;
        else if (a == "--ke-budget-csv" && i + 1 < argc)
            ke_budget_csv = argv[++i];
        else if (a == "--c2")
        {
            // Doc203 §13: Pr≈1e3, Ra≈1e4; keep L,rho,cp,beta,g,dT; solve nu,alpha,k.
            case_id = "C2";
            const Real Pr_t = Real(1.0e3);
            const Real Ra_t = Real(1.0e4);
            const Real nu_alpha = g * beta * dT * L * L * L / Ra_t;
            nu = std::sqrt(Pr_t * nu_alpha);
            const Real alpha_t = nu / Pr_t;
            k = alpha_t * rho0 * cp;
        }
        else if (a == "--dvd-exact")
        {
            // Doc204 §3.2: dvd_exact_ra1e5_pr071 — keep L,rho,cp,beta,g,dT; retarget nu,k.
            case_id = "C1_dvd_exact_ra1e5_pr071";
            const Real Pr_t = Real(0.71);
            const Real Ra_t = Real(1.0e5);
            const Real B = g * beta * dT * L * L * L;
            nu = std::sqrt(Pr_t * B / Ra_t);
            const Real alpha_t = nu / Pr_t;
            k = alpha_t * rho0 * cp;
        }
        else if (a == "--dens-reg")
        {
            enable_dens_reg = true;
            dens_reg_cli = true;
        }
        else if (a == "--no-dens-reg")
        {
            enable_dens_reg = false;
            dens_reg_cli = true;
        }
        else if (a == "--no-dens-sum")
            enable_dens_sum = false;
        else if (a == "--dens-norm")
            enable_dens_norm = true;
        else if (a == "--no-dens-norm")
            enable_dens_norm = false;
        else if (a == "--kernel-corr")
            enable_kernel_corr = true;
        else if (a == "--dens-mode" && i + 1 < argc)
            dens_mode_cli = std::atoi(argv[++i]);
        else
            sph_argv.push_back(argv[i]);
    }
    if (thermal_1d)
        case_id = "C1_thermal_1d";
    if (thermal_1d && fixed_T_linear)
    {
        std::cerr << "Doc209: refuse combining --thermal-1d and --fixed-T-linear\n";
        return 1;
    }
    if (fixed_T_linear)
        case_id = case_id + "_fixed_T_linear";
    const Real alpha = k / (rho0 * cp);
    const Real mu = rho0 * nu;
    const Real Ra = g * beta * dT * L * L * L / (nu * alpha);
    const Real Pr = nu / alpha;
    const Real U_char = alpha / L;
    // WCSPH: Ma≪1 vs estimated circulation; c0² ≫ gL.
    const Real U_est = (case_id == "C2") ? Real(20) * U_char : Real(70) * U_char;
    const Real c0 = std::max({Real(10) * U_est, Real(20) * std::sqrt(g * L),
                              Real(10) * std::sqrt(g * beta * dT * L)});
    const bool use_noriemann = (scheme == "noriemann");
    const bool use_central = (scheme == "central");
    const bool use_dissipative = (scheme == "dissipative");
    // Pure NoRiemann default: dens_sum+norm, no Internal reg (unless --dens-reg).
    if (use_noriemann && enable_dens_sum && !dens_reg_cli)
        enable_dens_reg = false;
    const Real dens_mode =
        dens_mode_cli >= 0 ? Real(dens_mode_cli) : (use_noriemann ? Real(3) : Real(1));
    const Real dp = L / n_div;
    const Real BW = Real(3) * dp;
    const Real end_time = tstar_end * L * L / alpha;

    BoundingBoxd bounds(Vec2d(-BW, -BW), Vec2d(L + BW, L + BW));
    SPHSystem sph_system(bounds, dp);
    int sph_argc = static_cast<int>(sph_argv.size());
    sph_system.handleCommandlineOptions(sph_argc, sph_argv.data());

    GeometricShapeBox fluid_shape(Transform(Vec2d(0.5 * L, 0.5 * L)), Vec2d(0.5 * L, 0.5 * L), "Fluid");
    FluidBody fluid(sph_system, fluid_shape);
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, c0);
    fluid.addMaterialProperty<Viscosity>(mu);
    fluid.generateParticles<BaseParticles, Lattice>();

    class WallBoundary : public ComplexShape
    {
      public:
        WallBoundary(const std::string &name, Real Len, Real wall_bw) : ComplexShape(name)
        {
            add<GeometricShapeBox>(Transform(Vec2d(0.5 * Len, 0.5 * Len)),
                                   Vec2d(0.5 * Len + wall_bw, 0.5 * Len + wall_bw));
            subtract<GeometricShapeBox>(Transform(Vec2d(0.5 * Len, 0.5 * Len)), Vec2d(0.5 * Len, 0.5 * Len));
        }
    };
    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall", L, BW));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall(fluid, {&wall});
    BaseParticles &fp = fluid.getBaseParticles();
    BaseParticles &wp = wall.getBaseParticles();

    auto regThermal = [](BaseParticles &p, Real Tinit, Real k0)
    {
        p.registerStateVariable<Real>(kT, Tinit);
        p.registerStateVariable<Real>(kDT, Real(0));
        p.registerStateVariable<Real>(kLap, Real(0));
        p.registerStateVariable<Real>(kCond, k0);
        p.addEvolvingVariable<Real>(kT);
        p.addEvolvingVariable<Real>(kDT);
    };
    regThermal(fp, T0, k);
    regThermal(wp, T0, k);
    if (thermal_1d || fixed_T_linear)
    {
        // thermal-1d: start near analytic steady. fixed-T: impose Eulerian linear drive from t=0.
        applyFixedLinearTemperature(fp, L, Th, Tc, T0);
    }
    fp.registerStateVariable<Real>("DensSumSigma0");
    fp.addEvolvingVariable<Real>("DensSumSigma0");

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cells(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cells(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> rel(fluid_inner, fluid_wall);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_pos(fluid);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        lin_corr(DynamicsArgs(fluid_inner, Real(0.5)), fluid_wall);

    // Default: NoKernelCorrection (doc77c default). --kernel-corr enables LinearCorrection.
    using A1NR = fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, NoRiemannSolverCK, NoKernelCorrectionCK>,
                                                     Contact<Wall, NoRiemannSolverCK, NoKernelCorrectionCK>>;
    using A2NR = fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK;
    using A1R = fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK;
    using A2R = fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK;
    using A1NRc = fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, NoRiemannSolverCK, LinearCorrectionCK>,
                                                      Contact<Wall, NoRiemannSolverCK, LinearCorrectionCK>>;
    using A2NRc = fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, NoRiemannSolverCK, LinearCorrectionCK>,
                                                      Contact<Wall, NoRiemannSolverCK, LinearCorrectionCK>>;
    using A1Rc = fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK;
    using A2Rc = fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK;
    using A1C = fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, OphelieCentralDensJumpCK, NoKernelCorrectionCK>,
                                                    Contact<Wall, OphelieCentralDensJumpCK, NoKernelCorrectionCK>>;
    using A2C = fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, OphelieCentralDensJumpCK, NoKernelCorrectionCK>,
                                                    Contact<Wall, OphelieCentralDensJumpCK, NoKernelCorrectionCK>>;
    using A1D = fluid_dynamics::AcousticStep1stHalf<
        Inner<OneLevel, DissipativeRiemannSolverCK, NoKernelCorrectionCK>,
        Contact<Wall, DissipativeRiemannSolverCK, NoKernelCorrectionCK>>;
    using A2D = fluid_dynamics::AcousticStep2ndHalf<
        Inner<OneLevel, DissipativeRiemannSolverCK, NoKernelCorrectionCK>,
        Contact<Wall, DissipativeRiemannSolverCK, NoKernelCorrectionCK>>;

    InteractionDynamicsCK<MainExecutionPolicy, A1NR> acoustic_1st_nr(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2NR> acoustic_2nd_nr(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A1R> acoustic_1st_r(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2R> acoustic_2nd_r(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A1NRc> acoustic_1st_nr_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2NRc> acoustic_2nd_nr_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A1Rc> acoustic_1st_r_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2Rc> acoustic_2nd_r_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A1C> acoustic_1st_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2C> acoustic_2nd_c(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A1D> acoustic_1st_d(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, A2D> acoustic_2nd_d(fluid_inner, fluid_wall);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK> viscous(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>, Contact<>>>
        density_sum(fluid_inner, fluid_wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, Internal>> density_reg(fluid);
    ParticleSortCK<MainExecutionPolicy> particle_sort(fluid);

    Gravity gravity(Vecd(0.0, -g)); // g-vector for Boussinesq only (no GravityForceCK on P1)
    StateDynamics<MainExecutionPolicy, BoussinesqBuoyancyForceCK> boussinesq(fluid, gravity, beta, Tref, kT);

    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> lap_inner(
        fluid_inner, kT, kCond, kLap, Real(0.01));
    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Contact<>>> lap_wall(
        fluid_wall, kT, kCond, kLap, Real(0.01), /*dirichlet_ghost=*/true);
    StateDynamics<MainExecutionPolicy, ApplyThermalDtCK> apply_th(fluid, rho0, cp, T0);

    const Real U_dt = Real(100) * U_char;
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> adv_dt(fluid, U_dt);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionViscousTimeStepCK> adv_visc_dt(fluid, U_dt);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> ac_dt(fluid);

    wall_normal.exec();
    fluid_cells.exec();
    wall_cells.exec();
    rel.exec();
    setWallThermal(wp, L, Th, Tc, T0, k, Real(1.0e6));

    // Doc105 P1: p_h absorbs ρ0 g; BoussinesqBuoyancyForceCK only; p'-momentum.
    // Construct acoustics first so BackgroundPressure exists.
    syncH<Real>(fp, "PerturbationMomentum");
    syncH<Real>(fp, "AcousticDensMode");
    syncH<Vecd>(fp, "ReferenceBodyAcceleration");
    syncH<Real>(fp, "BackgroundPressure");
    syncH<Vecd>(fp, "Position");
    {
        Real *pm = fp.getVariableDataByName<Real>("PerturbationMomentum");
        Real *dm = fp.getVariableDataByName<Real>("AcousticDensMode");
        Vecd *aref = fp.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
        Real *pbg = fp.getVariableDataByName<Real>("BackgroundPressure");
        Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            pm[i] = Real(1);
            dm[i] = dens_mode;
            aref[i] = Vecd::Zero();
            pbg[i] = rho0 * g * (L - pos[i][1]);
        }
        syncD<Real>(fp, "PerturbationMomentum");
        syncD<Real>(fp, "AcousticDensMode");
        syncD<Vecd>(fp, "ReferenceBodyAcceleration");
        syncD<Real>(fp, "BackgroundPressure");
    }

    auto run_acoustic_half = [&](bool first, Real dt)
    {
        if (use_central)
        {
            if (first)
                acoustic_1st_c.exec(dt);
            else
                acoustic_2nd_c.exec(dt);
        }
        else if (use_dissipative)
        {
            if (first)
                acoustic_1st_d.exec(dt);
            else
                acoustic_2nd_d.exec(dt);
        }
        else if (use_noriemann)
        {
            if (enable_kernel_corr)
            {
                if (first)
                    acoustic_1st_nr_c.exec(dt);
                else
                    acoustic_2nd_nr_c.exec(dt);
            }
            else
            {
                if (first)
                    acoustic_1st_nr.exec(dt);
                else
                    acoustic_2nd_nr.exec(dt);
            }
        }
        else if (enable_kernel_corr)
        {
            if (first)
                acoustic_1st_r_c.exec(dt);
            else
                acoustic_2nd_r_c.exec(dt);
        }
        else
        {
            if (first)
                acoustic_1st_r.exec(dt);
            else
                acoustic_2nd_r.exec(dt);
        }
    };
    auto run_acoustic = [&](Real dt)
    {
        run_acoustic_half(true, dt);
        run_acoustic_half(false, dt);
    };

    SingleVariable<Real> *sv_t = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    std::cout << std::scientific << std::setprecision(4)
              << "[" << case_id << "] scheme=" << scheme << " dens_sum=" << enable_dens_sum
              << " dens_reg=" << enable_dens_reg << " dens_norm=" << enable_dens_norm
              << " dens_mode=" << double(dens_mode) << " kcorr=" << enable_kernel_corr
              << " thermal_1d=" << (thermal_1d ? 1 : 0)
              << " fixed_T_linear=" << (fixed_T_linear ? 1 : 0)
              << " ke_budget=" << (ke_budget ? 1 : 0)
              << " L=" << double(L) << " dp=" << double(dp)
              << " Ra=" << double(Ra) << " Pr=" << double(Pr) << " Uchar=" << double(U_char)
              << " c0=" << double(c0) << " t*_end=" << double(tstar_end) << " t_end=" << double(end_time)
              << " n=" << fp.TotalRealParticles() << "\n";

    Real Umax_end = Real(0);
    int steps = 0;
    bool diverged = false;
    int high_u_streak = 0;
    bool dens_sigma0_ready = false;
    double dens_sigma0_mean = 0.0;
    double T_L2_rel = 0.0;
    int adv_count = 0;
    double ke_sum_d1 = 0.0, ke_sum_d2 = 0.0, ke_sum_ke0 = 0.0;
    size_t ke_nsamp = 0;
    std::ofstream ke_ofs;
    if (ke_budget)
    {
        ke_ofs.open(ke_budget_csv);
        ke_ofs << "tstar,adv,ac_sub,dt,ke0,ke1,ke2,dke_1st,dke_2nd\n";
        ke_ofs << std::setprecision(16);
    }
    if (thermal_1d)
    {
        // Doc204 §4.1: freeze hydro; only contact/inner Laplace + ApplyThermalDt.
        const Real dt_th_max = Real(0.1) * dp * dp / (alpha + TinyReal);
        while (sv_t->getValue() < end_time && !diverged)
        {
            Real dt_th = dt_th_max;
            const Real t_left = end_time - sv_t->getValue();
            if (dt_th > t_left)
                dt_th = t_left;
            if (!(dt_th > TinyReal))
                break;
            setWallThermal(wp, L, Th, Tc, T0, k, Real(1.0e6));
            syncH<Real>(fp, kLap);
            Real *lap = fp.getVariableDataByName<Real>(kLap);
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                lap[i] = Real(0);
            syncD<Real>(fp, kLap);
            lap_inner.exec();
            lap_wall.exec();
            apply_th.exec(dt_th);
            sv_t->incrementValue(dt_th);
            ++steps;
            if (steps % 500 == 0)
            {
                const Real tstar = sv_t->getValue() * alpha / (L * L);
                std::cout << "  step~" << steps << " t*=" << double(tstar) << " (thermal-1d frozen flow)\n";
            }
        }
        // Doc209 §5.1: volume-weighted L2 — √[Σ Vi(Ti−Tref)² / (ΔT² Σ Vi)]. Never divide by N twice.
        syncH<Vecd>(fp, "Position");
        syncH<Real>(fp, kT);
        syncH<Real>(fp, "VolumetricMeasure");
        Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        Real *T = fp.getVariableDataByName<Real>(kT);
        Real *Vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
        double se = 0, sv = 0;
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            const Real Tss = Th + (Tc - Th) * (pos[i][0] / (L + TinyReal));
            const double e = double(T[i] - Tss);
            const double w = double(Vol[i]);
            se += w * e * e;
            sv += w;
        }
        T_L2_rel = std::sqrt(se / (double(dT) * double(dT) * sv + 1e-30));
        std::cout << "[" << case_id << "] T_L2_rel_vs_linear=" << T_L2_rel
                  << " (Doc209_volume_weighted)\n";
    }
    while (!thermal_1d && sv_t->getValue() < end_time && !diverged)
    {
        if (enable_dens_sum)
        {
            density_sum.exec();
            if (enable_dens_norm && !dens_sigma0_ready)
            {
                syncH<Real>(fp, "DensitySummation");
                syncH<Real>(fp, "DensSumSigma0");
                const Real *rho_sum = fp.getVariableDataByName<Real>("DensitySummation");
                Real *sigma0 = fp.getVariableDataByName<Real>("DensSumSigma0");
                dens_sigma0_mean = 0.0;
                for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                {
                    sigma0[i] = std::max(rho_sum[i], TinyReal);
                    dens_sigma0_mean += double(sigma0[i]);
                }
                dens_sigma0_mean /= double(std::max<size_t>(fp.TotalRealParticles(), 1));
                syncD<Real>(fp, "DensSumSigma0");
                dens_sigma0_ready = true;
                std::cout << "[" << case_id << "] dens_sum_normalize sigma0_mean=" << dens_sigma0_mean
                          << " packing=" << (dens_sigma0_mean / (double(rho0) + 1e-30)) << "\n";
            }
            if (enable_dens_reg)
                density_reg.exec();
            if (enable_dens_norm)
            {
                syncH<Real>(fp, "Density");
                syncH<Real>(fp, "DensitySummation");
                syncH<Real>(fp, "DensSumSigma0");
                Real *rho = fp.getVariableDataByName<Real>("Density");
                const Real *rho_sum = fp.getVariableDataByName<Real>("DensitySummation");
                const Real *sigma0 = fp.getVariableDataByName<Real>("DensSumSigma0");
                for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                    rho[i] = rho0 * (rho_sum[i] / (sigma0[i] + TinyReal));
                syncD<Real>(fp, "Density");
            }
        }
        // Doc105/French: rewrite p_h = ρ0 g (H−y) every advection (particles moved).
        syncH<Vecd>(fp, "Position");
        syncH<Real>(fp, "BackgroundPressure");
        {
            Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
            Real *pbg = fp.getVariableDataByName<Real>("BackgroundPressure");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                pbg[i] = rho0 * g * (L - pos[i][1]);
            syncD<Real>(fp, "BackgroundPressure");
        }
        adv_setup.exec();
        zeroVec(fp, "Force");
        if (enable_kernel_corr)
            lin_corr.exec();
        const Real dt_adv = std::min(adv_dt.exec(), adv_visc_dt.exec());
        // ViscousForceWithWallCK is ForcePriorCK — do NOT overwrite ForcePrior afterwards.
        viscous.exec();
        boussinesq.exec();

        const Real tstar_now = sv_t->getValue() * alpha / (L * L);
        const bool sample_ke =
            ke_budget && (adv_count < 120 || (adv_count % 100 == 0) || tstar_now > tstar_end - Real(0.03));

        Real relax = Real(0);
        size_t ac_sub = 0;
        while (relax < dt_adv && ac_sub < 2000)
        {
            Real dt = ac_dt.exec();
            const Real dt_left = dt_adv - relax;
            if (dt > dt_left)
                dt = dt_left;
            if (!(dt > TinyReal) || !std::isfinite(double(dt)))
            {
                diverged = true;
                break;
            }
            if (sample_ke)
            {
                const double ke0 = computeParticleKE(fp);
                run_acoustic_half(true, dt);
                const double ke1 = computeParticleKE(fp);
                run_acoustic_half(false, dt);
                const double ke2 = computeParticleKE(fp);
                const double d1 = ke1 - ke0;
                const double d2 = ke2 - ke1;
                ke_sum_d1 += d1;
                ke_sum_d2 += d2;
                ke_sum_ke0 += ke0;
                ++ke_nsamp;
                const Real tstar_s = sv_t->getValue() * alpha / (L * L);
                ke_ofs << double(tstar_s) << "," << adv_count << "," << ac_sub << "," << double(dt) << ","
                       << ke0 << "," << ke1 << "," << ke2 << "," << d1 << "," << d2 << "\n";
            }
            else
            {
                run_acoustic(dt);
            }
            boussinesq.exec(); // doc77c: refresh buoyancy each acoustic substep
            relax += dt;
            sv_t->incrementValue(dt);
            ++steps;
            ++ac_sub;
        }
        ++adv_count;
        if (diverged)
            break;
        update_pos.exec();

        if (fixed_T_linear)
        {
            // Remap T from fixed Eulerian field after particles move; skip thermal PDE.
            applyFixedLinearTemperature(fp, L, Th, Tc, T0);
            setWallThermal(wp, L, Th, Tc, T0, k, Real(1.0e6));
        }
        else
        {
            Real t_th = Real(0);
            const Real dt_th_max = Real(0.1) * dp * dp / (alpha + TinyReal);
            size_t n_th = 0;
            while (t_th < relax && n_th < 10000)
            {
                const Real dt_th = std::min(dt_th_max, relax - t_th);
                setWallThermal(wp, L, Th, Tc, T0, k, Real(1.0e6));
                syncH<Real>(fp, kLap);
                Real *lap = fp.getVariableDataByName<Real>(kLap);
                for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                    lap[i] = Real(0);
                syncD<Real>(fp, kLap);
                lap_inner.exec();
                lap_wall.exec();
                apply_th.exec(dt_th);
                t_th += dt_th;
                ++n_th;
            }
        }

        // Skip sort: current ViscousForce/Boussinesq force arrays are not evolving.
        (void)particle_sort;
        fluid_cells.exec();
        wall_cells.exec();
        rel.exec();

        syncH<Vecd>(fp, "Velocity");
        syncH<Real>(fp, "Density");
        Vecd *vel_now = fp.getVariableDataByName<Vecd>("Velocity");
        Real *rho_now = fp.getVariableDataByName<Real>("Density");
        Real um = Real(0);
        Real rho_min = rho_now[0], rho_max = rho_now[0];
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            um = std::max(um, vel_now[i].norm());
            rho_min = std::min(rho_min, rho_now[i]);
            rho_max = std::max(rho_max, rho_now[i]);
        }
        const Real um_star = um / U_char;
        // Coarse mesh can spike briefly (GPU); require sustained runaway or hard NaN/∞.
        if (!std::isfinite(double(um_star)) || um_star > Real(2000))
        {
            std::cout << "[" << case_id << "] ABORT diverge Umax*=" << double(um_star) << " at t*="
                      << double(sv_t->getValue() * alpha / (L * L))
                      << " rho=[" << double(rho_min) << "," << double(rho_max) << "]\n";
            diverged = true;
            Umax_end = um;
            break;
        }
        if (um_star > Real(400))
            ++high_u_streak;
        else
            high_u_streak = 0;
        if (high_u_streak >= 8)
        {
            std::cout << "[" << case_id << "] ABORT sustained high Umax*=" << double(um_star) << " at t*="
                      << double(sv_t->getValue() * alpha / (L * L))
                      << " rho=[" << double(rho_min) << "," << double(rho_max) << "]\n";
            diverged = true;
            Umax_end = um;
            break;
        }
        if (steps % 200 == 0)
        {
            const Real tstar = sv_t->getValue() * alpha / (L * L);
            std::cout << "  step~" << steps << " t*=" << double(tstar) << " Umax*=" << double(um_star)
                      << " rho=[" << double(rho_min) << "," << double(rho_max) << "]\n";
        }
        Umax_end = um;
    }

    Real Vmax_hcline = Real(0), Umax_vcline = Real(0); // DVD Table V probes
    Real x_at_Vmax = Real(0), y_at_Umax = Real(0);
    Real Nu_grad_h = Real(0), Nu_grad_c = Real(0), Nu_contact_h = Real(0), Nu_contact_c = Real(0);
    Real Nu_fit_h = Real(0), Nu_fit_c = Real(0);
    double Q_hot = 0, Q_cold = 0;
    Real dxs_mean_hot = Real(0);
    Real x_first_mean_hot = Real(0);
    if (!diverged)
    {
        syncH<Vecd>(fp, "Velocity");
        syncH<Vecd>(fp, "Position");
        syncH<Real>(fp, kT);
        syncH<Real>(fp, "VolumetricMeasure");
        Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        Real *T = fp.getVariableDataByName<Real>(kT);
        Real *Vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
        Umax_end = Real(0);
        const Real grad_layer = Real(2) * dp;
        const Real fit_layer = Real(4) * dp; // Doc204 independent wall-normal fit window
        // Doc209 §5.2: geometric centerline IDW (no band particle-max, no 90% hot-side preference).
        const Real support = centerline_support_factor * dp;
        const int n_cl = 81;
        std::ofstream clf(centerline_csv);
        clf << std::setprecision(16)
            << "line,s_over_L,u_star,v_star,theta,nearest_dist_over_L,n_support,support_dp_factor\n";
        Real vmax_pos = Real(0), x_vmax = Real(0);
        Real umax_pos = Real(0), y_umax = Real(0);
        for (int k = 0; k < n_cl; ++k)
        {
            const Real s = (Real(k) + Real(0.5)) / Real(n_cl) * L;
            // horizontal midplane y=L/2: v(x)
            {
                auto sm = sampleCenterlineIDW(pos, vel, T, fp.TotalRealParticles(), L, Tc, dT, U_char, s,
                                              Real(0.5) * L, support);
                clf << "y_eq_half," << double(s / L) << "," << double(sm.u_star) << "," << double(sm.v_star)
                    << "," << double(sm.theta) << "," << double(sm.nearest_dist_over_L) << "," << sm.n_support
                    << "," << double(centerline_support_factor) << "\n";
                if (sm.n_support > 0 && sm.v_star > vmax_pos)
                {
                    vmax_pos = sm.v_star;
                    x_vmax = s / L;
                }
            }
            // vertical midplane x=L/2: u(y)
            {
                auto sm = sampleCenterlineIDW(pos, vel, T, fp.TotalRealParticles(), L, Tc, dT, U_char,
                                              Real(0.5) * L, s, support);
                clf << "x_eq_half," << double(s / L) << "," << double(sm.u_star) << "," << double(sm.v_star)
                    << "," << double(sm.theta) << "," << double(sm.nearest_dist_over_L) << "," << sm.n_support
                    << "," << double(centerline_support_factor) << "\n";
                if (sm.n_support > 0 && sm.u_star > umax_pos)
                {
                    umax_pos = sm.u_star;
                    y_umax = s / L;
                }
            }
        }
        clf.close();
        // sensitivity: support=1.5 dp
        Real vmax_alt = Real(0), umax_alt = Real(0);
        const Real support_alt = Real(1.5) * dp;
        for (int k = 0; k < n_cl; ++k)
        {
            const Real s = (Real(k) + Real(0.5)) / Real(n_cl) * L;
            auto smv = sampleCenterlineIDW(pos, vel, T, fp.TotalRealParticles(), L, Tc, dT, U_char, s,
                                           Real(0.5) * L, support_alt);
            auto smu = sampleCenterlineIDW(pos, vel, T, fp.TotalRealParticles(), L, Tc, dT, U_char,
                                           Real(0.5) * L, s, support_alt);
            vmax_alt = std::max(vmax_alt, smv.v_star);
            umax_alt = std::max(umax_alt, smu.u_star);
        }
        Vmax_hcline = vmax_pos;
        x_at_Vmax = x_vmax;
        Umax_vcline = umax_pos;
        y_at_Umax = y_umax;
        std::cout << "[" << case_id << "] centerline_IDW support_dp=" << double(centerline_support_factor)
                  << " Vmax*=" << double(Vmax_hcline) << "@x/L=" << double(x_at_Vmax)
                  << " Umax*=" << double(Umax_vcline) << "@y/L=" << double(y_at_Umax)
                  << " sens_1p5dp V*=" << double(vmax_alt) << " U*=" << double(umax_alt)
                  << " csv=" << centerline_csv << "\n";

        double Nu_h = 0, Nu_c = 0, wh = 0, wc = 0;
        // least-squares T = a + b * x_walldist  (vol-weighted)
        double sh_x = 0, sh_t = 0, sh_xx = 0, sh_xt = 0, sh_w = 0;
        double sc_x = 0, sc_t = 0, sc_xx = 0, sc_xt = 0, sc_w = 0;
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            const Real um = vel[i].norm();
            Umax_end = std::max(Umax_end, um);
            const Real x = pos[i][0], y = pos[i][1];
            (void)y;
            const Real theta = (T[i] - Tc) / (dT + TinyReal);
            // One-sided (Tw-Tp)/x: first-layer Dirichlet contact biases Tp→Tw → under-Nu.
            if (x > TinyReal && x < grad_layer)
            {
                const Real dxs = x / L;
                Nu_h += double(((Real(1) - theta) / (dxs + TinyReal)) * Vol[i]);
                wh += double(Vol[i]);
                dxs_mean_hot += dxs * Vol[i];
                x_first_mean_hot += (x / L) * Vol[i];
            }
            if (x < L - TinyReal && x > L - grad_layer)
            {
                const Real dxs = (L - x) / L;
                Nu_c += double((theta / (dxs + TinyReal)) * Vol[i]);
                wc += double(Vol[i]);
            }
            if (x > TinyReal && x < fit_layer)
            {
                const double w = double(Vol[i]);
                const double xw = double(x);
                const double tw = double(T[i]);
                sh_w += w;
                sh_x += w * xw;
                sh_t += w * tw;
                sh_xx += w * xw * xw;
                sh_xt += w * xw * tw;
            }
            if (x < L - TinyReal && x > L - fit_layer)
            {
                const double w = double(Vol[i]);
                const double xw = double(L - x); // distance from cold wall
                const double tw = double(T[i]);
                sc_w += w;
                sc_x += w * xw;
                sc_t += w * tw;
                sc_xx += w * xw * xw;
                sc_xt += w * xw * tw;
            }
        }
        Nu_grad_h = wh > TinyReal ? Real(Nu_h / wh) : Real(0);
        Nu_grad_c = wc > TinyReal ? Real(Nu_c / wc) : Real(0);
        if (wh > TinyReal)
        {
            dxs_mean_hot = Real(dxs_mean_hot / wh);
            x_first_mean_hot = Real(x_first_mean_hot / wh);
        }
        // Nu from dT/dn: hot wall T decreases with x → b_h<0 → Nu = -b L/ΔT
        // cold: T vs dist-from-cold; toward fluid T rises → b_c>0 → Nu = b L/ΔT
        auto slope = [](double sx, double st, double sxx, double sxt, double sw) -> double {
            const double det = sw * sxx - sx * sx;
            if (!(std::abs(det) > 1e-30) || !(sw > 0))
                return 0.0;
            return (sw * sxt - sx * st) / det;
        };
        const double bh = slope(sh_x, sh_t, sh_xx, sh_xt, sh_w);
        const double bc = slope(sc_x, sc_t, sc_xx, sc_xt, sc_w);
        Nu_fit_h = Real((-bh) * double(L) / (double(dT) + 1e-30));
        Nu_fit_c = Real((bc) * double(L) / (double(dT) + 1e-30));

        // Contact Nu: one thermal snapshot — Q = Σ(−lap_contact Vol), Nu = (Q/A)/(k ΔT/L)
        syncH<Real>(fp, kLap);
        Real *lap = fp.getVariableDataByName<Real>(kLap);
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            lap[i] = Real(0);
        syncD<Real>(fp, kLap);
        lap_inner.exec();
        std::vector<Real> lap_inner_buf(fp.TotalRealParticles());
        syncH<Real>(fp, kLap);
        lap = fp.getVariableDataByName<Real>(kLap);
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            lap_inner_buf[i] = lap[i];
        setWallThermal(wp, L, Th, Tc, T0, k, Real(1.0e6));
        lap_wall.exec();
        syncH<Real>(fp, kLap);
        syncH<Vecd>(fp, "Position");
        syncH<Real>(fp, "VolumetricMeasure");
        lap = fp.getVariableDataByName<Real>(kLap);
        pos = fp.getVariableDataByName<Vecd>("Position");
        Vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
        const Real wall_band = std::max(Real(2) * dp, Real(0.1) * L);
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            const Real lap_c = lap[i] - lap_inner_buf[i];
            const double dE = -double(lap_c * Vol[i]); // power into particle (Doc77c raw)
            const Real x = pos[i][0];
            if (x <= wall_band && x <= L - x)
                Q_hot += dE;
            else if (L - x <= wall_band)
                Q_cold += dE;
        }
        const Real A_wall = L; // 2D unit depth
        const Real q_scale = k * dT / (L + TinyReal);
        Nu_contact_h = Real(Q_hot / (double(A_wall) * double(q_scale) + 1e-30));
        Nu_contact_c = Real((-Q_cold) / (double(A_wall) * double(q_scale) + 1e-30));
    }
    const Real Umax_star = Umax_end / U_char;
    const Real Nu_avg = Real(0.5) * (Nu_grad_h + Nu_grad_c);
    const Real Nu_contact_avg = Real(0.5) * (Nu_contact_h + Nu_contact_c);
    const double q_imbalance =
        (std::abs(Q_hot) + std::abs(Q_cold) > 1e-30)
            ? std::abs(Q_hot + Q_cold) / (0.5 * (std::abs(Q_hot) + std::abs(Q_cold)))
            : 0.0;
    const bool energy_balance_pass =
        thermal_1d ? (q_imbalance < 0.05 && std::abs(double(Nu_contact_avg) - 1.0) < 0.05)
                   : (q_imbalance < 0.05);
    const Real Nu_fit_avg = Real(0.5) * (Nu_fit_h + Nu_fit_c);
    // Doc204 mid/fine gates vs DVD Table V (contact Nu only). smoke ≠ quantitative.
    const double err_Nu = std::abs(double(Nu_contact_avg) - 4.519) / 4.519;
    const double err_V = std::abs(double(Vmax_hcline) - 68.59) / 68.59;
    const double err_U = std::abs(double(Umax_vcline) - 34.73) / 34.73;
    const double err_x = std::abs(double(x_at_Vmax) - 0.066);
    const double err_y = std::abs(double(y_at_Umax) - 0.855);
    const bool reference_match =
        (!thermal_1d && !fixed_T_linear && case_id == "C1_dvd_exact_ra1e5_pr071" && !diverged &&
         err_Nu <= 0.05 && err_V <= 0.10 && err_U <= 0.10 && err_x <= 0.05 && err_y <= 0.05);
    // thermal-1d: Doc204 — L2≤2%, independent fit flux≤3%; contact≤5% on coarse n40.
    const bool quantitative_pass =
        thermal_1d ? (energy_balance_pass && T_L2_rel < 0.02 &&
                      std::abs(double(Nu_fit_avg) - 1.0) <= 0.03 &&
                      std::abs(double(Nu_contact_avg) - 1.0) <= 0.05)
                   : reference_match;
    // Smoke/consistency: circulation on, not runaway. C1 DVD compare only meaningful for air Ra~1e5.
    // thermal-1d: frozen u≈0; pass if Nu_contact near 1 and T_L2 small.
    const Real umin_gate = thermal_1d ? Real(0) : ((case_id == "C2") ? Real(0.5) : Real(1));
    const bool pass_smoke = thermal_1d
                                ? (!diverged && std::abs(double(Nu_contact_avg) - 1.0) < 0.15 && T_L2_rel < 0.05)
                                : (!diverged && Umax_star > umin_gate && Umax_star < Real(300));
    std::cout << std::scientific << std::setprecision(4)
              << "[" << case_id << "] done Umax*=" << double(Umax_star)
              << " Vmax_hcline*=" << double(Vmax_hcline) << " @x/L=" << double(x_at_Vmax)
              << " (dvd 68.59@0.066) Umax_vcline*=" << double(Umax_vcline)
              << " @y/L=" << double(y_at_Umax) << " (dvd 34.73@0.855)"
              << " Nu_grad_h=" << double(Nu_grad_h) << " Nu_grad_c=" << double(Nu_grad_c)
              << " Nu_grad_avg=" << double(Nu_avg)
              << " Nu_fit_h=" << double(Nu_fit_h) << " Nu_fit_c=" << double(Nu_fit_c)
              << " Nu_fit_avg=" << double(Nu_fit_avg)
              << " Nu_contact_h=" << double(Nu_contact_h) << " Nu_contact_c=" << double(Nu_contact_c)
              << " Nu_contact_avg=" << double(Nu_contact_avg)
              << (thermal_1d ? " (analytic Nu=1)" : " (dvd 4.519)")
              << " Q_hot=" << Q_hot << " Q_cold=" << Q_cold
              << " dxs_mean_hot=" << double(dxs_mean_hot) << " x_first_mean_hot=" << double(x_first_mean_hot)
              << " q_imbalance=" << q_imbalance
              << " -> " << (pass_smoke ? "SMOKE_PASS" : "SMOKE_FAIL")
              << " quant=" << (quantitative_pass ? "PASS" : "FAIL")
              << " ebal=" << (energy_balance_pass ? "PASS" : "FAIL") << "\n";

    std::ofstream ofs(out_json);
    const char *status_note = "C1_paper_table_PARTIAL";
    if (case_id == "C2")
        status_note = "CONSISTENCY_ONLY";
    else if (fixed_T_linear)
        status_note = "Doc209_§10.1_FIXED_T_LINEAR_HYDRO_DIAGNOSIS";
    else if (ke_budget)
        status_note = "Doc209_§10.3_ACOUSTIC_KE_BUDGET";
    else if (case_id == "C1_dvd_exact_ra1e5_pr071")
        status_note = "C1_dvd_exact_FOR_DVD_TABLE_V";
    else if (thermal_1d)
        status_note = "C1_thermal_1d_FLUX_AUDIT";
    if (ke_budget)
    {
        if (ke_ofs)
            ke_ofs.close();
        const double mean_d1 = ke_nsamp ? ke_sum_d1 / double(ke_nsamp) : 0.0;
        const double mean_d2 = ke_nsamp ? ke_sum_d2 / double(ke_nsamp) : 0.0;
        const double mean_ke0 = ke_nsamp ? ke_sum_ke0 / double(ke_nsamp) : 0.0;
        std::cout << std::scientific << std::setprecision(6)
                  << "[" << case_id << "] ke_budget samples=" << ke_nsamp
                  << " sum_dke_1st=" << ke_sum_d1 << " sum_dke_2nd=" << ke_sum_d2
                  << " mean_dke_1st=" << mean_d1 << " mean_dke_2nd=" << mean_d2
                  << " mean_ke0=" << mean_ke0
                  << " csv=" << ke_budget_csv << "\n";
    }
    ofs << std::setprecision(16) << "{\n  \"doc\": \"204\", \"case\": \"" << case_id << "\",\n"
        << "  \"status_note\": \"" << status_note << "\",\n"
        << "  \"scheme\": \"" << scheme << "+P1+BoussinesqCK+Viscous+ContactLaplace\",\n"
        << "  \"thermal_1d\": " << (thermal_1d ? "true" : "false") << ",\n"
        << "  \"fixed_T_linear\": " << (fixed_T_linear ? "true" : "false") << ",\n"
        << "  \"ke_budget\": " << (ke_budget ? "true" : "false") << ",\n"
        << "  \"ke_nsamp\": " << ke_nsamp << ",\n"
        << "  \"ke_sum_dke_1st\": " << ke_sum_d1 << ",\n"
        << "  \"ke_sum_dke_2nd\": " << ke_sum_d2 << ",\n"
        << "  \"ke_mean_dke_1st\": " << (ke_nsamp ? ke_sum_d1 / double(ke_nsamp) : 0.0) << ",\n"
        << "  \"ke_mean_dke_2nd\": " << (ke_nsamp ? ke_sum_d2 / double(ke_nsamp) : 0.0) << ",\n"
        << "  \"ke_mean_ke0\": " << (ke_nsamp ? ke_sum_ke0 / double(ke_nsamp) : 0.0) << ",\n"
        << "  \"dens_sum\": " << (enable_dens_sum ? "true" : "false")
        << ", \"dens_reg\": " << (enable_dens_reg ? "true" : "false")
        << ", \"dens_norm\": " << (enable_dens_norm ? "true" : "false")
        << ", \"kernel_corr\": " << (enable_kernel_corr ? "true" : "false") << ",\n"
        << "  \"Ra\": " << double(Ra) << ", \"Pr\": " << double(Pr) << ",\n"
        << "  \"n_div\": " << double(n_div) << ", \"tstar_end\": " << double(tstar_end) << ",\n"
        << "  \"Umax_star\": " << double(Umax_star)
        << ", \"Vmax_hcline_star\": " << double(Vmax_hcline)
        << ", \"x_at_Vmax_over_L\": " << double(x_at_Vmax)
        << ", \"Umax_vcline_star\": " << double(Umax_vcline)
        << ", \"y_at_Umax_over_L\": " << double(y_at_Umax)
        << ", \"Nu_grad_hot\": " << double(Nu_grad_h) << ", \"Nu_grad_cold\": " << double(Nu_grad_c)
        << ", \"Nu_grad_avg\": " << double(Nu_avg)
        << ", \"Nu_fit_hot\": " << double(Nu_fit_h) << ", \"Nu_fit_cold\": " << double(Nu_fit_c)
        << ", \"Nu_fit_avg\": " << double(Nu_fit_avg)
        << ", \"Nu_contact_hot\": " << double(Nu_contact_h) << ", \"Nu_contact_cold\": " << double(Nu_contact_c)
        << ", \"Nu_contact_avg\": " << double(Nu_contact_avg) << ",\n"
        << "  \"Q_hot\": " << Q_hot << ", \"Q_cold\": " << Q_cold
        << ", \"q_imbalance\": " << q_imbalance
        << ", \"dxs_mean_hot\": " << double(dxs_mean_hot)
        << ", \"x_first_mean_hot\": " << double(x_first_mean_hot)
        << ", \"T_L2_rel_vs_linear\": " << T_L2_rel << ",\n"
        << "  \"dvd_Nu\": 4.519, \"dvd_Vmax\": 68.59, \"dvd_Umax\": 34.73,\n"
        << "  \"err_Nu_contact\": " << err_Nu << ", \"err_Vmax\": " << err_V << ", \"err_Umax\": " << err_U
        << ", \"err_x_over_L\": " << err_x << ", \"err_y_over_L\": " << err_y << ",\n"
        << "  \"diverged\": " << (diverged ? "true" : "false")
        << ", \"smoke_pass\": " << (pass_smoke ? "true" : "false")
        << ", \"energy_balance_pass\": " << (energy_balance_pass ? "true" : "false")
        << ", \"reference_match\": " << (reference_match ? "true" : "false")
        << ", \"quantitative_pass\": " << (quantitative_pass ? "true" : "false") << "\n}\n";
    return pass_smoke ? 0 : 1;
}
