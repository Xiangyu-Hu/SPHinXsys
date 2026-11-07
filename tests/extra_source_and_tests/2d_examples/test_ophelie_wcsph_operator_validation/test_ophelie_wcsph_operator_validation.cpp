/**
 * Doc203 O0: central WCSPH operator validation.
 *   --case o0a : host paper-form vs NoRiemann mirror (no SPH)
 *   --case o0b : AcousticStep(NoRiemann) GPU/CK RHS vs independent host sum
 */
#include "ophelie_central_wcsph_pair_ref.h"
#include "sphinxsys.h"

// Minimal device sync (avoid linking extra_sources_2d / 3D-only units).
namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{
template <typename DataType>
inline void syncVariableToHost(BaseParticles &particles, const std::string &name)
{
#if SPHINXSYS_USE_SYCL
    particles.template getVariableByName<DataType>(name)->prepareForOutput(execution::par_device);
#else
    (void)particles;
    (void)name;
#endif
}
template <typename DataType>
inline void syncVariableToDevice(BaseParticles &particles, const std::string &name)
{
#if SPHINXSYS_USE_SYCL
    particles.template getVariableByName<DataType>(name)->finalizeLoadIn(execution::par_device);
#else
    (void)particles;
    (void)name;
#endif
}
} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using namespace SPH::electromagnetics::ophelie::hydro;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
struct Metrics
{
    double max_a_rel = 0;
    double max_a_abs = 0;
    double max_drho_rel = 0;
    double max_drho_abs = 0;
    double pair_force_rel = 0;
    size_t n_interior = 0;
    bool pass = false;
};

std::vector<PairState2D> makeLattice(int n, double dp, double rho0, bool perturb)
{
    std::vector<PairState2D> P;
    P.reserve(size_t(n * n));
    const double Vol = dp * dp;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            PairState2D p;
            p.x = (i + 0.5) * dp;
            p.y = (j + 0.5) * dp;
            if (perturb)
            {
                p.x += 0.02 * dp * std::sin(3.1 * p.x + 1.7 * p.y);
                p.y += 0.02 * dp * std::cos(2.3 * p.x - 1.1 * p.y);
            }
            p.rho = rho0;
            p.Vol = Vol;
            p.mass = rho0 * Vol;
            P.push_back(p);
        }
    return P;
}

bool isInterior(double x, double y, double L, double margin)
{
    return x > margin && x < L - margin && y > margin && y < L - margin;
}

Metrics runHostIdentity(const std::string &name, std::vector<PairState2D> P, double h, double L, double margin,
                        double tol, double a_char, double drho_char)
{
    Metrics m;
    double sumFx = 0, sumFy = 0, sumAbs = 0;
    for (size_t i = 0; i < P.size(); ++i)
    {
        const auto paper = centralRhsPaperForm(i, P, h);
        const auto mirror = centralRhsNoRiemannMirror(i, P, h);
        sumFx += paper.force_x;
        sumFy += paper.force_y;
        sumAbs += std::abs(paper.force_x) + std::abs(paper.force_y);
        if (!isInterior(P[i].x, P[i].y, L, margin))
            continue;
        ++m.n_interior;
        const double dax = std::abs(paper.ax - mirror.ax);
        const double day = std::abs(paper.ay - mirror.ay);
        const double ddr = std::abs(paper.drho_dt - mirror.drho_dt);
        m.max_a_abs = std::max(m.max_a_abs, std::max(dax, day));
        m.max_drho_abs = std::max(m.max_drho_abs, ddr);
        m.max_a_rel = std::max(m.max_a_rel, dax / a_char);
        m.max_a_rel = std::max(m.max_a_rel, day / a_char);
        m.max_drho_rel = std::max(m.max_drho_rel, ddr / drho_char);
    }
    m.pair_force_rel = sumAbs > 0 ? std::sqrt(sumFx * sumFx + sumFy * sumFy) / sumAbs : 0;
    m.pass = (m.max_a_rel < tol) && (m.max_drho_rel < tol);
    std::cout << std::scientific << std::setprecision(3) << "[O0a] " << name << " n_int=" << m.n_interior
              << " max_a_rel=" << m.max_a_rel << " max_drho_rel=" << m.max_drho_rel
              << " pair|ΣF|/Σ|F|=" << m.pair_force_rel << " -> " << (m.pass ? "PASS" : "FAIL") << "\n";
    return m;
}

int runO0a(const std::string &out_json)
{
    const int n = 21;
    const double dp = 1.0 / n;
    const double L = 1.0;
    const double h = 1.3 * dp;
    const double rho0 = 1000.0;
    const double margin = 2.0 * h;
    const double tol = 1.0e-10;
    const double p_ref = 1.0e5;
    const double a_char = std::max(1.0, p_ref / (rho0 * h));
    const double U_ref = 0.01;
    const double drho_char = std::max(1.0, rho0 * U_ref / h);

    auto applyConstP = [](std::vector<PairState2D> &P, double p0)
    {
        for (auto &q : P)
            q.p = p0;
    };
    auto applyLinearP = [](std::vector<PairState2D> &P, double g)
    {
        for (auto &q : P)
            q.p = g * q.y;
    };
    auto applySmoothU = [](std::vector<PairState2D> &P)
    {
        for (auto &q : P)
        {
            q.u = 0.01 * std::sin(2.0 * M_PI * q.x) * std::cos(2.0 * M_PI * q.y);
            q.v = -0.01 * std::cos(2.0 * M_PI * q.x) * std::sin(2.0 * M_PI * q.y);
        }
    };

    int fails = 0;
    std::vector<Metrics> all;
    {
        auto P = makeLattice(n, dp, rho0, false);
        applyConstP(P, p_ref);
        all.push_back(runHostIdentity("const_p_zero_u_regular", P, h, L, margin, tol, a_char, drho_char));
    }
    {
        auto P = makeLattice(n, dp, rho0, true);
        applyConstP(P, p_ref);
        all.push_back(runHostIdentity("const_p_zero_u_perturb", P, h, L, margin, tol, a_char, drho_char));
    }
    {
        auto P = makeLattice(n, dp, rho0, false);
        applyLinearP(P, rho0 * 9.81);
        all.push_back(runHostIdentity("linear_p_zero_u", P, h, L, margin, tol, a_char, drho_char));
    }
    {
        auto P = makeLattice(n, dp, rho0, false);
        applyConstP(P, p_ref);
        applySmoothU(P);
        all.push_back(runHostIdentity("const_p_smooth_u", P, h, L, margin, tol, a_char, drho_char));
    }
    {
        auto P = makeLattice(n, dp, rho0, false);
        applyLinearP(P, rho0 * 9.81);
        applySmoothU(P);
        all.push_back(runHostIdentity("linear_p_smooth_u", P, h, L, margin, tol, a_char, drho_char));
    }
    for (const auto &m : all)
        if (!m.pass)
            ++fails;

    std::ofstream ofs(out_json);
    ofs << std::setprecision(16) << "{\n  \"doc\": \"203\", \"case\": \"O0a\",\n  \"fails\": " << fails << "\n}\n";
    std::cout << "O0a summary fails=" << fails << "\n";
    return fails == 0 ? 0 : 1;
}

void hostZeroVec(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Vecd>(particles, name);
    Vecd *v = particles.getVariableDataByName<Vecd>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Vecd::Zero();
    syncVariableToDevice<Vecd>(particles, name);
}

void hostZeroReal(BaseParticles &particles, const std::string &name)
{
    syncVariableToHost<Real>(particles, name);
    Real *v = particles.getVariableDataByName<Real>(name);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        v[i] = Real(0);
    syncVariableToDevice<Real>(particles, name);
}

enum class FieldCase
{
    ConstPZeroU,
    LinearPZeroU,
    ConstPSmoothU
};

int runO0b(const std::string &out_json, int argc, char **argv)
{
    // Geometry: unit square fluid, Inner only (no wall). Interior band for full-support check.
    const Real L = Real(1);
    const Real dp = Real(1) / Real(21);
    const Real rho0 = Real(1000);
    const Real c0 = Real(10); // small c0 unused at dt=0; keeps DensityFromPressure mild
    const Real p_ref = Real(1.0e5);

    // Filter Doc203-only flags so SPHSystem boost::program_options does not reject them.
    std::vector<char *> sph_argv;
    sph_argv.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--case" || a == "--out")
        {
            if (i + 1 < argc)
                ++i;
            continue;
        }
        sph_argv.push_back(argv[i]);
    }
    int sph_argc = static_cast<int>(sph_argv.size());

    BoundingBoxd system_domain_bounds(Vec2d(-dp, -dp), Vec2d(L + dp, L + dp));
    SPHSystem sph_system(system_domain_bounds, dp);
    sph_system.handleCommandlineOptions(sph_argc, sph_argv.data());

    GeometricShapeBox fluid_shape(Transform(Vec2d(0.5 * L, 0.5 * L)), Vec2d(0.5 * L, 0.5 * L), "Fluid");
    FluidBody fluid(sph_system, fluid_shape);
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, c0);
    fluid.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    BaseParticles &particles = fluid.getBaseParticles();
    auto &mat = DynamicCast<WeaklyCompressibleFluid>(&fluid, fluid.getMatterMaterial());

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> cells(fluid);
    UpdateRelation<MainExecutionPolicy, Inner<>> rel(fluid_inner);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);

    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep1stHalf<
                              Inner<OneLevel, NoRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_1st(fluid_inner);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<
                              Inner<OneLevel, NoRiemannSolverCK, NoKernelCorrectionCK>>>
        acoustic_2nd(fluid_inner);

    cells.exec();
    rel.exec();
    adv_setup.exec();

    const Real h = fluid.getSPHAdaptation().ReferenceSmoothingLength();
    const Real margin = Real(2) * h;
    // Single precision (SYCL float) thresholds from Doc203 O0.
    const Real tol_rel = Real(1.0e-4);
    const Real a_char = std::max(Real(1), p_ref / (rho0 * h));
    const Real U_ref = Real(0.01);
    const Real drho_char = std::max(Real(1), rho0 * U_ref / h);

    auto applyField = [&](FieldCase fc)
    {
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Vecd>(particles, "Velocity");
        syncVariableToHost<Real>(particles, "Density");
        syncVariableToHost<Real>(particles, "Pressure");
        syncVariableToHost<Real>(particles, "BackgroundPressure");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<Real>(particles, "Mass");

        Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
        Real *rho = particles.getVariableDataByName<Real>("Density");
        Real *p = particles.getVariableDataByName<Real>("Pressure");
        Real *pbg = particles.getVariableDataByName<Real>("BackgroundPressure");
        Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        Real *mass = particles.getVariableDataByName<Real>("Mass");
        const size_t n = particles.TotalRealParticles();
        for (size_t i = 0; i < n; ++i)
        {
            pbg[i] = Real(0);
            Real p_des = p_ref;
            if (fc == FieldCase::LinearPZeroU)
                p_des = rho0 * Real(9.81) * pos[i][1];
            rho[i] = mat.DensityFromPressure(p_des);
            p[i] = p_des; // Init will reset from EOS; keep consistent
            Vol[i] = mass[i] / rho[i]; // AdvectionStepSetup semantics after ρ write
            vel[i] = Vecd::Zero();
            if (fc == FieldCase::ConstPSmoothU)
            {
                const Real x = pos[i][0], y = pos[i][1];
                vel[i][0] = U_ref * std::sin(Real(2) * Real(M_PI) * x) * std::cos(Real(2) * Real(M_PI) * y);
                vel[i][1] = -U_ref * std::cos(Real(2) * Real(M_PI) * x) * std::sin(Real(2) * Real(M_PI) * y);
            }
        }
        syncVariableToDevice<Vecd>(particles, "Velocity");
        syncVariableToDevice<Real>(particles, "Density");
        syncVariableToDevice<Real>(particles, "Pressure");
        syncVariableToDevice<Real>(particles, "BackgroundPressure");
        syncVariableToDevice<Real>(particles, "VolumetricMeasure");
        hostZeroVec(particles, "Force");
        hostZeroVec(particles, "ForcePrior");
        hostZeroReal(particles, "DensityChangeRate");
    };

    auto compareAfterKernels = [&](const std::string &name, bool check_force, bool check_drho) -> Metrics
    {
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Vecd>(particles, "Velocity");
        syncVariableToHost<Vecd>(particles, "Force");
        syncVariableToHost<Real>(particles, "Density");
        syncVariableToHost<Real>(particles, "Pressure");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<Real>(particles, "Mass");
        syncVariableToHost<Real>(particles, "DensityChangeRate");

        Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
        Vecd *force = particles.getVariableDataByName<Vecd>("Force");
        Real *rho = particles.getVariableDataByName<Real>("Density");
        Real *p = particles.getVariableDataByName<Real>("Pressure");
        Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        Real *mass = particles.getVariableDataByName<Real>("Mass");
        Real *drho = particles.getVariableDataByName<Real>("DensityChangeRate");
        const size_t n = particles.TotalRealParticles();

        std::vector<PairState2D> P(n);
        for (size_t i = 0; i < n; ++i)
        {
            P[i].x = double(pos[i][0]);
            P[i].y = double(pos[i][1]);
            P[i].u = double(vel[i][0]);
            P[i].v = double(vel[i][1]);
            P[i].p = double(p[i]);
            P[i].rho = double(rho[i]);
            P[i].Vol = double(Vol[i]);
            P[i].mass = double(mass[i]);
        }

        Metrics m;
        double sumFx = 0, sumFy = 0, sumAbs = 0;
        for (size_t i = 0; i < n; ++i)
        {
            sumFx += double(force[i][0]);
            sumFy += double(force[i][1]);
            sumAbs += std::abs(double(force[i][0])) + std::abs(double(force[i][1]));
            if (!isInterior(P[i].x, P[i].y, double(L), double(margin)))
                continue;
            ++m.n_interior;
            const auto ref = centralRhsPaperForm(i, P, double(h));
            if (check_force)
            {
                const double ax = double(force[i][0] / mass[i]);
                const double ay = double(force[i][1] / mass[i]);
                const double dax = std::abs(ax - ref.ax);
                const double day = std::abs(ay - ref.ay);
                m.max_a_abs = std::max(m.max_a_abs, std::max(dax, day));
                m.max_a_rel = std::max(m.max_a_rel, dax / double(a_char));
                m.max_a_rel = std::max(m.max_a_rel, day / double(a_char));
            }
            if (check_drho)
            {
                const double ddr = std::abs(double(drho[i]) - ref.drho_dt);
                m.max_drho_abs = std::max(m.max_drho_abs, ddr);
                m.max_drho_rel = std::max(m.max_drho_rel, ddr / double(drho_char));
            }
        }
        m.pair_force_rel = sumAbs > 0 ? std::sqrt(sumFx * sumFx + sumFy * sumFy) / sumAbs : 0;
        bool ok = true;
        if (check_force)
            ok = ok && (m.max_a_rel < double(tol_rel));
        if (check_drho)
            ok = ok && (m.max_drho_rel < double(tol_rel));
        m.pass = ok;
        std::cout << std::scientific << std::setprecision(3) << "[O0b] " << name << " n=" << n
                  << " n_int=" << m.n_interior << " h=" << double(h) << " max_a_rel=" << m.max_a_rel
                  << " max_a_abs=" << m.max_a_abs << " max_drho_rel=" << m.max_drho_rel
                  << " pair|ΣF|/Σ|F|=" << m.pair_force_rel << " -> " << (m.pass ? "PASS" : "FAIL") << "\n";
        return m;
    };

    struct Job
    {
        std::string name;
        FieldCase fc;
        bool force;
        bool dens;
    };
    const Job jobs[] = {
        {"const_p_zero_u_force", FieldCase::ConstPZeroU, true, false},
        {"linear_p_zero_u_force", FieldCase::LinearPZeroU, true, false},
        {"const_p_smooth_u_force", FieldCase::ConstPSmoothU, true, false},
        {"const_p_smooth_u_drho", FieldCase::ConstPSmoothU, false, true},
        {"linear_p_smooth_u_both", FieldCase::LinearPZeroU, true, true}, // dens~0; still checks
    };

    int fails = 0;
    for (const auto &job : jobs)
    {
        applyField(job.fc);
        // dt=0: Init refreshes p=EOS(ρ); Interact fills RHS; Update no state change.
        acoustic_1st.exec(Real(0));
        if (job.dens)
        {
            // 2nd half adds continuity into DensityChangeRate; overwrites Force with PJump(=0).
            acoustic_2nd.exec(Real(0));
        }
        // If both force+dens requested: force must be read before 2nd. Split such jobs.
        if (job.force && job.dens)
        {
            // Re-run: force after 1st only
            applyField(job.fc);
            acoustic_1st.exec(Real(0));
            Metrics mf = compareAfterKernels(job.name + "/force", true, false);
            if (!mf.pass)
                ++fails;
            acoustic_2nd.exec(Real(0));
            Metrics md = compareAfterKernels(job.name + "/drho", false, true);
            if (!md.pass)
                ++fails;
        }
        else
        {
            Metrics m = compareAfterKernels(job.name, job.force, job.dens);
            if (!m.pass)
                ++fails;
        }
    }

    std::ofstream ofs(out_json);
    ofs << std::setprecision(16) << "{\n  \"doc\": \"203\", \"case\": \"O0b\",\n"
        << "  \"n_particles\": " << particles.TotalRealParticles() << ",\n"
        << "  \"h\": " << double(h) << ", \"dp\": " << double(dp) << ",\n"
        << "  \"tol_rel\": " << double(tol_rel) << ", \"fails\": " << fails << ",\n"
        << "  \"solver\": \"AcousticStep NoRiemannSolverCK NoKernelCorrectionCK\"\n}\n";
    std::cout << "O0b summary fails=" << fails << " out=" << out_json << "\n";
    return fails == 0 ? 0 : 1;
}

/** Doc203 O2a: closed square, u=0, p'=0, NoRiemann+wall, no viscosity/density-reg/gravity. */
int runO2a(const std::string &out_json, int argc, char **argv)
{
    const Real L = Real(1);
    const Real dp = Real(1) / Real(20);
    const Real BW = Real(3) * dp;
    const Real rho0 = Real(1000);
    const Real c0 = Real(10);
    const Real t_cross = L / c0;
    const Real end_time = Real(10) * t_cross; // 10 acoustic crossings
    const Real U_ref_proxy = Real(1.0e-2);    // until C1 Uref; Doc203 ≤1% Uref
    const Real u_tol = Real(0.01) * U_ref_proxy;

    std::vector<char *> sph_argv;
    sph_argv.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--case" || a == "--out")
        {
            if (i + 1 < argc)
                ++i;
            continue;
        }
        sph_argv.push_back(argv[i]);
    }
    int sph_argc = static_cast<int>(sph_argv.size());

    BoundingBoxd bounds(Vec2d(-BW, -BW), Vec2d(L + BW, L + BW));
    SPHSystem sph_system(bounds, dp);
    sph_system.handleCommandlineOptions(sph_argc, sph_argv.data());

    GeometricShapeBox fluid_shape(Transform(Vec2d(0.5 * L, 0.5 * L)), Vec2d(0.5 * L, 0.5 * L), "Fluid");
    FluidBody fluid(sph_system, fluid_shape);
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, c0);
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
    BaseParticles &particles = fluid.getBaseParticles();

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cells(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cells(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> rel(fluid_inner, fluid_wall);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_pos(fluid);

    using Acoustic1stNoR =
        fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, NoRiemannSolverCK, NoKernelCorrectionCK>,
                                            Contact<Wall, NoRiemannSolverCK, NoKernelCorrectionCK>>;
    using Acoustic2ndNoR = fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK;
    InteractionDynamicsCK<MainExecutionPolicy, Acoustic1stNoR> acoustic_1st(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, Acoustic2ndNoR> acoustic_2nd(fluid_inner, fluid_wall);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> adv_dt(fluid, U_ref_proxy);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> ac_dt(fluid);

    wall_normal.exec();
    fluid_cells.exec();
    wall_cells.exec();
    rel.exec();

    // Rest IC: ρ=ρ0, u=0 (default). Force cleared each advection via setup.
    Real t = Real(0);
    int steps = 0;
    Real u_max_end = Real(0);
    Real rho_minmax_span = Real(0);
    while (t < end_time)
    {
        adv_setup.exec();
        hostZeroVec(particles, "ForcePrior");
        const Real dt_adv = adv_dt.exec();
        Real t_relax = Real(0);
        while (t_relax < dt_adv)
        {
            const Real dt = ac_dt.exec();
            acoustic_1st.exec(dt);
            acoustic_2nd.exec(dt);
            t_relax += dt;
            t += dt;
            ++steps;
        }
        update_pos.exec();
        if (steps % 50 == 0)
        {
            fluid_cells.exec();
            rel.exec();
        }
        else
        {
            fluid_cells.exec();
            rel.exec();
        }
    }

    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Density");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const size_t n = particles.TotalRealParticles();
    Real umax = Real(0);
    Real rmin = rho[0], rmax = rho[0];
    for (size_t i = 0; i < n; ++i)
    {
        umax = std::max(umax, vel[i].norm());
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
    }
    u_max_end = umax;
    rho_minmax_span = rmax - rmin;
    const bool pass = (u_max_end < u_tol);
    std::cout << std::scientific << std::setprecision(3)
              << "[O2a] n=" << n << " steps~" << steps << " t=" << double(t)
              << " Umax=" << double(u_max_end) << " tol=" << double(u_tol)
              << " dρ=" << double(rho_minmax_span)
              << " -> " << (pass ? "PASS" : "FAIL") << "\n";

    std::ofstream ofs(out_json);
    ofs << std::setprecision(16) << "{\n  \"doc\": \"203\", \"case\": \"O2a\",\n"
        << "  \"scheme\": \"NoRiemann+Wall NoKernelCorrection\",\n"
        << "  \"end_time\": " << double(end_time) << ", \"Umax\": " << double(u_max_end)
        << ", \"u_tol\": " << double(u_tol) << ", \"drho_span\": " << double(rho_minmax_span)
        << ", \"pass\": " << (pass ? "true" : "false") << "\n}\n";
    return pass ? 0 : 1;
}

/**
 * Doc203 O2b / Doc105 P1: stable strat. T=Tref+a*ζ.
 * p_h = ρ0 g (H−ζ) in BackgroundPressure (omit GravityForceCK);
 * p' = ½ ρ0 g β a (ζ²−H²) → ρ via EOS; PerturbationMomentum=1;
 * ForcePrior = Boussinesq (= −m β(T−Tref) g_vec). NoRiemann central candidate.
 */
int runO2b(const std::string &out_json, int argc, char **argv)
{
    const Real L = Real(1);
    const Real H = L;
    const Real dp = Real(1) / Real(40); // mid mesh: discrete ∇p' better balances buoyancy
    const Real BW = Real(3) * dp;
    const Real rho0 = Real(1000);
    const Real c0 = Real(50); // higher c0 → smaller WC hydro error on p'
    const Real g = Real(9.81);
    const Real beta = Real(1.0e-3);
    const Real Tref = Real(300);
    const Real a_strat = Real(10); // K/m; ΔT = a*H = 10 K, stable (hotter upward)
    const Real t_cross = L / c0;
    const Real end_time = Real(10) * t_cross;
    // Buoyancy velocity scale; Doc203 1% Uref — O2b uses U_b until C1 Uref exists.
    const Real U_b = std::sqrt(std::max(Real(0), g * beta * a_strat * H * H));
    const Real U_ref_proxy = std::max(Real(1.0e-2), U_b);
    const Real u_tol = Real(0.01) * U_ref_proxy;
    const std::string Tfield = "Temperature";
    const Real y_bottom = Real(0);

    std::vector<char *> sph_argv;
    sph_argv.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--case" || a == "--out")
        {
            if (i + 1 < argc)
                ++i;
            continue;
        }
        sph_argv.push_back(argv[i]);
    }
    int sph_argc = static_cast<int>(sph_argv.size());

    BoundingBoxd bounds(Vec2d(-BW, -BW), Vec2d(L + BW, L + BW));
    SPHSystem sph_system(bounds, dp);
    sph_system.handleCommandlineOptions(sph_argc, sph_argv.data());

    GeometricShapeBox fluid_shape(Transform(Vec2d(0.5 * L, 0.5 * H)), Vec2d(0.5 * L, 0.5 * H), "Fluid");
    FluidBody fluid(sph_system, fluid_shape);
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, c0);
    fluid.generateParticles<BaseParticles, Lattice>();

    class WallBoundary : public ComplexShape
    {
      public:
        WallBoundary(const std::string &name, Real Len, Real height, Real wall_bw) : ComplexShape(name)
        {
            add<GeometricShapeBox>(Transform(Vec2d(0.5 * Len, 0.5 * height)),
                                   Vec2d(0.5 * Len + wall_bw, 0.5 * height + wall_bw));
            subtract<GeometricShapeBox>(Transform(Vec2d(0.5 * Len, 0.5 * height)),
                                        Vec2d(0.5 * Len, 0.5 * height));
        }
    };
    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall", L, H, BW));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall(fluid, {&wall});
    BaseParticles &particles = fluid.getBaseParticles();
    particles.registerStateVariable<Real>(Tfield);
    particles.addEvolvingVariable<Real>(Tfield);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cells(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cells(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> rel(fluid_inner, fluid_wall);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_pos(fluid);

    // LinearCorrection improves discrete ∇p' vs body force near walls (Doc77a path).
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        lin_corr(DynamicsArgs(fluid_inner, Real(0.5)), fluid_wall);
    using Acoustic1stNoR =
        fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, NoRiemannSolverCK, LinearCorrectionCK>,
                                            Contact<Wall, NoRiemannSolverCK, LinearCorrectionCK>>;
    using Acoustic2ndNoR =
        fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, NoRiemannSolverCK, LinearCorrectionCK>,
                                            Contact<Wall, NoRiemannSolverCK, LinearCorrectionCK>>;
    InteractionDynamicsCK<MainExecutionPolicy, Acoustic1stNoR> acoustic_1st(fluid_inner, fluid_wall);
    InteractionDynamicsCK<MainExecutionPolicy, Acoustic2ndNoR> acoustic_2nd(fluid_inner, fluid_wall);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> adv_dt(fluid, U_ref_proxy);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> ac_dt(fluid);

    wall_normal.exec();
    fluid_cells.exec();
    wall_cells.exec();
    rel.exec();
    lin_corr.exec();
    adv_setup.exec();

    // Doc105 P1 fields (registered by AcousticStep).
    syncVariableToHost<Real>(particles, "PerturbationMomentum");
    syncVariableToHost<Real>(particles, "AcousticDensMode");
    syncVariableToHost<Vecd>(particles, "ReferenceBodyAcceleration");
    syncVariableToHost<Real>(particles, "BackgroundPressure");
    {
        Real *pm = particles.getVariableDataByName<Real>("PerturbationMomentum");
        Real *dm = particles.getVariableDataByName<Real>("AcousticDensMode");
        Vecd *aref = particles.getVariableDataByName<Vecd>("ReferenceBodyAcceleration");
        const size_t nn = particles.TotalRealParticles();
        for (size_t i = 0; i < nn; ++i)
        {
            pm[i] = Real(1);   // momentum on p'
            dm[i] = Real(3);   // omit dens p-jump (Doc201 strat); NoRiemann anyway 0
            aref[i] = Vecd::Zero();
        }
        syncVariableToDevice<Real>(particles, "PerturbationMomentum");
        syncVariableToDevice<Real>(particles, "AcousticDensMode");
        syncVariableToDevice<Vecd>(particles, "ReferenceBodyAcceleration");
    }

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, Tfield);
    syncVariableToHost<Vecd>(particles, "Velocity");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *p_bg = particles.getVariableDataByName<Real>("BackgroundPressure");
    Real *T = particles.getVariableDataByName<Real>(Tfield);
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0);
    Real mass_0 = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real zeta = std::min(H, std::max(Real(0), pos[i][1] - y_bottom));
        T[i] = Tref + a_strat * zeta;
        p_bg[i] = rho0 * g * (H - zeta); // base hydro only
        const Real p_prime = Real(0.5) * rho0 * g * beta * a_strat * (zeta * zeta - H * H);
        rho[i] = rho0 + p_prime * inv_c2;
        mass[i] = rho[i] * Vol[i]; // Doc87: Mass fixed after init
        vel[i] = Vecd::Zero();
        mass_0 += mass[i];
    }
    // Freeze Vol from initial lattice spacing for SPH volume; Mass carries ρ.
    syncVariableToDevice<Real>(particles, Tfield);
    syncVariableToDevice<Real>(particles, "Density");
    syncVariableToDevice<Real>(particles, "Mass");
    syncVariableToDevice<Real>(particles, "BackgroundPressure");
    syncVariableToDevice<Vecd>(particles, "Velocity");

    auto applyBoussinesqForcePrior = [&]()
    {
        syncVariableToHost<Real>(particles, Tfield);
        syncVariableToHost<Real>(particles, "Mass");
        syncVariableToHost<Vecd>(particles, "ForcePrior");
        Real *Th = particles.getVariableDataByName<Real>(Tfield);
        Real *m = particles.getVariableDataByName<Real>("Mass");
        Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
        const size_t nn = particles.TotalRealParticles();
        for (size_t i = 0; i < nn; ++i)
        {
            const Real factor = -m[i] * beta * (Th[i] - Tref);
            fp[i][0] = Real(0);
            fp[i][1] = factor * (-g);
        }
        syncVariableToDevice<Vecd>(particles, "ForcePrior");
    };

    Real t = Real(0);
    int steps = 0;
    Real u_max_end = Real(0);
    Real u_max_mid = Real(0);
    bool mid_sampled = false;
    Real mass_end = Real(0);
    while (t < end_time)
    {
        adv_setup.exec(); // Vol=m/ρ each advect — OK with fixed Mass
        applyBoussinesqForcePrior();
        const Real dt_adv = adv_dt.exec();
        Real t_relax = Real(0);
        while (t_relax < dt_adv)
        {
            const Real dt = ac_dt.exec();
            acoustic_1st.exec(dt);
            acoustic_2nd.exec(dt);
            t_relax += dt;
            t += dt;
            ++steps;
        }
        update_pos.exec();
        fluid_cells.exec();
        rel.exec();
        lin_corr.exec();
        if (!mid_sampled && t >= Real(0.5) * end_time)
        {
            syncVariableToHost<Vecd>(particles, "Velocity");
            Vecd *v = particles.getVariableDataByName<Vecd>("Velocity");
            Real um = Real(0);
            for (size_t i = 0; i < n; ++i)
                um = std::max(um, v[i].norm());
            u_max_mid = um;
            mid_sampled = true;
        }
    }

    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "Mass");
    vel = particles.getVariableDataByName<Vecd>("Velocity");
    rho = particles.getVariableDataByName<Real>("Density");
    mass = particles.getVariableDataByName<Real>("Mass");
    Real umax = Real(0);
    Real rmin = rho[0], rmax = rho[0];
    for (size_t i = 0; i < n; ++i)
    {
        umax = std::max(umax, vel[i].norm());
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
        mass_end += mass[i];
    }
    u_max_end = umax;
    const Real mass_rel = std::abs(mass_end - mass_0) / std::max(mass_0, Real(1e-30));
    const bool growing = (u_max_end > Real(2) * std::max(u_max_mid, Real(1e-12)) && u_max_end > u_tol);
    const bool pass = (u_max_end < u_tol) && !growing && (mass_rel < Real(1e-6));
    std::cout << std::scientific << std::setprecision(3)
              << "[O2b] n=" << n << " steps~" << steps << " t=" << double(t)
              << " Umax=" << double(u_max_end) << " Umid=" << double(u_max_mid)
              << " tol=" << double(u_tol) << " dρ=" << double(rmax - rmin)
              << " mass_rel=" << double(mass_rel)
              << " -> " << (pass ? "PASS" : "FAIL") << "\n";

    std::ofstream ofs(out_json);
    ofs << std::setprecision(16) << "{\n  \"doc\": \"203\", \"case\": \"O2b_P1\",\n"
        << "  \"scheme\": \"NoRiemann+Wall P1 bgp+p'+Boussinesq\",\n"
        << "  \"beta\": " << double(beta) << ", \"a_strat\": " << double(a_strat) << ",\n"
        << "  \"Umax\": " << double(u_max_end) << ", \"Umid\": " << double(u_max_mid)
        << ", \"u_tol\": " << double(u_tol) << ", \"mass_rel\": " << double(mass_rel)
        << ", \"pass\": " << (pass ? "true" : "false") << "\n}\n";
    return pass ? 0 : 1;
}

} // namespace

int main(int argc, char **argv)
{
    std::string which = "o0b";
    std::string out_json = "o0_report.json";
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--case" && i + 1 < argc)
            which = argv[++i];
        else if (a == "--out" && i + 1 < argc)
            out_json = argv[++i];
    }
    if (which == "o0a")
        return runO0a(out_json);
    if (which == "o0b")
        return runO0b(out_json, argc, argv);
    if (which == "o2a")
        return runO2a(out_json, argc, argv);
    if (which == "o2b")
        return runO2b(out_json, argc, argv);
    std::cerr << "Unknown --case " << which << " (o0a|o0b|o2a|o2b)\n";
    return 2;
}
