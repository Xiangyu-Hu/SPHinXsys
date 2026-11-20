/**
 * @file test_3d_ophelie_doc226_two_layer_couette.cpp
 * @brief Doc227 D: two-layer Couette on SYCL/CK GPU with ParticleViscosity.
 *
 * Periodic: host-built ghost particles with translated positions, temporarily
 * included in TotalRealParticles for CK CLL + neighbor build (classic
 * PeriodicCondition CLL injection is invisible to CK NeighborSearch).
 *
 * Bottom U=0, top U=(U,0,0); μ=μ1 below z=a, μ=μ2 above.
 */
#include "electromagnetic_ophelie_french_material_laws.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
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
Real analytic_u(Real z, Real H, Real a, Real U, Real mu1, Real mu2)
{
    const Real tau = U / (a / mu1 + (H - a) / mu2);
    if (z <= a)
        return tau * z / mu1;
    return tau * a / mu1 + tau * (z - a) / mu2;
}

Real analytic_tau(Real H, Real a, Real U, Real mu1, Real mu2)
{
    return U / (a / mu1 + (H - a) / mu2);
}

class FluidBox : public ComplexShape
{
  public:
    explicit FluidBox(const std::string &name, Real Lx, Real Ly, Real H) : ComplexShape(name)
    {
        const Vecd half(0.5 * Lx, 0.5 * Ly, 0.5 * H);
        add<GeometricShapeBox>(Transform(half), half);
    }
};

class WallPlates : public ComplexShape
{
  public:
    explicit WallPlates(const std::string &name, Real Lx, Real Ly, Real H, Real wall_thickness)
        : ComplexShape(name)
    {
        const Vecd half_xy(0.5 * Lx + wall_thickness, 0.5 * Ly + wall_thickness, 0.5 * wall_thickness);
        add<GeometricShapeBox>(Transform(Vecd(0.5 * Lx, 0.5 * Ly, -0.5 * wall_thickness)), half_xy);
        add<GeometricShapeBox>(Transform(Vecd(0.5 * Lx, 0.5 * Ly, H + 0.5 * wall_thickness)), half_xy);
    }
};

inline void hostSetTwoLayerMu(BaseParticles &particles, Real a, Real mu1, Real mu2)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        mu[i] = (pos[i][2] <= a) ? mu1 : mu2;
    syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
}

inline void hostSetWallVelocities(BaseParticles &wall, Real H, Real U)
{
    syncVariableToHost<Vecd>(wall, "Position");
    Vecd *vel = wall.registerStateVariableData<Vecd>("Velocity");
    Vecd *vel_ave = wall.registerStateVariableData<Vecd>("AverageVelocity");
    const Vecd *pos = wall.getVariableDataByName<Vecd>("Position");
    const size_t n = wall.TotalRealParticles();
    const Real z_mid = 0.5 * H;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd uw = (pos[i][2] > z_mid) ? Vecd(U, Real(0), Real(0)) : Vecd::Zero();
        vel[i] = uw;
        vel_ave[i] = uw;
    }
    syncVariableToDevice<Vecd>(wall, "Velocity");
    syncVariableToDevice<Vecd>(wall, "AverageVelocity");
}

inline void hostWrapPeriodicXY(BaseParticles &particles, Real Lx, Real Ly, UnsignedInt n_real)
{
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    for (UnsignedInt i = 0; i < n_real; ++i)
    {
        Real x = pos[i][0], y = pos[i][1];
        x -= Lx * std::floor(x / Lx);
        y -= Ly * std::floor(y / Ly);
        pos[i][0] = x;
        pos[i][1] = y;
    }
    syncVariableToDevice<Vecd>(particles, "Position");
}

/**
 * Build XY periodic ghosts for CK: copy near-boundary real particles to ghost
 * slots with ±L translation. Returns exclusive end index (n_real + n_ghost).
 * Ghosts must already be reserved in ParticlesBound.
 */
inline UnsignedInt hostFillPeriodicGhostsXY(BaseParticles &particles, Real Lx, Real Ly, Real cut,
                                             UnsignedInt n_real, UnsignedInt ghost0, UnsignedInt ghost_cap)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    UnsignedInt g = ghost0;
    const UnsignedInt ghost_end = ghost0 + ghost_cap;
    auto push = [&](UnsignedInt src, const Vecd &pimg) {
        if (g >= ghost_end)
        {
            std::cerr << "[doc227][D] FATAL ghost overflow g=" << g << " cap_end=" << ghost_end << "\n";
            std::exit(2);
        }
        particles.updateGhostParticle(g, src);
        pos[g] = pimg;
        ++g;
    };
    for (UnsignedInt i = 0; i < n_real; ++i)
    {
        const Real x = pos[i][0], y = pos[i][1];
        const bool xl = x < cut, xh = x > Lx - cut, yl = y < cut, yh = y > Ly - cut;
        if (xl)
            push(i, Vecd(x + Lx, y, pos[i][2]));
        if (xh)
            push(i, Vecd(x - Lx, y, pos[i][2]));
        if (yl)
            push(i, Vecd(x, y + Ly, pos[i][2]));
        if (yh)
            push(i, Vecd(x, y - Ly, pos[i][2]));
        // corners (both axes)
        if (xl && yl)
            push(i, Vecd(x + Lx, y + Ly, pos[i][2]));
        if (xl && yh)
            push(i, Vecd(x + Lx, y - Ly, pos[i][2]));
        if (xh && yl)
            push(i, Vecd(x - Lx, y + Ly, pos[i][2]));
        if (xh && yh)
            push(i, Vecd(x - Lx, y - Ly, pos[i][2]));
    }
    // updateGhostParticle copies host state; push fields CK interactions need on device.
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Vecd>(particles, "Velocity");
    syncVariableToDevice<Real>(particles, "Density");
    syncVariableToDevice<Real>(particles, "Mass");
    syncVariableToDevice<Real>(particles, "VolumetricMeasure");
    syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
    return g;
}

struct ProfileMetrics
{
    Real u_L2 = 0, u_L2_raw = 0, u_interface_err = 0;
    Real shear_lo_err = 0, shear_hi_err = 0;
    Real u_interface_num = 0, u_interface_ana = 0, u_max_abs = 0;
    Real uy_rms = 0, uz_rms = 0;
    size_t n_profile = 0, n_outlier = 0;
};

inline ProfileMetrics hostEvaluateProfile(BaseParticles &particles, Real H, Real a, Real U, Real mu1,
                                          Real mu2, Real h_sm, Real dp, UnsignedInt n_real)
{
    ProfileMetrics m;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real u_out = Real(2) * std::abs(U) + TinyReal;
    const Real v_out = Real(0.5) * std::abs(U) + TinyReal;
    std::vector<std::pair<Real, Real>> samples;
    Real s2 = 0, w = 0, sy2 = 0, sz2 = 0;
    for (UnsignedInt i = 0; i < n_real; ++i)
    {
        const Real z = pos[i][2];
        if (z < dp || z > H - dp)
            continue;
        const Real u = vel[i][0];
        const Real ua = analytic_u(z, H, a, U, mu1, mu2);
        samples.emplace_back(z, u);
        s2 += (u - ua) * (u - ua);
        w += ua * ua;
        sy2 += vel[i][1] * vel[i][1];
        sz2 += vel[i][2] * vel[i][2];
        ++m.n_profile;
        m.u_max_abs = std::max(m.u_max_abs, std::abs(u));
        if (std::abs(u) > u_out || std::abs(vel[i][1]) > v_out || std::abs(vel[i][2]) > v_out)
            ++m.n_outlier;
    }
    m.u_L2 = m.u_L2_raw =
        (w > TinyReal) ? std::sqrt(s2 / w) : std::sqrt(s2 / Real(std::max<size_t>(m.n_profile, 1)));
    m.uy_rms = std::sqrt(sy2 / Real(std::max<size_t>(m.n_profile, 1)));
    m.uz_rms = std::sqrt(sz2 / Real(std::max<size_t>(m.n_profile, 1)));
    Real ui_sum = 0;
    size_t ni = 0;
    for (const auto &p : samples)
    {
        if (std::abs(p.first - a) > dp)
            continue;
        ui_sum += p.second;
        ++ni;
    }
    m.u_interface_num = ni ? ui_sum / Real(ni) : Real(0);
    m.u_interface_ana = analytic_u(a, H, a, U, mu1, mu2);
    m.u_interface_err =
        std::abs(m.u_interface_num - m.u_interface_ana) / (std::abs(m.u_interface_ana) + TinyReal);
    const Real tau = analytic_tau(H, a, U, mu1, mu2);
    auto band = [&](Real zc, Real half_w) {
        Real sz = 0, su = 0, sz2b = 0, szu = 0;
        size_t c = 0;
        for (const auto &p : samples)
        {
            if (std::abs(p.first - zc) > half_w)
                continue;
            sz += p.first;
            su += p.second;
            sz2b += p.first * p.first;
            szu += p.first * p.second;
            ++c;
        }
        if (c < 4)
            return Real(0);
        const Real inv = Real(1) / Real(c);
        const Real zbar = sz * inv, ubar = su * inv;
        const Real varz = sz2b * inv - zbar * zbar;
        if (varz <= TinyReal)
            return Real(0);
        return (szu * inv - zbar * ubar) / varz;
    };
    const Real z_lo = a - Real(2) * h_sm - Real(2) * dp;
    const Real z_hi = a + Real(2) * h_sm + Real(2) * dp;
    m.shear_lo_err =
        std::abs(mu1 * band(std::max(z_lo, Real(3) * dp), Real(1.5) * dp) - tau) / (std::abs(tau) + TinyReal);
    m.shear_hi_err =
        std::abs(mu2 * band(std::min(z_hi, H - Real(3) * dp), Real(1.5) * dp) - tau) / (std::abs(tau) + TinyReal);
    return m;
}
} // namespace

int main(int argc, char *argv[])
{
    Real H = Real(0.02), Lx = H, U = Real(0.01), rho0 = Real(2750);
    Real mu1 = Real(4), mu2 = Real(40), a_frac = Real(0.5);
    int n_dp = 40;
    Real end_time = Real(0.2), c0 = Real(10);
    Real acoustic_cfl = Real(0.55), advection_cfl = Real(0.25);
    std::string summary_json;
    bool degenerate = false, doc227_periodic_precheck = false, doc227_wrap_only = false;
    // DensitySummationCK on this lattice is ~+2% high → acoustic blow-up with c0=10.
    // Default: continuity-only (skip summation); --doc227-density-summation to enable.
    bool use_density_summation = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::strncmp(argv[i], "--dp-div=", 9) == 0)
            n_dp = std::atoi(argv[i] + 9);
        else if (std::strncmp(argv[i], "--end-time=", 11) == 0)
            end_time = std::strtod(argv[i] + 11, nullptr);
        else if (std::strncmp(argv[i], "--mu2=", 6) == 0)
            mu2 = std::strtod(argv[i] + 6, nullptr);
        else if (std::strcmp(argv[i], "--degenerate-equal-mu") == 0)
            degenerate = true;
        else if (std::strncmp(argv[i], "--Lx-factor=", 12) == 0)
            Lx = H * std::strtod(argv[i] + 12, nullptr);
        else if (std::strncmp(argv[i], "--summary-json=", 15) == 0)
            summary_json = argv[i] + 15;
        else if (std::strncmp(argv[i], "--c0=", 5) == 0)
            c0 = std::strtod(argv[i] + 5, nullptr);
        else if (std::strncmp(argv[i], "--acoustic-cfl=", 15) == 0)
            acoustic_cfl = std::strtod(argv[i] + 15, nullptr);
        else if (std::strncmp(argv[i], "--advection-cfl=", 16) == 0)
            advection_cfl = std::strtod(argv[i] + 16, nullptr);
        else if (std::strcmp(argv[i], "--doc227-periodic-precheck") == 0)
            doc227_periodic_precheck = true;
        else if (std::strcmp(argv[i], "--doc227-wrap-only") == 0)
            doc227_wrap_only = true;
        else if (std::strcmp(argv[i], "--doc227-density-summation") == 0)
            use_density_summation = true;
    }
    if (degenerate)
        mu2 = mu1;

    const Real dp = H / Real(n_dp);
    const Real Ly = Real(8) * dp;
    const Real a = a_frac * H;
    const Real wall_thickness = Real(3) * dp;

    // XY margin must cover periodic ghosts (±cut ≈ 1.5 h); walls need Z margin only.
    const Real xy_margin = Real(6) * dp;
    BoundingBoxd bounds(Vecd(-xy_margin, -xy_margin, -Real(2) * wall_thickness),
                        Vecd(Lx + xy_margin, Ly + xy_margin, H + Real(2) * wall_thickness));
    SPHSystem sph_system(bounds, dp);

    // Ghost pool ~2× real particles (XY faces + corners).
    Ghost<ReserveSizeFactor> ghost_reserve(Real(2.0));
    FluidBody fluid(sph_system, makeShared<FluidBox>("Fluid", Lx, Ly, H));
    fluid.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, c0);
    fluid.addMaterialProperty<ParticleViscosity>(mu1);
    fluid.generateParticlesWithReserve<BaseParticles, Lattice>(ghost_reserve);
    ParticleViscosity &particle_mu = fluid.getMaterialProperty<ParticleViscosity>();
    Viscosity &visc_ref = fluid.getMaterialProperty<Viscosity>();

    SolidBody wall(sph_system, makeShared<WallPlates>("Wall", Lx, Ly, H, wall_thickness));
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall(fluid, {&wall});

    BaseParticles &fp = fluid.getBaseParticles();
    const UnsignedInt n_real = fp.TotalRealParticles();
    const UnsignedInt ghost0 = n_real;
    const UnsignedInt ghost_cap = fp.ParticlesBound() - ghost0;

    fp.registerStateVariable<Real>(kOphelieDynamicViscosityField, mu1);
    fp.addEvolvingVariable<Real>(kOphelieDynamicViscosityField);
    particle_mu.bindParticleViscosityField(fp, kOphelieDynamicViscosityField);
    hostSetTwoLayerMu(fp, a, mu1, mu2);
    visc_ref.setReferenceViscosity(std::max(mu1, mu2));
    {
        syncVariableToHost<Real>(fp, "Density");
        Real *rho = fp.getVariableDataByName<Real>("Density");
        for (UnsignedInt i = 0; i < n_real; ++i)
            rho[i] = rho0;
        syncVariableToDevice<Real>(fp, "Density");
    }

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(execution::par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(execution::par_host);
    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    // NormalFromBodyShapeCK is HostKernel-only (not device).
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall).exec();
    auto &fluid_cll = main_methods.addCellLinkedListDynamics(fluid);
    auto &wall_cll = main_methods.addCellLinkedListDynamics(wall);
    auto &update_rel = main_methods.addRelationDynamics(fluid_inner, fluid_wall);

    auto &adv_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(fluid);
    auto &update_pos = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(fluid);
    auto &ac1 =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(fluid_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_wall);
    auto &ac2 =
        main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK,
                                                    NoKernelCorrectionCK>(fluid_inner);
    auto &ac2_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(fluid_wall);
    ac2.addPostContactInteraction(ac2_wall);
    auto &dens =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(fluid_inner)
            .addPostContactInteraction(fluid_wall)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, Internal>(fluid);

    auto *fluid_viscous = &main_methods.addInteractionDynamicsWithUpdate<
        fluid_dynamics::ViscousForceCK, ParticleViscosity, NoKernelCorrectionCK>(fluid_inner);
    auto *fluid_viscous_wall = &main_methods.addInteractionDynamics<
        fluid_dynamics::ViscousForceCK, Wall, ParticleViscosity, NoKernelCorrectionCK>(fluid_wall);
    fluid_viscous->addPostContactInteraction(*fluid_viscous_wall);

    auto &advection_viscous_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionViscousTimeStepCK>(fluid, U, advection_cfl);
    auto &acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(fluid, acoustic_cfl);

    const Real h_actual = fluid.getSPHAdaptation().ReferenceSmoothingLength();
    const Real cut = Real(1.5) * h_actual;

    auto refresh_periodic_config = [&]() {
        hostWrapPeriodicXY(fp, Lx, Ly, n_real);
        UnsignedInt n_active = n_real;
        if (!doc227_wrap_only)
            n_active = hostFillPeriodicGhostsXY(fp, Lx, Ly, cut, n_real, ghost0, ghost_cap);
        // CK CLL/neighbors see ghosts as temporary "real" particles.
        fp.svTotalRealParticles()->setValue(n_active);
        fluid_cll.exec();
        wall_cll.exec();
        update_rel.exec();
        fp.svTotalRealParticles()->setValue(n_real);
    };

    hostSetWallVelocities(wall.getBaseParticles(), H, U);
    refresh_periodic_config();

    std::cout << "[doc227][D] Couette CK/GPU ParticleViscosity H=" << H << " Lx=" << Lx << " Ly=" << Ly
              << " dp=" << dp << " a=" << a << " mu1=" << mu1 << " mu2=" << mu2 << " U=" << U
              << " n_real=" << n_real << " ghost_cap=" << ghost_cap << " h=" << h_actual
              << " wrap_only=" << (doc227_wrap_only ? 1 : 0)
              << " dens_sum=" << (use_density_summation ? 1 : 0) << " end_time=" << end_time << std::endl;
    std::cout << "[doc227][D] periodic=CK_host_ghosts_XY (classic PeriodicCondition not used)\n"
              << std::flush;

    if (doc227_periodic_precheck)
    {
        syncVariableToHost<Vecd>(fp, "Position");
        syncVariableToHost<Vecd>(fp, "Velocity");
        Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        for (UnsignedInt i = 0; i < n_real; ++i)
            vel[i] = Vecd(U, Real(0), Real(0));
        syncVariableToDevice<Vecd>(fp, "Velocity");
        refresh_periodic_config();

        UnsignedInt *particle_offset =
            fluid_inner.dvParticleOffset()->DelegatedData(execution::par_host);
        // After relation build on device, sync offset to host for counting.
        fluid_inner.dvParticleOffset()->prepareForOutput(execution::par_ck);
        particle_offset = fluid_inner.dvParticleOffset()->Data();

        Real mean_int = 0, mean_edge_x = 0;
        size_t n_int = 0, n_edge_x = 0;
        for (UnsignedInt i = 0; i < n_real; ++i)
        {
            const Real x = pos[i][0];
            const bool near_x = (x < cut) || (x > Lx - cut);
            const size_t nn = particle_offset[i + 1] - particle_offset[i];
            if (near_x)
            {
                mean_edge_x += Real(nn);
                ++n_edge_x;
            }
            else
            {
                mean_int += Real(nn);
                ++n_int;
            }
        }
        mean_int = n_int ? mean_int / Real(n_int) : Real(0);
        mean_edge_x = n_edge_x ? mean_edge_x / Real(n_edge_x) : Real(0);
        const Real ratio = mean_edge_x / (mean_int + TinyReal);
        const bool pass = (!doc227_wrap_only) && (ratio >= Real(0.85)) && (mean_int > Real(5));
        std::cout << "[doc227][D2] precheck wrap_only=" << (doc227_wrap_only ? 1 : 0)
                  << " mean_nbr_interior=" << mean_int << " mean_nbr_edge_x=" << mean_edge_x
                  << " edge_x/int=" << ratio << " PERIODIC_NEIGHBORS_VERIFIED="
                  << (pass ? "PASS" : "FAIL") << std::endl
                  << std::flush;
        if (!summary_json.empty())
        {
            std::ofstream js(summary_json);
            js << std::scientific << std::setprecision(8);
            js << "{\n  \"case\": \"doc227_periodic_precheck\",\n"
               << "  \"path\": \"ck_gpu_host_ghosts\",\n"
               << "  \"wrap_only\": " << (doc227_wrap_only ? 1 : 0) << ",\n"
               << "  \"mean_nbr_interior\": " << mean_int << ",\n"
               << "  \"mean_nbr_edge_x\": " << mean_edge_x << ",\n"
               << "  \"edge_x_over_interior\": " << ratio << ",\n"
               << "  \"PERIODIC_NEIGHBORS_VERIFIED\": \"" << (pass ? "PASS" : "FAIL") << "\",\n"
               << "  \"pass\": " << (pass ? 1 : 0) << "\n}\n";
        }
        return pass ? 0 : 1;
    }

    // IC analytic
    {
        syncVariableToHost<Vecd>(fp, "Position");
        syncVariableToHost<Vecd>(fp, "Velocity");
        const Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        for (UnsignedInt i = 0; i < n_real; ++i)
            vel[i] = Vecd(analytic_u(pos[i][2], H, a, U, mu1, mu2), Real(0), Real(0));
        syncVariableToDevice<Vecd>(fp, "Velocity");
    }

    {
        const auto ic = hostEvaluateProfile(fp, H, a, U, mu1, mu2, h_actual, dp, n_real);
        std::cout << "[doc227][D] IC u_L2=" << ic.u_L2 << " u_max=" << ic.u_max_abs
                  << " n_out=" << ic.n_outlier << "\n"
                  << std::flush;
    }
    std::cout << "[doc227][D] dt_adv0=" << advection_viscous_time_step.exec()
              << " entering time loop\n"
              << std::flush;

    auto hostRhoStats = [&](const char *tag) {
        syncVariableToHost<Real>(fp, "Density");
        const Real *rho = fp.getVariableDataByName<Real>("Density");
        Real rmin = rho[0], rmax = rho[0], rsum = 0;
        for (UnsignedInt i = 0; i < n_real; ++i)
        {
            rmin = std::min(rmin, rho[i]);
            rmax = std::max(rmax, rho[i]);
            rsum += rho[i];
        }
        std::cout << "[doc227][D] " << tag << " rho min/avg/max=" << rmin << "/"
                  << (rsum / Real(n_real)) << "/" << rmax << std::endl
                  << std::flush;
    };

    size_t adv_steps = 0;
    int steps = 0;
    while (!time_stepper.isEndTime(end_time))
    {
        ++adv_steps;
        if (use_density_summation)
        {
            dens.exec();
            if (adv_steps <= 2)
                hostRhoStats("after_dens");
            // Ghosts were built with pre-dens state; refresh matching images.
            refresh_periodic_config();
        }
        adv_setup.exec();
        hostSetWallVelocities(wall.getBaseParticles(), H, U);
        fluid_viscous->exec();
        if (adv_steps <= 2)
        {
            syncVariableToHost<Vecd>(fp, "ForcePrior");
            const Vecd *fpv = fp.getVariableDataByName<Vecd>("ForcePrior");
            Real fmax = 0;
            for (UnsignedInt i = 0; i < n_real; ++i)
                fmax = std::max(fmax, fpv[i].norm());
            std::cout << "[doc227][D] after_viscous ForcePrior_max=" << fmax << std::endl
                      << std::flush;
            hostRhoStats("before_acoustic");
        }

        const Real Dt = advection_viscous_time_step.exec();
        if (!(Dt > TinyReal) || !std::isfinite(Dt))
        {
            std::cerr << "[doc227][D] FATAL bad Dt=" << Dt << "\n";
            return 2;
        }
        Real integrated = 0, dt_last = 0;
        int sub = 0;
        while (integrated + TinyReal < Dt && !time_stepper.isEndTime(end_time))
        {
            const Real dt_ac = acoustic_time_step.exec();
            Real dt = std::min(dt_ac, Dt - integrated);
            if (!(dt > Real(1.0e-12)))
            {
                std::cerr << "[doc227][D] FATAL collapsed dt\n";
                return 2;
            }
            time_stepper.incrementPhysicalTime(dt);
            ac1.exec(dt);
            ac2.exec(dt);
            integrated += dt;
            dt_last = dt;
            ++steps;
            ++sub;
        }
        update_pos.exec();
        // No particle_sort — previously associated with long CPU stalls at fine dp.
        refresh_periodic_config();
        hostSetTwoLayerMu(fp, a, mu1, mu2);

        if (adv_steps <= 5 || adv_steps % 20 == 0)
        {
            const auto mid = hostEvaluateProfile(fp, H, a, U, mu1, mu2, h_actual, dp, n_real);
            std::cout << "[doc227][D] t=" << time_stepper.getPhysicalTime() << " adv=" << adv_steps
                      << " steps=" << steps << " Dt=" << Dt << " dt_ac=" << dt_last
                      << " u_L2=" << mid.u_L2 << " u_int_err=" << mid.u_interface_err
                      << " n_out=" << mid.n_outlier << " u_max=" << mid.u_max_abs << std::endl
                      << std::flush;
        }
    }

    const auto met = hostEvaluateProfile(fp, H, a, U, mu1, mu2, h_actual, dp, n_real);
    const bool pass_L2 = met.u_L2 <= Real(0.05);
    const bool pass_ui = met.u_interface_err <= Real(0.05);
    const bool shear_ok = met.n_profile < 50 ||
                          (met.shear_lo_err <= Real(0.10) && met.shear_hi_err <= Real(0.10));
    const bool state_ok = met.n_outlier * 20 <= met.n_profile;
    const bool pass = pass_L2 && pass_ui && shear_ok && state_ok && !doc227_wrap_only;

    std::cout << "[doc227][D] SUMMARY n_profile=" << met.n_profile << " n_outlier=" << met.n_outlier
              << " u_L2=" << met.u_L2 << " u_int_err=" << met.u_interface_err
              << " shear_lo_err=" << met.shear_lo_err << " shear_hi_err=" << met.shear_hi_err
              << " pass=" << (pass ? 1 : 0) << std::endl
              << std::flush;

    if (!summary_json.empty())
    {
        std::ofstream js(summary_json);
        js << std::scientific << std::setprecision(8);
        js << "{\n  \"case\": \"doc227_two_layer_couette\",\n"
           << "  \"path\": \"ck_gpu_ParticleViscosity_host_ghosts\",\n"
           << "  \"H\": " << H << ", \"Lx\": " << Lx << ", \"dp\": " << dp << ", \"a\": " << a << ",\n"
           << "  \"mu1\": " << mu1 << ", \"mu2\": " << mu2 << ", \"U\": " << U << ",\n"
           << "  \"degenerate_equal_mu\": " << (degenerate ? 1 : 0) << ",\n"
           << "  \"wrap_only\": " << (doc227_wrap_only ? 1 : 0) << ",\n"
           << "  \"n_profile\": " << met.n_profile << ", \"n_outlier\": " << met.n_outlier << ",\n"
           << "  \"u_L2\": " << met.u_L2 << ", \"u_interface_err\": " << met.u_interface_err << ",\n"
           << "  \"shear_lo_err\": " << met.shear_lo_err << ", \"shear_hi_err\": " << met.shear_hi_err
           << ",\n"
           << "  \"VARIABLE_MU_OPERATOR_VERIFIED\": \"" << (pass ? "PASS" : "FAIL") << "\",\n"
           << "  \"pass\": " << (pass ? 1 : 0) << "\n}\n";
    }
    return pass ? 0 : 1;
}
