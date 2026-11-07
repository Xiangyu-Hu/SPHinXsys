/**
 * @file test_3d_ophelie_minimal_cylinder_convection.cpp
 * @brief Doc205: minimal closed-cylinder convection wiring (C1-like assembly).
 *
 * Fixed T(r)=Tref+40*((r/R)^2-1/2); Boussinesq only; no EM/heat equation.
 * Cases: A0/A1 (all no-slip on closed wall); B0/B1 (French dual-wall: acoustic closed +
 * viscous open-cup omit lid + AcousticWallFreeSlip on closed lid).
 */
#include "electromagnetic_ophelie_boussinesq.h"
#include "sphinxsys.h"

#include <cmath>
#include <filesystem>
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

Transform axisZCylinderTransform(const Vec3d &center)
{
    const Rotation3d rot(Eigen::Quaternion<Real>::FromTwoVectors(Vec3d::UnitX(), Vec3d::UnitZ()));
    return Transform(rot, center);
}

struct Cli
{
    std::string case_id = "A0"; // A0 A1 B0 B1
    Real dp = 0.015;
    Real end_time = 0.2;
    Real R = 0.25;
    Real H = 0.185;
    Real rho0 = 2750.0;
    Real mu = 4.0;
    Real beta = 0.0; // overridden by case
    Real t_ref = 1473.0;
    Real dT_radial = 40.0;
    Real c0 = 20.0;
    Real g = 9.81;
    bool use_bg_hydro = true; // C1-like P1 pending wall check; logged
    std::string out_dir = "docs/ophelie/reference/doc205_runs/scratch";
    std::string summary_json;
};

bool parseCli(int argc, char **argv, Cli &cli, std::vector<char *> &sph_argv)
{
    sph_argv.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto take = [&](const std::string &key, Real &dst) {
            const std::string p = key + "=";
            if (a.rfind(p, 0) == 0)
            {
                dst = static_cast<Real>(std::atof(a.c_str() + p.size()));
                return true;
            }
            return false;
        };
        auto take_s = [&](const std::string &key, std::string &dst) {
            const std::string p = key + "=";
            if (a.rfind(p, 0) == 0)
            {
                dst = a.substr(p.size());
                return true;
            }
            return false;
        };
        if (take_s("--case", cli.case_id) || take("--dp", cli.dp) || take("--end-time", cli.end_time) ||
            take("--R", cli.R) || take("--H", cli.H) || take("--rho0", cli.rho0) || take("--mu", cli.mu) ||
            take("--beta", cli.beta) || take("--t-ref", cli.t_ref) || take("--dT-radial", cli.dT_radial) ||
            take("--c0", cli.c0) || take("--g", cli.g) || take_s("--out-dir", cli.out_dir) ||
            take_s("--summary-json", cli.summary_json))
            continue;
        if (a == "--no-bg-hydro")
        {
            cli.use_bg_hydro = false;
            continue;
        }
        if (a == "--help" || a == "-h")
        {
            std::cout << "Doc205 minimal cylinder\n"
                      << "  --case=A0|A1|B0|B1 --dp= --end-time= --out-dir=\n";
            return false;
        }
        sph_argv.push_back(argv[i]);
    }
    if (cli.case_id == "A0" || cli.case_id == "B0")
        cli.beta = 0.0;
    else if (cli.case_id == "A1" || cli.case_id == "B1")
    {
        if (!(cli.beta > 0))
            cli.beta = Real(1.0e-5);
    }
    else
    {
        std::cerr << "[doc205] unknown --case=" << cli.case_id << "\n";
        return false;
    }
    if (cli.out_dir == "docs/ophelie/reference/doc205_runs/scratch" && !cli.summary_json.empty())
    {
        const auto p = std::filesystem::path(cli.summary_json).parent_path();
        if (!p.empty())
            cli.out_dir = p.string();
    }
    return true;
}

class FluidCylinder : public ComplexShape
{
  public:
    FluidCylinder(const std::string &name, Real R, Real H) : ComplexShape(name)
    {
        const Vec3d c(0, 0, Real(0.5) * H);
        add<GeometricShapeCylinder>(axisZCylinderTransform(c), R, Real(0.5) * H);
    }
};

class ClosedCylinderWall : public ComplexShape
{
  public:
    ClosedCylinderWall(const std::string &name, Real R, Real H, Real bw) : ComplexShape(name)
    {
        const Vec3d c(0, 0, Real(0.5) * H);
        add<GeometricShapeCylinder>(axisZCylinderTransform(c), R + bw, Real(0.5) * H + bw);
        subtract<GeometricShapeCylinder>(axisZCylinderTransform(c), R, Real(0.5) * H);
    }
};

/**
 * Side+bottom only (no lid) — Doc205 B viscous contact.
 * Matches French open-cup: outer cylinder ends at z=H (flush rim), no particles above melt.
 */
class OpenCupCylinderWall : public ComplexShape
{
  public:
    OpenCupCylinderWall(const std::string &name, Real R, Real H, Real bw) : ComplexShape(name)
    {
        const Vec3d c(0, 0, Real(0.5) * H);
        const Vec3d outer_c(0, 0, Real(0.5) * H - Real(0.5) * bw);
        const Real outer_half = Real(0.5) * H + Real(0.5) * bw;
        add<GeometricShapeCylinder>(axisZCylinderTransform(outer_c), R + bw, outer_half);
        subtract<GeometricShapeCylinder>(axisZCylinderTransform(c), R, Real(0.5) * H);
    }
};

inline Real TofR(Real r, Real R, Real Tref, Real dT)
{
    const Real rn = r / (R + TinyReal);
    return Tref + dT * (rn * rn - Real(0.5));
}

inline void applyFixedT(BaseParticles &p, Real R, Real Tref, Real dT)
{
    syncVariableToHost<Vecd>(p, "Position");
    syncVariableToHost<Real>(p, kT);
    Vecd *pos = p.getVariableDataByName<Vecd>("Position");
    Real *T = p.getVariableDataByName<Real>(kT);
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        const Real r = std::hypot(pos[i][0], pos[i][1]);
        T[i] = TofR(r, R, Tref, dT);
    }
    syncVariableToDevice<Real>(p, kT);
}

inline void writeBgHydro(BaseParticles &p, Real rho0, Real g, Real H)
{
    syncVariableToHost<Vecd>(p, "Position");
    syncVariableToHost<Real>(p, "BackgroundPressure");
    Vecd *pos = p.getVariableDataByName<Vecd>("Position");
    Real *pbg = p.getVariableDataByName<Real>("BackgroundPressure");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        const Real z = std::min(H, std::max(Real(0), pos[i][2]));
        pbg[i] = rho0 * g * (H - z);
    }
    syncVariableToDevice<Real>(p, "BackgroundPressure");
}

/** Doc205: analytic wall normals — avoid NormalFromBodyShapeCK hang on cylinders. */
inline void setAnalyticCylinderWallNormals(BaseParticles &wall, Real R, Real H)
{
    syncVariableToHost<Vecd>(wall, "Position");
    syncVariableToHost<Vecd>(wall, "NormalDirection");
    Vecd *pos = wall.getVariableDataByName<Vecd>("Position");
    Vecd *n = wall.getVariableDataByName<Vecd>("NormalDirection");
    const Real corner = Real(0.75) * (H > R ? R : H); // face ownership band ~O(dp) set by caller scale
    (void)corner;
    for (size_t i = 0; i < wall.TotalRealParticles(); ++i)
    {
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        const Real r = std::hypot(x, y);
        const Real dist_side = r - R;     // >0 outside fluid
        const Real dist_bot = -z;         // >0 below z=0
        const Real dist_top = z - H;      // >0 above z=H
        // Pick dominant exterior face; normal points toward fluid.
        if (dist_top >= dist_bot && dist_top >= dist_side)
            n[i] = Vecd(0, 0, -1);
        else if (dist_bot >= dist_side)
            n[i] = Vecd(0, 0, 1);
        else if (r > TinyReal)
            n[i] = Vecd(-x / r, -y / r, 0); // -ê_r
        else
            n[i] = Vecd(0, 0, 1);
    }
    syncVariableToDevice<Vecd>(wall, "NormalDirection");
}

inline void stage(const char *msg)
{
    std::cout << "[doc205][stage] " << msg << std::endl;
    std::cout.flush();
}

/** Tag lid on closed wall for acoustic free-slip dens ghost (French Doc201). */
inline size_t tagAcousticWallFreeSlipLid(BaseParticles &wall, Real H, Real bw)
{
    syncVariableToHost<Vecd>(wall, "Position");
    syncVariableToHost<Vecd>(wall, "NormalDirection");
    syncVariableToHost<Real>(wall, "AcousticWallFreeSlip");
    const Vecd *pos = wall.getVariableDataByName<Vecd>("Position");
    const Vecd *n = wall.getVariableDataByName<Vecd>("NormalDirection");
    Real *fs = wall.getVariableDataByName<Real>("AcousticWallFreeSlip");
    const Real z_cut = H - Real(0.25) * bw;
    size_t n_lid = 0;
    for (size_t i = 0; i < wall.TotalRealParticles(); ++i)
    {
        const bool is_lid = (pos[i][2] > z_cut) && (std::abs(n[i][2]) > Real(0.5));
        fs[i] = is_lid ? Real(1) : Real(0);
        if (is_lid)
            ++n_lid;
    }
    syncVariableToDevice<Real>(wall, "AcousticWallFreeSlip");
    return n_lid;
}

} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    std::vector<char *> sph_argv;
    if (!parseCli(argc, argv, cli, sph_argv))
        return 1;

    const Real dp = cli.dp;
    const Real BW = Real(3) * dp;
    const Real R = cli.R, H = cli.H;
    const bool lid_slip = (cli.case_id == "B0" || cli.case_id == "B1");

    BoundingBoxd bounds(Vec3d(-R - BW, -R - BW, -BW), Vec3d(R + BW, R + BW, H + BW));
    SPHSystem sph_system(bounds, dp);
    int sph_argc = static_cast<int>(sph_argv.size());
    sph_system.handleCommandlineOptions(sph_argc, sph_argv.data());

    stage("generate fluid lattice");
    FluidBody fluid(sph_system, makeShared<FluidCylinder>("Fluid", R, H));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(cli.rho0, cli.c0);
    fluid.addMaterialProperty<Viscosity>(cli.mu);
    fluid.generateParticles<BaseParticles, Lattice>();
    stage(("fluid n=" + std::to_string(fluid.getBaseParticles().TotalRealParticles())).c_str());

    stage("generate wall_closed lattice (acoustic + dens)");
    SolidBody wall_closed(sph_system, makeShared<ClosedCylinderWall>("WallClosed", R, H, BW));
    wall_closed.defineMatterMaterial<Solid>();
    wall_closed.generateParticles<BaseParticles, Lattice>();
    stage(("wall_closed n=" + std::to_string(wall_closed.getBaseParticles().TotalRealParticles())).c_str());

    // Viscous open-cup always generated; A uses closed for viscous, B uses noslip (omit lid).
    stage("generate wall_noslip lattice (open cup)");
    SolidBody wall_noslip(sph_system, makeShared<OpenCupCylinderWall>("WallNoSlip", R, H, BW));
    wall_noslip.defineMatterMaterial<Solid>();
    wall_noslip.generateParticles<BaseParticles, Lattice>();
    stage(("wall_noslip n=" + std::to_string(wall_noslip.getBaseParticles().TotalRealParticles())).c_str());

    Inner<> fluid_inner(fluid);
    Contact<> fluid_closed(fluid, {&wall_closed});
    Contact<> fluid_noslip(fluid, {&wall_noslip});
    Contact<> &fluid_viscous = lid_slip ? fluid_noslip : fluid_closed;
    BaseParticles &fp = fluid.getBaseParticles();
    BaseParticles &wp_closed = wall_closed.getBaseParticles();
    BaseParticles &wp_noslip = wall_noslip.getBaseParticles();

    fp.registerStateVariable<Real>(kT, cli.t_ref);
    fp.addEvolvingVariable<Real>(kT);
    fp.registerStateVariable<Real>("DensSumSigma0");
    fp.addEvolvingVariable<Real>("DensSumSigma0");

    stage("construct dynamics (acoustic/viscous/dens/boussinesq)");
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cells(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_closed_cells(wall_closed);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_noslip_cells(wall_noslip);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>, Contact<>> rel(fluid_inner, fluid_closed, fluid_noslip);
    // Construct once to register NormalDirection; DO NOT exec (hangs on cylinder shapes).
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_closed_normal_register(wall_closed);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_noslip_normal_register(wall_noslip);
    (void)wall_closed_normal_register;
    (void)wall_noslip_normal_register;
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> adv_setup(fluid);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> update_pos(fluid);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK> acoustic_1st(
        fluid_inner, fluid_closed);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK> acoustic_2nd(
        fluid_inner, fluid_closed);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK> viscous(fluid_inner,
                                                                                              fluid_viscous);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>, Contact<>>> density_sum(
        fluid_inner, fluid_closed);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, Internal>> density_reg(fluid);

    Gravity gravity(Vecd(0, 0, -cli.g));
    StateDynamics<MainExecutionPolicy, BoussinesqBuoyancyForceCK> boussinesq(fluid, gravity, cli.beta, cli.t_ref, kT);

    const Real U_dt = Real(0.05);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> adv_dt(fluid, U_dt);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionViscousTimeStepCK> adv_visc_dt(fluid, U_dt);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> ac_dt(fluid);

    stage("analytic wall normals");
    setAnalyticCylinderWallNormals(wp_closed, R, H);
    setAnalyticCylinderWallNormals(wp_noslip, R, H);
    size_t n_lid_fs = 0;
    if (lid_slip)
    {
        n_lid_fs = tagAcousticWallFreeSlipLid(wp_closed, H, BW);
        stage(("AcousticWallFreeSlip lid tags=" + std::to_string(n_lid_fs)).c_str());
    }
    stage("cell lists");
    fluid_cells.exec();
    wall_closed_cells.exec();
    wall_noslip_cells.exec();
    stage("neighbor relations");
    rel.exec();
    stage("relations done");

    applyFixedT(fp, R, cli.t_ref, cli.dT_radial);

    // P1 flags (C1-like): perturbation momentum + dens mode
    syncVariableToHost<Real>(fp, "PerturbationMomentum");
    syncVariableToHost<Real>(fp, "AcousticDensMode");
    {
        Real *pm = fp.getVariableDataByName<Real>("PerturbationMomentum");
        Real *dm = fp.getVariableDataByName<Real>("AcousticDensMode");
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            pm[i] = Real(1);
            dm[i] = Real(1);
        }
        syncVariableToDevice<Real>(fp, "PerturbationMomentum");
        syncVariableToDevice<Real>(fp, "AcousticDensMode");
    }
    if (cli.use_bg_hydro)
        writeBgHydro(fp, cli.rho0, cli.g, H);
    else
    {
        syncVariableToHost<Real>(fp, "BackgroundPressure");
        Real *pbg = fp.getVariableDataByName<Real>("BackgroundPressure");
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            pbg[i] = Real(0);
        syncVariableToDevice<Real>(fp, "BackgroundPressure");
    }

    // zero velocity / fix mass
    syncVariableToHost<Vecd>(fp, "Velocity");
    syncVariableToHost<Real>(fp, "Density");
    syncVariableToHost<Real>(fp, "Mass");
    syncVariableToHost<Real>(fp, "VolumetricMeasure");
    {
        Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        Real *rho = fp.getVariableDataByName<Real>("Density");
        Real *mass = fp.getVariableDataByName<Real>("Mass");
        const Real *vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            vel[i] = Vecd::Zero();
            rho[i] = cli.rho0;
            mass[i] = cli.rho0 * vol[i];
        }
        syncVariableToDevice<Vecd>(fp, "Velocity");
        syncVariableToDevice<Real>(fp, "Density");
        syncVariableToDevice<Real>(fp, "Mass");
    }

    SingleVariable<Real> *sv_t = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    std::cout << std::scientific << std::setprecision(4)
              << "[doc205] case=" << cli.case_id << " n_fluid=" << fp.TotalRealParticles()
              << " n_wall_closed=" << wp_closed.TotalRealParticles()
              << " n_wall_noslip=" << wp_noslip.TotalRealParticles() << " lid_slip=" << (lid_slip ? 1 : 0)
              << " n_lid_fs=" << n_lid_fs << " dp=" << double(dp) << " R=" << double(R)
              << " H=" << double(H) << " beta=" << double(cli.beta) << " mu=" << double(cli.mu)
              << " c0=" << double(cli.c0) << " bg_hydro=" << (cli.use_bg_hydro ? 1 : 0)
              << " end_time=" << double(cli.end_time) << "\n";

    bool dens_sigma0_ready = false;
    int steps = 0;
    Real Umax_end = 0;
    while (sv_t->getValue() < cli.end_time)
    {
        density_sum.exec();
        if (!dens_sigma0_ready)
        {
            syncVariableToHost<Real>(fp, "DensitySummation");
            syncVariableToHost<Real>(fp, "DensSumSigma0");
            const Real *rho_sum = fp.getVariableDataByName<Real>("DensitySummation");
            Real *sigma0 = fp.getVariableDataByName<Real>("DensSumSigma0");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                sigma0[i] = std::max(rho_sum[i], TinyReal);
            syncVariableToDevice<Real>(fp, "DensSumSigma0");
            dens_sigma0_ready = true;
        }
        density_reg.exec();
        {
            syncVariableToHost<Real>(fp, "Density");
            syncVariableToHost<Real>(fp, "DensitySummation");
            syncVariableToHost<Real>(fp, "DensSumSigma0");
            Real *rho = fp.getVariableDataByName<Real>("Density");
            const Real *rho_sum = fp.getVariableDataByName<Real>("DensitySummation");
            const Real *sigma0 = fp.getVariableDataByName<Real>("DensSumSigma0");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                rho[i] = cli.rho0 * (rho_sum[i] / (sigma0[i] + TinyReal));
            syncVariableToDevice<Real>(fp, "Density");
        }
        if (cli.use_bg_hydro)
            writeBgHydro(fp, cli.rho0, cli.g, H);

        adv_setup.exec();
        // zero Force then viscous (ForcePrior) + boussinesq
        {
            syncVariableToHost<Vecd>(fp, "Force");
            Vecd *F = fp.getVariableDataByName<Vecd>("Force");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
                F[i] = Vecd::Zero();
            syncVariableToDevice<Vecd>(fp, "Force");
        }
        applyFixedT(fp, R, cli.t_ref, cli.dT_radial);
        viscous.exec();
        boussinesq.exec();

        const Real dt_adv = std::min(adv_dt.exec(), adv_visc_dt.exec());
        Real relax = 0;
        size_t ac_sub = 0;
        while (relax < dt_adv && ac_sub < 2000)
        {
            Real dt = ac_dt.exec();
            const Real left = dt_adv - relax;
            if (dt > left)
                dt = left;
            if (!(dt > TinyReal))
                break;
            acoustic_1st.exec(dt);
            acoustic_2nd.exec(dt);
            boussinesq.exec();
            relax += dt;
            sv_t->incrementValue(dt);
            ++steps;
            ++ac_sub;
        }
        update_pos.exec();
        applyFixedT(fp, R, cli.t_ref, cli.dT_radial);
        fluid_cells.exec();
        wall_closed_cells.exec();
        wall_noslip_cells.exec();
        rel.exec();

        syncVariableToHost<Vecd>(fp, "Velocity");
        Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        Real um = 0;
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            um = std::max(um, vel[i].norm());
        Umax_end = um;
        if (steps % 200 == 0)
            std::cout << "  step~" << steps << " t=" << double(sv_t->getValue()) << " Umax=" << double(um) << "\n";

        // A0/B0: abort if runaway motion (quiet control fail)
        if ((cli.case_id == "A0" || cli.case_id == "B0") && um > Real(1.0e-2))
        {
            std::cout << "[doc205] QUIET_CONTROL_FAIL Umax=" << double(um) << " at t=" << double(sv_t->getValue())
                      << "\n";
            break;
        }
    }

    const bool quiet_ok =
        (cli.case_id == "A0" || cli.case_id == "B0") ? (Umax_end < Real(1.0e-2)) : true;

    // --- end diagnostics: U_rms, section mass-flux proxy, field CSV ---
    syncVariableToHost<Vecd>(fp, "Position");
    syncVariableToHost<Vecd>(fp, "Velocity");
    syncVariableToHost<Real>(fp, "Density");
    syncVariableToHost<Real>(fp, "VolumetricMeasure");
    syncVariableToHost<Real>(fp, kT);
    Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
    Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
    Real *rho = fp.getVariableDataByName<Real>("Density");
    Real *Vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
    Real *T = fp.getVariableDataByName<Real>(kT);

    double Urms = 0, K = 0;
    for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
    {
        const double u2 = double(vel[i].squaredNorm());
        Urms += u2;
        K += 0.5 * double(rho[i] * Vol[i]) * u2;
    }
    Urms = std::sqrt(Urms / double(std::max<size_t>(fp.TotalRealParticles(), 1)));
    const double U_buoy = std::sqrt(double(cli.g * cli.beta * cli.dT_radial * H) + 1e-30);

    // Disk slices at z/H = 0.25, 0.5, 0.75: mdot± = Σ ρ uz Vol / Δz_band  (band proxy)
    const Real band = Real(1.5) * dp;
    struct SliceFlux
    {
        double mplus = 0, mminus = 0;
        int n = 0, n_up = 0, n_dn = 0;
    };
    auto slice = [&](Real z_over_H) -> SliceFlux {
        SliceFlux s;
        const Real z0 = z_over_H * H;
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
        {
            if (std::abs(pos[i][2] - z0) > band)
                continue;
            const double mdot = double(rho[i] * vel[i][2] * Vol[i]) / double(band); // ~∫ρuz dA
            ++s.n;
            if (vel[i][2] > 0)
            {
                s.mplus += mdot;
                ++s.n_up;
            }
            else if (vel[i][2] < 0)
            {
                s.mminus += mdot;
                ++s.n_dn;
            }
        }
        return s;
    };
    const SliceFlux s25 = slice(Real(0.25)), s50 = slice(Real(0.5)), s75 = slice(Real(0.75));
    auto EQ = [](const SliceFlux &s) -> double {
        const double den = std::max(s.mplus, std::abs(s.mminus));
        if (den < 1e-12)
            return -1.0; // NO_RESOLVED_FLOW marker
        return std::abs(s.mplus + s.mminus) / den;
    };

    // core vs wall mean uz (r/R <0.4 vs >0.7) at mid height
    double uz_core = 0, uz_wall = 0, w_core = 0, w_wall = 0;
    for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
    {
        if (std::abs(pos[i][2] - Real(0.5) * H) > band)
            continue;
        const Real rn = std::hypot(pos[i][0], pos[i][1]) / R;
        const double w = double(Vol[i]);
        if (rn < Real(0.4))
        {
            uz_core += w * double(vel[i][2]);
            w_core += w;
        }
        else if (rn > Real(0.7))
        {
            uz_wall += w * double(vel[i][2]);
            w_wall += w;
        }
    }
    if (w_core > 0)
        uz_core /= w_core;
    if (w_wall > 0)
        uz_wall /= w_wall;

    std::error_code ec;
    std::filesystem::create_directories(cli.out_dir, ec);
    const std::string field_csv = cli.out_dir + "/field_end.csv";
    {
        std::ofstream csv(field_csv);
        csv << "x_m,y_m,z_m,ux,uy,uz,T_K,rho,vol\n" << std::setprecision(10);
        for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            csv << double(pos[i][0]) << ',' << double(pos[i][1]) << ',' << double(pos[i][2]) << ','
                << double(vel[i][0]) << ',' << double(vel[i][1]) << ',' << double(vel[i][2]) << ','
                << double(T[i]) << ',' << double(rho[i]) << ',' << double(Vol[i]) << '\n';
    }

    const bool buoyancy_response =
        (cli.beta > 0) && (Urms > Real(10) * Real(3.9e-6)); // vs A0 Umax floor ~4e-6
    const bool both_ways_mid = (s50.n_up > 10 && s50.n_dn > 10 && s50.mplus > 0 && s50.mminus < 0);
    const double eq50 = EQ(s50);

    std::cout << std::scientific << std::setprecision(4)
              << "[doc205] done case=" << cli.case_id << " t=" << double(sv_t->getValue())
              << " Umax=" << double(Umax_end) << " Urms=" << Urms << " U_buoy=" << U_buoy
              << " Urms/U_buoy=" << (Urms / U_buoy)
              << " uz_core=" << uz_core << " uz_wall=" << uz_wall
              << " m+_50=" << s50.mplus << " m-_50=" << s50.mminus << " EQ_50=" << eq50
              << " both_ways_mid=" << (both_ways_mid ? 1 : 0)
              << " quiet_ok=" << (quiet_ok ? 1 : 0) << " WIRING_OK=1"
              << " BUOYANCY_RESPONSE=" << (buoyancy_response ? 1 : 0)
              << " field=" << field_csv << "\n";

    if (!cli.summary_json.empty())
    {
        std::ofstream ofs(cli.summary_json);
        ofs << std::setprecision(16) << "{\n  \"doc\": \"205\", \"case\": \"" << cli.case_id << "\",\n"
            << "  \"WIRING_OK\": true,\n"
            << "  \"QUIET_CONTROL_OK\": " << (quiet_ok ? "true" : "false") << ",\n"
            << "  \"BUOYANCY_RESPONSE\": " << (buoyancy_response ? "true" : "false") << ",\n"
            << "  \"CLOSED_RETURN_CONFIRMED\": false,\n"
            << "  \"B_STATUS\": \""
            << (lid_slip ? "WIRED_VISCous_OMIT_PLUS_ACOUSTIC_FREESLIP" : "N/A") << "\",\n"
            << "  \"lid_slip\": " << (lid_slip ? "true" : "false") << ", \"n_lid_fs\": " << n_lid_fs << ",\n"
            << "  \"n_wall_closed\": " << wp_closed.TotalRealParticles()
            << ", \"n_wall_noslip\": " << wp_noslip.TotalRealParticles() << ",\n"
            << "  \"Umax\": " << double(Umax_end) << ", \"Urms\": " << Urms
            << ", \"U_buoy\": " << U_buoy << ", \"t_end\": " << double(sv_t->getValue()) << ",\n"
            << "  \"uz_core_mid\": " << uz_core << ", \"uz_wall_mid\": " << uz_wall << ",\n"
            << "  \"mdot_plus_z05\": " << s50.mplus << ", \"mdot_minus_z05\": " << s50.mminus
            << ", \"EQ_z05\": " << eq50 << ", \"both_ways_mid\": " << (both_ways_mid ? "true" : "false")
            << ",\n"
            << "  \"mdot_plus_z025\": " << s25.mplus << ", \"mdot_minus_z025\": " << s25.mminus
            << ", \"EQ_z025\": " << EQ(s25) << ",\n"
            << "  \"mdot_plus_z075\": " << s75.mplus << ", \"mdot_minus_z075\": " << s75.mminus
            << ", \"EQ_z075\": " << EQ(s75) << ",\n"
            << "  \"beta\": " << double(cli.beta) << ", \"dp\": " << double(dp) << ", \"n_fluid\": "
            << fp.TotalRealParticles() << ",\n"
            << "  \"use_bg_hydro\": " << (cli.use_bg_hydro ? "true" : "false")
            << ", \"field_csv\": \"" << field_csv << "\"\n}\n";
    }
    return quiet_ok ? 0 : 2;
}
