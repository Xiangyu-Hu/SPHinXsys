/**
 * @file test_3d_ophelie_viscosity_reference.cpp
 * @brief Doc180 P1/P2 + Doc200 B2: SYCL/CK ViscousForceCK vs production InteractionSplit
 *        damping, against analytic channel decay / two-layer Couette / top free-slip decay.
 *
 * Geometry: fluid slab between z=0 and z=H. Default: solid plates above/below (no-slip).
 * Doc200 top_free_slip: viscous contact on bottom plate only (French viscous_omit analog);
 * top has no viscous wall — positions frozen (impermeability stand-in for acoustic lid).
 * Lateral sides have no wall contact (wide aspect). Positions frozen — viscosity only.
 */
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_implicit_viscosity_split.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
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
    std::string problem = "constant_decay"; /**< constant_decay | layered_couette | top_free_slip */
    std::string viscosity_mode = "explicit"; /**< explicit | implicit_split */
    Real dp = 0.01; /**< Doc180: one refinement from 0.0125 if coarse fails ~6% */
    Real H = 0.1;
    Real W = 0.4;   /**< wide aspect; metrics exclude lateral incomplete-support band */
    Real rho0 = 1000.0;
    Real mu = 10.0;
    Real mu1 = 5.0;
    Real mu2 = 50.0; /**< ratio 10 */
    Real h1 = 0.05;  /**< lower layer thickness */
    Real U0 = 0.1;   /**< decay amplitude or Couette top speed */
    Real c0 = 20.0;
    Real t_star = 0.1; /**< ν t / H^2 for decay observation */
    Real dt = 0.0;     /**< 0 = auto from viscous CFL */
    Real end_time = 0.0;
    std::string out_dir = ".";
    std::string summary_json;
};

bool parseCli(int argc, char **argv, Cli &cli)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto take = [&](const char *k, Real &d) {
            const std::string p = std::string(k) + "=";
            if (a.rfind(p, 0) == 0)
            {
                d = static_cast<Real>(std::atof(a.c_str() + p.size()));
                return true;
            }
            return false;
        };
        auto take_s = [&](const char *k, std::string &d) {
            const std::string p = std::string(k) + "=";
            if (a.rfind(p, 0) == 0)
            {
                d = a.substr(p.size());
                return true;
            }
            return false;
        };
        if (take_s("--problem", cli.problem) || take_s("--viscosity-mode", cli.viscosity_mode) ||
            take_s("--out-dir", cli.out_dir) || take_s("--summary-json", cli.summary_json))
            continue;
        if (take("--dp", cli.dp) || take("--H", cli.H) || take("--W", cli.W) || take("--rho0", cli.rho0) ||
            take("--mu", cli.mu) || take("--mu1", cli.mu1) || take("--mu2", cli.mu2) || take("--h1", cli.h1) ||
            take("--U0", cli.U0) || take("--c0", cli.c0) || take("--t-star", cli.t_star) || take("--dt", cli.dt) ||
            take("--end-time", cli.end_time))
            continue;
        std::cerr << "unknown arg: " << a << std::endl;
        return false;
    }
    if (cli.problem != "constant_decay" && cli.problem != "layered_couette" &&
        cli.problem != "top_free_slip")
        return false;
    if (cli.viscosity_mode != "explicit" && cli.viscosity_mode != "implicit_split")
        return false;
    return true;
}

class FluidSlab : public ComplexShape
{
  public:
    FluidSlab(const std::string &name, Real W, Real H) : ComplexShape(name)
    {
        const Vecd half(0.5 * W, 0.5 * W, 0.5 * H);
        const Vecd center(0.5 * W, 0.5 * W, 0.5 * H);
        add<GeometricShapeBox>(Transform(center), half);
    }
};

class WallPlates : public ComplexShape
{
  public:
    WallPlates(const std::string &name, Real W, Real H, Real wall_t) : ComplexShape(name)
    {
        const Vecd half_xy(0.5 * W + wall_t, 0.5 * W + wall_t, 0.5 * wall_t);
        add<GeometricShapeBox>(Transform(Vecd(0.5 * W, 0.5 * W, -0.5 * wall_t)), half_xy);
        add<GeometricShapeBox>(Transform(Vecd(0.5 * W, 0.5 * W, H + 0.5 * wall_t)), half_xy);
    }
};

/** Doc200 B2: bottom plate only — viscous no-slip; top free of viscous contact. */
class WallBottomOnly : public ComplexShape
{
  public:
    WallBottomOnly(const std::string &name, Real W, Real wall_t) : ComplexShape(name)
    {
        const Vecd half_xy(0.5 * W + wall_t, 0.5 * W + wall_t, 0.5 * wall_t);
        add<GeometricShapeBox>(Transform(Vecd(0.5 * W, 0.5 * W, -0.5 * wall_t)), half_xy);
    }
};

Real analyticDecayUx(Real z, Real t, Real U0, Real H, Real nu)
{
    const Real k = Pi / H;
    return U0 * std::sin(k * z) * std::exp(-nu * k * k * t);
}

/** Doc200: u_x = U0 sin(π z /(2H)) exp[-ν (π/(2H))² t]; u=0 at z=0, du/dz=0 at z=H. */
Real analyticTopFreeSlipUx(Real z, Real t, Real U0, Real H, Real nu)
{
    const Real k = Pi / (Real(2) * H);
    return U0 * std::sin(k * z) * std::exp(-nu * k * k * t);
}

Real analyticCouetteUx(Real z, Real H, Real h1, Real mu1, Real mu2, Real U)
{
    const Real h2 = H - h1;
    const Real tau = U / (h1 / mu1 + h2 / mu2);
    if (z <= h1)
        return (tau / mu1) * z;
    return (tau / mu1) * h1 + (tau / mu2) * (z - h1);
}

struct Metrics
{
    Real amp_num = 0;
    Real amp_ref = 0;
    Real amp_rel_err = 1;
    Real L2_rel = 1;
    Real K = 0;
    Real uz_rms = 0;
    Real uz_maxabs = 0;
    Real uz_rms_over_U0 = 1;
    size_t n_interior = 0;
    Real classic_nb_mean = 0;
};

Metrics evaluateDecay(BaseParticles &particles, Real t, Real U0, Real H, Real W, Real nu, Real margin_z,
                      Real margin_xy)
{
    Metrics m;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double num = 0, den = 0, err2 = 0, ref2 = 0, K = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        // Exclude wall layers and open lateral edges (no periodic BC; infinite-plate metric).
        if (z < margin_z || z > H - margin_z || x < margin_xy || x > W - margin_xy || y < margin_xy ||
            y > W - margin_xy)
            continue;
        const Real uref = analyticDecayUx(z, t, U0, H, nu);
        const Real unum = vel[i][0];
        const Real w = vol[i];
        const Real s = std::sin(Pi * z / H);
        num += static_cast<double>(unum * s * w);
        den += static_cast<double>(s * s * w);
        err2 += static_cast<double>((unum - uref) * (unum - uref) * w);
        ref2 += static_cast<double>(uref * uref * w);
        K += 0.5 * static_cast<double>(unum * unum * w);
        ++m.n_interior;
    }
    m.amp_num = den > TinyReal ? static_cast<Real>(num / den) : 0;
    m.amp_ref = U0 * std::exp(-nu * (Pi / H) * (Pi / H) * t);
    m.amp_rel_err = std::abs(m.amp_num - m.amp_ref) /
                    std::max(std::abs(m.amp_ref), static_cast<Real>(TinyReal));
    m.L2_rel = static_cast<Real>(std::sqrt(err2) /
                                 std::max(std::sqrt(ref2), static_cast<double>(TinyReal)));
    m.K = static_cast<Real>(K);
    return m;
}

Metrics evaluateCouette(BaseParticles &particles, Real H, Real W, Real h1, Real mu1, Real mu2, Real U,
                        Real margin_z, Real margin_xy)
{
    Metrics m;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double err2 = 0, ref2 = 0, K = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        if (z < margin_z || z > H - margin_z || x < margin_xy || x > W - margin_xy || y < margin_xy ||
            y > W - margin_xy)
            continue;
        const Real uref = analyticCouetteUx(z, H, h1, mu1, mu2, U);
        const Real unum = vel[i][0];
        const Real w = vol[i];
        err2 += static_cast<double>((unum - uref) * (unum - uref) * w);
        ref2 += static_cast<double>(uref * uref * w);
        K += 0.5 * static_cast<double>(unum * unum * w);
        ++m.n_interior;
    }
    m.L2_rel = static_cast<Real>(std::sqrt(err2) /
                                 std::max(std::sqrt(ref2), static_cast<double>(TinyReal)));
    m.amp_rel_err = m.L2_rel;
    m.K = static_cast<Real>(K);
    return m;
}

/**
 * Doc200 B2 metrics: exclude bottom wall layer + open lateral band.
 * Top incomplete-support band also excluded (no viscous lid); report separately.
 */
Metrics evaluateTopFreeSlip(BaseParticles &particles, Real t, Real U0, Real H, Real W, Real nu,
                            Real margin_z_bot, Real margin_z_top, Real margin_xy)
{
    Metrics m;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real k = Pi / (Real(2) * H);
    double num = 0, den = 0, err2 = 0, ref2 = 0, K = 0, uz2 = 0, vol_sum = 0;
    Real uz_max = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0], y = pos[i][1], z = pos[i][2];
        if (z < margin_z_bot || z > H - margin_z_top || x < margin_xy || x > W - margin_xy ||
            y < margin_xy || y > W - margin_xy)
            continue;
        const Real uref = analyticTopFreeSlipUx(z, t, U0, H, nu);
        const Real unum = vel[i][0];
        const Real w = vol[i];
        const Real s = std::sin(k * z);
        num += static_cast<double>(unum * s * w);
        den += static_cast<double>(s * s * w);
        err2 += static_cast<double>((unum - uref) * (unum - uref) * w);
        ref2 += static_cast<double>(uref * uref * w);
        K += 0.5 * static_cast<double>(unum * unum * w);
        uz2 += static_cast<double>(vel[i][2] * vel[i][2] * w);
        vol_sum += static_cast<double>(w);
        uz_max = std::max(uz_max, std::abs(vel[i][2]));
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
    m.uz_maxabs = uz_max;
    m.uz_rms_over_U0 = m.uz_rms / std::max(std::abs(U0), static_cast<Real>(TinyReal));
    return m;
}
} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
        return 2;

    const Real wall_t = 3.0 * cli.dp;
    const Real nu = cli.mu / cli.rho0;
    const bool top_free = (cli.problem == "top_free_slip");
    if (!(cli.end_time > TinyReal))
    {
        if (cli.problem == "constant_decay" || top_free)
            cli.end_time = cli.t_star * cli.H * cli.H / nu;
        else
            cli.end_time = Real(2.0) * cli.H * cli.H / (cli.mu1 / cli.rho0); // settle ~ viscous time
    }
    if (!(cli.dt > TinyReal))
    {
        const Real mu_ref = cli.problem == "layered_couette" ? std::max(cli.mu1, cli.mu2) : cli.mu;
        // Explicit viscous CFL safety; implicit_split may still use this as baseline step.
        cli.dt = Real(0.1) * (cli.dp * cli.dp) * cli.rho0 / mu_ref;
    }

    std::filesystem::create_directories(cli.out_dir);
    const BoundingBoxd bounds(Vecd(-wall_t, -wall_t, -wall_t),
                              Vecd(cli.W + wall_t, cli.W + wall_t, cli.H + wall_t));
    SPHSystem sph_system(bounds, cli.dp);
    FluidBody fluid(sph_system, makeShared<FluidSlab>("Fluid", cli.W, cli.H));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(cli.rho0, cli.c0);
    fluid.addMaterialProperty<ParticleViscosity>(cli.mu);
    fluid.generateParticles<BaseParticles, Lattice>();
    ParticleViscosity &pvisc = fluid.getMaterialProperty<ParticleViscosity>();

    // Doc200 B2: bottom-only viscous wall (French viscous_omit analog). Else both plates.
    SharedPtr<Shape> wall_shape =
        top_free ? SharedPtr<Shape>(makeShared<WallBottomOnly>("WallBottom", cli.W, wall_t))
                 : SharedPtr<Shape>(makeShared<WallPlates>("Wall", cli.W, cli.H, wall_t));
    SolidBody wall(sph_system, wall_shape);
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> fluid_inner(fluid);
    Contact<> fluid_wall(fluid, {&wall});
    UniquePtr<InnerRelation> fluid_inner_classic;
    UniquePtr<ContactRelation> fluid_wall_classic;
    if (cli.viscosity_mode == "implicit_split")
    {
        fluid_inner_classic = makeUnique<InnerRelation>(fluid);
        fluid_wall_classic = makeUnique<ContactRelation>(fluid, RealBodyVector{&wall});
    }

    BaseParticles &fp = fluid.getBaseParticles();
    BaseParticles &wp = wall.getBaseParticles();
    // ViscousForceCK getVariable("Velocity"); AcousticStep normally registers it — do the same here.
    fp.registerStateVariable<Vecd>("Velocity");
    fp.registerStateVariable<Vecd>("ForcePrior");
    fp.registerStateVariable<Vecd>("Force");
    wp.registerStateVariable<Vecd>("Velocity");
    fp.registerStateVariable<Real>(kOphelieDynamicViscosityField, cli.mu);
    pvisc.bindParticleViscosityField(fp, kOphelieDynamicViscosityField);
    syncVariableToDevice<Real>(fp, kOphelieDynamicViscosityField);

    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal(wall);
    wall_normal.exec();

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fluid_cll(fluid);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cll(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> update_rel(fluid_inner, fluid_wall);

    using ViscousCK =
        fluid_dynamics::ViscousForceCK<Inner<WithUpdate, ParticleViscosity, NoKernelCorrectionCK>,
                                       Contact<Wall, ParticleViscosity, NoKernelCorrectionCK>>;
    InteractionDynamicsCK<MainExecutionPolicy, ViscousCK> viscous_ck(fluid_inner, fluid_wall);

    UniquePtr<InteractionSplit<DampingPairwiseInner<Vecd, ParticleViscosityHarmonicDampingRate>>> damp_inner;
    UniquePtr<InteractionSplit<DampingPairwiseContactWall<Vecd, ParticleViscosityOneSideDampingRate>>> damp_wall;
    if (cli.viscosity_mode == "implicit_split")
    {
        damp_inner =
            makeUnique<InteractionSplit<DampingPairwiseInner<Vecd, ParticleViscosityHarmonicDampingRate>>>(
                *fluid_inner_classic, "Velocity");
        damp_wall =
            makeUnique<InteractionSplit<DampingPairwiseContactWall<Vecd, ParticleViscosityOneSideDampingRate>>>(
                *fluid_wall_classic, "Velocity");
    }

    auto update_classic = [&]() {
        if (!fluid_inner_classic)
            return;
        syncVariableToHost<Vecd>(fp, "Position");
        syncVariableToHost<Vecd>(wp, "Position");
        fluid.updateCellLinkedList();
        wall.updateCellLinkedList();
        fluid_inner_classic->updateConfiguration();
        fluid_wall_classic->updateConfiguration();
    };

    // Init fields
    syncVariableToHost<Vecd>(fp, "Position");
    syncVariableToHost<Vecd>(fp, "Velocity");
    syncVariableToHost<Real>(fp, kOphelieDynamicViscosityField);
    Vecd *fpos = fp.getVariableDataByName<Vecd>("Position");
    Vecd *fvel = fp.getVariableDataByName<Vecd>("Velocity");
    Real *fmu = fp.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    const size_t nf = fp.TotalRealParticles();
    for (size_t i = 0; i < nf; ++i)
    {
        const Real z = fpos[i][2];
        if (cli.problem == "constant_decay")
        {
            fmu[i] = cli.mu;
            fvel[i] = Vecd(analyticDecayUx(z, Real(0), cli.U0, cli.H, nu), 0, 0);
        }
        else if (top_free)
        {
            fmu[i] = cli.mu;
            fvel[i] = Vecd(analyticTopFreeSlipUx(z, Real(0), cli.U0, cli.H, nu), 0, 0);
        }
        else
        {
            fmu[i] = (z < cli.h1) ? cli.mu1 : cli.mu2;
            fvel[i] = Vecd(analyticCouetteUx(z, cli.H, cli.h1, cli.mu1, cli.mu2, cli.U0), 0, 0);
            // Small perturbation for regression (Doc180): not only hold IC.
            fvel[i][0] += Real(0.02) * cli.U0 * std::sin(Real(2) * Pi * z / cli.H);
        }
    }
    syncVariableToDevice<Real>(fp, kOphelieDynamicViscosityField);
    syncVariableToDevice<Vecd>(fp, "Velocity");

    // Wall velocities
    syncVariableToHost<Vecd>(wp, "Position");
    syncVariableToHost<Vecd>(wp, "Velocity");
    Vecd *wpos = wp.getVariableDataByName<Vecd>("Position");
    Vecd *wvel = wp.getVariableDataByName<Vecd>("Velocity");
    const size_t nw = wp.TotalRealParticles();
    for (size_t i = 0; i < nw; ++i)
    {
        wvel[i] = Vecd::Zero();
        if (cli.problem == "layered_couette" && wpos[i][2] > cli.H)
            wvel[i] = Vecd(cli.U0, 0, 0);
    }
    syncVariableToDevice<Vecd>(wp, "Velocity");

    fluid_cll.exec();
    wall_cll.exec();
    update_rel.exec();
    update_classic();

    Real classic_nb_mean = 0;
    if (fluid_inner_classic)
    {
        size_t s = 0;
        for (size_t i = 0; i < nf; ++i)
            s += fluid_inner_classic->inner_configuration_[i].current_size_;
        classic_nb_mean = static_cast<Real>(s) / static_cast<Real>(nf);
    }

    std::ofstream hist(cli.out_dir + "/history.csv");
    hist << "step,time_s,amp_num,amp_ref,amp_rel_err,L2_rel,K,uz_rms_over_U0\n";

    Real t = 0;
    size_t step = 0;
    const Real margin_z = Real(2.5) * cli.dp;
    const Real margin_z_top = top_free ? Real(2.5) * cli.dp : margin_z; /**< incomplete support at free top */
    const Real margin_xy = Real(3.5) * cli.dp; /**< open lateral edges */
    Metrics last{};
    auto write_profile = [&](const std::string &path) {
        syncVariableToHost<Vecd>(fp, "Position");
        syncVariableToHost<Vecd>(fp, "Velocity");
        const Vecd *pos = fp.getVariableDataByName<Vecd>("Position");
        const Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
        std::ofstream pf(path);
        pf << "x,y,z,ux,uy,uz,ux_ref\n";
        for (size_t i = 0; i < nf; ++i)
        {
            const Real z = pos[i][2];
            Real uref = 0;
            if (cli.problem == "constant_decay")
                uref = analyticDecayUx(z, t, cli.U0, cli.H, nu);
            else if (top_free)
                uref = analyticTopFreeSlipUx(z, t, cli.U0, cli.H, nu);
            else
                uref = analyticCouetteUx(z, cli.H, cli.h1, cli.mu1, cli.mu2, cli.U0);
            pf << pos[i][0] << "," << pos[i][1] << "," << z << "," << vel[i][0] << "," << vel[i][1] << ","
               << vel[i][2] << "," << uref << "\n";
        }
    };
    if (top_free)
        std::cout << "[doc200][b2] top_free_slip viscous_wall=bottom_only acoustic_lid=frozen_positions "
                     "analytic=U0*sin(pi*z/(2H))*exp(-nu*(pi/(2H))^2*t)"
                  << std::endl;
    while (t + TinyReal < cli.end_time)
    {
        const Real dt_step = std::min(cli.dt, cli.end_time - t);
        if (cli.viscosity_mode == "explicit")
        {
            hostClearViscousForcePriorContribution(fp);
            viscous_ck.exec();
            syncVariableToHost<Vecd>(fp, "ViscousForce");
            syncVariableToHost<Vecd>(fp, "Velocity");
            syncVariableToHost<Real>(fp, "Mass");
            Vecd *vf = fp.getVariableDataByName<Vecd>("ViscousForce");
            Vecd *vel = fp.getVariableDataByName<Vecd>("Velocity");
            Real *mass = fp.getVariableDataByName<Real>("Mass");
            for (size_t i = 0; i < nf; ++i)
                vel[i] += vf[i] / std::max(mass[i], TinyReal) * dt_step;
            syncVariableToDevice<Vecd>(fp, "Velocity");
        }
        else
        {
            syncVariableToHost<Vecd>(fp, "Velocity");
            syncVariableToHost<Real>(fp, kOphelieDynamicViscosityField);
            syncVariableToHost<Real>(fp, "Mass");
            syncVariableToHost<Real>(fp, "VolumetricMeasure");
            syncVariableToHost<Vecd>(wp, "Velocity");
            syncVariableToHost<Real>(wp, "VolumetricMeasure");
            damp_inner->exec(dt_step);
            damp_wall->exec(dt_step);
            syncVariableToDevice<Vecd>(fp, "Velocity");
        }
        t += dt_step;
        ++step;
        if (step % 20 == 0 || t + TinyReal >= cli.end_time)
        {
            if (cli.problem == "constant_decay")
                last = evaluateDecay(fp, t, cli.U0, cli.H, cli.W, nu, margin_z, margin_xy);
            else if (top_free)
                last = evaluateTopFreeSlip(fp, t, cli.U0, cli.H, cli.W, nu, margin_z, margin_z_top, margin_xy);
            else
                last = evaluateCouette(fp, cli.H, cli.W, cli.h1, cli.mu1, cli.mu2, cli.U0, margin_z, margin_xy);
            last.classic_nb_mean = classic_nb_mean;
            hist << step << "," << t << "," << last.amp_num << "," << last.amp_ref << "," << last.amp_rel_err
                 << "," << last.L2_rel << "," << last.K << "," << last.uz_rms_over_U0 << "\n";
            hist.flush();
            std::cout << "[doc180][visc-ref] step=" << step << " t=" << t << " L2_rel=" << last.L2_rel
                      << " amp_rel=" << last.amp_rel_err << " uz/U0=" << last.uz_rms_over_U0
                      << " n_int=" << last.n_interior << std::endl;
        }
    }
    write_profile(cli.out_dir + "/profile_end.csv");

    const bool nb_ok = (cli.viscosity_mode == "explicit") || (classic_nb_mean > Real(1));
    // Doc200: profile L2≤5%, amp≤5%, uz_rms/U0≤1%.
    const bool top_ok =
        (!top_free) || (last.amp_rel_err <= Real(0.05) && last.uz_rms_over_U0 <= Real(0.01));
    const bool passed = nb_ok && last.L2_rel <= Real(0.05) &&
                        (cli.problem == "layered_couette" || last.amp_rel_err <= Real(0.05)) && top_ok;

    const std::string sj =
        cli.summary_json.empty() ? (cli.out_dir + "/summary.json") : cli.summary_json;
    std::ofstream js(sj);
    const Real t_star_out =
        (cli.problem == "constant_decay" || top_free) ? (nu * t) / (cli.H * cli.H) : 0.0;
    js << "{\n"
       << "  \"problem\": \"" << cli.problem << "\",\n"
       << "  \"viscosity_mode\": \"" << cli.viscosity_mode << "\",\n"
       << "  \"backend_note\": \"explicit=ViscousForceCK(SYCL); implicit=InteractionSplit+DampingPairwise(CPU) "
          "with Doc177 CK-match factors 0.5/2\",\n"
       << "  \"top_bc\": \"" << (top_free ? "flat_free_slip_viscous_omit_bottom_only" : "both_plates_noslip")
       << "\",\n"
       << "  \"dp\": " << cli.dp << ",\n"
       << "  \"H\": " << cli.H << ",\n"
       << "  \"end_time_s\": " << t << ",\n"
       << "  \"dt_s\": " << cli.dt << ",\n"
       << "  \"t_star\": " << t_star_out << ",\n"
       << "  \"n_steps\": " << step << ",\n"
       << "  \"classic_nb_mean\": " << classic_nb_mean << ",\n"
       << "  \"n_interior\": " << last.n_interior << ",\n"
       << "  \"amp_num\": " << last.amp_num << ",\n"
       << "  \"amp_ref\": " << last.amp_ref << ",\n"
       << "  \"amp_rel_err\": " << last.amp_rel_err << ",\n"
       << "  \"L2_rel\": " << last.L2_rel << ",\n"
       << "  \"uz_rms\": " << last.uz_rms << ",\n"
       << "  \"uz_maxabs\": " << last.uz_maxabs << ",\n"
       << "  \"uz_rms_over_U0\": " << last.uz_rms_over_U0 << ",\n"
       << "  \"K\": " << last.K << ",\n"
       << "  \"ck_match_inner\": " << kOphelieCkMatchInnerDampingFactor << ",\n"
       << "  \"ck_match_wall\": " << kOphelieCkMatchWallDampingFactor << ",\n"
       << "  \"viscosity_constant_passed\": "
       << (cli.problem == "constant_decay" && passed ? "true" : "false") << ",\n"
       << "  \"viscosity_variable_wall_passed\": "
       << (cli.problem == "layered_couette" && passed ? "true" : "false") << ",\n"
       << "  \"top_free_slip_passed\": " << (top_free && passed ? "true" : "false") << ",\n"
       << "  \"passed\": " << (passed ? "true" : "false") << "\n"
       << "}\n";
    std::cout << "[doc180][visc-ref] passed=" << passed << " L2_rel=" << last.L2_rel
              << " uz_rms/U0=" << last.uz_rms_over_U0 << " classic_nb_mean=" << classic_nb_mean
              << " summary=" << sj << std::endl;
    return passed ? 0 : 1;
}
