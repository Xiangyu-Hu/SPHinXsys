/**
 * @file spdb_programmed_temperature_sycl.cpp
 * @brief SYCL/GPU dual-slot programmed-temperature demo.
 * @details Uses ContactPenaltyForceCK (e_ij / ContactForce formula), NOT
 *          RepulsionForceCK (flush-face unsafe). Materials/thermal/gate are
 *          demo_assumption. Backend must be SPHINXSYS_USE_SYCL=ON.
 */
#include "sphinxsys.h"

#include "spdb_ck_axial.h"
#include "spdb_contact_penalty_ck.h"
#include "spdb_contact_scaled.h"
#include "spdb_dual_slot_metrics.h"
#include "spdb_ren2014_solid.h"
#include "spdb_temperature_program.h"
#include "spdb_thermal_expansion.h"
#include "spdb_velocity_decay_ck.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

using namespace SPH;

namespace
{
constexpr Real demo_E = 120.0e9;
constexpr Real demo_nu = 0.30;
constexpr Real demo_rho = 7700.0;
constexpr Real demo_alpha = 1.2e-5;
constexpr Real demo_T0_K = 293.15;
constexpr Real demo_T0_C = 20.0;
constexpr Real cylinder_radius = 0.020;
constexpr Real specimen_half_stack = 0.020;
constexpr Real nominal_pressure = 5.0e6;
constexpr Real mesh_volume_ref = 2.38170881694208e-5;
// Device light calm + sparse host InteractionSplit pairwise damp (test_3d_sliding / CPU SPDB).
constexpr Real velocity_decay_default = 5.0e4; // 1/s
constexpr Real physical_viscosity_default = 5.0e6;
constexpr Real damp_random_default = 0.5; // sliding uses 0.5; CPU SPDB used 0.2
constexpr int host_damp_stride_default = 10; // every N acoustic steps

Real nominal_force() { return nominal_pressure * Pi * cylinder_radius * cylinder_radius; }

BoundingBoxd domain_bounds(Real dp)
{
    return BoundingBoxd(Vec3d(-cylinder_radius - 2.0 * dp, -cylinder_radius - 2.0 * dp,
                              -specimen_half_stack - 2.0 * dp),
                        Vec3d(cylinder_radius + 2.0 * dp, cylinder_radius + 2.0 * dp,
                              specimen_half_stack + 2.0 * dp));
}

class LowerFromSTL : public ComplexShape
{
  public:
    explicit LowerFromSTL(const std::string &file) : ComplexShape("ProgramTempLower")
    {
        add<TriangleMeshShapeSTL>(file, Vec3d::Zero(), 1.0);
    }
};

class UpperFromSTL : public ComplexShape
{
  public:
    explicit UpperFromSTL(const std::string &file) : ComplexShape("ProgramTempUpper")
    {
        add<TriangleMeshShapeSTL>(file, Vec3d::Zero(), 1.0);
    }
};

void sync_ren2014_host(BaseParticles &p)
{
    // Public API: device -> host (synchronizeWithDevice is private).
    const ParallelDevicePolicy pol{};
    p.getVariableByName<Matd>("DeformationGradient")->prepareForOutput(pol);
    p.getVariableByName<Matd>("InversePlasticRightCauchyStrain")->prepareForOutput(pol);
    p.getVariableByName<Real>("EquivalentInelasticStrain")->prepareForOutput(pol);
    p.getVariableByName<Real>("ArrheniusRangeFlag")->prepareForOutput(pol);
    p.getVariableByName<Vecd>("Position")->prepareForOutput(pol);
    p.getVariableByName<Real>("Density")->prepareForOutput(pol);
    p.getVariableByName<Vecd>("Velocity")->prepareForOutput(pol);
}

void push_ren2014_device(BaseParticles &p)
{
    // Public API: host -> device (synchronizeToDevice is private).
    const ParallelDevicePolicy pol{};
    p.getVariableByName<Matd>("DeformationGradient")->finalizeLoadIn(pol);
    p.getVariableByName<Matd>("InversePlasticRightCauchyStrain")->finalizeLoadIn(pol);
    p.getVariableByName<Real>("EquivalentInelasticStrain")->finalizeLoadIn(pol);
    p.getVariableByName<Real>("ArrheniusRangeFlag")->finalizeLoadIn(pol);
    p.getVariableByName<Vecd>("Position")->finalizeLoadIn(pol);
    p.getVariableByName<Real>("Density")->finalizeLoadIn(pol);
    p.getVariableByName<Real>("Temperature")->finalizeLoadIn(pol);
}

Real mean_inelastic_host(BaseParticles &p)
{
    Real *eq = p.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real sum = 0.0;
    const size_t n = p.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        sum += eq[i];
    return n > 0 ? sum / static_cast<Real>(n) : 0.0;
}

Real particle_volume_sum(BaseParticles &p)
{
    Real *vol = p.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
        sum += vol[i];
    return sum;
}

enum class CaseMode
{
    GeometryOnly,
    F1,
    F2,
    F3
};

struct Options
{
    CaseMode mode = CaseMode::F3;
    Real dp = 0.001;
    Real end_time_s = 60.0;
    Real macro_dt = 10.0;
    Real max_macro_dt = 30.0;
    Real stiffness_scale = 1.0;
    Real relax_cap = 1.5e-3;
    Real velocity_decay = velocity_decay_default;
    Real physical_viscosity = physical_viscosity_default;
    Real damp_random = damp_random_default;
    int host_damp_stride = host_damp_stride_default;
    Real vmax_calm = 2.0e-3;
    Real vtp_every = 100.0; // write VTP every N physical seconds (and always at end)
    bool apply_load = true;
    bool enable_thermal = true;
    bool enable_inelastic = true;
    bool use_program_temperature = true;
};

Options parse_args(int argc, char *argv[])
{
    Options o;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto need = [&](Real &dst) {
            if (i + 1 < argc)
                dst = std::atof(argv[++i]);
        };
        if (a == "--geometry_only")
            o.mode = CaseMode::GeometryOnly;
        else if (a == "--F1")
        {
            o.mode = CaseMode::F1;
            o.apply_load = true;
            o.enable_thermal = false;
            o.enable_inelastic = false;
            o.use_program_temperature = false;
            o.end_time_s = 0.0;
        }
        else if (a == "--F2")
        {
            o.mode = CaseMode::F2;
            o.enable_inelastic = false;
        }
        else if (a == "--F3")
            o.mode = CaseMode::F3;
        else if (a == "--dp")
            need(o.dp);
        else if (a == "--end_time")
            need(o.end_time_s);
        else if (a == "--macro_dt")
            need(o.macro_dt);
        else if (a == "--max_macro_dt")
            need(o.max_macro_dt);
        else if (a == "--contact_stiffness_scale")
            need(o.stiffness_scale);
        else if (a == "--relax_cap")
            need(o.relax_cap);
        else if (a == "--physical_viscosity")
            need(o.physical_viscosity);
        else if (a == "--damp_random")
            need(o.damp_random);
        else if (a == "--host_damp_stride")
        {
            if (i + 1 < argc)
                o.host_damp_stride = std::atoi(argv[++i]);
        }
        else if (a == "--velocity_decay")
            need(o.velocity_decay);
        else if (a == "--vmax_calm")
            need(o.vmax_calm);
        else if (a == "--vtp_every")
            need(o.vtp_every);
    }
    return o;
}

/** Host-side J2 stress using Ren2014 Cauchy (plastic history aware). */
void fill_ren2014_equivalent_stress(BaseParticles &p, Real G0, Real K0)
{
    Matd *F = p.getVariableDataByName<Matd>("DeformationGradient");
    Matd *Cp = p.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *seq = p.getVariableDataByName<Real>("Ren2014EquivalentStress");
    const size_t n = p.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        Real s = 0.0;
        ren2014_cauchy_stress(F[i], Cp[i], G0, K0, &s, nullptr);
        seq[i] = s;
    }
}

const char *mode_name(CaseMode m)
{
    switch (m)
    {
    case CaseMode::GeometryOnly:
        return "geometry_only";
    case CaseMode::F1:
        return "F1_elastic_20C";
    case CaseMode::F2:
        return "F2_thermal_elastic";
    case CaseMode::F3:
        return "F3_thermal_arrhenius";
    }
    return "unknown";
}

std::vector<Real> build_landing_times(const spdb::TemperatureProgram &prog, Real end_s, Real macro_dt,
                                      Real max_dt)
{
    std::set<Real> marks{0.0, end_s};
    for (const auto &k : prog.knots())
        if (k.first <= end_s + 1.0e-12)
            marks.insert(k.first);
    for (Real t : prog.crossing_times_C(850.0))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    for (Real t : prog.crossing_times_C(900.0))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    // Extra landings approaching the prior abort band (~775 C / ~990 s).
    for (Real t : prog.crossing_times_C(700.0))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    for (Real t : prog.crossing_times_C(750.0))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    for (Real t : prog.crossing_times_C(800.0))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    std::vector<Real> sorted(marks.begin(), marks.end());
    std::vector<Real> out{0.0};
    for (size_t i = 1; i < sorted.size(); ++i)
    {
        Real a = sorted[i - 1];
        const Real b = sorted[i];
        while (b - a > 1.0e-12)
        {
            const Real T_a = prog.temperature_C(a);
            Real local_dt = macro_dt;
            if (T_a >= 700.0)
                local_dt = std::min(local_dt, Real(5.0));
            if (T_a >= 750.0)
                local_dt = std::min(local_dt, Real(2.5));
            if (T_a >= 800.0)
                local_dt = std::min(local_dt, Real(1.0));
            local_dt = std::min(local_dt, max_dt);
            if (b - a <= local_dt + 1.0e-12)
                break;
            a += local_dt;
            out.push_back(a);
        }
        out.push_back(b);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(),
                          [](Real x, Real y) { return std::abs(x - y) < 1.0e-9; }),
              out.end());
    return out;
}

void write_status(const char *phase, Real physical_s, Real end_s, Real T_C, Real gate, int finished)
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    char iso[64];
    std::strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&tt));
    std::ofstream out("hold_status.txt");
    out << std::setprecision(10);
    out << "updated_at=" << iso << "\nphase=" << phase << "\nphysical_s=" << physical_s
        << "\nend_s=" << end_s << "\nprogress_frac=" << (end_s > 0 ? physical_s / end_s : 0.0)
        << "\ntemperature_C=" << T_C << "\ngate=" << gate << "\nbackend=SYCL_GPU\nfinished=" << finished
        << "\n";
}

} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    std::cerr << "This binary requires SPHINXSYS_USE_SYCL=ON (GPU/SYCL build).\n";
    return 2;
#endif
    Options opt = parse_args(ac, av);
    try
    {
        spdb::TemperatureProgram prog;
        prog.load_csv("input/temperature_program.csv");

        const Real dp = opt.dp;
        const Real force_total = opt.apply_load ? nominal_force() : 0.0;
        const Real slab = 0.55 * dp;
        Real thermal_stretch = 1.0;
        Real physical_clock = 0.0;

        std::cout << "SPDB programmed_temperature SYCL\n"
                  << "  mode=" << mode_name(opt.mode) << "\n"
                  << "  backend=SPHINXSYS_USE_SYCL contact=ContactPenaltyForceCK(e_ij)\n"
                  << "  material=demo_assumption E=" << demo_E << " alpha=" << demo_alpha << "\n"
                  << "  NOT RepulsionForceCK; NOT measured 20Cr13 full-T prediction.\n"
                  << "  force_N=" << force_total << " dp=" << dp << " end_s=" << opt.end_time_s
                  << " macro_dt=" << opt.macro_dt << "\n"
                  << "  relax_cap=" << opt.relax_cap << " velocity_decay=" << opt.velocity_decay
                  << " vmax_calm=" << opt.vmax_calm << "\n"
                  << "  host_split_damp eta=" << opt.physical_viscosity
                  << " random=" << opt.damp_random << " stride=" << opt.host_damp_stride << "\n"
                  << "  damping=VelocityDecayCK + host InteractionSplit/DampingPairwiseInner\n"
                  << "  vtp fields=Temperature,EquivalentInelasticStrain,VonMisesStress,"
                     "VonMisesStrain,Displacement,Ren2014EquivalentStress"
                  << " every_s=" << opt.vtp_every << "\n";

        SPHSystem sph_system(domain_bounds(dp), dp);
        // Do not call handleCommandlineOptions(ac, av): Boost only accepts library
        // flags (--relax/--reload/...). Our --F1/--geometry_only/--dp are handled by parse_args.

        auto &lower_shape = sph_system.addShape<LowerFromSTL>("lower_m.stl");
        auto &upper_shape = sph_system.addShape<UpperFromSTL>("upper_m.stl");

        auto &lower_body = sph_system.addBody<SolidBody>(lower_shape);
        lower_body.defineAdaptationRatios(1.15, 1.0);
        ContactScaledRen2014 &lower_mat = lower_body.defineMatterMaterial<ContactScaledRen2014>(
            demo_rho, demo_E, demo_nu, demo_T0_K, false, opt.stiffness_scale);
        lower_body.defineBodyLevelSetShape().correctLevelSetSign();
        lower_body.generateParticles<BaseParticles, Lattice>();

        auto &upper_body = sph_system.addBody<SolidBody>(upper_shape);
        upper_body.defineAdaptationRatios(1.15, 1.0);
        ContactScaledRen2014 &upper_mat = upper_body.defineMatterMaterial<ContactScaledRen2014>(
            demo_rho, demo_E, demo_nu, demo_T0_K, false, opt.stiffness_scale);
        upper_body.defineBodyLevelSetShape().correctLevelSetSign();
        upper_body.generateParticles<BaseParticles, Lattice>();

        BaseParticles &lower_particles = lower_body.getBaseParticles();
        BaseParticles &upper_particles = upper_body.getBaseParticles();
        lower_particles.registerStateVariable<Real>("Temperature",
                                                    [&](size_t) -> Real { return demo_T0_K; });
        upper_particles.registerStateVariable<Real>("Temperature",
                                                    [&](size_t) -> Real { return demo_T0_K; });
        // Extra viz fields for ParaView (stress / strain / displacement).
        lower_particles.registerStateVariable<Real>("Ren2014EquivalentStress",
                                                    [&](size_t) -> Real { return 0.0; });
        upper_particles.registerStateVariable<Real>("Ren2014EquivalentStress",
                                                    [&](size_t) -> Real { return 0.0; });

        const Real vol_l = particle_volume_sum(lower_particles);
        const Real vol_u = particle_volume_sum(upper_particles);
        Vecd *pos_l0 = lower_particles.ParticlePositions();
        Vecd *pos_u0 = upper_particles.ParticlePositions();
        const spdb::DualSlotCavityMetrics cav_initial =
            spdb::measure_dual_slot_cavity(pos_l0, lower_particles.TotalRealParticles(), pos_u0,
                                           upper_particles.TotalRealParticles(), dp);
        std::cout << "  lower_n=" << lower_particles.TotalRealParticles()
                  << " upper_n=" << upper_particles.TotalRealParticles()
                  << " vol_rel_err_l=" << std::abs(vol_l - mesh_volume_ref) / mesh_volume_ref
                  << " cavity0_m=" << cav_initial.center_net_height
                  << " floor_n=" << cav_initial.floor_n << " roof_n=" << cav_initial.roof_n << "\n";

        if (opt.mode == CaseMode::GeometryOnly)
        {
            std::ofstream out("spdb_program_temp_sycl_geometry_metrics.txt");
            out << std::setprecision(10);
            out << "case=test_3d_spdb_programmed_temperature_sycl\nbackend=SYCL\n";
            out << "lower_particles=" << lower_particles.TotalRealParticles() << "\n";
            out << "upper_particles=" << upper_particles.TotalRealParticles() << "\n";
            out << "particle_volume_lower_m3=" << vol_l << "\n";
            out << "particle_volume_upper_m3=" << vol_u << "\n";
            out << "mesh_volume_ref_m3=" << mesh_volume_ref << "\n";
            out << "initial_cavity_net_height_m=" << cav_initial.center_net_height << "\n";
            out << "floor_n=" << cav_initial.floor_n << "\nroof_n=" << cav_initial.roof_n << "\n";
            write_status("geometry_only", 0, 0, demo_T0_C, 0, 1);
            return 0;
        }

        GeometricShapeBox bottom_shape(Transform(Vec3d(0.0, 0.0, -specimen_half_stack + 0.5 * dp)),
                                       Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab),
                                       "BottomAxialSupport");
        GeometricShapeBox top_shape(Transform(Vec3d(0.0, 0.0, specimen_half_stack - 0.5 * dp)),
                                    Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "TopAxialLoad");
        auto &bottom_part = lower_body.addBodyPart<BodyRegionByParticle>(bottom_shape);
        auto &top_part = upper_body.addBodyPart<BodyRegionByParticle>(top_shape);
        if (bottom_part.body_part_particles_.empty() ||
            (opt.apply_load && top_part.body_part_particles_.empty()))
        {
            std::cerr << "empty support/load layer\n";
            return 1;
        }
        const Real force_per =
            opt.apply_load ? -force_total / static_cast<Real>(top_part.body_part_particles_.size()) : 0.0;

        auto &lower_inner = sph_system.addInnerRelation(lower_body, ConfigType::Lagrangian);
        auto &upper_inner = sph_system.addInnerRelation(upper_body, ConfigType::Lagrangian);
        auto &lower_upper = sph_system.addContactRelation(lower_body, upper_body);
        auto &upper_lower = sph_system.addContactRelation(upper_body, lower_body);

        SPHSolver sph_solver(sph_system);
        auto &main_methods = sph_solver.getMainMethodContainer();

        ParticleDynamicsGroup update_cfg;
        update_cfg.add(&main_methods.addCellLinkedListDynamics(lower_body));
        update_cfg.add(&main_methods.addCellLinkedListDynamics(upper_body));
        update_cfg.add(&main_methods.addRelationDynamics(lower_inner));
        update_cfg.add(&main_methods.addRelationDynamics(upper_inner));
        update_cfg.add(&main_methods.addRelationDynamics(lower_upper));
        update_cfg.add(&main_methods.addRelationDynamics(upper_lower));

        auto &corr_l = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(lower_inner);
        auto &corr_u = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(upper_inner);
        auto &step1_l = main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration1stHalf, ContactScaledRen2014, NoKernelCorrectionCK>(
            lower_inner);
        auto &step2_l =
            main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(
                lower_inner);
        auto &step1_u = main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration1stHalf, ContactScaledRen2014, NoKernelCorrectionCK>(
            upper_inner);
        auto &step2_u =
            main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(
                upper_inner);

        auto &factor_l = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(lower_upper);
        auto &factor_u = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(upper_lower);
        auto &force_l =
            main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::ContactPenaltyForceCK>(
                lower_upper);
        auto &force_u =
            main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::ContactPenaltyForceCK>(
                upper_lower);

        auto &dt_l = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(lower_body, 0.2);
        auto &dt_u = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(upper_body, 0.2);
        auto &support = main_methods.addStateDynamics<AxialZSupportCK>(bottom_part);
        auto &decay_l =
            main_methods.addStateDynamics<VelocityDecayCK>(lower_body, opt.velocity_decay);
        auto &decay_u =
            main_methods.addStateDynamics<VelocityDecayCK>(upper_body, opt.velocity_decay);

        // Exact sliding / CPU SPDB damping on host (cell-split Gauss–Seidel).
        // Sync pos/vel -> host CLL + InnerRelation -> damp -> push vel. Sparse stride
        // so PCIe cost stays acceptable; device VelocityDecayCK remains every step.
        InnerRelation lower_inner_host(lower_body);
        InnerRelation upper_inner_host(upper_body);
        DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
            damp_l_host(opt.damp_random, lower_inner_host, "Velocity", opt.physical_viscosity);
        DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
            damp_u_host(opt.damp_random, upper_inner_host, "Velocity", opt.physical_viscosity);

        SingleVariable<Real> force_scale("ForceScale", opt.apply_load ? 1.0 : 0.0);
        auto &apply_force =
            main_methods.addStateDynamics<NominalAxialForceCK>(top_part, force_per, force_scale);

        auto &recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
        // Scalar/history fields already on particles.
        recorder.addToWrite<Real>(lower_body, "Temperature");
        recorder.addToWrite<Real>(upper_body, "Temperature");
        recorder.addToWrite<Real>(lower_body, "EquivalentInelasticStrain");
        recorder.addToWrite<Real>(upper_body, "EquivalentInelasticStrain");
        recorder.addToWrite<Real>(lower_body, "Ren2014EquivalentStress");
        recorder.addToWrite<Real>(upper_body, "Ren2014EquivalentStress");
        // Derived after device->host sync inside VtpCK::writeToFile (host SimpleDynamics).
        recorder.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(lower_body);
        recorder.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(upper_body);
        recorder.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(lower_body);
        recorder.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(upper_body);
        recorder.addDerivedVariableRecording<SimpleDynamics<Displacement>>(lower_body);
        recorder.addDerivedVariableRecording<SimpleDynamics<Displacement>>(upper_body);

        auto write_vtp = [&](size_t ite = std::numeric_limits<size_t>::max()) {
            sync_ren2014_host(lower_particles);
            sync_ren2014_host(upper_particles);
            fill_ren2014_equivalent_stress(lower_particles, lower_mat.ShearModulus(),
                                           lower_mat.BulkModulus());
            fill_ren2014_equivalent_stress(upper_particles, upper_mat.ShearModulus(),
                                           upper_mat.BulkModulus());
            // Push host-filled stress so VtpCK prepareForOutput does not wipe it with device zeros.
            {
                const ParallelDevicePolicy pol{};
                lower_particles.getVariableByName<Real>("Ren2014EquivalentStress")->finalizeLoadIn(pol);
                upper_particles.getVariableByName<Real>("Ren2014EquivalentStress")->finalizeLoadIn(pol);
            }
            if (ite == std::numeric_limits<size_t>::max())
                recorder.writeToFile();
            else
                recorder.writeToFile(ite);
        };

        update_cfg.exec();
        corr_l.exec();
        corr_u.exec();
        write_vtp(0);

        size_t acoustic_steps = 0;
        auto host_pairwise_damp = [&](Real dt) {
            if (!(dt > 0.0) || !std::isfinite(dt))
                return;
            const ParallelDevicePolicy pol{};
            lower_particles.getVariableByName<Vecd>("Position")->prepareForOutput(pol);
            upper_particles.getVariableByName<Vecd>("Position")->prepareForOutput(pol);
            lower_particles.getVariableByName<Vecd>("Velocity")->prepareForOutput(pol);
            upper_particles.getVariableByName<Vecd>("Velocity")->prepareForOutput(pol);
            // Classic host CLL (cell_index_lists_), independent of CK device CLL buffers.
            lower_body.updateCellLinkedList();
            upper_body.updateCellLinkedList();
            lower_inner_host.updateConfiguration();
            upper_inner_host.updateConfiguration();
            damp_l_host.exec(dt);
            damp_u_host.exec(dt);
            lower_particles.getVariableByName<Vecd>("Velocity")->finalizeLoadIn(pol);
            upper_particles.getVariableByName<Vecd>("Velocity")->finalizeLoadIn(pol);
        };

        auto one_acoustic = [&](Real dt) {
            if (opt.apply_load)
                apply_force.exec();
            factor_l.exec();
            force_l.exec();
            factor_u.exec();
            force_u.exec();
            step1_l.exec(dt);
            support.exec();
            decay_l.exec(dt);
            support.exec();
            step2_l.exec(dt);
            step1_u.exec(dt);
            decay_u.exec(dt);
            step2_u.exec(dt);
            update_cfg.exec();
            ++acoustic_steps;
            if (opt.host_damp_stride > 0 && (acoustic_steps % static_cast<size_t>(opt.host_damp_stride) == 0))
                host_pairwise_damp(dt);
        };

        auto measure_vmax = [&]() -> Real {
            sync_ren2014_host(lower_particles);
            sync_ren2014_host(upper_particles);
            Vecd *vel_l = lower_particles.getVariableDataByName<Vecd>("Velocity");
            Vecd *vel_u = upper_particles.getVariableDataByName<Vecd>("Velocity");
            Real vmax = 0.0;
            for (size_t i = 0; i < lower_particles.TotalRealParticles(); ++i)
                vmax = std::max(vmax, vel_l[i].norm());
            for (size_t i = 0; i < upper_particles.TotalRealParticles(); ++i)
                vmax = std::max(vmax, vel_u[i].norm());
            return vmax;
        };

        auto relax_to_balance = [&](Real cap, const char *label) -> int {
            Real pseudo = 0.0;
            size_t steps = 0;
            int calm = 0;
            Real vmax_last = 0.0;
            while (pseudo < cap)
            {
                const Real dt = std::min(dt_l.exec(), dt_u.exec());
                if (!(dt > 0.0) || !std::isfinite(dt))
                {
                    std::cerr << "relax_abort reason=bad_dt label=" << label << " pseudo=" << pseudo
                              << " steps=" << steps << " dt=" << dt << " vmax_last=" << vmax_last << "\n";
                    return 1;
                }
                one_acoustic(dt);
                pseudo += dt;
                ++steps;
                if (steps % 500 == 0)
                {
                    corr_l.exec();
                    corr_u.exec();
                    vmax_last = measure_vmax();
                    std::cout << label << " N=" << steps << " pseudo=" << pseudo << " vmax=" << vmax_last
                              << "\n";
                    if (pseudo > 2.5e-4 && vmax_last < opt.vmax_calm)
                    {
                        ++calm;
                        if (calm >= 3)
                            return 0;
                    }
                    else
                        calm = 0;
                    // Extra host-split damping when kinetic energy rises (no velocity wipe).
                    if (vmax_last > 1.0)
                    {
                        const Real dt_boost = std::min(dt_l.exec(), dt_u.exec());
                        for (int s = 0; s < 5; ++s)
                            host_pairwise_damp(dt_boost > 0.0 ? dt_boost : Real(1.0e-7));
                    }
                }
            }
            vmax_last = measure_vmax();
            std::cout << label << " cap_reached pseudo=" << pseudo << " vmax=" << vmax_last << "\n";
            // Best-effort at cap. Do not abort on residual vmax.
            return 0;
        };

        TickCount t0 = TickCount::now();
        if (relax_to_balance(opt.relax_cap, "elastic") != 0)
            return 1;
        sync_ren2014_host(lower_particles);
        sync_ren2014_host(upper_particles);
        Vecd *pos_l = lower_particles.ParticlePositions();
        Vecd *pos_u = upper_particles.ParticlePositions();
        const spdb::DualSlotCavityMetrics cav_elastic =
            spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                           upper_particles.TotalRealParticles(), dp);
        std::cout << "elastic_balance cavity_m=" << cav_elastic.center_net_height
                  << " floor_n=" << cav_elastic.floor_n << " roof_n=" << cav_elastic.roof_n
                  << " d_cavity_m=" << (cav_elastic.center_net_height - cav_initial.center_net_height)
                  << "\n";
        write_vtp();
        write_status("elastic_balance", 0, opt.end_time_s, demo_T0_C, 0, 0);

        if (opt.mode == CaseMode::F1 || !(opt.end_time_s > 0.0))
        {
            std::ofstream out("spdb_program_temp_sycl_F1_metrics.txt");
            out << std::setprecision(10);
            out << "case=test_3d_spdb_programmed_temperature_sycl\nmode=F1\nbackend=SYCL\n";
            out << "material_source=demo_assumption\n";
            out << "force_N=" << force_total << "\n";
            out << "initial_cavity_net_height_m=" << cav_initial.center_net_height << "\n";
            out << "cavity_net_height_m=" << cav_elastic.center_net_height << "\n";
            out << "delta_cavity_from_initial_m="
                << (cav_elastic.center_net_height - cav_initial.center_net_height) << "\n";
            out << "floor_n=" << cav_elastic.floor_n << "\nroof_n=" << cav_elastic.roof_n << "\n";
            out << "wall_time_s=" << (TickCount::now() - t0).seconds() << "\n";
            write_status("finished", 0, 0, demo_T0_C, 0, 1);
            return 0;
        }

        const auto landings =
            build_landing_times(prog, opt.end_time_s, opt.macro_dt, opt.max_macro_dt);
        int aborted = 0;
        Real next_vtp_t = opt.vtp_every > 0.0 ? opt.vtp_every : opt.end_time_s + 1.0;
        for (size_t k = 1; k < landings.size(); ++k)
        {
            const Real t_target = landings[k];
            const Real dt_phys = t_target - physical_clock;
            if (!(dt_phys > 0.0))
                continue;
            const Real T_C = opt.use_program_temperature ? prog.temperature_C(t_target) : demo_T0_C;
            const Real T_K = spdb::celsius_to_kelvin(T_C);
            const Real gate =
                opt.enable_inelastic ? spdb::nonelastic_gate_weight_C(T_C) : 0.0;
            const Real stretch_new =
                opt.enable_thermal ? spdb::thermal_stretch_from_C(T_C, demo_alpha, demo_T0_C)
                                   : thermal_stretch;
            const Real factor = stretch_new / thermal_stretch;

            sync_ren2014_host(lower_particles);
            sync_ren2014_host(upper_particles);
            if (opt.enable_thermal && std::abs(factor - 1.0) > 1.0e-15)
            {
                // Split large isotropic kicks so contact does not go unstable in one jump.
                const Real max_step = 2.0e-4; // ~α·17 K at demo alpha
                const int n_sub = std::max(1, static_cast<int>(std::ceil(std::abs(std::log(factor)) / max_step)));
                const Real factor_sub = std::pow(factor, Real(1.0) / Real(n_sub));
                for (int s = 0; s < n_sub; ++s)
                {
                    spdb::apply_isotropic_thermal_factor(lower_particles, factor_sub,
                                                        Vecd(0.0, 0.0, -specimen_half_stack));
                    spdb::apply_isotropic_thermal_factor(upper_particles, factor_sub, Vecd(0.0, 0.0, 0.0));
                }
                thermal_stretch = stretch_new;
            }
            spdb::assign_uniform_temperature(lower_particles, T_K);
            spdb::assign_uniform_temperature(upper_particles, T_K);

            if (opt.enable_inelastic && gate > 0.0)
            {
                const auto c_l = commit_ren2014_physical_increment(
                    lower_particles, lower_mat.ShearModulus(), T_K, dt_phys, gate);
                const auto c_u = commit_ren2014_physical_increment(
                    upper_particles, upper_mat.ShearModulus(), T_K, dt_phys, gate);
                if (c_l.accepted != 1 || c_u.accepted != 1)
                {
                    std::cerr << "inelastic_abort at t_target=" << t_target << " T_C=" << T_C << "\n";
                    aborted = 1;
                    break;
                }
            }
            // Do NOT zero velocities here: dumping KE under unbalanced contact caused early vmax blow-ups.
            push_ren2014_device(lower_particles);
            push_ren2014_device(upper_particles);

            // Calm thermal kick with a few host-split damp passes before acoustic post-relax.
            {
                const Real dt0 = std::min(dt_l.exec(), dt_u.exec());
                for (int s = 0; s < 8; ++s)
                    host_pairwise_damp(dt0 > 0.0 ? dt0 : Real(1.0e-7));
            }

            // Slightly more pseudo-time after thermal map updates than elastic preload.
            const Real post_cap = std::max(opt.relax_cap, Real(2.5e-3));
            if (relax_to_balance(post_cap, "post") != 0)
            {
                std::cerr << "relax_abort at t_target=" << t_target << " T_C=" << T_C << " gate=" << gate
                          << "\n";
                aborted = 1;
                break;
            }
            physical_clock = t_target;
            sync_ren2014_host(lower_particles);
            sync_ren2014_host(upper_particles);
            const Real cav =
                spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                               upper_particles.TotalRealParticles(), dp)
                    .center_net_height;
            std::cout << "accepted t=" << physical_clock << " T_C=" << T_C << " gate=" << gate
                      << " cavity_m=" << cav << " mean_inl=" << mean_inelastic_host(lower_particles)
                      << "\n";
            write_status("accepted", physical_clock, opt.end_time_s, T_C, gate, 0);
            if (opt.vtp_every > 0.0 && physical_clock + 1.0e-9 >= next_vtp_t)
            {
                write_vtp();
                while (next_vtp_t <= physical_clock + 1.0e-9)
                    next_vtp_t += opt.vtp_every;
            }
        }

        sync_ren2014_host(lower_particles);
        sync_ren2014_host(upper_particles);
        write_vtp(); // final state for stress/strain ParaView
        const Real wall = (TickCount::now() - t0).seconds();
        const auto cav_final =
            spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                           upper_particles.TotalRealParticles(), dp);
        char name[192];
        std::snprintf(name, sizeof(name), "spdb_program_temp_sycl_%s_dp%.4g_end%.4g_metrics.txt",
                      mode_name(opt.mode), dp, opt.end_time_s);
        std::ofstream out(name);
        out << std::setprecision(10);
        const Real T_end = prog.temperature_C(physical_clock);
        const Real gate_end =
            opt.enable_inelastic ? spdb::nonelastic_gate_weight_C(T_end) : 0.0;
        out << "case=test_3d_spdb_programmed_temperature_sycl\n";
        out << "mode=" << mode_name(opt.mode) << "\nbackend=SYCL_GPU\n";
        out << "contact=ContactPenaltyForceCK_e_ij\n";
        out << "material_source=demo_assumption\n";
        out << "physical_s=" << physical_clock << "\nend_s_requested=" << opt.end_time_s << "\n";
        out << "aborted=" << aborted << "\n";
        out << "temperature_C=" << T_end << "\ngate=" << gate_end << "\n";
        out << "initial_cavity_net_height_m=" << cav_initial.center_net_height << "\n";
        out << "elastic_cavity_net_height_m=" << cav_elastic.center_net_height << "\n";
        out << "final_cavity_net_height_m=" << cav_final.center_net_height << "\n";
        out << "delta_cavity_elastic_m="
            << (cav_elastic.center_net_height - cav_initial.center_net_height) << "\n";
        out << "delta_cavity_program_m="
            << (cav_final.center_net_height - cav_elastic.center_net_height) << "\n";
        out << "floor_n=" << cav_final.floor_n << "\nroof_n=" << cav_final.roof_n << "\n";
        out << "stack_height_m=" << cav_final.stack_height << "\n";
        out << "mean_inelastic_lower=" << mean_inelastic_host(lower_particles) << "\n";
        out << "mean_inelastic_upper=" << mean_inelastic_host(upper_particles) << "\n";
        out << "wall_time_s=" << wall << "\n";
        write_status(aborted ? "aborted" : "finished", physical_clock, opt.end_time_s, T_end, gate_end,
                     aborted ? 0 : 1);
        std::cout << "done mode=" << mode_name(opt.mode) << " physical_s=" << physical_clock
                  << " aborted=" << aborted << " wall_s=" << wall << "\n";
        return aborted;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "exception: " << ex.what() << "\n";
        return 1;
    }
}
