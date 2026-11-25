/**
 * @file spdb_programmed_temperature.cpp
 * @brief Dual-slot STL specimen with uniform programmed temperature + 5 MPa demo.
 * @details Extra-test only. Uses ContactForce + ContactScaledRen2014.
 *          Thermal expansion and 850–900 C gate are demo_assumption.
 *          NOT heat transfer, NOT electric/joule, NOT graphite, NOT welding.
 */
#include "sphinxsys.h"

#include "spdb_contact_scaled.h"
#include "spdb_dual_slot_metrics.h"
#include "spdb_ren2014_solid.h"
#include "spdb_temperature_program.h"
#include "spdb_thermal_expansion.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
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
// demo_assumption fallback from package demo_config.json — NOT Ren2014, NOT measured.
constexpr Real demo_E = 120.0e9;
constexpr Real demo_nu = 0.30;
constexpr Real demo_rho = 7700.0;
constexpr Real demo_alpha = 1.2e-5;
constexpr Real demo_T0_K = 293.15;
constexpr Real demo_T0_C = 20.0;
constexpr Real cylinder_radius = 0.020;
constexpr Real specimen_half_stack = 0.020;
constexpr Real nominal_pressure = 5.0e6;
constexpr Real physical_viscosity = 5.0e6;
constexpr Real mesh_volume_ref = 2.38170881694208e-5;

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

Real kinetic_energy(BaseParticles &particles)
{
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real sum = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        sum += 0.5 * mass[i] * vel[i].squaredNorm();
    return sum;
}

Real max_speed(BaseParticles &particles)
{
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real vmax = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        vmax = std::max(vmax, vel[i].norm());
    return vmax;
}

Real sum_field_z(BaseParticles &particles, const char *name)
{
    Vecd *field = particles.getVariableDataByName<Vecd>(name);
    Real sum = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        sum += field[i][2];
    return sum;
}

Real particle_volume_sum(BaseParticles &particles)
{
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        sum += vol[i];
    return sum;
}

Real mean_inelastic(BaseParticles &particles)
{
    Real *eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real sum = 0.0;
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
        sum += eq[i];
    return n > 0 ? sum / static_cast<Real>(n) : 0.0;
}

Real max_inelastic(BaseParticles &particles)
{
    Real *eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real m = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        m = std::max(m, eq[i]);
    return m;
}

Real min_inelastic(BaseParticles &particles)
{
    Real *eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real m = std::numeric_limits<Real>::max();
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        m = std::min(m, eq[i]);
    return particles.TotalRealParticles() > 0 ? m : 0.0;
}

struct InterfaceGap
{
    size_t pairs = 0;
    Real max_penetration = 0.0;
};

InterfaceGap measure_interface_gap(const IndexVector &lower_ids, Vecd *pos_l, const IndexVector &upper_ids,
                                   Vecd *pos_u, Real dp)
{
    InterfaceGap gap;
    for (size_t i : lower_ids)
    {
        Real best_r = std::numeric_limits<Real>::max();
        size_t j_best = upper_ids.empty() ? 0 : upper_ids.front();
        for (size_t j : upper_ids)
        {
            const Real dx = pos_l[i][0] - pos_u[j][0];
            const Real dy = pos_l[i][1] - pos_u[j][1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            if (r < best_r)
            {
                best_r = r;
                j_best = j;
            }
        }
        if (best_r > 0.6 * dp)
            continue;
        const Real surface = (pos_u[j_best][2] - pos_l[i][2]) - dp;
        gap.max_penetration = std::max(gap.max_penetration, std::max(Real(0.0), -surface));
        ++gap.pairs;
    }
    return gap;
}

struct Snap
{
    std::vector<Vecd> pos, vel;
    std::vector<Matd> deformation, inverse_plastic;
    std::vector<Real> density, inelastic, range_flag, temperature;
};

Snap capture(BaseParticles &p)
{
    const size_t n = p.TotalRealParticles();
    Snap s;
    Vecd *pos = p.ParticlePositions();
    Vecd *vel = p.getVariableDataByName<Vecd>("Velocity");
    Matd *F = p.getVariableDataByName<Matd>("DeformationGradient");
    Matd *Cp = p.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *rho = p.getVariableDataByName<Real>("Density");
    Real *eq = p.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *flag = p.getVariableDataByName<Real>("ArrheniusRangeFlag");
    Real *T = p.getVariableDataByName<Real>("Temperature");
    s.pos.assign(pos, pos + n);
    s.vel.assign(vel, vel + n);
    s.deformation.assign(F, F + n);
    s.inverse_plastic.assign(Cp, Cp + n);
    s.density.assign(rho, rho + n);
    s.inelastic.assign(eq, eq + n);
    s.range_flag.assign(flag, flag + n);
    s.temperature.assign(T, T + n);
    return s;
}

void restore(BaseParticles &p, const Snap &s)
{
    const size_t n = p.TotalRealParticles();
    Vecd *pos = p.ParticlePositions();
    Vecd *vel = p.getVariableDataByName<Vecd>("Velocity");
    Matd *F = p.getVariableDataByName<Matd>("DeformationGradient");
    Matd *Cp = p.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *rho = p.getVariableDataByName<Real>("Density");
    Real *eq = p.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *flag = p.getVariableDataByName<Real>("ArrheniusRangeFlag");
    Real *T = p.getVariableDataByName<Real>("Temperature");
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = s.pos[i];
        vel[i] = s.vel[i];
        F[i] = s.deformation[i];
        Cp[i] = s.inverse_plastic[i];
        rho[i] = s.density[i];
        eq[i] = s.inelastic[i];
        flag[i] = s.range_flag[i];
        T[i] = s.temperature[i];
    }
}

void write_status(const char *phase, Real physical_s, Real end_s, Real T_C, Real gate, Real ke,
                  Real max_pen, Real cavity, Real mean_inl, Real mean_inu, int aborted = 0, int finished = 0)
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    char iso[64];
    std::strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&tt));
    std::ofstream out("hold_status.txt");
    out << std::setprecision(10);
    out << "updated_at=" << iso << "\nphase=" << phase << "\nphysical_s=" << physical_s
        << "\nend_s=" << end_s << "\nprogress_frac=" << (end_s > 0 ? physical_s / end_s : 0.0)
        << "\ntemperature_C=" << T_C << "\ngate=" << gate << "\nke_J=" << ke
        << "\nmax_penetration_m=" << max_pen << "\ncavity_net_height_m=" << cavity
        << "\nmean_inelastic_lower=" << mean_inl << "\nmean_inelastic_upper=" << mean_inu
        << "\naborted=" << aborted << "\nfinished=" << finished << "\n";
    std::ofstream beat("hold_heartbeat");
    beat << iso << " t=" << physical_s << " T_C=" << T_C << " phase=" << phase << "\n";
}

enum class CaseMode
{
    GeometryOnly,
    TemperatureUnit,
    F0,
    F1,
    F2,
    F3
};

struct Options
{
    CaseMode mode = CaseMode::F3;
    Real dp = 0.001;
    Real end_time_s = 1800.0;
    Real macro_dt = 10.0;
    Real max_macro_dt = 30.0;
    Real stiffness_scale = 1.0;
    Real relax_cap = 1.5e-3;
    Real gate_off_C = 850.0;
    Real gate_on_C = 900.0;
    bool hard_gate_at_900 = false; // sensitivity: off below 900, on above
    bool enable_thermal = true;
    bool enable_inelastic = true;
    bool apply_load = true;
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
        else if (a == "--temperature_unit")
            o.mode = CaseMode::TemperatureUnit;
        else if (a == "--F0")
        {
            o.mode = CaseMode::F0;
            o.apply_load = false;
            o.enable_inelastic = false;
            o.enable_thermal = true;
            o.use_program_temperature = true;
        }
        else if (a == "--F1")
        {
            o.mode = CaseMode::F1;
            o.apply_load = true;
            o.enable_inelastic = false;
            o.enable_thermal = false;
            o.use_program_temperature = false;
            o.end_time_s = 0.0; // elastic balance only
        }
        else if (a == "--F2")
        {
            o.mode = CaseMode::F2;
            o.apply_load = true;
            o.enable_inelastic = false;
            o.enable_thermal = true;
            o.use_program_temperature = true;
        }
        else if (a == "--F3")
        {
            o.mode = CaseMode::F3;
            o.apply_load = true;
            o.enable_inelastic = true;
            o.enable_thermal = true;
            o.use_program_temperature = true;
        }
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
        else if (a == "--hard_gate_900")
            o.hard_gate_at_900 = true;
    }
    return o;
}

const char *mode_name(CaseMode m)
{
    switch (m)
    {
    case CaseMode::GeometryOnly:
        return "geometry_only";
    case CaseMode::TemperatureUnit:
        return "temperature_unit";
    case CaseMode::F0:
        return "F0_free_thermal";
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
                                      Real max_dt, Real gate_off_C, Real gate_on_C)
{
    std::set<Real> marks;
    marks.insert(0.0);
    marks.insert(end_s);
    for (const auto &k : prog.knots())
        if (k.first <= end_s + 1.0e-12)
            marks.insert(k.first);
    for (Real t : prog.crossing_times_C(gate_off_C))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    for (Real t : prog.crossing_times_C(gate_on_C))
        if (t > 0.0 && t < end_s)
            marks.insert(t);
    // Fill with macro steps but never exceed max_dt; land exactly on marks.
    std::vector<Real> sorted(marks.begin(), marks.end());
    std::vector<Real> out;
    out.push_back(0.0);
    for (size_t i = 1; i < sorted.size(); ++i)
    {
        Real a = sorted[i - 1];
        const Real b = sorted[i];
        while (b - a > max_dt + 1.0e-12)
        {
            Real step = std::min(macro_dt, max_dt);
            // Prefer landing so remaining is not a tiny sliver.
            if (b - (a + step) < 0.5 * macro_dt && b - (a + step) > 1.0e-12)
                step = 0.5 * (b - a);
            a += step;
            out.push_back(a);
        }
        out.push_back(b);
    }
    // Unique / sorted
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(),
                          [](Real x, Real y) { return std::abs(x - y) < 1.0e-9; }),
              out.end());
    return out;
}

Real gate_weight(const Options &opt, Real T_C)
{
    if (opt.hard_gate_at_900)
        return T_C >= 900.0 ? 1.0 : 0.0;
    return spdb::nonelastic_gate_weight_C(T_C, opt.gate_off_C, opt.gate_on_C);
}

int run_temperature_unit()
{
    spdb::TemperatureProgram prog;
    prog.load_csv("input/temperature_program.csv");
    if (!spdb::temperature_program_unit_ok(prog))
    {
        std::cerr << "temperature_program unit check FAILED\n";
        return 1;
    }
    const double g850 = spdb::nonelastic_gate_weight_C(850.0);
    const double g875 = spdb::nonelastic_gate_weight_C(875.0);
    const double g900 = spdb::nonelastic_gate_weight_C(900.0);
    if (g850 != 0.0 || std::abs(g875 - 0.5) > 1.0e-12 || g900 != 1.0)
    {
        std::cerr << "gate unit check FAILED g850=" << g850 << " g875=" << g875 << " g900=" << g900 << "\n";
        return 1;
    }
    const double eth = demo_alpha * (950.0 - 20.0);
    std::cout << "temperature_unit OK\n"
              << "  T(150)_C=" << prog.temperature_C(150.0) << "\n"
              << "  free_thermal_strain_950C=" << eth << " (expect 0.01116)\n"
              << "  10mm_growth_m=" << 0.010 * eth << " 20mm_growth_m=" << 0.020 * eth << "\n";
    return std::abs(eth - 0.01116) < 1.0e-12 ? 0 : 1;
}

int run_case(Options opt)
{
    if (opt.mode == CaseMode::TemperatureUnit)
        return run_temperature_unit();

    spdb::TemperatureProgram prog;
    prog.load_csv("input/temperature_program.csv");
    if (!spdb::temperature_program_unit_ok(prog))
    {
        std::cerr << "temperature program failed unit check\n";
        return 1;
    }

    const Real dp = opt.dp;
    const Real force_total = opt.apply_load ? nominal_force() : 0.0;
    const Real ramp_time = 2.5e-4;
    const Real slab = 0.55 * dp;
    Real thermal_stretch = 1.0;
    Real physical_clock = 0.0;

    std::cout << "SPDB programmed_temperature\n"
              << "  mode=" << mode_name(opt.mode) << "\n"
              << "  material=demo_assumption E_Pa=" << demo_E << " nu=" << demo_nu
              << " rho=" << demo_rho << " alpha=" << demo_alpha << "\n"
              << "  NOT measured 20Cr13 full-T properties; Ren2014 rate with demo gate.\n"
              << "  excluded=heat_equation,electric,joule,graphite,welding,phase_change\n"
              << "  force_N=" << force_total << " dp_m=" << dp << " end_s=" << opt.end_time_s << "\n"
              << "  stl=input/lower_m.stl + input/upper_m.stl scale=1 (already meters)\n";

    SPHSystem sph_system(domain_bounds(dp), dp);
    SolidBody lower_body(sph_system, makeShared<LowerFromSTL>("lower_m.stl"));
    lower_body.defineAdaptationRatios(1.15, 1.0);
    ContactScaledRen2014 &lower_mat = lower_body.defineMatterMaterial<ContactScaledRen2014>(
        demo_rho, demo_E, demo_nu, demo_T0_K, false, opt.stiffness_scale);
    lower_body.defineBodyLevelSetShape().correctLevelSetSign();
    lower_body.generateParticles<BaseParticles, Lattice>();

    SolidBody upper_body(sph_system, makeShared<UpperFromSTL>("upper_m.stl"));
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
    lower_particles.addVariableToWrite<Real>("Temperature");
    upper_particles.addVariableToWrite<Real>("Temperature");
    lower_particles.addVariableToWrite<Real>("EquivalentInelasticStrain");
    upper_particles.addVariableToWrite<Real>("EquivalentInelasticStrain");

    const Real vol_l = particle_volume_sum(lower_particles);
    const Real vol_u = particle_volume_sum(upper_particles);
    const Real vol_err_l = std::abs(vol_l - mesh_volume_ref) / mesh_volume_ref;
    const Real vol_err_u = std::abs(vol_u - mesh_volume_ref) / mesh_volume_ref;
    std::cout << "  lower_particles=" << lower_particles.TotalRealParticles()
              << " upper_particles=" << upper_particles.TotalRealParticles() << "\n"
              << "  particle_volume_lower_m3=" << vol_l << " rel_err=" << vol_err_l << "\n"
              << "  particle_volume_upper_m3=" << vol_u << " rel_err=" << vol_err_u << "\n"
              << "  mesh_volume_ref_m3=" << mesh_volume_ref << "\n";

    Vecd *pos_l = lower_particles.ParticlePositions();
    Vecd *pos_u = upper_particles.ParticlePositions();
    const spdb::DualSlotCavityMetrics cav0 =
        spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                       upper_particles.TotalRealParticles(), dp);
    std::cout << "  initial_cavity_net_height_m=" << cav0.center_net_height
              << " floor_n=" << cav0.floor_n << " roof_n=" << cav0.roof_n
              << " stack_height_m=" << cav0.stack_height << "\n";

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.writeToFile(0);

    if (opt.mode == CaseMode::GeometryOnly)
    {
        std::ofstream out("spdb_program_temp_geometry_metrics.txt");
        out << std::setprecision(10);
        out << "case=test_3d_spdb_programmed_temperature\nmode=geometry_only\n";
        out << "material_source=demo_assumption\n";
        out << "lower_particles=" << lower_particles.TotalRealParticles() << "\n";
        out << "upper_particles=" << upper_particles.TotalRealParticles() << "\n";
        out << "particle_volume_lower_m3=" << vol_l << "\n";
        out << "particle_volume_upper_m3=" << vol_u << "\n";
        out << "mesh_volume_ref_m3=" << mesh_volume_ref << "\n";
        out << "volume_rel_err_lower=" << vol_err_l << "\n";
        out << "volume_rel_err_upper=" << vol_err_u << "\n";
        out << "cavity_net_height_m=" << cav0.center_net_height << "\n";
        out << "cavity_width_x_m=" << cav0.width_x << "\n";
        out << "cavity_width_y_m=" << cav0.width_y << "\n";
        out << "cavity_volume_claimed=0\n";
        out << "stl_scale=1\n";
        std::cout << "geometry_only done. volume_rel_err_lower=" << vol_err_l
                  << " volume_rel_err_upper=" << vol_err_u << "\n";
        // Soft check: lattice volume often differs more than surface mesh; report only.
        return 0;
    }

    InnerRelation lower_inner(lower_body);
    InnerRelation upper_inner(upper_body);
    SurfaceContactRelation lower_upper(lower_body, {&upper_body});
    SurfaceContactRelation upper_lower(upper_body, {&lower_body});

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_l(lower_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_u(upper_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_l_1st(lower_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_l_2nd(lower_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_u_1st(upper_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_u_2nd(upper_inner);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor_l(lower_upper);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor_u(upper_lower);
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_l(lower_upper);
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_u(upper_lower);

    GeometricShapeBox bottom_shape(Transform(Vec3d(0.0, 0.0, -specimen_half_stack + 0.5 * dp)),
                                   Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "BottomAxialSupport");
    GeometricShapeBox top_shape(Transform(Vec3d(0.0, 0.0, specimen_half_stack - 0.5 * dp)),
                                Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "TopAxialLoad");
    BodyRegionByParticle bottom_part(lower_body, bottom_shape);
    BodyRegionByParticle top_part(upper_body, top_shape);
    if (bottom_part.body_part_particles_.empty())
    {
        std::cerr << "empty bottom support layer\n";
        return 1;
    }
    if (opt.apply_load && top_part.body_part_particles_.empty())
    {
        std::cerr << "empty top load layer\n";
        return 1;
    }

    SimpleDynamics<FixedInAxisDirection> constrain_bottom(bottom_part, Vec3d(1.0, 1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_com_l(lower_body, Vec3d(1.0, 1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_com_u(upper_body, Vec3d(1.0, 1.0, 0.0));
    const Vec3d force_per_particle =
        top_part.body_part_particles_.empty()
            ? Vec3d::Zero()
            : Vec3d(0.0, 0.0, -force_total / static_cast<Real>(top_part.body_part_particles_.size()));
    SimpleDynamics<solid_dynamics::ForceInBodyRegion> apply_force(top_part, force_per_particle, ramp_time);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        damp_l(0.2, lower_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        damp_u(0.2, upper_inner, "Velocity", physical_viscosity);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_l(lower_body, 0.2);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_u(upper_body, 0.2);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_l.exec();
    corrected_u.exec();

    IndexVector lower_face = spdb::dual_ids_near_z_outside_slot(
        pos_l, lower_particles.TotalRealParticles(), -0.5 * dp, 0.55 * dp);
    IndexVector upper_face = spdb::dual_ids_near_z_outside_slot(
        pos_u, upper_particles.TotalRealParticles(), 0.5 * dp, 0.55 * dp);
    if (lower_face.empty() || upper_face.empty())
    {
        std::cerr << "empty join-face ring lower_n=" << lower_face.size()
                  << " upper_n=" << upper_face.size() << "\n";
        return 1;
    }

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    auto one_acoustic = [&](Real dt) {
        if (opt.apply_load)
            apply_force.exec(dt);
        contact_factor_l.exec();
        contact_force_l.exec();
        contact_factor_u.exec();
        contact_force_u.exec();
        stress_l_1st.exec(dt);
        constrain_bottom.exec(dt);
        constrain_com_l.exec(dt);
        damp_l.exec(dt);
        constrain_bottom.exec(dt);
        constrain_com_l.exec(dt);
        stress_l_2nd.exec(dt);
        stress_u_1st.exec(dt);
        constrain_com_u.exec(dt);
        damp_u.exec(dt);
        constrain_com_u.exec(dt);
        stress_u_2nd.exec(dt);
        lower_body.updateCellLinkedList();
        upper_body.updateCellLinkedList();
        lower_inner.updateConfiguration();
        upper_inner.updateConfiguration();
        lower_upper.updateConfiguration();
        upper_lower.updateConfiguration();
    };

    auto relax_to_balance = [&](Real cap, const char *label) -> int {
        Real previous_pen = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
        int calm = 0;
        size_t steps = 0;
        Real pseudo = 0.0;
        while (pseudo < cap)
        {
            const Real dt = SMIN(dt_l.exec(), dt_u.exec());
            if (!(dt > 0.0))
                return 1;
            one_acoustic(dt);
            ++steps;
            pseudo += dt;
            if (std::strcmp(label, "elastic") == 0)
                physical_time = pseudo;
            if (max_speed(lower_particles) > 20.0 || max_speed(upper_particles) > 20.0)
            {
                std::cerr << label << " abort: speed\n";
                return 1;
            }
            if (steps % 500 == 0)
            {
                corrected_l.exec();
                corrected_u.exec();
                const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
                const Real ke = kinetic_energy(lower_particles) + kinetic_energy(upper_particles);
                const Real dpen = std::abs(gap.max_penetration - previous_pen);
                std::cout << label << " N=" << steps << " pseudo_s=" << pseudo
                          << " max_pen_m=" << gap.max_penetration << " ke_J=" << ke << "\n";
                if (pseudo > ramp_time && ke < 1.0e-6 && max_speed(lower_particles) < 2.0e-3 &&
                    max_speed(upper_particles) < 2.0e-3 && dpen < 1.0e-9)
                {
                    ++calm;
                    if (calm >= 3)
                        return 0;
                }
                else
                    calm = 0;
                previous_pen = gap.max_penetration;
            }
        }
        return 0;
    };

    std::ofstream csv("program_temperature_observables.csv");
    csv << "physical_s,temperature_C,temperature_K,gate,force_N,reaction_z_N,reaction_rel_err,"
           "cavity_net_height_m,min_clearance_m,width_x_m,width_y_m,lower_height_m,upper_height_m,"
           "stack_height_m,max_penetration_m,mean_inelastic_l,mean_inelastic_u,max_inelastic_l,"
           "max_inelastic_u,min_inelastic_l,T_particle_min_K,T_particle_max_K,thermal_stretch,"
           "wall_time_s\n";

    TickCount t_wall0 = TickCount::now();
    auto append_csv = [&](Real T_C, Real gate) {
        const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
        const spdb::DualSlotCavityMetrics cav =
            spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                           upper_particles.TotalRealParticles(), dp);
        Vecd *force_l = lower_particles.getVariableDataByName<Vecd>("Force");
        Vecd *prior_l = lower_particles.getVariableDataByName<Vecd>("ForcePrior");
        Real bottom_internal_z = 0.0;
        for (size_t i : bottom_part.body_part_particles_)
            bottom_internal_z += force_l[i][2] + prior_l[i][2];
        const Real reaction_z = -bottom_internal_z;
        const Real reaction_err =
            force_total > 0.0 ? std::abs(std::abs(reaction_z) - force_total) / force_total : 0.0;
        Real tmin = 0.0, tmax = 0.0;
        spdb::temperature_minmax(lower_particles, &tmin, &tmax);
        Real tmin_u = 0.0, tmax_u = 0.0;
        spdb::temperature_minmax(upper_particles, &tmin_u, &tmax_u);
        tmin = std::min(tmin, tmin_u);
        tmax = std::max(tmax, tmax_u);
        const Real wall = (TickCount::now() - t_wall0).seconds();
        csv << std::setprecision(10) << physical_clock << "," << T_C << "," << spdb::celsius_to_kelvin(T_C)
            << "," << gate << "," << force_total << "," << reaction_z << "," << reaction_err << ","
            << cav.center_net_height << "," << cav.min_clearance << "," << cav.width_x << "," << cav.width_y
            << "," << cav.lower_height << "," << cav.upper_height << "," << cav.stack_height << ","
            << gap.max_penetration << "," << mean_inelastic(lower_particles) << ","
            << mean_inelastic(upper_particles) << "," << max_inelastic(lower_particles) << ","
            << max_inelastic(upper_particles) << "," << min_inelastic(lower_particles) << "," << tmin << ","
            << tmax << "," << thermal_stretch << "," << wall << "\n";
        csv.flush();
    };

    // --- 20 C preload balance (or free geometry for F0) ---
    spdb::assign_uniform_temperature(lower_particles, demo_T0_K);
    spdb::assign_uniform_temperature(upper_particles, demo_T0_K);
    if (relax_to_balance(opt.relax_cap, "elastic") != 0)
        return 1;
    corrected_l.exec();
    corrected_u.exec();
    const InterfaceGap gap_elastic = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const spdb::DualSlotCavityMetrics cav_elastic =
        spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                       upper_particles.TotalRealParticles(), dp);
    std::cout << "elastic_balance max_pen_m=" << gap_elastic.max_penetration
              << " cavity_m=" << cav_elastic.center_net_height << "\n"
              << "  physical clock starts; loading-period inelastic omitted.\n";
    physical_clock = 0.0;
    append_csv(demo_T0_C, 0.0);
    write_states.writeToFile();
    write_status("elastic_balance", 0.0, opt.end_time_s, demo_T0_C, 0.0,
                 kinetic_energy(lower_particles) + kinetic_energy(upper_particles),
                 gap_elastic.max_penetration, cav_elastic.center_net_height, 0.0, 0.0);

    if (opt.mode == CaseMode::F1 || !(opt.end_time_s > 0.0))
    {
        std::ofstream out("spdb_program_temp_F1_metrics.txt");
        out << std::setprecision(10);
        out << "case=test_3d_spdb_programmed_temperature\nmode=F1\n";
        out << "material_source=demo_assumption\n";
        out << "temperature_C=20\nforce_N=" << force_total << "\n";
        out << "cavity_net_height_m=" << cav_elastic.center_net_height << "\n";
        out << "max_penetration_m=" << gap_elastic.max_penetration << "\n";
        out << "delta_cavity_from_initial_m=" << (cav_elastic.center_net_height - cav0.center_net_height)
            << "\n";
        write_status("finished", 0.0, 0.0, demo_T0_C, 0.0,
                     kinetic_energy(lower_particles) + kinetic_energy(upper_particles),
                     gap_elastic.max_penetration, cav_elastic.center_net_height, 0.0, 0.0, 0, 1);
        return 0;
    }

    const std::vector<Real> landings =
        build_landing_times(prog, opt.end_time_s, opt.macro_dt, opt.max_macro_dt, opt.gate_off_C, opt.gate_on_C);
    int aborted = 0;
    for (size_t k = 1; k < landings.size(); ++k)
    {
        const Real t_target = landings[k];
        const Real dt_phys = t_target - physical_clock;
        if (!(dt_phys > 0.0))
            continue;
        const Real T_C = opt.use_program_temperature ? prog.temperature_C(t_target) : demo_T0_C;
        const Real T_K = spdb::celsius_to_kelvin(T_C);
        const Real gate = opt.enable_inelastic ? gate_weight(opt, T_C) : 0.0;
        const Real stretch_new =
            opt.enable_thermal ? spdb::thermal_stretch_from_C(T_C, demo_alpha, demo_T0_C) : thermal_stretch;
        const Real factor = stretch_new / thermal_stretch;

        const Snap saved_l = capture(lower_particles);
        const Snap saved_u = capture(upper_particles);

        // Thermal: scale about supported bottom (lower) and join plane (upper).
        if (opt.enable_thermal && std::abs(factor - 1.0) > 1.0e-15)
        {
            spdb::apply_isotropic_thermal_factor(lower_particles, factor, Vecd(0.0, 0.0, -specimen_half_stack));
            spdb::apply_isotropic_thermal_factor(upper_particles, factor, Vecd(0.0, 0.0, 0.0));
            thermal_stretch = stretch_new;
        }
        spdb::assign_uniform_temperature(lower_particles, T_K);
        spdb::assign_uniform_temperature(upper_particles, T_K);

        if (opt.enable_inelastic && gate > 0.0)
        {
            const Ren2014PhysicalCommit c_l = commit_ren2014_physical_increment(
                lower_particles, lower_mat.ShearModulus(), T_K, dt_phys, gate);
            const Ren2014PhysicalCommit c_u = commit_ren2014_physical_increment(
                upper_particles, upper_mat.ShearModulus(), T_K, dt_phys, gate);
            if (c_l.accepted != 1 || c_u.accepted != 1)
            {
                restore(lower_particles, saved_l);
                restore(upper_particles, saved_u);
                thermal_stretch = spdb::thermal_stretch_from_C(
                    opt.use_program_temperature ? prog.temperature_C(physical_clock) : demo_T0_C, demo_alpha,
                    demo_T0_C);
                std::cerr << "commit rejected at t=" << physical_clock << " -> " << t_target << "\n";
                aborted = 1;
                break;
            }
        }

        physical_time = ramp_time;
        if (relax_to_balance(opt.relax_cap, "post") != 0)
        {
            restore(lower_particles, saved_l);
            restore(upper_particles, saved_u);
            aborted = 1;
            break;
        }
        physical_clock = t_target;
        append_csv(T_C, gate);
        write_status("accepted", physical_clock, opt.end_time_s, T_C, gate,
                     kinetic_energy(lower_particles) + kinetic_energy(upper_particles),
                     measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration,
                     spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                                    upper_particles.TotalRealParticles(), dp)
                         .center_net_height,
                     mean_inelastic(lower_particles), mean_inelastic(upper_particles));
        // Key frames: breakpoints and every ~150 s.
        bool key = false;
        for (const auto &kn : prog.knots())
            if (std::abs(kn.first - physical_clock) < 1.0e-9)
                key = true;
        if (static_cast<int>(std::round(physical_clock)) % 150 == 0)
            key = true;
        if (key)
            write_states.writeToFile();
        std::cout << "accepted t=" << physical_clock << " T_C=" << T_C << " gate=" << gate
                  << " stretch=" << thermal_stretch
                  << " cavity_m="
                  << spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                                    upper_particles.TotalRealParticles(), dp)
                         .center_net_height
                  << " mean_inl=" << mean_inelastic(lower_particles) << "\n";
    }

    const InterfaceGap gap_final = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const spdb::DualSlotCavityMetrics cav_final =
        spdb::measure_dual_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                       upper_particles.TotalRealParticles(), dp);
    const Real wall = (TickCount::now() - t_wall0).seconds();
    char metrics_name[192];
    std::snprintf(metrics_name, sizeof(metrics_name), "spdb_program_temp_%s_dp%.4g_end%.4g_metrics.txt",
                  mode_name(opt.mode), dp, opt.end_time_s);
    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_programmed_temperature\n";
    out << "mode=" << mode_name(opt.mode) << "\n";
    out << "material_source=demo_assumption\n";
    out << "material_note=NOT_measured_20Cr13_full_T; Ren2014_rate_with_demo_gate\n";
    out << "excluded_physics=heat_equation,electric,joule,graphite,welding,phase_change\n";
    out << "backend=CPU_build-cpu-icpx\n";
    out << "youngs_modulus_Pa=" << demo_E << "\nnu=" << demo_nu << "\nrho=" << demo_rho
        << "\nalpha_per_K=" << demo_alpha << "\n";
    out << "force_N=" << force_total << "\nnominal_pressure_Pa=" << nominal_pressure << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "physical_s=" << physical_clock << "\nend_s_requested=" << opt.end_time_s << "\n";
    out << "aborted=" << aborted << "\n";
    out << "initial_cavity_net_height_m=" << cav0.center_net_height << "\n";
    out << "elastic_cavity_net_height_m=" << cav_elastic.center_net_height << "\n";
    out << "final_cavity_net_height_m=" << cav_final.center_net_height << "\n";
    out << "delta_cavity_elastic_m=" << (cav_elastic.center_net_height - cav0.center_net_height) << "\n";
    out << "delta_cavity_program_m=" << (cav_final.center_net_height - cav_elastic.center_net_height) << "\n";
    out << "free_thermal_20mm_growth_expect_m=" << 0.020 * demo_alpha * (950.0 - 20.0) << "\n";
    out << "max_penetration_elastic_m=" << gap_elastic.max_penetration << "\n";
    out << "max_penetration_final_m=" << gap_final.max_penetration << "\n";
    out << "mean_inelastic_lower=" << mean_inelastic(lower_particles) << "\n";
    out << "mean_inelastic_upper=" << mean_inelastic(upper_particles) << "\n";
    out << "max_inelastic_lower=" << max_inelastic(lower_particles) << "\n";
    out << "max_inelastic_upper=" << max_inelastic(upper_particles) << "\n";
    out << "cavity_volume_claimed=0\n";
    out << "wall_time_s=" << wall << "\n";
    out << "observables_csv=program_temperature_observables.csv\n";
    write_status(aborted ? "aborted" : "finished", physical_clock, opt.end_time_s,
                 opt.use_program_temperature ? prog.temperature_C(physical_clock) : demo_T0_C,
                 opt.enable_inelastic ? gate_weight(opt, prog.temperature_C(physical_clock)) : 0.0,
                 kinetic_energy(lower_particles) + kinetic_energy(upper_particles), gap_final.max_penetration,
                 cav_final.center_net_height, mean_inelastic(lower_particles), mean_inelastic(upper_particles),
                 aborted, aborted ? 0 : 1);
    std::cout << "done mode=" << mode_name(opt.mode) << " physical_s=" << physical_clock
              << " aborted=" << aborted << " wall_s=" << wall << " metrics=" << metrics_name << "\n";
    return aborted;
}

} // namespace

int main(int ac, char *av[])
{
    Options opt = parse_args(ac, av);
    try
    {
        return run_case(opt);
    }
    catch (const std::exception &ex)
    {
        std::cerr << "exception: " << ex.what() << "\n";
        return 1;
    }
}
