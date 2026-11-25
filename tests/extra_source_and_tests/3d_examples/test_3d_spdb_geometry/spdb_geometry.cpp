/**
 * @file 	spdb_geometry.cpp
 * @brief 	SPDB P1: parametric geometry_only + no-slot elastic compression baseline.
 *
 * Case executable. Reusable sources live in
 * extra_src/shared (Arrhenius) and extra_src/for_3D_build (workpiece geometry).
 *
 * Usage:
 *   test_3d_spdb_geometry --geometry_only [--no_slot]
 *   test_3d_spdb_geometry --elastic_no_slot
 *   test_3d_spdb_geometry --elastic_axial
 *
 * Material is synthetic elastic verification data (source=synthetic).
 * Do NOT interpret results as 950 C 20Cr13 predictions.
 *   test_3d_spdb_geometry --material_point
 *
 * --material_point checks Ren 2014 hot-compression coefficients.
 * It is not a 950 C / 7200 s creep prediction.
 */
#include "spdb_arrhenius.h"
#include "spdb_ck_axial.h"
#include "spdb_contact_scaled.h"
#include "spdb_geometry.h"
#include "spdb_ren2014_solid.h"
#include "spdb_slot_hold.h"
#include "spdb_slot_metrics.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

using namespace SPH;

namespace
{
struct RunOptions
{
    bool geometry_only = false;
    bool elastic_no_slot = false;
    bool elastic_axial = false;
    bool elastic_contact = false;
    bool slot_contact = false;
    bool slot_hold = false;
    bool material_point = false;
    bool ck_hold = false;
    bool ck_contact = false;
    Real dp = 0.001;
    Real end_time = 0.0;
    Real contact_stiffness_scale = 1.0;
    Real physical_dt = 60.0;
    Real damping_factor = 0.125;
    Real contact_youngs = 0.0;
    int hold_steps = 1;
    SlotMode slot_mode = SlotMode::LowerOnly;
};

RunOptions parse_options(int ac, char *av[])
{
    RunOptions opt;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--geometry_only") == 0)
            opt.geometry_only = true;
        else if (std::strcmp(av[i], "--elastic_no_slot") == 0)
            opt.elastic_no_slot = true;
        else if (std::strcmp(av[i], "--elastic_axial") == 0)
            opt.elastic_axial = true;
        else if (std::strcmp(av[i], "--elastic_contact") == 0)
            opt.elastic_contact = true;
        else if (std::strcmp(av[i], "--slot_contact") == 0)
            opt.slot_contact = true;
        else if (std::strcmp(av[i], "--slot_hold") == 0)
            opt.slot_hold = true;
        else if (std::strcmp(av[i], "--material_point") == 0)
            opt.material_point = true;
        else if (std::strcmp(av[i], "--ck_hold") == 0)
            opt.ck_hold = true;
        else if (std::strcmp(av[i], "--ck_contact") == 0)
            opt.ck_contact = true;
        else if (std::strcmp(av[i], "--dp") == 0 && i + 1 < ac)
            opt.dp = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--end_time") == 0 && i + 1 < ac)
            opt.end_time = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--contact_stiffness_scale") == 0 && i + 1 < ac)
            opt.contact_stiffness_scale = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--contact_youngs") == 0 && i + 1 < ac)
            opt.contact_youngs = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--physical_dt") == 0 && i + 1 < ac)
            opt.physical_dt = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--damping_factor") == 0 && i + 1 < ac)
            opt.damping_factor = std::strtod(av[++i], nullptr);
        else if (std::strcmp(av[i], "--hold_steps") == 0 && i + 1 < ac)
            opt.hold_steps = std::atoi(av[++i]);
        else if (std::strcmp(av[i], "--no_slot") == 0)
            opt.slot_mode = SlotMode::None;
        else if (std::strcmp(av[i], "--both_slots") == 0)
            opt.slot_mode = SlotMode::Both;
        else if (std::strcmp(av[i], "--help") == 0)
        {
            std::cout << "SPDB geometry case options:\n"
                      << "  --geometry_only     Generate particles, report metrics, write VTK, exit\n"
                      << "  --no_slot           No blind slot (default for geometry is lower-only slot)\n"
                      << "  --both_slots        Upper and lower slots (not a first-round default run)\n"
                      << "  --elastic_no_slot   Short synthetic elastic compression (no slot)\n"
                      << "  --elastic_axial     No-slot uniaxial 10 MPa vs Hooke (synthetic E)\n"
                      << "  --elastic_contact   Two-body flush contact, no slot\n"
                      << "  --slot_contact      Lower blind-slot contact, demo E=100 GPa by default\n"
                      << "  --slot_hold         Slot contact + Ren2014 physical hold (ContactForce path)\n"
                      << "  --material_point    Ren 2014 Arrhenius material-point checks\n"
                      << "  --ck_hold           No-slot CK hold: elastic relax, then physical Δt\n"
                      << "  --ck_contact        Two no-slot cylinders, CK RepulsionForceCK (not for flush)\n"
                      << "  --physical_dt <s>   One accepted hold increment (default 60, not 7200)\n"
                      << "  --hold_steps <n>    How many accepted increments (default 1)\n"
                      << "  --damping_factor <x> StructureIntegration numerical damping (default 0.125)\n"
                      << "  --dp <m>            Particle spacing (default 0.001)\n"
                      << "  --end_time <s>      Override the mode's end / pseudo-time cap\n"
                      << "  --contact_stiffness_scale <x>\n"
                      << "                      Multiply rho*c^2 contact stiffness (default 1)\n"
                      << "  --contact_youngs <Pa>\n"
                      << "                      Elastic-contact modulus (default 2 GPa synthetic)\n";
            std::exit(0);
        }
    }
    if (!opt.geometry_only && !opt.elastic_no_slot && !opt.elastic_axial && !opt.elastic_contact &&
        !opt.slot_contact && !opt.slot_hold && !opt.material_point && !opt.ck_hold && !opt.ck_contact)
        opt.geometry_only = true; // safe default for first-round smoke
    if (opt.elastic_no_slot || opt.elastic_axial || opt.elastic_contact || opt.ck_hold || opt.ck_contact)
        opt.slot_mode = SlotMode::None;
    if (opt.slot_contact || opt.slot_hold)
        opt.slot_mode = SlotMode::LowerOnly;
    return opt;
}

Real mean_coordinate(const IndexVector &ids, Vecd *pos, int axis)
{
    Real sum = 0.0;
    for (size_t i : ids)
        sum += pos[i][axis];
    return ids.empty() ? 0.0 : sum / static_cast<Real>(ids.size());
}

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

//----------------------------------------------------------------------
// No-slot uniaxial compression against the infinitesimal Hooke solution.
// One cylinder, bottom layer fixed in z only, nominal force on the top
// layer, lateral surfaces free. Synthetic E = Youngs_modulus_axial.
//----------------------------------------------------------------------
int run_elastic_axial(Real dp, Real end_time)
{
    particle_spacing_ref = dp;
    resolution_ref = dp;
    const Real force_total = nominal_axial_force();
    const Real eps_z_analytical = -nominal_pressure / Youngs_modulus_axial;
    const Real eps_r_analytical = -poisson * eps_z_analytical;
    const Real ramp_time = 1.0e-3;
    const Real output_period = 1.0e-3;
    const Real interior_inset = std::max(4.0 * dp, 0.004);
    const Real core_inset = std::max(3.0 * dp, 0.003);
    const Real radial_halfband = std::max(2.0 * dp, 0.002);

    std::cout << "SPDB elastic_axial\n"
              << "  NOTE: synthetic LinearElasticSolid, E=" << Youngs_modulus_axial
              << " Pa, nu=" << poisson << ". NOT 950C 20Cr13.\n"
              << "  nominal_pressure_Pa=" << nominal_pressure
              << "  nominal_force_N=" << force_total << "\n"
              << "  analytical_eps_z=" << eps_z_analytical
              << "  analytical_eps_r=" << eps_r_analytical << "\n"
              << "  steel-like E=2e11 would shorten the 20 mm bar by about 1e-6 m"
              << " at this force; that displacement is not the quantity computed here.\n"
              << "  dp_m=" << dp << " end_time_s=" << end_time << "\n";

    SPHSystem sph_system(system_bounds_for(dp), dp);
    SolidBody lower_body(sph_system, makeShared<LowerWorkpiece>(SlotMode::None, "LowerWorkpiece"));
    lower_body.defineAdaptationRatios(1.15, 1.0);
    lower_body.defineMatterMaterial<LinearElasticSolid>(rho0_s, Youngs_modulus_axial, poisson);
    lower_body.generateParticles<BaseParticles, Lattice>();
    report_body_metrics("LowerWorkpiece", lower_body, theoretical_lower_volume(false));

    InnerRelation lower_inner(lower_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected(lower_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_1st(lower_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_2nd(lower_inner);

    const Real slab = 0.55 * particle_spacing_ref;
    GeometricShapeBox bottom_shape(
        Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * particle_spacing_ref)),
        Vec3d(cylinder_radius + particle_spacing_ref, cylinder_radius + particle_spacing_ref, slab),
        "BottomAxialSupport");
    GeometricShapeBox top_shape(
        Transform(Vec3d(0.0, 0.0, -0.5 * particle_spacing_ref)),
        Vec3d(cylinder_radius + particle_spacing_ref, cylinder_radius + particle_spacing_ref, slab),
        "TopAxialLoad");
    BodyRegionByParticle bottom_part(lower_body, bottom_shape);
    BodyRegionByParticle top_part(lower_body, top_shape);
    if (bottom_part.body_part_particles_.empty() || top_part.body_part_particles_.empty())
    {
        std::cerr << "elastic_axial: empty support or load layer\n";
        return 1;
    }

    // Keep vx, vy. Zero vz. Vec components of 0 are the constrained axes.
    SimpleDynamics<FixedInAxisDirection> constrain_bottom(bottom_part, Vec3d(1.0, 1.0, 0.0));
    // Remove rigid translation in x and y. Leave z free so the bar can shorten.
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_com(lower_body, Vec3d(1.0, 1.0, 0.0));

    const Vec3d force_per_particle(0.0, 0.0, -force_total / static_cast<Real>(top_part.body_part_particles_.size()));
    SimpleDynamics<solid_dynamics::ForceInBodyRegion> apply_force(top_part, force_per_particle, ramp_time);

    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        damp(0.2, lower_inner, "Velocity", physical_viscosity);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> acoustic_dt(lower_body, 0.2);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected.exec();

    BaseParticles &particles = lower_body.getBaseParticles();
    Vecd *pos = particles.ParticlePositions();
    std::vector<Vec3d> pos0(particles.TotalRealParticles());
    for (size_t i = 0; i < pos0.size(); ++i)
        pos0[i] = pos[i];

    const Real z_bottom0 = mean_coordinate(bottom_part.body_part_particles_, pos, 2);
    const Real z_top0 = mean_coordinate(top_part.body_part_particles_, pos, 2);
    const Real gauge0 = z_top0 - z_bottom0;

    IndexVector interior_bottom_ids;
    IndexVector interior_top_ids;
    for (size_t i = 0; i < pos0.size(); ++i)
    {
        if (std::abs(pos0[i][2] - (z_bottom0 + interior_inset)) < 0.5 * particle_spacing_ref)
            interior_bottom_ids.push_back(i);
        if (std::abs(pos0[i][2] - (z_top0 - interior_inset)) < 0.5 * particle_spacing_ref)
            interior_top_ids.push_back(i);
    }
    const Real z_ib0 = mean_coordinate(interior_bottom_ids, pos, 2);
    const Real z_it0 = mean_coordinate(interior_top_ids, pos, 2);
    const Real interior_gauge0 = z_it0 - z_ib0;

    std::vector<size_t> radial_ids;
    const Real mid_z = 0.5 * (z_bottom0 + z_top0);
    for (size_t i = 0; i < pos0.size(); ++i)
    {
        const Real r0 = std::sqrt(pos0[i][0] * pos0[i][0] + pos0[i][1] * pos0[i][1]);
        if (r0 > 0.75 * cylinder_radius && std::abs(pos0[i][2] - mid_z) < radial_halfband)
            radial_ids.push_back(i);
    }

    auto axial_strain = [&]() -> Real
    {
        const Real z_bottom = mean_coordinate(bottom_part.body_part_particles_, pos, 2);
        const Real z_top = mean_coordinate(top_part.body_part_particles_, pos, 2);
        return gauge0 != 0.0 ? ((z_top - z_bottom) - gauge0) / gauge0 : 0.0;
    };
    auto radial_strain = [&]() -> Real
    {
        if (radial_ids.empty())
            return 0.0;
        Real sum = 0.0;
        for (size_t i : radial_ids)
        {
            const Real r0 = std::sqrt(pos0[i][0] * pos0[i][0] + pos0[i][1] * pos0[i][1]);
            const Real r = std::sqrt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]);
            sum += (r - r0) / r0;
        }
        return sum / static_cast<Real>(radial_ids.size());
    };

    BodyStatesRecordingToVtp write_states(sph_system);
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int aborted = 0;
    Real dt = 0.0;
    write_states.writeToFile(0);
    TickCount t1 = TickCount::now();

    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period && physical_time < end_time)
        {
            apply_force.exec(dt);
            stress_1st.exec(dt);
            constrain_bottom.exec(dt);
            constrain_com.exec(dt);
            damp.exec(dt);
            constrain_bottom.exec(dt);
            constrain_com.exec(dt);
            stress_2nd.exec(dt);

            ite++;
            dt = acoustic_dt.exec();
            integration_time += dt;
            physical_time += dt;
            if (max_speed(particles) > 20.0)
            {
                std::cerr << "elastic_axial abort: max speed exceeded 20 m/s at t="
                          << physical_time << "\n";
                aborted = 1;
                break;
            }
            lower_body.updateCellLinkedList();
            lower_inner.updateConfiguration();
        }
        corrected.exec();
        std::cout << "N=" << ite << " Time: " << physical_time
                  << " dt: " << dt
                  << " eps_z: " << axial_strain()
                  << " eps_r: " << radial_strain()
                  << " ke_J: " << kinetic_energy(particles) << "\n";
        write_states.writeToFile();
        if (aborted)
            break;
    }

    TickCount t2 = TickCount::now();
    const Real eps_z = axial_strain();
    const Real eps_r = radial_strain();
    const Real z_bottom = mean_coordinate(bottom_part.body_part_particles_, pos, 2);
    const Real z_top = mean_coordinate(top_part.body_part_particles_, pos, 2);
    const Real z_ib = mean_coordinate(interior_bottom_ids, pos, 2);
    const Real z_it = mean_coordinate(interior_top_ids, pos, 2);
    const Real eps_z_interior = interior_gauge0 != 0.0 ? ((z_it - z_ib) - interior_gauge0) / interior_gauge0 : 0.0;

    Matd *F = particles.getVariableDataByName<Matd>("DeformationGradient");
    const Real lambda0 = Youngs_modulus_axial * poisson /
                         ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const Real shear = Youngs_modulus_axial / (2.0 * (1.0 + poisson));
    Real pk2_zz_sum = 0.0;
    size_t core_n = 0;
    for (size_t i = 0; i < pos0.size(); ++i)
    {
        const Real r0 = std::sqrt(pos0[i][0] * pos0[i][0] + pos0[i][1] * pos0[i][1]);
        if (r0 < 0.5 * cylinder_radius && pos0[i][2] > z_bottom0 + core_inset &&
            pos0[i][2] < z_top0 - core_inset)
        {
            const Matd strain = 0.5 * (F[i].transpose() + F[i]) - Matd::Identity();
            const Matd pk2 = lambda0 * strain.trace() * Matd::Identity() + 2.0 * shear * strain;
            pk2_zz_sum += pk2(2, 2);
            ++core_n;
        }
    }
    const Real pk2_zz = core_n > 0 ? pk2_zz_sum / static_cast<Real>(core_n) : 0.0;

    Vecd *force = particles.getVariableDataByName<Vecd>("Force");
    Vecd *force_prior = particles.getVariableDataByName<Vecd>("ForcePrior");
    Real applied_z = 0.0;
    for (size_t i : top_part.body_part_particles_)
        applied_z += force_prior[i][2];
    Real bottom_internal_z = 0.0;
    for (size_t i : bottom_part.body_part_particles_)
        bottom_internal_z += force[i][2] + force_prior[i][2];
    const Real reaction_z = -bottom_internal_z;

    const Real err_z = std::abs(eps_z - eps_z_analytical) / std::abs(eps_z_analytical);
    const Real err_z_interior = std::abs(eps_z_interior - eps_z_analytical) / std::abs(eps_z_analytical);
    const Real err_r = std::abs(eps_r - eps_r_analytical) / std::abs(eps_r_analytical);
    const Real err_stress = std::abs(pk2_zz - (-nominal_pressure)) / nominal_pressure;
    const Real err_force = std::abs(std::abs(reaction_z) - force_total) / force_total;
    char metrics_name[160];
    std::snprintf(metrics_name, sizeof(metrics_name), "spdb_elastic_axial_dp%.4g_metrics.txt", dp);

    std::cout << "elastic_axial wall time: " << (t2 - t1).seconds() << " s\n"
              << "  gauge0_m=" << gauge0
              << "  bottom_dz_m=" << (z_bottom - z_bottom0)
              << "  top_dz_m=" << (z_top - z_top0) << "\n"
              << "  eps_z=" << eps_z << " analytical=" << eps_z_analytical << " rel_err=" << err_z << "\n"
              << "  eps_z_interior=" << eps_z_interior << " rel_err=" << err_z_interior
              << " interior_n=" << interior_bottom_ids.size() << "/" << interior_top_ids.size() << "\n"
              << "  eps_r=" << eps_r << " analytical=" << eps_r_analytical << " rel_err=" << err_r << "\n"
              << "  core_pk2_zz_Pa=" << pk2_zz << " rel_err=" << err_stress << " core_n=" << core_n << "\n"
              << "  applied_force_z_N=" << applied_z
              << "  reaction_z_N=" << reaction_z << " rel_err=" << err_force << "\n"
              << "  metrics_file=" << metrics_name << "\n"
              << "  NOTE: synthetic modulus chosen so the shortening is observable at this spacing.\n";

    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "mode=elastic_axial\n";
    out << "units=SI\n";
    out << "material_source=synthetic\n";
    out << "material_note=NOT_950C_20Cr13_prediction\n";
    out << "constitutive=LinearElasticSolid_Integration1stHalfPK2\n";
    out << "youngs_modulus_Pa=" << Youngs_modulus_axial << "\n";
    out << "poisson=" << poisson << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "interior_inset_m=" << interior_inset << "\n";
    out << "core_inset_m=" << core_inset << "\n";
    out << "nominal_pressure_Pa=" << nominal_pressure << "\n";
    out << "nominal_force_N=" << force_total << "\n";
    out << "load_mode=constant_nominal_force\n";
    out << "ramp_time_s=" << ramp_time << "\n";
    out << "end_time_s=" << end_time << "\n";
    out << "physical_time_s=" << physical_time << "\n";
    out << "aborted=" << aborted << "\n";
    out << "iterations=" << ite << "\n";
    out << "bottom_particles=" << bottom_part.body_part_particles_.size() << "\n";
    out << "top_particles=" << top_part.body_part_particles_.size() << "\n";
    out << "radial_sample_particles=" << radial_ids.size() << "\n";
    out << "gauge0_m=" << gauge0 << "\n";
    out << "bottom_dz_m=" << (z_bottom - z_bottom0) << "\n";
    out << "top_dz_m=" << (z_top - z_top0) << "\n";
    out << "eps_z=" << eps_z << "\n";
    out << "eps_z_analytical=" << eps_z_analytical << "\n";
    out << "eps_z_rel_err=" << err_z << "\n";
    out << "eps_z_interior=" << eps_z_interior << "\n";
    out << "eps_z_interior_rel_err=" << err_z_interior << "\n";
    out << "interior_gauge0_m=" << interior_gauge0 << "\n";
    out << "interior_bottom_particles=" << interior_bottom_ids.size() << "\n";
    out << "interior_top_particles=" << interior_top_ids.size() << "\n";
    out << "eps_r=" << eps_r << "\n";
    out << "eps_r_analytical=" << eps_r_analytical << "\n";
    out << "eps_r_rel_err=" << err_r << "\n";
    out << "core_pk2_zz_Pa=" << pk2_zz << "\n";
    out << "core_pk2_zz_rel_err=" << err_stress << "\n";
    out << "core_particles=" << core_n << "\n";
    out << "applied_force_z_N=" << applied_z << "\n";
    out << "reaction_z_N=" << reaction_z << "\n";
    out << "reaction_rel_err=" << err_force << "\n";
    out << "kinetic_energy_J=" << kinetic_energy(particles) << "\n";
    out << "wall_time_s=" << (t2 - t1).seconds() << "\n";
    return aborted;
}

struct InterfaceGap
{
    size_t pairs = 0;
    size_t penetrating = 0;
    Real min_center_gap = 0.0;
    Real mean_center_gap = 0.0;
    Real min_surface_gap = 0.0;
    Real max_penetration = 0.0;
};

// Surface gap = center gap - dp. The lattice faces start dp/2 inside each body,
// so a flush geometric join has center gap = dp and surface gap = 0.
// A negative surface gap means those faces have crossed.
InterfaceGap measure_interface_gap(const IndexVector &lower_ids, Vecd *pos_l,
                                   const IndexVector &upper_ids, Vecd *pos_u, Real dp)
{
    InterfaceGap gap;
    if (lower_ids.empty() || upper_ids.empty())
        return gap;
    Real min_c = std::numeric_limits<Real>::max();
    Real min_s = std::numeric_limits<Real>::max();
    Real sum_c = 0.0;
    for (size_t i : lower_ids)
    {
        Real best_r = std::numeric_limits<Real>::max();
        size_t j_best = upper_ids.front();
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
        const Real center = pos_u[j_best][2] - pos_l[i][2];
        const Real surface = center - dp;
        min_c = std::min(min_c, center);
        min_s = std::min(min_s, surface);
        sum_c += center;
        gap.max_penetration = std::max(gap.max_penetration, std::max(Real(0.0), -surface));
        if (surface < 0.0)
            ++gap.penetrating;
        ++gap.pairs;
    }
    if (gap.pairs > 0)
    {
        gap.min_center_gap = min_c;
        gap.mean_center_gap = sum_c / static_cast<Real>(gap.pairs);
        gap.min_surface_gap = min_s;
    }
    return gap;
}

Real sum_field_z(BaseParticles &particles, const char *name)
{
    Vecd *field = particles.getVariableDataByName<Vecd>(name);
    Real sum = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        sum += field[i][2];
    return sum;
}

IndexVector ids_near_z(Vecd *pos, size_t n, Real z_target, Real tol)
{
    IndexVector ids;
    for (size_t i = 0; i < n; ++i)
        if (std::abs(pos[i][2] - z_target) <= tol)
            ids.push_back(i);
    return ids;
}

//----------------------------------------------------------------------
// Two flush cylinders, no slot. Same nominal force as the single-bar
// check. Lower bottom fixed in z, upper top force-controlled.
// Reports interface penetration and force balance. Synthetic E.
//----------------------------------------------------------------------
int run_elastic_contact(Real dp, Real end_time, Real stiffness_scale, Real youngs_modulus, SlotMode slot_mode)
{
    particle_spacing_ref = dp;
    resolution_ref = dp;
    const bool with_slot = (slot_mode == SlotMode::LowerOnly || slot_mode == SlotMode::Both);
    const Real force_total = nominal_axial_force();
    const Real eps_z_analytical = -nominal_pressure / youngs_modulus;
    const Real ramp_time = std::min(1.0e-3, 0.25 * end_time);
    const bool demo_modulus = std::abs(youngs_modulus - 100.0e9) / 100.0e9 < 1.0e-6;
    const Real output_period = std::min(1.0e-3, std::max(0.25 * end_time, 1.0e-6));
    const Real slab = 0.55 * dp;

    std::cout << "SPDB elastic_contact\n"
              << "  NOTE: LinearElasticSolid, E=" << youngs_modulus
              << " Pa, nu=" << poisson << ". "
              << (demo_modulus ? "demo_assumption, NOT a 950C 20Cr13 prediction.\n"
                               : "synthetic, NOT 950C 20Cr13.\n")
              << "  slot_mode=" << static_cast<int>(slot_mode)
              << "  contact=ContactForce (not RepulsionForceCK)\n"
              << "  nominal_pressure_Pa=" << nominal_pressure
              << "  nominal_force_N=" << force_total << "\n"
              << "  dp_m=" << dp << " end_time_s=" << end_time << "\n"
              << "  contact_stiffness_scale=" << stiffness_scale
              << " (multiplies rho*c^2; c_contact = c0*sqrt(scale))\n"
              << "  surface_gap = (z_upper - z_lower) - dp on join-face pairs outside the slot.\n"
              << "  surface_gap < 0 means the reconstructed faces have crossed.\n"
              << "  Cavity metrics use particle centers. Volume is not claimed.\n";

    SPHSystem sph_system(system_bounds_for(dp), dp);
    SolidBody lower_body(sph_system, makeShared<LowerWorkpiece>(slot_mode, "LowerWorkpiece"));
    lower_body.defineAdaptationRatios(1.15, 1.0);
    ContactScaledLinearElastic &lower_mat = lower_body.defineMatterMaterial<ContactScaledLinearElastic>(
        rho0_s, youngs_modulus, poisson, stiffness_scale);
    lower_body.generateParticles<BaseParticles, Lattice>();
    SolidBody upper_body(sph_system, makeShared<UpperWorkpiece>(slot_mode, "UpperWorkpiece"));
    upper_body.defineAdaptationRatios(1.15, 1.0);
    ContactScaledLinearElastic &upper_mat = upper_body.defineMatterMaterial<ContactScaledLinearElastic>(
        rho0_s, youngs_modulus, poisson, stiffness_scale);
    upper_body.generateParticles<BaseParticles, Lattice>();
    report_body_metrics("LowerWorkpiece", lower_body, theoretical_lower_volume(with_slot));
    report_body_metrics("UpperWorkpiece", upper_body, theoretical_upper_volume());

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

    GeometricShapeBox bottom_shape(
        Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab),
        "BottomAxialSupport");
    GeometricShapeBox top_shape(
        Transform(Vec3d(0.0, 0.0, cylinder_height - 0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab),
        "TopAxialLoad");
    BodyRegionByParticle bottom_part(lower_body, bottom_shape);
    BodyRegionByParticle top_part(upper_body, top_shape);
    if (bottom_part.body_part_particles_.empty() || top_part.body_part_particles_.empty())
    {
        std::cerr << "elastic_contact: empty support or load layer\n";
        return 1;
    }

    SimpleDynamics<FixedInAxisDirection> constrain_bottom(bottom_part, Vec3d(1.0, 1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_com_l(lower_body, Vec3d(1.0, 1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_com_u(upper_body, Vec3d(1.0, 1.0, 0.0));
    const Vec3d force_per_particle(0.0, 0.0, -force_total / static_cast<Real>(top_part.body_part_particles_.size()));
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

    BaseParticles &lower_particles = lower_body.getBaseParticles();
    BaseParticles &upper_particles = upper_body.getBaseParticles();
    Vecd *pos_l = lower_particles.ParticlePositions();
    Vecd *pos_u = upper_particles.ParticlePositions();
    IndexVector lower_face =
        with_slot ? spdb::ids_near_z_outside_slot(pos_l, lower_particles.TotalRealParticles(), -0.5 * dp, 0.5 * dp)
                  : ids_near_z(pos_l, lower_particles.TotalRealParticles(), -0.5 * dp, 0.5 * dp);
    IndexVector upper_face =
        with_slot ? spdb::ids_near_z_outside_slot(pos_u, upper_particles.TotalRealParticles(), 0.5 * dp, 0.5 * dp)
                  : ids_near_z(pos_u, upper_particles.TotalRealParticles(), 0.5 * dp, 0.5 * dp);
    if (lower_face.empty() || upper_face.empty())
    {
        std::cerr << "elastic_contact: empty interface layer\n";
        return 1;
    }

    const Real z_lb0 = mean_coordinate(bottom_part.body_part_particles_, pos_l, 2);
    const Real z_lf0 = mean_coordinate(lower_face, pos_l, 2);
    const Real z_uf0 = mean_coordinate(upper_face, pos_u, 2);
    const Real z_ut0 = mean_coordinate(top_part.body_part_particles_, pos_u, 2);
    const Real gauge_l0 = z_lf0 - z_lb0;
    const Real gauge_u0 = z_ut0 - z_uf0;
    const InterfaceGap gap0 = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const spdb::SlotCavityMetrics cavity0 =
        with_slot ? spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                             upper_particles.TotalRealParticles(), dp)
                  : spdb::SlotCavityMetrics{};

    auto strain_of = [](Real z_a, Real z_b, Real gauge0) -> Real
    {
        return gauge0 != 0.0 ? ((z_b - z_a) - gauge0) / gauge0 : 0.0;
    };

    BodyStatesRecordingToVtp write_states(sph_system);
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int aborted = 0;
    Real dt = 0.0;
    write_states.writeToFile(0);
    std::cout << "initial_pairs=" << gap0.pairs
              << " initial_min_surface_gap_m=" << gap0.min_surface_gap
              << " initial_mean_center_gap_m=" << gap0.mean_center_gap;
    if (with_slot)
        std::cout << " initial_cavity_net_height_m=" << cavity0.center_net_height
                  << " floor_n=" << cavity0.floor_n << " roof_n=" << cavity0.roof_n;
    std::cout << "\n";
    TickCount t1 = TickCount::now();

    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period && physical_time < end_time)
        {
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

            ite++;
            dt = SMIN(dt_l.exec(), dt_u.exec());
            integration_time += dt;
            physical_time += dt;
            if (max_speed(lower_particles) > 20.0 || max_speed(upper_particles) > 20.0)
            {
                std::cerr << "elastic_contact abort: max speed exceeded 20 m/s at t="
                          << physical_time << "\n";
                aborted = 1;
                break;
            }
            lower_body.updateCellLinkedList();
            upper_body.updateCellLinkedList();
            lower_inner.updateConfiguration();
            upper_inner.updateConfiguration();
            lower_upper.updateConfiguration();
            upper_lower.updateConfiguration();
        }
        corrected_l.exec();
        corrected_u.exec();
        const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
        const Real eps_l = strain_of(mean_coordinate(bottom_part.body_part_particles_, pos_l, 2),
                                     mean_coordinate(lower_face, pos_l, 2), gauge_l0);
        const Real eps_u = strain_of(mean_coordinate(upper_face, pos_u, 2),
                                     mean_coordinate(top_part.body_part_particles_, pos_u, 2), gauge_u0);
        std::cout << "N=" << ite << " Time: " << physical_time
                  << " dt: " << dt
                  << " eps_l: " << eps_l
                  << " eps_u: " << eps_u
                  << " min_surface_gap_m: " << gap.min_surface_gap
                  << " max_penetration_m: " << gap.max_penetration
                  << " ke_J: " << (kinetic_energy(lower_particles) + kinetic_energy(upper_particles));
        if (with_slot)
        {
            const spdb::SlotCavityMetrics cav =
                spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                          upper_particles.TotalRealParticles(), dp);
            std::cout << " cavity_net_height_m: " << cav.center_net_height
                      << " min_clearance_m: " << cav.min_clearance
                      << " top_wall_defl_m: " << cav.top_wall_deflection;
        }
        std::cout << "\n";
        write_states.writeToFile();
        if (aborted)
            break;
    }

    TickCount t2 = TickCount::now();
    const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const Real z_lb = mean_coordinate(bottom_part.body_part_particles_, pos_l, 2);
    const Real z_lf = mean_coordinate(lower_face, pos_l, 2);
    const Real z_uf = mean_coordinate(upper_face, pos_u, 2);
    const Real z_ut = mean_coordinate(top_part.body_part_particles_, pos_u, 2);
    const Real eps_l = strain_of(z_lb, z_lf, gauge_l0);
    const Real eps_u = strain_of(z_uf, z_ut, gauge_u0);
    const Real applied_z = sum_field_z(upper_particles, "ForceInBodyRegion");
    const Real contact_l_z = sum_field_z(lower_particles, "RepulsionForce");
    const Real contact_u_z = sum_field_z(upper_particles, "RepulsionForce");
    Vecd *force_l = lower_particles.getVariableDataByName<Vecd>("Force");
    Vecd *prior_l = lower_particles.getVariableDataByName<Vecd>("ForcePrior");
    Real bottom_internal_z = 0.0;
    for (size_t i : bottom_part.body_part_particles_)
        bottom_internal_z += force_l[i][2] + prior_l[i][2];
    const Real reaction_z = -bottom_internal_z;
    const Real err_l = std::abs(eps_l - eps_z_analytical) / std::abs(eps_z_analytical);
    const Real err_u = std::abs(eps_u - eps_z_analytical) / std::abs(eps_z_analytical);
    const Real contact_imbalance = std::abs(contact_l_z + contact_u_z) / force_total;
    const Real reaction_err = std::abs(std::abs(reaction_z) - force_total) / force_total;
    const Real contact_vs_load = std::abs(std::abs(contact_u_z) - force_total) / force_total;

    const spdb::SlotCavityMetrics cavity =
        with_slot ? spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                             upper_particles.TotalRealParticles(), dp)
                  : spdb::SlotCavityMetrics{};
    char metrics_name[192];
    if (with_slot)
        std::snprintf(metrics_name, sizeof(metrics_name), "spdb_slot_contact_dp%.4g_E%.4g_k%.4g_metrics.txt", dp,
                      youngs_modulus, stiffness_scale);
    else
        std::snprintf(metrics_name, sizeof(metrics_name), "spdb_elastic_contact_dp%.4g_E%.4g_k%.4g_metrics.txt", dp,
                      youngs_modulus, stiffness_scale);
    std::cout << "elastic_contact wall time: " << (t2 - t1).seconds() << " s\n"
              << "  eps_l=" << eps_l << " rel_err=" << err_l
              << "  eps_u=" << eps_u << " rel_err=" << err_u << "\n"
              << "  min_surface_gap_m=" << gap.min_surface_gap
              << "  max_penetration_m=" << gap.max_penetration
              << "  penetration_over_dp=" << gap.max_penetration / dp
              << "  penetrating_pairs=" << gap.penetrating << "/" << gap.pairs << "\n"
              << "  applied_z_N=" << applied_z
              << "  contact_lower_z_N=" << contact_l_z
              << "  contact_upper_z_N=" << contact_u_z
              << "  reaction_z_N=" << reaction_z << "\n"
              << "  contact_pair_imbalance=" << contact_imbalance
              << "  reaction_rel_err=" << reaction_err
              << "  contact_vs_load=" << contact_vs_load << "\n";
    if (with_slot)
        std::cout << "  cavity_net_height_m=" << cavity.center_net_height
                  << "  delta_net_height_m=" << (cavity.center_net_height - cavity0.center_net_height)
                  << "  min_clearance_m=" << cavity.min_clearance
                  << "  width_x_m=" << cavity.width_x
                  << "  width_y_m=" << cavity.width_y
                  << "  top_wall_defl_m=" << cavity.top_wall_deflection << "\n"
                  << "  NOTE: cavity change must not be read as creep; this run is elastic only.\n";
    std::cout << "  metrics_file=" << metrics_name << "\n";

    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "mode=" << (with_slot ? "slot_contact" : "elastic_contact") << "\n";
    out << "units=SI\n";
    out << "material_source=" << (demo_modulus ? "demo_assumption" : "synthetic") << "\n";
    out << "material_note=NOT_950C_20Cr13_prediction\n";
    out << "constitutive=LinearElasticSolid_Integration1stHalfPK2\n";
    out << "contact=ContactForce\n";
    out << "slot_mode=" << static_cast<int>(slot_mode) << "\n";
    out << "contact_stiffness_scale=" << stiffness_scale << "\n";
    out << "contact_stiffness_Pa=" << lower_mat.ContactStiffness() << "\n";
    out << "upper_contact_stiffness_Pa=" << upper_mat.ContactStiffness() << "\n";
    out << "youngs_modulus_Pa=" << youngs_modulus << "\n";
    out << "poisson=" << poisson << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "nominal_pressure_Pa=" << nominal_pressure << "\n";
    out << "nominal_force_N=" << force_total << "\n";
    out << "load_mode=constant_nominal_force\n";
    out << "friction=none\n";
    out << "ramp_time_s=" << ramp_time << "\n";
    out << "end_time_s=" << end_time << "\n";
    out << "physical_time_s=" << physical_time << "\n";
    out << "aborted=" << aborted << "\n";
    out << "iterations=" << ite << "\n";
    out << "lower_particles=" << lower_particles.TotalRealParticles() << "\n";
    out << "upper_particles=" << upper_particles.TotalRealParticles() << "\n";
    out << "lower_face_particles=" << lower_face.size() << "\n";
    out << "upper_face_particles=" << upper_face.size() << "\n";
    out << "gauge_lower0_m=" << gauge_l0 << "\n";
    out << "gauge_upper0_m=" << gauge_u0 << "\n";
    out << "lower_bottom_dz_m=" << (z_lb - z_lb0) << "\n";
    out << "lower_face_dz_m=" << (z_lf - z_lf0) << "\n";
    out << "upper_face_dz_m=" << (z_uf - z_uf0) << "\n";
    out << "upper_top_dz_m=" << (z_ut - z_ut0) << "\n";
    out << "eps_z_lower=" << eps_l << "\n";
    out << "eps_z_upper=" << eps_u << "\n";
    out << "eps_z_analytical=" << eps_z_analytical << "\n";
    out << "eps_z_lower_rel_err=" << err_l << "\n";
    out << "eps_z_upper_rel_err=" << err_u << "\n";
    out << "initial_pairs=" << gap0.pairs << "\n";
    out << "initial_mean_center_gap_m=" << gap0.mean_center_gap << "\n";
    out << "initial_min_surface_gap_m=" << gap0.min_surface_gap << "\n";
    out << "pairs=" << gap.pairs << "\n";
    out << "penetrating_pairs=" << gap.penetrating << "\n";
    out << "min_center_gap_m=" << gap.min_center_gap << "\n";
    out << "mean_center_gap_m=" << gap.mean_center_gap << "\n";
    out << "min_surface_gap_m=" << gap.min_surface_gap << "\n";
    out << "max_penetration_m=" << gap.max_penetration << "\n";
    out << "max_penetration_over_dp=" << (dp > 0.0 ? gap.max_penetration / dp : 0.0) << "\n";
    out << "cavity_volume_claimed=0\n";
    out << "initial_cavity_net_height_m=" << cavity0.center_net_height << "\n";
    out << "cavity_net_height_m=" << cavity.center_net_height << "\n";
    out << "delta_cavity_net_height_m=" << (cavity.center_net_height - cavity0.center_net_height) << "\n";
    out << "min_clearance_m=" << cavity.min_clearance << "\n";
    out << "cavity_width_x_m=" << cavity.width_x << "\n";
    out << "cavity_width_y_m=" << cavity.width_y << "\n";
    out << "top_wall_deflection_m=" << cavity.top_wall_deflection << "\n";
    out << "cavity_floor_particles=" << cavity.floor_n << "\n";
    out << "cavity_roof_particles=" << cavity.roof_n << "\n";
    out << "applied_force_z_N=" << applied_z << "\n";
    out << "contact_lower_z_N=" << contact_l_z << "\n";
    out << "contact_upper_z_N=" << contact_u_z << "\n";
    out << "contact_pair_imbalance=" << contact_imbalance << "\n";
    out << "reaction_z_N=" << reaction_z << "\n";
    out << "reaction_rel_err=" << reaction_err << "\n";
    out << "contact_vs_load=" << contact_vs_load << "\n";
    out << "kinetic_energy_J=" << (kinetic_energy(lower_particles) + kinetic_energy(upper_particles)) << "\n";
    out << "wall_time_s=" << (t2 - t1).seconds() << "\n";
    return aborted;
}

void write_metrics_file(const std::string &path, SolidBody &lower, SolidBody &upper,
                        SlotMode slot_mode, bool with_slot_theory)
{
    std::ofstream out(path);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "units=SI\n";
    out << "material_source=synthetic\n";
    out << "material_note=NOT_950C_20Cr13_prediction\n";
    out << "particle_spacing_m=" << particle_spacing_ref << "\n";
    out << "cylinder_radius_m=" << cylinder_radius << "\n";
    out << "cylinder_height_m=" << cylinder_height << "\n";
    out << "slot_mode=" << static_cast<int>(slot_mode) << "\n";
    out << "theoretical_slot_volume_m3=" << theoretical_slot_volume() << "\n";
    out << "theoretical_lower_volume_m3=" << theoretical_lower_volume(with_slot_theory) << "\n";
    out << "theoretical_upper_volume_m3=" << theoretical_upper_volume() << "\n";
    out << "lower_particles=" << lower.getBaseParticles().TotalRealParticles() << "\n";
    out << "upper_particles=" << upper.getBaseParticles().TotalRealParticles() << "\n";
    out << "lower_discrete_volume_m3=" << particle_volume_sum(lower.getBaseParticles()) << "\n";
    out << "upper_discrete_volume_m3=" << particle_volume_sum(upper.getBaseParticles()) << "\n";
}

struct CkHoldPoint
{
    Real eps_z = 0.0;
    Real eps_z_interior = 0.0;
    Real eps_r = 0.0;
    Real ke = 0.0;
    Real speed = 0.0;
    Real core_cauchy_zz = 0.0;
    Real core_sigma_eq = 0.0;
    Real core_inelastic = 0.0;
    Real max_inelastic = 0.0;
    size_t core_n = 0;
    Real reaction_z = 0.0;
    Real applied_z = 0.0;
};

struct MechanicalSnapshot
{
    std::vector<Vecd> pos;
    std::vector<Vecd> vel;
    std::vector<Matd> deformation;
    std::vector<Matd> deformation_rate;
    std::vector<Matd> inverse_plastic;
    std::vector<Real> density;
    std::vector<Real> inelastic;
    std::vector<Real> range_flag;
};

MechanicalSnapshot capture_mechanics(BaseParticles &particles)
{
    const size_t n = particles.TotalRealParticles();
    MechanicalSnapshot snap;
    Vecd *pos = particles.ParticlePositions();
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Matd *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    Matd *rate = particles.getVariableDataByName<Matd>("DeformationRate");
    Matd *inverse_plastic = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *density = particles.getVariableDataByName<Real>("Density");
    Real *inelastic = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *range_flag = particles.getVariableDataByName<Real>("ArrheniusRangeFlag");
    snap.pos.assign(pos, pos + n);
    snap.vel.assign(vel, vel + n);
    snap.deformation.assign(deformation, deformation + n);
    snap.deformation_rate.assign(rate, rate + n);
    snap.inverse_plastic.assign(inverse_plastic, inverse_plastic + n);
    snap.density.assign(density, density + n);
    snap.inelastic.assign(inelastic, inelastic + n);
    snap.range_flag.assign(range_flag, range_flag + n);
    return snap;
}

void restore_mechanics(BaseParticles &particles, const MechanicalSnapshot &snap)
{
    const size_t n = particles.TotalRealParticles();
    Vecd *pos = particles.ParticlePositions();
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Matd *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    Matd *rate = particles.getVariableDataByName<Matd>("DeformationRate");
    Matd *inverse_plastic = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *density = particles.getVariableDataByName<Real>("Density");
    Real *inelastic = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *range_flag = particles.getVariableDataByName<Real>("ArrheniusRangeFlag");
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = snap.pos[i];
        vel[i] = snap.vel[i];
        deformation[i] = snap.deformation[i];
        rate[i] = snap.deformation_rate[i];
        inverse_plastic[i] = snap.inverse_plastic[i];
        density[i] = snap.density[i];
        inelastic[i] = snap.inelastic[i];
        range_flag[i] = snap.range_flag[i];
    }
}

void remove_rigid_xy(BaseParticles &particles, const Vecd &com0)
{
    Vecd *pos = particles.ParticlePositions();
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    Real mass_sum = 0.0;
    Real mx = 0.0;
    Real my = 0.0;
    Real mvx = 0.0;
    Real mvy = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        mass_sum += mass[i];
        mx += mass[i] * pos[i][0];
        my += mass[i] * pos[i][1];
        mvx += mass[i] * vel[i][0];
        mvy += mass[i] * vel[i][1];
    }
    if (!(mass_sum > 0.0))
        return;
    const Real dx = mx / mass_sum - com0[0];
    const Real dy = my / mass_sum - com0[1];
    const Real dvx = mvx / mass_sum;
    const Real dvy = mvy / mass_sum;
    for (size_t i = 0; i < n; ++i)
    {
        pos[i][0] -= dx;
        pos[i][1] -= dy;
        vel[i][0] -= dvx;
        vel[i][1] -= dvy;
    }
}

//----------------------------------------------------------------------
// No-slot cylinder on StructureIntegration1stHalf.
// Pseudo-time relaxes the explicit solid with inelastic updates frozen.
// Each accepted hold step calls commit_ren2014_physical_increment once.
//----------------------------------------------------------------------
int run_ck_hold(Real dp, Real pseudo_cap, Real physical_dt, Real damping_factor, int hold_steps)
{
    const Real demo_youngs = 100.0e9;
    const Real demo_temperature = 1223.15;
    const Real force_total = nominal_axial_force();
    const Real eps_z_hooke = -nominal_pressure / demo_youngs;
    const Real ramp_pseudo = 0.25 * pseudo_cap;
    const Real interior_inset = std::max(4.0 * dp, 0.004);
    const Real core_inset = std::max(3.0 * dp, 0.003);
    const Real radial_halfband = std::max(2.0 * dp, 0.002);
    particle_spacing_ref = dp;

    std::cout << "SPDB ck_hold\n"
              << "  model=ren2014_x20cr13  E_Pa=" << demo_youngs
              << " source=demo_assumption  T_K=" << demo_temperature << "\n"
              << "  NOT a 950 C 20Cr13 process prediction and NOT a 7200 s result.\n"
              << "  inelastic_during_relaxation=false\n"
              << "  nominal_force_N=" << force_total
              << "  hooke_eps_z=" << eps_z_hooke << "\n"
              << "  dp_m=" << dp << " pseudo_cap_s=" << pseudo_cap
              << " physical_dt_s=" << physical_dt
              << " hold_steps=" << hold_steps
              << " damping_factor=" << damping_factor << "\n"
              << "  Loading to 10 MPa is numerical relaxation. Hold time starts at the elastic equilibrium.\n";

    SPHSystem sph_system(system_bounds_for(dp), dp);
    char output_folder[192];
    std::snprintf(output_folder, sizeof(output_folder), "./output_ck_hold_dp%.4g_damp%.4g", dp, damping_factor);
    IO::getEnvironment().resetOutputFolder(output_folder);
    auto &shape = sph_system.addShape<LowerWorkpiece>(SlotMode::None, "CkHoldCylinder");
    auto &body = sph_system.addBody<SolidBody>(shape);
    body.defineAdaptationRatios(1.15, 1.0);
    auto &material = body.defineMatterMaterial<Ren2014X20Cr13Solid>(rho0_s, demo_youngs, poisson, demo_temperature, false);
    body.generateParticles<BaseParticles, Lattice>();
    report_body_metrics("CkHoldCylinder", body, theoretical_lower_volume(false));

    auto &inner = sph_system.addInnerRelation(body, ConfigType::Lagrangian);
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();

    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(body));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(inner));
    auto &correction = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(inner);
    auto &step_1st = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, Ren2014X20Cr13Solid, LinearCorrectionCK>(inner, damping_factor);
    auto &step_2nd = main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(inner);
    auto &acoustic = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(body, 0.3);

    const Real slab = 0.55 * dp;
    auto &bottom_shape = sph_system.addShape<GeometricShapeBox>(
        Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "BottomAxialSupport");
    auto &top_shape = sph_system.addShape<GeometricShapeBox>(
        Transform(Vec3d(0.0, 0.0, -0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "TopAxialLoad");
    BodyRegionByParticle bottom_part(body, bottom_shape);
    BodyRegionByParticle top_part(body, top_shape);
    if (bottom_part.body_part_particles_.empty() || top_part.body_part_particles_.empty())
    {
        std::cerr << "ck_hold: empty support or load layer\n";
        return 1;
    }

    SingleVariable<Real> *force_scale = sph_system.registerSystemVariable<Real>("AxialForceScale", 0.0);
    const Real force_per_particle = -force_total / static_cast<Real>(top_part.body_part_particles_.size());
    auto &support = main_methods.addStateDynamics<AxialZSupportCK<BodyRegionByParticle>>(bottom_part);
    auto &apply_force = main_methods.addStateDynamics<NominalAxialForceCK<BodyRegionByParticle>>(
        top_part, force_per_particle, *force_scale);
    auto &recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);

    lagrangian_configuration.exec();
    correction.exec();

    BaseParticles &particles = body.getBaseParticles();
    Vecd *pos = particles.ParticlePositions();
    const size_t n_particles = particles.TotalRealParticles();
    std::vector<Vec3d> pos0(n_particles);
    Real mass_sum = 0.0;
    Vecd com0 = Vecd::Zero();
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    for (size_t i = 0; i < n_particles; ++i)
    {
        pos0[i] = pos[i];
        mass_sum += mass[i];
        com0 += mass[i] * pos[i];
    }
    if (mass_sum > 0.0)
        com0 /= mass_sum;

    const Real z_bottom0 = mean_coordinate(bottom_part.body_part_particles_, pos, 2);
    const Real z_top0 = mean_coordinate(top_part.body_part_particles_, pos, 2);
    const Real gauge0 = z_top0 - z_bottom0;
    IndexVector interior_bottom_ids;
    IndexVector interior_top_ids;
    std::vector<size_t> radial_ids;
    std::vector<size_t> core_ids;
    const Real mid_z = 0.5 * (z_bottom0 + z_top0);
    for (size_t i = 0; i < n_particles; ++i)
    {
        if (std::abs(pos0[i][2] - (z_bottom0 + interior_inset)) < 0.5 * dp)
            interior_bottom_ids.push_back(i);
        if (std::abs(pos0[i][2] - (z_top0 - interior_inset)) < 0.5 * dp)
            interior_top_ids.push_back(i);
        const Real r0 = std::sqrt(pos0[i][0] * pos0[i][0] + pos0[i][1] * pos0[i][1]);
        if (r0 > 0.75 * cylinder_radius && std::abs(pos0[i][2] - mid_z) < radial_halfband)
            radial_ids.push_back(i);
        if (r0 < 0.5 * cylinder_radius && pos0[i][2] > z_bottom0 + core_inset && pos0[i][2] < z_top0 - core_inset)
            core_ids.push_back(i);
    }
    const Real z_ib0 = mean_coordinate(interior_bottom_ids, pos, 2);
    const Real z_it0 = mean_coordinate(interior_top_ids, pos, 2);
    const Real interior_gauge0 = z_it0 - z_ib0;
    const Real shear = material.ShearModulus();
    const Real bulk = material.BulkModulus();

    auto sample = [&]() -> CkHoldPoint
    {
        CkHoldPoint point;
        const Real z_bottom = mean_coordinate(bottom_part.body_part_particles_, pos, 2);
        const Real z_top = mean_coordinate(top_part.body_part_particles_, pos, 2);
        point.eps_z = gauge0 != 0.0 ? ((z_top - z_bottom) - gauge0) / gauge0 : 0.0;
        const Real z_ib = mean_coordinate(interior_bottom_ids, pos, 2);
        const Real z_it = mean_coordinate(interior_top_ids, pos, 2);
        point.eps_z_interior = interior_gauge0 != 0.0 ? ((z_it - z_ib) - interior_gauge0) / interior_gauge0 : 0.0;
        if (!radial_ids.empty())
        {
            Real sum = 0.0;
            for (size_t i : radial_ids)
            {
                const Real r0 = std::sqrt(pos0[i][0] * pos0[i][0] + pos0[i][1] * pos0[i][1]);
                const Real r = std::sqrt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]);
                sum += (r - r0) / r0;
            }
            point.eps_r = sum / static_cast<Real>(radial_ids.size());
        }
        Matd *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
        Matd *inverse_plastic = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
        Real *inelastic = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
        Real zz_sum = 0.0;
        Real eq_sum = 0.0;
        Real p_sum = 0.0;
        point.core_n = core_ids.size();
        for (size_t i : core_ids)
        {
            Real sigma_eq = 0.0;
            const Matd sigma = ren2014_cauchy_stress(deformation[i], inverse_plastic[i], shear, bulk, &sigma_eq, nullptr);
            zz_sum += sigma(2, 2);
            eq_sum += sigma_eq;
            p_sum += inelastic[i];
        }
        if (point.core_n > 0)
        {
            point.core_cauchy_zz = zz_sum / static_cast<Real>(point.core_n);
            point.core_sigma_eq = eq_sum / static_cast<Real>(point.core_n);
            point.core_inelastic = p_sum / static_cast<Real>(point.core_n);
        }
        for (size_t i = 0; i < n_particles; ++i)
            point.max_inelastic = std::max(point.max_inelastic, inelastic[i]);
        point.ke = kinetic_energy(particles);
        point.speed = max_speed(particles);
        Vecd *force = particles.getVariableDataByName<Vecd>("Force");
        Vecd *force_prior = particles.getVariableDataByName<Vecd>("ForcePrior");
        for (size_t i : top_part.body_part_particles_)
            point.applied_z += force_prior[i][2];
        Real bottom_internal_z = 0.0;
        for (size_t i : bottom_part.body_part_particles_)
            bottom_internal_z += force[i][2] + force_prior[i][2];
        point.reaction_z = -bottom_internal_z;
        return point;
    };

    struct RelaxOutcome
    {
        int aborted = 0;
        int equilibrated = 0;
        size_t steps = 0;
        Real pseudo = 0.0;
    };

    auto relax = [&](Real cap, Real ramp_end, const char *label) -> RelaxOutcome
    {
        RelaxOutcome outcome;
        Real previous_eps = sample().eps_z;
        int calm_windows = 0;
        const size_t report_every = 500;
        while (outcome.pseudo < cap)
        {
            const Real dt = acoustic.exec();
            if (!(dt > 0.0) || dt > 1.0e-4)
            {
                std::cerr << "ck_hold " << label << " bad acoustic dt=" << dt << "\n";
                outcome.aborted = 1;
                break;
            }
            const Real scale = (ramp_end > 0.0 && outcome.pseudo < ramp_end) ? outcome.pseudo / ramp_end : 1.0;
            force_scale->setValue(scale);
            apply_force.exec(dt);
            step_1st.exec(dt);
            support.exec();
            remove_rigid_xy(particles, com0);
            step_2nd.exec(dt);
            support.exec();
            remove_rigid_xy(particles, com0);
            outcome.pseudo += dt;
            ++outcome.steps;
            if (sample().speed > 20.0)
            {
                std::cerr << "ck_hold " << label << " abort: max speed exceeded 20 m/s\n";
                outcome.aborted = 1;
                break;
            }
            if (outcome.steps % report_every == 0)
            {
                const CkHoldPoint point = sample();
                const Real deps = std::abs(point.eps_z - previous_eps);
                std::cout << label << " N=" << outcome.steps
                          << " pseudo_s=" << outcome.pseudo
                          << " hold_physical_s=unchanged"
                          << " dt=" << dt
                          << " scale=" << scale
                          << " eps_z=" << point.eps_z
                          << " ke_J=" << point.ke
                          << " vmax=" << point.speed << "\n";
                // Residual acoustic ringing at dp = 1 mm is about 2e-7 in strain.
                // That is smaller than the discretization error on eps_z ~ 1e-4.
                if (outcome.pseudo > ramp_end && point.ke < 1.0e-6 && point.speed < 2.0e-3 && deps < 5.0e-7)
                {
                    ++calm_windows;
                    if (calm_windows >= 3)
                    {
                        outcome.equilibrated = 1;
                        break;
                    }
                }
                else
                    calm_windows = 0;
                previous_eps = point.eps_z;
            }
        }
        return outcome;
    };

    recorder.writeToFile();
    TickCount t1 = TickCount::now();
    const RelaxOutcome elastic = relax(pseudo_cap, ramp_pseudo, "elastic");
    const CkHoldPoint at_elastic = sample();
    Real hold_time = 0.0;
    Real inelastic_after_commit = 0.0;
    Real scalar_inelastic = 0.0;
    Real inelastic_drift = 0.0;
    int commit_accepted = 0;
    int post_equilibrated = 0;
    RelaxOutcome post;
    int status = 0;
    if (elastic.aborted || !elastic.equilibrated)
        status = 1;

    if (status == 0 && hold_steps > 0)
    {
        for (int step = 0; step < hold_steps; ++step)
        {
            const MechanicalSnapshot saved = capture_mechanics(particles);
            const CkHoldPoint before = sample();
            const Ren2014PhysicalCommit committed =
                commit_ren2014_physical_increment(particles, shear, demo_temperature, physical_dt);
            if (!committed.accepted)
            {
                restore_mechanics(particles, saved);
                std::cerr << "ck_hold physical increment rejected. large_strain=" << committed.large_strain
                          << " rejected=" << committed.rejected << "\n";
                status = 2;
                break;
            }
            inelastic_after_commit = sample().core_inelastic;
            const spdb::EquivalentReturn scalar = spdb::equivalent_stress_backward_euler(
                before.core_sigma_eq, shear, physical_dt, demo_temperature, before.core_inelastic,
                spdb::StartupPolicy::EndpointExtend);
            scalar_inelastic = scalar.inelastic_eq;
            post = relax(pseudo_cap, 0.0, "post");
            const CkHoldPoint after = sample();
            inelastic_drift = std::abs(after.core_inelastic - inelastic_after_commit);
            if (post.aborted || !post.equilibrated || inelastic_drift > 1.0e-12)
            {
                restore_mechanics(particles, saved);
                std::cerr << "ck_hold post-increment relaxation was not accepted. drift=" << inelastic_drift << "\n";
                status = 3;
                break;
            }
            hold_time += physical_dt;
            commit_accepted = 1;
            post_equilibrated = 1;
            std::cout << "accepted hold_physical_s=" << hold_time
                      << " core_inelastic=" << after.core_inelastic
                      << " scalar_inelastic=" << scalar_inelastic << "\n";
        }
    }
    const CkHoldPoint final_point = sample();
    TickCount t2 = TickCount::now();
    recorder.writeToFile();

    char metrics_name[192];
    std::snprintf(metrics_name, sizeof(metrics_name), "spdb_ck_hold_dp%.4g_damp%.4g_metrics.txt", dp, damping_factor);
    const Real err_z = std::abs(eps_z_hooke) > 0.0 ? std::abs(at_elastic.eps_z - eps_z_hooke) / std::abs(eps_z_hooke) : 0.0;
    const Real err_z_interior =
        std::abs(eps_z_hooke) > 0.0 ? std::abs(at_elastic.eps_z_interior - eps_z_hooke) / std::abs(eps_z_hooke) : 0.0;
    const Real err_stress = nominal_pressure > 0.0 ? std::abs(at_elastic.core_cauchy_zz - (-nominal_pressure)) / nominal_pressure : 0.0;
    const Real err_force = force_total > 0.0 ? std::abs(std::abs(final_point.reaction_z) - force_total) / force_total : 0.0;

    std::cout << "ck_hold wall time: " << (t2 - t1).seconds() << " s\n"
              << "  elastic_equilibrated=" << elastic.equilibrated
              << " elastic_pseudo_s=" << elastic.pseudo
              << " elastic_steps=" << elastic.steps << "\n"
              << "  eps_z=" << at_elastic.eps_z << " hooke=" << eps_z_hooke << " rel_err=" << err_z << "\n"
              << "  eps_z_interior=" << at_elastic.eps_z_interior << " rel_err=" << err_z_interior << "\n"
              << "  core_cauchy_zz_Pa=" << at_elastic.core_cauchy_zz << " rel_err=" << err_stress
              << " core_n=" << at_elastic.core_n << "\n"
              << "  hold_physical_s=" << hold_time
              << " commit_accepted=" << commit_accepted
              << " post_equilibrated=" << post_equilibrated << "\n"
              << "  core_inelastic_after_commit=" << inelastic_after_commit
              << " scalar_at_core_mean_stress=" << scalar_inelastic
              << " inelastic_drift_during_second_relax=" << inelastic_drift << "\n"
              << "  final_eps_z=" << final_point.eps_z
              << " final_core_inelastic=" << final_point.core_inelastic << "\n"
              << "  reaction_rel_err=" << err_force
              << " status=" << status << "\n"
              << "  metrics_file=" << metrics_name << "\n";

    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "mode=ck_hold\n";
    out << "units=SI\n";
    out << "material_source=demo_assumption\n";
    out << "material_note=NOT_950C_20Cr13_process_prediction\n";
    out << "constitutive=Ren2014X20Cr13Solid_StructureIntegration1stHalf\n";
    out << "inelastic_during_relaxation=false\n";
    out << "physical_increment=commit_ren2014_physical_increment_once_per_accepted_step\n";
    out << "youngs_modulus_Pa=" << demo_youngs << "\n";
    out << "poisson=" << poisson << "\n";
    out << "temperature_K=" << demo_temperature << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "numerical_damping_factor=" << damping_factor << "\n";
    out << "acoustic_cfl=0.3\n";
    out << "pseudo_cap_s=" << pseudo_cap << "\n";
    out << "ramp_pseudo_s=" << ramp_pseudo << "\n";
    out << "elastic_pseudo_s=" << elastic.pseudo << "\n";
    out << "elastic_steps=" << elastic.steps << "\n";
    out << "elastic_equilibrated=" << elastic.equilibrated << "\n";
    out << "elastic_eps_z=" << at_elastic.eps_z << "\n";
    out << "hooke_eps_z=" << eps_z_hooke << "\n";
    out << "elastic_eps_z_rel_err=" << err_z << "\n";
    out << "elastic_eps_z_interior=" << at_elastic.eps_z_interior << "\n";
    out << "elastic_eps_z_interior_rel_err=" << err_z_interior << "\n";
    out << "elastic_eps_r=" << at_elastic.eps_r << "\n";
    out << "elastic_core_cauchy_zz_Pa=" << at_elastic.core_cauchy_zz << "\n";
    out << "elastic_core_cauchy_zz_rel_err=" << err_stress << "\n";
    out << "elastic_core_sigma_eq_Pa=" << at_elastic.core_sigma_eq << "\n";
    out << "core_n=" << at_elastic.core_n << "\n";
    out << "physical_dt_s=" << physical_dt << "\n";
    out << "hold_steps_requested=" << hold_steps << "\n";
    out << "hold_physical_s=" << hold_time << "\n";
    out << "hold_is_not_7200s=1\n";
    out << "commit_accepted=" << commit_accepted << "\n";
    out << "post_equilibrated=" << post_equilibrated << "\n";
    out << "post_pseudo_s=" << post.pseudo << "\n";
    out << "core_inelastic_after_commit=" << inelastic_after_commit << "\n";
    out << "scalar_inelastic_at_core_mean_stress=" << scalar_inelastic << "\n";
    out << "inelastic_drift_during_second_relax=" << inelastic_drift << "\n";
    out << "final_eps_z=" << final_point.eps_z << "\n";
    out << "final_core_inelastic=" << final_point.core_inelastic << "\n";
    out << "final_max_inelastic=" << final_point.max_inelastic << "\n";
    out << "applied_force_z_N=" << final_point.applied_z << "\n";
    out << "reaction_z_N=" << final_point.reaction_z << "\n";
    out << "reaction_rel_err=" << err_force << "\n";
    out << "kinetic_energy_J=" << final_point.ke << "\n";
    out << "status=" << status << "\n";
    out << "wall_time_s=" << (t2 - t1).seconds() << "\n";
    return status;
}

//----------------------------------------------------------------------
// Two flush cylinders on StructureIntegration1stHalf plus RepulsionForceCK.
// E = 100 GPa is the demo assumption. Penetration is measured here, not
// extrapolated from the E = 2 GPa synthetic contact runs.
//----------------------------------------------------------------------
int run_ck_contact(Real dp, Real pseudo_cap, Real physical_dt, Real damping_factor, Real stiffness_scale,
                   int hold_steps)
{
    const Real demo_youngs = 100.0e9;
    const Real demo_temperature = 1223.15;
    const Real force_total = nominal_axial_force();
    const Real eps_z_hooke = -nominal_pressure / demo_youngs;
    const Real ramp_pseudo = 0.25 * pseudo_cap;
    const Real slab = 0.55 * dp;
    particle_spacing_ref = dp;

    std::cout << "SPDB ck_contact\n"
              << "  model=ren2014_x20cr13  E_Pa=" << demo_youngs
              << " source=demo_assumption  T_K=" << demo_temperature << "\n"
              << "  NOT a 950 C 20Cr13 process prediction.\n"
              << "  contact=RepulsionForceCK  stiffness_scale=" << stiffness_scale << "\n"
              << "  nominal_force_N=" << force_total
              << "  hooke_eps_z_each_body=" << eps_z_hooke << "\n"
              << "  dp_m=" << dp << " pseudo_cap_s=" << pseudo_cap
              << " physical_dt_s=" << physical_dt
              << " hold_steps=" << hold_steps
              << " damping_factor=" << damping_factor << "\n";

    SPHSystem sph_system(system_bounds_for(dp), dp);
    char output_folder[192];
    std::snprintf(output_folder, sizeof(output_folder), "./output_ck_contact_dp%.4g_k%.4g", dp, stiffness_scale);
    IO::getEnvironment().resetOutputFolder(output_folder);

    auto &lower_shape = sph_system.addShape<LowerWorkpiece>(SlotMode::None, "CkContactLower");
    auto &upper_shape = sph_system.addShape<UpperWorkpiece>(SlotMode::None, "CkContactUpper");
    auto &lower = sph_system.addBody<SolidBody>(lower_shape);
    auto &upper = sph_system.addBody<SolidBody>(upper_shape);
    lower.defineAdaptationRatios(1.15, 1.0);
    upper.defineAdaptationRatios(1.15, 1.0);
    auto &lower_material = lower.defineMatterMaterial<ContactScaledRen2014>(
        rho0_s, demo_youngs, poisson, demo_temperature, false, stiffness_scale);
    auto &upper_material = upper.defineMatterMaterial<ContactScaledRen2014>(
        rho0_s, demo_youngs, poisson, demo_temperature, false, stiffness_scale);
    lower.generateParticles<BaseParticles, Lattice>();
    upper.generateParticles<BaseParticles, Lattice>();
    report_body_metrics("CkContactLower", lower, theoretical_lower_volume(false));
    report_body_metrics("CkContactUpper", upper, theoretical_upper_volume());

    auto &lower_inner = sph_system.addInnerRelation(lower, ConfigType::Lagrangian);
    auto &upper_inner = sph_system.addInnerRelation(upper, ConfigType::Lagrangian);
    auto &lower_upper = sph_system.addContactRelation(lower, upper);
    auto &upper_lower = sph_system.addContactRelation(upper, lower);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    ParticleDynamicsGroup update_contact;
    update_contact.add(&main_methods.addCellLinkedListDynamics(lower));
    update_contact.add(&main_methods.addCellLinkedListDynamics(upper));
    update_contact.add(&main_methods.addRelationDynamics(lower_upper));
    update_contact.add(&main_methods.addRelationDynamics(upper_lower));
    auto &inner_configuration = main_methods.addRelationDynamics(lower_inner, upper_inner);
    auto &correction_l = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(lower_inner);
    auto &correction_u = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(upper_inner);
    auto &step_l_1st = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, Ren2014X20Cr13Solid, LinearCorrectionCK>(lower_inner, damping_factor);
    auto &step_l_2nd = main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(lower_inner);
    auto &step_u_1st = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, Ren2014X20Cr13Solid, LinearCorrectionCK>(upper_inner, damping_factor);
    auto &step_u_2nd = main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(upper_inner);
    auto &acoustic_l = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(lower, 0.3);
    auto &acoustic_u = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(upper, 0.3);

    // Shape level-set normals are far too expensive at this spacing. The cylinders
    // are axis-aligned, so the outward normal is the nearest face of the reference box
    // in cylindrical coordinates. Small-strain hold does not rotate these normals.
    auto set_cylinder_normals = [](BaseParticles &particles, Real z_lo, Real z_hi, Real radius)
    {
        Vecd *pos = particles.ParticlePositions();
        Vecd *normal = particles.registerStateVariableData<Vecd>("NormalDirection");
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        {
            const Real x = pos[i][0];
            const Real y = pos[i][1];
            const Real z = pos[i][2];
            const Real radial = std::sqrt(x * x + y * y);
            const Real dist_top = std::abs(z_hi - z);
            const Real dist_bottom = std::abs(z - z_lo);
            const Real dist_side = std::abs(radius - radial);
            if (dist_top <= dist_bottom && dist_top <= dist_side)
                normal[i] = Vecd(0.0, 0.0, 1.0);
            else if (dist_bottom <= dist_side)
                normal[i] = Vecd(0.0, 0.0, -1.0);
            else if (radial > 0.0)
                normal[i] = Vecd(x / radial, y / radial, 0.0);
            else
                normal[i] = Vecd(0.0, 0.0, 1.0);
        }
    };
    set_cylinder_normals(lower.getBaseParticles(), -cylinder_height, 0.0, cylinder_radius);
    set_cylinder_normals(upper.getBaseParticles(), 0.0, cylinder_height, cylinder_radius);
    auto &factor_l = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(lower_upper);
    auto &factor_u = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(upper_lower);
    auto &repulsion_l = main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::RepulsionForceCK>(lower_upper);
    auto &repulsion_u = main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::RepulsionForceCK>(upper_lower);

    auto &bottom_shape = sph_system.addShape<GeometricShapeBox>(
        Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "BottomAxialSupport");
    auto &top_shape = sph_system.addShape<GeometricShapeBox>(
        Transform(Vec3d(0.0, 0.0, cylinder_height - 0.5 * dp)),
        Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "TopAxialLoad");
    BodyRegionByParticle bottom_part(lower, bottom_shape);
    BodyRegionByParticle top_part(upper, top_shape);
    if (bottom_part.body_part_particles_.empty() || top_part.body_part_particles_.empty())
    {
        std::cerr << "ck_contact: empty support or load layer\n";
        return 1;
    }
    SingleVariable<Real> *force_scale = sph_system.registerSystemVariable<Real>("AxialForceScale", 0.0);
    const Real force_per_particle = -force_total / static_cast<Real>(top_part.body_part_particles_.size());
    auto &support = main_methods.addStateDynamics<AxialZSupportCK<BodyRegionByParticle>>(bottom_part);
    auto &apply_force = main_methods.addStateDynamics<NominalAxialForceCK<BodyRegionByParticle>>(
        top_part, force_per_particle, *force_scale);

    std::cout << "ck_contact configuring contact\n" << std::flush;
    update_contact.exec();
    std::cout << "ck_contact configuring inner relations\n" << std::flush;
    inner_configuration.exec();
    std::cout << "ck_contact correction matrices\n" << std::flush;
    correction_l.exec();
    correction_u.exec();
    std::cout << "ck_contact configuration done\n" << std::flush;

    BaseParticles &lower_particles = lower.getBaseParticles();
    BaseParticles &upper_particles = upper.getBaseParticles();
    Vecd *pos_l = lower_particles.ParticlePositions();
    Vecd *pos_u = upper_particles.ParticlePositions();
    IndexVector lower_face = ids_near_z(pos_l, lower_particles.TotalRealParticles(), -0.5 * dp, 0.5 * dp);
    IndexVector upper_face = ids_near_z(pos_u, upper_particles.TotalRealParticles(), 0.5 * dp, 0.5 * dp);
    if (lower_face.empty() || upper_face.empty())
    {
        std::cerr << "ck_contact: empty interface layer\n";
        return 1;
    }

    auto com_of = [](BaseParticles &particles) -> Vecd
    {
        Vecd *pos = particles.ParticlePositions();
        Real *mass = particles.getVariableDataByName<Real>("Mass");
        Real mass_sum = 0.0;
        Vecd com = Vecd::Zero();
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        {
            mass_sum += mass[i];
            com += mass[i] * pos[i];
        }
        if (!(mass_sum > 0.0))
            return Vecd::Zero();
        com /= mass_sum;
        return com;
    };
    const Vecd com_l0 = com_of(lower_particles);
    const Vecd com_u0 = com_of(upper_particles);
    const Real z_lb0 = mean_coordinate(bottom_part.body_part_particles_, pos_l, 2);
    const Real z_lf0 = mean_coordinate(lower_face, pos_l, 2);
    const Real z_uf0 = mean_coordinate(upper_face, pos_u, 2);
    const Real z_ut0 = mean_coordinate(top_part.body_part_particles_, pos_u, 2);
    const Real gauge_l0 = z_lf0 - z_lb0;
    const Real gauge_u0 = z_ut0 - z_uf0;

    auto strain_of = [](Real z_a, Real z_b, Real gauge0) -> Real
    {
        return gauge0 != 0.0 ? ((z_b - z_a) - gauge0) / gauge0 : 0.0;
    };
    auto sample_eps = [&]() -> std::pair<Real, Real>
    {
        return {strain_of(mean_coordinate(bottom_part.body_part_particles_, pos_l, 2),
                          mean_coordinate(lower_face, pos_l, 2), gauge_l0),
                strain_of(mean_coordinate(upper_face, pos_u, 2),
                          mean_coordinate(top_part.body_part_particles_, pos_u, 2), gauge_u0)};
    };

    struct RelaxOutcome
    {
        int aborted = 0;
        int equilibrated = 0;
        size_t steps = 0;
        Real pseudo = 0.0;
    };
    auto relax = [&](Real cap, Real ramp_end, const char *label) -> RelaxOutcome
    {
        RelaxOutcome outcome;
        Real previous_eps = sample_eps().first;
        Real previous_penetration = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
        int calm_windows = 0;
        std::cout << label << " relax start\n" << std::flush;
        while (outcome.pseudo < cap)
        {
            const Real dt = std::min(acoustic_l.exec(), acoustic_u.exec());
            if (!(dt > 0.0) || dt > 1.0e-4)
            {
                std::cerr << "ck_contact " << label << " bad acoustic dt=" << dt << "\n";
                outcome.aborted = 1;
                break;
            }
            const Real scale = (ramp_end > 0.0 && outcome.pseudo < ramp_end) ? outcome.pseudo / ramp_end : 1.0;
            force_scale->setValue(scale);
            // Displacement per acoustic step is far below dp, so the contact
            // neighbor list is refreshed every 25 steps rather than every step.
            if (outcome.steps % 25 == 0)
                update_contact.exec();
            factor_l.exec();
            repulsion_l.exec();
            factor_u.exec();
            repulsion_u.exec();
            apply_force.exec(dt);
            step_l_1st.exec(dt);
            support.exec();
            remove_rigid_xy(lower_particles, com_l0);
            step_l_2nd.exec(dt);
            support.exec();
            remove_rigid_xy(lower_particles, com_l0);
            step_u_1st.exec(dt);
            remove_rigid_xy(upper_particles, com_u0);
            step_u_2nd.exec(dt);
            remove_rigid_xy(upper_particles, com_u0);
            outcome.pseudo += dt;
            ++outcome.steps;
            if (max_speed(lower_particles) > 20.0 || max_speed(upper_particles) > 20.0)
            {
                std::cerr << "ck_contact " << label << " abort: max speed exceeded 20 m/s\n";
                outcome.aborted = 1;
                break;
            }
            if (outcome.steps % 100 == 0)
            {
                const auto eps = sample_eps();
                const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
                const Real ke = kinetic_energy(lower_particles) + kinetic_energy(upper_particles);
                const Real speed = std::max(max_speed(lower_particles), max_speed(upper_particles));
                std::cout << label << " N=" << outcome.steps
                          << " pseudo_s=" << outcome.pseudo
                          << " dt=" << dt
                          << " scale=" << scale
                          << " eps_l=" << eps.first
                          << " eps_u=" << eps.second
                          << " max_penetration_m=" << gap.max_penetration
                          << " ke_J=" << ke
                          << " vmax=" << speed << "\n"
                          << std::flush;
                const bool calm = outcome.pseudo > ramp_end && ke < 1.0e-5 && speed < 2.0e-3 &&
                                  std::abs(eps.first - previous_eps) < 5.0e-7 &&
                                  std::abs(gap.max_penetration - previous_penetration) < 2.0e-6;
                calm_windows = calm ? calm_windows + 1 : 0;
                if (calm_windows >= 15)
                {
                    outcome.equilibrated = 1;
                    break;
                }
                previous_eps = eps.first;
                previous_penetration = gap.max_penetration;
            }
        }
        return outcome;
    };

    TickCount t1 = TickCount::now();
    const RelaxOutcome elastic = relax(pseudo_cap, ramp_pseudo, "elastic");
    const auto eps_elastic = sample_eps();
    const InterfaceGap gap_elastic = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    Real hold_time = 0.0;
    Real penetration_after_hold = gap_elastic.max_penetration;
    int status = (elastic.aborted || !elastic.equilibrated) ? 1 : 0;
    if (status == 0 && hold_steps > 0)
    {
        for (int step = 0; step < hold_steps; ++step)
        {
            const MechanicalSnapshot saved_l = capture_mechanics(lower_particles);
            const MechanicalSnapshot saved_u = capture_mechanics(upper_particles);
            const Real penetration_before = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
            const Ren2014PhysicalCommit committed_l =
                commit_ren2014_physical_increment(lower_particles, lower_material.ShearModulus(), demo_temperature, physical_dt);
            const Ren2014PhysicalCommit committed_u =
                commit_ren2014_physical_increment(upper_particles, upper_material.ShearModulus(), demo_temperature, physical_dt);
            if (!committed_l.accepted || !committed_u.accepted)
            {
                restore_mechanics(lower_particles, saved_l);
                restore_mechanics(upper_particles, saved_u);
                std::cerr << "ck_contact physical increment rejected\n";
                status = 2;
                break;
            }
            const Real inelastic_l = committed_l.mean_inelastic;
            const Real inelastic_u = committed_u.mean_inelastic;
            const RelaxOutcome post = relax(pseudo_cap, 0.0, "post");
            const Real inelastic_l_after = [&]()
            {
                Real *field = lower_particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
                Real sum = 0.0;
                const size_t n = lower_particles.TotalRealParticles();
                for (size_t i = 0; i < n; ++i)
                    sum += field[i];
                return n > 0 ? sum / static_cast<Real>(n) : 0.0;
            }();
            const Real inelastic_u_after = [&]()
            {
                Real *field = upper_particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
                Real sum = 0.0;
                const size_t n = upper_particles.TotalRealParticles();
                for (size_t i = 0; i < n; ++i)
                    sum += field[i];
                return n > 0 ? sum / static_cast<Real>(n) : 0.0;
            }();
            if (post.aborted || !post.equilibrated || std::abs(inelastic_l_after - inelastic_l) > 1.0e-12 ||
                std::abs(inelastic_u_after - inelastic_u) > 1.0e-12)
            {
                restore_mechanics(lower_particles, saved_l);
                restore_mechanics(upper_particles, saved_u);
                std::cerr << "ck_contact post-increment relaxation was not accepted\n";
                status = 3;
                break;
            }
            hold_time += physical_dt;
            penetration_after_hold = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
            std::cout << "accepted hold_physical_s=" << hold_time
                      << " penetration_before_m=" << penetration_before
                      << " penetration_after_m=" << penetration_after_hold
                      << " mean_inelastic_l=" << inelastic_l
                      << " mean_inelastic_u=" << inelastic_u << "\n";
        }
    }

    const auto eps_final = sample_eps();
    const InterfaceGap gap_final = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const Real contact_l_z = sum_field_z(lower_particles, ("RepulsionForce" + lower_upper.Name()).c_str());
    const Real contact_u_z = sum_field_z(upper_particles, ("RepulsionForce" + upper_lower.Name()).c_str());
    Vecd *force_l = lower_particles.getVariableDataByName<Vecd>("Force");
    Vecd *prior_l = lower_particles.getVariableDataByName<Vecd>("ForcePrior");
    Real bottom_internal_z = 0.0;
    for (size_t i : bottom_part.body_part_particles_)
        bottom_internal_z += force_l[i][2] + prior_l[i][2];
    const Real reaction_z = -bottom_internal_z;
    const Real applied_z = [&]()
    {
        Vecd *prior = upper_particles.getVariableDataByName<Vecd>("ForcePrior");
        Real sum = 0.0;
        for (size_t i : top_part.body_part_particles_)
            sum += prior[i][2];
        return sum;
    }();
    TickCount t2 = TickCount::now();
    const Real material_shortening = (z_lf0 - mean_coordinate(lower_face, pos_l, 2)) +
                                     (mean_coordinate(top_part.body_part_particles_, pos_u, 2) - z_ut0);
    const Real end_displacement = mean_coordinate(top_part.body_part_particles_, pos_u, 2) - z_ut0;
    const Real contact_imbalance = force_total > 0.0 ? std::abs(contact_l_z + contact_u_z) / force_total : 0.0;
    const Real reaction_err = force_total > 0.0 ? std::abs(std::abs(reaction_z) - force_total) / force_total : 0.0;

    char metrics_name[192];
    std::snprintf(metrics_name, sizeof(metrics_name), "spdb_ck_contact_dp%.4g_k%.4g_metrics.txt", dp, stiffness_scale);
    std::cout << "ck_contact wall time: " << (t2 - t1).seconds() << " s\n"
              << "  elastic_equilibrated=" << elastic.equilibrated
              << " elastic_pseudo_s=" << elastic.pseudo
              << " elastic_steps=" << elastic.steps << "\n"
              << "  eps_l=" << eps_elastic.first << " eps_u=" << eps_elastic.second
              << " hooke_each=" << eps_z_hooke << "\n"
              << "  max_penetration_m=" << gap_elastic.max_penetration
              << " penetration_over_dp=" << gap_elastic.max_penetration / dp
              << " pairs=" << gap_elastic.penetrating << "/" << gap_elastic.pairs << "\n"
              << "  hold_physical_s=" << hold_time
              << " penetration_after_hold_m=" << penetration_after_hold
              << " penetration_drift_m=" << (penetration_after_hold - gap_elastic.max_penetration) << "\n"
              << "  material_shortening_m=" << material_shortening
              << " top_displacement_m=" << end_displacement << "\n"
              << "  contact_lower_z_N=" << contact_l_z
              << " contact_upper_z_N=" << contact_u_z
              << " reaction_z_N=" << reaction_z
              << " reaction_rel_err=" << reaction_err
              << " contact_pair_imbalance=" << contact_imbalance << "\n"
              << "  contact_stiffness_Pa=" << lower_material.ContactStiffness()
              << " status=" << status << "\n"
              << "  metrics_file=" << metrics_name << "\n";

    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "mode=ck_contact\n";
    out << "units=SI\n";
    out << "material_source=demo_assumption\n";
    out << "material_note=NOT_950C_20Cr13_process_prediction\n";
    out << "constitutive=ContactScaledRen2014_StructureIntegration1stHalf_RepulsionForceCK\n";
    out << "youngs_modulus_Pa=" << demo_youngs << "\n";
    out << "poisson=" << poisson << "\n";
    out << "temperature_K=" << demo_temperature << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "contact_stiffness_scale=" << stiffness_scale << "\n";
    out << "contact_stiffness_Pa=" << lower_material.ContactStiffness() << "\n";
    out << "upper_contact_stiffness_Pa=" << upper_material.ContactStiffness() << "\n";
    out << "numerical_damping_factor=" << damping_factor << "\n";
    out << "pseudo_cap_s=" << pseudo_cap << "\n";
    out << "elastic_pseudo_s=" << elastic.pseudo << "\n";
    out << "elastic_steps=" << elastic.steps << "\n";
    out << "elastic_equilibrated=" << elastic.equilibrated << "\n";
    out << "elastic_eps_l=" << eps_elastic.first << "\n";
    out << "elastic_eps_u=" << eps_elastic.second << "\n";
    out << "hooke_eps_z_each=" << eps_z_hooke << "\n";
    out << "max_penetration_m=" << gap_elastic.max_penetration << "\n";
    out << "max_penetration_over_dp=" << (dp > 0.0 ? gap_elastic.max_penetration / dp : 0.0) << "\n";
    out << "min_surface_gap_m=" << gap_elastic.min_surface_gap << "\n";
    out << "penetrating_pairs=" << gap_elastic.penetrating << "\n";
    out << "interface_pairs=" << gap_elastic.pairs << "\n";
    out << "hold_physical_s=" << hold_time << "\n";
    out << "penetration_after_hold_m=" << penetration_after_hold << "\n";
    out << "penetration_drift_m=" << (penetration_after_hold - gap_elastic.max_penetration) << "\n";
    out << "final_eps_l=" << eps_final.first << "\n";
    out << "final_eps_u=" << eps_final.second << "\n";
    out << "final_max_penetration_m=" << gap_final.max_penetration << "\n";
    out << "material_shortening_m=" << material_shortening << "\n";
    out << "top_displacement_m=" << end_displacement << "\n";
    out << "applied_force_z_N=" << applied_z << "\n";
    out << "contact_lower_z_N=" << contact_l_z << "\n";
    out << "contact_upper_z_N=" << contact_u_z << "\n";
    out << "reaction_z_N=" << reaction_z << "\n";
    out << "reaction_rel_err=" << reaction_err << "\n";
    out << "contact_pair_imbalance=" << contact_imbalance << "\n";
    out << "status=" << status << "\n";
    out << "wall_time_s=" << (t2 - t1).seconds() << "\n";
    return status;
}
} // namespace

int main(int ac, char *av[])
{
    RunOptions opt = parse_options(ac, av);
    if (opt.material_point)
        return spdb::run_ren2014_material_point();
    const int mode_count = static_cast<int>(opt.elastic_axial) + static_cast<int>(opt.elastic_contact) +
                           static_cast<int>(opt.slot_contact) + static_cast<int>(opt.slot_hold) +
                           static_cast<int>(opt.ck_hold) + static_cast<int>(opt.ck_contact) +
                           static_cast<int>(opt.material_point) + static_cast<int>(opt.elastic_no_slot) +
                           static_cast<int>(opt.geometry_only);
    if (mode_count > 1)
    {
        std::cerr << "Choose only one run mode\n";
        return 2;
    }
    if (opt.dp <= 0.0)
    {
        std::cerr << "--dp must be positive\n";
        return 2;
    }
    if (opt.ck_hold)
    {
        if (!(opt.physical_dt > 0.0) || !(opt.damping_factor > 0.0) || opt.hold_steps < 0)
        {
            std::cerr << "--physical_dt and --damping_factor must be positive; --hold_steps must be >= 0\n";
            return 2;
        }
        return run_ck_hold(opt.dp, opt.end_time > 0.0 ? opt.end_time : 3.0e-3, opt.physical_dt, opt.damping_factor,
                           opt.hold_steps);
    }
    if (opt.ck_contact)
    {
        if (!(opt.physical_dt > 0.0) || !(opt.damping_factor > 0.0) || !(opt.contact_stiffness_scale > 0.0) ||
            opt.hold_steps < 0)
        {
            std::cerr << "--physical_dt, --damping_factor and --contact_stiffness_scale must be positive\n";
            return 2;
        }
        return run_ck_contact(opt.dp, opt.end_time > 0.0 ? opt.end_time : 3.0e-3, opt.physical_dt, opt.damping_factor,
                              opt.contact_stiffness_scale, opt.hold_steps);
    }
    if (opt.elastic_axial)
        return run_elastic_axial(opt.dp, opt.end_time > 0.0 ? opt.end_time : 8.0e-3);
    if (opt.slot_contact)
    {
        if (!(opt.contact_stiffness_scale > 0.0))
        {
            std::cerr << "--contact_stiffness_scale must be positive\n";
            return 2;
        }
        const Real youngs = opt.contact_youngs > 0.0 ? opt.contact_youngs : 100.0e9;
        return run_elastic_contact(opt.dp, opt.end_time > 0.0 ? opt.end_time : 1.5e-3, opt.contact_stiffness_scale,
                                   youngs, SlotMode::LowerOnly);
    }
    if (opt.slot_hold)
    {
        if (!(opt.physical_dt > 0.0) || !(opt.contact_stiffness_scale > 0.0) || opt.hold_steps < 0)
        {
            std::cerr << "--physical_dt and --contact_stiffness_scale must be positive; --hold_steps >= 0\n";
            return 2;
        }
        return run_slot_hold(opt.dp, opt.end_time > 0.0 ? opt.end_time : 1.5e-3, opt.physical_dt,
                             opt.contact_stiffness_scale, opt.hold_steps);
    }
    if (opt.elastic_contact)
    {
        if (!(opt.contact_stiffness_scale > 0.0))
        {
            std::cerr << "--contact_stiffness_scale must be positive\n";
            return 2;
        }
        const Real youngs = opt.contact_youngs > 0.0 ? opt.contact_youngs : Youngs_modulus_axial;
        return run_elastic_contact(opt.dp, opt.end_time > 0.0 ? opt.end_time : 4.0e-3, opt.contact_stiffness_scale,
                                   youngs, SlotMode::None);
    }

    std::cout << "SPDB first-round case\n"
              << "  geometry_only=" << opt.geometry_only
              << " elastic_no_slot=" << opt.elastic_no_slot
              << " slot_mode=" << static_cast<int>(opt.slot_mode) << "\n"
              << "  NOTE: synthetic elastic material — not 950C 20Cr13.\n";

    SPHSystem sph_system(system_domain_bounds, resolution_ref);

    SolidBody lower_body(sph_system, makeShared<LowerWorkpiece>(opt.slot_mode, "LowerWorkpiece"));
    lower_body.defineAdaptationRatios(1.15, 1.0);
    lower_body.defineMatterMaterial<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
    lower_body.generateParticles<BaseParticles, Lattice>();

    SolidBody upper_body(sph_system, makeShared<UpperWorkpiece>(opt.slot_mode, "UpperWorkpiece"));
    upper_body.defineAdaptationRatios(1.15, 1.0);
    upper_body.defineMatterMaterial<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
    upper_body.generateParticles<BaseParticles, Lattice>();

    const bool lower_has_slot = (opt.slot_mode == SlotMode::LowerOnly || opt.slot_mode == SlotMode::Both);
    report_body_metrics("LowerWorkpiece", lower_body, theoretical_lower_volume(lower_has_slot));
    report_body_metrics("UpperWorkpiece", upper_body,
                        opt.slot_mode == SlotMode::Both
                            ? theoretical_full_cylinder_volume() - theoretical_slot_volume()
                            : theoretical_upper_volume());

    BodyStatesRecordingToVtp write_states(sph_system);
    write_metrics_file("spdb_geometry_metrics.txt", lower_body, upper_body, opt.slot_mode, lower_has_slot);

    if (opt.geometry_only)
    {
        sph_system.initializeSystemCellLinkedLists();
        write_states.writeToFile(0);
        std::cout << "geometry_only complete. Wrote VTK and spdb_geometry_metrics.txt\n";
        return 0;
    }

    //----------------------------------------------------------------------
    // Elastic compression baseline (no slot): synthetic material.
    // Bottom of lower fixed; light axial gravity on upper; contact + damping.
    //----------------------------------------------------------------------
    InnerRelation lower_inner(lower_body);
    InnerRelation upper_inner(upper_body);
    SurfaceContactRelation lower_upper_contact(lower_body, {&upper_body});
    SurfaceContactRelation upper_lower_contact(upper_body, {&lower_body});

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_lower(lower_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_upper(upper_inner);

    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_lower_1st(lower_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_lower_2nd(lower_inner);
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_upper_1st(upper_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_upper_2nd(upper_inner);

    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor_lower(lower_upper_contact);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor_upper(upper_lower_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_lower(lower_upper_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_upper(upper_lower_contact);

    // Fix a thin bottom slab of the lower cylinder (end support).
    GeometricShapeBox lower_holder_shape(Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * particle_spacing_ref)),
                                         Vec3d(cylinder_radius, cylinder_radius, particle_spacing_ref),
                                         "LowerHolder");
    BodyRegionByParticle lower_holder(lower_body, lower_holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constrain_lower(lower_holder);

    // Light compressive body force on upper workpiece (synthetic demo load).
    Gravity gravity(Vec3d(0.0, 0.0, -50.0)); // m/s^2, demonstration only
    SimpleDynamics<GravityForce<Gravity>> apply_gravity(upper_body, gravity);

    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        damp_lower(0.2, lower_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        damp_upper(0.2, upper_inner, "Velocity", physical_viscosity);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_lower(lower_body, 0.2);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> dt_upper(upper_body, 0.2);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_lower.exec();
    corrected_upper.exec();
    apply_gravity.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 2.0e-5; // short smoke / baseline run
    Real output_period = end_time / 5.0;
    Real dt = 0.0;

    write_states.writeToFile(0);
    TickCount t1 = TickCount::now();

    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 50 == 0)
            {
                std::cout << "N=" << ite << " Time: " << physical_time << " dt: " << dt << "\n";
            }

            contact_factor_lower.exec();
            contact_force_lower.exec();
            contact_factor_upper.exec();
            contact_force_upper.exec();

            stress_lower_1st.exec(dt);
            constrain_lower.exec(dt);
            damp_lower.exec(dt);
            constrain_lower.exec(dt);
            stress_lower_2nd.exec(dt);

            stress_upper_1st.exec(dt);
            damp_upper.exec(dt);
            stress_upper_2nd.exec(dt);

            ite++;
            dt = SMIN(dt_lower.exec(), dt_upper.exec());
            integration_time += dt;
            physical_time += dt;

            lower_body.updateCellLinkedList();
            upper_body.updateCellLinkedList();
            lower_upper_contact.updateConfiguration();
            upper_lower_contact.updateConfiguration();
        }
        write_states.writeToFile();
    }

    TickCount t2 = TickCount::now();
    std::cout << "elastic_no_slot wall time: " << (t2 - t1).seconds() << " s\n";
    std::cout << "NOTE: synthetic material / short smoke run — not a validated hold-pressure result.\n";
    write_metrics_file("spdb_elastic_no_slot_metrics.txt", lower_body, upper_body, opt.slot_mode, false);
    return 0;
}
