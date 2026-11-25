/**
 * @file spdb_slot_hold.cpp
 * @brief Blind-slot hold: ContactForce + Ren2014 physical increments.
 * @details Extra-test only. Acoustic relaxation keeps inelastic frozen.
 *          Each accepted physical step calls commit_ren2014_physical_increment once
 *          on both bodies. Contact is ContactForce, not RepulsionForceCK.
 */
#include "spdb_slot_hold.h"

#include "spdb_contact_scaled.h"
#include "spdb_geometry.h"
#include "spdb_ren2014_solid.h"
#include "spdb_slot_metrics.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

using namespace SPH;

namespace
{
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

Real sum_field_z(BaseParticles &particles, const char *name)
{
    Vecd *field = particles.getVariableDataByName<Vecd>(name);
    Real sum = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        sum += field[i][2];
    return sum;
}

struct InterfaceGap
{
    size_t pairs = 0;
    size_t penetrating = 0;
    Real min_surface_gap = 0.0;
    Real max_penetration = 0.0;
};

InterfaceGap measure_interface_gap(const IndexVector &lower_ids, Vecd *pos_l, const IndexVector &upper_ids,
                                   Vecd *pos_u, Real dp)
{
    InterfaceGap gap;
    Real min_s = std::numeric_limits<Real>::max();
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
        min_s = std::min(min_s, surface);
        gap.max_penetration = std::max(gap.max_penetration, std::max(Real(0.0), -surface));
        if (surface < 0.0)
            ++gap.penetrating;
        ++gap.pairs;
    }
    if (gap.pairs > 0)
        gap.min_surface_gap = min_s;
    return gap;
}

struct HoldSnapshot
{
    std::vector<Vecd> pos, vel;
    std::vector<Matd> deformation, inverse_plastic;
    std::vector<Real> density, inelastic, range_flag;
};

HoldSnapshot capture(BaseParticles &particles)
{
    const size_t n = particles.TotalRealParticles();
    HoldSnapshot s;
    Vecd *pos = particles.ParticlePositions();
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Matd *F = particles.getVariableDataByName<Matd>("DeformationGradient");
    Matd *Cp = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *flag = particles.getVariableDataByName<Real>("ArrheniusRangeFlag");
    s.pos.assign(pos, pos + n);
    s.vel.assign(vel, vel + n);
    s.deformation.assign(F, F + n);
    s.inverse_plastic.assign(Cp, Cp + n);
    s.density.assign(rho, rho + n);
    s.inelastic.assign(eq, eq + n);
    s.range_flag.assign(flag, flag + n);
    return s;
}

void restore(BaseParticles &particles, const HoldSnapshot &s)
{
    const size_t n = particles.TotalRealParticles();
    Vecd *pos = particles.ParticlePositions();
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Matd *F = particles.getVariableDataByName<Matd>("DeformationGradient");
    Matd *Cp = particles.getVariableDataByName<Matd>("InversePlasticRightCauchyStrain");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *eq = particles.getVariableDataByName<Real>("EquivalentInelasticStrain");
    Real *flag = particles.getVariableDataByName<Real>("ArrheniusRangeFlag");
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = s.pos[i];
        vel[i] = s.vel[i];
        F[i] = s.deformation[i];
        Cp[i] = s.inverse_plastic[i];
        rho[i] = s.density[i];
        eq[i] = s.inelastic[i];
        flag[i] = s.range_flag[i];
    }
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

// Sidecar for disconnected monitoring. Overwrite each update so a short
// `cat hold_status.txt` shows the latest hold clock and phase.
void write_hold_status(const char *phase, Real hold_physical, int hold_step, int hold_steps,
                       Real physical_dt, Real pseudo, Real ke, Real max_penetration,
                       Real cavity_net_height, Real mean_inelastic_l, Real mean_inelastic_u,
                       int aborted = 0, int finished = 0)
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    char iso[64];
    std::strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&tt));

    std::ofstream out("hold_status.txt");
    out << std::setprecision(10);
    out << "updated_at=" << iso << "\n";
    out << "phase=" << phase << "\n";
    out << "hold_physical_s=" << hold_physical << "\n";
    out << "hold_step=" << hold_step << "\n";
    out << "hold_steps_requested=" << hold_steps << "\n";
    out << "physical_dt_s=" << physical_dt << "\n";
    out << "target_hold_s=" << physical_dt * static_cast<Real>(hold_steps) << "\n";
    out << "progress_frac="
        << (hold_steps > 0 ? hold_physical / (physical_dt * static_cast<Real>(hold_steps)) : 0.0) << "\n";
    out << "relax_pseudo_s=" << pseudo << "\n";
    out << "ke_J=" << ke << "\n";
    out << "max_penetration_m=" << max_penetration << "\n";
    out << "cavity_net_height_m=" << cavity_net_height << "\n";
    out << "mean_inelastic_lower=" << mean_inelastic_l << "\n";
    out << "mean_inelastic_upper=" << mean_inelastic_u << "\n";
    out << "aborted=" << aborted << "\n";
    out << "finished=" << finished << "\n";
    out << "hold_is_7200s=" << (std::abs(physical_dt * hold_steps - 7200.0) < 1.0e-6 ? 1 : 0) << "\n";
    out.flush();

    std::ofstream beat("hold_heartbeat");
    beat << iso << " hold_s=" << hold_physical << " phase=" << phase << "\n";
}
} // namespace

int SPH::run_slot_hold(Real dp, Real end_time, Real physical_dt, Real stiffness_scale, int hold_steps)
{
    particle_spacing_ref = dp;
    resolution_ref = dp;
    const Real demo_youngs = 100.0e9;
    const Real demo_temperature = 1223.15;
    const Real force_total = nominal_axial_force();
    const Real ramp_time = std::min(1.0e-3, 0.25 * end_time);
    const Real output_period = std::min(1.0e-3, std::max(0.25 * end_time, 1.0e-6));
    const Real slab = 0.55 * dp;
    const Real relax_cap = std::max(end_time, 1.5e-3);

    std::cout << "SPDB slot_hold\n"
              << "  model=ren2014_x20cr13 E_Pa=" << demo_youngs
              << " source=demo_assumption T_K=" << demo_temperature << "\n"
              << "  NOT a 950 C 20Cr13 process prediction and NOT a 7200 s result.\n"
              << "  contact=ContactForce  inelastic_during_relaxation=false\n"
              << "  slot=lower_blind_cube  friction=none\n"
              << "  nominal_force_N=" << force_total
              << "  dp_m=" << dp
              << "  physical_dt_s=" << physical_dt
              << "  hold_steps=" << hold_steps
              << "  stiffness_scale=" << stiffness_scale << "\n"
              << "  Cavity metrics use particle centers. Volume is not claimed.\n"
              << "  Penetration must not be painted as creep.\n";

    SPHSystem sph_system(system_bounds_for(dp), dp);
    SolidBody lower_body(sph_system, makeShared<LowerWorkpiece>(SlotMode::LowerOnly, "SlotHoldLower"));
    lower_body.defineAdaptationRatios(1.15, 1.0);
    ContactScaledRen2014 &lower_mat = lower_body.defineMatterMaterial<ContactScaledRen2014>(
        rho0_s, demo_youngs, poisson, demo_temperature, false, stiffness_scale);
    lower_body.generateParticles<BaseParticles, Lattice>();
    SolidBody upper_body(sph_system, makeShared<UpperWorkpiece>(SlotMode::LowerOnly, "SlotHoldUpper"));
    upper_body.defineAdaptationRatios(1.15, 1.0);
    ContactScaledRen2014 &upper_mat = upper_body.defineMatterMaterial<ContactScaledRen2014>(
        rho0_s, demo_youngs, poisson, demo_temperature, false, stiffness_scale);
    upper_body.generateParticles<BaseParticles, Lattice>();
    std::cout << "  lower_particles=" << lower_body.getBaseParticles().TotalRealParticles()
              << " upper_particles=" << upper_body.getBaseParticles().TotalRealParticles() << "\n";

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

    GeometricShapeBox bottom_shape(Transform(Vec3d(0.0, 0.0, -cylinder_height + 0.5 * dp)),
                                   Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "BottomAxialSupport");
    GeometricShapeBox top_shape(Transform(Vec3d(0.0, 0.0, cylinder_height - 0.5 * dp)),
                                Vec3d(cylinder_radius + dp, cylinder_radius + dp, slab), "TopAxialLoad");
    BodyRegionByParticle bottom_part(lower_body, bottom_shape);
    BodyRegionByParticle top_part(upper_body, top_shape);
    if (bottom_part.body_part_particles_.empty() || top_part.body_part_particles_.empty())
    {
        std::cerr << "slot_hold: empty support or load layer\n";
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
        spdb::ids_near_z_outside_slot(pos_l, lower_particles.TotalRealParticles(), -0.5 * dp, 0.5 * dp);
    IndexVector upper_face =
        spdb::ids_near_z_outside_slot(pos_u, upper_particles.TotalRealParticles(), 0.5 * dp, 0.5 * dp);
    if (lower_face.empty() || upper_face.empty())
    {
        std::cerr << "slot_hold: empty join-face ring\n";
        return 1;
    }

    const spdb::SlotCavityMetrics cavity0 =
        spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                  upper_particles.TotalRealParticles(), dp);
    BodyStatesRecordingToVtp write_states(sph_system);
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    write_states.writeToFile(0);
    std::cout << "initial_cavity_net_height_m=" << cavity0.center_net_height
              << " floor_n=" << cavity0.floor_n << " roof_n=" << cavity0.roof_n << "\n";

    auto one_acoustic = [&](Real dt)
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
        lower_body.updateCellLinkedList();
        upper_body.updateCellLinkedList();
        lower_inner.updateConfiguration();
        upper_inner.updateConfiguration();
        lower_upper.updateConfiguration();
        upper_lower.updateConfiguration();
    };

    Real hold_physical = 0.0;
    int hold_step_done = 0;
    auto relax_to_balance = [&](Real cap, const char *label) -> int
    {
        Real previous_pen = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
        int calm = 0;
        size_t steps = 0;
        Real pseudo = 0.0;
        Real dt = 0.0;
        while (pseudo < cap)
        {
            dt = SMIN(dt_l.exec(), dt_u.exec());
            if (!(dt > 0.0))
                return 1;
            one_acoustic(dt);
            ++steps;
            pseudo += dt;
            // Acoustic clock only during relax; hold clock is hold_physical.
            if (std::strcmp(label, "elastic") == 0)
                physical_time = pseudo;
            if (max_speed(lower_particles) > 20.0 || max_speed(upper_particles) > 20.0)
            {
                std::cerr << "slot_hold " << label << " abort: max speed > 20 m/s\n";
                write_hold_status("aborted_speed", hold_physical, hold_step_done, hold_steps, physical_dt, pseudo,
                                  kinetic_energy(lower_particles) + kinetic_energy(upper_particles), previous_pen,
                                  0.0, mean_inelastic(lower_particles), mean_inelastic(upper_particles), 1, 0);
                return 1;
            }
            if (steps % 500 == 0)
            {
                corrected_l.exec();
                corrected_u.exec();
                const InterfaceGap gap = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
                const Real ke = kinetic_energy(lower_particles) + kinetic_energy(upper_particles);
                const Real dpen = std::abs(gap.max_penetration - previous_pen);
                const spdb::SlotCavityMetrics cav =
                    spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                              upper_particles.TotalRealParticles(), dp);
                std::cout << label << " N=" << steps << " pseudo_s=" << pseudo << " dt=" << dt
                          << " max_penetration_m=" << gap.max_penetration << " ke_J=" << ke << "\n";
                write_hold_status(label, hold_physical, hold_step_done, hold_steps, physical_dt, pseudo, ke,
                                  gap.max_penetration, cav.center_net_height, mean_inelastic(lower_particles),
                                  mean_inelastic(upper_particles));
                if (pseudo > ramp_time && ke < 1.0e-6 &&
                    max_speed(lower_particles) < 2.0e-3 && max_speed(upper_particles) < 2.0e-3 && dpen < 1.0e-9)
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

    TickCount t1 = TickCount::now();
    if (relax_to_balance(relax_cap, "elastic") != 0)
        return 1;
    corrected_l.exec();
    corrected_u.exec();

    const InterfaceGap gap_elastic = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const spdb::SlotCavityMetrics cavity_elastic =
        spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                  upper_particles.TotalRealParticles(), dp);
    const Real inelastic_l0 = mean_inelastic(lower_particles);
    const Real inelastic_u0 = mean_inelastic(upper_particles);
    std::cout << "elastic_balance max_penetration_m=" << gap_elastic.max_penetration
              << " cavity_net_height_m=" << cavity_elastic.center_net_height
              << " top_wall_defl_m=" << cavity_elastic.top_wall_deflection << "\n"
              << "  hold clock starts here; loading-period inelastic was omitted.\n";
    write_hold_status("elastic_balance", 0.0, 0, hold_steps, physical_dt, 0.0,
                      kinetic_energy(lower_particles) + kinetic_energy(upper_particles),
                      gap_elastic.max_penetration, cavity_elastic.center_net_height, inelastic_l0, inelastic_u0);

    int aborted = 0;
    for (int step = 0; step < hold_steps; ++step)
    {
        const HoldSnapshot saved_l = capture(lower_particles);
        const HoldSnapshot saved_u = capture(upper_particles);
        const Real pen_before = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
        const spdb::SlotCavityMetrics cav_before =
            spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                      upper_particles.TotalRealParticles(), dp);

        const Ren2014PhysicalCommit commit_l =
            commit_ren2014_physical_increment(lower_particles, lower_mat.ShearModulus(), demo_temperature, physical_dt);
        const Ren2014PhysicalCommit commit_u =
            commit_ren2014_physical_increment(upper_particles, upper_mat.ShearModulus(), demo_temperature, physical_dt);
        if (commit_l.rejected || commit_u.rejected)
        {
            restore(lower_particles, saved_l);
            restore(upper_particles, saved_u);
            std::cerr << "slot_hold physical commit rejected at hold_s=" << hold_physical << "\n";
            aborted = 1;
            break;
        }

        // Keep PhysicalTime at/above ramp so ForceInBodyRegion stays at full load.
        physical_time = ramp_time;
        if (relax_to_balance(relax_cap, "post") != 0)
        {
            restore(lower_particles, saved_l);
            restore(upper_particles, saved_u);
            aborted = 1;
            break;
        }
        hold_physical += physical_dt;
        ++hold_step_done;
        const Real pen_after = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp).max_penetration;
        const spdb::SlotCavityMetrics cav_after =
            spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                      upper_particles.TotalRealParticles(), dp);
        const Real inl = mean_inelastic(lower_particles);
        const Real inu = mean_inelastic(upper_particles);
        std::cout << "accepted hold_physical_s=" << hold_physical
                  << " mean_inelastic_l=" << inl
                  << " mean_inelastic_u=" << inu
                  << " max_penetration_m=" << pen_after
                  << " d_penetration_m=" << (pen_after - pen_before)
                  << " cavity_net_height_m=" << cav_after.center_net_height
                  << " d_cavity_m=" << (cav_after.center_net_height - cav_before.center_net_height)
                  << " top_wall_defl_m=" << cav_after.top_wall_deflection << "\n";
        write_hold_status("accepted", hold_physical, hold_step_done, hold_steps, physical_dt, 0.0,
                          kinetic_energy(lower_particles) + kinetic_energy(upper_particles), pen_after,
                          cav_after.center_net_height, inl, inu);
        write_states.writeToFile();
    }

    TickCount t2 = TickCount::now();
    const InterfaceGap gap_final = measure_interface_gap(lower_face, pos_l, upper_face, pos_u, dp);
    const spdb::SlotCavityMetrics cavity_final =
        spdb::measure_slot_cavity(pos_l, lower_particles.TotalRealParticles(), pos_u,
                                  upper_particles.TotalRealParticles(), dp);
    const Real applied_z = sum_field_z(upper_particles, "ForceInBodyRegion");
    const Real contact_l_z = sum_field_z(lower_particles, "RepulsionForce");
    const Real contact_u_z = sum_field_z(upper_particles, "RepulsionForce");
    Vecd *force_l = lower_particles.getVariableDataByName<Vecd>("Force");
    Vecd *prior_l = lower_particles.getVariableDataByName<Vecd>("ForcePrior");
    Real bottom_internal_z = 0.0;
    for (size_t i : bottom_part.body_part_particles_)
        bottom_internal_z += force_l[i][2] + prior_l[i][2];
    const Real reaction_z = -bottom_internal_z;
    const Real reaction_err = std::abs(std::abs(reaction_z) - force_total) / force_total;

    char metrics_name[192];
    std::snprintf(metrics_name, sizeof(metrics_name), "spdb_slot_hold_dp%.4g_dt%.4g_n%d_k%.4g_metrics.txt", dp,
                  physical_dt, hold_steps, stiffness_scale);
    std::cout << "slot_hold wall time: " << (t2 - t1).seconds() << " s\n"
              << "  hold_physical_s=" << hold_physical << " aborted=" << aborted << "\n"
              << "  max_penetration_elastic_m=" << gap_elastic.max_penetration
              << "  max_penetration_final_m=" << gap_final.max_penetration
              << "  d_penetration_hold_m=" << (gap_final.max_penetration - gap_elastic.max_penetration) << "\n"
              << "  cavity_net_height_elastic_m=" << cavity_elastic.center_net_height
              << "  cavity_net_height_final_m=" << cavity_final.center_net_height
              << "  d_cavity_hold_m=" << (cavity_final.center_net_height - cavity_elastic.center_net_height) << "\n"
              << "  mean_inelastic_l=" << mean_inelastic(lower_particles)
              << "  mean_inelastic_u=" << mean_inelastic(upper_particles)
              << "  (elastic_baseline_l=" << inelastic_l0 << " u=" << inelastic_u0 << ")\n"
              << "  reaction_rel_err=" << reaction_err << "\n"
              << "  metrics_file=" << metrics_name << "\n";

    std::ofstream out(metrics_name);
    out << std::setprecision(10);
    out << "case=test_3d_spdb_geometry\n";
    out << "mode=slot_hold\n";
    out << "units=SI\n";
    out << "material_source=demo_assumption\n";
    out << "material_note=NOT_950C_20Cr13_process_prediction\n";
    out << "constitutive=ContactScaledRen2014_Integration1stHalfPK2\n";
    out << "contact=ContactForce\n";
    out << "inelastic_during_relaxation=false\n";
    out << "physical_increment=commit_ren2014_physical_increment_once_per_accepted_step\n";
    out << "youngs_modulus_Pa=" << demo_youngs << "\n";
    out << "temperature_K=" << demo_temperature << "\n";
    out << "particle_spacing_m=" << dp << "\n";
    out << "contact_stiffness_scale=" << stiffness_scale << "\n";
    out << "contact_stiffness_Pa=" << lower_mat.ContactStiffness() << "\n";
    out << "physical_dt_s=" << physical_dt << "\n";
    out << "hold_steps_requested=" << hold_steps << "\n";
    out << "hold_physical_s=" << hold_physical << "\n";
    out << "hold_is_7200s=" << (std::abs(hold_physical - 7200.0) < 1.0e-6 ? 1 : 0) << "\n";
    out << "hold_is_not_7200s=" << (std::abs(hold_physical - 7200.0) < 1.0e-6 ? 0 : 1) << "\n";
    out << "aborted=" << aborted << "\n";
    write_hold_status(aborted ? "aborted" : "finished", hold_physical, hold_step_done, hold_steps, physical_dt, 0.0,
                      kinetic_energy(lower_particles) + kinetic_energy(upper_particles), gap_final.max_penetration,
                      cavity_final.center_net_height, mean_inelastic(lower_particles),
                      mean_inelastic(upper_particles), aborted, aborted ? 0 : 1);
    out << "cavity_volume_claimed=0\n";
    out << "initial_cavity_net_height_m=" << cavity0.center_net_height << "\n";
    out << "elastic_cavity_net_height_m=" << cavity_elastic.center_net_height << "\n";
    out << "final_cavity_net_height_m=" << cavity_final.center_net_height << "\n";
    out << "delta_cavity_elastic_m=" << (cavity_elastic.center_net_height - cavity0.center_net_height) << "\n";
    out << "delta_cavity_hold_m=" << (cavity_final.center_net_height - cavity_elastic.center_net_height) << "\n";
    out << "min_clearance_final_m=" << cavity_final.min_clearance << "\n";
    out << "cavity_width_x_m=" << cavity_final.width_x << "\n";
    out << "cavity_width_y_m=" << cavity_final.width_y << "\n";
    out << "top_wall_deflection_elastic_m=" << cavity_elastic.top_wall_deflection << "\n";
    out << "top_wall_deflection_final_m=" << cavity_final.top_wall_deflection << "\n";
    out << "max_penetration_elastic_m=" << gap_elastic.max_penetration << "\n";
    out << "max_penetration_final_m=" << gap_final.max_penetration << "\n";
    out << "delta_penetration_hold_m=" << (gap_final.max_penetration - gap_elastic.max_penetration) << "\n";
    out << "mean_inelastic_lower=" << mean_inelastic(lower_particles) << "\n";
    out << "mean_inelastic_upper=" << mean_inelastic(upper_particles) << "\n";
    out << "max_inelastic_lower=" << max_inelastic(lower_particles) << "\n";
    out << "max_inelastic_upper=" << max_inelastic(upper_particles) << "\n";
    out << "applied_force_z_N=" << applied_z << "\n";
    out << "contact_lower_z_N=" << contact_l_z << "\n";
    out << "contact_upper_z_N=" << contact_u_z << "\n";
    out << "reaction_z_N=" << reaction_z << "\n";
    out << "reaction_rel_err=" << reaction_err << "\n";
    out << "wall_time_s=" << (t2 - t1).seconds() << "\n";
    return aborted;
}
