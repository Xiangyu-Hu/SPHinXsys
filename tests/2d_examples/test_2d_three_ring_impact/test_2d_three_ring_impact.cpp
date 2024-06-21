/**
 * @file 	Three ring impact.cpp
 * @brief 	This is the case file for the test of dynamic contacts between shell and solid.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

void relax_solid(RealBody &body, InnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 200 == 0)
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
}

class Ring : public MultiPolygonShape
{
  public:
    explicit Ring(const std::string &shape_name, const Vec2d &center, Real radius_inner, Real radius_outer) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, radius_outer, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(center, radius_inner, 100, ShapeBooleanOps::sub);
    }
};

namespace SPH
{
class ShellRing;
template <>
class ParticleGenerator<ShellRing> : public ParticleGenerator<Surface>
{
    const Vec2d center_;
    const Real mid_srf_radius_;
    const Real dp_;
    const Real thickness_;

  public:
    ParticleGenerator(SPHBody &sph_body, const Vec2d &center, Real mid_srf_radius, Real dp, Real thickness)
        : ParticleGenerator<Surface>(sph_body),
          center_(center),
          mid_srf_radius_(mid_srf_radius),
          dp_(dp),
          thickness_(thickness){};
    void initializeGeometricVariables() override
    {
        const auto number_of_particles = int(2 * Pi * mid_srf_radius_ / dp_);
        const Real dtheta = 2 * Pi / Real(number_of_particles);
        for (int n = 0; n < number_of_particles; n++)
        {
            const Vec2d center_to_pos = mid_srf_radius_ * Vec2d(cos(n * dtheta), sin(n * dtheta));
            initializePositionAndVolumetricMeasure(center_ + center_to_pos, dp_);
            initializeSurfaceProperties(center_to_pos.normalized(), thickness_);
        }
    }
};
} // namespace SPH

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

void three_ring_impact(int resolution_factor_l, int resolution_factor_m, int resolution_factor_s)
{
    // end time
    const Real end_time = 1.0;

    // main geometric parameters
    // inner and outer diameter of the large, medium and small ring
    const Real diameter_inner_l = 26;
    const Real diameter_outer_l = 30;
    const Real diameter_inner_m = 10;
    const Real diameter_outer_m = 12;
    const Real diameter_inner_s = 8;
    const Real diameter_outer_s = 10;

    // mid surface diameter and thickness
    const Real mid_srf_diameter_m = 0.5 * (diameter_inner_m + diameter_outer_m);
    const Real mid_srf_diameter_s = 0.5 * (diameter_inner_s + diameter_outer_s);
    const Real thickness_l = 0.5 * (diameter_outer_l - diameter_inner_l);
    const Real thickness_m = 0.5 * (diameter_outer_m - diameter_inner_m);
    const Real thickness_s = 0.5 * (diameter_outer_s - diameter_inner_s);

    // center position
    const Vec2d center_l(0, 0);
    const Vec2d center_m = 0.5 * Vec2d(-7.9, 7.9);
    const Vec2d center_s = 0.5 * Vec2d(7.9, -8.5);

    // resolution
    const Real dp_ref = thickness_l / 4.0; // 4 particles in thickness direction
    const Real dp_l = dp_ref / Real(resolution_factor_l);
    const Real dp_m = dp_ref / Real(resolution_factor_m);
    const Real dp_s = dp_ref / Real(resolution_factor_s);

    // material
    const Real rho_l = 1.0;
    const Real rho_m = 0.01;
    const Real rho_s = 0.1;

    const Real youngs_modulus_l = 288e3;
    const Real youngs_modulus_m = 2250;
    const Real youngs_modulus_s = 10e3;

    const Real possion_ratio = 0.125;

    const Real physical_viscosity_l = get_physical_viscosity_general(rho_l, youngs_modulus_l, thickness_l);
    const Real physical_viscosity_m = get_physical_viscosity_general(rho_m, youngs_modulus_m, thickness_m);
    const Real physical_viscosity_s = get_physical_viscosity_general(rho_s, youngs_modulus_s, thickness_s);

    auto material_l = makeShared<NeoHookeanSolid>(rho_l, youngs_modulus_l, possion_ratio);
    auto material_m = makeShared<NeoHookeanSolid>(rho_m, youngs_modulus_m, possion_ratio);
    auto material_s = makeShared<NeoHookeanSolid>(rho_s, youngs_modulus_s, possion_ratio);

    // Bounding box
    BoundingBox bb_system(center_l - 0.5 * diameter_outer_l * Vec2d::Ones(),
                          center_l + 0.5 * diameter_outer_l * Vec2d::Ones());

    // System
    SPHSystem system(bb_system, dp_l);
    IOEnvironment io_environment(system);

    // Body
    SolidBody ring_l_body(system, makeShared<Ring>("RingLarge", center_l, 0.5 * diameter_inner_l, 0.5 * diameter_outer_l));
    ring_l_body.defineBodyLevelSetShape();
    ring_l_body.assignMaterial(material_l.get());
    ring_l_body.generateParticles<BaseParticles, Lattice>();
    auto particles_l = &ring_l_body.getBaseParticles();

    SolidBody ring_m_body(system, makeShared<DefaultShape>("RingMedium"));
    ring_m_body.defineAdaptationRatios(1.15, dp_l / dp_m);
    ring_m_body.assignMaterial(material_m.get());
    ring_m_body.generateParticles<SurfaceParticles, ShellRing>(center_m, 0.5 * mid_srf_diameter_m, dp_m, thickness_m);
    auto particles_m = &ring_m_body.getBaseParticles();

    SolidBody ring_s_body(system, makeShared<DefaultShape>("RingSmall"));
    ring_s_body.defineAdaptationRatios(1.15, dp_l / dp_s);
    ring_s_body.assignMaterial(material_s.get());
    ring_s_body.generateParticles<SurfaceParticles, ShellRing>(center_s, 0.5 * mid_srf_diameter_s, dp_s, thickness_s);
    auto particles_s = &ring_s_body.getBaseParticles();

    // Inner relation
    InnerRelation ring_l_inner(ring_l_body);
    InnerRelation ring_m_inner(ring_m_body);
    InnerRelation ring_s_inner(ring_s_body);

    // relaxation
    relax_solid(ring_l_body, ring_l_inner);

    // Methods
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_l(ring_l_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_l(ring_l_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_l(ring_l_inner);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec2d>>>
        velocity_damping_l(0.2, ring_l_inner, "Velocity", physical_viscosity_l);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size_l(ring_l_body);

    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration_m(ring_m_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half_m(ring_m_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half_m(ring_m_inner);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec2d>>>
        velocity_damping_m(0.2, ring_m_inner, "Velocity", physical_viscosity_m);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec2d>>>
        rotation_damping_m(0.2, ring_m_inner, "AngularVelocity", physical_viscosity_m);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_m(ring_m_body);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size_m(ring_m_body);

    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration_s(ring_s_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half_s(ring_s_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half_s(ring_s_inner);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec2d>>>
        velocity_damping_s(0.2, ring_s_inner, "Velocity", physical_viscosity_s);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec2d>>>
        rotation_damping_s(0.2, ring_s_inner, "AngularVelocity", physical_viscosity_s);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_s(ring_s_body);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size_s(ring_s_body);

    // self contact
    ShellSelfContactRelation self_contact_m(ring_m_body);

    // Contact relation
    // shell-shell
    SurfaceContactRelationFromShellToShell contact_s_to_m(ring_s_body, {&ring_m_body}, {true});
    SurfaceContactRelationFromShellToShell contact_m_to_s(ring_m_body, {&ring_s_body}, {true});
    // shell-solid
    SurfaceContactRelationToShell contact_m_to_l(ring_m_body, {&ring_l_body});
    // solid-shell
    SurfaceContactRelationFromShell contact_l_to_m(ring_l_body, {&ring_m_body}, {true});

    // Inner relation of curvature
    ShellInnerRelationWithContactKernel curvature_inner_m_with_s_kernel(ring_m_body, ring_s_body);
    ShellInnerRelationWithContactKernel curvature_inner_s_with_m_kernel(ring_s_body, ring_m_body);

    // Contact method
    // curvature update
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> curvature_m(ring_m_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> average_curvature_m_with_s_kernel(curvature_inner_m_with_s_kernel);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> average_curvature_s_with_m_kernel(curvature_inner_s_with_m_kernel);

    // density
    // self contact
    InteractionDynamics<solid_dynamics::ShellSelfContactDensitySummation> self_contact_density_m(self_contact_m);
    // shell-shell
    InteractionDynamics<solid_dynamics::ContactDensitySummationFromShell> contact_density_s_to_m(contact_s_to_m);
    InteractionDynamics<solid_dynamics::ContactDensitySummationFromShell> contact_density_m_to_s(contact_m_to_s);
    // shell-solid
    InteractionDynamics<solid_dynamics::ContactDensitySummation> contact_density_m_to_l(contact_m_to_l);
    // solid-shell
    InteractionDynamics<solid_dynamics::ContactDensitySummationFromShell> contact_density_l_to_m(contact_l_to_m);

    // force
    // self contact
    InteractionWithUpdate<solid_dynamics::SelfContactForce> self_contact_forces_m(self_contact_m);
    // shell-shell
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_forces_s_to_m(contact_s_to_m);
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_forces_m_to_s(contact_m_to_s);
    // shell-solid
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_m_to_l(contact_m_to_l);
    // solid-shell
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_l_to_m(contact_l_to_m);

    // Inital condition
    auto vel_ic = [&, vel = Vec2d(-30, 30)]()
    {
        auto &vel_ = *particles_s->getVariableByName<Vec2d>("Velocity");
        particle_for(
            execution::ParallelPolicy(),
            IndexRange(0, particles_s->total_real_particles_),
            [&](size_t index_i)
            {
                vel_[index_i] = vel;
            });
    };

    // Boundary condition
    auto fix_id_l = [&]()
    {
        IndexVector fixed_particles_id;
        for (size_t i = 0; i < particles_l->total_real_particles_; ++i)
        {
            Real radius = (particles_l->ParticlePositions()[i] - center_l).norm();
            if (radius > 0.5 * diameter_outer_l - 0.7 * dp_l)
                fixed_particles_id.push_back(i);
        }
        return fixed_particles_id;
    }();
    auto fix_bc = [&]()
    {
        auto &vel_ = *particles_l->getVariableByName<Vec2d>("Velocity");
        particle_for(
            execution::ParallelPolicy(),
            fix_id_l,
            [&](size_t index_i)
            {
                vel_[index_i] = Vec2d::Zero();
            });
    };

    // Check
    auto check_nan = [&](BaseParticles *particles)
    {
        const auto &pos_ = particles->ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    corrected_configuration_l.exec();
    corrected_configuration_m.exec();
    corrected_configuration_s.exec();

    curvature_m.compute_initial_curvature();
    average_curvature_m_with_s_kernel.exec();
    average_curvature_s_with_m_kernel.exec();
    contact_s_to_m.updateConfiguration();
    contact_m_to_s.updateConfiguration();

    vel_ic();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min({computing_time_step_size_l.exec(), computing_time_step_size_m.exec(), computing_time_step_size_s.exec()});
    auto run_simulation = [&]()
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                contact_density_s_to_m.exec();
                contact_forces_s_to_m.exec();

                contact_density_m_to_s.exec();
                contact_forces_m_to_s.exec();

                contact_density_l_to_m.exec();
                contact_force_l_to_m.exec();

                contact_density_m_to_l.exec();
                contact_force_m_to_l.exec();

                self_contact_density_m.exec();
                self_contact_forces_m.exec();

                dt = std::min({computing_time_step_size_l.exec(), computing_time_step_size_m.exec(), computing_time_step_size_s.exec()});
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                stress_relaxation_first_half_l.exec(dt);
                stress_relaxation_first_half_m.exec(dt);
                stress_relaxation_first_half_s.exec(dt);

                fix_bc();

                velocity_damping_l.exec(dt);
                velocity_damping_m.exec(dt);
                rotation_damping_m.exec(dt);
                rotation_damping_s.exec(dt);
                velocity_damping_s.exec(dt);

                fix_bc();

                stress_relaxation_second_half_l.exec(dt);
                stress_relaxation_second_half_m.exec(dt);
                stress_relaxation_second_half_s.exec(dt);

                ring_l_body.updateCellLinkedList();
                ring_m_body.updateCellLinkedList();
                ring_s_body.updateCellLinkedList();
                update_normal_m.exec();
                update_normal_s.exec();
                curvature_inner_m_with_s_kernel.updateConfiguration();
                curvature_inner_s_with_m_kernel.updateConfiguration();
                average_curvature_m_with_s_kernel.exec();
                average_curvature_s_with_m_kernel.exec();
                curvature_m.exec();
                self_contact_m.updateConfiguration();
                contact_s_to_m.updateConfiguration();
                contact_m_to_s.updateConfiguration();
                contact_m_to_l.updateConfiguration();
                contact_l_to_m.updateConfiguration();

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                { // checking if any position has become nan
                    check_nan(particles_l);
                    check_nan(particles_m);
                    check_nan(particles_s);
                }
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        vtp_output.writeToFile(ite_output + 1);
        exit(0);
    }
}
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    three_ring_impact(2, 2, 2);
}
