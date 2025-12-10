/**
 * @file 	hydrostatic_fsi.cpp
 * @brief 	shell deformation due to hydrostatic pressure under gravity.
 * @details This is the one of the basic test cases
 * for understanding SPH method for fluid-shell-interaction (FSI) simulation.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

namespace SPH
{
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::vector<Vecd> &shape, const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	wall body shape definition.
//----------------------------------------------------------------------
class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real DH;
    Real DL;
    Real particle_spacing_gate;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real DH, Real DL, Real particle_spacing_gate)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          DH(DH), DL(DL), particle_spacing_gate(particle_spacing_gate) {};
    void prepareGeometricData() override
    {
        const auto particle_number_wall = int(DH / particle_spacing_gate);
        for (int i = 0; i < particle_number_wall; i++)
        {
            Real x1 = -0.5 * particle_spacing_gate;
            Real x2 = DL + 0.5 * particle_spacing_gate;
            Real y = 0.5 * particle_spacing_gate + Real(i) * particle_spacing_gate;
            const Vec2d normal_direction_1(1.0, 0.0);
            const Vec2d normal_direction_2(-1.0, 0.0);
            addPositionAndVolumetricMeasure(Vecd(x1, y), particle_spacing_gate);
            addSurfaceProperties(normal_direction_1, particle_spacing_gate);
            addPositionAndVolumetricMeasure(Vecd(x2, y), particle_spacing_gate);
            addSurfaceProperties(normal_direction_2, particle_spacing_gate);
        }
    }
};
//----------------------------------------------------------------------
//	gate body shape definition.
//----------------------------------------------------------------------
class Gate;
template <>
class ParticleGenerator<SurfaceParticles, Gate> : public ParticleGenerator<SurfaceParticles>
{
    Real DL;
    Real BW;
    Real particle_spacing_gate;
    Real Gate_thickness;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real DL, Real BW, Real particle_spacing_gate, Real Gate_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          DL(DL), BW(BW), particle_spacing_gate(particle_spacing_gate), Gate_thickness(Gate_thickness) {};
    void prepareGeometricData() override
    {
        const auto particle_number_gate = int((DL + 2 * BW) / particle_spacing_gate);
        // generate particles for the elastic gate
        for (int i = 0; i < particle_number_gate; i++)
        {
            Real x = -BW + 0.5 * particle_spacing_gate + Real(i) * particle_spacing_gate;
            Real y = -0.5 * particle_spacing_gate;
            addPositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_gate);
            const Vec2d normal_direction(0, 1.0);
            addSurfaceProperties(normal_direction, Gate_thickness);
        }
    }
};
} // namespace SPH

void hydrostatic_fsi(const Real particle_spacing_gate, const Real particle_spacing_ref)
{
    //----------------------------------------------------------------------
    //	Basic geometry parameters and numerical setup.
    //----------------------------------------------------------------------
    const Real DL = 1.0;              /**< Tank length. */
    const Real DH = 2.1;              /**< Tank height. */
    const Real Dam_L = 1.0;           /**< Water block width. */
    const Real Dam_H = 2.0;           /**< Water block height. */
    const Real Gate_thickness = 0.05; /**< Width of the gate. */
    const Real BW = particle_spacing_ref * 4.0;
    const BoundingBoxd system_domain_bounds(Vec2d(-BW, -std::max(particle_spacing_gate, Gate_thickness)), Vec2d(DL + BW, DH + Gate_thickness));
    // observer location
    const StdVec<Vecd> observation_location = {Vecd(0.5 * Dam_L, -0.5 * particle_spacing_gate)};
    //----------------------------------------------------------------------
    //	Define the corner point of water block geometry.
    //----------------------------------------------------------------------
    const Vec2d DamP_lb(0.0, 0.0);     /**< Left bottom. */
    const Vec2d DamP_lt(0.0, Dam_H);   /**< Left top. */
    const Vec2d DamP_rt(Dam_L, Dam_H); /**< Right top. */
    const Vec2d DamP_rb(Dam_L, 0.0);   /**< Right bottom. */
    //----------------------------------------------------------------------
    //	Geometry definition.
    //----------------------------------------------------------------------
    auto createWaterBlockShape = [&]()
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(DamP_lb);
        water_block_shape.push_back(DamP_lt);
        water_block_shape.push_back(DamP_rt);
        water_block_shape.push_back(DamP_rb);
        water_block_shape.push_back(DamP_lb);

        return water_block_shape;
    };
    //----------------------------------------------------------------------
    //	Define the geometry for gate constrain.
    //----------------------------------------------------------------------
    Vec2d ConstrainLP_lb(-BW, -particle_spacing_ref);
    Vec2d ConstrainLP_lt(-BW, 0.0);
    Vec2d ConstrainLP_rt(0.0, 0.0);
    Vec2d ConstrainLP_rb(0.0, -particle_spacing_ref);
    Vec2d ConstrainRP_lb(Dam_L, -particle_spacing_ref);
    Vec2d ConstrainRP_lt(Dam_L, 0.0);
    Vec2d ConstrainRP_rt(Dam_L + BW, 0.0);
    Vec2d ConstrainRP_rb(Dam_L + BW, -particle_spacing_ref);
    //----------------------------------------------------------------------
    //	create Gate constrain shape
    //----------------------------------------------------------------------
    auto createGateConstrainShape = [&]()
    {
        // geometry
        std::vector<Vecd> gate_constraint_shape_left;
        gate_constraint_shape_left.push_back(ConstrainLP_lb);
        gate_constraint_shape_left.push_back(ConstrainLP_lt);
        gate_constraint_shape_left.push_back(ConstrainLP_rt);
        gate_constraint_shape_left.push_back(ConstrainLP_rb);
        gate_constraint_shape_left.push_back(ConstrainLP_lb);

        std::vector<Vecd> gate_constraint_shape_right;
        gate_constraint_shape_right.push_back(ConstrainRP_lb);
        gate_constraint_shape_right.push_back(ConstrainRP_lt);
        gate_constraint_shape_right.push_back(ConstrainRP_rt);
        gate_constraint_shape_right.push_back(ConstrainRP_rb);
        gate_constraint_shape_right.push_back(ConstrainRP_lb);

        MultiPolygon multi_polygon;
        multi_polygon.addAPolygon(gate_constraint_shape_left, ShapeBooleanOps::add);
        multi_polygon.addAPolygon(gate_constraint_shape_right, ShapeBooleanOps::add);
        return multi_polygon;
    };
    //----------------------------------------------------------------------
    //	Material properties of the fluid.
    //----------------------------------------------------------------------
    const Real rho0_f = 1000.0;                       /**< Reference density of fluid. */
    const Real gravity_g = 9.81;                      /**< Value of gravity. */
    const Real U_ref = 2.0 * sqrt(Dam_H * gravity_g); /**< Characteristic velocity. */
    const Real c_f = 10.0 * U_ref;                    /**< Reference sound speed. */
    const Real Re = 0.1;                              /**< Reynolds number. */
    const Real mu_f = rho0_f * U_ref * DL / Re;       /**< Dynamics viscosity. */
    //----------------------------------------------------------------------
    //	Material properties of the elastic gate.
    //----------------------------------------------------------------------
    const Real rho0_s = 2700.0; /**< Reference solid density. */
    const Real poisson = 0.495; /**< Poisson ratio. */
    const Real Ae = 6.75e10;    /**< Normalized Youngs Modulus. */
    const Real Youngs_modulus = Ae;
    const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * Gate_thickness * Gate_thickness;
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>(createWaterBlockShape(), "WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_gate);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<SurfaceParticles, WallBoundary>(DH, DL, particle_spacing_gate);

    SolidBody gate(sph_system, makeShared<DefaultShape>("Gate"));
    gate.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_gate);
    gate.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<SurfaceParticles, Gate>(DL, BW, particle_spacing_gate, Gate_thickness);
    //----------------------------------------------------------------------
    //	Particle and body creation of gate observer.
    //----------------------------------------------------------------------
    ObserverBody gate_observer(sph_system, "Observer");
    gate_observer.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_gate);
    gate_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation gate_inner(gate);
    // shell normal should point from fluid to shell
    // normal corrector set to true if shell normal is currently pointing from shell to fluid
    ContactRelationFromShellToFluid water_block_contact(water_block, {&wall_boundary, &gate}, {true, true});
    ContactRelationFromFluidToShell gate_contact(gate, {&water_block}, {true});
    ContactRelation gate_observer_contact(gate_observer, {&gate});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(gate, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> gate_corrected_configuration(gate_inner);

    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> gate_stress_relaxation_first_half(gate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> gate_stress_relaxation_second_half(gate_inner);

    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> gate_computing_time_step_size(gate);
    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> gate_constraint(gate_constraint_part);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> gate_update_normal(gate);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> gate_curvature(shell_curvature_inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> gate_position_damping(0.2, gate_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> gate_rotation_damping(0.2, gate_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define fluid methods which are used in this case.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_fluid_density(water_block_inner, water_block_contact);

    /** Compute time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(water_block, U_ref);
    /** Compute time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec2d, FixedDampingRate>>>
        fluid_damping(0.2, DynamicsArgs(water_block_inner, "Velocity", mu_f), DynamicsArgs(water_block_contact, "Velocity", mu_f));
    //----------------------------------------------------------------------
    //	Define fsi methods which are used in this case.
    //----------------------------------------------------------------------
    /** Compute the average velocity of gate. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
    /** Compute the force exerted on elastic gate due to fluid pressure. */
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> fluid_pressure_force_on_gate(gate_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system);
    write_real_body_states_to_vtp.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states_to_vtp.addToWrite<Real>(gate, "Average1stPrincipleCurvature");
    write_real_body_states_to_vtp.addToWrite<Real>(gate, "Average2ndPrincipleCurvature");
    write_real_body_states_to_vtp.addToWrite<Vecd>(gate, "PressureForceFromFluid");
    write_real_body_states_to_vtp.addDerivedVariableRecording<SimpleDynamics<Displacement>>(gate);
    /** Output the observed displacement of gate center. */
    ObservedQuantityRecording<Vecd> write_beam_tip_displacement("Displacement", gate_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing linear reproducing configuration for the insert body. */
    gate_corrected_configuration.exec();
    gate_curvature.exec();
    /** update fluid-shell contact*/
    water_block_contact.updateConfiguration();
    gate_contact.updateConfiguration();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.2; /**< End time. */
    Real output_interval = end_time / 100.0;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                fluid_damping.exec(dt);
                /** Fluid relaxation and force computation. */
                pressure_relaxation.exec(dt);
                fluid_pressure_force_on_gate.exec();
                density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    gate_stress_relaxation_first_half.exec(dt_s);
                    gate_constraint.exec();
                    gate_position_damping.exec(dt_s);
                    gate_rotation_damping.exec(dt_s);
                    gate_constraint.exec();
                    gate_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    dt_s = gate_computing_time_step_size.exec();
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update normal direction on elastic body. */
            gate_update_normal.exec();
            /** Update curvature. */
            gate.updateCellLinkedList();
            shell_curvature_inner.updateConfiguration();
            gate_curvature.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList(); // water particle motion is small

            /** one need update configuration after periodic condition. */
            water_block_complex.updateConfiguration();
            gate_contact.updateConfiguration();

            /** Output the observed data. */
            write_beam_tip_displacement.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        write_real_body_states_to_vtp.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    // analytical solution
    const Real p = rho0_f * gravity_g * Dam_H;
    const Real I = 1 / 12.0 * 1.0 * std::pow(Gate_thickness, 3);
    const Real max_disp_analytical = p * std::pow(Dam_L, 4) / 384.0 / Youngs_modulus / I;

    // Compare with analytical solution
    const Real max_disp = std::abs(write_beam_tip_displacement.getObservedQuantity()[0][1]);
    const Real error = std::abs((max_disp_analytical - max_disp) / max_disp_analytical) * 100.0;

    std::cout << "Analytical displacement: " << max_disp_analytical
              << "\t Displacement: " << max_disp
              << "\t Error: " << error << "%"
              << std::endl;

    // gtest
    EXPECT_NEAR(max_disp_analytical, max_disp, max_disp_analytical * 15e-2);
}

TEST(hydrostatic_fsi, dp_2)
{ // for CI
    const Real Gate_thickness = 0.05;
    const Real particle_spacing_gate = Gate_thickness / 2.0;
    const Real particle_spacing_ref = particle_spacing_gate;
    hydrostatic_fsi(particle_spacing_gate, particle_spacing_ref);
}

TEST(DISABLED_hydrostatic_fsi, dp_4)
{ // for CI
    const Real Gate_thickness = 0.05;
    const Real particle_spacing_gate = Gate_thickness / 4.0;
    const Real particle_spacing_ref = particle_spacing_gate;
    hydrostatic_fsi(particle_spacing_gate, particle_spacing_ref);
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
