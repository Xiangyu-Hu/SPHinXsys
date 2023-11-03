/**
 * @file 	channel_flow_shell.cpp
 * @brief 	This is a test of fluid-shell interaction.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h"

#include "case.h" //	case file to setup the test case
using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    std::cout << "System established..." << std::endl;
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    std::cout << "handleCommandlineOptions..." << std::endl;
    IOEnvironment io_environment(sph_system);
    std::cout << "io_environment..." << std::endl;
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    std::cout << "Body generation starts..." << std::endl;
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    std::cout << "Fluid generation finished..." << std::endl;

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
    std::cout << "Wall generation finished..." << std::endl;

    SolidBody valve(sph_system, makeShared<DefaultShape>("Valve"));
    valve.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_valve);
    valve.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    valve.generateParticles<ValveParticleGenerator>();
    valve.addBodyStateForRecording<Vecd>("PseudoNormal");
    std::cout << "Valve generation finished..." << std::endl;
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation valve_inner(valve);
    // Curvature calculation
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> valve_curvature(valve_inner);
    /** Must construct ShellContactRelation first! Otherwise ComplexRelation creates ContactRelation*/
    ShellContactRelation water_valve_contact(water_block, {&valve});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});

    ComplexRelation water_valve_complex(water_block_inner, water_valve_contact);
    ComplexRelation water_wall_complex(water_block_inner, water_wall_contact);

    ContactRelation valve_water_contact(valve, {&water_block});
    std::cout << "Contact relation finished..." << std::endl;
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** wall norm */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_valve_contact, water_wall_complex);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration1stHalfRiemann> pressure_relaxation(water_wall_contact, water_valve_complex);
    Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration2ndHalf> density_relaxation(water_wall_contact, water_valve_complex);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_valve_contact, water_wall_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithShellandWall> viscous_acceleration(water_wall_contact, water_valve_complex);
    /** Inflow boundary condition. */
    BodyAlignedBoxByCell inflow_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, water_block.getBodyShapeBounds(), xAxis);
    std::cout << "Fluid algs finished..." << std::endl;
    /** Algorithms for solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> valve_time_step_size(valve);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> valve_corrected_configuration(valve_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> valve_stress_relaxation_first(valve_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> valve_stress_relaxation_second(valve_inner);
    std::cout << "Solid algs finished..." << std::endl;
    /** FSI */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_valve(valve_water_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdate> fluid_force_on_valve_update(valve_water_contact, viscous_force_on_valve);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(valve);
    std::cout << "Fsi algs finished..." << std::endl;
    /** constraint and damping */
    BoundaryGeometry valve_boundary_geometry(valve, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> valve_constrain(valve_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        baffle_position_damping(0.2, valve_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        baffle_rotation_damping(0.2, valve_inner, "AngularVelocity", physical_viscosity);
    std::cout << "Constraint and damping finished..." << std::endl;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    valve.addBodyStateForRecording<Real>("MeanCurvature");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    std::cout << "Cell link finished..." << std::endl;
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    std::cout << "Periodic condition finished..." << std::endl;
    /** initialize configurations for all bodies. */
    // sph_system.initializeSystemConfigurations();
    water_block_inner.updateConfiguration();
    std::cout << "water inner finished..." << std::endl;
    valve_inner.updateConfiguration();
    std::cout << "valve inner finished..." << std::endl;
    water_wall_contact.updateConfiguration();
    std::cout << "water wall contact finished..." << std::endl;
    valve_water_contact.updateConfiguration();
    std::cout << "valve water contact finished..." << std::endl;
    water_valve_contact.updateConfiguration();
    std::cout << "water valve contact finished..." << std::endl;
    /**wall normal*/
    wall_boundary_normal_direction.exec();
    std::cout << "Wall normal finished..." << std::endl;
    /** initial curvature*/
    valve_corrected_configuration.exec();
    valve_curvature.compute_initial_curvature();
    water_valve_contact.updateConfiguration();
    std::cout << "Curvature finished..." << std::endl;

    // Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness(water_block_inner, water_valve_contact);
    check_kernel_completeness.exec();
    water_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_block.addBodyStateForRecording<Real>("TotalKernel");
    water_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.5;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_correction.exec();

            viscous_force_on_valve.exec();

            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /**FSI for pressure force*/
                fluid_force_on_valve_update.exec();
                /** velocity */
                parabolic_inflow.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics time stepping. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(valve_time_step_size.exec(), dt - dt_s_sum);
                    valve_stress_relaxation_first.exec(dt_s);
                    valve_constrain.exec();
                    valve_stress_relaxation_second.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt
                          << "	dt / dt_s = " << inner_ite_dt_s
                          << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            // update curvature before water shell contact update
            valve_curvature.exec();
            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            valve.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();

            water_block_inner.updateConfiguration();
            water_wall_contact.updateConfiguration();
            water_valve_contact.updateConfiguration();

            valve_water_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        check_kernel_completeness.exec();
        /** write run-time observation into file */
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
