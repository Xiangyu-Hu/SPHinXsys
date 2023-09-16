/**
 * @file 	stfb.cpp
 * @brief 	This is the case file for 2D still floaing body.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
#include "stfb.h" //header for this case

int main(int ac, char *av[])
{
    std::cout << "Mass " << StructureMass << " str_sup " << FlStA << " rho_s " << rho_s << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<Real>("VolumetricMeasure");

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody structure(system, makeShared<TransformShape<GeometricShapeBox>>(
                                    Transform(structure_translation), structure_halfsize, "Structure"));
    structure.defineParticlesAndMaterial<SolidParticles, Solid>(rho_s);
    structure.generateParticles<ParticleGeneratorLattice>();

    ObserverBody observer(system, "Observer");
    observer.defineAdaptationRatios(1.15, 2.0);
    observer.generateParticles<ObserverParticleGenerator>(
        StdVec<Vecd>{obs});
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation structure_inner(structure);
    ComplexRelation water_block_complex(water_block_inner, {&wall_boundary, &structure});
    ContactRelation structure_contact(structure, {&water_block});
    ContactRelation observer_contact_with_water(observer, {&water_block});
    ContactRelation observer_contact_with_structure(observer, {&structure});
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<OffsetInitialPosition> structure_offset_position(structure, offset);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> str_normal(structure);
    /** corrected strong configuration. */
    InteractionWithUpdate<CorrectedConfigurationInner> str_corrected_conf(structure_inner);
    /** Time step initialization, add gravity. */
    SimpleDynamics<TimeStepInitialization> initialize_time_step_to_fluid(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    /** time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** pressure relaxation using Verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Fluid force on structure. */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(structure_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_solid(structure_contact, viscous_force_on_solid);
    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    /** set up the multi body system. */
    SimTK::MultibodySystem MBsystem;
    /** the bodies or matter of the system. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** the forces of the system. */
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** mass properties of the fixed spot. */
    StructureSystemForSimbody structure_multibody(structure, makeShared<TransformShape<GeometricShapeBox>>(
                                                                 Transform(structure_translation), structure_halfsize, "Structure"));
    /** Mass properties of the constrained spot.
     * SimTK::MassProperties(mass, center of mass, inertia)
     */
    SimTK::Body::Rigid structure_info(*structure_multibody.body_part_mass_properties_);
    /**
     * @brief  ** Create a %Planar mobilizer between an existing parent (inboard) body P
     *	and a new child (outboard) body B created by copying the given \a bodyInfo
     *	into a privately-owned Body within the constructed %MobilizedBody object.
     *	Specify the mobilizer frames F fixed to parent P and M fixed to child B.
     * @param[in] inboard(SimTKVec3) Defines the location of the joint point relative to the parent body.
     * @param[in] outboard(SimTKVec3) Defines the body's origin location to the joint point.
     * @note	The body's origin location can be the mass center, the the center of mass should be SimTKVec3(0)
     * 			in SimTK::MassProperties(mass, com, inertia)
     */
    SimTK::MobilizedBody::Planar structure_mob(matter.Ground(), SimTK::Transform(SimTKVec3(G[0], G[1], 0.0)), structure_info, SimTK::Transform(SimTKVec3(0.0, 0.0, 0.0)));
    /**
     * @details Add gravity to mb body.
     * @param[in,out] forces, The subsystem to which this force should be added.
     * @param[in]     matter, The subsystem containing the bodies that will be affected.
     * @param[in]    gravity, The default gravity vector v, interpreted as v=g*d where g=|\a gravity| is
     *				a positive scalar and d is the "down" direction unit vector d=\a gravity/g.
     * @param[in]  zeroHeight This is an optional specification of the default value for the height
     *				up the gravity vector that is considered to be "zero" for purposes of
     *				calculating the gravitational potential energy. The default is
     *				\a zeroHeight == 0, i.e., a body's potential energy is defined to be zero
     *				when the height of its mass center is the same as the height of the Ground
     *				origin. The zero height will have the value specified here unless
     *				explicitly changed within a particular State use the setZeroHeight()
     *			method.
     * @par Force Each body B that has not been explicitly excluded will experience a force
     *		fb = mb*g*d, applied to its center of mass, where mb is the mass of body B.
     * @par Potential Energy
     *		Gravitational potential energy for a body B is mb*g*hb where hb is the height of
     *		body B's mass center over an arbitrary "zero" height hz (default is hz=0),
     *		measured along the "up" direction -d. If pb is the Ground frame vector giving
     *		the position of body B's mass center, its height over or under hz is
     *		hb=pb*(-d) - hz. Note that this is a signed quantity so the potential energy is
     *		also signed. 0.475
     */
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(0.0, -gravity_g, 0.0), 0.0);
    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_structure(structure_multibody, MBsystem, structure_mob, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_on_structure(structure_multibody, MBsystem, structure_mob, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    BodyRegionByCell wave_probe_buffer(water_block, makeShared<TransformShape<GeometricShapeBox>>(
                                                        Transform(gauge_translation), gauge_halfsize, "FreeSurfaceGauge"));
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>> wave_gauge(io_environment, wave_probe_buffer);
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_observer_position(observer_contact_with_structure, "Position", "Position");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_str_displacement("Position", io_environment, observer_contact_with_structure);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    structure_offset_position.exec();
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    str_normal.exec();
    str_corrected_conf.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    write_str_displacement.writeToFile(0);
    wave_gauge.writeToFile(0);
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    int number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = total_physical_time;
    Real output_interval = end_time / 100;
    Real dt = 0.0;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_interval)
        {
            initialize_time_step_to_fluid.exec();

            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            /** Viscous force exerting on structure. */
            viscous_force_on_solid.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = get_fluid_time_step_size.exec();

                pressure_relaxation.exec(dt);
                fluid_force_on_solid.exec();
                density_relaxation.exec(dt);
                /** coupled rigid body dynamics. */
                if (total_time >= relax_time)
                {
                    SimTK::State &state_for_update = integ.updAdvancedState();
                    force_on_bodies.clearAllBodyForces(state_for_update);
                    force_on_bodies.setOneBodyForce(state_for_update, structure_mob, force_on_structure.exec());
                    integ.stepBy(dt);
                    constraint_on_structure.exec();
                }
                interpolation_observer_position.exec();

                relaxation_time += dt;
                integral_time += dt;
                total_time += dt;
                if (total_time >= relax_time)
                    GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << "	Total Time = " << total_time
                          << "	Physical Time = " << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            water_block.updateCellLinkedListWithParticleSort(100);
            wall_boundary.updateCellLinkedList();
            structure.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            structure_contact.updateConfiguration();
            observer_contact_with_water.updateConfiguration();

            if (total_time >= relax_time)
            {
                write_str_displacement.writeToFile(number_of_iterations);
                wave_gauge.writeToFile(number_of_iterations);
            }
        }

        TickCount t2 = TickCount::now();
        if (total_time >= relax_time)
            write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_str_displacement.generateDataBase(0.001);
        wave_gauge.generateDataBase(0.001);
    }
    else
    {
        write_str_displacement.testResult();
        wave_gauge.testResult();
    }

    return 0;
}
