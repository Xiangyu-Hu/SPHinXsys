/**
 * @file 	nonlinear_wave_fsi.h
 * @brief 	This is the 2d case file for wave impact with tension leg moored floating structure.
 * @author  Nicolò Salis
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
#include "nonlinear_wave_fsi.h" //header for this case

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody structure(sph_system, makeShared<FloatingStructure>("Structure"));
    structure.defineMaterial<Solid>(rho_s);
    structure.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer(sph_system, "Observer");
    observer.defineAdaptationRatios(1.15, 2.0);
    observer.generateParticles<ObserverParticles>(StdVec<Vecd>{obs});
    //---------------------------------------------------------
    // PRESSURE PROBES
    //---------------------------------------------------------
    ObserverBody fp2(sph_system, "FluidObserver2");
    Real fp2x = 12.466;
    Real fp2y = 0.968;
    StdVec<Vecd> fp2l = {Vecd(fp2x, fp2y)};
    fp2.defineAdaptationRatios(1.15, 2.0);
    fp2.generateParticles<ObserverParticles>(fp2l);
    ObserverBody fp3(sph_system, "FluidObserver3");
    Real fp3x = 12.466;
    Real fp3y = 1.013;
    StdVec<Vecd> fp3l = {Vecd(fp3x, fp3y)};
    fp3.defineAdaptationRatios(1.15, 2.0);
    fp3.generateParticles<ObserverParticles>(fp3l);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation structure_inner(structure);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &structure});
    ContactRelation structure_contact(structure, {&water_block});
    ContactRelation observer_contact_with_water(observer, {&water_block});
    ContactRelation observer_contact_with_structure(observer, {&structure});
    ContactRelation fp2_contact_s(fp2, {&structure});
    ContactRelation fp3_contact_s(fp3, {&structure});
    ContactRelation fp2_contact_w(fp2, {&water_block});
    ContactRelation fp3_contact_w(fp3, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    BodyRegionByParticle wave_maker(wall_boundary, makeShared<MultiPolygonShape>(createWaveMakerShape()));
    SimpleDynamics<WaveMaking> wave_making(wave_maker);

    SimpleDynamics<OffsetInitialPosition> structure_offset_position(structure, offset);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> structure_normal_direction(structure);
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(ConstructorArgs(water_block_inner, 0.1), water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Fluid force on structure. */
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(structure_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> fluid_force_on_structure(structure_contact);
    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    /** set up the multi body system. */
    SimTK::MultibodySystem MBsystem;
    /** the bodies or matter of the system. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** the forces of the system. */
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SimTK::CableTrackerSubsystem cables(MBsystem);
    /** mass properties of the fixed spot. */
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1, SimTK::Vec3(0), SimTK::UnitInertia(1)));
    /** mass properties of the structure. */
    StructureSystemForSimbody structure_multibody(structure, makeShared<MultiPolygonShape>(createStructureShape(), "StructureMultiBody"));
    /** Mass properties of the constrained spot.
     * SimTK::MassProperties(mass, center of mass, inertia)
     */
    SimTK::Body::Rigid structure_info(*structure_multibody.body_part_mass_properties_);

    SimTK::MobilizedBody::Planar tethered_spot(matter.Ground(), SimTK::Transform(SimTK::Vec3(G[0], G[1], 0.0)),
                                               structure_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spotA(matter.Ground(), SimTK::Transform(SimTK::Vec3(tethering_pointA[0], tethering_pointA[1], 0.0)),
                                           fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /*---------------------------------------------------------------------------*/
    SimTK::MobilizedBody::Weld fixed_spotB(matter.Ground(), SimTK::Transform(SimTK::Vec3(tethering_pointB[0], tethering_pointB[1], 0.0)),
                                           fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));

    // A SEASIDE PILLAR; B PORTSIDE PILLAR
    Vecd disp_cable_endA = cable_endA - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineA(cables, fixed_spotA, SimTK::Vec3(0.0, 0.0, 0.0), tethered_spot,
                                     SimTK::Vec3(disp_cable_endA[0], disp_cable_endA[1], 0.0));
    /*-----------------------------------------------------------------------------*/
    Vecd disp_cable_endB = cable_endB - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineB(cables, fixed_spotB, SimTK::Vec3(0.0, 0.0, 0.0), tethered_spot,
                                     SimTK::Vec3(disp_cable_endB[0], disp_cable_endB[1], 0.0));

    Real lengthA = G[1] + disp_cable_endA[1];
    Real lengthB = G[1] + disp_cable_endB[1];

    SimTK::CableSpring tethering_springA(forces, tethering_lineA, 3.163E5, lengthA, 2.);
    SimTK::CableSpring tethering_springB(forces, tethering_lineB, 3.163E5, lengthB, 2.);
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, -gravity_g, 0.0), 0.0);
    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    //----------------------------------------------------------------------
    std::cout << "MASS CENTER GOAL"
              << " " << G[0] << " " << G[1] << std::endl;
    std::cout << "MASS CENTER SIMBODY"
              << " " << structure_multibody.initial_mass_center_[0] << " "
              << structure_multibody.initial_mass_center_[1] << std::endl;
    std::cout << "INERTIA " << Ix << " " << Iy << " " << Iz << std::endl;
    //----------------------------------------------------------------------
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_structure(structure_multibody, MBsystem, tethered_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_on_structure(structure_multibody, MBsystem, tethered_spot, integ);

    //----------------------------------------------------------------------
    //	SimBody Output
    //----------------------------------------------------------------------
    WriteSimBodyCableData write_cable_A(sph_system, integ, tethering_springA, "A");
    WriteSimBodyCableData write_cable_B(sph_system, integ, tethering_springB, "B");
    WriteSimBodyPlanarData write_planar(sph_system, integ, tethered_spot);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    /** WaveProbe. */
    BodyRegionByCell wave_probe_buffer(water_block, makeShared<MultiPolygonShape>(createWaveGauge(), "WaveGauge"));
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>>
        wave_gauge(wave_probe_buffer, "FreeSurfaceHeight");
    /** StructureMovement. */
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_observer_position(observer_contact_with_structure, "Position", "Position");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_str_displacement("Position", observer_contact_with_structure);
    /** StructurePressureProbes. Position is updated with structure movement. */
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_fp2_position(fp2_contact_s, "Position", "Position");
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_fp3_position(fp3_contact_s, "Position", "Position");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_pressure_fp2("Pressure", fp2_contact_w);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_pressure_fp3("Pressure", fp3_contact_w);
    RestartIO restart_io(sph_system);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    structure_offset_position.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    structure_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedListWithParticleSort(100);
        wall_boundary.updateCellLinkedList();
        structure.updateCellLinkedList();
        water_block_complex.updateConfiguration();
        structure_contact.updateConfiguration();
        observer_contact_with_water.updateConfiguration();
        fp2_contact_w.updateConfiguration();
        fp3_contact_w.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 1000;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = total_physical_time;
    Real output_interval = end_time / 25;
    Real dt = 0.0;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(number_of_iterations);
    write_str_displacement.writeToFile(number_of_iterations);
    wave_gauge.writeToFile(number_of_iterations);
    write_recorded_pressure_fp2.writeToFile(number_of_iterations);
    write_recorded_pressure_fp3.writeToFile(number_of_iterations);
    write_cable_A.writeToFile(number_of_iterations);
    write_cable_B.writeToFile(number_of_iterations);
    write_planar.writeToFile(number_of_iterations);

    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            corrected_configuration_fluid.exec();
            viscous_force.exec();
            /** Viscous force exerting on structure. */
            viscous_force_on_solid.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = get_fluid_time_step_size.exec();

                pressure_relaxation.exec(dt);
                fluid_force_on_structure.exec();
                density_relaxation.exec(dt);
                /** coupled rigid body dynamics. */
                if (total_time >= relax_time)
                {
                    SimTK::State &state_for_update = integ.updAdvancedState();
                    force_on_bodies.clearAllBodyForces(state_for_update);
                    force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_structure.exec());
                    integ.stepBy(dt);
                    constraint_on_structure.exec();
                    wave_making.exec(dt);
                }
                interpolation_observer_position.exec();
                interpolation_fp2_position.exec();
                interpolation_fp3_position.exec();

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
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            water_block.updateCellLinkedListWithParticleSort(100);
            wall_boundary.updateCellLinkedList();
            structure.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            structure_contact.updateConfiguration();
            observer_contact_with_water.updateConfiguration();

            fp2_contact_w.updateConfiguration();
            fp3_contact_w.updateConfiguration();
            if (total_time >= relax_time)
            {
                write_str_displacement.writeToFile(number_of_iterations);

                wave_gauge.writeToFile(number_of_iterations);

                write_recorded_pressure_fp2.writeToFile(number_of_iterations);
                write_recorded_pressure_fp3.writeToFile(number_of_iterations);

                write_cable_A.writeToFile(number_of_iterations);
                write_cable_B.writeToFile(number_of_iterations);

                write_planar.writeToFile(number_of_iterations);
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

    if (sph_system.GenerateRegressionData())
    {
        write_str_displacement.generateDataBase(1.0e-3);
        write_recorded_pressure_fp2.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_str_displacement.testResult();
        write_recorded_pressure_fp2.testResult();
    }

    return 0;
}
