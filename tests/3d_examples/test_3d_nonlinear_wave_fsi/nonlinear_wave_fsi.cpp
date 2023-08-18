/**
 * @file 	nonlinear_wave_fsi.h
 * @brief 	This is the 3d case file for wave impact with tension leg moored floating structure.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
#include "io_simbody_cable.h"   //output for cable data
#include "io_simbody_free.h"    //output for free body data
#include "nonlinear_wave_fsi.h" //header for this case

int main(int ac, char *av[])
{
    std::cout << "Mass " << StructureMass << " Volume " << StructureVol << " rho_str " << StructureDensity << std::endl;

    //----------------------------------------------------------------------
    //	Build up the Structure Relax SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system_fit(system_domain_bounds, particle_spacing_structure);
    system_fit.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment_fit(system_fit);
    SolidBody structure_fit(system_fit, makeShared<FloatingStructure>("Structure_Fit"));
    structure_fit.defineAdaptation<ParticleRefinementNearSurface>(1.3, 0.7, 3);
    structure_fit.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment_fit);
    structure_fit.defineParticlesAndMaterial<SolidParticles, Solid>(StructureDensity);
    structure_fit.generateParticles<ParticleGeneratorMultiResolution>();
    structure_fit.addBodyStateForRecording<Real>("SmoothingLengthRatio");

    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    AdaptiveInnerRelation structure_adaptive_inner(structure_fit);

    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(structure_fit);
    /** A  Physics relaxation step. */
    relax_dynamics::RelaxationStepInner relaxation_step_inner(structure_adaptive_inner, true);
    SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(structure_fit);
    /** Write the particle reload files. */
    ReloadParticleIO write_particle_reload_files(io_environment_fit, structure_fit);

    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    update_smoothing_length_ratio.exec();
    structure_fit.updateCellLinkedList();

    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 100)
    {
        update_smoothing_length_ratio.exec();
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 10 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
        }
    }
    std::cout << "The physics relaxation process of imported model finish !" << std::endl;

    write_particle_reload_files.writeToFile(0);

    //----------------------------------------------------------------------
    //	Build up the Main environment of a SPHSystem with global controls.
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
    water_block.addBodyStateForRecording<Vecd>("Position");
    water_block.addBodyStateForRecording<Real>("Pressure");

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody structure(system, makeShared<FloatingStructure>("Structure"));
    structure.defineParticlesAndMaterial<SolidParticles, Solid>(StructureDensity);
    structure.generateParticles<ParticleGeneratorReload>(io_environment, "Structure_Fit");

    ObserverBody observer(system, "Observer");
    observer.defineAdaptationRatios(h, 2.0);
    observer.generateParticles<ObserverParticleGenerator>(
        StdVec<Vecd>{obs});

    ObserverBody WMobserver(system, "WMObserver");
    WMobserver.defineAdaptationRatios(h, 2.0);
    Vecd WMpos0 = Vecd(0.0, -Maker_width / 2, HWM / 2);
    WMobserver.generateParticles<ObserverParticleGenerator>(
        StdVec<Vecd>{WMpos0});

    //---------------------------------------------------------
    // PRESSURE PROBES
    //---------------------------------------------------------
   
    ObserverBody fp1(system, "FP1");
    Real fp1x = Strx - 0.12;
    Real fp1y = Stry;
    Real fp1z = 1.013;
    StdVec<Vecd> fp1l = {Vecd(fp1x, fp1y, fp1z)};
    fp1.generateParticles<ObserverParticleGenerator>(fp1l);
    
    ObserverBody bp1(system, "BP1");
    Real bp1x = Strx - 0.295;
    Real bp1y = Stry + .035;
    Real bp1z = 0.933;
    StdVec<Vecd> bp1l = {Vecd(bp1x, bp1y, bp1z)};
    bp1.generateParticles<ObserverParticleGenerator>(bp1l);

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
    ContactRelation WMobserver_contact_with_water(WMobserver, {&water_block});
    ContactRelation WMobserver_contact_with_wall(WMobserver, {&wall_boundary});


    ContactRelation fp1_contact_s(fp1, {&structure});
    ContactRelation bp1_contact_s(bp1, {&structure});

    ContactRelation fp1_contact_w(fp1, {&water_block});
    ContactRelation bp1_contact_w(bp1, {&water_block});

    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<OffsetInitialPosition> structure_offset_position(structure, offset);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> structure_normal_direction(structure);
    /** Time step initialization, add gravity. */
    SimpleDynamics<TimeStepInitialization> initialize_time_step_to_fluid(water_block, makeShared<Gravity>(Vecd(0.0, 0.0, -gravity_g)));
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    /** time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** pressure relaxation using Verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannCorrectWithWall> pressure_relaxation(water_block_complex);
    // Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
    /** corrected strong configuration. */
    InteractionWithUpdate<CorrectedConfigurationComplex> corrected_configuration_fluid(water_block_complex, 2, 0.4);
    InteractionWithUpdate<CorrectedConfigurationInner> structure_corrected_configuration(structure_inner);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Damp waves */
    Vecd translation_damping(0.5 * DW, 9.5 , 0.5 * HWM);
    Vecd damping(0.5 * DW, 0.5, 0.5 * HWM);
    BodyRegionByCell damping_buffer(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_damping), damping));
    SimpleDynamics<fluid_dynamics::DampingBoundaryCondition> damping_wave(damping_buffer);
    /** Fluid force on structure. */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(structure_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_structure(structure_contact, viscous_force_on_solid);
    /** constrain region of the part of wall boundary. */
    BodyRegionByParticle wave_maker(wall_boundary, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_wave_maker), wave_maker_shape));
    SimpleDynamics<WaveMaking> wave_making(wave_maker);

    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    std::cout << "Volume " << StructureVol << std::endl;
    std::cout << "MASS " << StructureDensity * StructureVol << std::endl;
    std::cout << "MASS CENTER " << G << std::endl;
    std::cout << "INERTIA " << Ix << " " << Iy << " " << Iz << std::endl;
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
    StructureSystemForSimbody structure_multibody(structure, makeShared<TriangleMeshShapeSTL>(stl_structure_path, translation_str, StructureScale));
    SimTK::Body::Rigid structure_info(*structure_multibody.body_part_mass_properties_);
    /** Free mobilizer */
    SimTK::MobilizedBody::Free tethered_struct(matter.Ground(), SimTK::Transform(SimTK::Vec3(translation_str[0], translation_str[1], translation_str[2])), structure_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /** Mobility of the fixed spot. */
    /*---------------------------------------------------------------------------*/
    SimTK::MobilizedBody::Weld fixed_spotAR(matter.Ground(), SimTK::Transform(SimTK::Vec3(ground_tethering_AR[0], ground_tethering_AR[1], ground_tethering_AR[2])),
                                            fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /*---------------------------------------------------------------------------*/
    SimTK::MobilizedBody::Weld fixed_spotAL(matter.Ground(), SimTK::Transform(SimTK::Vec3(ground_tethering_AL[0], ground_tethering_AL[1], ground_tethering_AL[2])),
                                            fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /*---------------------------------------------------------------------------*/
    SimTK::MobilizedBody::Weld fixed_spotBR(matter.Ground(), SimTK::Transform(SimTK::Vec3(ground_tethering_BR[0], ground_tethering_BR[1], ground_tethering_BR[2])),
                                            fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /*---------------------------------------------------------------------------*/
    SimTK::MobilizedBody::Weld fixed_spotBL(matter.Ground(), SimTK::Transform(SimTK::Vec3(ground_tethering_BL[0], ground_tethering_BL[1], ground_tethering_BL[2])),
                                            fixed_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    /*---------------------------------------------------------------------------*/
    // A SEASIDE PILLARS; B PORTSIDE PILLARS
    /*-----------------------------------------------------------------------------*/
    Vecd disp_cable_endAR = structure_tethering_AR - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineAR(cables, fixed_spotAR, SimTK::Vec3(0.0, 0.0, 0.0), tethered_struct, SimTK::Vec3(disp_cable_endAR[0], disp_cable_endAR[1], disp_cable_endAR[2]));
    /*-----------------------------------------------------------------------------*/
    Vecd disp_cable_endAL = structure_tethering_AL - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineAL(cables, fixed_spotAL, SimTK::Vec3(0.0, 0.0, 0.0), tethered_struct, SimTK::Vec3(disp_cable_endAL[0], disp_cable_endAL[1], disp_cable_endAL[2]));
    /*-----------------------------------------------------------------------------*/
    Vecd disp_cable_endBR = structure_tethering_BR - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineBR(cables, fixed_spotBR, SimTK::Vec3(0.0, 0.0, 0.0), tethered_struct, SimTK::Vec3(disp_cable_endBR[0], disp_cable_endBR[1], disp_cable_endBR[2]));
    /*-----------------------------------------------------------------------------*/
    Vecd disp_cable_endBL = structure_tethering_BL - structure_multibody.initial_mass_center_;
    SimTK::CablePath tethering_lineBL(cables, fixed_spotBL, SimTK::Vec3(0.0, 0.0, 0.0), tethered_struct, SimTK::Vec3(disp_cable_endBL[0], disp_cable_endBL[1], disp_cable_endBL[2]));

    Real lengthAR = G[2] + disp_cable_endAR[2];
    Real lengthAL = G[2] + disp_cable_endAL[2];
    Real lengthBR = G[2] + disp_cable_endBR[2];
    Real lengthBL = G[2] + disp_cable_endBL[2];

    SimTK::CableSpring tethering_springAR(forces, tethering_lineAR, 3.163E5, lengthAR, 1.);
    SimTK::CableSpring tethering_springAL(forces, tethering_lineAL, 3.163E5, lengthAL, 1.);
    SimTK::CableSpring tethering_springBR(forces, tethering_lineBR, 3.163E5, lengthBR, 1.);
    SimTK::CableSpring tethering_springBL(forces, tethering_lineBL, 3.163E5, lengthBL, 1.);
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, 0.0, -gravity_g), 0.0);
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
        force_on_structure(structure_multibody, MBsystem, tethered_struct, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_on_structure(structure_multibody, MBsystem, tethered_struct, integ);

    //----------------------------------------------------------------------
    //	Cable SimBody Output
    //----------------------------------------------------------------------
    WriteSimBodyCableData write_cable_AR(io_environment, integ, tethering_springAR, "AR");
    WriteSimBodyCableData write_cable_AL(io_environment, integ, tethering_springAL, "AL");
    WriteSimBodyCableData write_cable_BR(io_environment, integ, tethering_springBR, "BR");
    WriteSimBodyCableData write_cable_BL(io_environment, integ, tethering_springBL, "BL");
    WriteSimBodyFreeRotationMatrix write_free_body_rotation(io_environment, integ, tethered_struct);
    WriteSimBodyVelocity write_free_body_velocity(io_environment, integ, tethered_struct);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    /** WaveProbes. */
    BodyRegionByCell wave_probe_buffer(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_WGauge), WGaugeDim));
    ReducedQuantityRecording<ReduceDynamics<FreeSurfaceHeightZ>> wave_gauge(io_environment, wave_probe_buffer);

    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_observer_position(observer_contact_with_structure, "Position", "Position");

    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_str_displacement("Position", io_environment, observer_contact_with_structure);

    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_WMobserver_position(WMobserver_contact_with_wall, "Position", "Position");
    ObservedQuantityRecording<Vecd>
        write_WM_displacement("Position", io_environment, WMobserver_contact_with_wall);

    
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_fp1_position(fp1_contact_s, "Position", "Position");
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_bp1_position(bp1_contact_s, "Position", "Position");
    
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_pressure_fp1("Pressure", io_environment, fp1_contact_w);
    ObservedQuantityRecording<Real>
        write_recorded_pressure_bp1("Pressure", io_environment, bp1_contact_w);
    
    RestartIO restart_io(io_environment, system.real_bodies_);

    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = total_physical_time;
    Real output_interval = end_time / 100;
    Real dt = 0.0;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    structure_offset_position.exec();
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    structure_normal_direction.exec();
    structure_corrected_configuration.exec();

    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
        water_block.updateCellLinkedListWithParticleSort(100);
        wall_boundary.updateCellLinkedList();
        structure.updateCellLinkedList();
        water_block_complex.updateConfiguration();
        structure_contact.updateConfiguration();
        observer_contact_with_water.updateConfiguration();
        WMobserver_contact_with_water.updateConfiguration();

        fp1_contact_w.updateConfiguration();
        bp1_contact_w.updateConfiguration();
    }

    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(number_of_iterations);
    write_str_displacement.writeToFile(number_of_iterations);
    write_WM_displacement.writeToFile(number_of_iterations);
    wave_gauge.writeToFile(number_of_iterations);

    write_recorded_pressure_fp1.writeToFile(number_of_iterations);
    write_recorded_pressure_bp1.writeToFile(number_of_iterations);

    write_cable_AR.writeToFile(number_of_iterations);
    write_cable_AL.writeToFile(number_of_iterations);
    write_cable_BR.writeToFile(number_of_iterations);
    write_cable_BL.writeToFile(number_of_iterations);

    write_free_body_rotation.writeToFile(number_of_iterations);
    write_free_body_velocity.writeToFile(number_of_iterations);

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
            corrected_configuration_fluid.exec();
            viscous_acceleration.exec();
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
                    force_on_bodies.setOneBodyForce(state_for_update, tethered_struct, force_on_structure.exec());
                    integ.stepBy(dt);
                    constraint_on_structure.exec();
                    wave_making.exec(dt);
                }
                interpolation_observer_position.exec();
                interpolation_fp1_position.exec();
                interpolation_bp1_position.exec();

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
            damping_wave.exec(Dt);
            water_block.updateCellLinkedListWithParticleSort(100);
            wall_boundary.updateCellLinkedList();
            structure.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            structure_contact.updateConfiguration();
            observer_contact_with_water.updateConfiguration();
            WMobserver_contact_with_water.updateConfiguration();

            fp1_contact_w.updateConfiguration();
            bp1_contact_w.updateConfiguration();

            if (total_time >= relax_time)
            {
                write_str_displacement.writeToFile(number_of_iterations);
                write_WM_displacement.writeToFile(number_of_iterations);
                wave_gauge.writeToFile(number_of_iterations);

                write_recorded_pressure_fp1.writeToFile(number_of_iterations);
                write_recorded_pressure_bp1.writeToFile(number_of_iterations);

                write_cable_AR.writeToFile(number_of_iterations);
                write_cable_AL.writeToFile(number_of_iterations);
                write_cable_BR.writeToFile(number_of_iterations);
                write_cable_BL.writeToFile(number_of_iterations);

                write_free_body_rotation.writeToFile(number_of_iterations);
                write_free_body_velocity.writeToFile(number_of_iterations);
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
        write_str_displacement.generateDataBase(1.0e-3);
        write_recorded_pressure_fp1.generateDataBase(1.0e-3);
    }
    else if (system.RestartStep() == 0)
    {
        write_str_displacement.testResult();
        write_recorded_pressure_fp1.testResult();
    }

    return 0;
}
