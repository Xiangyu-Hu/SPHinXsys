/**
 * @file 	stfb_ck.cpp
 * @brief 	This is the test for 2D still floating body using compute kernel.
 * @author   Nicolò Salis and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 10.0; /**< TOTAL SIMULATION TIME*/
Real DL = 3.0;                   /**< Tank length. */
Real DH = 4.0;                   /**< Tank height. */
Real WH = 2.0;                   /**< Water block height. */
Real L = 1.0;                    /**< Base of the floating body. */
Real particle_spacing_ref = L / 20;
Real BW = particle_spacing_ref * 4.0; /**< Extending width for BCs. */
BoundingBoxd system_domain_bounds(Vec2d(-DL - BW, -DH - BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                    /**< Reference density of fluid. */
Real gravity_g = 9.81;                   /**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                   /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * DL, 0.5 * WH);
Vec2d water_block_translation = Vec2d(0.0, -0.5 * WH);
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(0.0, 0.0);
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = Vec2d(0.0, 0.0);
//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 700;
/**< Area of the solid structure*/
Real FlStA = L * L;
/**< Density of the solid structure*/
Real rho_s = StructureMass / FlStA;
/* Equilibrium position of the solid structure*/
Real H = -(rho_s / rho0_f * L - L / 2); /**< initial placement of float body */

Real bcmx = 0;
Real bcmy = H + 0;
Vec2d G(bcmx, bcmy);
Real Ix = L * L * L * L / 3;
Real Iy = L * L * L * L / 3;
Real Iz = StructureMass / 12 * (L * L + L * L);

/** Structure observer position*/
Vec2d obs = G;
/** Structure definition*/
Vec2d structure_halfsize = Vec2d(0.5 * L, 0.5 * L);
Vec2d structure_translation = Vec2d(0.0, H);
//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------
class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    StructureSystemForSimbody(SPHBody &sph_body, Shape &shape)
        : SolidBodyPartForSimbody(sph_body, shape)
    {
        // Vec2d mass_center(G[0], G[1]);
        // initial_mass_center_ = SimTKVec3(mass_center[0], mass_center[1], 0.0);
        body_part_mass_properties_ =
            mass_properties_keeper_
                .createPtr<SimTK::MassProperties>(StructureMass, SimTKVec3(0.0), SimTK::UnitInertia(Ix, Iy, Iz));
    }
};
//----------------------------------------------------------------------
//	Dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(water_block_translation), water_block_halfsize);
        subtract<GeometricShapeBox>(Transform(structure_translation), structure_halfsize);
    }
};
//----------------------------------------------------------------------
//	create measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vec2d gauge_halfsize = Vec2d(0.5 * h, 0.5 * DH);
Vec2d gauge_translation = Vec2d(DL / 3, 0.5 * DH);
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "Mass " << StructureMass << " str_sup " << FlStA << " rho_s " << rho_s << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    GeometricShapeBox structure_shape(Transform(structure_translation), structure_halfsize, "Structure");
    SolidBody structure(sph_system, structure_shape);
    structure.defineMaterial<Solid>(rho_s);
    structure.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer(sph_system, "Observer");
    observer.defineAdaptationRatios(1.15, 2.0);
    observer.generateParticles<ObserverParticles>(StdVec<Vecd>{obs});
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    GeometricShapeBox wave_probe_buffer_shape(Transform(gauge_translation), gauge_halfsize, "FreeSurfaceGauge");
    BodyRegionByCell wave_probe_buffer(water_block, wave_probe_buffer_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_block_contact(water_block, {&wall_boundary, &structure});
    Contact<> structure_contact(structure, {&water_block});
    Contact<> observer_contact(observer, {&structure}, ConfigType::Lagrangian);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary);  // run on CPU
    auto &structure_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(structure); // run on CPU

    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &structure_cell_linked_list = main_methods.addCellLinkedListDynamics(structure);
    auto &water_block_update_complex_relation = main_methods.addRelationDynamics(water_block_inner, water_block_contact);
    auto &structure_update_contact_relation = main_methods.addRelationDynamics(structure_contact);
    auto &observer_update_contact_relation = main_methods.addRelationDynamics(observer_contact);
    auto &particle_sort = main_methods.addSortDynamics(water_block);

    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(water_block, gravity);
    auto &water_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);

    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact);

    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner);
    auto &fluid_acoustic_step_2nd_half_with_wall =
        main_methods.addInteractionDynamics<
            fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact);
    fluid_acoustic_step_2nd_half.addPostContactInteraction(fluid_acoustic_step_2nd_half_with_wall);

    auto &fluid_density_regularization =
        main_methods.addInteractionDynamicsWithUpdate<
                        fluid_dynamics::DensityRegularization, FreeSurface, AllParticles>(water_block_inner)
            .addPostContactInteraction(water_block_contact);

    auto &fluid_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_f);
    auto &fluid_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);

    auto &fluid_viscous_force = main_methods.addInteractionDynamicsWithUpdate<
        fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(water_block_inner);
    auto &fluid_viscous_force_from_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(water_block_contact);
    fluid_viscous_force.addPostContactInteraction(fluid_viscous_force_from_wall);

    auto &viscous_force_on_structure =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::ViscousForceFromFluid, std::remove_reference_t<decltype(fluid_viscous_force_from_wall)>>(structure_contact);
    auto &pressure_force_on_structure =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::PressureForceFromFluid, std::remove_reference_t<decltype(fluid_acoustic_step_2nd_half_with_wall)>>(structure_contact);
    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);                // the bodies or matter of the system
    SimTK::GeneralForceSubsystem forces(MBsystem);                 // the forces of the system
    SimbodyStateEngine simbody_state_engine(sph_system, MBsystem); // the state engine of the system

    StructureSystemForSimbody structure_multibody(structure, structure_shape);
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
    SimTK::MobilizedBody::Planar structure_mob(matter.Ground(), SimTK::Transform(SimTKVec3(G[0], G[1], 0.0)),
                                               structure_info, SimTK::Transform(SimTKVec3(0.0, 0.0, 0.0)));
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
    SimTK::State simbody_state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    auto &force_on_structure = main_methods.addReduceDynamics<
        solid_dynamics::TotalForceOnBodyPartForSimBodyCK>(structure_multibody, MBsystem, structure_mob, integ);
    auto &constraint_on_structure = main_methods.addStateDynamics<
        solid_dynamics::ConstraintBodyPartBySimBodyCK>(structure_multibody, MBsystem, structure_mob, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);
    auto &write_real_body_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &wave_gauge = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, UpperFrontInAxisDirectionCK>(wave_probe_buffer, "FreeSurfaceHeight");
    auto &write_structure_position = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>("Position", observer_contact);

    SingularVariable<SimTK::SpatialVec> sv_action_on_structure("ActionOnStructure", SimTK::SpatialVec(0));
    SingularVariableRecording<SimTK::SpatialVec> action_on_structure_recording(sph_system, &sv_action_on_structure);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();
    structure_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	From here, methods run on device.
    //----------------------------------------------------------------------
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Construct static cell linked list and Lagrangian configurations.
    //----------------------------------------------------------------------
    wall_cell_linked_list.exec();
    structure_cell_linked_list.exec();
    observer_update_contact_relation.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    Real restart_time = 0.0;
    size_t restart_step = sph_system.RestartStep();
    if (restart_step != 0)
    {
        restart_time = restart_io.readRestartFiles(restart_step);
        structure_cell_linked_list.exec();

        simbody_state_engine.readStateFromXml(restart_step, simbody_state);
        simbody_state.setTime(Real(restart_time));
    }
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(simbody_state);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper time_stepper(sph_system, total_physical_time, restart_time);
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    auto &trigger_FSI = time_stepper.addTriggerByPhysicalTime(1.0);
    size_t advection_steps = restart_step;
    int screening_interval = 100;
    int restart_interval = 1000;
    int observation_interval = screening_interval / 2;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 100.0);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    structure_update_contact_relation.exec();

    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    fluid_viscous_force.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    write_structure_position.writeToFile();
    wave_gauge.writeToFile();
    /** statistics for computing time. */
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_FSI;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        if (trigger_FSI())
        {
            pressure_force_on_structure.exec();
            SimTK::State &state_for_update = integ.updAdvancedState();
            force_on_bodies.clearAllBodyForces(state_for_update);
            sv_action_on_structure.setValue(force_on_structure.exec());
            force_on_bodies.setOneBodyForce(state_for_update, structure_mob, sv_action_on_structure.getValue());
            integ.stepBy(acoustic_dt);
            constraint_on_structure.exec();
        }
        interval_FSI += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            water_update_particle_position.exec();

            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Physical Time = " << time_stepper.getPhysicalTime()
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (trigger_FSI() && advection_steps % observation_interval == 0)
            {
                write_structure_position.writeToFile(advection_steps);
                wave_gauge.writeToFile(advection_steps);
                action_on_structure_recording.writeToFile(advection_steps);
            }

            if (advection_steps % restart_interval == 0)
            {
                restart_io.writeToFile(advection_steps);
                simbody_state_engine.writeStateToXml(advection_steps, integ);
            }

            if (trigger_FSI() && state_recording())
            {
                write_real_body_states.writeToFile();
            }

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (advection_steps % 100)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            structure_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            structure_update_contact_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_viscous_force.exec();
            if (trigger_FSI())
            {
                viscous_force_on_structure.exec();
            };
            interval_advection_step += TickCount::now() - time_instance;
        }
    }

    //----------------------------------------------------------------------
    // Summary for wall time used for real computations.
    //----------------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_FSI = "
              << interval_FSI.seconds() << "\n";
    //----------------------------------------------------------------------
    // Post-run regression test to ensure that the case is validated
    //----------------------------------------------------------------------

    if (sph_system.GenerateRegressionData())
    {
        write_structure_position.generateDataBase(0.001);
        wave_gauge.generateDataBase(0.001);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_structure_position.testResult();
        wave_gauge.testResult();
    }

    return 0;
}
