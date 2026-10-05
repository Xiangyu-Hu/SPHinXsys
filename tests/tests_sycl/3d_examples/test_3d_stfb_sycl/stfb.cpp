/**
 * @file 	 stfb.cpp
 * @brief 	 This is the case file for 3D still floating body using computing kernel.
 * @author   Nicolò Salis and Xiangyu Hu
 */
#include "simbody_system.h"
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 5.0; /**< TOTAL SIMULATION TIME*/
Real DW = 3.0;                  /**< Water length. */
Real DL = 3.0;                  /**< Tank length. */
Real DH = 2.5;                  /**< Tank height. */
Real WH = 2.0;                  /**< Water block height. */
Real L = 1.0;                   /**< Base of the floating body. */
Real particle_spacing_ref = L / 10;
Real BW = particle_spacing_ref * 4.0;          /**< Extending width for BCs. */
Real Maker_width = particle_spacing_ref * 4.0; /**< Width of the wave_maker. */
BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DW + BW, DL + BW, DH + BW));
Vecd offset = Vecd::Zero();
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                  /**< Reference density of fluid. */
Real gravity_g = 9.81;                 /**< Value of gravity. */
Real U_f = 2.0 * sqrt(WH * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                 /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 700;
/**< Area of the solid structure*/
Real FlStA = L * L * L;
/**< Density of the solid structure*/
Real rho_s = StructureMass / FlStA;
/* Equilibrium position of the solid structure*/
Real H = WH - (rho_s / rho0_f * L - L / 2); /**< Strart placement of Flt Body*/

Real bcmx = DL / 2;
Real bcmy = DL / 2;
Real bcmz = H;
Vecd G(bcmx, bcmy, bcmz);
Real Ix = StructureMass / 12 * (L * L + L * L);
Real Iy = StructureMass / 12 * (L * L + L * L);
Real Iz = StructureMass / 12 * (L * L + L * L);
/** Structure observer position */
Vecd obs = G;
/** Geometry definition. */
Vecd halfsize_structure(0.5 * L, 0.5 * L, 0.5 * L);
Vecd structure_pos(G[0], G[1], G[2]);
Transform translation_str(structure_pos);
//------------------------------------------------------------------------------
// Geometric shapes used in the case
//------------------------------------------------------------------------------
class FloatingStructure : public ComplexShape
{
  public:
    explicit FloatingStructure(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_str), halfsize_structure);
    }
};
//----------------------------------------------------------------------
//	Water block
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        Vecd halfsize_water(0.5 * DW, 0.5 * DL, 0.5 * WH);
        Vecd water_pos(0.5 * DW, 0.5 * DL, 0.5 * WH);
        Transform translation_water(water_pos);
        add<GeometricShapeBox>(Transform(translation_water), halfsize_water);
        subtract<GeometricShapeBox>(Transform(translation_str), halfsize_structure);
    }
};
//----------------------------------------------------------------------
//	Wall geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wall_outer(0.5 * DW + BW, 0.5 * DL + BW, 0.5 * DH + BW);
        Vecd wall_outer_pos(0.5 * DW, 0.5 * DL, 0.5 * DH);
        Transform translation_wall_outer(wall_outer_pos);
        add<GeometricShapeBox>(Transform(translation_wall_outer), halfsize_wall_outer);

        Vecd halfsize_wall_inner(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Vecd wall_inner_pos(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Transform translation_wall_inner(wall_inner_pos);
        subtract<GeometricShapeBox>(Transform(translation_wall_inner), halfsize_wall_inner);
    }
};
//----------------------------------------------------------------------
//	create measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vecd FS_gaugeDim(0.5 * h, 0.5 * h, 0.5 * DH);
Vecd FS_gauge(DW / 3, DL / 3, 0.5 * DH);
Transform translation_FS_gauge(FS_gauge);
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
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody structure(sph_system, makeShared<FloatingStructure>("Structure"));
    structure.defineMatterMaterial<Solid>(rho_s);
    structure.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer(sph_system, "Observer");
    observer.defineAdaptationRatios(1.15, 2.0);
    observer.generateParticles<ObserverParticles>(StdVec<Vecd>{obs});

    TriangleMeshShapeBrick structure_mesh(halfsize_structure, 1, structure_pos);
    ObserverBody structure_proxy(sph_system, "StructureObserver");
    structure_proxy.generateParticles<ObserverParticles>(structure_mesh);
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    GeometricShapeBox wave_probe_buffer_shape(Transform(translation_FS_gauge), FS_gauge);
    BodyRegionByCell wave_probe_buffer(water_block, wave_probe_buffer_shape);

    BodySurfaceLayer structure_surface(structure, 1.0);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> fluid_wall_contact(water_block, wall_boundary);
    Contact<> fluid_structure_contact(water_block, structure);
    Contact<> structure_contact(structure, water_block);
    Contact<> observer_contact(observer, structure, ConfigType::Lagrangian);
    Contact<Relation<SPHBody, BodyPartByParticle>> structure_proxy_contact(
        structure_proxy, structure_surface, ConfigType::Lagrangian);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> structure_cell_linked_list(structure);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>, Contact<>>
        water_block_update_complex_relation(water_block_inner, fluid_wall_contact, fluid_structure_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>>
        structure_update_contact_relation(structure_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>>
        observer_update_contact_relation(observer_contact);
    UpdateRelation<MainExecutionPolicy, Contact<Relation<SPHBody, BodyPartByParticle>>>
        structure_proxy_update_contact_relation(structure_proxy_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);

    Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(water_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> structure_boundary_normal_direction(structure);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_update_particle_position(water_block);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        fluid_acoustic_step_1st_half(water_block_inner);
    fluid_acoustic_step_1st_half
        .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_wall_contact)
        .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_structure_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>
        fluid_acoustic_step_2nd_half(water_block_inner);

    fluid_acoustic_step_2nd_half
        .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_wall_contact)
        .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_structure_contact)
        .addGeneralPostInteraction<FSI::PressureForceFromFluid, WithUpdate, AcousticRiemannSolverCK, NoKernelCorrectionCK>(
            structure_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::CompressionSummation<Inner<>, Contact<>, Contact<>>>
        fluid_density_summation(water_block_inner, fluid_wall_contact, fluid_structure_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, WeaklyCompressibleFluid, FreeSurface>>
        fluid_density_regularization(water_block);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceCK<Inner<WithUpdate, Viscosity, NoKernelCorrectionCK>>>
        fluid_viscous_force(water_block_inner);
    fluid_viscous_force
        .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(fluid_wall_contact)
        .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(fluid_structure_contact)
        .addGeneralPostInteraction<FSI::ViscousForceFromFluid, WithUpdate, Viscosity, NoKernelCorrectionCK>(
            structure_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_block, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>> fluid_acoustic_time_step(water_block);

    ArbitraryDynamicsSequence<
        StateDynamics<MainExecutionPolicy, solid_dynamics::UpdateDisplacementFromPosition>,
        InteractionDynamicsCK<MainExecutionPolicy, Interpolation<Contact<Vecd, Relation<SPHBody, BodyPartByParticle>>>>,
        InteractionDynamicsCK<MainExecutionPolicy, Interpolation<Contact<Vecd, Relation<SPHBody, BodyPartByParticle>>>>,
        StateDynamics<MainExecutionPolicy, solid_dynamics::UpdatePositionFromDisplacement>>
        update_structure_proxy_states(
            structure,
            DynamicsArgs(structure_proxy_contact, std::string("Velocity")),
            DynamicsArgs(structure_proxy_contact, std::string("Displacement")),
            structure_proxy);
    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    SimbodySystem simbody_system;
    GeometricShapeBox fix_spot_shape(Transform(translation_str), halfsize_structure);
    SolidBodyPartForSimbodyCK structure_multibody(structure, fix_spot_shape);
    std::string structure_name = simbody_system.createRigidBody(structure_multibody);
    SimTK::MobilizedBody &structure_mob = simbody_system.createMobilizedBody(structure_name, "Ground", "Planar");
    simbody_system.addUniformGravity(Vec3d(0.0, 0.0, -gravity_g));
    simbody_system.initializeStateForIntegrator();
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamicsCK<MainExecutionPolicy, solid_dynamics::TotalForceOnBodyPartForSimBodyCK>
        force_on_structure(structure_multibody, simbody_system);
    StateDynamics<MainExecutionPolicy, solid_dynamics::ConstraintBodyPartBySimBodyCK>
        constraint_on_structure(structure_multibody, simbody_system);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_real_body_states(sph_system);
    BodyStatesRecordingToTriangleMeshVtp write_structure_surface(structure_proxy, structure_mesh);
    write_structure_surface.addToWrite<Vecd>(structure_proxy, "Velocity");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<
        MainExecutionPolicy, UpperFrontInAxisDirectionCK<BodyRegionByCell>>>
        wave_gauge(wave_probe_buffer, "FreeSurfaceHeight");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd>>
        write_structure_position(observer_contact, "Position");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();
    structure_boundary_normal_direction.exec();
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    structure_cell_linked_list.exec();

    water_block_update_complex_relation.exec();
    structure_update_contact_relation.exec();
    observer_update_contact_relation.exec();
    structure_proxy_update_contact_relation.exec();
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = total_physical_time;
    Real output_interval = end_time / 200;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    write_structure_position.writeToFile(number_of_iterations);
    wave_gauge.writeToFile(number_of_iterations);
    write_structure_surface.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_interval)
        {
            fluid_density_summation.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            fluid_viscous_force.exec();

            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);

                if (total_time >= relax_time) // coupled rigid body dynamics.
                {
                    TorqueAndForce torque_and_force = force_on_structure.exec();
                    simbody_system.updateForceOnBody(structure_mob, torque_and_force);
                    simbody_system.stepSimbodySystemBy(acoustic_dt);
                    constraint_on_structure.exec();
                }

                relaxation_time += acoustic_dt;
                integral_time += acoustic_dt;
                total_time += acoustic_dt;
                if (total_time >= relax_time)
                    sv_physical_time->incrementValue(acoustic_dt);
            }
            water_update_particle_position.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << "	Total Time = " << total_time
                          << "	Physical Time = " << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";
            }
            number_of_iterations++;
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            structure_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            structure_update_contact_relation.exec();

            if (total_time >= relax_time)
            {
                write_structure_position.writeToFile(number_of_iterations);
                wave_gauge.writeToFile(number_of_iterations);
            }
        }

        TickCount t2 = TickCount::now();
        if (total_time >= relax_time)
        {
            write_real_body_states.writeToFile();
            update_structure_proxy_states.exec();
            write_structure_surface.writeToFile();
        }
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_structure_position.generateDataBase(1e-3);
        wave_gauge.generateDataBase(1e-3);
    }
    else
    {
        write_structure_position.testResult();
        wave_gauge.testResult();
    }

    return 0;
}
