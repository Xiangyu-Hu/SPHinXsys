/**
 * @file 	2d_eulerian_flow_around_cylinder_LG.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder coupling with the Laguerre Gauss kernel.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;                        /**< Channel length. */
Real DH = 10.0;                        /**< Channel height. */
Real resolution_ref = 1.0 / 4.0;       /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Vec2d cylinder_center(4.0, DH / 2.0);  /**< Location of the cylinder center. */
Real cylinder_radius = 1.0;            /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	User defined boundary condition for far field
//----------------------------------------------------------------------
class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
  public:
    explicit FarFieldBoundary(BaseInnerRelation &inner_relation)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation)
    {
        rho_farfield_ = rho0_f;
        sound_speed_ = c_f;
        vel_farfield_ = Vecd(U_f, 0.0);
    };
    virtual ~FarFieldBoundary() {};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(true);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    ComplexShape water_block_shape("WaterBlock");
    water_block_shape.add<GeometricShapeBox>(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge), "OuterBoundary");
    water_block_shape.subtract<GeometricShapeBall>(cylinder_center, cylinder_radius);
    FluidBody water_block(sph_system, water_block_shape);
    water_block.getSPHAdaptation().resetKernel<KernelLaguerreGauss>();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);

    GeometricShapeBall cylinder_shape(cylinder_center, cylinder_radius, "Cylinder");
    SolidBody cylinder(sph_system, cylinder_shape);
    cylinder.defineAdaptationRatios(1.3, 2.0);
    cylinder.getSPHAdaptation().resetKernel<KernelLaguerreGauss>();
    cylinder.defineMaterial<Solid>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        LevelSetShape *outer_level_set_shape = water_block.defineComponentLevelSetShape("OuterBoundary", 2.0)
                                                   ->writeLevelSet(sph_system);
        water_block.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_water_block_outer_surface(water_block, "OuterBoundary");

        LevelSetShape *cylinder_level_set_shape = cylinder.defineBodyLevelSetShape(2.0)
                                                      ->writeLevelSet(sph_system);
        cylinder.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_cylinder_surface(cylinder);

        Inner<> cylinder_inner(cylinder);
        Inner<> water_block_inner(water_block);
        Contact<> water_block_contact(water_block, {&cylinder});
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(sph_system);
        auto &main_methods = sph_solver.addParticleMethodContainer(par_host);
        auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

        auto &cylinder_cell_linked_list = main_methods.addCellLinkedListDynamics(cylinder);
        auto &water_block_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
        auto &cylinder_update_inner_relation = main_methods.addRelationDynamics(cylinder_inner);
        auto &water_block_update_complex_relation = main_methods.addRelationDynamics(water_block_inner, water_block_contact);

        auto &random_cylinder_particles = main_methods.addStateDynamics<RandomizeParticlePositionCK>(cylinder);
        auto &random_water_block_particles = main_methods.addStateDynamics<RandomizeParticlePositionCK>(water_block);

        auto &cylinder_relaxation_residual =
            main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(cylinder_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(cylinder, *cylinder_level_set_shape);
        auto &cylinder_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(cylinder);
        auto &cylinder_update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(cylinder);
        auto &cylinder_level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_cylinder_surface);

        auto &water_block_relaxation_residual =
            main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(water_block_inner)
                .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_block_contact)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(water_block, *outer_level_set_shape);
        auto &water_block_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(water_block);
        auto &water_block_update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(water_block);
        auto &water_block_level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_water_block_outer_surface);
        //----------------------------------------------------------------------
        //	Run on CPU after relaxation finished and output results.
        //----------------------------------------------------------------------
        auto &cylinder_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(cylinder);
        auto &water_block_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(water_block);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
        body_state_recorder.addToWrite<Vecd>(cylinder, "NormalDirection");
        body_state_recorder.addToWrite<Vecd>(water_block, "NormalDirection");
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&cylinder, &water_block});
        write_particle_reload_files.addToReload<Vecd>(cylinder, "NormalDirection");
        write_particle_reload_files.addToReload<Vecd>(water_block, "NormalDirection");
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        random_cylinder_particles.exec();
        random_water_block_particles.exec();
        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        body_state_recorder.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            cylinder_cell_linked_list.exec();
            cylinder_update_inner_relation.exec();

            cylinder_relaxation_residual.exec();
            Real cylinder_relaxation_step = cylinder_relaxation_scaling.exec();
            cylinder_update_particle_position.exec(cylinder_relaxation_step);
            cylinder_level_set_bounding.exec();

            water_block_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            water_block_relaxation_residual.exec();
            Real water_block_relaxation_step = water_block_relaxation_scaling.exec();
            water_block_update_particle_position.exec(water_block_relaxation_step);
            water_block_level_set_bounding.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;
        cylinder_normal_direction.exec();
        water_block_normal_direction.exec();
        write_particle_reload_files.writeToFile();
        return 0;
    }
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    water_block.getSPHAdaptation().resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");
    cylinder.getSPHAdaptation().resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //	Note that the same relation should be defined only once.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ContactRelation water_block_contact(water_block, {&cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_wall_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(water_block_inner, water_block_contact);
    InteractionDynamics<SmearedSurfaceIndication> smeared_surface(water_block_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> cylinder_kernel_correction_matrix(cylinder_inner, cylinder_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> water_block_kernel_correction_matrix(water_block_inner, water_block_contact);
    InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block, 0.5);
    InteractionWithUpdate<FarFieldBoundary> variable_reset_in_boundary_condition(water_block_inner);
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    write_real_body_states.addToWrite<Vecd>(water_block, "Velocity");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>>
        write_total_viscous_force_from_fluid(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_pressure_force_from_fluid_body(cylinder, "PressureForceFromFluid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    surface_indicator.exec();
    smeared_surface.exec();
    variable_reset_in_boundary_condition.exec();
    cylinder_kernel_correction_matrix.exec();
    water_block_kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 80.0;
    Real output_interval = 5.0; /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            viscous_force.exec();
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);

            integration_time += dt;
            physical_time += dt;
            variable_reset_in_boundary_condition.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        viscous_force_from_fluid.exec();
        pressure_force_from_fluid.exec();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        write_total_pressure_force_from_fluid_body.writeToFile(number_of_iterations);

        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    write_total_viscous_force_from_fluid.testResult();

    return 0;
}
