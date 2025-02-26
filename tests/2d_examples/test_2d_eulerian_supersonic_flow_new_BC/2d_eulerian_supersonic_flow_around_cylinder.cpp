/**
 * @file 	2d_eulerian_supersonic_flow_around_cylinder.cpp
 * @brief 	This is the compressible test for Eulerian supersonic flow around a cylinder
            with the FVM boundary algorithm, i.e., zero-order consistency, in the SPHinXsys.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "2d_eulerian_supersonic_flow_around_cylinder.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(false);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(false);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("FluidBlock"));
    fluid_block.getSPHAdaptation().resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    fluid_block.defineBodyLevelSetShape();
    fluid_block.defineMaterial<CompressibleFluid>(rho_reference, heat_capacity_ratio);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? fluid_block.generateParticlesWithReserve<BaseParticles, Reload>(ghost_boundary, fluid_block.getName())
        : fluid_block.generateParticlesWithReserve<BaseParticles, Lattice>(ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //	Note that the same relation should be defined only once.
    //----------------------------------------------------------------------
    InnerRelation fluid_block_inner(fluid_block);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(fluid_block);
        BodyStatesRecordingToVtp write_real_body_states(fluid_block);
        ReloadParticleIO write_real_body_particle_reload_files(fluid_block);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(fluid_block_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;

        write_real_body_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<FreeSurfaceIndication<Inner<>>> surface_indicator(fluid_block_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfHLLCWithLimiterRiemann> pressure_relaxation(fluid_block_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfHLLCWithLimiterRiemann> density_and_energy_relaxation(fluid_block_inner);
    SimpleDynamics<SupersonicFlowInitialCondition> initial_condition(fluid_block);
    ReduceDynamics<fluid_dynamics::EulerianCompressibleAcousticTimeStepSize> get_fluid_time_step_size(fluid_block, 0.1);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> kernel_correction_matrix(fluid_block_inner);
    InteractionDynamics<KernelGradientCorrectionInner> kernel_gradient_update(fluid_block_inner);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    surface_indicator.exec();
    //----------------------------------------------------------------------
    //	Prepare the ghost particles and its configuration and Boundary conditions setup.
    //	Ghost kernel gradient update is to make the boundary particles strictly achieve zero-order consistency.
    //  strictly achieve zero-order consistency
    //----------------------------------------------------------------------
    GhostCreationInESPH ghost_creation_for_boundary_condition(fluid_block_inner, ghost_boundary);
    SupersonicFlowBoundaryConditionSetup boundary_condition_setup(fluid_block_inner, ghost_creation_for_boundary_condition);
    InteractionWithUpdate<GhostKernelGradientUpdate> ghost_kernel_gradient_update(fluid_block_inner);
    boundary_condition_setup.resetBoundaryConditions();
    kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    ghost_kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MaximumSpeed>> write_maximum_speed(fluid_block);
    write_real_body_states.addToWrite<int>(fluid_block, "Indicator");
    write_real_body_states.addToWrite<Vecd>(fluid_block, "Velocity");
    write_real_body_states.addToWrite<Real>(fluid_block, "Pressure");
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    Real end_time = 40.0;
    Real output_interval = 1.0; /**< time stamps for output. */
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
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_and_energy_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();

            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
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
        write_maximum_speed.generateDataBase(1.0e-2);
    }
    else
    {
        write_maximum_speed.testResult();
    }

    return 0;
}
