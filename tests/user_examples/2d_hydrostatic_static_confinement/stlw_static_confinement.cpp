/**
 * @file 	stlw.cpp
 * @brief 	This is the case file for 2D still water.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
#include "stlw_static_confinement.h" //header for this case
#include "level_set_confinement.h"

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<TransformShape<GeometricShapeBox>>(
                                      Transform(water_block_translation), water_block_halfsize, "Structure"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<Real>("VolumetricMeasure");

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    //ComplexRelation water_block_complex(water_block_inner, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    //SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Time step initialization, add gravity. */
    SimpleDynamics<TimeStepInitialization> initialize_time_step_to_fluid(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);
    /** time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** pressure relaxation using Verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> density_relaxation(water_block_inner);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_block_inner);
    /** Define the confinement condition for wall. */
    NearShapeSurface near_surface_wall(water_block, makeShared<WallBoundary>("Wall"));
    near_surface_wall.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinementGeneral confinement_condition_wall(near_surface_wall);
    /** Push back the static confinement conditiont to corresponding dynamics. */
    update_density_by_summation.post_processes_.push_back(&confinement_condition_wall.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_wall.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.surface_bounding_);
    viscous_acceleration.post_processes_.push_back(&confinement_condition_wall.viscous_acceleration_);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    BodyRegionByCell wave_probe_buffer(water_block, makeShared<TransformShape<GeometricShapeBox>>(
                                                        Transform(gauge_translation), gauge_halfsize, "FreeSurfaceGauge"));
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>>
        wave_gauge(io_environment, wave_probe_buffer, "FreeSurfaceHeight");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    //wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
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

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = get_fluid_time_step_size.exec();

                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);

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
            //wall_boundary.updateCellLinkedList();
            water_block_inner.updateConfiguration();

            if (total_time >= relax_time)
            {
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

    if (system.GenerateRegressionData())
    {
        wave_gauge.generateDataBase(0.1);
    }
    else
    {
        wave_gauge.testResult();
    }

    return 0;
}
