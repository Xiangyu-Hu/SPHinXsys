/**
 * @file twisting_column.cpp
 * @brief This is an example of solid with classic neohookean model
 * to demonstrate the robustness of the formulation with Kirchhoff stress decomposition.
 * @author Chi Zhang  and Xiangyu Hu
 * @ref 	DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"

#include "twisting_column.h" /**< Case setup for this example. */

using namespace SPH;

int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    /** create a body with corresponding material, particles and reaction model. */
    SolidBody column(system, makeShared<Column>("Column"));
    column.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    column.generateParticles<ParticleGeneratorLattice>();
    /** Define Observer. */
    ObserverBody my_observer(system, "MyObserver");
    my_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<InitialCondition> initial_condition(column);
    /** Corrected configuration. */
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(column_inner);
    /** Time step size calculation. We use CFL = 0.5 due to the very large twisting speed. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(column, 0.5);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half(column_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(column_inner);
    /** Constrain the holder. */
    BodyRegionByParticle holder(column, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_holder), halfsize_holder, "Holder"));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_velocity("Velocity", io_environment, my_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", io_environment, my_observer_contact);
    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    initial_condition.exec();
    corrected_configuration.exec();
    write_states.writeToFile(0);
    write_displacement.writeToFile(0);
    write_velocity.writeToFile(0);
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    int ite = 0;
    Real end_time = 0.5;
    Real output_period = end_time / 250.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
            write_displacement.writeToFile(ite);
            write_velocity.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_displacement.generateDataBase(0.005);
        write_velocity.generateDataBase(0.005);
    }
    else
    {
        write_displacement.testResult();
        write_velocity.testResult();
    }

    return 0;
}
