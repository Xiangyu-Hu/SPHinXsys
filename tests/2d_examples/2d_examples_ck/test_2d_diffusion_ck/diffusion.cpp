/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test using splitting solver.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys_ck.h" //SPHinXsys Library
#include <gtest/gtest.h>
using namespace SPH; // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 2.0;
Real H = 0.4;
Real resolution_ref = H / 40.0;
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 1.0e-4;
Real end_time = 1.0; // end time of the simulation
//----------------------------------------------------------------------
// Define extra classes which are used in the main program.
// These classes are defined under the namespace of SPH.
//----------------------------------------------------------------------
class InitialDistribution : public ReturnFunction<Real>
{
  public:
    InitialDistribution(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser){};

        Real operator()(const Vecd &position)
        {
            if (position[0] >= 0.45 && position[0] <= 0.55)
            {
                return 1.0;
            }
            if (position[0] >= 1.0)
            {
                return exp(-0.25 * ((position[0] - 1.5) * (position[0] - 1.5)) / diffusion_coeff);
            }
            return 0.0;
        };
    };
};
//----------------------------------------------------------------------
//	Google test items
//----------------------------------------------------------------------
Real approximated_value(1.0);
Real expected_value = 1.0 / sqrt(end_time + 1.0);
TEST(Diffusion, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_value), 5.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_value << std::endl;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox diffusion_block(BoundingBox(Vecd(0.0, 0.0), Vecd(L, H)), "DiffusionBlock");
    SolidBody diffusion_body(sph_system, diffusion_block);
    diffusion_body.defineClosure<Solid, IsotropicDiffusion>(
        Solid(), ConstructArgs(diffusion_species_name, diffusion_coeff));
    diffusion_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of fluid observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    StdVec<Vecd> observation_location = {Vecd(1.5, 0.2)};
    temperature_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par);
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
    auto &diffusion_body_cell_linked_list = main_methods.addCellLinkedListDynamics(diffusion_body);
    auto &diffusion_body_update_inner_relation = main_methods.addRelationDynamics(diffusion_body_inner);
    auto &observer_update_contact_relation = main_methods.addRelationDynamics(observer_contact);

    auto &diffusion_body_linear_correction_matrix =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(diffusion_body_inner);
    auto &setup_diffusion_initial_condition =
        main_methods.addStateDynamics<
            InitialCondition<SPHBody, InitialDistribution>>(diffusion_body, diffusion_species_name);

    auto &diffusion_relaxation =
        main_methods.addInteractionDynamics<
            Dissipation, Splitting, IsotropicDiffusion>(diffusion_body_inner, diffusion_species_name);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(diffusion_body, diffusion_species_name);
    auto &observe_temperature = main_methods.addIODynamics<
        ObservedQuantityRecording, Real, RestoringCorrection>(diffusion_species_name, observer_contact);
    auto &reduce_total_species = main_methods.addReduceDynamics<QuantitySum<Real>>(diffusion_body, diffusion_species_name);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(end_time);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t time_steps = 0;
    Real output_peroid = 0.1;
    Real observe_interval = 0.1 * output_peroid;
    auto &state_recording = time_stepper.addTriggerByInterval(output_peroid);
    auto &observation_recording = time_stepper.addTriggerByInterval(observe_interval);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    diffusion_body_cell_linked_list.exec();
    diffusion_body_update_inner_relation.exec();
    observer_update_contact_relation.exec();

    diffusion_body_linear_correction_matrix.exec();
    setup_diffusion_initial_condition.exec();
    Real initial_total_species = reduce_total_species.exec();
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval_computing;
    TimeInterval interval_output;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    observe_temperature.writeToFile(time_steps);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real diffusion_dt = time_stepper.incrementPhysicalTime(0.5);
        diffusion_relaxation.exec(diffusion_dt);
        time_steps += 1;
        interval_computing += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        time_instance = TickCount::now();
        if (state_recording())
        {
            body_state_recorder.writeToFile();
        }

        if (observation_recording())
        {
            observe_temperature.writeToFile(time_steps);
        }

        if (time_steps % 10 == 0)
        {
            std::cout << "N=" << time_steps << " Time: "
                      << time_stepper.getPhysicalTime() << "	dt: " << diffusion_dt << "\n";
            std::cout << "Initial total species: " << initial_total_species
                      << "	Present total species: " << reduce_total_species.exec() << "\n";
        }
        interval_output += TickCount::now() - time_instance;
    }

    TimeInterval tt = TickCount::now() - t1 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    approximated_value = *observe_temperature.getObservedQuantity();
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
