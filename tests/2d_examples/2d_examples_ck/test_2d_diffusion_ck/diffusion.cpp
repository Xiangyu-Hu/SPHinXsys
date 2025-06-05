/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test using splitting solver.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys_ck.h" //SPHinXsys Library
#include <gtest/gtest.h>
using namespace SPH; // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 2.0;
Real H = 0.4;
Real resolution_ref = H / 80.0;
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
std::string phi_pair_wise = "PhiPairWise";
std::string phi_projection = "PhiProjection";
std::string phi_explicit = "PhiExplicit";
const Real diffusion_coeff = 1.0e-4;
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
Real expected_value = 1.0 / sqrt(end_time + 1.0);

Real approximated_pair_wise(1.0);
TEST(DiffusionPairWise, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_pair_wise), 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_pair_wise << std::endl;
};

Real approximated_projection(1.0);
TEST(DiffusionProjection, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_projection), 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_projection << std::endl;
};

Real approximated_explicit(1.0);
TEST(DiffusionExplicit, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_explicit), 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_explicit << std::endl;
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
        Solid(), ConstructArgs(phi_pair_wise, diffusion_coeff));
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

    auto &initial_condition_pair_wise =
        main_methods.addStateDynamics<
            InitialCondition<SPHBody, InitialDistribution>>(diffusion_body, phi_pair_wise);
    auto &diffusion_relaxation_pair_wise =
        main_methods.addInteractionDynamics<
            PairwiseDissipation, Splitting, IsotropicDiffusion>(diffusion_body_inner, phi_pair_wise);

    auto &initial_condition_projection =
        main_methods.addStateDynamics<
            InitialCondition<SPHBody, InitialDistribution>>(diffusion_body, phi_projection);
    using MainExecutionPolicy = execution::ParallelPolicy; // define execution policy for this case
    ImplicitDissipation<MainExecutionPolicy, Inner<IsotropicDiffusion>>
        diffusion_relaxation_projection(diffusion_body_inner, phi_projection, 1.0e-3);

    IsotropicDiffusion isotropic_diffusion_explict(phi_explicit, diffusion_coeff);
    auto &initial_condition_explicit =
        main_methods.addStateDynamics<
            InitialCondition<SPHBody, InitialDistribution>>(diffusion_body, phi_explicit);
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, &isotropic_diffusion_explict);
    auto &diffusion_relaxation_explicit =
        main_methods.addRungeKuttaSequence<
            InteractionDynamicsCK,
            DiffusionRelaxationCK<Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>>,
            DiffusionRelaxationCK<Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>>>(
            diffusion_body_inner, &isotropic_diffusion_explict);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_pair_wise);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_projection);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_explicit);

    auto &observe_temperature_pair_wise = main_methods.addIODynamics<
        ObservedQuantityRecording, Real>(phi_pair_wise, observer_contact);
    auto &observe_temperature_projection = main_methods.addIODynamics<
        ObservedQuantityRecording, Real>(phi_projection, observer_contact);
    auto &observe_temperature_explicit = main_methods.addIODynamics<
        ObservedQuantityRecording, Real>(phi_explicit, observer_contact);
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

    initial_condition_pair_wise.exec();
    initial_condition_projection.exec();
    initial_condition_explicit.exec();
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount time_instance;
    TimeInterval interval_pair_wise;
    TimeInterval interval_projection;
    TimeInterval interval_explicit;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    observe_temperature_pair_wise.writeToFile(time_steps);
    observe_temperature_projection.writeToFile(time_steps);
    observe_temperature_explicit.writeToFile(time_steps);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        Real diffusion_dt = time_stepper.incrementPhysicalTime(1.0 / 16.0);

        time_instance = TickCount::now();
        diffusion_relaxation_pair_wise.exec(diffusion_dt);
        interval_pair_wise += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        time_stepper.integrateMatchedTimeInterval(diffusion_relaxation_projection, diffusion_dt);
        interval_projection += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        time_stepper.integrateMatchedTimeInterval(diffusion_relaxation_explicit, diffusion_dt, get_time_step_size);
        interval_explicit += TickCount::now() - time_instance;
        time_steps += 1;

        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (state_recording())
        {
            body_state_recorder.writeToFile();
        }

        if (observation_recording())
        {
            observe_temperature_pair_wise.writeToFile(time_steps);
            observe_temperature_projection.writeToFile(time_steps);
            observe_temperature_explicit.writeToFile(time_steps);
        }

        if (time_steps % 1 == 0)
        {
            std::cout << "N=" << time_steps << " Time: "
                      << time_stepper.getPhysicalTime() << "	dt: " << diffusion_dt << "\n";
        }
    }

    TimeInterval tt = interval_pair_wise + interval_projection + interval_explicit;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Pair-wise diffusion time: " << interval_pair_wise.seconds() << " seconds." << std::endl;
    std::cout << "Projection diffusion time: " << interval_projection.seconds() << " seconds." << std::endl;
    std::cout << "Explicit diffusion time: " << interval_explicit.seconds() << " seconds." << std::endl;

    approximated_pair_wise = *observe_temperature_pair_wise.getObservedQuantity();
    approximated_projection = *observe_temperature_projection.getObservedQuantity();
    approximated_explicit = *observe_temperature_explicit.getObservedQuantity();

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
