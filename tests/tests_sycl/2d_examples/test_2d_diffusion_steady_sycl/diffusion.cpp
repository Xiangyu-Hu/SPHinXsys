/**
 * @file diffusion_steady.cpp
 * @brief This is the example to test the solution of steady dissipation system
 * uing pseudo-time integration.
 * @author Xiangyu Hu
 */
#include "sphinxsys_sycl.h" //SPHinXsys Library
#include <gtest/gtest.h>
using namespace SPH; // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 50.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
std::string phi_pairwise = "PhiPairWise";
std::string phi_implicit = "PhiImplicit";
std::string phi_explicit = "PhiExplicit";
const Real diffusion_coeff = 1.0;
//----------------------------------------------------------------------
//	Parameters for initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 300.0;
Real high_temperature = 300.0;
Real low_temperature = 400.0;
Real heat_source = 1000.0;
//----------------------------------------------------------------------
// Define extra classes which are used in the main program.
// These classes are defined under the namespace of SPH.
//----------------------------------------------------------------------
class InitialProfile : public ReturnFunction<Real>
{
    Real low_temperature_;
    Real high_temperature_;

  public:
    InitialProfile()
        : low_temperature_(low_temperature),
          high_temperature_(high_temperature) {};

    Real operator()(const Vecd &position)
    {
        if (position[1] < 0.0)
        {
            return low_temperature_;
        }
        if (position[1] > 1.0)
        {
            return high_temperature_;
        }
        return 0.0;
    };
};
//----------------------------------------------------------------------
//	Google test items
//----------------------------------------------------------------------
Real expected_value = 350.0;

Real approximated_pairwise(1.0);
TEST(DiffusionPairWise, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_pairwise) / expected_value, 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_pairwise << std::endl;
};

Real approximated_implicit(1.0);
TEST(DiffusionImplicit, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_implicit) / expected_value, 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_implicit << std::endl;
};

Real approximated_explicit(1.0);
TEST(DiffusionExplicit, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_explicit) / expected_value, 1.0e-2);
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
        Solid(), ConstructArgs(phi_implicit, diffusion_coeff));
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    ComplexShape wall_complex_shape("WallBoundary");
    wall_complex_shape.add<GeometricShapeBox>(BoundingBox(Vecd(0.4 * L, -BW), Vecd(0.6 * L, 0.0)));
    wall_complex_shape.add<GeometricShapeBox>(BoundingBox(Vecd(0.4 * L, H), Vecd(0.6 * L, H + BW)));
    SolidBody wall_boundary(sph_system, wall_complex_shape);
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> diffusion_body_contact(diffusion_body, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_device);
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
    auto &wall_boundary_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &diffusion_body_update_relation = main_methods.addRelationDynamics(diffusion_body_inner, diffusion_body_contact);

    auto &initial_condition_pairwise =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, ConstantValue<Real>>>(diffusion_body, phi_pairwise, initial_temperature);
    auto &boundary_condition_pairwise =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, SpatialDistribution<InitialProfile>>>(wall_boundary, phi_pairwise);
    auto &diffusion_relaxation_pairwise =
        main_methods.addInteractionDynamics<ProjectionDissipation, Splitting, IsotropicDiffusion>(diffusion_body_inner, phi_pairwise)
            .addPreContactInteraction<Splitting, Dirichlet<IsotropicDiffusion>>(diffusion_body_contact, phi_pairwise);
    auto &pairwise_average_species = main_methods.addReduceDynamics<QuantityAverage<Real>>(diffusion_body, phi_pairwise);

    auto &initial_condition_implicit =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, ConstantValue<Real>>>(diffusion_body, phi_implicit, initial_temperature);
    auto &boundary_condition_implicit =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, SpatialDistribution<InitialProfile>>>(wall_boundary, phi_implicit);
    using MainExecutionPolicy = execution::ParallelDevicePolicy; // define execution policy for this case
    ImplicitDissipation<MainExecutionPolicy, Inner<IsotropicDiffusion>>
        diffusion_relaxation_implicit(diffusion_body_inner, phi_implicit, 1.0e-6);
    diffusion_relaxation_implicit.addPostContactInteraction<Dirichlet<IsotropicDiffusion>>(diffusion_body_contact, phi_implicit);
    auto &implicit_average_species = main_methods.addReduceDynamics<QuantityAverage<Real>>(diffusion_body, phi_implicit);

    IsotropicDiffusion isotropic_diffusion_explict(phi_explicit, diffusion_coeff);
    auto &initial_condition_explicit =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, ConstantValue<Real>>>(diffusion_body, phi_explicit, initial_temperature);
    auto &boundary_condition_explicit =
        main_methods.addStateDynamics<VariableAssignment<SPHBody, SpatialDistribution<InitialProfile>>>(wall_boundary, phi_explicit);
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, &isotropic_diffusion_explict);
    auto &diffusion_relaxation_explicit =
        main_methods.addRungeKuttaSequence<
            InteractionDynamicsCK,
            DiffusionRelaxationCK<Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                                  Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>,
            DiffusionRelaxationCK<Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                                  Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>>(
            DynamicsArgs(diffusion_body_inner, &isotropic_diffusion_explict),
            DynamicsArgs(diffusion_body_contact, &isotropic_diffusion_explict));
    auto &explicit_average_species = main_methods.addReduceDynamics<QuantityAverage<Real>>(diffusion_body, phi_explicit);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_implicit);
    body_state_recorder.addToWrite<Real>(wall_boundary, phi_pairwise);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_pairwise);
    body_state_recorder.addToWrite<Real>(wall_boundary, phi_implicit);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_explicit);
    body_state_recorder.addToWrite<Real>(wall_boundary, phi_explicit);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(5.0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t time_steps = 0;
    Real output_peroid = 0.5;
    auto &state_recording = time_stepper.addTriggerByInterval(output_peroid);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    diffusion_body_cell_linked_list.exec();
    wall_boundary_cell_linked_list.exec();
    diffusion_body_update_relation.exec();

    initial_condition_pairwise.exec();
    boundary_condition_pairwise.exec();
    Real initial_pairwise_average_species = pairwise_average_species.exec();

    initial_condition_implicit.exec();
    boundary_condition_implicit.exec();
    diffusion_relaxation_implicit.initializeImplicitDissipation();
    Real initial_implicit_average_species = implicit_average_species.exec();

    initial_condition_explicit.exec();
    boundary_condition_explicit.exec();
    Real initial_explicit_average_species = explicit_average_species.exec();
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount time_instance;
    TimeInterval interval_pairwise;
    TimeInterval interval_implicit;
    TimeInterval interval_explicit;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        Real diffusion_dt = time_stepper.incrementPhysicalTime(1.0 / 100.0);

        time_instance = TickCount::now();
        time_stepper.integrateMatchedTimeInterval(diffusion_relaxation_pairwise, diffusion_dt, 8);
        interval_pairwise += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        time_stepper.integrateMatchedTimeInterval(diffusion_relaxation_implicit, diffusion_dt);
        interval_implicit += TickCount::now() - time_instance;

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

        if (time_steps % 100 == 0)
        {
            std::cout << "N=" << time_steps << " Time: "
                      << time_stepper.getPhysicalTime() << "	dt: " << diffusion_dt << "\n";
            std::cout << "Initial pairwise average species: " << initial_pairwise_average_species
                      << "	Present pairwise average species: " << pairwise_average_species.exec() << "\n";
            std::cout << "Initial implicit average species: " << initial_implicit_average_species
                      << "	Present implicit average species: " << implicit_average_species.exec() << "\n";
            std::cout << "Initial explicit average species: " << initial_explicit_average_species
                      << "	Present explicit average species: " << explicit_average_species.exec() << "\n";
        }
    }

    TimeInterval tt = interval_pairwise + interval_implicit + interval_explicit;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "PairWise diffusion computing time: " << interval_pairwise.seconds() << " seconds." << std::endl;
    std::cout << "Implicit diffusion computing time: " << interval_implicit.seconds() << " seconds." << std::endl;
    std::cout << "Explicit diffusion computing time: " << interval_explicit.seconds() << " seconds." << std::endl;

    approximated_pairwise = pairwise_average_species.exec();
    approximated_implicit = implicit_average_species.exec();
    approximated_explicit = explicit_average_species.exec();

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
