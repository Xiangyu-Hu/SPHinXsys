/**
 * @file diffusion_steady.cpp
 * @brief This is the example to test the solution of steady dissipation system
 * uing pseudo-time integration.
 * @author Xiangyu Hu
 */
#include "sphinxsys_ck.h" //SPHinXsys Library
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
std::string phi_implicit = "PhiImplicit";
const Real diffusion_coeff = 1.0;
//----------------------------------------------------------------------
//	Parameters for initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 300.0;
Real zero_residue = 0.0;
const Real high_temperature = 300.0;
const Real low_temperature = 400.0;
Real heat_source = 1000.0;
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
            if (position[1] < 0.0)
            {
                return low_temperature;
            }
            if (position[1] > 1.0)
            {
                return high_temperature;
            }
            return 0.0;
        };
    };
};
//----------------------------------------------------------------------
//	Google test items
//----------------------------------------------------------------------
Real expected_value = 587.88;

Real approximated_implicit(1.0);
TEST(DiffusionPairWise, Error)
{
    EXPECT_LT(ABS(expected_value - approximated_implicit), 1.0e-2);
    std::cout << "Reference Value: " << expected_value << " and "
              << "Predicted Value: " << approximated_implicit << std::endl;
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
    auto &wall_boundary_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &diffusion_body_update_relation = main_methods.addRelationDynamics(diffusion_body_inner, diffusion_body_contact);

    auto &initial_condition_implicit =
        main_methods.addStateDynamics<InitialCondition<SPHBody, UniformDistribution<Real>>>(diffusion_body, phi_implicit, initial_temperature);
    auto &boundary_condition_implicit =
        main_methods.addStateDynamics<InitialCondition<SPHBody, InitialDistribution>>(wall_boundary, phi_implicit);
    auto &zero_residue_boundary_condition =
        main_methods.addStateDynamics<InitialCondition<SPHBody, UniformDistribution<Real>>>(wall_boundary, "Residue" + phi_implicit, zero_residue);
    using MainExecutionPolicy = execution::ParallelPolicy; // define execution policy for this case
    ImplicitDissipation<MainExecutionPolicy, Inner<IsotropicDiffusion>>
        diffusion_relaxation_implicit(diffusion_body_inner, phi_implicit, 1.0e-6);
    diffusion_relaxation_implicit.addPostContactInteraction<Dirichlet<IsotropicDiffusion>>(diffusion_body_contact, phi_implicit);
    auto &reduce_average_species = main_methods.addReduceDynamics<QuantityAverage<Real>>(diffusion_body, phi_implicit);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(diffusion_body, phi_implicit);
    body_state_recorder.addToWrite<Real>(wall_boundary, phi_implicit);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(1.0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t time_steps = 0;
    Real output_peroid = 0.1;
    auto &state_recording = time_stepper.addTriggerByInterval(output_peroid);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    diffusion_body_cell_linked_list.exec();
    wall_boundary_cell_linked_list.exec();
    diffusion_body_update_relation.exec();

    initial_condition_implicit.exec();
    boundary_condition_implicit.exec();
    zero_residue_boundary_condition.exec();
    Real initial_average_species = reduce_average_species.exec();
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount time_instance;
    TimeInterval interval_implicit;
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
        time_stepper.integrateMatchedTimeInterval(diffusion_relaxation_implicit, diffusion_dt);
        interval_implicit += TickCount::now() - time_instance;
        time_steps += 1;

        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (state_recording())
        {
            body_state_recorder.writeToFile();
        }

        if (time_steps % 1 == 0)
        {
            std::cout << "N=" << time_steps << " Time: "
                      << time_stepper.getPhysicalTime() << "	dt: " << diffusion_dt << "\n";
            std::cout << "Initial average species: " << initial_average_species
                      << "	Present average species: " << reduce_average_species.exec() << "\n";
        }
    }

    TimeInterval tt = interval_implicit;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Pair-wise diffusion time: " << interval_implicit.seconds() << " seconds." << std::endl;

    approximated_implicit = reduce_average_species.exec();

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
