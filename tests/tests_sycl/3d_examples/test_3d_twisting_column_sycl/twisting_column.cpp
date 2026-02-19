/**
 * @file twisting_column.cpp
 * @brief This is an example of solid with classic Neohookean model
 * to demonstrate the robustness of the formulation with Kirchhoff stress decomposition.
 * @author Chi Zhang and Xiangyu Hu
 * @ref DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real total_physical_time = 0.5; /**< TOTAL SIMULATION TIME*/
Real PL = 6.0;                  /**< X-direction size. */
Real PH = 1.0;                  /**< Y-direction size. */
Real PW = 1.0;                  /**< Z-direction size. */
Real particle_spacing_ref = PH / 10.0;
/** You can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real BW = particle_spacing_ref * 0.0; /**< no wall boundary in this case. */
Real SL = particle_spacing_ref * 1.0; /**< Length of the holder is one layer particle. */
Vecd halfsize_column(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_column(0.5 * (PL - SL), 0.0, 0.0);
Vecd halfsize_holder(0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
Vecd translation_holder(-0.5 * (SL + BW), 0.0, 0.0);
Vec3d domain_lower_bound(-SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
Vec3d domain_upper_bound(PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
StdVec<Vecd> observation_location = {Vecd(PL, 0.0, 0.0)};
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 1100.0; /**< Reference density. */
Real poisson = 0.45;  /**< Poisson ratio. */
Real Youngs_modulus = 1.7e7;
Real angular_0 = -300.0;
//------------------------------------------------------------------------------
// define a velocity profile for initial condition
//------------------------------------------------------------------------------
class VelocityProfile : public ReturnFunction<Vecd>
{
    Real angular_0_ = angular_0;
    Real PL_ = PL;

  public:
    VelocityProfile() {};

    Vecd operator()(const Vec3d &position)
    {
        Real x = position[0];
        Real y = position[1];
        Real z = position[2];
        Real angular_velocity = angular_0_ * sin((M_PI * x) / (2.0 * PL_));
        Real local_radius = sqrt(pow(y, 2) + pow(z, 2));
        Real angular = atan2(y, z);

        if (x > 0.0)
        {
            return Vecd(0.0, angular_velocity * local_radius * cos(angular),
                        -angular_velocity * local_radius * sin(angular));
        }
        else
        {
            return Vecd::Zero();
        }
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    // Build up an SPHSystem and IO environment.
    // Please the make sure the global domain bounds are correctly defined.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Setup geometry first.
    //----------------------------------------------------------------------
    auto &column_shape = sph_system.addShape<ComplexShape>("Column");
    column_shape.add<GeometricShapeBox>(Transform(translation_column), halfsize_column);
    auto &holder_shape = sph_system.addShape<GeometricShapeBox>(Transform(translation_holder), halfsize_holder, "Holder");
    column_shape.add<GeometricShapeBox>(holder_shape);
    //----------------------------------------------------------------------
    // Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    auto &column = sph_system.addBody<SolidBody>(column_shape);
    column.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    column.generateParticles<BaseParticles, Lattice>();
    BodyRegionByParticle holder(column, holder_shape);

    auto &my_observer = sph_system.addBody<ObserverBody>("MyObserver");
    my_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    auto &column_inner = sph_system.addInnerRelation(column, ConfigType::Lagrangian);
    auto &my_observer_contact = sph_system.addContactRelation(my_observer, column, ConfigType::Lagrangian);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
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
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<VelocityProfile>>(column, "Velocity").exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(column));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(column_inner));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(my_observer_contact));

    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);
    auto &column_acoustic_step_1st_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, LinearCorrectionCK>(column_inner);
    auto &column_acoustic_step_2nd_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(column_inner);

    auto &column_acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(column, 0.5);
    auto &constraint_holder = main_methods.addStateDynamics<FixBodyPartConstraintCK>(holder);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_displacement = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>("Position", my_observer_contact);
    auto &write_velocity = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>("Velocity", my_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    size_t acoustic_steps = 1;
    int screening_interval = 100;
    int observation_interval = screening_interval;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 250.0);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    lagrangian_configuration.exec();
    column_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    write_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_acoustic_step;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    write_displacement.writeToFile(acoustic_steps);
    write_velocity.writeToFile(acoustic_steps);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(column_acoustic_time_step);
        column_acoustic_step_1st_half.exec(acoustic_dt);
        constraint_holder.exec();
        column_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        /** Output body state during the simulation according output_interval. */
        time_instance = TickCount::now();
        /** screen output, write body observables and restart files  */
        if (acoustic_steps == 1 || acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (acoustic_steps % observation_interval == 0)
        {
            write_displacement.writeToFile(acoustic_steps);
            write_velocity.writeToFile(acoustic_steps);
        }

        if (state_recording())
        {
            body_state_recorder.writeToFile();
        }
        interval_output += TickCount::now() - time_instance;

        acoustic_steps++;
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
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
