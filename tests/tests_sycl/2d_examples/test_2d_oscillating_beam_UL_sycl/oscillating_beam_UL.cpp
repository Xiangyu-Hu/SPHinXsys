/* ----------------------------------------------------------------------------*
 *                    SPHinXsys: 2D oscillation beam                           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases for understanding SPH method for    *
 * solid simulation based on updated Lagrangian method.                        *
 * A generalized hourglass control method is used here.                         *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu                              *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real total_physical_time = 1.0; /**< TOTAL SIMULATION TIME*/
Real PL = 0.2;                  // beam length
Real PH = 0.02;                 // for thick plate
Real SL = 0.06;                 // depth of the insert
Real global_resolution = PH / 10;
Real BW = global_resolution * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
                                  Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;       // Poisson ratio
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
Real U_ref = vf * c0 * (M * (cos(kl) - cosh(kl)) - N * (sin(kl) - sinh(kl))) / Q;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
GeometricShapeBox beam_base_shape(
    BoundingBoxd(Vec2d(-SL - BW, -PH / 2 - BW), Vec2d(0.0, PH / 2 + BW)), "BeamBase");
GeometricShapeBox beam_column(
    BoundingBoxd(Vec2d(-SL, -PH / 2), Vec2d(PL, PH / 2)), "BeamColumn");
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//------------------------------------------------------------------------------
// define a linear profile for initial condition
//------------------------------------------------------------------------------
class LinearProfile : public ReturnFunction<Vecd>
{
    Real vf_ = vf;
    Real c0_ = c0;
    Real kl_ = kl;
    Real M_ = M;
    Real N_ = N;
    Real Q_ = Q;
    Real PL_ = PL;

  public:
    LinearProfile() {};

    Vecd operator()(const Vec2d &position)
    {
        Real x = position[0] / PL_;
        Vecd result = Vec2d::Zero();
        if (x > 0.0)
        {
            result[1] = vf_ * c0_ * (M_ * (cos(kl_ * x) - cosh(kl_ * x)) - N_ * (sin(kl_ * x) - sinh(kl_ * x))) / Q_;
        }
        return result;
    }
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &beam_shape = sph_system.addShape<ComplexShape>("BeamBody");
    beam_shape.add(&beam_base_shape);
    beam_shape.add(&beam_column);
    auto &beam = sph_system.addBody<RealBody>(beam_shape);
    beam.defineMaterial<GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
    beam.generateParticles<BaseParticles, Lattice>();
    BodyRegionByParticle beam_base(beam, beam_base_shape);

    ObserverBody beam_observer(sph_system, "BeamObserver");
    beam_observer.getSPHAdaptation().resetAdaptationRatios(1.15, 2.0);
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    /**body relation topology */
    auto &beam_inner = sph_system.addInnerRelation(beam);
    auto &beam_observer_contact = sph_system.addContactRelation(beam_observer, beam);
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
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<LinearProfile>>(beam, "Velocity").exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup update_beam_configuration;
    update_beam_configuration.add(&main_methods.addCellLinkedListDynamics(beam));
    update_beam_configuration.add(&main_methods.addRelationDynamics(beam_inner));
    auto &beam_observer_contact_relation = main_methods.addRelationDynamics(beam_observer_contact);

    auto &beam_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(beam);
    auto &beam_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(beam);
    auto &beam_base_constraint = main_methods.addStateDynamics<ConstantConstraintCK, Vecd>(beam_base, "Velocity", Vec2d::Zero());
    auto &beam_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(beam_inner);

    ParticleDynamicsGroup column_shear_force;
    column_shear_force.add(&main_methods.addInteractionDynamics<LinearGradient, Vecd>(beam_inner, "Velocity"));
    column_shear_force.add(&main_methods.addInteractionDynamicsOneLevel<continuum_dynamics::ShearIntegration, GeneralContinuum>(beam_inner));

    auto &beam_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep1stHalf, NoRiemannSolverCK, NoKernelCorrectionCK>(beam_inner);
    auto &beam_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf, DissipativeRiemannSolverCK, NoKernelCorrectionCK>(beam_inner);

    auto &beam_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(beam, U_ref, 0.2);
    auto &beam_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(beam, 0.4);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(beam, "Pressure");
    body_state_recorder.addToWrite<Real>(beam, "Density");
    body_state_recorder.addDerivedVariableToWrite<continuum_dynamics::VonMisesStressCK>(beam);
    auto &record_beam_mechanical_energy = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, TotalKineticEnergyCK>(beam);
    auto &beam_observer_position = main_methods.addObserveRecorder<Vecd>("Position", beam_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(beam_advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep() + 1;
    int screening_interval = 100;
    int observation_interval = screening_interval * 2;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 100.0);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    update_beam_configuration.exec();

    beam_advection_step_setup.exec();
    beam_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    record_beam_mechanical_energy.writeToFile(advection_steps);
    beam_observer_position.writeToFile(advection_steps);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
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
        Real acoustic_dt = time_stepper.incrementPhysicalTime(beam_acoustic_time_step);
        column_shear_force.exec(acoustic_dt);
        beam_acoustic_step_1st_half.exec(acoustic_dt);
        beam_base_constraint.exec();
        beam_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(beam_advection_time_step))
        {
            advection_steps++;
            beam_update_particle_position.exec();

            /** Output body state during the simulation according output_interval. */
            time_instance = TickCount::now();
            /** screen output, write body observables and restart files  */
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                record_beam_mechanical_energy.writeToFile(advection_steps);
                beam_observer_contact_relation.exec();
                beam_observer_position.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            update_beam_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            beam_advection_step_setup.exec();
            beam_linear_correction_matrix.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        record_beam_mechanical_energy.generateDataBase(1.0e-3);
    }
    else
    {
        record_beam_mechanical_energy.testResult();
    }
    return 0;
}