/**
 * @file 	two_phase_dambreak_ck.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	sphinxsim@einsimo.com
 */
#include "near_wall_boundary.hpp"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 20.0;    /**< TOTAL SIMULATION TIME*/
Real DL = 5.3;                      /**< Tank length. */
Real DH = 2.0;                      /**< Tank height. */
Real LL = 2.0;                      /**< Liquid column length. */
Real LH = 1.0;                      /**< Liquid column height. */
Real particle_spacing_ref = 0.05;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// Observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of water. */
Real rho0_a = 0.001;                     /**< Reference density of air. */
Real gravity_g = 1.0;                    /**< Gravity force of fluid. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
namespace SPH
{
namespace fluid_dynamics
{

template <typename...>
class BuoyancyForce;

template <typename... Parameters>
class BuoyancyForce<Contact<WithUpdate, Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using BaseDynamicsType = Interaction<Contact<Parameters...>>;

  public:
    BuoyancyForce(Contact<Parameters...> &water_air_relation, const Gravity &gravity, Real threshold = 0.3)
        : BaseDynamicsType(water_air_relation), ForcePriorCK(this->particles_, "BuoyancyForce"),
          gravity_(gravity), threshold_by_dimensions_(threshold * Dimensions),
          dv_mass(this->particles_->template getVariableByName<Real>("Mass")),
          dv_isolation_factor_(nullptr), // from contact body
          dv_buoyancy_force_(this->getCurrentForce())
    {
        dv_isolation_factor_ = this->contact_particles_->template getVariableByName<Real>("IsolationFactor");
    };
    virtual ~BuoyancyForce() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser),
              gravity_(encloser.gravity_),
              threshold_by_dimensions_(encloser.threshold_by_dimensions_),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              mass_(encloser.dv_mass->DelegatedData(ex_policy)),
              isolation_factor_(encloser.dv_isolation_factor_->DelegatedData(ex_policy)),
              buoyancy_force_(encloser.dv_buoyancy_force_->DelegatedData(ex_policy)) {};

        void interact(size_t index_i, Real dt = 0.0)
        {
            Vecd buoyancy_i = Vecd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                UnsignedInt index_j = this->neighbor_index_[n];
                Real W_ij_Vi = this->W_ij(index_i, index_j) * Vol_[index_i];

                // Reverse gravity when air particle isolation factor less than threshold
                // Additional mass is a assume for only fraction of exclused volume
                // so that the summation is only contributed on water particles
                // Also, the original gravity should be canceled
                Vecd buoyancy_j = -SMAX(1.0 - isolation_factor_[index_j] / threshold_by_dimensions_, 0.0) *
                                  2.0 * gravity_.InducedAcceleration(Vecd::Zero(), 0.0);
                buoyancy_i += W_ij_Vi * buoyancy_j;
            }

            buoyancy_force_[index_i] = buoyancy_i * mass_[index_i];
        };

      protected:
        Gravity gravity_;
        Real threshold_by_dimensions_;
        Real *Vol_, *mass_, *isolation_factor_;
        Vecd *buoyancy_force_;
    };

  protected:
    const Gravity gravity_;
    Real threshold_by_dimensions_;
    DiscreteVariable<Real> *dv_mass, *dv_isolation_factor_;
    DiscreteVariable<Vecd> *dv_buoyancy_force_;
};
} // namespace fluid_dynamics
} // namespace SPH
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Geometric shapes.
    //----------------------------------------------------------------------
    GeometricShapeBox water_body_shape(BoundingBoxd(Vecd(0.0, 0.0), Vecd(LL, LH)), "WaterBody");
    GeometricShapeBox inner_wall_shape(BoundingBoxd(Vecd(0.0, 0.0), Vecd(DL, DH)), "InnerWall");
    GeometricShapeBox outer_wall_shape(BoundingBoxd(Vecd(-BW, -BW), Vecd(DL + BW, DH + BW)), "OuterWall");

    ComplexShape air_body_shape("AirBody");
    air_body_shape.add(&inner_wall_shape);
    air_body_shape.subtract(&water_body_shape);

    ComplexShape wall_boundary_shape("WallBoundary");
    wall_boundary_shape.add(&outer_wall_shape);
    wall_boundary_shape.subtract(&inner_wall_shape);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, water_body_shape);
    water_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_body.generateParticles<BaseParticles, Lattice>();

    FluidBody air_body(sph_system, air_body_shape);
    air_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_a, c_f);
    air_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, wall_boundary_shape);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody water_observer(sph_system, "FluidObserver");
    water_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &water_inner = sph_system.addInnerRelation(water_body);
    auto &water_air_contact = sph_system.addContactRelation(water_body, air_body);
    auto &water_wall_contact = sph_system.addContactRelation(water_body, wall_boundary);
    auto &air_inner = sph_system.addInnerRelation(air_body);
    auto &air_water_contact = sph_system.addContactRelation(air_body, water_body);
    auto &air_wall_contact = sph_system.addContactRelation(air_body, wall_boundary);
    auto &water_observer_contact = sph_system.addContactRelation(water_observer, water_body);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    auto &host_methods = sph_solver.getHostMethodContainer();
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
    host_methods.addStateDynamics<NormalFromSubShapeAndOpCK>(wall_boundary, "InnerWall").exec(); // host method able to run immediately

    StdVec<FluidBody *> fluid_bodies = {&water_body, &air_body};

    ParticleDynamicsGroup update_fluid_cell_linked_list = main_methods.addCellLinkedListDynamics(fluid_bodies);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    ParticleDynamicsGroup update_fluid_relation;
    update_fluid_relation.add(&main_methods.addRelationDynamics(water_inner, water_air_contact, water_wall_contact));
    update_fluid_relation.add(&main_methods.addRelationDynamics(air_inner, air_water_contact, air_wall_contact));
    ParticleDynamicsGroup update_fluid_configuration = update_fluid_cell_linked_list + update_fluid_relation;
    auto &water_observer_contact_relation = main_methods.addRelationDynamics(water_observer_contact);
    ParticleDynamicsGroup fluid_particle_sorting = main_methods.addSortDynamics(fluid_bodies);

    Gravity gravity(Vecd(0.0, -gravity_g));
    ParticleDynamicsGroup constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(fluid_bodies, gravity);
    ParticleDynamicsGroup advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(fluid_bodies);
    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(fluid_bodies);

    ParticleDynamicsGroup fluid_linear_correction_matrix;
    fluid_linear_correction_matrix.add(
        &main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(water_inner)
             .addPostContactInteraction(water_air_contact)
             .addPostContactInteraction(water_wall_contact));
    fluid_linear_correction_matrix.add(
        &main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(air_inner)
             .addPostContactInteraction(air_water_contact)
             .addPostContactInteraction(air_wall_contact));

    ParticleDynamicsGroup acoustic_step_1st_half;
    acoustic_step_1st_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact));
    acoustic_step_1st_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(air_wall_contact));

    ParticleDynamicsGroup acoustic_step_2nd_half;
    acoustic_step_2nd_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact));
    acoustic_step_2nd_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(air_wall_contact));

    ParticleDynamicsGroup density_regularization;
    density_regularization.add(
        &main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(water_inner)
             .addPostContactInteraction(water_wall_contact)
             .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeSurface>(water_body));
    density_regularization.add(
        &main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(air_inner)
             .addPostContactInteraction(air_water_contact)
             .addPostContactInteraction(air_wall_contact)
             .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, Internal>(air_body));

    ParticleDynamicsGroup transport_correction;
    transport_correction.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(air_wall_contact)
             .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter>(air_body));
    transport_correction.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(water_wall_contact)
             .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter>(water_body));

    auto &air_near_wall_bounding = main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::BoundingFromWallCK>(air_wall_contact);

    ParticleDynamicsGroup buoyancy_force;
    buoyancy_force.add(&main_methods.addInteractionDynamics<fluid_dynamics::IsolationFactor>(air_inner));
    buoyancy_force.add(&main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::BuoyancyForce>(water_air_contact, gravity, 0.3));

    ReduceDynamicsGroup advection_time_step = main_methods.addReduceDynamics<ReduceMin<Real>, fluid_dynamics::AdvectionTimeStepCK>(fluid_bodies, U_ref);
    ReduceDynamicsGroup acoustic_time_step = main_methods.addReduceDynamics<ReduceMin<Real>, fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(fluid_bodies);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(air_body, "NormalDirection");
    body_state_recorder.addToWrite<Real>(air_body, "SignedDistance");
    body_state_recorder.addToWrite<Real>(air_body, "IsolationFactor");
    body_state_recorder.addToWrite<Vecd>(water_body, "BuoyancyForce");
    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);
    auto &record_water_mechanical_energy = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, TotalMechanicalEnergyCK>(water_body, gravity);
    auto &water_observer_pressure = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Real>( water_observer_contact, "Pressure");
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        time_stepper.setPhysicalTime(restart_io.readRestartFiles(sph_system.RestartStep()));
    }
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep();
    int screening_interval = 100;
    int observation_interval = screening_interval * 2;
    int restart_output_interval = screening_interval * 10;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 100);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_cell_linked_list.exec();
    update_fluid_configuration.exec();
    water_observer_contact_relation.exec();

    constant_gravity.exec();
    density_regularization.exec();
    advection_step_setup.exec();
    fluid_linear_correction_matrix.exec();
    buoyancy_force.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    record_water_mechanical_energy.writeToFile(advection_steps);
    water_observer_pressure.writeToFile(advection_steps);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    TickCount t0 = TickCount::now();
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime(total_physical_time))
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(advection_time_step))
        {
            advection_steps++;
            update_particle_position.exec();
            air_near_wall_bounding.exec();

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
                record_water_mechanical_energy.writeToFile(advection_steps);
                water_observer_contact_relation.exec();
                water_observer_pressure.writeToFile(advection_steps);
            }

            if (advection_steps % restart_output_interval == 0)
            {
                restart_io.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (advection_steps % 100)
            {
                fluid_particle_sorting.exec();
            }
            update_fluid_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            density_regularization.exec();
            advection_step_setup.exec();
            transport_correction.exec();
            fluid_linear_correction_matrix.exec();
            buoyancy_force.exec();
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
    //----------------------------------------------------------------------
    // Post-run regression test to ensure that the case is validated
    //----------------------------------------------------------------------
    if (sph_system.GenerateRegressionData())
    {
        record_water_mechanical_energy.generateDataBase(1.0e-3);
        water_observer_pressure.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        record_water_mechanical_energy.testResult();
        water_observer_pressure.testResult();
    }

    return 0;
}
