/**
 * @file	rotation_patch.cpp
 * @brief	2D rotation patch example.
 * @author	Yaru Ren, Chi Zhang and Xiangyu Hu
 */

#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;                        /**< Liquid column length. */
Real LH = 1.0;                        /**< Liquid column height. */
Real particle_spacing_ref = LL / 100; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;   /**< Extending width for boundary conditions. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(LL + BW, LH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of fluid. */
Real U_max = 1.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */

Vec2d DamP_lb(-LL / 2, -LH / 2); /**< Left bottom. */
Vec2d DamP_lt(-LL / 2, LH / 2);  /**< Left top. */
Vec2d DamP_rt(LL / 2, LH / 2);   /**< Right top. */
Vec2d DamP_rb(LL / 2, -LH / 2);  /**< Right bottom. */

//----------------------------------------------------------------------
//	Geometry definition.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(DamP_lb);
    water_block_shape.push_back(DamP_lt);
    water_block_shape.push_back(DamP_rt);
    water_block_shape.push_back(DamP_rb);
    water_block_shape.push_back(DamP_lb);

    return water_block_shape;
}
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

/**
 * application dependent initial velocity
 */
class InitialVelocity
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    InitialVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          fluid_particles_(dynamic_cast<BaseParticles *>(&sph_body.getBaseParticles())),
          p_(fluid_particles_->getVariableDataByName<Real>("Pressure")),
          rho_(fluid_particles_->getVariableDataByName<Real>("Density")) {};

    void update(size_t index_i, Real dt)
    {
        Real omega = 1.0;
        /** initial velocity profile */
        vel_[index_i][0] = omega * pos_[index_i][1];
        vel_[index_i][1] = -omega * pos_[index_i][0];
        /** initial pressure */
        for (size_t m = 0; m != 100; ++m)
        {
            for (size_t n = 0; n != 100; ++n)
            {
                if (m % 2 == 1 && n % 2 == 1)
                {
                    Real x_star = pos_[index_i][0] + LL / 2;
                    Real y_star = pos_[index_i][1] + LL / 2;
                    Real coefficient1 = m * n * PI * PI * (pow((m * PI / LL), 2) + pow((n * PI / LL), 2));
                    p_[index_i] += rho0_f * (-32 * omega * omega) / coefficient1 * sin(m * PI * x_star / LL) * sin(n * PI * y_star / LL);
                }
            }
        }
        rho_[index_i] = p_[index_i] / pow(c_f, 2) + rho0_f;
    }

  protected:
    BaseParticles *fluid_particles_;
    Real *p_, *rho_;
};
//----------------------------------------------------------------------
//	wave gauge
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-LL, -LL), Vec2d(LL + 10 * BW, LH + 10 * BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    // Using relaxed particle distribution if needed
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(0.0, 0.0)};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationInner> free_surface_indicator(water_body_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_fluid(water_body_inner, 0.5);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionInnerRiemann> fluid_pressure_relaxation_correct(water_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> fluid_density_relaxation(water_body_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationInnerFreeStream> update_density_by_summation(water_body_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<NoLimiter, BulkParticles>> transport_velocity_correction(water_body_inner);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);
    SimpleDynamics<InitialVelocity> initial_condition(water_block);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "DensitySummation");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Matd>(water_block, "LinearGradientCorrectionMatrix");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    ObservedQuantityRecording<Real> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    if (sph_system.RestartStep() != 0)
    {
        physical_time = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_body_inner.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 8.0;
    Real output_interval = 0.04;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First state recording before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();
            free_surface_indicator.exec();
            corrected_configuration_fluid.exec();
            update_density_by_summation.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                Real acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_pressure_relaxation_correct.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                inner_ite_dt++;
                physical_time += acoustic_dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";

                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);

                write_water_kinetic_energy.writeToFile(number_of_iterations);
                write_recorded_water_pressure.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_body_inner.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        body_states_recording.writeToFile();

        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_water_kinetic_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_water_kinetic_energy.testResult();
    }

    return 0;
};
