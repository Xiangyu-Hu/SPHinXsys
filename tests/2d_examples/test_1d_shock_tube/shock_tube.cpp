/**
 * @file shock_tube.cpp
 * @brief This is a test to show the standard Sod shock tube case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and simulation setup.
//----------------------------------------------------------------------
Real DL = 5.0;                           /**< Tube length. */
Real particle_spacing_ref = 1.0 / 200.0; /**< Initial reference particle spacing. */
Real DH = particle_spacing_ref * 4;      /**< Tube height. */
BoundingBoxd system_domain_bounds(Vec2d(-2.0 / 5.0 * DL, 0.0), Vec2d(3.0 / 5.0 * DL, DH));
Real rho0_l = 1.0;              /**< initial density of left state. */
Real rho0_r = 0.125;            /**< initial density of right state. */
Vecd velocity_l = Vecd::Zero(); /**< initial velocity of left state. */
Vecd velocity_r = Vecd::Zero(); /**< initial velocity of right state. */
Real p_l = 1.0;                 /**< initial pressure of left state. */
Real p_r = 0.1;                 /**< initial pressure of right state. */
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real heat_capacity_ratio = 1.4; /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Cases-dependent geometry
//----------------------------------------------------------------------
class WaveBlock : public MultiPolygonShape
{
  public:
    explicit WaveBlock(const std::string &body_name)
        : MultiPolygonShape(body_name)
    {
        std::vector<Vecd> waves_block_shape{
            Vecd(-2.0 / 5.0 * DL, 0.0), Vecd(-2.0 / 5.0 * DL, DH), Vecd(3.0 / 5.0 * DL, DH),
            Vecd(3.0 / 5.0 * DL, 0.0), Vecd(-2.0 / 5.0 * DL, 0.0)};
        multi_polygon_.addAPolygon(waves_block_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class ShockTubeInitialCondition : public fluid_dynamics::CompressibleFluidInitialCondition
{
  public:
    explicit ShockTubeInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::CompressibleFluidInitialCondition(sph_body) {};
    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][0] < DL / 10.0)
        {
            // initial left state pressure,momentum and energy profile
            rho_[index_i] = rho0_l;
            mass_[index_i] = rho_[index_i] * Vol_[index_i];
            p_[index_i] = p_l;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i] = velocity_l;
            mom_[index_i] = mass_[index_i] * vel_[index_i];
            E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
        }
        if (pos_[index_i][0] > DL / 10.0)
        {
            // initial right state pressure,momentum and energy profile
            rho_[index_i] = rho0_r;
            mass_[index_i] = rho_[index_i] * Vol_[index_i];
            p_[index_i] = p_r;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i] = velocity_r;
            mom_[index_i] = mass_[index_i] * vel_[index_i];
            E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
        }
    }

  protected:
    Real gamma_ = heat_capacity_ratio;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Create body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody wave_body(sph_system, makeShared<WaveBlock>("WaveBody"));
    wave_body.defineMaterial<CompressibleFluid>(rho0_l, heat_capacity_ratio);
    wave_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The inner relation defines the particle configuration for particles within a body.
    //	The contact relation defines the particle configuration between the bodies.
    //----------------------------------------------------------------------
    InnerRelation wave_body_inner(wave_body);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfHLLCRiemann> pressure_relaxation(wave_body_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfHLLCRiemann> density_and_energy_relaxation(wave_body_inner);

    SimpleDynamics<ShockTubeInitialCondition> waves_initial_condition(wave_body);
    PeriodicAlongAxis periodic_along_y(wave_body.getSPHBodyBounds(), yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(wave_body, periodic_along_y);
    ReduceDynamics<fluid_dynamics::EulerianCompressibleAcousticTimeStepSize> get_wave_time_step_size(wave_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> kernel_correction_matrix(wave_body_inner);
    InteractionDynamics<KernelGradientCorrectionInner> kernel_gradient_update(wave_body_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations of the simulation.
    //	Regression tests are also defined here.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(wave_body, "TotalEnergy");
    body_states_recording.addToWrite<Real>(wave_body, "Density");
    RegressionTestEnsembleAverage<ReducedQuantityRecording<MaximumSpeed>>
        write_maximum_speed(wave_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    waves_initial_condition.exec();
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.2;
    Real output_interval = 0.01; /**< Time stamps for output of body states. */
    //----------------------------------------------------------------------
    // Output the start states of bodies.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Statistics for computing CPU time.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        //	Integrate time (loop) until the next output time.
        while (integration_time < output_interval)
        {
            Real dt = get_wave_time_step_size.exec();
            // Dynamics including pressure and density and energy relaxation.
            integration_time += dt;
            pressure_relaxation.exec(dt);
            density_and_energy_relaxation.exec(dt);
            physical_time += dt;

            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_maximum_speed.generateDataBase(1.0e-3, 1.0e-3);
    }
    else
    {
        write_maximum_speed.testResult();
    }

    return 0;
}
