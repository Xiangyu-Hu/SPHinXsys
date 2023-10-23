/**
 * @file 	eulerian_taylor_green.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_taylor_green vortex flow example.
 * @author 	Chi Zhang, Zhentong Wang and Xiangyu Hu
 */
#include "general_eulerian_fluid_dynamics.hpp" // eulerian classes for compressible fluid only.
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;                    /**< box length. */
Real DH = 1.0;                    /**< box height. */
Real resolution_ref = 1.0 / 25.0; /**< Global reference resolution. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                  /**< Reference density of fluid. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Reference sound speed. */
Real Re = 100;                      /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> water_body_shape;
        water_body_shape.push_back(Vecd(0.0, 0.0));
        water_body_shape.push_back(Vecd(0.0, DH));
        water_body_shape.push_back(Vecd(DL, DH));
        water_body_shape.push_back(Vecd(DL, 0.0));
        water_body_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class TaylorGreenInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit TaylorGreenInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
          rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
          mom_(*particles_->getVariableByName<Vecd>("Momentum")){};

    void update(size_t index_i, Real dt)
    {
        /** initial momentum and energy profile */
        rho_[index_i] = rho0_f;
        vel_[index_i][0] = -cos(2.0 * Pi * pos_[index_i][0]) *
            sin(2.0 * Pi * pos_[index_i][1]);
        vel_[index_i][1] = sin(2.0 * Pi * pos_[index_i][0]) *
            cos(2.0 * Pi * pos_[index_i][1]);
        mom_[index_i] = rho_[index_i] * vel_[index_i];
    }

protected:
    StdLargeVec<Vecd>& pos_, & vel_;
    StdLargeVec<Real>& rho_, & p_;
    StdLargeVec<Vecd>& mom_;

};
using namespace SPH; //	Namespace cite here.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); //Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);       //Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);// handle command line arguemnts.
#endif
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
        : water_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
       //	Define body relation map.
       //----------------------------------------------------------------------
        InnerRelation water_body_inner(water_body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_body);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_water_body_to_vtp(io_environment, { &water_body });
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, { &water_body });

        /* Relaxation method: including based on the 0th and 1st order consistency. */
        InteractionWithUpdate<KernelCorrectionMatrixInner> kernel_correction_inner(water_body_inner);
        relax_dynamics::RelaxationStepInnerImplicit relaxation_step_inner(water_body_inner);

        PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
        PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
        //----------------------------------------------------------------------  
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_water_body_particles.exec(0.15);
        sph_system.initializeSystemCellLinkedLists();
        periodic_condition_x.update_cell_linked_list_.exec();
        periodic_condition_y.update_cell_linked_list_.exec();
        sph_system.initializeSystemConfigurations();
        write_water_body_to_vtp.writeToFile(0);

        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        TickCount t1 = TickCount::now();

        int ite = 0; //iteration step for the total relaxation step.

        Real last_zero_maximum_residual = 1;
        Real last_zero_average_residual = 1;
        Real last_first_maximum_residual = 1;
        Real last_first_average_residual = 1;

        Real current_zero_maximum_residual = 1;
        Real current_zero_average_residual = 1;
        Real current_first_maximum_residaul = 1;
        Real current_first_average_residaul = 1;

        GlobalStaticVariables::physical_time_ = ite;
        /* The procedure to obtain uniform particle distribution that satisfies the 0ht order consistency. */
        while (current_zero_maximum_residual > 0.0001)
        {
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            water_body.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            water_body_inner.updateConfiguration();


            ite++;

            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
                std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << ": average = " << current_zero_average_residual << std::endl;
                write_water_body_to_vtp.writeToFile(ite);
            }
        }
        std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << ": average = " << current_zero_average_residual << std::endl;

        ite++;
        write_water_body_to_vtp.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);

        TickCount t2 = TickCount::now();
        TickCount::interval_t tt;
        tt = t2 - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_body);
    InnerRelation water_body_correct_inner(water_body);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<fluid_dynamics::ICEIntegration1stHalfHLLERiemann> pressure_relaxation(water_body_inner);
    InteractionWithUpdate<fluid_dynamics::ICEIntegration2ndHalfHLLERiemann> density_and_energy_relaxation(water_body_inner);
    SimpleDynamics<TaylorGreenInitialCondition> initial_condition(water_body);
    SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_body_correct_inner);
    InteractionWithUpdate<KernelCorrectionMatrixInner> kernel_correction_matrix(water_body_inner);
    InteractionWithUpdate<KernelCorrectionMatrixInner> kernel_correction_matrix_eu(water_body_correct_inner);
    InteractionDynamics<KernelGradientCorrectionInner> kernel_gradient_update(kernel_correction_matrix_eu);
    water_body.addBodyStateForRecording<Real>("Pressure");
    water_body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    /** Output the mechanical energy of fluid body. */
    RegressionTestEnsembleAverage<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_total_mechanical_energy(io_environment, water_body);
    /** Output the maximum speed of the fluid body. */
    RegressionTestEnsembleAverage<ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>>>
        write_maximum_speed(io_environment, water_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    kernel_correction_matrix.exec();
    kernel_correction_matrix_eu.exec();
    kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 5.0;
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    /** Output the mechanical energy of fluid. */
    write_total_mechanical_energy.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force. */
            time_step_initialization.exec();
            Real dt = get_fluid_time_step_size.exec();
            viscous_acceleration.exec();
            /** Dynamics including pressure relaxation. */
            integration_time += dt;
            pressure_relaxation.exec(dt);
            density_and_energy_relaxation.exec(dt);
            GlobalStaticVariables::physical_time_ += dt;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_total_mechanical_energy.writeToFile(number_of_iterations);
        write_maximum_speed.writeToFile(number_of_iterations);
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    write_total_mechanical_energy.testResult();
    write_maximum_speed.testResult();

    return 0;
}
