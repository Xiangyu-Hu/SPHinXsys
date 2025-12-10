/**
 * @file 	taylor_green.cpp
 * @brief 	2D taylor_green vortex flow example.
 * @details This is the one of the basic test cases.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // SPHinXsys namespace.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;                     /**< box length. */
Real DH = 1.0;                     /**< box height. */
Real resolution_ref = 1.0 / 100.0; /**< Global reference resolution. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                  /**< Reference density of fluid. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Reference sound speed. */
Real Re = 100;                      /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Fluid body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class TaylorGreenInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit TaylorGreenInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        vel_[index_i][0] = -cos(2.0 * Pi * pos_[index_i][0]) *
                           sin(2.0 * Pi * pos_[index_i][1]);
        vel_[index_i][1] = sin(2.0 * Pi * pos_[index_i][0]) *
                           cos(2.0 * Pi * pos_[index_i][1]);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d::Zero(), Vec2d(DL, DH));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    // Using relaxed particle distribution if needed
    sph_system.ReloadParticles()
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<TaylorGreenInitialCondition> initial_condition(water_block);
    /** Pressure relaxation algorithm by using Verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous.
     * The other reason is that we are using transport velocity formulation,
     * which will also introduce numerical dissipation slightly. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerNoRiemann> density_relaxation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationInner> update_density_by_summation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<TruncatedLinear, AllParticles>> transport_velocity_correction(water_block_inner);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicAlongAxis periodic_along_y(water_block.getSPHBodyBounds(), yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, periodic_along_x);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_block, periodic_along_y);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    ReloadParticleIO write_particle_reload_files(water_block);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>>
        write_total_kinetic_energy(water_block);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MaximumSpeed>>
        write_maximum_speed(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    initial_condition.exec();
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 5.0;
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;              /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_total_kinetic_energy.writeToFile(0);
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                // avoid possible smaller acoustic time step size for viscous flow
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                relaxation_time += dt;
                integration_time += dt;
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            water_block.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            water_block_inner.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_total_kinetic_energy.writeToFile(number_of_iterations);
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

    write_particle_reload_files.writeToFile();

    if (sph_system.GenerateRegressionData())
    {
        write_total_kinetic_energy.generateDataBase(1.0e-3);
        write_maximum_speed.generateDataBase(1.0e-3);
    }
    else if (!sph_system.ReloadParticles())
    {
        write_total_kinetic_energy.testResult();
        write_maximum_speed.testResult();
    }

    return 0;
}
