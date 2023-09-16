/**
 * @file 	throat.cpp
 * @brief 	2D in a channel with a throat.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for non-Newtonian low Reynolds number flows.
 *			The choice of maximum flow speed, sound speed and time step size follows
 *			Morris et al. Modeling Low Reynolds Number Incompressible Flows Using SPH.
 *			Journal of Computational Physics, Volume 136, 1997, 214-226.
 *			https://doi.org/10.1006/jcph.1997.5776
 *			Note that as we use implicit time stepping for the viscous term,
 *			the time step size does not need to follow the viscous time step criteria
 *			and is the same of that for pressure and density relaxations.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 4.0;                  // channel height
Real DT = 1.0;                  // throat height
Real DL = 24.0;                 // channel length
Real resolution_ref = 0.1;      // particle spacing
Real BW = resolution_ref * 4.0; // boundary width
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real gravity_g = 1.0; /**< Gravity force of fluid. */
Real Re = 0.001;      /**< Reynolds number defined in the channel */
// obtain viscosity according planar Poiseuille flow solution in the channel
Real mu_f = rho0_f * sqrt(0.5 * rho0_f * pow(0.5 * DH, 3) * gravity_g / Re);
// maximum flow velocity in the channel
Real U_c = 0.5 * pow(0.5 * DH, 2) * gravity_g * rho0_f / mu_f;
//	predicted overall maximum velocity for this case is in the throat according to incompressible condition
Real U_f = U_c * DH / DT;
// For low Reynolds number flow the weakly compressible formulation need to
// consider viscosity for artificial sound speed.
Real c_f = 10.0 * SMAX(U_f, sqrt(mu_f / rho0_f * U_f / DT));
Real mu_p_f = 0.6 * mu_f;
Real lambda_f = 10.0;
//----------------------------------------------------------------------
//	Fluid body cases-dependent geometries.
//----------------------------------------------------------------------
class FluidBlock : public MultiPolygonShape
{
  public:
    explicit FluidBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> pnts;
        pnts.push_back(Vecd(-0.5 * DL, -0.5 * DH));
        pnts.push_back(Vecd(-0.5 * DL, 0.5 * DH));
        pnts.push_back(Vecd(-DL / 6.0, 0.5 * DH));
        pnts.push_back(Vecd(-DL / 6.0, -0.5 * DH));
        pnts.push_back(Vecd(-0.5 * DL, -0.5 * DH));

        std::vector<Vecd> pnts1;
        pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5 * DT));
        pnts1.push_back(Vecd(-DL / 6.0 - BW, 0.5 * DT));
        pnts1.push_back(Vecd(DL / 6.0 + BW, 0.5 * DT));
        pnts1.push_back(Vecd(DL / 6.0 + BW, -0.5 * DT));
        pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5 * DT));

        std::vector<Vecd> pnts2;
        pnts2.push_back(Vecd(DL / 6.0, -0.5 * DH));
        pnts2.push_back(Vecd(DL / 6.0, 0.5 * DH));
        pnts2.push_back(Vecd(0.5 * DL, 0.5 * DH));
        pnts2.push_back(Vecd(0.5 * DL, -0.5 * DH));
        pnts2.push_back(Vecd(DL / 6.0, -0.5 * DH));

        multi_polygon_.addAPolygon(pnts, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(pnts1, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(pnts2, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Cases-dependent wall boundary geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> pnts3;
        pnts3.push_back(Vecd(-0.5 * DL - BW, -0.5 * DH - BW));
        pnts3.push_back(Vecd(-0.5 * DL - BW, 0.5 * DH + BW));
        pnts3.push_back(Vecd(0.5 * DL + BW, 0.5 * DH + BW));
        pnts3.push_back(Vecd(0.5 * DL + BW, -0.5 * DH - BW));
        pnts3.push_back(Vecd(-0.5 * DL - BW, -0.5 * DH - BW));

        std::vector<Vecd> pnts;
        pnts.push_back(Vecd(-0.5 * DL - 2.0 * BW, -0.5 * DH));
        pnts.push_back(Vecd(-0.5 * DL - 2.0 * BW, 0.5 * DH));
        pnts.push_back(Vecd(-DL / 6.0, 0.5 * DH));
        pnts.push_back(Vecd(-DL / 6.0, -0.5 * DH));
        pnts.push_back(Vecd(-0.5 * DL - 2.0 * BW, -0.5 * DH));

        std::vector<Vecd> pnts1;
        pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5 * DT));
        pnts1.push_back(Vecd(-DL / 6.0 - BW, 0.5 * DT));
        pnts1.push_back(Vecd(DL / 6.0 + BW, 0.5 * DT));
        pnts1.push_back(Vecd(DL / 6.0 + BW, -0.5 * DT));
        pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5 * DT));

        std::vector<Vecd> pnts2;
        pnts2.push_back(Vecd(DL / 6.0, -0.5 * DH));
        pnts2.push_back(Vecd(DL / 6.0, 0.5 * DH));
        pnts2.push_back(Vecd(0.5 * DL + 2.0 * BW, 0.5 * DH));
        pnts2.push_back(Vecd(0.5 * DL + 2.0 * BW, -0.5 * DH));
        pnts2.push_back(Vecd(DL / 6.0, -0.5 * DH));

        multi_polygon_.addAPolygon(pnts3, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(pnts, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(pnts1, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(pnts2, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-0.5 * DL - BW, -0.5 * DH - BW),
                                     Vec2d(0.5 * DL + BW, 0.5 * DH + BW));
    SPHSystem system(system_domain_bounds, resolution_ref);
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(system, makeShared<FluidBlock>("FluidBody"));
    fluid_block.defineParticlesAndMaterial<BaseParticles, Oldroyd_B_Fluid>(rho0_f, c_f, mu_f, lambda_f, mu_p_f);
    fluid_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd::Zero()};
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation fluid_block_inner(fluid_block);
    ComplexRelation fluid_block_complex(fluid_block_inner, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&fluid_block});
    //-------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-------------------------------------------------------------------
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingGhostParticles periodic_condition(fluid_block, fluid_block.getBodyShapeBounds(), xAxis);
    // evaluation of density by summation approach
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_block_complex);
    // time step size without considering sound wave speed and viscosity
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSizeForImplicitViscosity> get_fluid_advection_time_step_size(fluid_block, U_f);
    // time step size with considering sound wave speed
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(fluid_block);
    // pressure relaxation using verlet time stepping
    Dynamics1Level<fluid_dynamics::Oldroyd_BIntegration1stHalfWithWall> pressure_relaxation(fluid_block_complex);
    pressure_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);
    Dynamics1Level<fluid_dynamics::Oldroyd_BIntegration2ndHalfWithWall> density_relaxation(fluid_block_complex);
    density_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);
    // define external force
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(fluid_block, makeShared<Gravity>(Vecd(gravity_g, 0.0)));
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_block_complex);
    // computing viscous effect implicitly and with update velocity directly other than viscous acceleration
    InteractionSplit<DampingPairwiseWithWall<Vec2d, DampingPairwiseInner>>
        implicit_viscous_damping(fluid_block_complex, "Velocity", mu_f);
    // impose transport velocity
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_block_complex);
    // computing vorticity in the flow
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_block_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_fluid_mechanical_energy(io_environment, fluid_block);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_fluid_pressure("Pressure", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    // initial periodic boundary condition
    periodic_condition.ghost_creation_.exec();
    system.initializeSystemConfigurations();
    // prepare quantities will be used once only
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 10;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 20.0;
    // time step size for ouput file
    Real output_interval = end_time / 20.0;
    Real dt = 0.0; // default acoustic time step sizes
    // statistics for computing time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < output_interval)
        {

            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                implicit_viscous_damping.exec(dt);
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != system.RestartStep())
                {
                    write_fluid_mechanical_energy.writeToFile(number_of_iterations);
                    write_recorded_fluid_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            // water block configuration and periodic condition
            periodic_condition.bounding_.exec();
            fluid_block.updateCellLinkedListWithParticleSort(100);
            periodic_condition.ghost_creation_.exec();
            fluid_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_fluid_mechanical_energy.generateDataBase(1.0e-2);
        write_recorded_fluid_pressure.generateDataBase(1.0e-2);
    }
    else
    {
        write_fluid_mechanical_energy.testResult();
        write_recorded_fluid_pressure.testResult();
    }

    return 0;
}
