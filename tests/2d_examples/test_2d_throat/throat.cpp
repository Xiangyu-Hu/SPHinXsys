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
Real global_resolution = 0.1;      // particle spacing
Real BW = global_resolution * 4.0; // boundary width
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
    BoundingBoxd system_domain_bounds(Vec2d(-0.5 * DL - BW, -0.5 * DH - BW),
                                     Vec2d(0.5 * DL + BW, 0.5 * DH + BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("FluidBody"));
    fluid_block.defineClosure<WeaklyCompressibleFluid, OldroydBViscosity>(
        ConstructArgs(rho0_f, c_f), ConstructArgs(mu_f, lambda_f, mu_p_f));
    Ghost<PeriodicAlongAxis> ghost_along_x(fluid_block.getSPHBodyBounds(), xAxis);
    fluid_block.generateParticlesWithReserve<BaseParticles, Lattice>(ghost_along_x);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd::Zero()};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation fluid_block_inner(fluid_block);
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation fluid_block_contact(fluid_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&fluid_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation fluid_block_complex(fluid_block_inner, fluid_block_contact);
    //-------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-------------------------------------------------------------------
    Gravity gravity(Vecd(gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(fluid_block, gravity);
    InteractionDynamics<NormalDirectionFromParticles> wall_boundary_normal_direction(wall_boundary_inner);
    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(fluid_block_contact);

    Dynamics1Level<fluid_dynamics::Oldroyd_BIntegration1stHalfWithWall> pressure_relaxation(fluid_block_inner, fluid_block_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<NoKernelCorrection>> update_velocity_gradient(fluid_block_inner, fluid_block_contact);
    Dynamics1Level<fluid_dynamics::Oldroyd_BIntegration2ndHalfWithWall> density_relaxation(fluid_block_inner, fluid_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_block_inner, fluid_block_contact);
    InteractionSplit<DampingPairwiseWithWall<Vec2d, FixedDampingRate>> implicit_viscous_damping(
        DynamicsArgs(fluid_block_inner, "Velocity", mu_f), DynamicsArgs(fluid_block_contact, "Velocity", mu_f));
    InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_block_inner, fluid_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(fluid_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(fluid_block);
    PeriodicConditionUsingGhostParticles periodic_condition(fluid_block, ghost_along_x);
    pressure_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);
    density_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_block_inner);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(fluid_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>>
        write_fluid_kinetic_energy(fluid_block);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_fluid_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    // initial periodic boundary condition
    periodic_condition.ghost_creation_.exec();
    sph_system.initializeSystemConfigurations();
    // prepare quantities will be used once only
    wall_boundary_normal_direction.exec();
    distance_to_wall.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 10;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 40.0;
    // time step size for ouput file
    Real output_interval = end_time / 40.0;
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
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < output_interval)
        {

            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            transport_velocity_correction.exec();
            distance_to_wall.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                implicit_viscous_damping.exec(dt);
                pressure_relaxation.exec(dt);
                update_velocity_gradient.exec();
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_fluid_kinetic_energy.writeToFile(number_of_iterations);
                    write_recorded_fluid_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            // water block configuration and periodic condition
            periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            fluid_block.updateCellLinkedList();
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

    if (sph_system.GenerateRegressionData())
    {
        write_fluid_kinetic_energy.generateDataBase(1.0e-2);
        write_recorded_fluid_pressure.generateDataBase(1.0e-2);
    }
    else
    {
        write_fluid_kinetic_energy.testResult();
        write_recorded_fluid_pressure.testResult();
    }

    return 0;
}
