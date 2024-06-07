/**
 * @file 	dambreak_split_merge.cpp
 * @brief 	2D dambreak example with split and merge particles.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 * @version 0.1
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Tank length. */
Real DH = 5.366;                    /**< Tank height. */
Real LL = 2.0;                      /**< Liquid column length. */
Real LH = 1.0;                      /**< Liquid column height. */
Real particle_spacing_ref = 0.04;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity force of fluid. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d water_block_translation = water_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
/**
 * application of refinement area
 */
MultiPolygon createRefinementArea()
{
    /** Geometry definition. */
    std::vector<Vecd> refinement_area;
    refinement_area.push_back(Vecd(4.0, -0.1));
    refinement_area.push_back(Vecd(4.0, 6.0));
    refinement_area.push_back(Vecd(6.0, 6.0));
    refinement_area.push_back(Vecd(6.0, -0.1));
    refinement_area.push_back(Vecd(4.0, -0.1));
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(refinement_area, ShapeBooleanOps::add);
    return multi_polygon;
};

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    /** I/O environment. */
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    TransformShape<GeometricShapeBox> initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineAdaptation<ParticleSplitAndMerge>(1.3, 1.0, 1);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    ParticleBuffer<ReserveSizeFactor> particle_split_buffer(20.0);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(particle_split_buffer);
    water_block.addBodyStateForRecording<Real>("SmoothingLengthRatio");
    water_block.addBodyStateForRecording<Real>("VolumetricMeasure");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<BaseParticles, Observer>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    AdaptiveInnerRelation water_inner(water_block);
    AdaptiveContactRelation water_contact(water_block, {&wall_boundary});
    AdaptiveContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_complex(water_inner, water_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_inner, water_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_inner, water_contact);

    MultiPolygonShape split_merge_region(createRefinementArea());
    InteractionWithUpdate<SplitWithMinimumDensityErrorWithWall, SequencedPolicy> particle_split_(water_inner, water_contact, split_merge_region);
    InteractionDynamics<MergeWithMinimumDensityErrorWithWall, SequencedPolicy> particle_merge_(water_inner, water_contact, split_merge_region);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplexAdaptive> fluid_density_by_summation(water_inner, water_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real End_Time = 20.0; /**< End time. */
    Real D_Time = 0.1;    /**< Time stamps for output of body states. */
    Real dt = 0.0;        /**< Default acoustic time step sizes. */
    int refinement_interval = 1;
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
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    write_recorded_water_pressure.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();
            fluid_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                fluid_pressure_relaxation.exec(dt);
                fluid_density_relaxation.exec(dt);
                dt = fluid_acoustic_time_step.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != 0)
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    write_recorded_water_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            /** Particle Refinement */
            if (number_of_iterations % refinement_interval == 0)
            {
                particle_split_.exec();
                particle_merge_.exec();
            }
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
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

    return 0;
};
