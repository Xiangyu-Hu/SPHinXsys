/**
 * @file 	filling_tank.cpp
 * @brief 	2D example to show that a tank is filled by emitter.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding how emitter inflow is working.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;              /**< Tank length. */
Real DH = 5.366;              /**< Tank height. */
Real resolution_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4; /**< Extending width for wall boundary. */
Real LL = 2.0 * BW;           /**< Inflow region length. */
Real LH = 0.125;              /**< Inflows region height. */
Real inlet_height = 1.0;      /**< Inflow location height */
Real inlet_distance = -BW;    /**< Inflow location distance */
Vec2d inlet_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d inlet_translation = Vec2d(inlet_distance, inlet_height) + inlet_halfsize;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
Real rho0_f = 1.0;                                      /**< Reference density of fluid. */
Real gravity_g = 1.0;                                   /**< Gravity force of fluid. */
Real U_f = 2.0 * sqrt(gravity_g * (inlet_height + LH)); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                  /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometries
//----------------------------------------------------------------------
/** create a outer wall polygon. */
std::vector<Vecd> CreateOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, -BW));

    return outer_wall_shape;
}
/** create a inner wall polygon. */
std::vector<Vecd> CreateInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, DH));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, 0.0));

    return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Case-dependent wall boundary
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(CreateOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(CreateInnerWallShape(), ShapeBooleanOps::sub);
        multi_polygon_.addABox(Transform(inlet_translation), inlet_halfsize, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
  public:
    InletInflowCondition(BodyAlignedBoxByParticle &aligned_box_part)
        : EmitterInflowCondition(aligned_box_part) {}

  protected:
    virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
    {
        return Vec2d(2.0, 0.0);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    /** Build up a SPHSystem */
    SPHSystem system(system_domain_bounds, resolution_ref);
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(system, makeShared<TransformShape<GeometricShapeBox>>(
                                     Transform(inlet_translation), inlet_halfsize, "WaterBody"));
    water_body.sph_adaptation_->resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_body.generateParticles<ParticleGeneratorLattice>();
    /**note that, as particle sort is activated (by default) for fluid particles,
     * the output occasionally does not reflect the real free surface indication due to sorting. */
    SolidBody wall(system, makeShared<WallBoundary>("Wall"));
    wall.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    ComplexRelation water_body_complex(water_body, {&wall});
    ContactRelation fluid_observer_contact_relation(fluid_observer, {&water_body});
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_body_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_body_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_body_complex);
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
        indicate_free_surface(water_body_complex);
    water_body.addBodyStateForRecording<Real>("PositionDivergence"); // for debug
    water_body.addBodyStateForRecording<int>("Indicator");    // for debug

    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_body, gravity_ptr);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    BodyAlignedBoxByParticle emitter(
        water_body, makeShared<AlignedBoxShape>(Transform(inlet_translation), inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, 350, 0);

    //----------------------------------------------------------------------
    //	File Output
    //----------------------------------------------------------------------
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp body_states_recording(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_water_mechanical_energy(io_environment, water_body, gravity_ptr);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact_relation);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    indicate_free_surface.exec();
    //----------------------------------------------------------------------
    //	Time stepping control parameters.
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 30.0;
    Real output_interval = 0.1;
    Real dt = 0.0; /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force and gravity. */
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);
                inflow_condition.exec();
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** inflow emitter injection*/
            emitter_injection.exec();
            /** Update cell linked list and configuration. */

            water_body.updateCellLinkedListWithParticleSort(100);
            water_body_complex.updateConfiguration();
            fluid_observer_contact_relation.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_water_mechanical_energy.writeToFile(number_of_iterations);
        indicate_free_surface.exec();
        body_states_recording.writeToFile();
        write_recorded_water_pressure.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    write_water_mechanical_energy.testResult();
    write_recorded_water_pressure.testResult();

    return 0;
}
