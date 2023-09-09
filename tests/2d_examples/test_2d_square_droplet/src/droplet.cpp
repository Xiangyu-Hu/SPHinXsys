/**
 * @file 	droplet.cpp
 * @brief 	A square droplet deforms to circle due to surface tension.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;                         /**< Tank length. */
Real DH = 2.0;                         /**< Tank height. */
Real LL = 1.0;                         /**< Liquid column length. */
Real LH = 1.0;                         /**< Liquid column height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of water. */
Real rho0_a = 0.001;     /**< Reference density of air. */
Real U_ref = 1.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
Real mu_f = 0.2;         /**< Water viscosity. */
Real mu_a = 0.002;       /**< Air viscosity. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.5, 0.5));
    water_block_shape.push_back(Vecd(0.5, LH + 0.5));
    water_block_shape.push_back(Vecd(LL + 0.5, LH + 0.5));
    water_block_shape.push_back(Vecd(LL + 0.5, 0.5));
    water_block_shape.push_back(Vecd(0.5, 0.5));
    return water_block_shape;
}
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, -BW));

    return outer_wall_shape;
}
std::vector<Vecd> createInnerWallShape()
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
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
// Air body shape definition.
//----------------------------------------------------------------------/**
class AirBlock : public MultiPolygonShape
{
  public:
    explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_shape = createOuterWallShape();
        std::vector<Vecd> inner_shape = createInnerWallShape();
        multi_polygon_.addAPolygon(outer_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_shape, ShapeBooleanOps::sub);
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
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<int>("Indicator");

    FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
    air_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
    air_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_block, {&air_block});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation air_water_complex(air_block, {&water_block});
    ContactRelation air_wall_contact(air_block, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block);
    SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
        update_water_density_by_summation(water_wall_contact, water_air_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
        update_air_density_by_summation(air_wall_contact, air_water_complex);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_wall_contact, air_water_complex, 2.0e-3);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
        water_transport_correction(water_air_complex, 2.0e-3);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_ref);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);
    /** Pressure relaxation for water by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
        water_pressure_relaxation(water_wall_contact, water_air_complex);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
        water_density_relaxation(water_wall_contact, water_air_complex);
    /** Extend Pressure relaxation is used for air. */
    Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
        air_pressure_relaxation(air_wall_contact, air_water_complex, 2.0);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
        air_density_relaxation(air_wall_contact, air_water_complex);
    /** Viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall>
        air_viscous_acceleration(air_wall_contact, air_water_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall>
        water_viscous_acceleration(water_wall_contact, water_air_complex);
    /** Surface tension. */
    InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationInner>
        surface_detection(water_air_complex.getInnerRelation());
    InteractionDynamics<fluid_dynamics::ColorFunctionGradientInner>
        color_gradient(water_air_complex.getInnerRelation());
    InteractionDynamics<fluid_dynamics::ColorFunctionGradientInterpolationInner>
        color_gradient_interpolation(water_air_complex.getInnerRelation());
    InteractionDynamics<fluid_dynamics::SurfaceTensionAccelerationInner>
        surface_tension_acceleration(water_air_complex.getInnerRelation(), 1.0);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 5.0;
    Real output_interval = 0.02; /**< Time stamps for output of body states. */
    Real dt = 0.0;               /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
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
            time_instance = TickCount::now();
            initialize_a_water_step.exec();
            initialize_a_air_step.exec();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_f, Dt_a);

            update_water_density_by_summation.exec();
            update_air_density_by_summation.exec();
            air_transport_correction.exec();
            water_transport_correction.exec();

            air_viscous_acceleration.exec();
            water_viscous_acceleration.exec();

            surface_detection.exec();
            color_gradient.exec();
            color_gradient_interpolation.exec();
            surface_tension_acceleration.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            water_block.updateCellLinkedListWithParticleSort(100);
            water_air_complex.updateConfiguration();
            water_wall_contact.updateConfiguration();

            air_block.updateCellLinkedListWithParticleSort(100);
            air_water_complex.updateConfiguration();
            air_wall_contact.updateConfiguration();

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
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";


    return 0;
}
