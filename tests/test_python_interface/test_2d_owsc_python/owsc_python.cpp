/**
 * @file     owsc_python.cpp
 * @brief    This is the test of wave interaction with Oscillating Wave Surge Converter (OWSC), for DRL.
 * @author   Mai Ye, Chi Zhang and Xiangyu Hu
 */
#include "owsc_python.h"
#include "custom_io_observation.h"
#include "custom_io_simbody.h"
#include "sphinxsys.h"
#include <pybind11/pybind11.h>

using namespace SPH;
namespace py = pybind11;

class SphBasicSystemSetting : public SphBasicGeometrySetting
{
  protected:
    BoundingBoxd system_domain_bounds;
    SPHSystem sph_system;

  public:
    SphBasicSystemSetting(int parallel_env, int episode_env)
        : system_domain_bounds(Vec2d(-DL_Extra - BW, -BW), Vec2d(DL + BW, DH + BW)),
          sph_system(system_domain_bounds, particle_spacing_ref)
    {
        sph_system.getIOEnvironment().appendOutputFolder(
            "env_" + std::to_string(parallel_env) + "_episode_" + std::to_string(episode_env));
    }
};

class SphFlapReloadEnvironment : public SphBasicSystemSetting
{
  protected:
    FluidBody water_block;
    SolidBody wall_boundary, flap;
    ObserverBody flap_observer, wave_velocity_observer;

  public:
    SphFlapReloadEnvironment(int parallel_env, int episode_env)
        : SphBasicSystemSetting(parallel_env, episode_env),
          water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
          wall_boundary(sph_system, makeShared<WallBoundary>("Wall")),
          flap(sph_system, makeShared<Flap>("Flap")),
          flap_observer(sph_system, "FlapObserver"),
          wave_velocity_observer(sph_system, "WaveVelocityObserver")
    {
        water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
        water_block.generateParticles<BaseParticles, Lattice>();

        wall_boundary.defineMaterial<Solid>();
        wall_boundary.generateParticles<BaseParticles, Lattice>();

        flap.defineMaterial<Solid>(rho0_s);
        flap.generateParticles<BaseParticles, Lattice>();

        flap_observer.generateParticles<ObserverParticles>(createFlapObserver());
        wave_velocity_observer.generateParticles<ObserverParticles>(createWaveVelocityObserver());
    }
};

class SimbodyEnvironment : public SphFlapReloadEnvironment
{
  protected:
    /** set up the multi body system. */
    SimTK::MultibodySystem MBsystem;
    /** the bodies or matter of the system. */
    SimTK::SimbodyMatterSubsystem matter;
    /** the forces of the system. */
    SimTK::GeneralForceSubsystem forces;
    /** mass properties of the fixed spot. */
    FlapSystemForSimbody flap_multibody;
    SimTK::Body::Rigid pin_spot_info;
    SimTK::MobilizedBody::Pin pin_spot;
    SimTK::Force::UniformGravity sim_gravity;
    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies;
    SimTK::Force::MobilityLinearDamper linear_damper;
    SimTK::State state;
    SimTK::RungeKuttaMersonIntegrator integ;

  public:
    SimbodyEnvironment(int parallel_env, int episode_env)
        : SphFlapReloadEnvironment(parallel_env, episode_env),
          matter(MBsystem),
          forces(MBsystem),
          flap_multibody(flap, makeShared<MultiPolygonShape>(createFlapSimbodyConstrainShape(), "FlapMultiBody")),
          pin_spot_info(*flap_multibody.body_part_mass_properties_),
          pin_spot(matter.Ground(), SimTK::Transform(SimTK::Vec3(7.92, 0.315, 0.0)), pin_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0))),
          sim_gravity(forces, matter, SimTK::Vec3(0.0, -gravity_g, 0.0), 0.0),
          force_on_bodies(forces, matter),
          linear_damper(forces, pin_spot, SimTK::MobilizerUIndex(0), 20), integ(MBsystem)
    {
        pin_spot.setDefaultAngle(0);
        state = MBsystem.realizeTopology();
        integ.setAccuracy(1e-3);
        integ.setAllowInterpolation(false);
        integ.initialize(state);
    }
};

class SphOWSC : public SimbodyEnvironment
{
  protected:
    SPHSystem &sph_system_;
    InnerRelation water_block_inner, flap_inner;
    ContactRelation water_block_contact, flap_contact, flap_observer_contact_with_water, flap_observer_contact_with_flap,
        wave_velocity_observer_contact_with_water;
    ComplexRelation water_block_complex;
    //----------------------------------------------------------------------
    //	    Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    Gravity gravity;
    SimpleDynamics<GravityForce<Gravity>> constant_gravity;
    SimpleDynamics<OffsetInitialPosition> flap_offset_position;
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction, flap_normal_direction;

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation;
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation;
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force;

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size;
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size;
    BodyRegionByCell damping_buffer;
    SimpleDynamics<fluid_dynamics::DampingBoundaryCondition> damping_wave;
    BodyRegionByParticle wave_maker;
    SimpleDynamics<WaveMaking> wave_making;

    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid;
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid;

    ParticleSorting<ParallelPolicy> particle_sorting;
    //----------------------------------------------------------------------
    //	    Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_spot_flap;
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_spot_flap;
    //----------------------------------------------------------------------
    //	    Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states;
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid;
    /** Velocity probe. */
    ObservedQuantityRecording<Vecd> wave_velocity_probe;
    ObservedQuantityRecording<Vecd> wave_velocity_on_flap_probe;
    /** Observer for flap position. */
    ObservedQuantityRecording<Vecd> flap_position_probe;
    // Interpolate the particle position in flap to move the observer accordingly.
    // Seems not used? TODO: observe displacement more accurate.
    InteractionDynamics<InterpolatingAQuantity<Vecd>> interpolation_flap_velocity_observer_position;
    WriteSimBodyPinDataExtended write_flap_pin_data;
    /** WaveProbes. */
    BodyRegionByCell wave_probe_buffer_no_0, wave_probe_buffer_no_1;
    ExtendedReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_0, wave_probe_1;
    //----------------------------------------------------------------------
    //	    Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    Real dt = 0.0;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    Real output_interval = 0.1;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;

  public:
    explicit SphOWSC(int parallel_env, int episode_env)
        : SimbodyEnvironment(parallel_env, episode_env),
          sph_system_(sph_system),
          water_block_inner(water_block),
          flap_inner(flap),
          water_block_contact(water_block, {&wall_boundary, &flap}),
          flap_contact(flap, {&water_block}),
          flap_observer_contact_with_water(flap_observer, {&water_block}),
          flap_observer_contact_with_flap(flap_observer, {&flap}),
          wave_velocity_observer_contact_with_water(wave_velocity_observer, {&water_block}),
          water_block_complex(water_block_inner, water_block_contact),

          gravity(Vecd(0.0, -gravity_g)),
          constant_gravity(water_block, gravity),
          flap_offset_position(flap, offset),
          wall_boundary_normal_direction(wall_boundary),
          flap_normal_direction(flap),

          pressure_relaxation(water_block_inner, water_block_contact),
          density_relaxation(water_block_inner, water_block_contact),
          update_density_by_summation(water_block_inner, water_block_contact),
          viscous_force(water_block_inner, water_block_contact),

          get_fluid_advection_time_step_size(water_block, U_f),
          get_fluid_time_step_size(water_block),
          damping_buffer(water_block, makeShared<MultiPolygonShape>(createDampingBufferShape())),
          damping_wave(damping_buffer),
          wave_maker(wall_boundary, makeShared<MultiPolygonShape>(createWaveMakerShape())),
          wave_making(wave_maker),

          viscous_force_from_fluid(flap_contact),
          pressure_force_from_fluid(flap_contact),

          particle_sorting(water_block),

          force_on_spot_flap(flap_multibody, MBsystem, pin_spot, integ),
          constraint_spot_flap(flap_multibody, MBsystem, pin_spot, integ),
          write_real_body_states(sph_system),
          write_total_viscous_force_from_fluid(flap, "ViscousForceFromFluid"),
          wave_velocity_probe("Velocity", wave_velocity_observer_contact_with_water),
          wave_velocity_on_flap_probe("Velocity", flap_observer_contact_with_water),
          flap_position_probe("Position", flap_observer_contact_with_flap),
          interpolation_flap_velocity_observer_position(flap_observer_contact_with_flap, "Position", "Position"),
          write_flap_pin_data(sph_system, integ, pin_spot),
          wave_probe_buffer_no_0(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(3.0), "WaveProbe_03")),
          wave_probe_buffer_no_1(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(5.0), "WaveProbe_05")),
          wave_probe_0(wave_probe_buffer_no_0, "FreeSurfaceHeight"),
          wave_probe_1(wave_probe_buffer_no_1, "FreeSurfaceHeight")
    {
        physical_time = 0.0;
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        flap_offset_position.exec();
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        wall_boundary_normal_direction.exec();
        flap_normal_direction.exec();
        constant_gravity.exec();
        //----------------------------------------------------------------------
        //	First output before the main loop.
        //----------------------------------------------------------------------
        write_real_body_states.writeToFile(0);
        write_total_viscous_force_from_fluid.writeToFile(0);
        wave_velocity_probe.writeToFile(0);
        wave_velocity_on_flap_probe.writeToFile(0);
        flap_position_probe.writeToFile(0);
        write_flap_pin_data.writeToFile(0);
        wave_probe_0.writeToFile(0);
        wave_probe_1.writeToFile(0);
    }

    virtual ~SphOWSC() {};
    //----------------------------------------------------------------------
    //	    For ctest.
    //----------------------------------------------------------------------
    int cmakeTest()
    {
        return 1;
    }
    //----------------------------------------------------------------------
    //	    Get flap angle and angle rate.
    //----------------------------------------------------------------------
    Real getFlapAngle()
    {
        return write_flap_pin_data.getAngleToPython(number_of_iterations);
    };

    Real getFlapAngleRate()
    {
        return write_flap_pin_data.getAngleRateToPython(number_of_iterations);
    };
    //----------------------------------------------------------------------
    //  Get wave height.
    //----------------------------------------------------------------------
    Real getWaveHeight(size_t number)
    {
        if (number == 0)
        {
            return wave_probe_0.getReducedQuantity();
        }
        else if (number == 1)
        {
            return wave_probe_1.getReducedQuantity();
        }
        return -1.0;
    }
    //----------------------------------------------------------------------
    //  Get wave velocity in front of the flap.
    //----------------------------------------------------------------------
    Real getWaveVelocity(int number, int direction)
    {
        return wave_velocity_probe.getObservedQuantity()[number][direction];
    };
    //----------------------------------------------------------------------
    //	     Get wave velocity on the flap.
    //----------------------------------------------------------------------
    Real getWaveVelocityOnFlap(int number, int direction)
    {
        return wave_velocity_on_flap_probe.getObservedQuantity()[number][direction];
    };
    //----------------------------------------------------------------------
    //	    Get flap position.
    //----------------------------------------------------------------------
    Real getFlapPositon(int number, int direction)
    {
        return flap_position_probe.getObservedQuantity()[number][direction];
    };
    //----------------------------------------------------------------------
    //	    Main loop of time stepping starts here. && For changing damping coefficient.
    //----------------------------------------------------------------------
    void runCase(Real pause_time_from_python, Real dampling_coefficient_from_python)
    {
        while (physical_time < pause_time_from_python)
        {
            Real integral_time = 0.0;
            while (integral_time < output_interval)
            {
                Real Dt = get_fluid_advection_time_step_size.exec();
                update_density_by_summation.exec();
                viscous_force.exec();
                /** Viscous force exerting on flap. */
                viscous_force_from_fluid.exec();

                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    pressure_relaxation.exec(dt);
                    pressure_force_from_fluid.exec();
                    density_relaxation.exec(dt);
                    /** coupled rigid body dynamics. */
                    if (total_time >= relax_time)
                    {
                        SimTK::State &state_for_update = integ.updAdvancedState();
                        linear_damper.setDamping(state_for_update, dampling_coefficient_from_python);
                        force_on_bodies.clearAllBodyForces(state_for_update);
                        force_on_bodies.setOneBodyForce(state_for_update, pin_spot, force_on_spot_flap.exec());
                        integ.stepBy(dt);
                        constraint_spot_flap.exec();
                        wave_making.exec(dt);
                    }
                    interpolation_flap_velocity_observer_position.exec();

                    dt = get_fluid_time_step_size.exec();
                    relaxation_time += dt;
                    integral_time += dt;
                    total_time += dt;
                    if (total_time >= relax_time)
                        physical_time += dt;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                              << "	Total Time = " << total_time
                              << "	Physical Time = " << physical_time
                              << "	Dt = " << Dt << "	dt = " << dt << "\n";
                }
                number_of_iterations++;
                damping_wave.exec(Dt);
                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                water_block.updateCellLinkedList();
                wall_boundary.updateCellLinkedList();
                flap.updateCellLinkedList();
                water_block_complex.updateConfiguration();
                flap_contact.updateConfiguration();
                flap_observer_contact_with_water.updateConfiguration();
                wave_velocity_observer_contact_with_water.updateConfiguration();
                if (total_time >= relax_time)
                {
                    write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
                    wave_velocity_probe.writeToFile(number_of_iterations);
                    wave_velocity_on_flap_probe.writeToFile(number_of_iterations);
                    flap_position_probe.writeToFile(number_of_iterations);
                    write_flap_pin_data.writeToFile(physical_time);
                    wave_probe_0.writeToFile(number_of_iterations);
                    wave_probe_1.writeToFile(number_of_iterations);
                }
            }

            TickCount t2 = TickCount::now();
            if (total_time >= relax_time)
                write_real_body_states.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;

        // This section is used for CMake testing (cmake test).
        // During reinforcement learning training, this part can be commented out.
        if (sph_system.GenerateRegressionData())
        {
            write_total_viscous_force_from_fluid.generateDataBase(1.0e-3);
        }
        else
        {
            write_total_viscous_force_from_fluid.testResult();
        }
    };
};

PYBIND11_MODULE(test_2d_owsc_python, m)
{
    py::class_<SphOWSC>(m, "owsc_from_sph_cpp")
        .def(py::init<const int &, const int &>())
        .def("cmake_test", &SphOWSC::cmakeTest)
        .def("get_flap_angle", &SphOWSC::getFlapAngle)
        .def("get_flap_angle_rate", &SphOWSC::getFlapAngleRate)
        .def("get_wave_height", &SphOWSC::getWaveHeight)
        .def("get_wave_velocity", &SphOWSC::getWaveVelocity)
        .def("get_wave_velocity_on_flap", &SphOWSC::getWaveVelocityOnFlap)
        .def("get_flap_position", &SphOWSC::getFlapPositon)
        .def("run_case", &SphOWSC::runCase);
}
