/**
 * @file    2d_flow_stream_around_fish.cpp
 * @brief   fish swimming driven by active muscles
 * @author  Yaru Ren and Xiangyu Hu
 */
#include <pybind11/pybind11.h>
#include "2d_flow_stream_around_fish.h"
#include "sphinxsys.h"
namespace py = pybind11;
using namespace SPH;


class SphFishEnvironment : public SphBasicGeometrySetting
{
  protected:
    BoundingBox system_domain_bounds;
    SPHSystem system;
    IOEnvironment io_environment;
    FluidBody water_block;
    SolidBody fish_body;
    ObserverBody fish_observer;

  public:
    SphFishEnvironment(int episode_env):
        system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW)),
        system(system_domain_bounds, particle_spacing_ref),
        io_environment(system, episode_env),
        water_block(system, makeShared<WaterBlock>("WaterBody")),
        fish_body(system, makeShared<FishBody>("FishBody")),
        fish_observer(system, "FishObserver")
    {
        /** Tag for computation start with relaxed body fitted particles distribution. */
        system.setRunParticleRelaxation(false);
        system.setReloadParticles(true);

        water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        water_block.generateParticles<ParticleGeneratorLattice>();

        fish_body.defineAdaptationRatios(1.15, 2.0);
        fish_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        fish_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
        (!system.RunParticleRelaxation() && system.ReloadParticles())
            ? fish_body.generateParticles<ParticleGeneratorReload>(io_environment, fish_body.getName())
            : fish_body.generateParticles<ParticleGeneratorLattice>();

        fish_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    }
};

class SphFishSimulation : public SphFishEnvironment
{
  protected:
    InnerRelation fish_inner, water_block_inner;
    ComplexRelation water_block_complex;
    ContactRelation fish_contact, fish_observer_contact_water, fish_observer_contact_solid;
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step;
    BodyAlignedBoxByParticle emitter;
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection;
    BodyAlignedBoxByCell emitter_buffer;
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition;
    BodyAlignedBoxByCell disposer;
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion;
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex> free_stream_surface_indicator;
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density;
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size;
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size;
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint;
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation;
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration;
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction;
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity;
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> fish_body_normal_direction;
    InteractionWithUpdate<CorrectedConfigurationInner> fish_body_corrected_configuration;
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid;
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluidRiemann> fluid_force_on_fish_update;
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration;
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<FishMaterialInitialization> composite_material_id;
    SimpleDynamics<ImposingActiveStrain> imposing_active_strain;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> fish_body_computing_time_step_size;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> fish_body_stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> fish_body_stress_relaxation_second_half;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> fish_body_update_normal;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states;
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>> write_total_viscous_force_on_insert_body;
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>> write_total_force_on_insert_body;
    ObservedQuantityRecording<Real> pressure_probe;
    ObservedQuantityRecording<Vecd> fish_velocity_probe;
    ObservedQuantityRecording<Vecd> fish_position_probe;
    InteractionDynamics<InterpolatingAQuantity<Vecd>> interpolation_particle_position;
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real D_Time = 0.01;  /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

  public:
    explicit SphFishSimulation(int episode_simulation) :  
        SphFishEnvironment(episode_simulation),
        fish_inner(fish_body), water_block_inner(water_block),
        water_block_complex(water_block, {&fish_body}),
        fish_contact(fish_body, {&water_block}), fish_observer_contact_water(fish_observer, { &water_block }), fish_observer_contact_solid(fish_observer, { &fish_body }),
        initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero())),
        emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize)),
        emitter_inflow_injection(emitter, 10, 0),
        emitter_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize)),
        emitter_buffer_inflow_condition(emitter_buffer),
        disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize)),
        disposer_outflow_deletion(disposer, 0),
        free_stream_surface_indicator(water_block_complex),
        update_fluid_density(water_block_complex),
        get_fluid_advection_time_step_size(water_block, U_f),
        get_fluid_time_step_size(water_block),
        velocity_boundary_condition_constraint(water_block),
        pressure_relaxation(water_block_complex),
        density_relaxation(water_block_complex),
        viscous_acceleration(water_block_complex),
        transport_velocity_correction(water_block_complex),
        compute_vorticity(water_block_inner),
        fish_body_normal_direction(fish_body),
        fish_body_corrected_configuration(fish_inner),
        viscous_force_on_solid(fish_contact),
        fluid_force_on_fish_update(fish_contact, viscous_force_on_solid),
        average_velocity_and_acceleration(fish_body),
        composite_material_id(fish_body),
        imposing_active_strain(fish_body),
        fish_body_computing_time_step_size(fish_body),
        fish_body_stress_relaxation_first_half(fish_inner),
        fish_body_stress_relaxation_second_half(fish_inner),
        fish_body_update_normal(fish_body),
        write_real_body_states(io_environment, system.real_bodies_),
        write_total_viscous_force_on_insert_body(io_environment, viscous_force_on_solid, "TotalViscousForceOnSolid"),
        write_total_force_on_insert_body(io_environment, fluid_force_on_fish_update, "TotalForceOnSolid"),
        pressure_probe("Pressure", io_environment, fish_observer_contact_water),
        fish_velocity_probe("Velocity", io_environment, fish_observer_contact_solid),
        fish_position_probe("Position", io_environment, fish_observer_contact_solid),
        interpolation_particle_position(fish_observer_contact_solid, "Position", "Position")
    {
        GlobalStaticVariables::physical_time_ = 0.0;
        /** We can output a method-specific particle data for debug */
        water_block.addBodyStateForRecording<Real>("Pressure");
        water_block.addBodyStateForRecording<int>("Indicator");

        /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
        pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);

        fish_body.addBodyStateForRecording<Real>("Density");
        fish_body.addBodyStateForRecording<int>("MaterialID");
        fish_body.addBodyStateForRecording<Matd>("ActiveStrain");

        system.initializeSystemCellLinkedLists();
        system.initializeSystemConfigurations();
        fish_body_normal_direction.exec();
        fish_body_corrected_configuration.exec();
        composite_material_id.exec();

        //----------------------------------------------------------------------
        //	First output before the main loop.
        //----------------------------------------------------------------------
        write_real_body_states.writeToFile();
        pressure_probe.writeToFile();
        fish_velocity_probe.writeToFile();
        fish_position_probe.writeToFile();
    }

    virtual ~SphFishSimulation(){};

    int cmakeTest()
    {
        return 1;
    }
    //----------------------------------------------------------------------
    //	Action from python environment.
    //----------------------------------------------------------------------
    void setLambdafromPython(Real lambda_from_python)
    {
        lambda_ = lambda_from_python;
    }

    void setFreqfromPython(Real freq_from_python)
    {
        frequency_ = freq_from_python;
    }
    //----------------------------------------------------------------------
    //	head position for python environment.
    //----------------------------------------------------------------------
    Real getFishHeadPositionX()
    {
        return fish_position_probe.getCurrentPositionX(18);
    }

    Real getFishHeadPositionY()
    {
        return fish_position_probe.getCurrentPositionY(18);
    }
    //----------------------------------------------------------------------
    //	observation for python environment.
    //----------------------------------------------------------------------
    Real getPressurePoint1()
    {
        return pressure_probe.getCurrentPressure(0);
    }

    Real getPressurePoint2()
    {
        return pressure_probe.getCurrentPressure(1);
    }

    Real getPressurePoint3()
    {
        return pressure_probe.getCurrentPressure(2);
    }

    Real getPressurePoint4()
    {
        return pressure_probe.getCurrentPressure(3);
    }

    Real getPressurePoint5()
    {
        return pressure_probe.getCurrentPressure(4);
    }

    Real getPressurePoint6()
    {
        return pressure_probe.getCurrentPressure(5);
    }

    Real getPressurePoint7()
    {
        return pressure_probe.getCurrentPressure(6);
    }

    Real getPressurePoint8()
    {
        return pressure_probe.getCurrentPressure(7);
    }

    Real getPressurePoint9()
    {
        return pressure_probe.getCurrentPressure(8);
    }

    Real getPressurePoint11()
    {
        return pressure_probe.getCurrentPressure(9);
    }

    Real getPressurePoint12()
    {
        return pressure_probe.getCurrentPressure(10);
    }

    Real getPressurePoint13()
    {
        return pressure_probe.getCurrentPressure(11);
    }

    Real getPressurePoint14()
    {
        return pressure_probe.getCurrentPressure(12);
    }

    Real getPressurePoint15()
    {
        return pressure_probe.getCurrentPressure(13);
    }

    Real getPressurePoint16()
    {
        return pressure_probe.getCurrentPressure(14);
    }

    Real getPressurePoint17()
    {
        return pressure_probe.getCurrentPressure(15);
    }

    Real getPressurePoint18()
    {
        return pressure_probe.getCurrentPressure(16);
    }

    Real getPressurePoint19()
    {
        return pressure_probe.getCurrentPressure(17);
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    void runCase(int episode, Real pause_time_from_python)
    {
        while (GlobalStaticVariables::physical_time_ < pause_time_from_python)
        {
            Real integration_time = 0.0;

            /** Integrate time (loop) until the next output time. */
            while (integration_time < D_Time)
            {
                initialize_a_fluid_step.exec();
                Real Dt = get_fluid_advection_time_step_size.exec();
                free_stream_surface_indicator.exec();
                update_fluid_density.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                /** FSI for viscous force. */
                viscous_force_on_solid.exec();
                /** Update normal direction on elastic body.*/
                fish_body_update_normal.exec();
                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt = get_fluid_time_step_size.exec();
                    /** Fluid pressure relaxation, first half. */
                    pressure_relaxation.exec(dt);
                    /** FSI for fluid force on solid body. */
                    fluid_force_on_fish_update.exec();
                    /** Fluid pressure relaxation, second half. */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    inner_ite_dt_s = 0;
                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        Real dt_s = SMIN(fish_body_computing_time_step_size.exec(), dt - dt_s_sum);
                        imposing_active_strain.exec();
                        fish_body_stress_relaxation_first_half.exec(dt_s);
                        fish_body_stress_relaxation_second_half.exec(dt_s);
                        dt_s_sum += dt_s;
                        inner_ite_dt_s++;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);
                    interpolation_particle_position.exec();

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                    emitter_buffer_inflow_condition.exec(dt);
                    inner_ite_dt++;
                }
                pressure_probe.writeToFile(number_of_iterations);
                fish_velocity_probe.writeToFile(number_of_iterations);
                write_total_viscous_force_on_insert_body.writeToFile(number_of_iterations);

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                            << GlobalStaticVariables::physical_time_
                            << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
                }
                number_of_iterations++;

                /** Water block configuration and periodic condition. */
                emitter_inflow_injection.exec();
                disposer_outflow_deletion.exec();

                water_block.updateCellLinkedListWithParticleSort(100);
                fish_body.updateCellLinkedList();
                /** one need update configuration after periodic condition. */
                water_block_complex.updateConfiguration();
                /** one need update configuration after periodic condition. */
                fish_contact.updateConfiguration();
                fish_observer_contact_water.updateConfiguration();
                //fish_observer_contact_solid.updateConfiguration();
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            compute_vorticity.exec();
            if (episode % 10 == 1 || episode >= 1000)
                write_real_body_states.writeToFile();

            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
    }

};


PYBIND11_MODULE(test_2d_flow_stream_around_fish_pybind, m)
{
	py::class_<SphFishSimulation>(m, "from_sph_cpp")
		.def(py::init<const int&>())
        .def("CmakeTest", &SphFishSimulation::cmakeTest)
        .def("SetLambda", &SphFishSimulation::setLambdafromPython)
        .def("SetFreq", &SphFishSimulation::setFreqfromPython)
        .def("GetFishHeadPositionX", &SphFishSimulation::getFishHeadPositionX)
        .def("GetFishHeadPositionY", &SphFishSimulation::getFishHeadPositionY)
        .def("GetPressurePoint1", &SphFishSimulation::getPressurePoint1)
        .def("GetPressurePoint2", &SphFishSimulation::getPressurePoint2)
        .def("GetPressurePoint3", &SphFishSimulation::getPressurePoint3)
        .def("GetPressurePoint4", &SphFishSimulation::getPressurePoint4)
        .def("GetPressurePoint5", &SphFishSimulation::getPressurePoint5)
        .def("GetPressurePoint6", &SphFishSimulation::getPressurePoint6)
        .def("GetPressurePoint7", &SphFishSimulation::getPressurePoint7)
        .def("GetPressurePoint8", &SphFishSimulation::getPressurePoint8)
        .def("GetPressurePoint9", &SphFishSimulation::getPressurePoint9)
        .def("GetPressurePoint11", &SphFishSimulation::getPressurePoint11)
        .def("GetPressurePoint12", &SphFishSimulation::getPressurePoint12)
        .def("GetPressurePoint13", &SphFishSimulation::getPressurePoint13)
        .def("GetPressurePoint14", &SphFishSimulation::getPressurePoint14)
        .def("GetPressurePoint15", &SphFishSimulation::getPressurePoint15)
        .def("GetPressurePoint16", &SphFishSimulation::getPressurePoint16)
        .def("GetPressurePoint17", &SphFishSimulation::getPressurePoint17)
        .def("GetPressurePoint18", &SphFishSimulation::getPressurePoint18)
        .def("GetPressurePoint19", &SphFishSimulation::getPressurePoint19)
		.def("RunCase", &SphFishSimulation::runCase);
}