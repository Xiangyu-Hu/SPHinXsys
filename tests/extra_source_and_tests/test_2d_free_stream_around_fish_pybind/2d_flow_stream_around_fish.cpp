/**
 * @file 2d_flow_stream_around_fish.cpp
 * @brief fish swimming driven by active muscles
 * @author Mai Ye, Yaru Ren and Xiangyu Hu
 */
#include <pybind11/pybind11.h>
#include "2d_flow_stream_around_fish.h"
#include "sphinxsys.h"

using namespace SPH;
namespace py = pybind11;

class SphBasicSystemSetting : public SphBasicGeometrySetting
{
protected:
    BoundingBox system_domain_bounds;
    SPHSystem sph_system;
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer;

  public:
    SphBasicSystemSetting(int episode_env) :
        system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW)),
        sph_system(system_domain_bounds, particle_spacing_ref), inlet_particle_buffer(0.5) {
        sph_system.setIOEnvironment();
    }
};

class SphFishRelaxationEnvironment : public SphBasicSystemSetting
{
protected:
    FluidBody water_block;
    SolidBody fish_body;
    ObserverBody fish_observer;
    

public:
    SphFishRelaxationEnvironment(int episode_env) :
        SphBasicSystemSetting(episode_env),
        water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
        fish_body(sph_system, makeShared<FishBody>("FishBody")),
        fish_observer(sph_system, "FishObserver")
    {

        water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

        fish_body.defineAdaptationRatios(1.15, 2.0);
        fish_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
        fish_body.defineMaterial<FishBodyComposite>();
        fish_body.generateParticles<BaseParticles, Lattice>();
        fish_observer.generateParticles<BaseParticles, Observer>(observation_locations);
    }
};

class SphFishRelaxation : public SphFishRelaxationEnvironment
{
protected:
    InnerRelation fish_inner, water_block_inner;
    ContactRelation water_block_contact, fish_contact, fish_observer_contact_water, fish_observer_contact_fish;
    ComplexRelation water_block_complex;
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_fish_body_particles;
    BodyStatesRecordingToVtp write_fish_body;
    ReloadParticleIO write_fish_particle_reload_files;
    relax_dynamics::RelaxationStepInner relaxation_step_inner_fish;
    int ite_p = 0;

public:
    explicit SphFishRelaxation(int episode_env) :
        SphFishRelaxationEnvironment(episode_env),
        fish_inner(fish_body),
        water_block_inner(water_block),
        water_block_contact(water_block, { &fish_body }),
        fish_contact(fish_body, { &water_block }),
        fish_observer_contact_water(fish_observer, { &water_block }),
        fish_observer_contact_fish(fish_observer, { &fish_body }),
        water_block_complex(water_block_inner, water_block_contact),
        random_fish_body_particles(fish_body),
        write_fish_body(fish_body),
        write_fish_particle_reload_files(fish_body),
        relaxation_step_inner_fish(fish_inner)
    {
        random_fish_body_particles.exec(0.25);
        relaxation_step_inner_fish.SurfaceBounding().exec();
        write_fish_body.writeToFile();
        while (ite_p < 1000)
        {
            relaxation_step_inner_fish.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the fish body N = " << ite_p << "\n";
                write_fish_body.writeToFile(ite_p);
            }
        }
        write_fish_particle_reload_files.writeToFile();

        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    }
};

class SphFishReloadTrainEnvironment : public SphBasicSystemSetting
{
protected:
    FluidBody water_block;
    SolidBody fish_body;
    ObserverBody fish_observer;

public:
    SphFishReloadTrainEnvironment(int episode_env) :
        SphBasicSystemSetting(episode_env),
        water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
        fish_body(sph_system, makeShared<FishBody>("FishBody")),
        fish_observer(sph_system, "FishObserver")
    {
        water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

        fish_body.defineAdaptationRatios(1.15, 2.0);
        fish_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
        fish_body.defineMaterial<FishBodyComposite>();
        fish_body.generateParticles<BaseParticles, Reload>(fish_body.getName());

        fish_observer.generateParticles<BaseParticles, Observer>(observation_locations);
    }
};

class SphFishReloadTrain : public SphFishReloadTrainEnvironment
{
protected:
    InnerRelation fish_inner, water_block_inner;
    ContactRelation water_block_contact, fish_contact, fish_observer_contact_water, fish_observer_contact_fish;
    ComplexRelation water_block_complex;
    
    SimpleDynamics<NormalDirectionFromBodyShape> fish_body_normal_direction;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> fish_body_corrected_configuration;
    SimpleDynamics<FishMaterialInitialization> composite_material_id_fish;
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> fish_body_stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> fish_body_stress_relaxation_second_half;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> fish_body_computing_time_step_size;
    SimpleDynamics<ImposingActiveStrain> imposing_active_strain;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> fish_body_update_normal;
    TimeDependentAcceleration time_dependent_acceleration;
    SimpleDynamics<GravityForce> apply_gravity_force;
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation;
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density;
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration;
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction;
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size;
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size;
    BodyAlignedBoxByParticle emitter;
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection;
    BodyAlignedBoxByCell emitter_buffer, disposer;
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition;
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion;
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint;
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity;
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_fish;
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> fluid_force_on_fish_update;

    solid_dynamics::AverageVelocityAndAcceleration average_fish_velocity_and_acceleration;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states;
    RestartIO restart_io;
    ObservedQuantityRecording<Real> fish_pressure_probe;
    ObservedQuantityRecording<Vecd> fish_velocity_probe;
    ObservedQuantityRecording<Vecd> fish_position_probe;
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int write_episode_vtp = 10;
    Real D_Time = 0.01;  /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

public:
    explicit SphFishReloadTrain(int episode_env) :
        SphFishReloadTrainEnvironment(episode_env),
        fish_inner(fish_body), water_block_inner(water_block),
        water_block_contact(water_block, { &fish_body }),
        fish_contact(fish_body, { &water_block }),
        fish_observer_contact_water(fish_observer, { &water_block }),
        fish_observer_contact_fish(fish_observer, { &fish_body }),
        water_block_complex(water_block_inner, water_block_contact),
        time_dependent_acceleration(Vec2d::Zero()),
        apply_gravity_force(water_block, time_dependent_acceleration),
        emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize)),
        emitter_inflow_injection(emitter, inlet_particle_buffer, xAxis),
        emitter_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize)),
        emitter_buffer_inflow_condition(emitter_buffer),
        disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize)),
        disposer_outflow_deletion(disposer, 0),
        free_stream_surface_indicator(water_block_inner, water_block_contact),
        update_fluid_density(water_block_inner, water_block_contact),
        get_fluid_advection_time_step_size(water_block, U_f),
        get_fluid_time_step_size(water_block),
        velocity_boundary_condition_constraint(water_block),
        pressure_relaxation(water_block_inner, water_block_contact),
        density_relaxation(water_block_inner, water_block_contact),
        viscous_acceleration(water_block_inner, water_block_contact),
        transport_velocity_correction(water_block_inner, water_block_contact),
        compute_vorticity(water_block_inner),
        fish_body_normal_direction(fish_body),
        fish_body_corrected_configuration(fish_inner),
        viscous_force_on_fish(fish_contact),
        fluid_force_on_fish_update(fish_contact),
        average_fish_velocity_and_acceleration(fish_body),
        composite_material_id_fish(fish_body),
        imposing_active_strain(fish_body),
        fish_body_computing_time_step_size(fish_body),
        fish_body_stress_relaxation_first_half(fish_inner),
        fish_body_stress_relaxation_second_half(fish_inner),
        fish_body_update_normal(fish_body),
        write_real_body_states(sph_system),
        restart_io(sph_system),
        fish_pressure_probe("Pressure", fish_observer_contact_water),
        fish_velocity_probe("Velocity", fish_observer_contact_fish),
        fish_position_probe("Position", fish_observer_contact_fish)
    {
        GlobalStaticVariables::physical_time_ = 0.0;
        write_real_body_states.addVariableRecording<Real>(water_block, "Pressure");
        write_real_body_states.addVariableRecording<int>(water_block, "Indicator");
        write_real_body_states.addVariableRecording<Real>(fish_body, "Density");
        write_real_body_states.addVariableRecording<int>(fish_body, "MaterialID");
        write_real_body_states.addVariableRecording<Matd>(fish_body, "ActiveStrain");
        /** initialize cell linked lists for all bodies. */
        sph_system.initializeSystemCellLinkedLists();
        /** initialize configurations for all bodies. */
        sph_system.initializeSystemConfigurations();
        /** computing surface normal direction for the fish. */
        fish_body_normal_direction.exec();
        /** computing linear reproducing configuration for the fish. */
        fish_body_corrected_configuration.exec();
        /** initialize material ids for the fish. */
        composite_material_id_fish.exec();

        write_real_body_states.writeToFile();
        fish_pressure_probe.writeToFile();
        fish_velocity_probe.writeToFile();
        fish_position_probe.writeToFile();
    }

    virtual ~SphFishReloadTrain() {};

    //----------------------------------------------------------------------
    //	Action from python environment.
    //----------------------------------------------------------------------


    void setFreqfromPython(Real freq_from_python)
    {
        imposing_active_strain.setFreqfromPython(freq_from_python);
    }

    //----------------------------------------------------------------------
    //      	observation for python environment.
    //----------------------------------------------------------------------
    Real getFishPressure(int number)
    {
        StdLargeVec<Real> *observedQuantity = 
            fish_pressure_probe.getObservedQuantity();
        return (*observedQuantity)[number];
    }

    Real getFishPositionX(int number)
    {
        StdLargeVec<Vecd> *observedQuantity = fish_position_probe.getObservedQuantity();
        return (*observedQuantity)[number][0];
    }

    Real getFishPositionY(int number)
    {
        StdLargeVec<Vecd> *observedQuantity = fish_position_probe.getObservedQuantity();
        return (*observedQuantity)[number][1];
    }

    void runCase(int episode, Real pause_time_from_python)
    {
        //----------------------------------------------------------------------
        //	Main loop starts here.
        //----------------------------------------------------------------------
        while (GlobalStaticVariables::physical_time_ < pause_time_from_python)
        {
            Real integration_time = 0.0;

            /** Integrate time (loop) until the next output time. */
            while (integration_time < D_Time)
            {
                apply_gravity_force.exec();
                Real Dt = get_fluid_advection_time_step_size.exec();
                free_stream_surface_indicator.exec();
                update_fluid_density.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                /** FSI for viscous force. */
                viscous_force_on_fish.exec();
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
                    average_fish_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        Real dt_s = SMIN(fish_body_computing_time_step_size.exec(), dt - dt_s_sum);
                        imposing_active_strain.exec();
                        fish_body_stress_relaxation_first_half.exec(dt_s);
                        fish_body_stress_relaxation_second_half.exec(dt_s);
                        dt_s_sum += dt_s;
                        inner_ite_dt_s++;
                    }
                    average_fish_velocity_and_acceleration.update_averages_.exec(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                    emitter_buffer_inflow_condition.exec(dt);
                    inner_ite_dt++;
                }

                fish_pressure_probe.writeToFile(number_of_iterations);
                fish_velocity_probe.writeToFile(number_of_iterations);
                fish_position_probe.writeToFile(number_of_iterations);


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
                /** one need update configuration after periodic condition. */
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            compute_vorticity.exec();
            //write_real_body_states.writeToFile(number_of_iterations);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
    }
};


PYBIND11_MODULE(test_2d_free_stream_around_fish_pybind, m)
{
    py::class_<SphFishRelaxation>(m, "from_sph_relaxation")
        .def(py::init<const int&>());

    py::class_<SphFishReloadTrain>(m, "from_sph_reload_and_train")
        .def(py::init<const int&>())
        .def("SetFreq", &SphFishReloadTrain::setFreqfromPython)
        .def("GetFishPressurePoint", &SphFishReloadTrain::getFishPressure)
        .def("GetFishPositionX", &SphFishReloadTrain::getFishPositionX)
        .def("GetFishPositionY", &SphFishReloadTrain::getFishPositionY)
        .def("RunCase", &SphFishReloadTrain::runCase);
}
