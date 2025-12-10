/**
 * @file	dambreak_python.cpp
 * @brief	2D dambreak example with pybind11 for python.
 * @details	This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"         //SPHinXsys Library.
#include <pybind11/pybind11.h> //pybind11 Library.
namespace py = pybind11;
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
class Parameter
{
  protected:
    Real DL = 5.366;                    /**< Water tank length. */
    Real DH = 5.366;                    /**< Water tank height. */
    Real LL = 2.0;                      /**< Water column length. */
    Real LH = 1.0;                      /**< Water column height. */
    Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
    Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
    //----------------------------------------------------------------------
    //	Material parameters.
    //----------------------------------------------------------------------
    Real rho0_f = 1.0;                       /**< Reference density of fluid. */
    Real gravity_g = 1.0;                    /**< Gravity. */
    Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
    Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
    //----------------------------------------------------------------------
    //	Geometric shapes used in this case.
    //----------------------------------------------------------------------
    Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
    Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
    Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
    Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
    Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
    Vec2d inner_wall_translation = inner_wall_halfsize;
};
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape, public Parameter
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//  Define system, geometry, material, particles and all other things.
//----------------------------------------------------------------------
class PreSettingCase : public Parameter
{
  protected:
    BoundingBoxd system_domain_bounds;
    SPHSystem sph_system;
    FluidBody water_block;
    SolidBody wall_boundary;
    StdVec<Vecd> observation_location;
    ObserverBody fluid_observer;

  public:
    PreSettingCase()
        : system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW)),
          sph_system(system_domain_bounds, particle_spacing_ref),
          water_block(sph_system, makeShared<GeometricShapeBox>(
                                      Transform(water_block_translation), water_block_halfsize, "WaterBody")),
          wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary")),
          observation_location({Vecd(DL, 0.2)}),
          fluid_observer(sph_system, "FluidObserver")
    {
        //----------------------------------------------------------------------
        //	Creating bodies with corresponding materials and particles.
        //----------------------------------------------------------------------
        water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
        water_block.generateParticles<BaseParticles, Lattice>();

        wall_boundary.defineMaterial<Solid>();
        wall_boundary.generateParticles<BaseParticles, Lattice>();

        fluid_observer.generateParticles<ObserverParticles>(observation_location);
    }
};
//----------------------------------------------------------------------
//  Define environment.
//----------------------------------------------------------------------
class Environment : public PreSettingCase
{
  protected:
    SPHSystem &sph_system_;
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner;
    ContactRelation water_wall_contact;
    ComplexRelation water_block_complex;
    ContactRelation fluid_observer_contact;
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity;
    SimpleDynamics<GravityForce<Gravity>> constant_gravity;
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction;

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation;
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation;

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> fluid_advection_time_step;
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step;
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting<ParallelPolicy> particle_sorting;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording;
    RestartIO restart_io;
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_water_mechanical_energy;
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real output_interval = 0.1;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

  public:
    explicit Environment(int set_restart_step)
        : PreSettingCase(),
          sph_system_(sph_system),
          water_block_inner(water_block),
          water_wall_contact(water_block, {&wall_boundary}),
          water_block_complex(water_block_inner, water_wall_contact),
          fluid_observer_contact(fluid_observer, {&water_block}),
          gravity(Vecd(0.0, -gravity_g)),
          constant_gravity(water_block, gravity),
          wall_boundary_normal_direction(wall_boundary),
          fluid_pressure_relaxation(water_block_inner, water_wall_contact),
          fluid_density_relaxation(water_block_inner, water_wall_contact),
          fluid_density_by_summation(water_block_inner, water_wall_contact),
          fluid_advection_time_step(water_block, U_ref),
          fluid_acoustic_time_step(water_block),
          particle_sorting(water_block),
          body_states_recording(sph_system),
          restart_io(sph_system),
          write_water_mechanical_energy(water_block, gravity),
          write_recorded_water_pressure("Pressure", fluid_observer_contact)
    {
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        wall_boundary_normal_direction.exec();
        constant_gravity.exec();
        /** set restart step. */
        sph_system.setRestartStep(set_restart_step);
        //----------------------------------------------------------------------
        //	First output before the main loop.
        //----------------------------------------------------------------------
        body_states_recording.writeToFile();
        write_water_mechanical_energy.writeToFile(sph_system.RestartStep());
        write_recorded_water_pressure.writeToFile(sph_system.RestartStep());
    }

    virtual ~Environment() {};
    //----------------------------------------------------------------------
    //	For ctest.
    //----------------------------------------------------------------------
    int cmakeTest()
    {
        return 1;
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    void runCase(Real End_time)
    {
        Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
        size_t number_of_iterations = 0;
        while (physical_time < End_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                /** outer loop for dual-time criteria time-stepping. */
                time_instance = TickCount::now();
                Real advection_dt = fluid_advection_time_step.exec();
                fluid_density_by_summation.exec();
                interval_computing_time_step += TickCount::now() - time_instance;

                time_instance = TickCount::now();
                Real relaxation_time = 0.0;
                Real acoustic_dt = 0.0;
                while (relaxation_time < advection_dt)
                {
                    /** inner loop for dual-time criteria time-stepping.  */
                    acoustic_dt = fluid_acoustic_time_step.exec();
                    fluid_pressure_relaxation.exec(acoustic_dt);
                    fluid_density_relaxation.exec(acoustic_dt);
                    relaxation_time += acoustic_dt;
                    integration_time += acoustic_dt;
                    physical_time += acoustic_dt;
                }
                interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

                /** screen output, write body reduced values and restart files  */
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time
                              << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                    if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                    {
                        write_water_mechanical_energy.writeToFile(number_of_iterations);
                        write_recorded_water_pressure.writeToFile(number_of_iterations);
                    }
                    if (number_of_iterations % restart_output_interval == 0)
                        restart_io.writeToFile(number_of_iterations);
                }
                number_of_iterations++;

                /** Update cell linked list and configuration. */
                time_instance = TickCount::now();
                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                water_block.updateCellLinkedList();
                water_block_complex.updateConfiguration();
                fluid_observer_contact.updateConfiguration();
                interval_updating_configuration += TickCount::now() - time_instance;
            }

            body_states_recording.writeToFile();
            TickCount t2 = TickCount::now();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TickCount::interval_t tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds()
                  << " seconds." << std::endl;
        std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
                  << interval_computing_time_step.seconds() << "\n";
        std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
                  << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
        std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
                  << interval_updating_configuration.seconds() << "\n";

        if (sph_system.GenerateRegressionData())
        {
            write_water_mechanical_energy.generateDataBase(1.0e-3);
            write_recorded_water_pressure.generateDataBase(1.0e-3);
        }

        else if (sph_system.RestartStep() == 0)
        {
            write_water_mechanical_energy.testResult();
            write_recorded_water_pressure.testResult();
        }
    }
};
//----------------------------------------------------------------------
//	Use pybind11 to expose.
//----------------------------------------------------------------------
/** test_2d_dambreak_python should be same with the project name */
PYBIND11_MODULE(test_2d_dambreak_python, m)
{
    py::class_<Environment>(m, "dambreak_from_sph_cpp")
        .def(py::init<const int &>())
        .def("CmakeTest", &Environment::cmakeTest)
        .def("RunCase", &Environment::runCase);
}
