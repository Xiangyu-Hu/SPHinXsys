/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
#include <pybind11/pybind11.h> //pybind11 Library.
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
class BasicParameters
{
  public:
    BasicParameters()= default;
    ~BasicParameters()= default;
    void setSimDoamin(Real lower_x, Real lower_y, Real lower_z, Real upper_x, Real upper_y, Real upper_z)
    {
        sim_Domain_Lower = Vecd(lower_x, lower_y, lower_z);
        sim_Domain_Upper = Vecd(upper_x, upper_y, upper_z);
    }
//    void setSimDoaminLower(Real lower_x, Real lower_y, Real lower_z)
//    {
//        sim_Domain_Lower = Vecd(lower_x, lower_y, lower_z);
//    }
    void setSimDoaminLower(const Vecd &lower)
    {
        sim_Domain_Lower = lower;
    }
//    void setSimDoaminUpper(Real upper_x, Real upper_y, Real upper_z)
//    {
//        sim_Domain_Upper = Vecd(upper_x, upper_y, upper_z);
//    }
    void setSimDoaminUpper(const Vecd &upper)
    {
        sim_Domain_Upper = upper;
    }
    void setParticleSpacingRef(Real particle_spacing)
    {
        particle_spacing_ref = particle_spacing;
    }
    void setGravity(Real gravity)
    {
        gravity_g = gravity;
    }
    void setEndTime(Real time)
    {
        end_time = time;
    }
    void setRho(Real rho)
    {
        rho0_f = rho;
    }
    void setVelRef(Real vel_ref)
    {
        U_ref = vel_ref;
        c_f = 10 * U_ref;

        //    Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
    }
    void setWaterFilePath(const string &water_file_path)
    {
        water_block_file_path = water_file_path;
    }
    void setRigidFilePath(const string &rigid_file_path)
    {
        rigid_block_file_path = rigid_file_path;
    }
    void setElasticFilePath(const string &elastic_file_path)
    {
        elastic_block_file_path = elastic_file_path;
    }
    void setWaterFileName(const std::vector<string> &water_name)
    {
        water_block_file_name = water_name;
    }
    void setRigidFileName(const std::vector<string> &rigid_name)
    {
        rigid_block_file_name = rigid_name;
    }
    void setElasticFileName(const std::vector<string> &elastic_name)
    {
        elastic_block_file_name = elastic_name;
    }
    void pushbackWaterFileName(const string &water_name)
    {
        water_block_file_name.push_back(water_name);
    }
    void pushbackRigidFileName(const string &rigid_name)
    {
        rigid_block_file_name.push_back(rigid_name);
    }
    void pushbackElasticFileName(const string &elastic_name)
    {
        elastic_block_file_name.push_back(elastic_name);
    }

    Vecd getSimDoaminLower() const
    {
        return sim_Domain_Lower;
    }
    Vecd getSimDoaminUpper() const
    {
        return sim_Domain_Upper;
    }
    Real getParticleSpacingRef() const
    {
        return particle_spacing_ref;
    }
    Real getGravity() const
    {
        return gravity_g;
    }
    Real getEndTime() const
    {
        return end_time;
    }
    Real getRho() const
    {
        return rho0_f;
    }
    Real getVelRef() const
    {
        return U_ref;
    }
    Real getSoundRef() const
    {
        return c_f;
    }
    string getWaterFilePath() const
    {
        return water_block_file_path;
    }
    string getRigidFilePath() const
    {
        return rigid_block_file_path;
    }
    string getElasticFilePath() const
    {
        return elastic_block_file_path;
    }
    std::vector<string> getWaterFileName() const
    {
        return water_block_file_name;
    }
    std::vector<string> getRigidFileName() const
    {
        return rigid_block_file_name;
    }
    std::vector<string> getElasticFileName() const
    {
        return elastic_block_file_name;
    }
    string getTotalWaterFilePath(int index = 0) const
    {
        return water_block_file_path + water_block_file_name[index] + ".stl";
    }
    string getTotalRigidFilePath(int index = 0) const
    {
        return rigid_block_file_path + rigid_block_file_name[index] + ".stl";
    }
    string getTotalElasticFilePath(int index = 0) const
    {
        return elastic_block_file_path + elastic_block_file_name[index] + ".stl";
    }

  private:
    Vecd sim_Domain_Lower;                                       /**< Simulation domain lower. */
    Vecd sim_Domain_Upper;                                       /**< Simulation domain upper. */

    Real particle_spacing_ref = 0.05;                            /**< Initial reference particle spacing. */
    Real gravity_g = 1.0;                                        /**< Gravity. */
    Real end_time = 20;                                          /**< Simulation end time. */

    Real rho0_f = 1.0;                                           /**< Reference density of fluid. */
    Real U_ref = 2.0 * sqrt(gravity_g * 1);                   /**< Characteristic velocity. */
    Real c_f = 10.0 * U_ref;                                     /**< Reference sound speed. */

    std::string water_block_file_path = "./Mesh/Water/";         /**< Water block stl files path. */
    std::string rigid_block_file_path = "./Mesh/Rigid/";         /**< Rigid block stl files path. */
    std::string elastic_block_file_path = "./Mesh/Elastic/";     /**< Elastic block stl files path. */
    std::vector<string> water_block_file_name;                       /**< Water block stl files name. */
    std::vector<string> rigid_block_file_name;                       /**< Rigid block stl files name. */
    std::vector<string> elastic_block_file_name;                     /**< Elastic block stl files name. */
};

//	define the water block shape
class WaterBlockFromSTL : public ComplexShape
{
  public:
    explicit WaterBlockFromSTL(const std::string &shape_name,
                               const std::string &file_path = "",
                               const Vecd &translation = Vecd(0, 0, 0),
                               const Real &scaling = 1
                               ) : ComplexShape(shape_name)
    {
//        if(!filepath.empty())
//        {
//            add<TriangleMeshShapeSTL>(filepath, translation, scaling);
//        }
        add<TriangleMeshShapeSTL>(file_path, translation, scaling);
    }
};
//	define the static solid wall boundary shape
class RigidBlockFromSTL : public ComplexShape
{
  public:
    explicit RigidBlockFromSTL(const std::string &shape_name,
                               const std::string &file_path = "",
                               const Vecd &translation = Vecd(0, 0, 0),
                               const Real &scaling = 1
                               ) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(file_path, translation, scaling);

//        if(!(file_path_inner.empty() || file_path_outer.empty()))
//        {
//            add<TriangleMeshShapeSTL>(file_path_outer, translation, scaling);
//            subtract<TriangleMeshShapeSTL>(file_path_inner, translation, scaling);
//        }
    }
};

//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        // add observation points
        positions_.push_back(Vecd(1, 0.01, 0.5 * 0.5));
        positions_.push_back(Vecd(1, 0.1, 0.5 * 0.5));
        positions_.push_back(Vecd(1, 0.2, 0.5 * 0.5));
        positions_.push_back(Vecd(1, 0.24, 0.5 * 0.5));
        positions_.push_back(Vecd(1, 0.252, 0.5 * 0.5));
        positions_.push_back(Vecd(1, 0.266, 0.5 * 0.5));
    }
};

//----------------------------------------------------------------------
//  Define system, geometry, material, particles and all other things.
//----------------------------------------------------------------------
class SimulationConfig
{
  public:
    SimulationConfig(BasicParameters parameters)
        : basic_parameters(parameters),
          system_domain_bounds(parameters.getSimDoaminLower(), parameters.getSimDoaminUpper()),
          sph_system(system_domain_bounds, parameters.getParticleSpacingRef()),
          io_environment(sph_system),
          fluid_block(sph_system, makeShared<WaterBlockFromSTL>("WaterBody", parameters.getTotalWaterFilePath())),
          rigid_block(sph_system, makeShared<RigidBlockFromSTL>("Wall", parameters.getTotalRigidFilePath())),
          fluid_observer(sph_system, "FluidObserver")
    {
        //----------------------------------------------------------------------
        //	Creating bodies with corresponding materials and particles.
        //----------------------------------------------------------------------
        fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(parameters.getRho(), parameters.getSoundRef());
        fluid_block.generateParticles<ParticleGeneratorLattice>();

        rigid_block.defineParticlesAndMaterial<SolidParticles, Solid>();
        rigid_block.generateParticles<ParticleGeneratorLattice>();
        rigid_block.addBodyStateForRecording<Vec3d>("NormalDirection");

        fluid_observer.generateParticles<WaterObserverParticleGenerator>();
    }

    BoundingBox& getBoundingBox()
    {
        return system_domain_bounds;
    }
    SPHSystem& getSPHSystem()
    {
        return sph_system;
    }
    IOEnvironment& getIOEnvironment()
    {
        return io_environment;
    }
    FluidBody& getFluidBody()
    {
        return fluid_block;
    }
    SolidBody& getSolidBody()
    {
        return rigid_block;
    }
    ObserverBody& getObserverBody()
    {
        return fluid_observer;
    }
    BasicParameters& getBasicParameters()
    {
        return basic_parameters;
    }

  private:
    BasicParameters basic_parameters;
    BoundingBox system_domain_bounds;
    SPHSystem sph_system;
    IOEnvironment io_environment;
    FluidBody fluid_block;
    SolidBody rigid_block;
    ObserverBody fluid_observer;
};

class Simulator
{
  protected:
    SimulationConfig simulation_config;
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
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation;
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation;
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation;
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction;
    SharedPtr<Gravity> gravity_ptr;
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization;
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step;
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step;
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
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

  public:
    explicit Simulator(BasicParameters parameters)
        : simulation_config(parameters),
          water_block_inner(simulation_config.getFluidBody()),
          water_wall_contact(simulation_config.getFluidBody(), {&simulation_config.getSolidBody()}),
          water_block_complex(water_block_inner, water_wall_contact),
          fluid_observer_contact(simulation_config.getObserverBody(), {&simulation_config.getFluidBody()}),
          fluid_pressure_relaxation(water_block_inner, water_wall_contact),
          fluid_density_relaxation(water_block_inner, water_wall_contact),
          fluid_density_by_summation(water_block_inner, water_wall_contact),
          wall_boundary_normal_direction(simulation_config.getSolidBody()),
          gravity_ptr(makeShared<Gravity>(Vec3d(0.0, -parameters.getGravity(), 0.0))),
          fluid_step_initialization(simulation_config.getFluidBody(), gravity_ptr),
          fluid_advection_time_step(simulation_config.getFluidBody(), parameters.getVelRef()),
          fluid_acoustic_time_step(simulation_config.getFluidBody()),
          body_states_recording(simulation_config.getIOEnvironment(), simulation_config.getSPHSystem().real_bodies_),
          restart_io(simulation_config.getIOEnvironment(), simulation_config.getSPHSystem().real_bodies_),
          write_water_mechanical_energy(simulation_config.getIOEnvironment(), simulation_config.getFluidBody(), gravity_ptr),
          write_recorded_water_pressure("Pressure", simulation_config.getIOEnvironment(), fluid_observer_contact)
    {
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        simulation_config.getSPHSystem().initializeSystemCellLinkedLists();
        simulation_config.getSPHSystem().initializeSystemConfigurations();
        wall_boundary_normal_direction.exec();
        //----------------------------------------------------------------------
        //	First output before the main loop.
        //----------------------------------------------------------------------
        body_states_recording.writeToFile(0);
        write_water_mechanical_energy.writeToFile(0);
    }

    virtual ~Simulator(){};
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
    void run()
    {
        /** Set restart number of iterations. */
        size_t number_of_iterations = simulation_config.getSPHSystem().RestartStep();
        Real output_interval = simulation_config.getBasicParameters().getEndTime() / 160.0;
        Real dt = 0.0;
        while (GlobalStaticVariables::physical_time_ < simulation_config.getBasicParameters().getEndTime())
        {
            Real integration_time = 0.0;
            while (integration_time < output_interval)
            {
                fluid_step_initialization.exec();
                Real Dt = fluid_advection_time_step.exec();
                fluid_density_by_summation.exec();

                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    fluid_pressure_relaxation.exec(dt);
                    fluid_density_relaxation.exec(dt);
                    dt = fluid_acoustic_time_step.exec();
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

                simulation_config.getFluidBody().updateCellLinkedListWithParticleSort(100);
                water_block_complex.updateConfiguration();
                fluid_observer_contact.updateConfiguration();
                write_recorded_water_pressure.writeToFile(number_of_iterations);
            }

            write_water_mechanical_energy.writeToFile(number_of_iterations);

            TickCount t2 = TickCount::now();
            body_states_recording.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

        if (simulation_config.getSPHSystem().GenerateRegressionData())
        {
            write_water_mechanical_energy.generateDataBase(1.0e-3);
            write_recorded_water_pressure.generateDataBase(1.0e-3);
        }
        else
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
PYBIND11_MODULE(test_3d_blender, m)
{
    py::class_<Simulator> Simulator(m, "Simulator");
    Simulator.def(py::init<BasicParameters>())
        .def("run", &Simulator::run);

    py::class_<BasicParameters> BasicParameters(m, "BasicParameters");
    BasicParameters.def(py::init<>())
        .def_property("sim_domain_lower", &BasicParameters::getSimDoaminLower, &BasicParameters::setSimDoaminLower)
        .def_property("sim_domain_upper", &BasicParameters::getSimDoaminUpper, &BasicParameters::setSimDoaminUpper)
        .def_property("particle_spacing_ref", &BasicParameters::getParticleSpacingRef, &BasicParameters::setParticleSpacingRef)
        .def_property("gravity_g", &BasicParameters::getGravity, &BasicParameters::setGravity)
        .def_property("end_time", &BasicParameters::getEndTime, &BasicParameters::setEndTime)
        .def_property("rho0_f", &BasicParameters::getRho, &BasicParameters::setRho)
        .def_property("U_ref", &BasicParameters::getVelRef, &BasicParameters::setVelRef)
        .def_property_readonly("c_f", &BasicParameters::getSoundRef)
        .def_property("water_block_file_path", &BasicParameters::getWaterFilePath, &BasicParameters::setWaterFilePath)
        .def_property("rigid_block_file_path", &BasicParameters::getRigidFilePath, &BasicParameters::setRigidFilePath)
        .def_property("elastic_block_file_path", &BasicParameters::getElasticFilePath, &BasicParameters::setElasticFilePath)
        .def_property("water_block_file_name", &BasicParameters::getWaterFileName, &BasicParameters::setWaterFileName)
        .def_property("rigid_block_file_name", &BasicParameters::getRigidFileName, &BasicParameters::setRigidFileName)
        .def_property("elastic_block_file_name", &BasicParameters::getElasticFileName, &BasicParameters::setElasticFileName)
        .def("push_back_water", &BasicParameters::pushbackWaterFileName)
        .def("push_back_rigid", &BasicParameters::pushbackRigidFileName)
        .def("push_back_elastic", &BasicParameters::pushbackElasticFileName);
}
