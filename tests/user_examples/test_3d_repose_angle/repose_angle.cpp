#include "all_continuum.h"
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;
// general parameters for geometry
Real radius = 0.1;                                         // liquid length
Real height = 0.1;                                         // liquid height
Real resolution_ref = radius / 10;                         // particle spacing
Real BW = resolution_ref * 4;                              // boundary width
Real DL = 2 * radius * (1 + 1.24 * height / radius) + 0.1; // tank length
Real DH = height + 0.02;                                   // tank height
Real DW = DL;                                              // tank width
// for material properties
Real rho0_s = 2600;           /**< Reference density of soil. */
Real gravity_g = 9.8;         /**< Gravity force of soil. */
Real Youngs_modulus = 5.98e6; // reference Youngs modulus
Real poisson = 0.3;           // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3 * (1 - 2 * poisson)));
Real friction_angle = 30 * Pi / 180;
/** Define the soil body. */
Real inner_circle_radius = radius;
int resolution(20);
class SoilBlock : public ComplexShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(DL / 2, 0.5 * height, DW / 2);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), inner_circle_radius,
                                       0.5 * height, resolution, translation_column);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd outer_wall_halfsize = Vecd(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd outer_wall_translation = Vecd(-BW, -BW, -BW) + outer_wall_halfsize;
        Vecd inner_wall_halfsize = Vecd(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Vecd inner_wall_translation = inner_wall_halfsize;
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class SoilInitialCondition : public continuum_dynamics::ContinuumInitialCondition
{
  public:
    explicit SoilInitialCondition(RealBody &granular_column)
        : continuum_dynamics::ContinuumInitialCondition(granular_column){};

  protected:
    void update(size_t index_i, Real dt)
    {
        /** initial stress */
        Real y = pos_[index_i][1];
        Real gama = 1 - sin(friction_angle);
        Real stress_yy = -rho0_s * gravity_g * y;
        stress_tensor_3D_[index_i](1, 1) = stress_yy;
        stress_tensor_3D_[index_i](0, 0) = stress_yy * gama;
        stress_tensor_3D_[index_i](2, 2) = stress_yy * gama;
    };
};
// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    RealBody soil_block(sph_system, makeShared<SoilBlock>("GranularBody"));
    soil_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    soil_block.defineParticlesAndMaterial<PlasticContinuumParticles, PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? soil_block.generateParticles<ParticleGeneratorReload>(io_environment, soil_block.getName())
        : soil_block.generateParticles<ParticleGeneratorLattice>();
    soil_block.addBodyStateForRecording<Real>("Pressure");
    soil_block.addBodyStateForRecording<Real>("Density");
    soil_block.addBodyStateForRecording<Real>("VerticalStress");
    soil_block.addBodyStateForRecording<Real>("AccDeviatoricPlasticStrain");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    ComplexRelation soil_block_complex(soil_block, {&wall_boundary});
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    // run particle relaxation
    if (sph_system.RunParticleRelaxation())
    {
        /**
         * @brief 	Methods used for particle relaxation.
         */
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_column_particles(soil_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_column_to_vtp(io_environment, soil_block);
        /** Write the particle reload files. */

        ReloadParticleIO write_particle_reload_files(io_environment, soil_block);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(soil_block_complex.getInnerRelation());
        /**
         * @brief 	Particle relaxation starts here.
         */
        random_column_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        body_states_recording.writeToFile(0.0);

        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
                write_column_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0.0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<SoilInitialCondition> soil_initial_condition(soil_block);
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<TimeStepInitialization> soil_step_initialization(soil_block, gravity_ptr);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> soil_acoustic_time_step(soil_block, 0.1);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> soil_density_by_summation(soil_block_complex);
    InteractionDynamics<continuum_dynamics::StressDiffusion> stress_diffusion(soil_block_complex.getInnerRelation());
    Dynamics1Level<continuum_dynamics::StressRelaxation1stHalfRiemannWithWall> granular_stress_relaxation_1st(soil_block_complex);
    Dynamics1Level<continuum_dynamics::StressRelaxation2ndHalfRiemannWithWall> granular_stress_relaxation_2nd(soil_block_complex);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    RestartIO restart_io(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_soil_mechanical_energy(io_environment, soil_block, gravity_ptr);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    soil_initial_condition.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        soil_block.updateCellLinkedList();
        soil_block_complex.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 500;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real End_Time = 0.5;         /**< End time. */
    Real D_Time = End_Time / 25; /**< Time stamps for output of body states. */
    Real Dt = 0.1 * D_Time;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_soil_stress_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_soil_mechanical_energy.writeToFile(number_of_iterations);
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
            soil_step_initialization.exec();

            soil_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = soil_acoustic_time_step.exec();

                granular_stress_relaxation_1st.exec(dt);
                stress_diffusion.exec();
                granular_stress_relaxation_2nd.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                interval_computing_soil_stress_relaxation += TickCount::now() - time_instance;

                /** screen output, write body reduced values and restart files  */
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4) << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << std::scientific << "	dt = " << dt << "\n";

                    if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                    {
                        write_soil_mechanical_energy.writeToFile(number_of_iterations);
                    }
                    if (number_of_iterations % restart_output_interval == 0)
                        restart_io.writeToFile(number_of_iterations);
                }
                number_of_iterations++;
                soil_block.updateCellLinkedList();
                soil_block_complex.updateConfiguration();
            }
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
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
    std::cout << std::fixed << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_soil_stress_relaxation = "
              << interval_computing_soil_stress_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    std::cout << "total time steps = " << number_of_iterations << "\n";

    // sph_system.GenerateRegressionData() = true;
    if (sph_system.GenerateRegressionData())
    {
        write_soil_mechanical_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_soil_mechanical_energy.testResult();
    }

    return 0;
}
