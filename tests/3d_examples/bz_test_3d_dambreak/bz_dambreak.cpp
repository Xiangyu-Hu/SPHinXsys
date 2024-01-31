/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation.                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry
Real DL = 1.610;                // tank length
Real DH = 1.610;                // tank height
Real DW = 0.15;                 // tank width
Real LL = 0.6;                  // liquid length
Real LH = 0.3;                  // liquid height
Real LW = 0.15;                 // liquid width
Real resolution_ref = LH / 40 ; // particle spacing
Real BW = resolution_ref * 4;   // boundary width

// for material properties of the fluid
Real rho0_f = 1000;
Real gravity_g = 9.81;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

//	define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Vecd center_water(DL - 0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_water(center_water);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_wall(halfsize_inner);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
    }
};

//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        // add observation points
        positions_.push_back(Vecd(0.0, 0.003, 0.5 * DW));
        positions_.push_back(Vecd(0.0, 0.03, 0.5 * DW));
        positions_.push_back(Vecd(0.0, 0.08, 0.5 * DW));
    }
};

//WaveGauges 
Real h = 1.3 * resolution_ref;
Vecd halfsize_No1(h, 0.3, h);
Vecd translation_No1(1.310, 0.3, 0.0);
Vecd halfsize_No2(h, 0.3, h);
Vecd translation_No2(0.496, 0.3, 0.0);
Vecd halfsize_No3(h, 0.3, h);
Vecd translation_No3(0.2485, 0.3, 0.0);

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<Real>("Pressure");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<WaterObserverParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0));
    //InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_matrix(water_block_complex, 0.5);
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannConsistencyWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    water_block.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt write_water_block_states(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);

    /** WaveProbes. */
    BodyRegionByCell wave_probe_buffer_no_1(water_block, makeShared<TransformShape<GeometricShapeBox>>
        (Transform(translation_No1), halfsize_No1, "WaveProbe_01"));
    ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
        wave_probe_1(io_environment, wave_probe_buffer_no_1);

    BodyRegionByCell wave_probe_buffer_no_2(water_block, makeShared<TransformShape<GeometricShapeBox>>
        (Transform(translation_No2), halfsize_No2, "WaveProbe_02"));
    ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
        wave_probe_2(io_environment, wave_probe_buffer_no_2);

    BodyRegionByCell wave_probe_buffer_no_3(water_block, makeShared<TransformShape<GeometricShapeBox>>
        (Transform(translation_No3), halfsize_No3, "WaveProbe_03"));
    ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
        wave_probe_3(io_environment, wave_probe_buffer_no_3);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //kernel_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 2.0;
    Real output_interval = end_time / 50.0;
    Real dt = 0.0; // default acoustic time step sizes
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_water_block_states.writeToFile(0);
    write_water_mechanical_energy.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            //kernel_correction_matrix.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
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

            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            write_recorded_water_pressure.writeToFile(number_of_iterations); //recording the pressure at probe 1,2,3;
            wave_probe_1.writeToFile(); //recording the wave height at probe 1;
            wave_probe_2.writeToFile(); //recording the wave height at probe 2;
            wave_probe_3.writeToFile(); //recording the wave height at probe 3;
            write_water_mechanical_energy.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        write_water_block_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        //write_water_mechanical_energy.generateDataBase(1.0e-3);
        //write_recorded_water_pressure.generateDataBase(1.0e-3);
    }
    else
    {
        //write_water_mechanical_energy.testResult();
        //write_recorded_water_pressure.testResult();
    }


    return 0;
}
