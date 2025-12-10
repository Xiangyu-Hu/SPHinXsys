/**
 * @file 	column_collapse_ck.cpp
 * @brief 	2D column collapse.
 * @details Column collapse using computing kernels.
 * @author Shuang Li, Xiangyu Hu and Shuaihao Zhang
 */
#include "sphinxsys.h"
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.5;                       /**< Tank length. */
Real DH = 0.15;                      /**< Tank height. */
Real LL = 0.2;                       /**< Soil column length. */
Real LH = 0.1;                       /**< Soil column height. */
Real particle_spacing_ref = LH / 50; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;  /**< Extending width for boundary conditions. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
Real rho0_s = 2040;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 5.84e6;                                             // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // sound speed
Real friction_angle = 21.9 * Pi / 180;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d soil_block_translation = soil_block_halfsize;
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
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
std::vector<Vecd> soil_shape{
    Vecd(0, 0), Vecd(0, LH), Vecd(LL, LH), Vecd(LL, 0), Vecd(0, 0)};

class Soil : public MultiPolygonShape
{
  public:
    explicit Soil(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox initial_soil_block(Transform(soil_block_translation), soil_block_halfsize, "GranularBody");
    RealBody soil_block(sph_system, initial_soil_block);
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle);
    soil_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> soil_cell_linked_list(soil_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);

    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------

    Inner<> soil_block_inner(soil_block);
    Contact<> soil_block_contact(soil_block, {&wall_boundary});

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> soil_block_update_complex_relation(soil_block_inner, soil_block_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(soil_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(soil_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> soil_advection_step_setup(soil_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> soil_update_particle_position(soil_block);

    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep1stHalfWithWallRiemannCK>
        soil_acoustic_step_1st_half(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep2ndHalfWithWallRiemannCK>
        soil_acoustic_step_2nd_half(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        soil_density_regularization(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::StressDiffusionInnerCK> stress_diffusion(soil_block_inner);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> soil_acoustic_time_step(soil_block, 0.4);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(soil_block, "Density");
    StateDynamics<MainExecutionPolicy, continuum_dynamics::VerticalStressCK> vertical_stress(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "VerticalStress");
    StateDynamics<MainExecutionPolicy, continuum_dynamics::AccDeviatoricPlasticStrainCK> accumulated_deviatoric_plastic_strain(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "AccDeviatoricPlasticStrain");
    RestartIOCK<MainExecutionPolicy> restart_io(sph_system);

    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>>
        write_mechanical_energy(soil_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    soil_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    soil_block_update_complex_relation.exec();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real End_Time = 0.8;         /**< End time. */
    Real D_Time = End_Time / 40; /**< Time stamps for output of body states. */
    Real Dt = 0.1 * D_Time;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_acoustic_steps;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            /** outer loop for dual-time criteria time-stepping. */
            soil_density_regularization.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                soil_advection_step_setup.exec();
                Real dt = soil_acoustic_time_step.exec();
                stress_diffusion.exec();
                soil_acoustic_step_1st_half.exec(dt);
                soil_acoustic_step_2nd_half.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);

                interval_acoustic_steps += TickCount::now() - time_instance;

                /** screen output, write body observables and restart files  */
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4) << "	Time = "
                              << sv_physical_time->getValue()
                              << std::scientific << "	dt = " << dt << "\n";

                    if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                    {
                        write_mechanical_energy.writeToFile(number_of_iterations);
                    }
                    if (number_of_iterations % restart_output_interval == 0)
                        restart_io.writeToFile(number_of_iterations);
                }
                soil_update_particle_position.exec();
                number_of_iterations++;
                /** Update cell linked list and configuration. */
                time_instance = TickCount::now();
                soil_cell_linked_list.exec();
                soil_block_update_complex_relation.exec();
                interval_updating_configuration += TickCount::now() - time_instance;
            }
        }
        vertical_stress.exec();
        accumulated_deviatoric_plastic_strain.exec();
        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
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
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_mechanical_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_mechanical_energy.testResult();
    }

    return 0;
};
