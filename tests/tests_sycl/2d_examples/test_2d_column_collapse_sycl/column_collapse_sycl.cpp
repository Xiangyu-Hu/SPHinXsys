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
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d soil_block_translation = soil_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
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
    auto &initial_soil_block = sph_system.addShape<GeometricShapeBox>(
        Transform(soil_block_translation), soil_block_halfsize, "GranularBody");
    auto &soil_block = sph_system.addBody<RealBody>(initial_soil_block);
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle);
    soil_block.generateParticles<BaseParticles, Lattice>();

    auto &wall_shape = sph_system.addShape<ComplexShape>("WallBoundary");
    wall_shape.add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
    wall_shape.subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    auto &wall_boundary = sph_system.addBody<SolidBody>(wall_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &soil_block_inner = sph_system.addInnerRelation(soil_block);
    auto &soil_block_contact = sph_system.addContactRelation(soil_block, wall_boundary);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary).exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);

    ParticleDynamicsGroup soil_update_configuration;
    soil_update_configuration.add(&main_methods.addCellLinkedListDynamics(soil_block));
    soil_update_configuration.add(&main_methods.addRelationDynamics(soil_block_inner, soil_block_contact));
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(soil_block, gravity);
    auto &soil_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(soil_block);
    auto &soil_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(soil_block);

    auto &soil_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
                        continuum_dynamics::PlasticAcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(soil_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(soil_block_contact);
    auto &soil_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
                        continuum_dynamics::PlasticAcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(soil_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(soil_block_contact);
    auto &soil_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(soil_block_inner)
            .addPostContactInteraction(soil_block_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(soil_block);
    auto &stress_diffusion = main_methods.addInteractionDynamics<continuum_dynamics::StressDiffusionCK>(soil_block_inner);

    auto &soil_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(soil_block, 0.4);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_state_recorder.addToWrite<Real>(soil_block, "Density");
    body_state_recorder.addDerivedVariableToWrite<continuum_dynamics::VerticalStressCK>(soil_block);
    body_state_recorder.addDerivedVariableToWrite<continuum_dynamics::AccDeviatoricPlasticStrainCK>(soil_block);
    auto &record_water_mechanical_energy = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, TotalMechanicalEnergyCK>(soil_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    constant_gravity.exec();
    wall_cell_linked_list.exec();
    soil_update_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    int observation_sample_interval = screen_output_interval * 2;
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
    body_state_recorder.writeToFile();
    record_water_mechanical_energy.writeToFile(number_of_iterations);
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
                        record_water_mechanical_energy.writeToFile(number_of_iterations);
                    }
                }
                soil_update_particle_position.exec();
                number_of_iterations++;
                /** Update cell linked list and configuration. */
                time_instance = TickCount::now();
                soil_update_configuration.exec();
                interval_updating_configuration += TickCount::now() - time_instance;
            }
        }
        body_state_recorder.writeToFile();
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
        record_water_mechanical_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        record_water_mechanical_energy.testResult();
    }

    return 0;
};
