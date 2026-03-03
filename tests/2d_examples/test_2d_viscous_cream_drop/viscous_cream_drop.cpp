/**
 * @file 	viscous_cream_drop.cpp
 * @brief   Viscous cream adhering to the platform falls due to gravity.
 * @author 	Liezhao Wu, Xiaojing Tang and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real global_resolution = 0.005;    /* reference resolution. */
Real DL = 0.5;                  /* platform length. */
Real DH = 2.0;                  /* height. */
Real BW = global_resolution * 5.0; /* width for boundary. */
Real cream_radius = global_resolution * 20.0;
Vec2d cream_center(0.0, -cream_radius);
BoundingBoxd system_domain_bounds(Vec2d(-DL, -DH), Vec2d(DL, BW));
Real sqrt_3 = sqrt(3);
// observer location
StdVec<Vecd> observation_location = {cream_center};
//----------------------------------------------------------------------
//	Global parameters on material properties.
//----------------------------------------------------------------------
Real gravity_g = 9.8;
// shaving cream material
Real rho0_s = 77.7;                                                                                     /*  density  */
Real Bulk_modulus = 1.09e5;                                                                             /*  bulk modulus */
Real Shear_modulus = 290.0;                                                                             /*  shear modulus */
Real yield_stress = 31.9;                                                                               /*  yield stress */
Real viscosity = 27.2;                                                                                  /* viscosity  */
Real Herschel_Bulkley_power = 0.22;                                                                     /*   Herschel_Bulkley_power. */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /*   Young's modulus  */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /*  Poisson's ratio */
//----------------------------------------------------------------------
//	Geometric shapes.
//----------------------------------------------------------------------
std::vector<Vecd> createPlatformShape()
{
    std::vector<Vecd> platform_shape;
    platform_shape.push_back(Vecd(-0.5 * DL, 0.0));
    platform_shape.push_back(Vecd(-0.5 * DL, BW));
    platform_shape.push_back(Vecd(0.5 * DL, BW));
    platform_shape.push_back(Vecd(0.5 * DL, 0.0));
    platform_shape.push_back(Vecd(-0.5 * DL, 0.0));

    return platform_shape;
}
std::vector<Vecd> createCreamUpperShape()
{
    std::vector<Vecd> cream_upper_shape;
    cream_upper_shape.push_back(Vecd(-sqrt_3 * cream_radius, 0.0));
    cream_upper_shape.push_back(Vecd(sqrt_3 * cream_radius, 0.0));
    cream_upper_shape.push_back(Vecd(sqrt_3 * cream_radius / 2.0, -3.0 * cream_radius / 2.0));
    cream_upper_shape.push_back(Vecd(-sqrt_3 * cream_radius / 2.0, -3.0 * cream_radius / 2.0));
    cream_upper_shape.push_back(Vecd(-sqrt_3 * cream_radius, 0.0));

    return cream_upper_shape;
}
class Cream : public MultiPolygonShape
{
  public:
    explicit Cream(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createPlatformShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createCreamUpperShape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(cream_center, cream_radius, 100, ShapeBooleanOps::add);
    }
};
MultiPolygon createPlatformConstraint()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createPlatformShape(), ShapeBooleanOps::add);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.setGenerateRegressionData(false);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody cream(sph_system, makeShared<Cream>("Cream"));
    cream.defineBodyLevelSetShape();
    cream.defineMaterial<ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                              yield_stress, viscosity, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cream.generateParticles<BaseParticles, Reload>(cream.getName())
        : cream.generateParticles<BaseParticles, Lattice>();

    ObserverBody cream_observer(sph_system, "CreamObserver");
    cream_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation cream_inner(cream);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> cream_random_particles(cream);
        RelaxationStepInner cream_relaxation_step_inner(cream_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_cream_state(sph_system);
        ReloadParticleIO write_particle_reload_files(cream);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        cream_random_particles.exec(0.25);
        write_cream_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            cream_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_cream_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of cream particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation cream_inner(cream);
    ContactRelation cream_observer_contact(cream_observer, {&cream});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(cream, gravity);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> cream_corrected_configuration(cream_inner);

    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf> cream_stress_relaxation_first_half(cream_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> cream_stress_relaxation_second_half(cream_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> cream_get_time_step_size(cream, 0.2);
    BodyRegionByParticle platform(cream, makeShared<MultiPolygonShape>(createPlatformConstraint()));
    SimpleDynamics<FixBodyPartConstraint> platform_constraint(platform);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        cream_displacement_recording("Position", cream_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    cream_corrected_configuration.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    cream_displacement_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 0.75;
    Real end_time = T0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real output_interval = 0.01 * T0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            if (ite % screen_output_interval == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";

                if (ite != 0 && ite % observation_sample_interval == 0)
                {
                    cream_displacement_recording.writeToFile(ite);
                }
            }
            cream_stress_relaxation_first_half.exec(dt);
            platform_constraint.exec(dt);
            cream_stress_relaxation_second_half.exec(dt);

            ite++;
            dt = cream_get_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        cream_displacement_recording.generateDataBase(0.05);
    }
    else
    {
        cream_displacement_recording.testResult();
    }

    return 0;
}
