/**
 * @file twisting_column.cpp
 * @brief This is an example of solid with classic Neohookean model
 * to demonstrate the robustness of the formulation with Kirchhoff stress decomposition.
 * @author Chi Zhang and Xiangyu Hu
 * @ref DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 6.0; /**< X-direction size. */
Real PH = 1.0; /**< Y-direction size. */
Real PW = 1.0; /**< Z-direction size. */
Real particle_spacing_ref = PH / 10.0;
/** You can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real BW = particle_spacing_ref * 0.0; /**< no wall boundary in this case. */
Real SL = particle_spacing_ref * 1.0; /**< Length of the holder is one layer particle. */
Vecd halfsize_column(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_column(0.5 * (PL - SL), 0.0, 0.0);
Vecd halfsize_holder(0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
Vecd translation_holder(-0.5 * (SL + BW), 0.0, 0.0);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 1100.0; /**< Reference density. */
Real poisson = 0.45;  /**< Poisson ratio. */
Real Youngs_modulus = 1.7e7;
Real angular_0 = -400.0;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
class Column : public ComplexShape
{
  public:
    explicit Column(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_column), halfsize_column);
        add<GeometricShapeBox>(Transform(translation_holder), halfsize_holder);
    }
};
//----------------------------------------------------------------------
//	Case dependent initial condition.
//----------------------------------------------------------------------
class InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        Real x = pos_[index_i][0];
        Real y = pos_[index_i][1];
        Real z = pos_[index_i][2];
        Real angular_velocity = angular_0 * sin((M_PI * x) / (2.0 * PL));
        Real local_radius = sqrt(pow(y, 2) + pow(z, 2));
        Real angular = atan2(y, z);

        if (x > 0.0)
        {
            vel_[index_i][1] = angular_velocity * local_radius * cos(angular);
            vel_[index_i][2] = -angular_velocity * local_radius * sin(angular);
        }
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    // Build up an SPHSystem and IO environment.
    // Please the make sure the global domain bounds are correctly defined.
    //----------------------------------------------------------------------
    Vec3d domain_lower_bound(-SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
    Vec3d domain_upper_bound(PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    // Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    SolidBody column(sph_system, makeShared<Column>("Column"));
    column.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    column.generateParticles<BaseParticles, Lattice>();

    ObserverBody my_observer(sph_system, "MyObserver");
    StdVec<Vecd> observation_location = {Vecd(PL, 0.0, 0.0)};
    my_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(column_inner);
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half(column_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(column_inner);

    /** Time step size calculation. We use CFL = 0.5 due to the very large twisting speed. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(column, 0.5);
    SimpleDynamics<InitialCondition> initial_condition(column);

    GeometricShapeBox holder_shape(Transform(translation_holder), halfsize_holder, "Holder");
    BodyRegionByParticle holder(column, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_velocity("Velocity", my_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", my_observer_contact);
    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    corrected_configuration.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    Real end_time = 0.5;
    Real output_period = end_time / 250.0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_displacement.writeToFile(number_of_iterations);
    write_velocity.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (number_of_iterations % 100 == 0)
            {
                std::cout << "N=" << number_of_iterations << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            number_of_iterations++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;
            write_displacement.writeToFile(number_of_iterations);
            write_velocity.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_displacement.generateDataBase(0.005);
        write_velocity.generateDataBase(0.005);
    }
    else
    {
        write_displacement.testResult();
        write_velocity.testResult();
    }

    return 0;
}
