/**
 * @file muscle_activation.cpp
 * @brief This is the first example of electro activation of myocardium
 * @author Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real PL = 1.0;                   /**< Length of the myocardium body. */
Real PH = 1.0;                   /**< Thickness of the myocardium body. */
Real PW = 1.0;                   /**< Width of the myocardium body. */
Real resolution_ref = PH / 25.0; /**< Initial particle spacing. */
Real SL = 4.0 * resolution_ref;  /**< Extension for holder. */
Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_s = 1.0;
std::array<Real, 4> a0 = {0.059, 0.0, 0.0, 0.0};
std::array<Real, 4> b0 = {8.023, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real reference_voltage = 30.0;
Real linear_active_stress_factor = -0.5;
/** reference stress or bulk modulus to achieve weakly compressible condition */
Real bulk_modulus = 30.0 * reference_voltage * fabs(linear_active_stress_factor);
//----------------------------------------------------------------------
//	Case dependent muscle activation history.
//----------------------------------------------------------------------
class MyocardiumActivation
    : public active_muscle_dynamics::MuscleActivation
{
  public:
    explicit MyocardiumActivation(SPHBody &sph_body)
        : active_muscle_dynamics::MuscleActivation(sph_body),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {};

    void update(size_t index_i, Real dt)
    {
        Real voltage = pos0_[index_i][0] <= 0 ? 0.0 : reference_voltage * pos0_[index_i][0] / PL;
        active_contraction_stress_[index_i] += *physical_time_ <= 1.0
                                                   ? linear_active_stress_factor * voltage * dt
                                                   : 0.0;
    };

  protected:
    Real *physical_time_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-SL, -SL, -SL), Vecd(PL + SL, PH + SL, PW + SL));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox muscle_body_shape(Transform(translation_myocardium), halfsize_myocardium, "MyocardiumMuscleBody");
    SolidBody muscle_body(sph_system, muscle_body_shape);
    muscle_body.defineMaterial<ActiveMuscle<Muscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    muscle_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation muscle_body_inner(muscle_body);
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(muscle_body_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(muscle_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(muscle_body_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(muscle_body);
    SimpleDynamics<MyocardiumActivation> myocardium_activation(muscle_body);
    BodyRegionByParticle holder(muscle_body, makeShared<GeometricShapeBox>(Transform(translation_holder), halfsize_holder));
    SimpleDynamics<FixedInAxisDirection> constrain_holder(holder, Vecd(0.0, 1.0, 1.0));
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 1.2;
    Real output_period = end_time / 60.0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }
            myocardium_activation.exec(dt);
            stress_relaxation_first_half.exec(dt);
            constrain_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;
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

    return 0;
}
