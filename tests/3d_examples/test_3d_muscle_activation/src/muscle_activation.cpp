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
Real a0[4] = {0.059, 0.0, 0.0, 0.0};
Real b0[4] = {8.023, 0.0, 0.0, 0.0};
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
        : active_muscle_dynamics::MuscleActivation(sph_body){};

    void update(size_t index_i, Real dt)
    {
        Real voltage = pos0_[index_i][0] <= 0 ? 0.0 : reference_voltage * pos0_[index_i][0] / PL;
        active_contraction_stress_[index_i] += GlobalStaticVariables::physical_time_ <= 1.0
                                                   ? linear_active_stress_factor * voltage * dt
                                                   : 0.0;
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-SL, -SL, -SL), Vecd(PL + SL, PH + SL, PW + SL));
    SPHSystem system(system_domain_bounds, resolution_ref);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    SolidBody myocardium_muscle_body(system, makeShared<TransformShape<GeometricShapeBox>>(
                                                 Transform(translation_myocardium), halfsize_myocardium, "MyocardiumMuscleBody"));
    myocardium_muscle_body.defineParticlesAndMaterial<
        ElasticSolidParticles, ActiveMuscle<Muscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    myocardium_muscle_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation myocardium_muscle_body_inner(myocardium_muscle_body);
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(myocardium_muscle_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(myocardium_muscle_body_inner);
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(myocardium_muscle_body_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(myocardium_muscle_body);
    SimpleDynamics<MyocardiumActivation> myocardium_activation(myocardium_muscle_body);
    BodyRegionByParticle holder(myocardium_muscle_body, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_holder), halfsize_holder));
    SimpleDynamics<solid_dynamics::FixedInAxisDirection> constrain_holder(holder, Vecd(0.0, 1.0, 1.0));
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            myocardium_activation.exec(dt);
            stress_relaxation_first_half.exec(dt);
            constrain_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
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
