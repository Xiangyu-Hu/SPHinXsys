/**
 * @file test_2d_anisotropic_beam.cpp
 * @brief This is a test cases using anisotropic kernel for simulating solid.
 * Particle space is anisotropic in different directions of the beam.
 * @author Xiaojing Tang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real total_physical_time = 1.0;                                  /**< TOTAL SIMULATION TIME*/
Real PL = 0.2;                                                   // beam length
Real PH = 0.02;                                                  // beam width
Real SL = 0.02;                                                  // constrained length
int y_num = 10;                                                  // particle number in y direction
Real anisotropic_ratio = 4.0;                                    // anisotropic ratio, also dp_x / dp_y
Real particle_spacing_ref = PH / y_num;                          // particle spacing in y direction
Real maximum_spacing = anisotropic_ratio * particle_spacing_ref; // large particle spacing, also the particle spacing in x direction
Real Total_PL = PL + SL;                                         // total length
int x_num = Total_PL / maximum_spacing;                          // particle number in x direction
//   anisotropic parameters
Vec2d scaling_vector = Vec2d(anisotropic_ratio, 1.0); // scaling_vector for defining the anisotropic kernel
Vec2d orientation_vector = Vec2d::UnitX();            // orientation vector for defining the anisotropic kernel
Real scaling_factor = 1.0 / anisotropic_ratio;        // scaling factor to calculate the time step
Real BW = particle_spacing_ref * 4;                   // boundary width, at least three particles
/** Domain bounds of the sph_system. */
BoundingBoxd system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;       // Poisson ratio
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometric shapes used in the sph_system.
//----------------------------------------------------------------------
GeometricShapeBox beam_shape(BoundingBoxd(Vecd(-SL, -PH / 2), Vecd(PL, PH / 2)), "BeamBody");
GeometricShapeBox beam_base_shape(BoundingBoxd(Vecd(-SL, -PH / 2), Vecd(0.0, PH / 2)), "BeamBase");
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	Adaptation used in the sph_system.
//----------------------------------------------------------------------
PrescribedAnisotropy y_refinement(scaling_vector, orientation_vector, particle_spacing_ref, 1.15, 1.0);
namespace SPH
{
//----------------------------------------------------------------------
//	particle generation considering the anisotropic resolution
//----------------------------------------------------------------------
template <>
class ParticleGenerator<BaseParticles, UserDefined> : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles) {};

    virtual void prepareGeometricData() override
    {
        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < y_num; j++)
            {
                Real x = -SL + (i + 0.5) * maximum_spacing;
                Real y = -PH / 2 + (j + 0.5) * particle_spacing_ref;
                addPositionAndVolumetricMeasure(Vecd(x, y), (particle_spacing_ref * maximum_spacing));
            }
        }
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class LinearProfile : public ReturnFunction<Vecd>
{
    Real vf_ = vf;
    Real kl_ = kl;
    Real M_ = M;
    Real N_ = N;
    Real Q_ = Q;
    Real PL_ = PL;
    Real c0_;

  public:
    LinearProfile(Real c0) : c0_(c0) {};

    Vecd operator()(const Vec2d &position)
    {
        Real x = position[0] / PL_;
        Vecd result = Vec2d::Zero();
        if (x > 0.0)
        {
            result[1] = vf_ * c0_ * (M_ * (cos(kl_ * x) - cosh(kl_ * x)) - N_ * (sin(kl_ * x) - sinh(kl_ * x))) / Q_;
        }
        return result;
    }
};
} // namespace SPH
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, maximum_spacing);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &beam_body = sph_system.addAdaptiveBody<RealBody, PrescribedAnisotropy>(y_refinement, beam_shape);
    auto *beam_material = beam_body.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, UserDefined>();
    auto &beam_base = beam_body.addBodyPart<BodyRegionByParticle>(beam_base_shape);

    auto &beam_observer = sph_system.addAdaptiveBody<ObserverBody, PrescribedAnisotropy>(y_refinement, "BeamObserver");
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    auto &beam_body_inner = sph_system.addInnerRelation(beam_body, ConfigType::Lagrangian);
    auto &my_observer_contact = sph_system.addContactRelation(beam_observer, beam_body, ConfigType::Lagrangian);
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
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<LinearProfile>>(
                    beam_body, "Velocity", beam_material->ReferenceSoundSpeed())
        .exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(beam_body));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(beam_body_inner));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(my_observer_contact));

    auto &linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(beam_body_inner);
    auto &constraint_holder = main_methods.addStateDynamics<ConstantConstraintCK, Vecd>(beam_base, "Velocity", Vec2d::Zero());

    auto &acoustic_step_1st_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, NoKernelCorrectionCK>(beam_body_inner);
    auto &acoustic_step_2nd_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(beam_body_inner);

    auto &acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(beam_body, 0.2);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &write_real_body_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    write_real_body_states.addToWrite<Real>(beam_body, "Density");
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    size_t acoustic_steps = 1;
    int screening_interval = 100;
    int observation_interval = screening_interval;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 100.0);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_acoustic_step;
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    lagrangian_configuration.exec();
    linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        constraint_holder.exec();
        acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        /** Output body state during the simulation according output_interval. */
        time_instance = TickCount::now();
        /** screen output, write body observables and restart files  */
        if (acoustic_steps == 1 || acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (state_recording())
        {
            write_real_body_states.writeToFile();
        }
        interval_output += TickCount::now() - time_instance;

        acoustic_steps++;
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";

    return 0;
}
