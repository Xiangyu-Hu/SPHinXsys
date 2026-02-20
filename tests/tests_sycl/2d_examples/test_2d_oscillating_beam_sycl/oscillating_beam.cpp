/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D oscillation beam example-one body version           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases, also the first case for            *
 * understanding SPH method for solid simulation.                              *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;  // beam length
Real PH = 0.02; // for thick plate; =0.01 for thin plate
Real SL = 0.06; // depth of the insert
// reference particle spacing
Real global_resolution = PH / 10.0;
Real BW = global_resolution * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
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
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a beam base shape
GeometricShapeBox beam_base_shape(
    BoundingBoxd(Vec2d(-SL - BW, -PH / 2 - BW), Vec2d(0.0, PH / 2 + BW)), "BeamBase");
GeometricShapeBox beam_column(
    BoundingBoxd(Vec2d(-SL, -PH / 2), Vec2d(PL, PH / 2)), "BeamColumn");
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class LinearProfile : public ReturnFunction<Vecd>
{
    Real vf_ = vf;
    Real c0_;
    Real kl_ = kl;
    Real M_ = M;
    Real N_ = N;
    Real Q_ = Q;
    Real PL_ = PL;

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
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
       //	Creating body, materials and particles.
       //----------------------------------------------------------------------
    ComplexShape beam_shape("BeamBody");
    beam_shape.add(&beam_base_shape);
    beam_shape.add(&beam_column);
    SolidBody beam_body(sph_system, beam_shape);
    auto *beam_material = beam_body.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, Lattice>();
    ComplexShape beam_constrain_shape("BeamConstrain");
    beam_constrain_shape.add(&beam_base_shape);
    beam_constrain_shape.subtract(&beam_column);
    BodyRegionByParticle beam_base(beam_body, beam_constrain_shape);

    ObserverBody beam_observer(sph_system, "BeamObserver");
    beam_observer.defineAdaptationRatios(1.15, 2.0);
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    auto &beam_body_inner = sph_system.addInnerRelation(beam_body, ConfigType::Lagrangian);
    auto &beam_observer_contact = sph_system.addContactRelation(beam_observer, beam_body, ConfigType::Lagrangian);
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
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(beam_body));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(beam_body_inner));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(beam_observer_contact));

    auto &beam_corrected_configuration = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(beam_body_inner);
    auto &stress_relaxation_first_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, LinearCorrectionCK>(beam_body_inner);
    auto &stress_relaxation_second_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(beam_body_inner);

    auto &computing_time_step_size = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(beam_body);
    auto &constraint_beam_base = main_methods.addStateDynamics<FixBodyPartConstraintCK>(beam_base);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_beam_states(sph_system);
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    lagrangian_configuration.exec();
    beam_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup computing time-step controls.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 1.0;
    Real end_time = T0;
    // time step size for output file
    Real output_interval = 0.01 * T0;
    Real Dt = 0.1 * output_interval; /**< Time period for data observing */
    Real dt = 0.0;                   // default acoustic time step sizes

    // statistics for computing time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    write_beam_states.writeToFile();

    // computation loop starts
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < output_interval)
        {

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                stress_relaxation_first_half.exec(dt);
                constraint_beam_base.exec();
                stress_relaxation_second_half.exec(dt);

                ite++;
                dt = computing_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }
            }
        }

        TickCount t2 = TickCount::now();
        write_beam_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
