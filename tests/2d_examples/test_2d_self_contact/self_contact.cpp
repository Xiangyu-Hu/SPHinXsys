/**
 * @file 	self_contact.cpp
 * @brief 	This is the case file for the test of dynamic self contact.
 * @author   Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;  // beam length
Real PH = 0.01; // for thick plate; 0.01 for thin plate
Real SL = 0.04; // depth of the insert
Real resolution_ref = PH / 10.0;
Real BW = resolution_ref * 4; // boundary width, at least three particles
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 1.0e5; // reference Youngs modulus
Real poisson = 0.45;         // Poisson ratio
//----------------------------------------------------------------------
//	Global parameters for initial condition
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.15;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometric elements used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> createBeamBaseShape()
{
    // geometry
    std::vector<Vecd> beam_base_shape;
    beam_base_shape.push_back(Vecd(-SL - BW, -PH / 2 - BW));
    beam_base_shape.push_back(Vecd(-SL - BW, PH / 2 + BW));
    beam_base_shape.push_back(Vecd(0.0, PH / 2 + BW));
    beam_base_shape.push_back(Vecd(0.0, -PH / 2 - BW));
    beam_base_shape.push_back(Vecd(-SL - BW, -PH / 2 - BW));

    return beam_base_shape;
}
std::vector<Vecd> createBeamShape()
{
    std::vector<Vecd> beam_shape;
    beam_shape.push_back(Vecd(-SL, -PH / 2));
    beam_shape.push_back(Vecd(-SL, PH / 2));
    beam_shape.push_back(Vecd(PL, PH / 2));
    beam_shape.push_back(Vecd(PL, -PH / 2));
    beam_shape.push_back(Vecd(-SL, -PH / 2));

    return beam_shape;
}
//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createBeamBaseShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Body part shape usually for impose constraints.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createBeamBaseShape(), ShapeBooleanOps::add);
    multi_polygon.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class BeamInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit BeamInitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        Real x = pos_[index_i][0] / PL;
        if (x > 0.0)
        {
            vel_[index_i][1] = vf * particles_->elastic_solid_.ReferenceSoundSpeed() / Q *
                               (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x)));
        }
    };
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0));
    SPHSystem system(system_domain_bounds, resolution_ref);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody beam_body(system, makeShared<Beam>("BeamBody"));
    beam_body.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<ParticleGeneratorLattice>();

    ObserverBody beam_observer(system, "BeamObserver");
    beam_observer.defineAdaptationRatios(1.15, 2.0);
    StdVec<Vecd> beam_observation_location = {Vecd(PL, 0.0)};
    beam_observer.generateParticles<ObserverParticleGenerator>(beam_observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    SelfSurfaceContactRelation beam_self_contact(beam_body);
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
    SimpleDynamics<TimeStepInitialization> reset_prior_acceleration(beam_body);
    InteractionWithUpdate<CorrectedConfigurationInner> beam_corrected_configuration(beam_body_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> beam_self_contact_density(beam_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> beam_self_contact_forces(beam_self_contact);
    BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_beam_base(beam_base);
    //-----------------------------------------------------------------------------
    //	outputs
    //-----------------------------------------------------------------------------
    IOEnvironment io_environment(system);
    beam_body.addBodyStateForRecording<Real>("SelfContactDensity");
    BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
    //-----------------------------------------------------------------------------
    //	Setup particle configuration and initial conditions
    //-----------------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    beam_initial_velocity.exec();
    beam_corrected_configuration.exec();
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    // starting time zero
    GlobalStaticVariables::physical_time_ = 0.0;
    write_beam_states.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);

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

    // computation loop starts
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < output_interval)
        {

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {

                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                reset_prior_acceleration.exec();
                beam_self_contact_density.exec();
                beam_self_contact_forces.exec();
                beam_body.updateCellLinkedList();
                beam_self_contact.updateConfiguration();

                stress_relaxation_first_half.exec(dt);
                constraint_beam_base.exec();
                stress_relaxation_second_half.exec(dt);

                ite++;
                dt = computing_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
        }

        write_beam_tip_displacement.writeToFile(ite);

        TickCount t2 = TickCount::now();
        write_beam_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_beam_tip_displacement.generateDataBase(1.0e-2);
    }
    else
    {
        write_beam_tip_displacement.testResult();
    }

    return 0;
}
