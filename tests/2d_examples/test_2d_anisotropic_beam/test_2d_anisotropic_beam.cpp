/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D anisotropic beam example-one body version           *
 * ----------------------------------------------------------------------------*
 * This is a test cases using ASPH method with anisotropic kerne for            *
 *  simulating solid. Particle space is anisotropic in different directions in this beam.  *
 * @author	Xiaojing Tang
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;  // beam length
Real PH = 0.02; // beam width
Real SL = 0.02; // constrained length

//   particle spacings and particle numbers
int y_num = 10;                                      // particle number in y direction
Real ratio_ = 4.0;                                   // anisotropic ratio, also dp_x / dp_y
Real resolution_ref = PH / y_num;                    // particle spacing in y direction
Real resolution_ref_large = ratio_ * resolution_ref; // large particle spacing, also the particle spacing in x direction
Real Total_PL = PL + SL;                             // total length
int x_num = Total_PL / resolution_ref_large;         // particle number in x direction
//   anisotropic parameters
Vec2d scaling_vector = Vec2d(1.0, 1.0 / ratio_); // scaling_vector for defining the anisotropic kernel
Real scaling_factor = 1.0 / ratio_;              // scaling factor to calculate the time step

//----------------------------------------------------------------------
//	particle generation considering the anisotropic resolution
//----------------------------------------------------------------------
class AnisotropicParticleGenerator : public ParticleGenerator
{
  public:
 AnisotropicParticleGenerator(SPHBody &sph_body) : ParticleGenerator(sph_body){};

    virtual void initializeGeometricVariables() override
    {
        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < y_num; j++)
            {
                Real x = -SL + (i + 0.5) * resolution_ref_large;
                Real y = -PH / 2 + (j + 0.5) * resolution_ref;
                initializePositionAndVolumetricMeasure(Vecd(x, y), (resolution_ref * resolution_ref_large));
            }
        }
    }
};

Real BW = resolution_ref * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
                                 Vec2d(PL + 3.0 * BW, PL / 2.0));
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
std::vector<Vecd> beam_base_shape{
    Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(0.0, PH / 2),
    Vecd(0.0, -PH / 2), Vecd(-SL, -PH / 2)};
// a beam shape
std::vector<Vecd> beam_shape{
    Vecd(0.0, -PH / 2), Vecd(0.0, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(0.0, -PH / 2)};
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
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
            vel_[index_i][1] = vf * particles_->elastic_solid_.ReferenceSoundSpeed() *
                               (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
        }
    };
};
//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::sub);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	calculate correction matrix B to keep the accuracy
//----------------------------------------------------------------------

class AnisotropicCorrectConfiguration : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    AnisotropicCorrectConfiguration(BaseInnerRelation &inner_relation, int beta = 0, Real alpha = Real(0))
        : LocalDynamics(inner_relation.getSPHBody()),
          GeneralDataDelegateInner(inner_relation),
          beta_(beta), alpha_(alpha),
          B_(*particles_->getVariableByName<Matd>("CorrectionMatrix")),
          pos_(particles_->pos_)
    {
        particles_->registerVariable(show_neighbor_, "ShowingNeighbor", Real(0.0));
    }
    virtual ~AnisotropicCorrectConfiguration(){};

  protected:
  protected:
    int beta_;
    Real alpha_;
    StdLargeVec<Matd> &B_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> show_neighbor_;

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Matd local_configuration = Eps * Matd::Identity();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
            Vecd e_ij = inner_neighborhood.e_ij_[n];
            if (index_i == 67)
            {
                show_neighbor_[index_j] = 1.0;
            };

            Vecd gradW_ijV_j = dW_ijV_j * e_ij;
            Vecd r_ji = pos_[index_i] - pos_[index_j];
            local_configuration -= r_ji * gradW_ijV_j.transpose();
        }
        B_[index_i] = local_configuration;
    };

    void update(size_t index_i, Real dt)
    {
        Real det_sqr = pow(B_[index_i].determinant(), beta_);
        Matd inverse = B_[index_i].inverse();
        B_[index_i] = (det_sqr * inverse + alpha_ * Matd::Identity()) / (alpha_ + det_sqr);
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
    SPHSystem system(system_domain_bounds, resolution_ref_large);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
       //	Creating body, materials and particles.
       //----------------------------------------------------------------------
    SolidBody beam_body(system, makeShared<Beam>("BeamBody"));
    beam_body.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    beam_body.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<AnisotropicParticleGenerator>();

    ObserverBody beam_observer(system, "BeamObserver");
    beam_observer.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    beam_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
    // corrected strong configuration 
    InteractionWithUpdate<AnisotropicCorrectConfiguration> beam_corrected_configuration(beam_body_inner);
    beam_body.addBodyStateForRecording<Real>("ShowingNeighbor");
    // time step size calculation
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
    // stress relaxation for the beam
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);
    // clamping a solid body part.  
    BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_beam_base(beam_base);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    beam_initial_velocity.exec();
    beam_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup computing time-step controls.
    //----------------------------------------------------------------------
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
    write_beam_states.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);

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
                stress_relaxation_first_half.exec(dt);
                constraint_beam_base.exec();
                stress_relaxation_second_half.exec(dt);

                ite++;
                dt = scaling_factor * computing_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }
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

    if (system.GenerateRegressionData())
    {
        write_beam_tip_displacement.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
    }
    else
    {
        write_beam_tip_displacement.testResult();
    }

    return 0;
}
