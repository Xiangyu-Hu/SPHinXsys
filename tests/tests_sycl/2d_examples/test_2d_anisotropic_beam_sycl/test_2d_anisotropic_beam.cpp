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
Real PL = 0.2;                                       // beam length
Real PH = 0.02;                                      // beam width
Real SL = 0.02;                                      // constrained length
int y_num = 10;                                      // particle number in y direction
Real ratio_ = 4.0;                                   // anisotropic ratio, also dp_x / dp_y
Real global_resolution = PH / y_num;                    // particle spacing in y direction
Real global_resolution_large = ratio_ * global_resolution; // large particle spacing, also the particle spacing in x direction
Real Total_PL = PL + SL;                             // total length
int x_num = Total_PL / global_resolution_large;         // particle number in x direction
//   anisotropic parameters
Vec2d scaling_vector = Vec2d(1.0, 1.0 / ratio_); // scaling_vector for defining the anisotropic kernel
Real scaling_factor = 1.0 / ratio_;              // scaling factor to calculate the time step
Real BW = global_resolution * 4;                    // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
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
namespace SPH
{
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
//	particle generation considering the anisotropic resolution
//----------------------------------------------------------------------
template <>
class ParticleGenerator<BaseParticles, Beam> : public ParticleGenerator<BaseParticles>
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
                Real x = -SL + (i + 0.5) * global_resolution_large;
                Real y = -PH / 2 + (j + 0.5) * global_resolution;
                addPositionAndVolumetricMeasure(Vecd(x, y), (global_resolution * global_resolution_large));
            }
        }
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
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body),
          elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_->getBaseMaterial())) {};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        Real x = pos_[index_i][0] / PL;
        if (x > 0.0)
        {
            vel_[index_i][1] = vf * elastic_solid_.ReferenceSoundSpeed() *
                               (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
        }
    };

  protected:
    ElasticSolid &elastic_solid_;
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

class AnisotropicCorrectConfiguration : public LocalDynamics, public DataDelegateInner
{
  public:
    AnisotropicCorrectConfiguration(BaseInnerRelation &inner_relation, int beta = 0, Real alpha = Real(0))
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          beta_(beta), alpha_(alpha), Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          B_(particles_->registerStateVariableData<Matd>("LinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          show_neighbor_(particles_->registerStateVariableData<Real>("ShowingNeighbor", Real(0.0))) {};
    virtual ~AnisotropicCorrectConfiguration() {};

  protected:
  protected:
    int beta_;
    Real alpha_;
    Real *Vol_;
    Matd *B_;
    Vecd *pos_;
    Real *show_neighbor_;

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Matd local_configuration = Eps * Matd::Identity();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
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
} // namespace SPH
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, global_resolution_large);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody beam_body(system, makeShared<Beam>("BeamBody"));
    beam_body.getSPHAdaptation().resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    beam_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, Beam>();

    ObserverBody beam_observer(system, "BeamObserver");
    beam_observer.getSPHAdaptation().resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    InteractionWithUpdate<AnisotropicCorrectConfiguration> beam_corrected_configuration(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(beam_body);
    SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
    BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
    SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_beam_states(beam_body);
    write_beam_states.addToWrite<Real>(beam_body, "ShowingNeighbor");
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_beam_tip_displacement("Position", beam_observer_contact);
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
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
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
    write_beam_tip_displacement.writeToFile(number_of_iterations);
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

                number_of_iterations++;
                dt = scaling_factor * computing_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                if (number_of_iterations % 100 == 0)
                {
                    std::cout << "N=" << number_of_iterations << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }
            }
        }

        TickCount t2 = TickCount::now();
        write_beam_tip_displacement.writeToFile(number_of_iterations);
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
