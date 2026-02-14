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
Real PL = 0.2;                                             // beam length
Real PH = 0.02;                                            // beam width
Real SL = 0.02;                                            // constrained length
int y_num = 10;                                            // particle number in y direction
Real ratio_ = 4.0;                                         // anisotropic ratio, also dp_x / dp_y
Real global_resolution = PH / y_num;                       // particle spacing in y direction
Real global_resolution_large = ratio_ * global_resolution; // large particle spacing, also the particle spacing in x direction
Real Total_PL = PL + SL;                                   // total length
int x_num = Total_PL / global_resolution_large;            // particle number in x direction
//   anisotropic parameters
Vec2d scaling_vector = Vec2d(1.0, 1.0 / ratio_); // scaling_vector for defining the anisotropic kernel
Vec2d orientation_vector = Vec2d::UnitX();       // orientation vector for defining the anisotropic kernel
Real scaling_factor = 1.0 / ratio_;              // scaling factor to calculate the time step
Real BW = global_resolution * 4;                 // boundary width, at least three particles
/** Domain bounds of the sph_system. */
BoundingBoxd system_domain_bounds(
    Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0));
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
AnisotropicAdaptation y_refinement(scaling_vector, orientation_vector, global_resolution, 1.15, 1.0);
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
    SPHSystem sph_system(system_domain_bounds, global_resolution_large);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &beam_body = sph_system.addAdaptiveBody<RealBody, AnisotropicAdaptation>(y_refinement, beam_shape);
    beam_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, UserDefined>();
    //    auto &beam_base = beam_body.addBodyPart<BodyRegionByParticle>(beam_base_shape);

    auto &beam_observer = sph_system.addAdaptiveBody<ObserverBody, AnisotropicAdaptation>(y_refinement, "BeamObserver");
    beam_observer.generateParticles<ObserverParticles>(observation_location);
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
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(beam_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &write_real_body_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    write_real_body_states.addToWrite<Real>(beam_body, "Density");
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    update_cell_linked_list.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    return 0;
}
