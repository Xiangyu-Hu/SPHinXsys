#ifndef CONTINUUM_DYNAMICS_INNER_H
#define CONTINUUM_DYNAMICS_INNER_H
#include "all_fluid_dynamics.h"
#include "constraint_dynamics.h"
#include "continuum_particles.h"
#include "general_continuum.h"

namespace SPH
{
namespace continuum_dynamics
{
typedef DataDelegateSimple<ContinuumParticles> ContinuumDataSimple;
typedef DataDelegateInner<ContinuumParticles> ContinuumDataInner;

typedef DataDelegateSimple<PlasticContinuumParticles> PlasticContinuumDataSimple;
typedef DataDelegateInner<PlasticContinuumParticles> PlasticContinuumDataInner;

class ContinuumInitialCondition : public LocalDynamics, public PlasticContinuumDataSimple
{
  public:
    explicit ContinuumInitialCondition(SPHBody &sph_body);
    virtual ~ContinuumInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Mat3d> &stress_tensor_3D_;
};

class ContinuumAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
{
  public:
    explicit ContinuumAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.5);
    virtual ~ContinuumAcousticTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
};

/**
 * @class ArtificialStressAcceleration
 * @brief Implemented according to the literature:
 * Gray, J.P., Monaghan, J.J. and Swift, R., 2001. SPH elastic dynamics.
 * Computer methods in applied mechanics and engineering, 190(49-50), pp.6641-6662.
 */
class ArtificialStressAcceleration : public LocalDynamics, public ContinuumDataInner
{
  public:
    explicit ArtificialStressAcceleration(BaseInnerRelation &inner_relation, Real epsilon, Real exponent);
    virtual ~ArtificialStressAcceleration(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_, reference_spacing_;
    Real epsilon_, exponent_;
    StdLargeVec<Matd> &shear_stress_, artificial_stress_;
    StdLargeVec<Real> &p_, &rho_;
    StdLargeVec<Vecd> &acc_prior_;

    Matd getArtificialStress(const Matd &stress_tensor_i, const Real &rho_i);
};

/**
 * @class BaseShearStressIntegration
 * @brief Evolution of shear stress
 */
template <class DataDelegationType>
class BaseShearStressIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseShearStressIntegration(BaseRelationType &base_relation);
    virtual ~BaseShearStressIntegration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Matd> &velocity_gradient_, &shear_stress_;
    StdLargeVec<Real> &p_, &von_mises_stress_;
};

/**
 * @class ShearStressIntegration
 */
class ShearStressIntegration : public BaseShearStressIntegration<ContinuumDataInner>
{
  public:
    explicit ShearStressIntegration(BaseInnerRelation &inner_relation);
    virtual ~ShearStressIntegration(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
};

/**
 * @class BaseShearStressAcceleration
 */
class BaseShearStressAcceleration : public LocalDynamics, public ContinuumDataInner
{
  public:
    explicit BaseShearStressAcceleration(BaseInnerRelation &inner_relation);
    virtual ~BaseShearStressAcceleration(){};

  protected:
    StdLargeVec<Matd> &shear_stress_;
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &acc_prior_;
};

/**
 * @class ShearStressAcceleration
 */
class ShearStressAcceleration : public BaseShearStressAcceleration
{
  public:
    explicit ShearStressAcceleration(BaseInnerRelation &inner_relation)
        : BaseShearStressAcceleration(inner_relation){};
    virtual ~ShearStressAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShearStressAccelerationWithWall
 */
class ShearStressAccelerationWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>
{
  public:
    explicit ShearStressAccelerationWithWall(BaseContactRelation &wall_relation)
        : fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>(wall_relation){};
    virtual ~ShearStressAccelerationWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class PlasticShearStressIntegration
 */
class PlasticShearStressIntegration : public BaseShearStressIntegration<PlasticContinuumDataInner>
{
  public:
    explicit PlasticShearStressIntegration(BaseInnerRelation &inner_relation);
    virtual ~PlasticShearStressIntegration(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    StdLargeVec<Mat3d> &stress_tensor_3D_, &strain_tensor_3D_, &stress_rate_3D_, &strain_rate_3D_;
    StdLargeVec<Mat3d> &elastic_strain_tensor_3D_, &elastic_strain_rate_3D_;
    Real E_, nu_;
    StdLargeVec<Matd> &shear_stress_;
};

/**
 * @class ShearAccelerationIntegration
 * @brief This designed for hourglass free formulation
 */
class ShearAccelerationIntegration : public LocalDynamics, public ContinuumDataInner
{
  public:
    explicit ShearAccelerationIntegration(BaseInnerRelation &inner_relation);
    virtual ~ShearAccelerationIntegration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
    Real G_, smoothing_length_;
    StdLargeVec<Vecd> &vel_, &acc_prior_, acc_shear_;
    StdLargeVec<Real> &rho_;
};

/**
 * @class BaseMotionConstraint
 */
template <class DynamicsIdentifier>
class BaseMotionConstraint : public BaseLocalDynamics<DynamicsIdentifier>, public ContinuumDataSimple
{
  public:
    explicit BaseMotionConstraint(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier), ContinuumDataSimple(identifier.getSPHBody()),
          pos_(particles_->pos_), pos0_(particles_->pos0_),
          n_(particles_->n_), n0_(particles_->n0_),
          vel_(particles_->vel_), acc_(particles_->acc_){};

    virtual ~BaseMotionConstraint(){};

  protected:
    StdLargeVec<Vecd> &pos_, &pos0_;
    StdLargeVec<Vecd> &n_, &n0_;
    StdLargeVec<Vecd> &vel_, &acc_;
};
/**@class FixConstraint
 * @brief Constraint with zero velocity.
 */
template <class DynamicsIdentifier>
class FixConstraint : public BaseMotionConstraint<DynamicsIdentifier>
{
  public:
    explicit FixConstraint(DynamicsIdentifier &identifier)
        : BaseMotionConstraint<DynamicsIdentifier>(identifier){};
    virtual ~FixConstraint(){};

    void update(size_t index_i, Real dt = 0.0) { this->vel_[index_i] = Vecd::Zero(); };
};
using FixBodyConstraint = FixConstraint<SPHBody>;
using FixBodyPartConstraint = FixConstraint<BodyPartByParticle>;

/**
 * @class FixedInAxisDirection
 * @brief Constrain the velocity of a solid body part.
 */
class FixedInAxisDirection : public BaseMotionConstraint<BodyPartByParticle>
{
  public:
    FixedInAxisDirection(BodyPartByParticle &body_part, Vecd constrained_axises = Vecd::Zero());
    virtual ~FixedInAxisDirection(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd constrain_matrix_;
};

/**
 * @class ConstrainSolidBodyMassCenter
 * @brief Constrain the mass center of a solid body.
 */
class ConstrainSolidBodyMassCenter : public LocalDynamics, public ContinuumDataSimple
{
  private:
    Real total_mass_;
    Matd correction_matrix_;
    Vecd velocity_correction_;
    StdLargeVec<Vecd> &vel_;
    ReduceDynamics<QuantityMoment<Vecd>> compute_total_momentum_;

  protected:
    virtual void setupDynamics(Real dt = 0.0) override;

  public:
    explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd::Ones());
    virtual ~ConstrainSolidBodyMassCenter(){};

    void update(size_t index_i, Real dt = 0.0);
};

class BaseRelaxationPlastic : public LocalDynamics, public PlasticContinuumDataInner
{
  public:
    explicit BaseRelaxationPlastic(BaseInnerRelation &inner_relation);
    virtual ~BaseRelaxationPlastic(){};

  protected:
    PlasticContinuum &plastic_continuum_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
    StdLargeVec<Mat3d> &stress_tensor_3D_, &strain_tensor_3D_, &stress_rate_3D_, &strain_rate_3D_;
    StdLargeVec<Mat3d> &elastic_strain_tensor_3D_, &elastic_strain_rate_3D_;
};

/**
 * @class StressDiffusion
 */
class StressDiffusion : public BaseRelaxationPlastic
{
  public:
    explicit StressDiffusion(BaseInnerRelation &inner_relation, SharedPtr<Gravity> gravity_ptr, int axis);
    virtual ~StressDiffusion(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real zeta_ = 0.1;
    int axis_;
    Real rho0_, gravity_, smoothing_length_;
    Real phi_, diffusion_coeff_;
};

class ShearStressAccelerationWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>
{
  public:
    explicit ShearStressAccelerationWithWall(BaseContactRelation &wall_relation)
        : fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>(wall_relation){};
    virtual ~ShearStressAccelerationWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_DYNAMICS_INNER_H
