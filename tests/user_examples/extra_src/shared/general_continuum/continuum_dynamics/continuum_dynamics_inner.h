#ifndef CONTINUUM_DYNAMICS_INNER_H
#define CONTINUUM_DYNAMICS_INNER_H
#include "constraint_dynamics.h"
#include "fluid_dynamics_inner.hpp"
#include "general_continuum.h"
#include "continuum_particles.h"

namespace SPH
{
namespace continuum_dynamics
{
typedef DataDelegateSimple<ContinuumParticles> ContinuumDataSimple;
typedef DataDelegateInner<ContinuumParticles> ContinuumDataInner;

class ContinuumInitialCondition : public LocalDynamics, public ContinuumDataSimple
{
  public:
    explicit ContinuumInitialCondition(SPHBody &sph_body);
    virtual ~ContinuumInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
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
 * @class BaseIntegration
 * @brief Pure abstract base class for all fluid relaxation schemes
 */
class BaseRelaxation : public LocalDynamics, public ContinuumDataInner
{
  public:
    explicit BaseRelaxation(BaseInnerRelation &inner_relation);
    virtual ~BaseRelaxation(){};

  protected:
    GeneralContinuum &continuum_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
};

/**
 * @class BaseArtificialStressRelaxation
 */
class BaseArtificialStressRelaxation : public BaseRelaxation
{
  public:
    explicit BaseArtificialStressRelaxation(BaseInnerRelation &inner_relation, Real epsilon = 0.3);
    virtual ~BaseArtificialStressRelaxation(){};
    Matd repulsiveForce(Matd stress_tensor_i, Real rho_i);

  protected:
    Real smoothing_length_, reference_spacing_, epsilon_;
};

/**
 * @class ArtificialNormalStressRelaxation
 */
class ArtificialNormalShearStressRelaxation : public BaseArtificialStressRelaxation
{
  public:
    explicit ArtificialNormalShearStressRelaxation(BaseInnerRelation &inner_relation, Real exponent = 4);
    virtual ~ArtificialNormalShearStressRelaxation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &shear_stress_;
    StdLargeVec<Vecd> &acc_shear_;
    Real exponent_;
};

/**
 * @class ShearStressRelaxation1stHalf
 */
class ShearStressRelaxation1stHalf : public BaseRelaxation
{
  public:
    explicit ShearStressRelaxation1stHalf(BaseInnerRelation &inner_relation);
    virtual ~ShearStressRelaxation1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &shear_stress_, &shear_stress_rate_;
    StdLargeVec<Vecd> &acc_shear_;
};
/**
 * @class ShearStressRelaxation2ndHalf
 */
class ShearStressRelaxation2ndHalf : public BaseRelaxation
{
  public:
    explicit ShearStressRelaxation2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~ShearStressRelaxation2ndHalf(){};
    // void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &shear_stress_, &shear_stress_rate_, &velocity_gradient_, &strain_tensor_, &strain_tensor_rate_;
    StdLargeVec<Real> &von_mises_stress_;
};

//=================================================================================================//
//===================================Non-hourglass formulation=====================================//
//=================================================================================================//
/**
    * @class Integration1stHalf
    * @brief Pressure relaxation scheme with the mostly used Riemann solver.
*/
template <class FluidDynamicsType>
class BaseIntegration1stHalf : public FluidDynamicsType
{
public:
    explicit BaseIntegration1stHalf(BaseInnerRelation& inner_relation);
    virtual ~BaseIntegration1stHalf() {};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:
    StdLargeVec<Vecd>& acc_shear_;
};
using Integration1stHalf = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalf>;
using Integration1stHalfRiemann = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalfRiemann>;

/**
* @class ShearAccelerationRelaxation
*/
class ShearAccelerationRelaxation : public BaseRelaxation
{
public:
    explicit ShearAccelerationRelaxation(BaseInnerRelation& inner_relation);
    virtual ~ShearAccelerationRelaxation() {};
    void interaction(size_t index_i, Real dt = 0.0);
protected:
    Real G_, smoothing_length_;
    StdLargeVec<Matd>& shear_stress_, & B_;
    StdLargeVec<Vecd>& acc_shear_;
};

/**
 * @class AngularConservativeShearAccelerationRelaxation
 */
class AngularConservativeShearAccelerationRelaxation : public ShearAccelerationRelaxation
{
public:
    explicit AngularConservativeShearAccelerationRelaxation(BaseInnerRelation& inner_relation)
        : ShearAccelerationRelaxation(inner_relation) {};
    virtual ~AngularConservativeShearAccelerationRelaxation() {};

    void interaction(size_t index_i, Real dt = 0.0);
    //void update(size_t index_i, Real dt = 0.0);

};

/**
* @class ShearStressRelaxation
*/
class ShearStressRelaxation : public BaseRelaxation
{
public:

    explicit ShearStressRelaxation(BaseInnerRelation& inner_relation);
    virtual ~ShearStressRelaxation() {};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
protected:
    StdLargeVec<Matd>& shear_stress_, & shear_stress_rate_, & velocity_gradient_, & strain_tensor_, & strain_tensor_rate_;
    StdLargeVec<Real>& von_mises_stress_, & von_mises_strain_, & Vol_;
    StdLargeVec<Matd>& B_;
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
    FixedInAxisDirection(BodyPartByParticle& body_part, Vecd constrained_axises = Vecd::Zero());
    virtual ~FixedInAxisDirection() {};
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
    StdLargeVec<Vecd>& vel_;
    ReduceDynamics<QuantityMoment<Vecd>> compute_total_momentum_;

protected:
    virtual void setupDynamics(Real dt = 0.0) override;

public:
    explicit ConstrainSolidBodyMassCenter(SPHBody& sph_body, Vecd constrain_direction = Vecd::Ones());
    virtual ~ConstrainSolidBodyMassCenter() {};

    void update(size_t index_i, Real dt = 0.0);
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_DYNAMICS_INNER_H