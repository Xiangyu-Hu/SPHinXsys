#ifndef CONTINUUM_DYNAMICS_INNER_H
#define CONTINUUM_DYNAMICS_INNER_H
#include "constraint_dynamics.h"
#include "fluid_dynamics_inner.hpp"
#include "general_continuum.h"
#include "continuum_particles.h"
#include "riemann_solver_extra.h"

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
    StdLargeVec<Mat3d>& stress_tensor_3D_;
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

//=================================================================================================//
//================================================Plastic==========================================//
//=================================================================================================//

class BaseRelaxationPlastic : public LocalDynamics, public PlasticContinuumDataInner
{
public:
    explicit BaseRelaxationPlastic(BaseInnerRelation& inner_relation);
    virtual ~BaseRelaxationPlastic() {};
    Matd reduceTensor(Mat3d tensor_3d);
    Mat3d increaseTensor(Matd tensor_2d);

protected:
    PlasticContinuum& plastic_continuum_;
    StdLargeVec<Real>& rho_, & p_, & drho_dt_;
    StdLargeVec<Vecd>& pos_, & vel_, & acc_, & acc_prior_;
    StdLargeVec<Mat3d>& stress_tensor_3D_, & strain_tensor_3D_, & stress_rate_3D_, & strain_rate_3D_;
    StdLargeVec<Mat3d>& elastic_strain_tensor_3D_, & elastic_strain_rate_3D_;
};

//=================================================================================================//
 //===============================BaseStressRelaxation1stHalf======================================//
 //=================================================================================================//
 /**
  * @class BaseIntegration1stHalf
  * @brief Template class for pressure relaxation scheme with the Riemann solver
  * as template variable
  */
template <class RiemannSolverType>
class BaseStressRelaxation1stHalf : public BaseRelaxationPlastic
{
public:
    explicit BaseStressRelaxation1stHalf(BaseInnerRelation& inner_relation);
    virtual ~BaseStressRelaxation1stHalf() {};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:
    virtual Vecd computeNonConservativeAcceleration(size_t index_i);

    StdLargeVec<Matd>& velocity_gradient_;
};
using StressRelaxation1stHalf = BaseStressRelaxation1stHalf<NoRiemannSolver>;
using StressRelaxation1stHalfRiemann = BaseStressRelaxation1stHalf<AcousticRiemannSolverExtra>;
using StressRelaxation1stHalfDissipativeRiemann = BaseStressRelaxation1stHalf<DissipativeRiemannSolverExtra>;

//=================================================================================================//
//===============================BaseStressRelaxation2ndHalf======================================//
//=================================================================================================//
/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseStressRelaxation2ndHalf : public BaseRelaxationPlastic
{
public:
    explicit BaseStressRelaxation2ndHalf(BaseInnerRelation& inner_relation);
    virtual ~BaseStressRelaxation2ndHalf() {};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:

    StdLargeVec<Matd>& velocity_gradient_;
    StdLargeVec<Real>& Vol_, & mass_, & von_mises_stress_;
    StdLargeVec<Real>& acc_deviatoric_plastic_strain_, & vertical_stress_;
    Real E_, nu_;
};
using StressRelaxation2ndHalf = BaseStressRelaxation2ndHalf<NoRiemannSolver>;
using StressRelaxation2ndHalfRiemann = BaseStressRelaxation2ndHalf<AcousticRiemannSolverExtra>;
using StressRelaxation2ndHalfDissipativeRiemann = BaseStressRelaxation2ndHalf<DissipativeRiemannSolverExtra>;
/**
* @class StressDiffusion
*/
class StressDiffusion : public BaseRelaxationPlastic
{
public:
    explicit StressDiffusion(BaseInnerRelation& inner_relation);
    virtual ~StressDiffusion() {};
    void interaction(size_t index_i, Real dt = 0.0);
protected:
    Real zeta_ = 0.1, fai_; //diffusion coeficient
    Real smoothing_length_, sound_speed_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_DYNAMICS_INNER_H