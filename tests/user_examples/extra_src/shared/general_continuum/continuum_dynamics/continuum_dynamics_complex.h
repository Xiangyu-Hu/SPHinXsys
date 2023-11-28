#ifndef CONTINUUM_DYNAMICS_COMPLEX_H
#define CONTINUUM_DYNAMICS_COMPLEX_H

#include "base_fluid_dynamics.h"
#include "continuum_dynamics_inner.hpp"
#include "continuum_particles.h"
#include "general_continuum.h"

namespace SPH
{
namespace continuum_dynamics
{

typedef DataDelegateContact<ContinuumParticles, SolidParticles, DataDelegateEmptyBase>
    ContinuumWallData;
typedef DataDelegateContact<ContinuumParticles, BaseParticles, DataDelegateEmptyBase>
    ContinuumContactData;
typedef DataDelegateContact<ContinuumParticles, SolidParticles> FSIContactData;

typedef DataDelegateContact<PlasticContinuumParticles, SolidParticles, DataDelegateEmptyBase>
    PlasticContinuumWallData;
typedef DataDelegateContact<PlasticContinuumParticles, BaseParticles, DataDelegateEmptyBase>
    PlasticContinuumContactData;
typedef DataDelegateContact<PlasticContinuumParticles, SolidParticles> PlasticFSIContactData;

/**
 * @class InteractionWithWall
 * @brief Base class adding interaction with wall to general relaxation process
 */
template <class BaseIntegrationType>
class InteractionWithWall : public BaseIntegrationType, public FSIContactData
{
  public:
    template <class BaseBodyRelationType, typename... Args>
    InteractionWithWall(BaseBodyRelationType &base_body_relation, BaseContactRelation &wall_contact_relation, Args &&...args);
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<Real> wall_inv_rho0_;
    StdVec<StdLargeVec<Real> *> wall_mass_;
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_force_ave_, wall_acc_ave_, wall_n_;
};

/**
 * @class BaseShearStressRelaxation1stHalfType
 */
template <class BaseShearStressRelaxation1stHalfType>
class BaseShearStressRelaxation1stHalfWithWall : public InteractionWithWall<BaseShearStressRelaxation1stHalfType>
{
  public:
    template <typename... Args>
    BaseShearStressRelaxation1stHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseShearStressRelaxation1stHalfType>(std::forward<Args>(args)...){};
    virtual ~BaseShearStressRelaxation1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    // virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
};

using ShearStressRelaxation1stHalfWithWall = BaseShearStressRelaxation1stHalfWithWall<ShearStressRelaxation1stHalf>;

/**
 * @class BaseShearStressRelaxation2ndHalfType
 */
template <class BaseShearStressRelaxation2ndHalfType>
class BaseShearStressRelaxation2ndHalfWithWall : public InteractionWithWall<BaseShearStressRelaxation2ndHalfType>
{
  public:
    template <typename... Args>
    BaseShearStressRelaxation2ndHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseShearStressRelaxation2ndHalfType>(std::forward<Args>(args)...){};
    virtual ~BaseShearStressRelaxation2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

using ShearStressRelaxation2ndHalfWithWall = BaseShearStressRelaxation2ndHalfWithWall<ShearStressRelaxation2ndHalf>;

template <class BaseStressRelaxation1stHalfType>
class BaseStressRelaxation1stHalfWithWall : public InteractionWithWall<BaseStressRelaxation1stHalfType>
{
  public:
    template <typename... Args>
    BaseStressRelaxation1stHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseStressRelaxation1stHalfType>(std::forward<Args>(args)...){};
    virtual ~BaseStressRelaxation1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    virtual Vecd computeNonConservativeForce(size_t index_i) override;
};
using StressRelaxation1stHalfWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalf>;
using StressRelaxation1stHalfRiemannWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalfRiemann>;
using StressRelaxation1stHalfDissipativeRiemannWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalfDissipativeRiemann>;

template <class BaseStressRelaxation2ndHalfType>
class BaseStressRelaxation2ndHalfWithWall : public InteractionWithWall<BaseStressRelaxation2ndHalfType>
{
  public:
    template <typename... Args>
    BaseStressRelaxation2ndHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseStressRelaxation2ndHalfType>(std::forward<Args>(args)...){};
    virtual ~BaseStressRelaxation2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using StressRelaxation2ndHalfWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalf>;
using StressRelaxation2ndHalfRiemannWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalfRiemann>;
using StressRelaxation2ndHalfDissipativeRiemannWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalfDissipativeRiemann>;

/**
 * @class StressDiffusionWithWall
 */
template <class BaseStressDiffusionType>
class BaseStressDiffusionWithWall : public InteractionWithWall<BaseStressDiffusionType>
{
  public:
    template <typename... Args>
    BaseStressDiffusionWithWall(Args &&...args)
        : InteractionWithWall<BaseStressDiffusionType>(std::forward<Args>(args)...){};
    virtual ~BaseStressDiffusionWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using StressDiffusionWithWall = BaseStressDiffusionWithWall<StressDiffusion>;
} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_COMPLEX_H
