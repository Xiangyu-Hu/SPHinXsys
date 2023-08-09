#ifndef CONTINUUM_DYNAMICS_COMPLEX_H
#define CONTINUUM_DYNAMICS_COMPLEX_H

#include "fluid_dynamics_complex.h"
#include "continuum_dynamics_inner.hpp"
#include "general_continuum.h"
#include "continuum_particles.h"

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
 * @class BaseShearStressRelaxation1stHalfType
 */
template <class BaseShearStressRelaxation1stHalfType>
class BaseShearStressRelaxation1stHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation1stHalfType>
{
  public:
    template <typename... Args>
    BaseShearStressRelaxation1stHalfWithWall(Args &&...args)
        : fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation1stHalfType>(std::forward<Args>(args)...){};
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
class BaseShearStressRelaxation2ndHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation2ndHalfType>
{
  public:
    template <typename... Args>
    BaseShearStressRelaxation2ndHalfWithWall(Args &&...args)
        : fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation2ndHalfType>(std::forward<Args>(args)...){};
    virtual ~BaseShearStressRelaxation2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

using ShearStressRelaxation2ndHalfWithWall = BaseShearStressRelaxation2ndHalfWithWall<ShearStressRelaxation2ndHalf>;

template <class BaseStressRelaxation1stHalfType>
class BaseStressRelaxation1stHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseStressRelaxation1stHalfType>
{
public:
    template <typename... Args>
    BaseStressRelaxation1stHalfWithWall(Args &&...args)
        : fluid_dynamics::InteractionWithWall<BaseStressRelaxation1stHalfType>(std::forward<Args>(args)...) {};
    virtual ~BaseStressRelaxation1stHalfWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
};
using StressRelaxation1stHalfWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalf>;
using StressRelaxation1stHalfRiemannWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalfRiemann>;
using StressRelaxation1stHalfDissipativeRiemannfWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalfDissipativeRiemann>;

template <class BaseStressRelaxation2ndHalfType>
class BaseStressRelaxation2ndHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseStressRelaxation2ndHalfType>
{
public:
    template <typename... Args>
    BaseStressRelaxation2ndHalfWithWall(Args &&...args)
        : fluid_dynamics::InteractionWithWall<BaseStressRelaxation2ndHalfType>(std::forward<Args>(args)...) {};
    virtual ~BaseStressRelaxation2ndHalfWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
using StressRelaxation2ndHalfWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalf>;
using StressRelaxation2ndHalfRiemannWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalfRiemann>;
using StressRelaxation2ndHalfDissipativeRiemannWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalfDissipativeRiemann>;

/**
* @class StressDiffusionWithWall
*/
template <class BaseStressDiffusionType>
class BaseStressDiffusionWithWall : public fluid_dynamics::InteractionWithWall<BaseStressDiffusionType>
{
public:
    template <typename... Args>
    BaseStressDiffusionWithWall(Args &&...args)
        : fluid_dynamics::InteractionWithWall<BaseStressDiffusionType>(std::forward<Args>(args)...) {};
    virtual ~BaseStressDiffusionWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
using StressDiffusionWithWall = BaseStressDiffusionWithWall<StressDiffusion>;
} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_COMPLEX_H
