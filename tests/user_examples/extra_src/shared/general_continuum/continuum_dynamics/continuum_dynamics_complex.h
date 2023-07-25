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
} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_COMPLEX_H
