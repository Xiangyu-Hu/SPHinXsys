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
typedef DataDelegateContact<ContinuumParticles, SolidParticles> FSIContactData;
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
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_force_ave_, wall_n_;
};

/**
 * @class Plastic:BaseStressRelaxationWithWall
 */
template <class BaseStressRelaxation1stHalfType>
class BaseStressRelaxation1stHalfWithWall : public InteractionWithWall<BaseStressRelaxation1stHalfType>
{
public:
    template <typename... Args>
    BaseStressRelaxation1stHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseStressRelaxation1stHalfType>(std::forward<Args>(args)...) {};
    virtual ~BaseStressRelaxation1stHalfWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual Vecd computeNonConservativeForce(size_t index_i) override;
};
using StressRelaxation1stHalfWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalf>;
using StressRelaxation1stHalfRiemannWithWall = BaseStressRelaxation1stHalfWithWall<StressRelaxation1stHalfRiemann>;

template <class BaseStressRelaxation2ndHalfType>
class BaseStressRelaxation2ndHalfWithWall : public InteractionWithWall<BaseStressRelaxation2ndHalfType>
{
public:
    template <typename... Args>
    BaseStressRelaxation2ndHalfWithWall(Args &&...args)
        : InteractionWithWall<BaseStressRelaxation2ndHalfType>(std::forward<Args>(args)...) {};
    virtual ~BaseStressRelaxation2ndHalfWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
using StressRelaxation2ndHalfWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalf>;
using StressRelaxation2ndHalfRiemannWithWall = BaseStressRelaxation2ndHalfWithWall<StressRelaxation2ndHalfRiemann>;

} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_COMPLEX_H
