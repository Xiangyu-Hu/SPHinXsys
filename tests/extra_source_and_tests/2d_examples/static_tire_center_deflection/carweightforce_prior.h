#ifndef CARWEIGHTFORCE_PRIOR_H
#define CARWEIGHTFORCE_PRIOR_H

#include "base_general_dynamics.h"

namespace SPH
{

template <class DynamicsIdentifier>
class BaseForcePrior1 : public BaseLocalDynamics<DynamicsIdentifier>
{
  protected:
    Vecd *force_prior_, *current_force_, *previous_force_;

  public:
    BaseForcePrior1(DynamicsIdentifier &identifier, const std::string &force_name);
    virtual ~BaseForcePrior1(){};
    void update(size_t index_i, Real dt = 0.0);
};
 using ForcePrior2 = BaseForcePrior1<BodyPartByParticle>;

template <class CarweightForcetype>
class CarweightForce : public ForcePrior2
{
  protected:
    const CarweightForcetype gravity_;
    Vecd *pos_;
    Real *mass_;
    Real *physical_time_;

  public:
    CarweightForce(BodyPartByParticle &sph_body, const CarweightForcetype &gravity);
    virtual ~CarweightForce(){};
    void update(size_t index_i, Real dt = 0.0);
};
} 
#endif 
