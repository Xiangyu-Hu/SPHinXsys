/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPDX-License-Identifier: Apache-2.0
 * ------------------------------------------------------------------------- */
/**
 * @file  subscale_buoyancy.h
 * @brief Implementation of buoyancy force for fluid dynamics.
 * @author sphinxsim@einsimo.com
 */

#ifndef SUBSCALE_BUOYANCY_H
#define SUBSCALE_BUOYANCY_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{
template <typename...>
class SubscaleBuoyancy;

template <typename... Parameters>
class SubscaleBuoyancy<Contact<WithUpdate, Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using BaseDynamicsType = Interaction<Contact<Parameters...>>;

  public:
    SubscaleBuoyancy(Contact<Parameters...> &water_air_relation, const Gravity &gravity, Real threshold = 0.3);
    virtual ~SubscaleBuoyancy() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Gravity gravity_;
        Real threshold_by_dimensions_;
        Real *Vol_, *mass_, *isolation_factor_;
        Vecd *subscale_buoyancy_;
    };

  protected:
    const Gravity gravity_;
    Real threshold_by_dimensions_;
    DiscreteVariable<Real> *dv_mass, *dv_isolation_factor_;
    DiscreteVariable<Vecd> *dv_subscale_buoyancy_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // SUBSCALE_BUOYANCY_H
