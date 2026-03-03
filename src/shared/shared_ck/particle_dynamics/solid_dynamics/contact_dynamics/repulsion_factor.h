/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file repulsion_factor.h
 * @brief TBD.
 * @details TBD.
 * @author Xiangyu Hu
 */

#ifndef REPULSION_FACTOR_H
#define REPULSION_FACTOR_H

#include "interaction_ck.hpp"

namespace SPH
{
namespace solid_dynamics
{

template <typename...>
class RepulsionFactor;

template <typename... Parameters>
class RepulsionFactor<Base, Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteractionType = Interaction<Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit RepulsionFactor(DynamicsIdentifier &identifier, const std::string &factor_name);
    virtual ~RepulsionFactor() {};

  protected:
    DiscreteVariable<Real> *dv_repulsion_factor_;
};

template <typename... Parameters>
class RepulsionFactor<Contact<Parameters...>>
    : public RepulsionFactor<Base, Contact<Parameters...>>
{
    using BaseInteractionType = RepulsionFactor<Base, Contact<Parameters...>>;

  public:
    explicit RepulsionFactor(Contact<Parameters...> &contact_relation);
    virtual ~RepulsionFactor() {};

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *repulsion_factor_;
        Real contact_inv_rho0_;
        Real *contact_mass_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FACTOR_H
