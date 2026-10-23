/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPDX-License-Identifier: Apache-2.0
 * ------------------------------------------------------------------------- */
/**
 * @file  near_wall_boundary.h
 * @brief Shared boundary condition implementations for near-wall interactions  .
 * @author sphinxsim@einsimo.com
 */

#ifndef NEAR_WALL_BOUNDARY_H
#define NEAR_WALL_BOUNDARY_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{
template <typename...>
class DistanceFromWallCK;

template <typename... Parameters>
class DistanceFromWallCK<Contact<WithInitialization, Parameters...>>
    : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit DistanceFromWallCK(Contact<Parameters...> &wall_contact_relation);
    virtual ~DistanceFromWallCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        Real distance_default_;
        Real *phi_;
        Vecd *n_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real spacing_ref_square_;
        Vecd *pos_;
        Vecd *n_, *wall_pos_, *wall_n_;
        Real *phi_, *wall_phi_;
    };

  protected:
    Real spacing_ref_, distance_default_, h_ref_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_n_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_wall_pos_, *dv_wall_n_;
    DiscreteVariable<Real> *dv_wall_phi_;
};

template <typename...>
class BoundingFromWallCK;

template <typename... Parameters>
class BoundingFromWallCK<Contact<OneLevel, Parameters...>>
    : public DistanceFromWallCK<Contact<WithInitialization, Parameters...>>
{
    using BaseInteraction = DistanceFromWallCK<Contact<WithInitialization, Parameters...>>;

  public:
    explicit BoundingFromWallCK(Contact<Parameters...> &wall_contact_relation);
    virtual ~BoundingFromWallCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real spacing_ref_, distance_min_;
        Real *phi_;
        Vecd *pos_, *n_;
    };

  protected:
    Real distance_min_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // NEAR_WALL_BOUNDARY_H
