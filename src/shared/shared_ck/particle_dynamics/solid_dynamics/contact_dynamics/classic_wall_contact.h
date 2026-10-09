// SPDX-License-Identifier: Apache-2.0
#ifndef CLASSIC_WALL_CONTACT_H
#define CLASSIC_WALL_CONTACT_H

#include "force_prior_ck.hpp"
#include "interaction_ck.hpp"
#include "base_body_part.h"

namespace SPH
{
namespace solid_dynamics
{
template <typename...>
class ClassicWallContactForceCK;

/**
 * The ContactFactorSummation / ContactForceFromWall law for one fixed wall.
 * Supply a BodySurfaceLayer constructed with its default three-spacing thickness
 * before deformation, and a contact relation for the full source body.
 * Both bodies must use ordinary volume particles and uniform, isotropic Wendland C2
 * adaptations; shell contact is not supported. Several independent wall instances
 * would not reproduce the classic multi-wall pressure sum, and are rejected for
 * the same source body.
 */
template <typename... Parameters>
class ClassicWallContactForceCK<Contact<WithUpdate, Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using BaseInteractionType = Interaction<Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    ClassicWallContactForceCK(DynamicsIdentifier &identifier, BodySurfaceLayer &surface);

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real inv_h_, factor_W_, factor_dW_, offset_W_, cutoff_, stiffness_;
        Real *Vol_, *contact_Vol_, *contact_factor_;
        Vecd *contact_force_;
        int *surface_mask_;
    };

  protected:
    Real inv_h_, factor_W_, factor_dW_, offset_W_, cutoff_, stiffness_;
    DiscreteVariable<Real> *dv_contact_factor_;
    DiscreteVariable<int> *dv_surface_mask_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // CLASSIC_WALL_CONTACT_H
