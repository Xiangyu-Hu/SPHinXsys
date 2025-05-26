#ifndef SURFACE_INDICATION_CK_H
#define SURFACE_INDICATION_CK_H

#include "base_fluid_dynamics.h"
#include "base_general_dynamics.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FreeSurfaceIndicationCK
 * @brief Free-surface indication with configurable relationship types.
 *
 * This template is specialized for different combinations of
 * "Base" + "Inner" / "Contact" relations to handle free-surface
 * detection and updating in SPH simulations.
 */
template <typename... RelationTypes>
class FreeSurfaceIndicationCK;

//=================================================================================================//
// Base relation version
//=================================================================================================//
/**
 * @class FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
 * @brief Basic free-surface indication logic, storing and computing
 *        positional divergence, surface indicators, etc.
 */
template <template <typename...> class RelationType, typename... Parameters>
class FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit FreeSurfaceIndicationCK(BaseRelationType &base_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    //------------------------------------------------------------------------------------------//
    /**
     * @class InteractKernel
     * @brief Base interaction kernel that calculates positional divergence,
     *        used in detecting free-surface particles.
     */
    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        int *indicator_;
        Real *pos_div_;
        Real *Vol_;
        Real threshold_by_dimensions_;
        Real smoothing_length_;
    };

  protected:
    DiscreteVariable<int> *dv_indicator_;
    DiscreteVariable<Real> *dv_pos_div_;
    DiscreteVariable<Real> *dv_Vol_;
    Real dv_threshold_by_dimensions_;
    Real dv_smoothing_length_;
};

//=================================================================================================//
// Inner relation version with "WithUpdate"
//=================================================================================================//
/**
 * @class FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>
 * @brief Extends the base free-surface indication for inner relations that also require
 *        an updating step (e.g., "WithUpdate").
 */
template <typename... Parameters>
class FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Inner<Parameters...> &inner_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    //------------------------------------------------------------------------------------------//
    /**
     * @class InteractKernel
     * @brief Interaction kernel for detecting free surface in an inner relation.
     */
    class InteractKernel
        : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0);

        /// Pointer to the previously stored surface indicator.
        int *previous_surface_indicator_;
    };

    //------------------------------------------------------------------------------------------//
    /**
     * @class UpdateKernel
     * @brief Post-processing/update kernel that modifies the surface indicator
     *        based on the newly computed positional divergences and neighbor information.
     */
    class UpdateKernel
        : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>> &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        int *previous_surface_indicator_;
        FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>> *outer_;
    };

  protected:
    DiscreteVariable<int> *dv_previous_surface_indicator_;
};
//=================================================================================================//
// Contact relation version
//=================================================================================================//
/**
 * @class FreeSurfaceIndicationCK<Contact<Parameters...>>
 * @brief Extends the base free-surface indication for contact relations.
 */
template <typename... Parameters>
class FreeSurfaceIndicationCK<Contact<Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Contact<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Contact<Parameters...> &contact_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    //------------------------------------------------------------------------------------------//
    /**
     * @class InteractKernel
     * @brief Interaction kernel for a contact relation, combining data from multiple contact bodies.
     */
    class InteractKernel
        : public FreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Contact<Parameters...>> &encloser,
                       size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

//------------------------------------------------------------------------------------------//
// Common type alias for complex free-surface indication (inner + contact).
using FreeSurfaceIndicationComplexSpatialTemporalCK = FreeSurfaceIndicationCK<Inner<WithUpdate>, Contact<>>;
//------------------------------------------------------------------------------------------//

} // namespace fluid_dynamics
} // namespace SPH

#endif // SURFACE_INDICATION_CK_H
