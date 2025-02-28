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
    DiscreteVariable<int> *dv_previous_surface_indicator_;
    Real dv_threshold_by_dimensions_;
    Real dv_smoothing_length_;
};

//=================================================================================================//
// Specialization for Inner<WithUpdate, Internal, Parameters...>
//=================================================================================================//
template <typename... Parameters>
class FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    /**
     * @class InteractKernel
     * @brief Interaction kernel for detecting free surface in an inner relation.
     */
    class InteractKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0);

        /// Pointer to the previously stored surface indicator.
        int *previous_surface_indicator_;
    };

    /**
     * @class UpdateKernel
     * @brief Update kernel that modifies the surface indicator based on the newly computed
     *        positional divergences and neighbor information.
     */
    class UpdateKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>> &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        int *previous_surface_indicator_;
    };

  protected:
    DiscreteVariable<int> *dv_previous_surface_indicator_;
};

//=================================================================================================//
// Specialization for Inner<Base, Internal, Parameters...>
//=================================================================================================//
/**
 * @brief Specialization for inner relations without update.
 *
 * This version omits any update kernel, so that the update step is not executed.
 */
template <typename... Parameters>
class FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    // Provide only an InteractKernel (no UpdateKernel)
    class InteractKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0);
    };
};

//=================================================================================================//
// Specialization for Inner<WithUpdate, SpatialTemporal, Parameters...>
//=================================================================================================//
template <typename... Parameters>
class FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    /**
     * @class InteractKernel
     * @brief Interaction kernel for detecting free surface in a spatial-temporal inner relation.
     */
    class InteractKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0);

        /// Pointer to the previously stored surface indicator.
        int *previous_surface_indicator_;
        int *is_near_surface_indicator_;
    };

    /**
     * @class UpdateKernel
     * @brief Update kernel that modifies the surface indicator based on the newly computed
     *        positional divergences and neighbor information.
     */
    class UpdateKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>> &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        int *previous_surface_indicator_;
        int *is_near_surface_indicator_;
    };

  protected:
    DiscreteVariable<int> *dv_previous_surface_indicator_;
    DiscreteVariable<int> *dv_is_near_previous_surface_indicator_;
    DiscreteVariable<int> *dv_is_near_surface_indicator_;
};

//=================================================================================================//
// Specialization for Contact<Parameters...>
//=================================================================================================//
template <typename... Parameters>
class FreeSurfaceIndicationCK<Contact<Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Contact<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    /**
     * @class InteractKernel
     * @brief Interaction kernel for a contact relation, combining data from multiple contact bodies.
     */
    class InteractKernel : public FreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel
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
// Common type aliases for free-surface indication.
using FreeSurfaceIndicationInnerCK = FreeSurfaceIndicationCK<Inner<WithUpdate, FreeSurface>>;
using FreeSurfaceIndicationInnerSpatialTemporalCK = FreeSurfaceIndicationCK<Inner<WithUpdate, FreeSurface>, Inner<WithUpdate, SpatialTemporal>>;
using FreeSurfaceIndicationComplexCK = FreeSurfaceIndicationCK<Inner<WithUpdate, Internal>, Contact<>>;
using FreeSurfaceIndicationComplexSpatialTemporalCK = FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal>, Contact<>>;

class SurfaceIndicationByAlignedBoxCK : public BaseLocalDynamics<AlignedBoxPartByCell>
{

  public:
    SurfaceIndicationByAlignedBoxCK(AlignedBoxPartByCell &aligned_box_part)
        : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
          sv_aligned_box_(aligned_box_part.svAlignedBox()),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_indicator_(particles_->template registerStateVariableOnly<int>("Indicator")),
          previous_surface_indicator_(particles_->template registerStateVariableOnly<int>("PreviousSurfaceIndicator"))
    {
    }
    virtual ~SurfaceIndicationByAlignedBoxCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0); // only works in sequenced policy

      protected:
        AlignedBox *aligned_box_;
        Vecd *pos_;
        int *indicator_;
    };

  protected:
    SingularVariable<AlignedBox> *sv_aligned_box_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_indicator_;
    DiscreteVariable<int> *previous_surface_indicator_;
};
} // namespace fluid_dynamics
} // namespace SPH

#endif // SURFACE_INDICATION_CK_H
