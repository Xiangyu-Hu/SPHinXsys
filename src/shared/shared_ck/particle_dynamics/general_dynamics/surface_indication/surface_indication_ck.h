#ifndef SURFACE_INDICATION_CK_H
#define SURFACE_INDICATION_CK_H

#include "base_general_dynamics.h"
#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{


template <typename...>
class IndicatorCK;
/**
 * @class FreeSurfaceIndicationCK
 * @brief A free-surface indication class for specialized relation types.
 */

template <typename... RelationTypes>
class FreeSurfaceIndicationCK;

/**
 * @class IndicatorCK
 * @brief This class template can be specialized for different relation types. 
 * Currently, this is not useful as this class was firstly construct for scenario with free surface + spatial.
 * Can be extend further once other type of surface indication is required. 
 * 
 */
template <>
class IndicatorCK<>
{
public:
    IndicatorCK(BaseParticles *particles) {}

    class ComputingKernel
    {
    public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        IndicatorCK<> &encloser,
                        ComputingKernelType &computing_kernel)
            : threshold_by_dimensions_(0.75 * Dimensions) {}

        Real operator()(Real &surface_indicator) { return surface_indicator; }

    protected:
        Real threshold_by_dimensions_;
    };
};

template<>
class IndicatorCK<Internal> : public IndicatorCK<> {
public:
    IndicatorCK(BaseParticles *particles) : IndicatorCK<>(particles) {}
    using ComputingKernel = typename IndicatorCK<>::ComputingKernel;
};

template <template <typename...> class RelationType, typename... Parameters>
class FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
public:
    template <class BaseRelationType>
    explicit FreeSurfaceIndicationCK(BaseRelationType &base_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

        virtual void interact(size_t index_i, Real dt = 0.0) {}

    protected:
        int *indicator_;
        Real *pos_div_, *Vol_;
        Real threshold_by_dimensions_;
        Real smoothing_length_;
    };

protected:
    DiscreteVariable<int> *dv_indicator_;
    DiscreteVariable<Real> *dv_pos_div_, *dv_Vol_;
    Real dv_threshold_by_dimensions_;
    Real dv_smoothing_length_;
};

template <class FlowType, typename... Parameters>
class FreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
    using IndicationKernel = typename IndicatorCK<FlowType>::ComputingKernel;

public:
    explicit FreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~FreeSurfaceIndicationCK() {}

    class InteractKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       FreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0) override;
        int *previous_surface_indicator_;
    };

    class UpdateKernel : public FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     FreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser);

        void update(size_t index_i, Real dt = 0.0);

    protected:
        IndicationKernel indication_;
        int *previous_surface_indicator_;
        FreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>* outer_;
    };

protected:
    IndicatorCK<FlowType> indicator_method_;
    DiscreteVariable<int> *dv_previous_surface_indicator_;
};

using FreeSurfaceIndicationInnerCK = FreeSurfaceIndicationCK<Inner<WithUpdate, Internal>>;

template <typename... Parameters>
class FreeSurfaceIndicationCK<Contact<Parameters...>>
    : public FreeSurfaceIndicationCK<Base, Contact<Parameters...>>
{
  public:
    explicit FreeSurfaceIndicationCK(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~FreeSurfaceIndicationCK(){};

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

using FreeSurfaceIndicationComplexCK = FreeSurfaceIndicationCK<Inner<WithUpdate, Internal>, Contact<>>;


} // namespace fluid_dynamics
} // namespace SPH

#endif // SURFACE_INDICATION_CK_H
