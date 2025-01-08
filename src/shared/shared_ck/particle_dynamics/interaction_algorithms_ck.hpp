#ifndef INTERACTION_ALGORITHMS_CK_HPP
#define INTERACTION_ALGORITHMS_CK_HPP

#include "interaction_algorithms_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Inner<Parameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : InteractionType<Inner<Parameters...>>(std::forward<Args>(args)...),
      kernel_implementation_(*this)
{
    this->registerComputingKernel(&kernel_implementation_);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Inner<Parameters...>>>::
    runInteraction(Real dt)
{
    InteractKernel *interact_kernel = kernel_implementation_.getComputingKernel();
    particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(this->identifier_),
                 [=](size_t i)
                 { interact_kernel->interact(i, dt); });
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Contact<Parameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : InteractionType<Contact<Parameters...>>(std::forward<Args>(args)...)
{
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        contact_kernel_implementation_.push_back(
            contact_kernel_implementation_ptrs_
                .template createPtr<KernelImplementation>(*this));
        this->registerComputingKernel(contact_kernel_implementation_.back(), k);
    }
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Contact<Parameters...>>>::
    runInteraction(Real dt)
{
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        InteractKernel *interact_kernel =
            contact_kernel_implementation_[k]->getComputingKernel(k);

        particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(this->identifier_),
                     [=](size_t i)
                     { interact_kernel->interact(i, dt); });
    }
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... Parameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<Parameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : InteractionDynamicsCK<
          ExecutionPolicy, Base,
          InteractionType<RelationType<Parameters...>>>(std::forward<Args>(args)...),
      InteractionDynamicsCK<Base>(), BaseDynamics<void>() {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<Parameters...>>>::
    exec(Real dt)
{
    this->setUpdated(this->identifier_.getSPHBody());
    this->setupDynamics(dt);
    InteractionDynamicsCK<Base>::runAllSteps(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<Parameters...>>>::
    runInteractionStep(Real dt)
{
    this->runInteraction(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<WithUpdate, OtherParameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : InteractionDynamicsCK<
          ExecutionPolicy, Base, InteractionType<RelationType<WithUpdate, OtherParameters...>>>(
          std::forward<Args>(args)...),
      InteractionDynamicsCK<WithUpdate>(),
      BaseDynamics<void>(),
      kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<WithUpdate, OtherParameters...>>>::
    exec(Real dt)
{
    this->setUpdated(this->identifier_.getSPHBody());
    this->setupDynamics(dt);
    InteractionDynamicsCK<WithUpdate>::runAllSteps(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<WithUpdate, OtherParameters...>>>::
    runInteractionStep(Real dt)
{
    this->runInteraction(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<WithUpdate, OtherParameters...>>>::
    runUpdateStep(Real dt)
{
    UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
    particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(this->identifier_),
                 [=](size_t i)
                 { update_kernel->update(i, dt); });
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<RelationType<OneLevel, OtherParameters...>>>(
          std::forward<Args>(args)...),
      InteractionDynamicsCK<OneLevel>(), BaseDynamics<void>(),
      initialize_kernel_implementation_(*this), update_kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>::
    exec(Real dt)
{
    this->setUpdated(this->identifier_.getSPHBody());
    this->setupDynamics(dt);
    InteractionDynamicsCK<OneLevel>::runAllSteps(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>::
    runInteractionStep(Real dt)
{
    this->runInteraction(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>::
    runInitializationStep(Real dt)
{
    InitializeKernel *initialize_kernel = initialize_kernel_implementation_.getComputingKernel();
    particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(this->identifier_),
                 [=](size_t i)
                 { initialize_kernel->initialize(i, dt); });
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>::
    runUpdateStep(Real dt)
{
    UpdateKernel *update_kernel = update_kernel_implementation_.getComputingKernel();
    particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(this->identifier_),
                 [=](size_t i)
                 { update_kernel->update(i, dt); });
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          class FirstInteraction, class... Others>
template <class FirstParameterSet, typename... OtherParameterSets>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction, Others...>>::
    InteractionDynamicsCK(FirstParameterSet &&first_parameter_set,
                          OtherParameterSets &&...other_parameter_sets)
    : InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction>>(first_parameter_set),
      other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType,
          class FirstInteraction, class... Others>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction, Others...>>::
    runInteractionStep(Real dt)
{
    InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction>>::runInteractionStep(dt);
    other_interactions_.runInteractionStep(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_ALGORITHMS_CK_HPP
