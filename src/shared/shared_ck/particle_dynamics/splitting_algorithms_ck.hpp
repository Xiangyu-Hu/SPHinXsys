#ifndef SPLITTING_ALGORITHMS_CK_HPP
#define SPLITTING_ALGORITHMS_CK_HPP

#include "splitting_algorithms_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<Inner<Splitting, Parameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : LocalDynamicsType(std::forward<Args>(args)...),
      InteractionDynamicsCK<ExecutionPolicy, InteractionType<Base>>(),
      kernel_implementation_(*this)
{
    this->registerComputingKernel(&kernel_implementation_);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Inner<Splitting, Parameters...>>>::
    runInteraction(Real dt)
{
    InteractKernel *interact_kernel = kernel_implementation_.getComputingKernel();
    particle_for(
        LoopRangeCK<ExecutionPolicy, Splitting>(
            DynamicCast<CellLinkedList>(this, this->identifier_.getCellLinkedList())),
        [=](size_t i)
        { interact_kernel->interact(i, 0.5 * dt); }); // half time step for splitting
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Inner<Splitting, Parameters...>>>::
    exec(Real dt)
{
    this->setUpdated(this->identifier_.getSPHBody());
    this->setupDynamics(dt);
    InteractionDynamicsCK<Base>::runAllSteps(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Inner<Splitting, Parameters...>>>::
    runInteractionStep(Real dt)
{
    this->runInteraction(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
template <typename... Args>
InteractionDynamicsCK<ExecutionPolicy, InteractionType<Contact<Splitting, Parameters...>>>::
    InteractionDynamicsCK(Args &&...args)
    : LocalDynamicsType(std::forward<Args>(args)...),
      InteractionDynamicsCK<ExecutionPolicy, InteractionType<Base>>()
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
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Contact<Splitting, Parameters...>>>::
    runInteraction(Real dt)
{
    Real split_dt = dt * 0.5; // half time step for splitting
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        InteractKernel *interact_kernel =
            contact_kernel_implementation_[k]->getComputingKernel(k);

        particle_for(
            LoopRangeCK<ExecutionPolicy, Splitting>(
                DynamicCast<CellLinkedList>(this, this->identifier_.getCellLinkedList())),
            [=](size_t i)
            { interact_kernel->interact(i, split_dt); });
    }
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Contact<Splitting, Parameters...>>>::
    exec(Real dt)
{
    this->setUpdated(this->identifier_.getSPHBody());
    this->setupDynamics(dt);
    InteractionDynamicsCK<Base>::runAllSteps(dt);
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
void InteractionDynamicsCK<ExecutionPolicy, InteractionType<Contact<Splitting, Parameters...>>>::
    runInteractionStep(Real dt)
{
    this->runInteraction(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // SPLITTING_ALGORITHMS_CK_HPP
