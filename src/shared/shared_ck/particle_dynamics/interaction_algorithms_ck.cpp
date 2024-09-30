#include "interaction_algorithms_ck.h"

namespace SPH
{
//=================================================================================================//
void InteractionDynamicsCK<Base>::runAllSteps(Real dt)
{
    for (size_t k = 0; k < this->pre_processes_.size(); ++k)
        this->pre_processes_[k]->exec(dt);

    runInteractionStep(dt);

    for (size_t k = 0; k < this->post_processes_.size(); ++k)
        this->post_processes_[k]->exec(dt);
}
//=================================================================================================//
void InteractionDynamicsCK<WithUpdate>::runAllSteps(Real dt)
{
    InteractionDynamicsCK<Base>::runAllSteps(dt);
    runUpdateStep(dt);
}
//=================================================================================================//
void InteractionDynamicsCK<WithInitialization>::runAllSteps(Real dt)
{
    runInitializationStep(dt);
    InteractionDynamicsCK<Base>::runAllSteps(dt);
}
//=================================================================================================//
void InteractionDynamicsCK<OneLevel>::runAllSteps(Real dt)
{
    runInitializationStep(dt);
    InteractionDynamicsCK<Base>::runAllSteps(dt);
    runUpdateStep(dt);
}
//=================================================================================================//
} // namespace SPH
