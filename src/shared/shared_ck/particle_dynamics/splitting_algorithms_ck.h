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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	splitting_algorithms_ck.h
 * @brief 	TBD
 * @author	Xiangyu Hu
 */

#ifndef SPLITTING_ALGORITHMS_CK_H
#define SPLITTING_ALGORITHMS_CK_H

#include "interaction_algorithms_ck.h"

namespace SPH
{
template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType<Inner<Splitting, Parameters...>>>
    : public InteractionType<Inner<Splitting, Parameters...>>,
      public InteractionDynamicsCK<ExecutionPolicy, InteractionType<Base>>
{
    using LocalDynamicsType = InteractionType<Inner<Splitting, Parameters...>>;
    using Identifier = typename LocalDynamicsType::Identifier;
    using InteractKernel = typename LocalDynamicsType::InteractKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};
    virtual void exec(Real dt = 0.0) override;
    virtual void runInteractionStep(Real dt = 0.0) override;

  protected:
    void runInteraction(Real dt);
};
} // namespace SPH
#endif // SPLITTING_ALGORITHMS_CK_H
