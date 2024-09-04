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
 * @file 	complex_algorithms_ck.h
 * @brief 	This is the classes for algorithms particle dynamics .
 * @detail	TBD
 * @author	Xiangyu Hu
 */

#ifndef COMPLEX_ALGORITHMS_CK_H
#define COMPLEX_ALGORITHMS_CK_H

#include "dynamics_algorithms_ck.h"
#include "interaction_dynamics_algorithms_ck.h"

namespace SPH
{
template <typename... T>
class SequencedCombination;

template <class ExecutionPolicy, typename... CommonParameters,
          template <typename... LocalDynamicsType> class AlgorithmType,
          template <typename... InteractionTypes> class LocalDynamicsName>
class SequencedCombination<AlgorithmType<ExecutionPolicy, LocalDynamicsName<>, CommonParameters...>> : public BaseDynamics<void>
{
  public:
    SequencedCombination() : BaseDynamics<void>(){};
    virtual void exec(Real dt = 0.0) override{};
};
template <class ExecutionPolicy, typename... CommonParameters,
          template <typename... LocalDynamicsType> class AlgorithmType,
          template <typename... InteractionTypes> class LocalDynamicsName,
          class FirstInteraction, class... OtherInteractions>
class SequencedCombination<AlgorithmType<ExecutionPolicy, LocalDynamicsName<FirstInteraction, OtherInteractions...>, CommonParameters...>>
    : public AlgorithmType<ExecutionPolicy, LocalDynamicsName<FirstInteraction, CommonParameters...>>
{
  protected:
    SequencedCombination<AlgorithmType<ExecutionPolicy, LocalDynamicsName<OtherInteractions...>, CommonParameters...>> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit SequencedCombination(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : AlgorithmType<ExecutionPolicy, LocalDynamicsName<FirstInteraction, CommonParameters...>>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    virtual void exec(Real dt = 0.0) override
    {
        AlgorithmType<ExecutionPolicy, LocalDynamicsName<FirstInteraction, CommonParameters...>>::exec(dt);
        other_interactions_.exec(dt);
    };
};

template <typename... T>
class ComplexInteraction;

template <class ExecutionPolicy, typename... CommonParameters,
          template <typename... InteractionTypes> class LocalInteractionName>
class ComplexInteraction<InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<>, CommonParameters...>>
{
  public:
    ComplexInteraction(){};
    void runMainStep(Real dt){};
};

template <class ExecutionPolicy, typename... CommonParameters,
          template <typename... InteractionTypes> class LocalInteractionName,
          class FirstInteraction, class... OtherInteractions>
class ComplexInteraction<InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<FirstInteraction, OtherInteractions...>, CommonParameters...>>
    : public InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<FirstInteraction, CommonParameters...>>
{
  protected:
    ComplexInteraction<InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<OtherInteractions...>, CommonParameters...>> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit ComplexInteraction(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<FirstInteraction, CommonParameters...>>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        InteractionDynamicsCK<ExecutionPolicy, LocalInteractionName<FirstInteraction, CommonParameters...>>::exec(dt);
        other_interactions_.exec(dt);
    };
};
} // namespace SPH
#endif // COMPLEX_ALGORITHMS_CK_H
