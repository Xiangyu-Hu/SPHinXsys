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

#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"

namespace SPH
{
template <typename... T>
class DynamicsSequence;  // same DynamicsType with the different constructor arguments

template <class ExecutionPolicy,
          template <typename... LocalDynamicsType> class DynamicsType>
class DynamicsSequence<DynamicsType<ExecutionPolicy>> : public BaseDynamics<void>
{
  public:
    DynamicsSequence() : BaseDynamics<void>() {};
    virtual void exec(Real dt = 0.0) override {};
};

template <class ExecutionPolicy,
          template <typename... LocalDynamicsType> class DynamicsType,
          class FirstLocalDynamics, class... OtherLocalDynamics>
class DynamicsSequence<DynamicsType<ExecutionPolicy, FirstLocalDynamics, OtherLocalDynamics...>>
    : public DynamicsType<ExecutionPolicy, FirstLocalDynamics>
{
  protected:
    DynamicsSequence<DynamicsType<ExecutionPolicy, OtherLocalDynamics...>> other_dynamics_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit DynamicsSequence(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : DynamicsType<ExecutionPolicy, FirstLocalDynamics>(first_parameter_set),
          other_dynamics_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    virtual void exec(Real dt = 0.0) override
    {
        DynamicsType<ExecutionPolicy, FirstLocalDynamics>::exec(dt);
        other_dynamics_.exec(dt);
    };
};

template <typename... T>
class ArbitraryDynamicsSequence; // different DynamicsType with the different constructor arguments

template <>
class ArbitraryDynamicsSequence<> : public BaseDynamics<void>
{
  public:
    ArbitraryDynamicsSequence() : BaseDynamics<void>() {};
    virtual void exec(Real dt = 0.0) override {};
};

template <class FistDynamicsType, class... OtherDynamicsTypes>
class ArbitraryDynamicsSequence<FistDynamicsType, OtherDynamicsTypes...> : public FistDynamicsType
{
  protected:
    ArbitraryDynamicsSequence<OtherDynamicsTypes...> other_dynamics_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit ArbitraryDynamicsSequence(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : FistDynamicsType(first_parameter_set),
          other_dynamics_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    virtual void exec(Real dt = 0.0) override
    {
        FistDynamicsType::exec(dt);
        other_dynamics_.exec(dt);
    };
};

template <typename... T>
class RungeKuttaSequence; // same DynamicsType with the same constructor arguments

template <class ExecutionPolicy,
          template <typename... LocalDynamicsType> class DynamicsType>
class RungeKuttaSequence<DynamicsType<ExecutionPolicy>> : public BaseDynamics<void>
{
  public:
    template <typename... Args>
    RungeKuttaSequence(Args &&...args) : BaseDynamics<void>(){};
    virtual void exec(Real dt = 0.0) override {};
};

template <class ExecutionPolicy,
          template <typename... LocalDynamicsType> class DynamicsType,
          class FirstLocalDynamics, class... OtherLocalDynamics>
class RungeKuttaSequence<DynamicsType<ExecutionPolicy, FirstLocalDynamics, OtherLocalDynamics...>>
    : public DynamicsType<ExecutionPolicy, FirstLocalDynamics>
{
  protected:
    RungeKuttaSequence<DynamicsType<ExecutionPolicy, OtherLocalDynamics...>> other_dynamics_;

  public:
    template <typename... Args>
    explicit RungeKuttaSequence(Args &&...args)
        : DynamicsType<ExecutionPolicy, FirstLocalDynamics>(std::forward<Args>(args)...),
          other_dynamics_(std::forward<Args>(args)...){};

    virtual void exec(Real dt = 0.0) override
    {
        DynamicsType<ExecutionPolicy, FirstLocalDynamics>::exec(dt);
        other_dynamics_.exec(dt);
    };
};
} // namespace SPH
#endif // COMPLEX_ALGORITHMS_CK_H
