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
 * @file 	particle_dynamics_dissipation.h
 * @brief A quantity damping by operator splitting schemes.
 * These methods modify the quantity directly.
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 * Note that, currently, these classes works only in single resolution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DISSIPATION_H
#define PARTICLE_DYNAMICS_DISSIPATION_H

#include "all_particle_dynamics.h"

namespace SPH
{
/* Base class to indicate the concept of operator splitting */
class OperatorSplitting
{
};

template <typename VariableType>
struct ErrorAndParameters
{
    VariableType error_;
    Real a_, c_;
    ErrorAndParameters()
        : error_(ZeroData<VariableType>::value),
          a_(0), c_(0){};
};

class FixedDampingRate
{
  public:
    FixedDampingRate(BaseParticles *particles, Real eta, Real c = 1.0)
        : damping_rate_(particles->registerSingleVariable<Real>("DampingRate", eta)),
          specific_capacity_(particles->registerSingleVariable<Real>("SpecificCapacity", c)),
          mass_(particles->getVariableDataByName<Real>("Mass")){};
    virtual ~FixedDampingRate(){};

    Real DampingRate(size_t i, size_t j) { return *damping_rate_; };
    Real DampingRate(size_t i) { return *damping_rate_; };
    Real Capacity(size_t i) { return (*specific_capacity_) * (*mass_)[i]; };

  protected:
    Real *damping_rate_;
    Real *specific_capacity_;
    StdLargeVec<Real> *mass_;
};

template <typename... InteractionTypes>
class Damping;

template <typename VariableType, typename DampingRateType, class DataDelegationType>
class Damping<Base, VariableType, DampingRateType, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public OperatorSplitting
{
  public:
    template <class BaseRelationType, typename... Args>
    explicit Damping(BaseRelationType &base_relation, const std::string &variable_name, Args &&...args);
    template <typename BodyRelationType, typename... Args, size_t... Is>
    Damping(ConstructorArgs<BodyRelationType, Args...> parameters, std::index_sequence<Is...>)
        : Damping(parameters.body_relation_, std::get<Is>(parameters.others_)...){};
    template <class BodyRelationType, typename... Args>
    explicit Damping(ConstructorArgs<BodyRelationType, Args...> parameters)
        : Damping(parameters, std::make_index_sequence<std::tuple_size_v<decltype(parameters.others_)>>{}){};
    virtual ~Damping(){};

  protected:
    typedef VariableType DampingVariable;
    std::string variable_name_;
    DampingRateType damping_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<VariableType> &variable_;
};

/**
 * @class Projection
 * @brief projection-based operator splitting method solving small system
 * around each particle using linear projection.
 * Note that here inner interaction and that with boundary are derived from the inner interaction.
 * This is because the error and parameters are computed based on both.
 */
class Projection;
template <typename VariableType, typename DampingRateType>
class Damping<Inner<Projection>, VariableType, DampingRateType>
    : public Damping<Base, VariableType, DampingRateType, DataDelegateInner>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Base, VariableType, DampingRateType, DataDelegateInner>(std::forward<Args>(args)...){};
    virtual ~Damping(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    void updateStates(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters);
};
template <typename VariableType, typename DampingRateType>
using DampingProjectionInner = Damping<Inner<Projection>, VariableType, DampingRateType>;

/**
 * @class Pairwise
 * @brief Pairwise-based operator splitting method solving for each pair of particles.
 */
class Pairwise;
template <typename VariableType, typename DampingRateType>
class Damping<Inner<Pairwise>, VariableType, DampingRateType>
    : public Damping<Base, VariableType, DampingRateType, DataDelegateInner>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Base, VariableType, DampingRateType, DataDelegateInner>(std::forward<Args>(args)...){};
    virtual ~Damping(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
template <typename VariableType, typename DampingRateType>
using DampingPairwiseInner = Damping<Inner<Pairwise>, VariableType, DampingRateType>;

template <typename VariableType, typename DampingRateType>
class Damping<Contact<Pairwise>, VariableType, DampingRateType>
    : public Damping<Base, VariableType, DampingRateType, DataDelegateContact>
{
  public:
    template <typename... Args>
    Damping(Args &&...args);
    virtual ~Damping(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<VariableType> *> contact_variable_;
};
template <typename VariableType, typename DampingRateType>
class Damping<Contact<Pairwise, Wall>, VariableType, DampingRateType>
    : public Damping<Contact<Pairwise>, VariableType, DampingRateType>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Contact<Pairwise>, VariableType, DampingRateType>(std::forward<Args>(args)...){};
    virtual ~Damping(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
template <typename VariableType, typename DampingRateType>
using DampingPairwiseContactWall = Damping<Contact<Pairwise, Wall>, VariableType, DampingRateType>;

template <typename VariableType, typename DampingRateType>
using DampingPairwiseWithWall = ComplexInteraction<Damping<Inner<Pairwise>, Contact<Pairwise, Wall>>,
                                                   VariableType, DampingRateType>;

/**
 * @class DampingWithRandomChoice
 * @brief A random choice method for obtaining static equilibrium state
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */
template <class DampingAlgorithmType>
class DampingWithRandomChoice : public DampingAlgorithmType
{
  protected:
    Real random_ratio_;
    bool RandomChoice();

  public:
    template <typename... Args>
    DampingWithRandomChoice(Real random_ratio, Args &&...args);
    virtual ~DampingWithRandomChoice(){};

    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // PARTICLE_DYNAMICS_DISSIPATION_H