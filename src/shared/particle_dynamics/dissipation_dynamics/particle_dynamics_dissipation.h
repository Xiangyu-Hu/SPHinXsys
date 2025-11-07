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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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

template <typename DataType>
struct ErrorAndParameters
{
    DataType error_;
    Real a_, c_;
    ErrorAndParameters()
        : error_(ZeroData<DataType>::value),
          a_(0), c_(0) {};
};

class FixedDampingRate
{
  public:
    FixedDampingRate(BaseParticles *particles, Real eta, Real c = 1.0)
        : damping_rate_(particles->registerSingularVariable<Real>("DampingRate", eta)->Data()),
          specific_capacity_(particles->registerSingularVariable<Real>("SpecificCapacity", c)->Data()),
          mass_(particles->getVariableDataByName<Real>("Mass")) {};
    virtual ~FixedDampingRate() {};

    Real DampingRate(size_t i, size_t j) { return *damping_rate_; };
    Real DampingRate(size_t i) { return *damping_rate_; };
    Real Capacity(size_t i) { return (*specific_capacity_) * mass_[i]; };

  protected:
    Real *damping_rate_;
    Real *specific_capacity_;
    Real *mass_;
};

template <typename... InteractionTypes>
class Damping;

template <typename DataType, typename DampingRateType, class DataDelegationType>
class Damping<Base, DataType, DampingRateType, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType, typename... Args>
    explicit Damping(BaseRelationType &base_relation, const std::string &name, Args &&...args);
    template <typename BodyRelationType, typename... Args, size_t... Is>
    Damping(DynamicsArgs<BodyRelationType, Args...> parameters, std::index_sequence<Is...>)
        : Damping(parameters.identifier_, std::get<Is>(parameters.others_)...){};
    template <class BodyRelationType, typename... Args>
    explicit Damping(DynamicsArgs<BodyRelationType, Args...> parameters)
        : Damping(parameters, std::make_index_sequence<std::tuple_size_v<decltype(parameters.others_)>>{}){};
    virtual ~Damping() {};

  protected:
    typedef DataType DampingVariable;
    std::string name_;
    DampingRateType damping_;
    Real *Vol_;
    DataType *data_field_;
};

/**
 * @class Projection
 * @brief projection-based operator splitting method solving small system
 * around each particle using linear projection.
 * Note that here inner interaction and that with boundary are derived from the inner interaction.
 * This is because the error and parameters are computed based on both.
 */
class Projection;
template <typename DataType, typename DampingRateType>
class Damping<Inner<Projection>, DataType, DampingRateType>
    : public Damping<Base, DataType, DampingRateType, DataDelegateInner>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Base, DataType, DampingRateType, DataDelegateInner>(std::forward<Args>(args)...){};
    virtual ~Damping() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ErrorAndParameters<DataType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    void updateStates(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters);
};
template <typename DataType, typename DampingRateType>
using DampingProjectionInner = Damping<Inner<Projection>, DataType, DampingRateType>;

/**
 * @class Pairwise
 * @brief Pairwise-based operator splitting method solving for each pair of particles.
 */
class Pairwise;
template <typename DataType, typename DampingRateType>
class Damping<Inner<Pairwise>, DataType, DampingRateType>
    : public Damping<Base, DataType, DampingRateType, DataDelegateInner>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Base, DataType, DampingRateType, DataDelegateInner>(std::forward<Args>(args)...){};
    virtual ~Damping() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
template <typename DataType, typename DampingRateType>
using DampingPairwiseInner = Damping<Inner<Pairwise>, DataType, DampingRateType>;

template <typename DataType, typename DampingRateType>
class Damping<Contact<Pairwise>, DataType, DampingRateType>
    : public Damping<Base, DataType, DampingRateType, DataDelegateContact>
{
  public:
    template <typename... Args>
    Damping(Args &&...args);
    virtual ~Damping() {};

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<DataType *> contact_data_field_;
};
template <typename DataType, typename DampingRateType>
class Damping<Contact<Pairwise, Wall>, DataType, DampingRateType>
    : public Damping<Contact<Pairwise>, DataType, DampingRateType>
{
  public:
    template <typename... Args>
    Damping(Args &&...args)
        : Damping<Contact<Pairwise>, DataType, DampingRateType>(std::forward<Args>(args)...){};
    virtual ~Damping() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
template <typename DataType, typename DampingRateType>
using DampingPairwiseContactWall = Damping<Contact<Pairwise, Wall>, DataType, DampingRateType>;

template <typename DataType, typename DampingRateType>
using DampingPairwiseWithWall = ComplexInteraction<Damping<Inner<Pairwise>, Contact<Pairwise, Wall>>,
                                                   DataType, DampingRateType>;

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
    virtual ~DampingWithRandomChoice() {};

    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // PARTICLE_DYNAMICS_DISSIPATION_H