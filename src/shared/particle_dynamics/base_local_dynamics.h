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
 * @file    base_local_dynamics.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_LOCAL_DYNAMICS_H
#define BASE_LOCAL_DYNAMICS_H

#include "base_data_package.h"
#include "base_particle_dynamics.h"
#include "sph_data_containers.h"

namespace SPH
{
//----------------------------------------------------------------------
// Interaction types for particle dynamics
//----------------------------------------------------------------------
class Base;            /**< Indicating base class for a method */
class Inner;           /**< Inner interaction: interaction within a body*/
class InnerAdaptive;   /**< Inner interaction with adaptive resolution */
class BaseInner;       /**< Base inner interaction */
class Contact;         /**< Contact interaction: interaction between a body with one or several another bodies */
class ContactAdaptive; /**< Contact interaction with adaptive resolution */
class BaseContact;     /**< Base contact interaction*/
class ContactBoundary; /**< Contact interaction with boundary */
class ContactWall;     /**< Contact interaction with wall boundary */
template <typename InteractionType>
class Extended; /**< An extened method of an interaction type */
//----------------------------------------------------------------------
// Particle group scope functors
//----------------------------------------------------------------------
class AllParticles
{
  public:
    explicit AllParticles(BaseParticles *base_particles){};
    bool operator()(size_t index_i)
    {
        return true;
    };
};

template <int INDICATOR>
class IndicatedParticles
{
    StdLargeVec<int> &indicator_;

  public:
    explicit IndicatedParticles(BaseParticles *base_particles)
        : indicator_(*base_particles->getVariableByName<int>("Indicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] == INDICATOR;
    };
};

using BulkParticles = IndicatedParticles<0>;

template <int INDICATOR>
class NotIndicatedParticles
{
    StdLargeVec<int> &indicator_;

  public:
    explicit NotIndicatedParticles(BaseParticles *base_particles)
        : indicator_(*base_particles->getVariableByName<int>("Indicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] != INDICATOR;
    };
};

template <typename DataType>
class PairAverageInner
{
    StdLargeVec<DataType> &variable_;

  public:
    explicit PairAverageInner(StdLargeVec<DataType> &variable)
        : variable_(variable){};
    DataType operator()(size_t index_i, size_t index_j)
    {
        return 0.5 * (variable_[index_i] + variable_[index_j]);
    };
};

template <typename DataType>
class PairAverageContact
{
    StdLargeVec<DataType> &inner_variable_;
    StdLargeVec<DataType> &contact_variable_;

  public:
    PairAverageContact(StdLargeVec<DataType> &inner_variable,
                       StdLargeVec<DataType> &contact_variable)
        : inner_variable_(inner_variable), contact_variable_(contact_variable){};
    DataType operator()(size_t index_i, size_t index_j)
    {
        return 0.5 * (inner_variable_[index_i] + contact_variable_[index_j]);
    };
};

class NoKernelCorrection
{
  public:
    NoKernelCorrection(BaseParticles *particles){};
    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};

class KernelCorrection
{
  public:
    KernelCorrection(BaseParticles *particles)
        : B_(*particles->getVariableByName<Matd>("KernelCorrectionMatrix")){};

    Matd operator()(size_t index_i)
    {
        return B_[index_i];
    };

  protected:
    StdLargeVec<Matd> &B_;
};

class SingleResolution
{
  public:
    SingleResolution(BaseParticles *particles){};
    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};

class AdaptiveResolution
{
  public:
    AdaptiveResolution(BaseParticles *particles)
        : h_ratio_(*particles->getVariableByName<Real>("SmoothingLengthRatio")){};

    Real operator()(size_t index_i)
    {
        return h_ratio_[index_i];
    };

  protected:
    StdLargeVec<Real> &h_ratio_;
};
//----------------------------------------------------------------------
// Particle reduce functors
//----------------------------------------------------------------------
template <class ReturnType>
struct ReduceSum
{
    ReturnType operator()(const ReturnType &x, const ReturnType &y) const { return x + y; };
};

struct ReduceMax
{
    Real operator()(Real x, Real y) const { return SMAX(x, y); };
};

struct ReduceMin
{
    Real operator()(Real x, Real y) const { return SMIN(x, y); };
};

struct ReduceOR
{
    bool operator()(bool x, bool y) const { return x || y; };
};

struct ReduceAND
{
    bool operator()(bool x, bool y) const { return x && y; };
};

struct ReduceLowerBound
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd lower_bound;
        for (int i = 0; i < lower_bound.size(); ++i)
            lower_bound[i] = SMIN(x[i], y[i]);
        return lower_bound;
    };
};
struct ReduceUpperBound
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd upper_bound;
        for (int i = 0; i < upper_bound.size(); ++i)
            upper_bound[i] = SMAX(x[i], y[i]);
        return upper_bound;
    };
};

/**
 * @class BaseLocalDynamics
 * @brief The base class for all local particle dynamics.
 */
template <class DynamicsIdentifier>
class BaseLocalDynamics
{
  public:
    explicit BaseLocalDynamics(DynamicsIdentifier &identifier)
        : identifier_(identifier), sph_body_(identifier.getSPHBody()){};
    virtual ~BaseLocalDynamics(){};
    SPHBody &getSPHBody() { return sph_body_; };
    DynamicsIdentifier &getDynamicsIdentifier() { return identifier_; };
    virtual void setupDynamics(Real dt = 0.0){}; // setup global parameters
  protected:
    DynamicsIdentifier &identifier_;
    SPHBody &sph_body_;
};
using LocalDynamics = BaseLocalDynamics<SPHBody>;

/**
 * @class BaseLocalDynamicsReduce
 * @brief The base class for all local particle dynamics for reducing.
 */
template <typename ReturnType, typename Operation, class DynamicsIdentifier>
class BaseLocalDynamicsReduce : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    BaseLocalDynamicsReduce(DynamicsIdentifier &identifier, ReturnType reference)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier), reference_(reference),
          quantity_name_("ReducedQuantity"){};
    virtual ~BaseLocalDynamicsReduce(){};

    using ReduceReturnType = ReturnType;
    ReturnType Reference() { return reference_; };
    std::string QuantityName() { return quantity_name_; };
    Operation &getOperation() { return operation_; };
    virtual ReturnType outputResult(ReturnType reduced_value) { return reduced_value; }

  protected:
    ReturnType reference_;
    Operation operation_;
    std::string quantity_name_;
};
template <typename ReturnType, typename Operation>
using LocalDynamicsReduce = BaseLocalDynamicsReduce<ReturnType, Operation, SPHBody>;

/**
 * @class Average
 * @brief Derives class for computing particle-wise averages.
 */
template <class ReduceSumType>
class Average : public ReduceSumType
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    Average(DynamicsIdentifier &identifier, Args &&...args)
        : ReduceSumType(identifier, std::forward<Args>(args)...){};
    virtual ~Average(){};
    using ReturnType = typename ReduceSumType::ReduceReturnType;

    virtual ReturnType outputResult(ReturnType reduced_value)
    {
        ReturnType sum = ReduceSumType::outputResult(reduced_value);
        return sum / Real(this->getDynamicsIdentifier().SizeOfLoopRange());
    }
};

/**
 * @class ConstructorArgs
 * @brief Class template argument deduction (CTAD) for constructor arguments.
 */
template <typename BodyRelationType, typename... OtherArgs>
struct ConstructorArgs
{
    BodyRelationType &body_relation_;
    std::tuple<OtherArgs...> others_;
    ConstructorArgs(BodyRelationType &body_relation, OtherArgs... other_args)
        : body_relation_(body_relation), others_(other_args...){};
};

/**
 * @class ComplexInteraction
 * @brief A class that integrates multiple local dynamics.
 * Typically, it includes an inner interaction and one or
 * several contact interaction and boundary conditions.
 */
template <typename... T>
class ComplexInteraction;

template <typename... ControlTypes, template <typename... InteractionTypes> class LocalDynamicsName>
class ComplexInteraction<LocalDynamicsName<>, ControlTypes...>
{
  public:
    ComplexInteraction(){};

    void interaction(size_t index_i, Real dt = 0.0){};
};

template <typename... ControlTypes, template <typename... InteractionTypes> class LocalDynamicsName,
          class FirstInteraction, class... OtherInteractions>
class ComplexInteraction<LocalDynamicsName<FirstInteraction, OtherInteractions...>, ControlTypes...>
    : public LocalDynamicsName<FirstInteraction, ControlTypes...>
{
  protected:
    ComplexInteraction<LocalDynamicsName<OtherInteractions...>, ControlTypes...> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit ComplexInteraction(FirstParameterSet &&first_parameter_set,
                                OtherParameterSets &&...other_parameter_sets)
        : LocalDynamicsName<FirstInteraction, ControlTypes...>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        LocalDynamicsName<FirstInteraction, ControlTypes...>::interaction(index_i, dt);
        other_interactions_.interaction(index_i, dt);
    };
};

/**
 * @class OldComplexInteraction
 * @brief A class that integrates multiple local dynamics.
 * Typically, it includes an inner interaction and one or
 * several contact interaction ad boundary conditions.
 */
template <typename... InteractionType>
class OldComplexInteraction;

template <>
class OldComplexInteraction<>
{
  public:
    OldComplexInteraction(){};

    void interaction(size_t index_i, Real dt = 0.0){};
};

template <class FirstInteraction, class... OtherInteractions>
class OldComplexInteraction<FirstInteraction, OtherInteractions...> : public FirstInteraction
{
  protected:
    OldComplexInteraction<OtherInteractions...> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit OldComplexInteraction(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : FirstInteraction(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        FirstInteraction::interaction(index_i, dt);
        other_interactions_.interaction(index_i, dt);
    };
};
} // namespace SPH
#endif // BASE_LOCAL_DYNAMICS_H
