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
#include "execution_policy.h"
#include "sphinxsys_containers.h"

namespace SPH
{
/**
 * @class BaseLocalDynamics
 * @brief The base class for all local particle dynamics.
 * @details The basic design idea is define local dynamics for local particle operations.
 * We split a general local dynamics into two parts in respect of functionality:
 * one is the action on singular data, which is carried within the function setupDynamics,
 * the other is the action on discrete variables, which will be carried out
 * by the computing kernel. In the scenarios of offloading computing,
 * the first function is generally carried on the host and the other on computing device.
 * We also split the local dynamics into two part in respect of memory management.
 */
template <class DynamicsIdentifier>
class BaseLocalDynamics
{
    UniquePtrsKeeper<BaseVariable> constant_entity_ptrs_;

  public:
    explicit BaseLocalDynamics(DynamicsIdentifier &identifier)
        : identifier_(identifier), sph_system_(identifier.getSPHSystem()),
          sph_body_(identifier.getSPHBody()),
          particles_(&sph_body_.getBaseParticles()){};
    virtual ~BaseLocalDynamics(){};
    DynamicsIdentifier &getDynamicsIdentifier() { return identifier_; };
    SPHBody &getSPHBody() { return sph_body_; };
    BaseParticles *getParticles() { return particles_; };
    virtual void setupDynamics(Real dt = 0.0){}; // setup global parameters
    void registerComputingKernel(Implementation<Base> *implementation)
    {
        sph_body_.registerComputingKernel(implementation);
    };

  protected:
    DynamicsIdentifier &identifier_;
    SPHSystem &sph_system_;
    SPHBody &sph_body_;
    BaseParticles *particles_;
};
using LocalDynamics = BaseLocalDynamics<SPHBody>;

/**
 * @class BaseLocalDynamicsReduce
 * @brief The base class for all local particle dynamics for reducing.
 */
template <typename Operation, class DynamicsIdentifier>
class BaseLocalDynamicsReduce : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    explicit BaseLocalDynamicsReduce(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          quantity_name_("ReducedQuantity"){};
    virtual ~BaseLocalDynamicsReduce(){};

    using ReturnType = decltype(Operation::reference_);
    ReturnType Reference() { return operation_.reference_; };
    std::string QuantityName() { return quantity_name_; };
    Operation &getOperation() { return operation_; };
    virtual ReturnType outputResult(ReturnType reduced_value) { return reduced_value; }

  protected:
    Operation operation_;
    std::string quantity_name_;
};
template <typename Operation>
using LocalDynamicsReduce = BaseLocalDynamicsReduce<Operation, SPHBody>;

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
    using ReturnType = typename ReduceSumType::ReturnType;

    virtual ReturnType outputResult(ReturnType reduced_value)
    {
        ReturnType sum = ReduceSumType::outputResult(reduced_value);
        return sum / Real(this->getDynamicsIdentifier().SizeOfLoopRange());
    }
};

/**
 * @class ConstructorArgs
 * @brief Class template argument deduction (CTAD) for constructor arguments.
 * @details Note that the form "XXX" is not std::string type, so we need to use
 * std::string("XXX") to convert it to std::string type.
 */
template <typename BodyRelationType, typename... OtherArgs>
struct ConstructorArgs
{
    BodyRelationType &body_relation_;
    std::tuple<OtherArgs...> others_;
    SPHBody &getSPHBody() { return body_relation_.getSPHBody(); };
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

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName>
class ComplexInteraction<LocalDynamicsName<>, CommonParameters...>
{
  public:
    ComplexInteraction(){};

    void interaction(size_t index_i, Real dt = 0.0){};
};

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName,
          class FirstInteraction, class... OtherInteractions>
class ComplexInteraction<LocalDynamicsName<FirstInteraction, OtherInteractions...>, CommonParameters...>
    : public LocalDynamicsName<FirstInteraction, CommonParameters...>
{
  protected:
    ComplexInteraction<LocalDynamicsName<OtherInteractions...>, CommonParameters...> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit ComplexInteraction(FirstParameterSet &&first_parameter_set,
                                OtherParameterSets &&...other_parameter_sets)
        : LocalDynamicsName<FirstInteraction, CommonParameters...>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        LocalDynamicsName<FirstInteraction, CommonParameters...>::interaction(index_i, dt);
        other_interactions_.interaction(index_i, dt);
    };
};
} // namespace SPH
#endif // BASE_LOCAL_DYNAMICS_H
