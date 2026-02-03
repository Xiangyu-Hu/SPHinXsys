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
 * @file    base_local_dynamics.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_LOCAL_DYNAMICS_H
#define BASE_LOCAL_DYNAMICS_H

#include "base_data_type_package.h"
#include "base_particle_dynamics.h"
#include "execution_policy.h"
#include "io_log.h"
#include "reduce_functors.h"
#include "sphinxsys_containers.h"

#include <type_traits>

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
  public:
    explicit BaseLocalDynamics(DynamicsIdentifier &identifier)
        : identifier_(&identifier), sph_system_(&identifier.getSPHSystem()),
          sph_body_(&identifier.getSPHBody()),
          sph_adaptation_(&sph_body_->getSPHAdaptation()),
          particles_(&sph_body_->getBaseParticles()),
          logger_(Log::get()) {};
    virtual ~BaseLocalDynamics() {};
    using Identifier = typename DynamicsIdentifier::BaseIdentifier;
    SPHBody &getSPHBody() { return *sph_body_; };
    BaseParticles &getBaseParticles() { return *particles_; };
    SPHAdaptation &getSPHAdaptation() { return *sph_adaptation_; };
    virtual void setupDynamics(Real dt = 0.0) {}; // setup global parameters

    class FinishDynamics
    {
      public:
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser){};
        void operator()() {};
    };

  protected:
    DynamicsIdentifier *identifier_;
    SPHSystem *sph_system_;
    SPHBody *sph_body_;
    SPHAdaptation *sph_adaptation_;
    BaseParticles *particles_;
    std::shared_ptr<spdlog::logger> logger_;
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
    typedef Operation OperationType;
    using ReturnType = typename Operation::ReturnType;
    explicit BaseLocalDynamicsReduce(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          reference_(ReduceReference<Operation>::value),
          quantity_name_("ReducedQuantity") {};
    virtual ~BaseLocalDynamicsReduce() {};

    ReturnType Reference() { return reference_; };
    std::string QuantityName() { return quantity_name_; };
    Operation &getOperation() { return operation_; };
    virtual ReturnType outputResult(ReturnType reduced_value) { return reduced_value; }

    class FinishDynamics
    {
      public:
        using OutputType = ReturnType;
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser){};
        ReturnType Result(ReturnType reduced_value)
        {
            return reduced_value;
        }
    };

  protected:
    Operation operation_;
    ReturnType reference_;
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
    virtual ~Average() {};
    using ReturnType = typename ReduceSumType::ReturnType;

    virtual ReturnType outputResult(ReturnType reduced_value)
    {
        ReturnType sum = ReduceSumType::outputResult(reduced_value);
        return sum / Real(this->identifier_->SizeOfLoopRange());
    }
};

/**
 * @class DynamicsArgs
 * @brief Class template argument deduction (CTAD) for constructing interaction dynamics.
 * @details Note that the form "XXX" is not std::string type, so we need to use
 * std::string("XXX") to convert it to std::string type.
 * Only the DynamicsIdentifier parameter is reference,
 * the other parameters should not use it, use pointer
 * instead.
 */
template <typename DynamicsIdentifier, typename... OtherArgs>
struct DynamicsArgs
{
    DynamicsIdentifier &identifier_;
    std::tuple<OtherArgs...> others_;
    SPHBody &getSPHBody() { return identifier_.getSPHBody(); };
    DynamicsArgs(DynamicsIdentifier &identifier, OtherArgs... other_args)
        : identifier_(identifier), others_(other_args...) {};
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
    ComplexInteraction() {};

    void interaction(size_t index_i, Real dt = 0.0) {};
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
