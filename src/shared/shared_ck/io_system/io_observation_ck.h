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
 * @file 	io_observation_ck.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef IO_OBSERVATION_CK_H
#define IO_OBSERVATION_CK_H

#include "io_observation.h"

#include "execution_policy.h"
#include "interpolation_dynamics.hpp"
#include "subdomain_fan_out.h"

namespace SPH
{
template <class ExecutionPolicy, typename DataType, typename... Parameters>
class ObservedQuantityRecording<ExecutionPolicy, DataType, Parameters...>
    : public BaseQuantityRecording
{
  protected:
    SPHBody &observer_;
    BaseParticles &base_particles_;
    BaseParticles &contact_particles_;
    ObservingQuantityCK<ExecutionPolicy, DataType, Parameters...> observation_method_;
    DiscreteVariable<DataType> *dv_interpolated_quantities_;
    size_t number_of_observe_;

  public:
    DataType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    template <typename... RelationParameters, typename... Args>
    ObservedQuantityRecording(Contact<RelationParameters...> &contact_relation, Args &&...args)
        : BaseQuantityRecording(
              contact_relation.getSPHBody().getSPHSystem(), contact_relation.getSPHBody().Name()),
          observer_(contact_relation.getSPHBody()),
          base_particles_(observer_.getBaseParticles()),
          contact_particles_(contact_relation.getContactParticles()),
          observation_method_(contact_relation, std::forward<Args>(args)...),
          dv_interpolated_quantities_(observation_method_.dvInterpolatedQuantities()),
          number_of_observe_(base_particles_.TotalRealParticles())
    {
        setFullPath(dv_interpolated_quantities_->Name());
    };
    virtual ~ObservedQuantityRecording() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        if (!header_written_)
        {
            std::ofstream out_file(filefullpath_output_.c_str(), std::ios::out);
            out_file << "run_time" << "   ";
            for (size_t i = 0; i != number_of_observe_; ++i)
            {
                std::string quantity_name_i = quantity_name_ + "[" + std::to_string(i) + "]";
                plt_engine_.writeAQuantityHeader(
                    out_file, dv_interpolated_quantities_->getValueWithScalingRef(i), quantity_name_i);
            }
            out_file << "\n";
            out_file.close();
            header_written_ = true;
        }
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_->getValueWithScalingRef() << "   ";
        observation_method_.exec();
        collectInterpolatedQuantities(ExecutionPolicy{});
        for (size_t i = 0; i != number_of_observe_; ++i)
        {
            plt_engine_.writeAQuantity(
                out_file, dv_interpolated_quantities_->getValueWithScalingRef(i));
        }
        out_file << "\n";
        out_file.close();
    };

    DataType *getObservedQuantity()
    {
        return this->dv_interpolated_quantities_->Data();
    };

  protected:
    template <class Policy>
    void collectInterpolatedQuantities(const Policy &ex_policy)
    {
        dv_interpolated_quantities_->prepareForOutput(ex_policy);
    };
    /** Decomposed run of the observed body. The observer is replicated: every subdomain
     *  interpolated every observer particle from its own owned and halo particles, which
     *  is complete only for the observer particles inside its slab. Each value is
     *  therefore taken from the subdomain owning the observation point. */
    template <class PolicyType>
    void collectInterpolatedQuantities(const DecomposedExecution<PolicyType> &ex_policy)
    {
        SubdomainExchangeInterface *exchange = contact_particles_.getSubdomainExchange();
        if (exchange == nullptr)
        { // observed body not decomposed, hence replicated as a whole
            dv_interpolated_quantities_->prepareForOutput(PolicyType{});
            return;
        }
        const Vecd *position = base_particles_.dvParticlePosition()->Data();
        const UnsignedInt width = dv_interpolated_quantities_->getWidth();
        DataType *host_data = dv_interpolated_quantities_->Data();
        StdVec<DataType> replica_copy(number_of_observe_ * width);
        for (int subdomain_id = 0; subdomain_id < execution::numberOfSubdomains(); ++subdomain_id)
        {
            execution::SubdomainScope scope(subdomain_id);
            const DataType *replica = dv_interpolated_quantities_->DelegatedData(ex_policy);
            execution::copyBetweenSubdomains(ex_policy, subdomain_id, replica_copy.data(), replica, replica_copy.size());
            for (size_t i = 0; i != number_of_observe_; ++i)
            {
                if (exchange->subdomainOf(position[i]) == subdomain_id)
                {
                    for (UnsignedInt entry = 0; entry < width; ++entry)
                        host_data[i * width + entry] = replica_copy[i * width + entry];
                }
            }
        }
    };

  public:

    size_t NumberOfObservedQuantity()
    {
        return number_of_observe_;
    };

    DiscreteVariable<DataType> &getObservedVariable()
    {
        return *dv_interpolated_quantities_;
    };
};

template <class ExecutionPolicy, class LocalReduceMethodType>
class ReducedQuantityRecording<ExecutionPolicy, LocalReduceMethodType> : public BaseQuantityRecording
{
  protected:
    ReduceDynamicsCK<ExecutionPolicy, LocalReduceMethodType> reduce_method_;

  public:
    /*< deduce variable type from reduce method. */
    using VariableType = typename LocalReduceMethodType::FinishDynamics::OutputType;
    VariableType type_indicator_; /*< this is an indicator to identify the variable type. */
    VariableType reduced_quantity_;

  public:
    template <class DynamicsIdentifier, typename... Args>
    ReducedQuantityRecording(DynamicsIdentifier &identifier, Args &&...args)
        : BaseQuantityRecording(identifier.getSPHBody().getSPHSystem(),
                                identifier.Name()),
          reduce_method_(identifier, std::forward<Args>(args)...),
          reduced_quantity_(ZeroData<VariableType>::value)
    {
        quantity_name_ = reduce_method_.QuantityName();
        setFullPath(quantity_name_);
    };
    virtual ~ReducedQuantityRecording() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        if (!header_written_)
        {
            std::ofstream out_file(filefullpath_output_.c_str(), std::ios::out);
            out_file << "\"run_time\"" << "   ";
            plt_engine_.writeAQuantityHeader(out_file, reduced_quantity_, quantity_name_);
            out_file << "\n";
            out_file.close();
            header_written_ = true;
        }
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_->getValue() << "   ";
        reduced_quantity_ = reduce_method_.exec();
        plt_engine_.writeAQuantity(out_file, reduced_quantity_);
        out_file << "\n";
        out_file.close();
        header_written_ = true;
    };

    VariableType *getObservedQuantity()
    {
        return &reduced_quantity_;
    };

    size_t NumberOfObservedQuantity()
    {
        return 1;
    };
};
} // namespace SPH
#endif // IO_OBSERVATION_CK_H
