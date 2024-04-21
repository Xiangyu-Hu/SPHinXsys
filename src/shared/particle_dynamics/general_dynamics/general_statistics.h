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
 * @file general_statistics.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef GENERAL_STATISTICS_H
#define GENERAL_STATISTICS_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class DataType>
struct SampleAdd
{
    using Sample = std::<DataType, size_t>;
    Sample reference;
    SampleAdd() : reference(ZeroData<ReturnType>::value, 0){};
    Sample operator()(const Sample &x, const Sample &y) const { return x + y; };
};

/**
 * @class SampleSummation
 * @brief Compute the summation of a discrete variable and the total number samples
 */
template <typename DataType, class DynamicsIdentifier>
class SampleSummation : public BaseLocalDynamicsReduce<SampleAdd<DataType>, DynamicsIdentifier>,
                        public GeneralDataDelegateSimple
{
  protected:
    using Sample = std::<DataType, size_t>;
    StdLargeVec<DataType> &variable_;

  public:
    using SampleDataType = DataType;

    explicit SampleSummation(SPHBody &sph_body, const std::string &variable_name)
        : BaseLocalDynamicsReduce<SampleAdd<DataType>, DynamicsIdentifier>(sph_body),
          GeneralDataDelegateSimple(sph_body),
          variable_(*this->particles_->template getVariableByName<DataType>(variable_name))
    {
        this->quantity_name_ = "Total" + variable_name + "AndSamples";
    };
    virtual ~QuantitySummation(){};

    Sample reduce(size_t index_i, Real dt = 0.0)
    {
        return Sample(variable_[index_i], 1);
    };
};

template <typename SampleSummationType, class ExecutionPolicy = ParallelPolicy>
class SampleMean : public BaseDynamics<typename SampleSummationType::SampleDataType>
{
  protected:
    using DataType = typename SampleSummationType::SampleDataType;
    ReduceDynamics<SampleSummationType, ExecutionPolicy> sample_summation_;

  public:
    template <class DynamicsIdentifier, typename... Args>
    explicit SampleMean(DynamicsIdentifier &identifier, Args &&...args)
        : BaseDynamics<DataType>(identifier.getSPHBody()),
          sample_summation_(identifier, std::forward<Args>(args)...){};

    virtual DataType exec(Real dt = 0.0) override
    {
        auto total = sample_summation_.exec();
        return total.first / Real(total.second);
    };
};
} // namespace SPH
#endif // GENERAL_STATISTICS_H
