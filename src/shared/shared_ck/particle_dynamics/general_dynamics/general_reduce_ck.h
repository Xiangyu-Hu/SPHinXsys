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
 * @file 	general_reduce_ck.h
 * @brief TBD
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_REDUCE_CK_H
#define GENERAL_REDUCE_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
class TotalKineticEnergyCK
    : public LocalDynamicsReduce<ReduceSum<Real>>
{

  public:
    explicit TotalKineticEnergyCK(SPHBody &sph_body);
    virtual ~TotalKineticEnergyCK() {};

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Real reduce(size_t index_i, Real dt = 0.0)
        {
            return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
        };

      protected:
        Real *mass_;
        Vecd *vel_;
    };

  protected:
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Vecd> *dv_vel_;
};

class TotalMechanicalEnergyCK : public TotalKineticEnergyCK
{

  public:
    explicit TotalMechanicalEnergyCK(SPHBody &sph_body, Gravity &gravity);
    virtual ~TotalMechanicalEnergyCK() {};

    class ReduceKernel : public TotalKineticEnergyCK::ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            return TotalKineticEnergyCK::ReduceKernel::reduce(index_i, dt) +
                   mass_[index_i] * gravity_.getPotential(pos_[index_i]);
        };

      protected:
        Gravity gravity_;
        Vecd *pos_;
    };

  protected:
    Gravity &gravity_;
    DiscreteVariable<Vecd> *dv_pos_;
};

template <typename DataType, class DynamicsIdentifier = SPHBody>
class QuantitySum : public BaseLocalDynamicsReduce<ReduceSum<DataType>, DynamicsIdentifier>
{
  public:
    QuantitySum(DynamicsIdentifier &identifier, const std::string &variable_name);
    virtual ~QuantitySum() {};

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        DataType reduce(size_t index_i, Real dt = 0.0)
        {
            return variable_[index_i];
        };

      protected:
        DataType *variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
};

template <typename DataType, class DynamicsIdentifier = SPHBody>
class QuantityAverage : public BaseLocalDynamicsReduce<ReduceSum<std::pair<DataType, Real>>, DynamicsIdentifier>
{
    using ReduceReturnType = std::pair<DataType, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, DynamicsIdentifier>;

  public:
    QuantityAverage(DynamicsIdentifier &identifier, const std::string &variable_name);
    virtual ~QuantityAverage() {};

    class FinishDynamics
    {
      public:
        using OutputType = DataType;
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser){};
        OutputType Result(const ReduceReturnType &reduced_value)
        {
            return reduced_value.first / reduced_value.second;
        }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            return ReduceReturnType(variable_[index_i], Real(1));
        };

      protected:
        DataType *variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
};

template <class DynamicsIdentifier>
class UpperFrontInAxisDirectionCK : public BaseLocalDynamicsReduce<ReduceMax, DynamicsIdentifier>
{
  public:
    UpperFrontInAxisDirectionCK(DynamicsIdentifier &identifier, const std::string &name, int axis = lastAxis);
    virtual ~UpperFrontInAxisDirectionCK() {};

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            return pos_[index_i][axis_];
        };

      protected:
        int axis_;
        Vecd *pos_;
    };

  protected:
    int axis_;
    DiscreteVariable<Vecd> *dv_pos_;
};

} // namespace SPH
#endif // GENERAL_REDUCE_CK_H
