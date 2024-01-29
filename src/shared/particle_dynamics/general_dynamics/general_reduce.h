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
 * @file 	general_reduce.h
 * @brief TBD
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_REDUCE_H
#define GENERAL_REDUCE_H

#include "base_general_dynamics.h"
#include <limits>

namespace SPH
{
/**
 * @class MaximumNorm
 * @brief  obtained the maximum norm of a variable
 */
template <typename DataType>
class MaximumNorm : public LocalDynamicsReduce<Real, ReduceMax>,
                    public GeneralDataDelegateSimple
{
  public:
    MaximumNorm(SPHBody &sph_body, const std::string &variable_name)
        : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
          GeneralDataDelegateSimple(sph_body),
          variable_(*particles_->getVariableByName<DataType>(variable_name)){};
    virtual ~MaximumNorm(){};
    virtual Real outputResult(Real reduced_value) override { return std::sqrt(reduced_value); }
    Real reduce(size_t index_i, Real dt = 0.0) { return getSquaredNorm(variable_[index_i]); };

  protected:
    StdLargeVec<DataType> &variable_;

    template <typename Datatype>
    Real getSquaredNorm(const Datatype &variable) { return variable.squaredNorm(); };

    Real getSquaredNorm(const Real &variable) { return variable * variable; };
};

/**
 * @class VelocityBoundCheck
 * @brief  check whether particle velocity within a given bound
 */
class VelocityBoundCheck : public LocalDynamicsReduce<bool, ReduceOR>,
                           public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &vel_;
    Real velocity_bound_;

  public:
    VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound);
    virtual ~VelocityBoundCheck(){};

    bool reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class 	UpperFrontInAxisDirection
 * @brief 	Get the upper front in an axis direction for a body or body part
 */
template <class DynamicsIdentifier>
class UpperFrontInAxisDirection : public BaseLocalDynamicsReduce<Real, ReduceMax, DynamicsIdentifier>,
                                  public GeneralDataDelegateSimple
{
  protected:
    int axis_;
    StdLargeVec<Vecd> &pos_;

  public:
    explicit UpperFrontInAxisDirection(DynamicsIdentifier &identifier, std::string name, int axis = lastAxis)
        : BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>(identifier, MinReal),
          GeneralDataDelegateSimple(identifier.getSPHBody()), axis_(axis), pos_(particles_->pos_)
    {
        this->quantity_name_ = name;
    }
    virtual ~UpperFrontInAxisDirection(){};

    Real reduce(size_t index_i, Real dt = 0.0) { return pos_[index_i][axis_]; };
};

/**
 * @class MaximumSpeed
 * @brief Get the maximum particle speed in a SPH body
 */
class MaximumSpeed : public LocalDynamicsReduce<Real, ReduceMax>,
                     public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &vel_;

  public:
    explicit MaximumSpeed(SPHBody &sph_body);
    virtual ~MaximumSpeed(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class	PositionLowerBound
 * @brief	the lower bound of a body by reduced particle positions.
 * 			TODO: a test using this method
 */
class PositionLowerBound : public LocalDynamicsReduce<Vecd, ReduceLowerBound>,
                           public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;

  public:
    explicit PositionLowerBound(SPHBody &sph_body);
    virtual ~PositionLowerBound(){};

    Vecd reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class	PositionUpperBound
 * @brief	the upper bound of a body by reduced particle positions.
 * 			TODO: a test using this method
 */
class PositionUpperBound : public LocalDynamicsReduce<Vecd, ReduceUpperBound>,
                           public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;

  public:
    explicit PositionUpperBound(SPHBody &sph_body);
    virtual ~PositionUpperBound(){};

    Vecd reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class QuantitySummation
 * @brief Compute the summation of  a particle variable in a body
 */
template <typename VariableType>
class QuantitySummation : public LocalDynamicsReduce<VariableType, ReduceSum<VariableType>>,
                          public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<VariableType> &variable_;

  public:
    explicit QuantitySummation(SPHBody &sph_body, const std::string &variable_name)
        : LocalDynamicsReduce<VariableType, ReduceSum<VariableType>>(sph_body, ZeroData<VariableType>::value),
          GeneralDataDelegateSimple(sph_body),
          variable_(*this->particles_->template getVariableByName<VariableType>(variable_name))
    {
        this->quantity_name_ = "Total" + variable_name;
    };
    virtual ~QuantitySummation(){};

    VariableType reduce(size_t index_i, Real dt = 0.0)
    {
        return variable_[index_i];
    };
};

/**
 * @class QuantityMoment
 * @brief Compute the moment of a body
 */
template <typename VariableType>
class QuantityMoment : public QuantitySummation<VariableType>
{
  protected:
    StdLargeVec<Real> &mass_;

  public:
    explicit QuantityMoment(SPHBody &sph_body, const std::string &variable_name)
        : QuantitySummation<VariableType>(sph_body, variable_name),
          mass_(this->particles_->mass_)
    {
        this->quantity_name_ = variable_name + "Moment";
    };
    virtual ~QuantityMoment(){};

    VariableType reduce(size_t index_i, Real dt = 0.0)
    {
        return mass_[index_i] * this->variable_[index_i];
    };
};

class TotalKineticEnergy
    : public LocalDynamicsReduce<Real, ReduceSum<Real>>,
      public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Real> &mass_;
    StdLargeVec<Vecd> &vel_;

  public:
    explicit TotalKineticEnergy(SPHBody &sph_body);
    virtual ~TotalKineticEnergy(){};
    Real reduce(size_t index_i, Real dt = 0.0);
};

class TotalMechanicalEnergy : public TotalKineticEnergy
{
  protected:
    Gravity &gravity_;
    StdLargeVec<Vecd> &pos_;

  public:
    explicit TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity);
    virtual ~TotalMechanicalEnergy(){};
    Real reduce(size_t index_i, Real dt = 0.0);
};

} // namespace SPH
#endif // GENERAL_REDUCE_H
