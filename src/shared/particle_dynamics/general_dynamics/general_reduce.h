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
 * @brief  obtained the maximum norm of a input
 */
template <typename InputFunction>
class MaximumNorm : public LocalDynamicsReduce<Real, ReduceMax>,
                    public GeneralDataDelegateSimple
{
  public:
    template <typename... Args>
    MaximumNorm(SPHBody &sph_body, Args &&...args)
        : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
          GeneralDataDelegateSimple(sph_body),
          input_function_(this->getParticles(), std::forward<Args>(args)...){};
    virtual ~MaximumNorm(){};
    virtual Real outputResult(Real reduced_value) override { return std::sqrt(reduced_value); }
    Real reduce(size_t index_i, Real dt = 0.0) { return getSquaredNorm(input_function_(index_i)); };

  protected:
    InputFunction input_function_;

    template <typename Datatype>
    Real getSquaredNorm(const Datatype &input) { return input.squaredNorm(); };

    Real getSquaredNorm(const Real &input) { return input * input; };
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
        : BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>(identifier, Real(MinRealNumber)),
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
 * @brief Compute the summation of  a particle input in a body
 */
template <typename inputType>
class QuantitySummation : public LocalDynamicsReduce<inputType, ReduceSum<inputType>>,
                          public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<inputType> &input_;

  public:
    explicit QuantitySummation(SPHBody &sph_body, const std::string &input_name)
        : LocalDynamicsReduce<inputType, ReduceSum<inputType>>(sph_body, ZeroData<inputType>::value),
          GeneralDataDelegateSimple(sph_body),
          input_(*this->particles_->template getinputByName<inputType>(input_name))
    {
        this->quantity_name_ = input_name + "Summation";
    };
    virtual ~QuantitySummation(){};

    inputType reduce(size_t index_i, Real dt = 0.0)
    {
        return input_[index_i];
    };
};

/**
 * @class QuantityMoment
 * @brief Compute the moment of a body
 */
template <typename inputType>
class QuantityMoment : public QuantitySummation<inputType>
{
  protected:
    StdLargeVec<Real> &mass_;

  public:
    explicit QuantityMoment(SPHBody &sph_body, const std::string &input_name)
        : QuantitySummation<inputType>(sph_body, input_name),
          mass_(this->particles_->mass_)
    {
        this->quantity_name_ = input_name + "Moment";
    };
    virtual ~QuantityMoment(){};

    inputType reduce(size_t index_i, Real dt = 0.0)
    {
        return mass_[index_i] * this->input_[index_i];
    };
};

/**
 * @class TotalMechanicalEnergy
 * @brief Compute the total mechanical (kinematic and potential) energy
 */
class TotalMechanicalEnergy
    : public LocalDynamicsReduce<Real, ReduceSum<Real>>,
      public GeneralDataDelegateSimple
{
  private:
    SharedPtrKeeper<Gravity> gravity_ptr_keeper_;

  protected:
    StdLargeVec<Real> &mass_;
    StdLargeVec<Vecd> &vel_, &pos_;
    Gravity *gravity_;

  public:
    explicit TotalMechanicalEnergy(SPHBody &sph_body, SharedPtr<Gravity> = makeShared<Gravity>(Vecd::Zero()));
    virtual ~TotalMechanicalEnergy(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_REDUCE_H
