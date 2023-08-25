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
 * @file 	general_dynamics.h
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_DYNAMICS_H
#define GENERAL_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.h"
#include "external_force.h"

#include <limits>

namespace SPH
{
typedef DataDelegateSimple<BaseParticles> GeneralDataDelegateSimple;
typedef DataDelegateInner<BaseParticles> GeneralDataDelegateInner;
typedef DataDelegateContact<BaseParticles, BaseParticles> GeneralDataDelegateContact;
typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
    GeneralDataDelegateContactOnly;

/**
 * @class BaseTimeStepInitialization
 * @brief base class for time step initialization.
 */
class BaseTimeStepInitialization : public LocalDynamics
{
  private:
    SharedPtrKeeper<Gravity> gravity_ptr_keeper_;

  protected:
    Gravity *gravity_;

  public:
    BaseTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> &gravity_ptr)
        : LocalDynamics(sph_body), gravity_(gravity_ptr_keeper_.assignPtr(gravity_ptr)){};
    virtual ~BaseTimeStepInitialization(){};
};

/**
 * @class TimeStepInitialization
 * @brief initialize a time step for a body.
 */
class TimeStepInitialization
    : public BaseTimeStepInitialization,
      public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_, &acc_prior_;

  public:
    TimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
    virtual ~TimeStepInitialization(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class RandomizeParticlePosition
 * @brief Randomize the initial particle position
 */
class RandomizeParticlePosition
    : public LocalDynamics,
      public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;
    Real randomize_scale_;

  public:
    explicit RandomizeParticlePosition(SPHBody &sph_body);
    virtual ~RandomizeParticlePosition(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ParticleSmoothing
 * @brief computing smoothed variable field by averaging with neighbors
 */
template <typename VariableType>
class ParticleSmoothing : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    explicit ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name)
        : LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
          W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)),
          smoothed_(*particles_->template getVariableByName<VariableType>(variable_name))
    {
        particles_->registerVariable(temp_, variable_name + "_temp");
    }

    virtual ~ParticleSmoothing(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real weight = W0_;
        VariableType summation = W0_ * smoothed_[index_i];
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            summation += inner_neighborhood.W_ij_[n] * smoothed_[index_j];
            weight += inner_neighborhood.W_ij_[n];
        }
        temp_[index_i] = summation / (weight + TinyReal);
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        smoothed_[index_i] = temp_[index_i];
    };

  protected:
    const Real W0_;
    StdLargeVec<VariableType> &smoothed_, temp_;
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
 * @class 	UpperFrontInXDirection
 * @brief 	Get the upper front In X Direction for a SPH body
 *			TODO: a test using this method
 */
class UpperFrontInXDirection : public LocalDynamicsReduce<Real, ReduceMax>,
                               public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;

  public:
    explicit UpperFrontInXDirection(SPHBody &sph_body);
    virtual ~UpperFrontInXDirection(){};

    Real reduce(size_t index_i, Real dt = 0.0);
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
        this->quantity_name_ = variable_name + "Summation";
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
    TotalMechanicalEnergy(SPHBody &sph_body, SharedPtr<Gravity> = makeShared<Gravity>(Vecd::Zero()));
    virtual ~TotalMechanicalEnergy(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_DYNAMICS_H