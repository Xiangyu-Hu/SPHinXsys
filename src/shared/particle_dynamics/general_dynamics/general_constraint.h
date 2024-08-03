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
 * @file general_constraint.h
 * @brief Particles are constrained on their position according to
 * different criteria.
 * @author	Xiangyu Hu
 */

#ifndef GENERAL_CONSTRAINT_H
#define GENERAL_CONSTRAINT_H

#include "base_general_dynamics.h"

namespace SPH
{

template <class DynamicsIdentifier, typename DataType>
class ConstantConstraint : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    ConstantConstraint(DynamicsIdentifier &identifier, const std::string &variable_name,
                       DataType constrained_value)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          variable_data_field_(this->particles_->template getVariableDataByName<DataType>(variable_name)),
          constrained_value_(constrained_value){};
    virtual ~ConstantConstraint(){};
    void update(size_t index_i, Real dt = 0.0)
    {
        variable_data_field_[index_i] = constrained_value_;
    };

  protected:
    DataType *variable_data_field_;
    DataType constrained_value_;
};

class LevelSetShape;

/**
 * @class ShapeSurfaceBounding
 * @brief constrain surface particles by
 * map constrained particles to geometry face and
 * r = r + phi * norm (vector distance to face)
 */
class ShapeSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    ShapeSurfaceBounding(NearShapeSurface &body_part);
    virtual ~ShapeSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    LevelSetShape *level_set_shape_;
    Real constrained_distance_;
};

/**
 * @class	MotionConstraint
 * @brief	Base class for constraining with prescribed motion.
 * 			Exact motion function will be defined in derive class.
 * 			Note that, we do not impose acceleration, so that this constraint
 * 			must be imposed after updating particle velocity by forces
 * 			and before updating particle position.
 */
template <class DynamicsIdentifier>
class MotionConstraint : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    explicit MotionConstraint(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          pos0_(this->particles_->template registerStateVariableFrom<Vecd>("InitialPosition", "Position")),
          vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")){};

    virtual ~MotionConstraint(){};

  protected:
    Vecd *pos_, *pos0_, *vel_;
};
using BodyPartMotionConstraint = MotionConstraint<BodyPartByParticle>;

/**@class FixConstraint
 * @brief Constraint with zero velocity.
 */
template <class DynamicsIdentifier>
class FixConstraint : public MotionConstraint<DynamicsIdentifier>
{
  public:
    explicit FixConstraint(DynamicsIdentifier &identifier)
        : MotionConstraint<DynamicsIdentifier>(identifier){};
    virtual ~FixConstraint(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        this->pos_[index_i] = this->pos0_[index_i];
        this->vel_[index_i] = Vecd::Zero();
    };
};
using FixBodyConstraint = FixConstraint<SPHBody>;
using FixBodyPartConstraint = FixConstraint<BodyPartByParticle>;

/**
 * @class FixedInAxisDirection
 * @brief Constrain the velocity of a solid body part.
 */
class FixedInAxisDirection : public MotionConstraint<BodyPartByParticle>
{
  public:
    explicit FixedInAxisDirection(BodyPartByParticle &body_part, Vecd constrained_axises = Vecd::Zero())
        : MotionConstraint<BodyPartByParticle>(body_part), constrain_matrix_(Matd::Identity())
    {
        for (int k = 0; k != Dimensions; ++k)
            constrain_matrix_(k, k) = constrained_axises[k];
    };
    virtual ~FixedInAxisDirection(){};
    void update(size_t index_i, Real dt = 0.0)
    {
        this->vel_[index_i] = constrain_matrix_ * this->vel_[index_i];
    };

  protected:
    Matd constrain_matrix_;
};
} // namespace SPH
#endif // GENERAL_CONSTRAINT_H
