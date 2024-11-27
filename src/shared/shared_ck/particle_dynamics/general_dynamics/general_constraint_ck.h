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
 * @file general_constraint_ck.h
 * @brief Particles are constrained on their position according to
 * different rules.
 * @author	Xiangyu Hu
 */

#ifndef GENERAL_CONSTRAINT_CK_H
#define GENERAL_CONSTRAINT_CK_H

#include "base_general_dynamics.h"

namespace SPH
{

template <class DynamicsIdentifier, class ConstraintType>
class MotionConstraintCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    explicit MotionConstraintCK(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_pos_(this->particles_->template getVariableDataByNameOnly<Vecd>("Position")),
          dv_pos0_(this->particles_->template registerStateVariableOnlyFrom<Vecd>("InitialPosition", "Position")),
          dv_vel_(this->particles_->template getVariableDataByNameOnly<Vecd>("Velocity")){};

    virtual ~MotionConstraintCK(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     MotionConstraintCK<DynamicsIdentifier, ConstraintType> &encloser)
            : constraint_(encloser.dv_pos_->DelegatedDataField(ex_policy),
                          encloser.dv_pos0_->DelegatedDataField(ex_policy),
                          encloser.dv_vel_->DelegatedDataField(ex_policy)){};
        void update(size_t index_i, Real dt = 0.0)
        {
            constraint_(index_i, dt);
        };

      protected:
        ConstraintType constraint_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos0_, *dv_vel_;
};

class FixConstraintCK
{
    Vecd *pos_, *pos0_, *vel_;

  public:
    FixConstraintCK(Vecd *pos, Vecd *pos0, Vecd *vel)
        : pos_(pos), pos0_(pos0), vel_(vel){};
    void operator()(size_t index_i, Real dt = 0.0)
    {
        pos_[index_i] = pos0_[index_i];
        vel_[index_i] = Vecd::Zero();
    };
};
using FixBodyPartConstraintCK = MotionConstraintCK<BodyPartByParticle, FixConstraintCK>;

} // namespace SPH
#endif // GENERAL_CONSTRAINT_CK_H
