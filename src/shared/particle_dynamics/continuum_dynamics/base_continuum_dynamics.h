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
 * @file base_continuum_dynamics.h
 * @brief Collection of headers and types used by all continuum dynamics classes.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#ifndef BASE_CONTINUUM_DYNAMICS_H
#define BASE_CONTINUUM_DYNAMICS_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace continuum_dynamics
{
/**
 * @class InteractionWithWall
 * @brief Base class adding interaction with wall to general relaxation process
 */
template <template <typename...> class BaseInteractionType>
class InteractionWithWall : public BaseInteractionType<DataDelegateContact>
{
  public:
    explicit InteractionWithWall(BaseContactRelation &wall_contact_relation)
        : BaseInteractionType<DataDelegateContact>(wall_contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            Solid &solid_material = DynamicCast<Solid>(this, this->contact_particles_[k]->getBaseMaterial());
            wall_vel_ave_.push_back(solid_material.AverageVelocity(this->contact_particles_[k]));
            wall_acc_ave_.push_back(solid_material.AverageAcceleration(this->contact_particles_[k]));
            wall_n_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
            wall_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<Vecd *> wall_vel_ave_, wall_acc_ave_, wall_n_;
    StdVec<Real *> wall_Vol_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // BASE_CONTINUUM_DYNAMICS_H
