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
 * @file 	free_stream_boundary.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef FREE_STREAM_BOUNDARY_H
#define FREE_STREAM_BOUNDARY_H

#include "fluid_boundary_state.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class   FreeStreamCondition
 * @brief   modify the velocity of free surface particles with far-field velocity
 *          TargetVelocity gives the velocity profile along the free-stream direction,
 *          i.e. x direction in local frame.
 */
template <typename TargetVelocity>
class FreeStreamCondition : public LocalDynamics
{
  protected:
    Transform transform_;
    Real rho0_;
    Real *rho_sum_;
    Vecd *pos_, *vel_;
    int *indicator_;
    TargetVelocity target_velocity;
    Real *physical_time_;

  public:
    explicit FreeStreamCondition(SPHBody &sph_body, const Transform &transform = Transform())
        : LocalDynamics(sph_body),
          transform_(transform), rho0_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceDensity()),
          rho_sum_(particles_->getVariableDataByName<Real>("DensitySummation")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          indicator_(particles_->getVariableDataByName<int>("Indicator")),
          target_velocity(*this),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {};
    virtual ~FreeStreamCondition() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i] == 1)
        {
            Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
            Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
            Real frame_u_stream_direction = frame_velocity[0];
            Real u_freestream = target_velocity(frame_position, frame_velocity, *physical_time_)[0];
            frame_velocity[0] = u_freestream + (frame_u_stream_direction - u_freestream) *
                                                   SMIN(rho_sum_[index_i], rho0_) / rho0_;
            vel_[index_i] = transform_.xformFrameVecToBase(frame_velocity);
        }
    };
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FREE_STREAM_BOUNDARY_H