#include "constraint_dynamics.h"
#include <numeric>

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
SpringConstrain::SpringConstrain(BodyPartByParticle &body_part, Real stiffness)
    : MotionConstraint<BodyPartByParticle>(body_part),
      stiffness_(stiffness * Vecd::Ones()),
      mass_(*particles_->getVariableDataByName<Real>("Mass")) {}
//=================================================================================================//
Vecd SpringConstrain::getAcceleration(Vecd &disp, Real mass)
{
    Vecd spring_force = Vecd::Zero();
    for (int i = 0; i < disp.size(); i++)
    {
        spring_force[i] = -stiffness_[i] * disp[i] / mass;
    }
    return spring_force;
}
//=================================================================================================//
void SpringConstrain::update(size_t index_i, Real dt)
{
    Vecd displacement = pos_[index_i] - pos0_[index_i];
    vel_[index_i] += dt * getAcceleration(displacement, mass_[index_i]);
}
//=================================================================================================//
PositionSolidBody::
    PositionSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd pos_end_center)
    : MotionConstraint<SPHBody>(sph_body),
      start_time_(start_time), end_time_(end_time), pos_end_center_(pos_end_center)
{
    BoundingBox bounds = sph_body.getSPHBodyBounds();
    pos_0_center_ = (bounds.first_ + bounds.second_) * 0.5;
    translation_ = pos_end_center_ - pos_0_center_;
}
//=================================================================================================//
Vecd PositionSolidBody::getDisplacement(size_t index_i, Real dt)
{
    // displacement from the initial position
    Vecd pos_final = pos0_[index_i] + translation_;
    return (pos_final - pos_[index_i]) * dt / (end_time_ - GlobalStaticVariables::physical_time_);
}
//=================================================================================================//
void PositionSolidBody::update(size_t index_i, Real dt)
{
    // only apply in the defined time period
    if (GlobalStaticVariables::physical_time_ >= start_time_ &&
        GlobalStaticVariables::physical_time_ <= end_time_)
    {
        pos_[index_i] = pos_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
        vel_[index_i] = Vecd::Zero();
    }
}
//=================================================================================================//
PositionScaleSolidBody::
    PositionScaleSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Real end_scale)
    : MotionConstraint<SPHBody>(sph_body),
      start_time_(start_time), end_time_(end_time), end_scale_(end_scale)
{
    BoundingBox bounds = sph_body.getSPHBodyBounds();
    pos_0_center_ = (bounds.first_ + bounds.second_) * 0.5;
}
//=================================================================================================//
Vecd PositionScaleSolidBody::getDisplacement(size_t index_i, Real dt)
{
    // displacement from the initial position
    Vecd pos_final = pos_0_center_ + end_scale_ * (pos0_[index_i] - pos_0_center_);
    return (pos_final - pos_[index_i]) * dt / (end_time_ - GlobalStaticVariables::physical_time_);
}
//=================================================================================================//
void PositionScaleSolidBody::update(size_t index_i, Real dt)
{
    // only apply in the defined time period
    if (GlobalStaticVariables::physical_time_ >= start_time_ &&
        GlobalStaticVariables::physical_time_ <= end_time_)
    {
        pos_[index_i] = pos_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
        vel_[index_i] = Vecd::Zero();
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
