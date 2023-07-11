#include "loading_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
ImposeExternalForce::ImposeExternalForce(SPHBody &sph_body)
    : LocalDynamics(sph_body), SolidDataSimple(sph_body), pos0_(particles_->pos0_), vel_(particles_->vel_) {}
//=================================================================================================//
void ImposeExternalForce::update(size_t index_i, Real dt)
{
    vel_[index_i] += dt * getAcceleration(pos0_[index_i]);
}
//=================================================================================================//
SpringDamperConstraintParticleWise::SpringDamperConstraintParticleWise(SPHBody &sph_body, Vecd stiffness, Real damping_ratio)
    : LocalDynamics(sph_body), SolidDataSimple(sph_body), pos_(particles_->pos_), pos0_(particles_->pos0_),
      vel_(particles_->vel_), acc_prior_(particles_->acc_prior_)
{
    // scale stiffness and damping by mass here, so it's not necessary in each iteration
    stiffness_ = stiffness / std::accumulate(&particles_->mass_[0], &particles_->mass_[particles_->total_real_particles_], 0.0);
    damping_coeff_ = stiffness_ * damping_ratio;
}
//=================================================================================================//
Vecd SpringDamperConstraintParticleWise::getSpringForce(size_t index_i, Vecd &disp)
{
    Vecd spring_force = Vecd::Zero();
    for (int i = 0; i < disp.size(); i++)
    {
        spring_force[i] = -stiffness_[i] * disp[i];
    }
    return spring_force;
}
//=================================================================================================//
Vecd SpringDamperConstraintParticleWise::getDampingForce(size_t index_i)
{
    Vecd damping_force = Vecd::Zero();
    for (int i = 0; i < vel_[index_i].size(); i++)
    {
        damping_force[i] = -damping_coeff_[i] * vel_[index_i][i];
    }
    return damping_force;
}
//=================================================================================================//
void SpringDamperConstraintParticleWise::update(size_t index_i, Real dt)
{
    Vecd delta_x = pos_[index_i] - pos0_[index_i];
    acc_prior_[index_i] += getSpringForce(index_i, delta_x);
    acc_prior_[index_i] += getDampingForce(index_i);
}
//=================================================================================================//
SpringNormalOnSurfaceParticles::SpringNormalOnSurfaceParticles(SPHBody &sph_body, bool outer_surface,
                                                               Vecd source_point, Real stiffness, Real damping_ratio)
    : LocalDynamics(sph_body), SolidDataSimple(sph_body), pos_(particles_->pos_),
      pos0_(particles_->pos0_), n_(particles_->n_), n0_(particles_->n0_), vel_(particles_->vel_),
      acc_prior_(particles_->acc_prior_), mass_(particles_->mass_),
      apply_spring_force_to_particle_(StdLargeVec<bool>(pos0_.size(), false))
{
    BodySurface surface_layer(sph_body);

    for (size_t particle_i : surface_layer.body_part_particles_)
    {
        Vecd vector_to_particle = source_point - pos0_[particle_i];
        Vecd normal = n0_[particle_i];
        Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
        // if outer surface, the normals close an angle greater than 90°
        // if the angle is greater than 90°, we apply the spring force to the surface particle
        Real epsilon = 1e-6; // to ignore exactly perpendicular surfaces
        if (outer_surface && cos_theta < -epsilon)
        {
            apply_spring_force_to_particle_[particle_i] = true;
        }
        // if not outer surface, it's inner surface, meaning the normals close an angle smaller than 90°
        // if the angle is less than 90°, we apply the spring force to the surface particle
        if (!outer_surface && cos_theta > epsilon)
        {
            apply_spring_force_to_particle_[particle_i] = true;
        }
    }
    // scale stiffness and damping by area here, so it's not necessary in each iteration
    // we take the area of the first particle, assuming they are uniform
    Real area = pow(particles_->Vol_[0], 2.0 / 3.0);
    stiffness_ = stiffness * area;
    damping_coeff_ = stiffness_ * damping_ratio;
}
//=================================================================================================//
Vecd SpringNormalOnSurfaceParticles::getSpringForce(size_t index_i, Vecd disp)
{
    Vecd normal = particles_->n0_[index_i];
    Vecd normal_disp = getVectorProjectionOfVector(disp, normal);
    Vecd spring_force_vector = -stiffness_ * normal_disp;

    return spring_force_vector;
}
//=================================================================================================//
Vecd SpringNormalOnSurfaceParticles::getDampingForce(size_t index_i)
{
    Vecd normal = particles_->n0_[index_i];
    Vecd velocity_n = vel_[index_i];
    Vecd normal_vel = getVectorProjectionOfVector(velocity_n, normal);
    Vecd damping_force_vector = -damping_coeff_ * normal_vel;

    return damping_force_vector;
}
//=================================================================================================//
void SpringNormalOnSurfaceParticles::update(size_t index_i, Real dt)
{
    if (apply_spring_force_to_particle_[index_i])
    {
        Vecd delta_x = pos_[index_i] - pos0_[index_i];
        acc_prior_[index_i] += getSpringForce(index_i, delta_x) / mass_[index_i];
        acc_prior_[index_i] += getDampingForce(index_i) / mass_[index_i];
    }
}
//=================================================================================================//
SpringOnSurfaceParticles::SpringOnSurfaceParticles(SPHBody &sph_body, Real stiffness, Real damping_ratio)
    : LocalDynamics(sph_body), SolidDataSimple(sph_body), pos_(particles_->pos_), pos0_(particles_->pos0_),
      vel_(particles_->vel_), acc_prior_(particles_->acc_prior_), mass_(particles_->mass_),
      apply_spring_force_to_particle_(StdLargeVec<bool>(pos0_.size(), false))
{
    BodySurface surface_layer(sph_body);
    // select which particles the spring is applied to
    // if the particle is in the surface layer, the force is applied
    for (size_t particle_i : surface_layer.body_part_particles_)
        apply_spring_force_to_particle_[particle_i] = true;

    // scale stiffness and damping by area here, so it's not necessary in each iteration
    // we take the area of the first particle, assuming they are uniform
    Real area = pow(particles_->Vol_[0], 2.0 / 3.0);
    stiffness_ = stiffness * area;
    damping_coeff_ = stiffness_ * damping_ratio;
}
//=================================================================================================//
void SpringOnSurfaceParticles::update(size_t index_i, Real dt)
{
    try
    {
        if (apply_spring_force_to_particle_[index_i])
        {
            acc_prior_[index_i] += -stiffness_ * (pos_[index_i] - pos0_[index_i]) / mass_[index_i];
            acc_prior_[index_i] += -damping_coeff_ * vel_[index_i] / mass_[index_i];
        }
    }
    catch (std::out_of_range &e)
    {
        throw std::runtime_error(std::string("SpringOnSurfaceParticles::Update: particle index out of bounds") + std::to_string(index_i));
    }
}
//=================================================================================================//
AccelerationForBodyPartInBoundingBox::AccelerationForBodyPartInBoundingBox(SPHBody &sph_body, BoundingBox &bounding_box, Vecd acceleration)
    : LocalDynamics(sph_body), SolidDataSimple(sph_body), pos_(particles_->pos_),
      acc_prior_(particles_->acc_prior_), bounding_box_(bounding_box), acceleration_(acceleration) {}
//=================================================================================================//
void AccelerationForBodyPartInBoundingBox::update(size_t index_i, Real dt)
{
    if (bounding_box_.checkContain(pos_[index_i]))
    {
        acc_prior_[index_i] += acceleration_;
    }
}
//=================================================================================================//
ForceInBodyRegion::ForceInBodyRegion(BodyPartByParticle &body_part, Vecd force, Real end_time)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), SolidDataSimple(sph_body_),
      pos0_(particles_->pos0_), acc_prior_(particles_->acc_prior_), acceleration_(Vecd::Zero()), end_time_(end_time)
{
    Real total_mass_in_region(0);
    for (size_t particle_i : body_part.body_part_particles_)
        total_mass_in_region += particles_->mass_[particle_i];
    acceleration_ = force / total_mass_in_region;
}
//=================================================================================================//
void ForceInBodyRegion::update(size_t index_i, Real dt)
{
    Real time_factor = SMIN(GlobalStaticVariables::physical_time_ / end_time_, Real(1.0));
    acc_prior_[index_i] = acceleration_ * time_factor;
}
//=================================================================================================//
SurfacePressureFromSource::SurfacePressureFromSource(BodyPartByParticle &body_part, Vecd source_point,
                                                     StdVec<std::array<Real, 2>> pressure_over_time)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), SolidDataSimple(sph_body_),
      pos0_(particles_->pos0_), n_(particles_->n_), acc_prior_(particles_->acc_prior_),
      mass_(particles_->mass_), pressure_over_time_(pressure_over_time),
      apply_pressure_to_particle_(StdLargeVec<bool>(pos0_.size(), false))
{
    BodySurface surface_layer(sph_body_);

    for (size_t particle_i : surface_layer.body_part_particles_)
    {
        Vecd vector_to_particle = source_point - particles_->pos0_[particle_i];
        Vecd normal = particles_->n0_[particle_i];
        Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
        // if the angle is less than 90°, we apply the pressure to the surface particle
        // ignore exactly perpendicular surfaces
        if (cos_theta > 1e-6)
        {
            apply_pressure_to_particle_[particle_i] = true;
        }
    }
}
//=================================================================================================//
Real SurfacePressureFromSource::getPressure()
{
    // check if we have reached the max time, if yes, return the last pressure
    bool max_time_reached = GlobalStaticVariables::physical_time_ > pressure_over_time_[pressure_over_time_.size() - 1][0];
    if (max_time_reached)
        return pressure_over_time_[pressure_over_time_.size() - 1][1];

    int interval = 0;
    for (size_t i = 0; i < pressure_over_time_.size(); i++)
    {
        if (GlobalStaticVariables::physical_time_ < pressure_over_time_[i][0])
        {
            interval = i;
            break;
        }
    }
    // interval has to be at least 1
    if (interval < 1)
        throw std::runtime_error(std::string("SurfacePressureFromSource::getPressure(): pressure_over_time input not correct, should start with {0.0, 0.0}"));
    // scale the pressure to the current time
    Real t_0 = pressure_over_time_[interval - 1][0];
    Real t_1 = pressure_over_time_[interval][0];
    Real p_0 = pressure_over_time_[interval - 1][1];
    Real p_1 = pressure_over_time_[interval][1];

    return p_0 + (p_1 - p_0) * (GlobalStaticVariables::physical_time_ - t_0) / (t_1 - t_0);
}
//=================================================================================================//
void SurfacePressureFromSource::update(size_t index_i, Real dt)
{
    if (apply_pressure_to_particle_[index_i])
    {
        Real area = pow(particles_->Vol_[index_i], 2.0 / 3.0);
        Real acc_from_pressure = getPressure() * area / mass_[index_i];
        // vector is made by multiplying it with the surface normal
        // add the acceleration to the particle
        acc_prior_[index_i] += (-1.0) * n_[index_i] * acc_from_pressure;
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
