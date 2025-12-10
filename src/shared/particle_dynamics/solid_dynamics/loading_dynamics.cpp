#include "loading_dynamics.h"
#include "base_general_dynamics.h"

#include <numeric>

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
SpringDamperConstraintParticleWise::
    SpringDamperConstraintParticleWise(SPHBody &sph_body, Vecd stiffness, Real damping_ratio)
    : LoadingForce(sph_body, "SpringDamperConstraintForce"),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      mass_(particles_->getVariableDataByName<Real>("Mass"))
{
    // scale stiffness and damping by mass here, so it's not necessary in each iteration
    stiffness_ = stiffness / std::accumulate(mass_, mass_ + particles_->TotalRealParticles(), Real(0.0));
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
    loading_force_[index_i] = getSpringForce(index_i, delta_x) * mass_[index_i] +
                              getDampingForce(index_i) * mass_[index_i];
    LoadingForce::update(index_i, dt);
}
//=================================================================================================//
SpringNormalOnSurfaceParticles::
    SpringNormalOnSurfaceParticles(SPHBody &sph_body, bool outer_surface,
                                   Vecd source_point, Real stiffness, Real damping_ratio)
    : LoadingForce(sph_body, "NormalSpringForceOnSurface"),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      n0_(particles_->registerStateVariableDataFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      is_spring_force_applied_(particles_->addUniqueDiscreteVariableData<bool>(
          "isSpringForceApplied", particles_->ParticlesBound(), false))
{
    BodySurface surface_layer(sph_body);

    for (size_t index_i : surface_layer.body_part_particles_)
    {
        Vecd vector_to_particle = source_point - pos0_[index_i];
        Vecd normal = n0_[index_i];
        Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
        // if outer surface, the normals close an angle greater than 90°
        // if the angle is greater than 90°, we apply the spring force to the surface particle
        Real epsilon = 1e-6; // to ignore exactly perpendicular surfaces
        if (outer_surface && cos_theta < -epsilon)
        {
            is_spring_force_applied_[index_i] = true;
        }
        // if not outer surface, it's inner surface, meaning the normals close an angle smaller than 90°
        // if the angle is less than 90°, we apply the spring force to the surface particle
        if (!outer_surface && cos_theta > epsilon)
        {
            is_spring_force_applied_[index_i] = true;
        }
    }
    // scale stiffness and damping by area here, so it's not necessary in each iteration
    // we take the area of the first particle, assuming they are uniform
    Real area = pow(Vol_[0], 2.0 / 3.0);
    stiffness_ = stiffness * area;
    damping_coeff_ = stiffness_ * damping_ratio;
}
//=================================================================================================//
Vecd SpringNormalOnSurfaceParticles::getSpringForce(size_t index_i, Vecd disp)
{
    Vecd normal = n0_[index_i];
    Vecd normal_disp = getVectorProjectionOfVector(disp, normal);
    Vecd spring_force_vector = -stiffness_ * normal_disp;

    return spring_force_vector;
}
//=================================================================================================//
Vecd SpringNormalOnSurfaceParticles::getDampingForce(size_t index_i)
{
    Vecd normal = n0_[index_i];
    Vecd velocity_n = vel_[index_i];
    Vecd normal_vel = getVectorProjectionOfVector(velocity_n, normal);
    Vecd damping_force_vector = -damping_coeff_ * normal_vel;

    return damping_force_vector;
}
//=================================================================================================//
void SpringNormalOnSurfaceParticles::update(size_t index_i, Real dt)
{
    if (is_spring_force_applied_[index_i])
    {
        Vecd delta_x = pos_[index_i] - pos0_[index_i];
        loading_force_[index_i] = getSpringForce(index_i, delta_x) +
                                  getDampingForce(index_i);
        LoadingForce::update(index_i, dt);
    }
}
//=================================================================================================//
SpringOnSurfaceParticles::
    SpringOnSurfaceParticles(SPHBody &sph_body, Real stiffness, Real damping_ratio)
    : LoadingForce(sph_body, "SpringForceOnSurface"),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      is_spring_force_applied_(particles_->addUniqueDiscreteVariableData<bool>(
          "isSpringForceApplied", particles_->ParticlesBound(), false))
{
    BodySurface surface_layer(sph_body);
    // select which particles the spring is applied to
    // if the particle is in the surface layer, the force is applied
    for (size_t index_i : surface_layer.body_part_particles_)
        is_spring_force_applied_[index_i] = true;

    // scale stiffness and damping by area here, so it's not necessary in each iteration
    // we take the area of the first particle, assuming they are uniform
    Real area = pow(Vol_[0], 2.0 / 3.0);
    stiffness_ = stiffness * area;
    damping_coeff_ = stiffness_ * damping_ratio;
}
//=================================================================================================//
void SpringOnSurfaceParticles::update(size_t index_i, Real dt)
{
    try
    {
        if (is_spring_force_applied_[index_i])
        {
            loading_force_[index_i] = -stiffness_ * (pos_[index_i] - pos0_[index_i]) -
                                      damping_coeff_ * vel_[index_i];
            ForcePrior::update(index_i, dt);
        }
    }
    catch (std::out_of_range &e)
    {
        throw std::runtime_error(std::string("SpringOnSurfaceParticles::Update: particle index out of bounds") + std::to_string(index_i));
    }
}
//=================================================================================================//
ExternalForceInBoundingBox::
    ExternalForceInBoundingBox(SPHBody &sph_body, BoundingBoxd &bounding_box, Vecd acceleration)
    : LoadingForce(sph_body, "ExternalForceInBoundingBox"),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      bounding_box_(bounding_box), acceleration_(acceleration) {}
//=================================================================================================//
void ExternalForceInBoundingBox::update(size_t index_i, Real dt)
{
    if (bounding_box_.checkContain(pos_[index_i]))
    {
        loading_force_[index_i] = acceleration_ * mass_[index_i];
        LoadingForce::update(index_i, dt);
    }
}
//=================================================================================================//
ForceInBodyRegion::
    ForceInBodyRegion(BodyPartByParticle &body_part, Vecd force, Real end_time)
    : BaseLoadingForce<BodyPartByParticle>(body_part, "ForceInBodyRegion"),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
      force_vector_(Vecd::Zero()), end_time_(end_time),
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))
{
    Real total_mass_in_region(0);
    for (size_t index_i : body_part.body_part_particles_)
        total_mass_in_region += mass_[index_i];
    force_vector_ = force;
}
//=================================================================================================//
void ForceInBodyRegion::update(size_t index_i, Real dt)
{
    Real time_factor = SMIN(*physical_time_ / end_time_, Real(1.0));
    loading_force_[index_i] = force_vector_ * time_factor;
    BaseLoadingForce<BodyPartByParticle>::update(index_i, dt);
}
//=================================================================================================//
SurfacePressureFromSource::
    SurfacePressureFromSource(BodyPartByParticle &body_part, Vecd source_point,
                              StdVec<std::array<Real, 2>> pressure_over_time)
    : BaseLoadingForce<BodyPartByParticle>(body_part, "SurfacePressureForce"),
      pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      pressure_over_time_(pressure_over_time),
      is_pressure_applied_(particles_->addUniqueDiscreteVariableData<bool>(
          "isPressureApplied", particles_->ParticlesBound(), false)),
      physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))
{
    BodySurface surface_layer(body_part.getSPHBody());

    for (size_t index_i : surface_layer.body_part_particles_)
    {
        Vecd vector_to_particle = source_point - pos0_[index_i];
        Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, n_[index_i]);
        // if the angle is less than 90°, we apply the pressure to the surface particle
        // ignore exactly perpendicular surfaces
        if (cos_theta > 1e-6)
        {
            is_pressure_applied_[index_i] = true;
        }
    }
}
//=================================================================================================//
Real SurfacePressureFromSource::getPressure()
{
    // check if we have reached the max time, if yes, return the last pressure
    bool max_time_reached = *physical_time_ > pressure_over_time_[pressure_over_time_.size() - 1][0];
    if (max_time_reached)
        return pressure_over_time_[pressure_over_time_.size() - 1][1];

    int interval = 0;
    for (size_t i = 0; i < pressure_over_time_.size(); i++)
    {
        if (*physical_time_ < pressure_over_time_[i][0])
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

    return p_0 + (p_1 - p_0) * (*physical_time_ - t_0) / (t_1 - t_0);
}
//=================================================================================================//
void SurfacePressureFromSource::update(size_t index_i, Real dt)
{
    if (is_pressure_applied_[index_i])
    {
        Real area = pow(Vol_[index_i], 2.0 / 3.0);
        Real acc_from_pressure = getPressure() * area / mass_[index_i];
        // vector is made by multiplying it with the surface normal
        // add the force to the particle
        loading_force_[index_i] = mass_[index_i] * (-1.0) * n_[index_i] * acc_from_pressure;
        BaseLoadingForce<BodyPartByParticle>::update(index_i, dt);
    }
}
//=================================================================================================//
PressureForceOnShell::PressureForceOnShell(SPHBody &sph_body, Real pressure)
    : LoadingForce(sph_body, "PressureForceOnShell"),
      pressure_(pressure),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")) {}
//=================================================================================================//
void PressureForceOnShell::update(size_t index_i, Real dt)
{
    loading_force_[index_i] = -pressure_ * Vol_[index_i] * n_[index_i];
    LoadingForce::update(index_i, dt);
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
