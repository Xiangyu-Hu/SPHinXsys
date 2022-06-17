/**
 * @file 	loading_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "loading_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		ImposeExternalForce::
			ImposeExternalForce(SolidBody &solid_body, SolidBodyPartForSimbody &body_part)
			: PartSimpleDynamicsByParticle(solid_body, body_part), SolidDataSimple(solid_body),
			  pos_0_(particles_->pos_0_), vel_n_(particles_->vel_n_),
			  vel_ave_(particles_->vel_ave_) {}
		//=================================================================================================//
		void ImposeExternalForce::Update(size_t index_i, Real dt)
		{
			Vecd induced_acceleration = getAcceleration(pos_0_[index_i]);
			vel_n_[index_i] += induced_acceleration * dt;
			vel_ave_[index_i] = vel_n_[index_i];
		}
		//=================================================================================================//
		SpringDamperConstraintParticleWise::
			SpringDamperConstraintParticleWise(SolidBody &solid_body, Vecd stiffness, Real damping_ratio)
			: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_)
		{
			// scale stiffness and damping by mass here, so it's not necessary in each iteration
			stiffness_ = stiffness / std::accumulate(particles_->mass_.begin(), particles_->mass_.end(), 0.0);
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Vecd SpringDamperConstraintParticleWise::getSpringForce(size_t index_i, Vecd &disp)
		{
			Vecd spring_force(0);
			for (int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i];
			}
			return spring_force;
		}
		//=================================================================================================//
		Vecd SpringDamperConstraintParticleWise::getDampingForce(size_t index_i)
		{
			Vecd damping_force(0);
			for (int i = 0; i < vel_n_[index_i].size(); i++)
			{
				damping_force[i] = -damping_coeff_[i] * vel_n_[index_i][i];
			}
			return damping_force;
		}
		//=================================================================================================//
		void SpringDamperConstraintParticleWise::Update(size_t index_i, Real dt)
		{
			Vecd delta_x = pos_n_[index_i] - pos_0_[index_i];
			dvel_dt_prior_[index_i] += getSpringForce(index_i, delta_x);
			dvel_dt_prior_[index_i] += getDampingForce(index_i);
		}
		//=================================================================================================//
		SpringNormalOnSurfaceParticles::
			SpringNormalOnSurfaceParticles(SolidBody &solid_body, BodyPartByParticle &body_part, bool outer_surface, Vecd source_point, Real stiffness, Real damping_ratio)
			: PartSimpleDynamicsByParticle(solid_body, body_part), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  n_(particles_->n_),
			  n_0_(particles_->n_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  apply_spring_force_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(solid_body);
			// select which particles the spring is applied to
			for (size_t particle_i : surface_layer.body_part_particles_)
			{
				// vector to the source point from the particle
				Vecd vector_to_particle = source_point - pos_0_[particle_i];
				// normal of the particle
				Vecd normal = n_0_[particle_i];
				// get the cos of the angle between the vector and the normal
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
			Real area = std::pow(particles_->Vol_[0], 2.0 / 3.0);
			stiffness_ = stiffness * area;
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Vecd SpringNormalOnSurfaceParticles::getSpringForce(size_t index_i, Vecd disp)
		{
			// normal of the particle
			Vecd normal = particles_->n_0_[index_i];
			// get the normal portion of the displacement, which is parallel to the normal of particles, meaning it is the normal vector * scalar
			Vecd normal_disp = getVectorProjectionOfVector(disp, normal);

			Vecd spring_force_vector = -stiffness_ * normal_disp;

			return spring_force_vector;
		}
		//=================================================================================================//
		Vecd SpringNormalOnSurfaceParticles::getDampingForce(size_t index_i)
		{
			// normal of the particle
			Vecd normal = particles_->n_0_[index_i];
			// velocity of the particle
			Vecd velocity_n = vel_n_[index_i];
			// get the normal portion of the velocity, which is parallel to the normal of particles, meaning it is the normal vector * scalar
			Vecd normal_vel = getVectorProjectionOfVector(velocity_n, normal);

			Vecd damping_force_vector = -damping_coeff_ * normal_vel;

			return damping_force_vector;
		}
		//=================================================================================================//
		void SpringNormalOnSurfaceParticles::Update(size_t index_i, Real dt)
		{
			try
			{
				if (apply_spring_force_to_particle_[index_i])
				{
					Vecd delta_x = pos_n_[index_i] - pos_0_[index_i];
					dvel_dt_prior_[index_i] += getSpringForce(index_i, delta_x) / mass_[index_i];
					dvel_dt_prior_[index_i] += getDampingForce(index_i) / mass_[index_i];
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("SpringNormalOnSurfaceParticles::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		SpringOnSurfaceParticles::
			SpringOnSurfaceParticles(SolidBody &body, Real stiffness, Real damping_ratio)
			: ParticleDynamicsSimple(body), SolidDataSimple(body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  apply_spring_force_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(body);
			// select which particles the spring is applied to
			// if the particle is in the surface layer, the force is applied
			for (size_t particle_i : surface_layer.body_part_particles_)
				apply_spring_force_to_particle_[particle_i] = true;

			// scale stiffness and damping by area here, so it's not necessary in each iteration
			// we take the area of the first particle, assuming they are uniform
			Real area = std::pow(particles_->Vol_[0], 2.0 / 3.0);
			stiffness_ = stiffness * area;
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void SpringOnSurfaceParticles::Update(size_t index_i, Real dt)
		{
			try
			{
				if (apply_spring_force_to_particle_[index_i])
				{
					dvel_dt_prior_[index_i] += -stiffness_ * (pos_n_[index_i] - pos_0_[index_i]) / mass_[index_i];
					dvel_dt_prior_[index_i] += -damping_coeff_ * vel_n_[index_i] / mass_[index_i];
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("SpringOnSurfaceParticles::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		AccelerationForBodyPartInBoundingBox::
			AccelerationForBodyPartInBoundingBox(SolidBody &solid_body, BoundingBox &bounding_box, Vecd acceleration)
			: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  bounding_box_(bounding_box),
			  acceleration_(acceleration)
		{
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void AccelerationForBodyPartInBoundingBox::Update(size_t index_i, Real dt)
		{
			try
			{
				if (bounding_box_.checkContain(pos_n_[index_i]))
				{
					dvel_dt_prior_[index_i] += acceleration_;
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("AccelerationForBodyPartInBoundingBox::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		ForceInBodyRegion::
			ForceInBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part, Vecd force, Real end_time)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_0_(particles_->pos_0_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  acceleration_(0),
			  end_time_(end_time)
		{
			// calculate acceleration: force / total mass
			acceleration_ = force / std::accumulate(particles_->mass_.begin(), particles_->mass_.end(), 0.0);
			// set ghost particles to zero
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void ForceInBodyRegion::Update(size_t index_i, Real dt)
		{
			try
			{
				Real time_factor = SMIN(GlobalStaticVariables::physical_time_ / end_time_, 1.0);
				dvel_dt_prior_[index_i] = acceleration_ * time_factor;
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("ForceInBodyRegion::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		SurfacePressureFromSource::
			SurfacePressureFromSource(SPHBody &sph_body, BodyPartByParticle &body_part, Vecd source_point, StdVec<std::array<Real, 2>> pressure_over_time)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_0_(particles_->pos_0_),
			  n_(particles_->n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  pressure_over_time_(pressure_over_time),
			  // set apply_pressure_to_particle_ to false for each particle
			  apply_pressure_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(sph_body);
			// select which particles the pressure is applied to
			for (size_t particle_i : surface_layer.body_part_particles_)
			{
				// vector to the source point from the particle
				Vecd vector_to_particle = source_point - particles_->pos_0_[particle_i];
				// normal of the particle
				Vecd normal = particles_->n_0_[particle_i];
				// get the cos of the angle between the vector and the normal
				Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
				// if the angle is less than 90°, we apply the pressure to the surface particle
				// ignore exactly perpendicular surfaces
				if (cos_theta > 1e-6) 
				{
					apply_pressure_to_particle_[particle_i] = true;
				}
			}

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Real SurfacePressureFromSource::getPressure()
		{
			// check if we have reached the max time, if yes, return the last pressure
			bool max_time_reached = GlobalStaticVariables::physical_time_ > pressure_over_time_[pressure_over_time_.size() - 1][0];
			if (max_time_reached)
				return pressure_over_time_[pressure_over_time_.size() - 1][1];

			int interval = 0;
			// find out the interval
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
		void SurfacePressureFromSource::Update(size_t index_i, Real dt)
		{
			try
			{
				if (apply_pressure_to_particle_[index_i])
				{
					// get the surface area of the particle, assuming it has a cubic volume
					// acceleration is particle force / particle mass
					Real area = std::pow(particles_->Vol_[index_i], 2.0 / 3.0);
					Real acc_from_pressure = getPressure() * area / mass_[index_i];
					// vector is made by multiplying it with the surface normal
					// add the acceleration to the particle
					dvel_dt_prior_[index_i] += (-1.0) * n_[index_i] * acc_from_pressure;
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("SurfacePressureFromSource::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
	}
}
