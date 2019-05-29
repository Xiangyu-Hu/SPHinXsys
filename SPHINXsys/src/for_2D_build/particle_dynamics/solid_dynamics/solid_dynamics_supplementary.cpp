#include "solid_dynamics.h"
#include "solid_body.h"
#include "solid_body_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "base_data_package.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "mesh_cell_linked_list.h"
#include "weakly_compressible_fluid_particles.h"
#include "weakly_compressible_fluid.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=========================================================================================//
		void UpdateElasticNormalDirection::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//deformation tensor in 2D
			Real F00 = elastic_data_i.F_(0, 0);
			Real F01 = elastic_data_i.F_(0, 1);
			Real F10 = elastic_data_i.F_(1, 0);
			Real F11 = elastic_data_i.F_(1, 1);

			//polar decomposition
			Mat2d R = Mat2d(F00 + F11, F01 - F10, F10 - F01, F00 + F11);
			R = R / sqrt(pow(F00 + F11, 2) + pow(F01 - F10, 2));

			solid_data_i.n_ = R * solid_data_i.n_0_;
		}
		//=========================================================================================//
		void ConstrianSoildBodyPartBySimBody
			::ConstraintAParticle(size_t index_particle_i,
				Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vec3 rr, pos, vel, acc;
			rr(0) = solid_data_i.pos_0_[0] - initial_mobod_origin_location_[0];
			rr(1) = solid_data_i.pos_0_[1] - initial_mobod_origin_location_[1];
			rr(2) = 0.0;
			mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
			/** this is how we calculate the particle position in after transform of MBbody.
			 * const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
			 * const SimTK::Vec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
			 * const SimTK::Vec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
			 * base_particle_data_i.pos_n_ = (p_GB + r).getSubVec<2>(0);
			 */
			base_particle_data_i.pos_n_ = pos.getSubVec<2>(0);
			base_particle_data_i.vel_n_ = vel.getSubVec<2>(0);
			base_particle_data_i.dvel_dt_ = acc.getSubVec<2>(0);
		}
		//=========================================================================================//
		SpatialVec ForceOnSolidBodyPartForSimBody::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i 	= particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vec3 force_from_particle(0);
			force_from_particle.updSubVec<2>(0) = solid_data_i.force_from_fluid_;
			Vec3 displacement(0);
			displacement.updSubVec<2>(0) = base_particle_data_i.pos_n_ 
				- current_mobod_origin_location_.getSubVec<2>(0);
			Vec3 torque_from_particle = cross(displacement, force_from_particle);

			return SpatialVec(torque_from_particle, force_from_particle);
		}
		//=========================================================================================//
		SpatialVec ForceOnElasticBodyPartForSimBody
			::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vec3 force_from_particle(0);
			force_from_particle.updSubVec<2>(0) = elastic_data_i.mass_
				*base_particle_data_i.dvel_dt_;
			Vec3 displacement(0);
			displacement.updSubVec<2>(0) = base_particle_data_i.pos_n_ 
				- current_mobod_origin_location_.getSubVec<2>(0);
			Vec3 torque_from_particle = cross(displacement, force_from_particle);

			return SpatialVec(torque_from_particle, force_from_particle);
		}
		//=================================================================================================//	
	}
}
