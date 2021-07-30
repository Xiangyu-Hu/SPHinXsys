/**
 * @file 	solid_dynamics_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "all_solid_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "neighbor_relation.h"
#include "base_kernel.h"
#include "base_data_package.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "cell_linked_list.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=========================================================================================//
		void UpdateElasticNormalDirection::Update(size_t index_i, Real dt)
		{
			Matd& F = F_[index_i];
			//deformation tensor in 2D
			Real F00 = F(0, 0);
			Real F01 = F(0, 1);
			Real F10 = F(1, 0);
			Real F11 = F(1, 1);

			//polar decomposition
			Mat2d R = Mat2d(F00 + F11, F01 - F10, F10 - F01, F00 + F11);
			R = R / sqrt(pow(F00 + F11, 2) + pow(F01 - F10, 2));

			n_[index_i] = R * n_0_[index_i];
		}
		//=========================================================================================//
		void ConstrainSolidBodyPartBySimBody::Update(size_t index_i, Real dt)
		{
			Vecd& pos_0_i = pos_0_[index_i];
			Vec3 rr, pos, vel, acc;
			rr(0) = pos_0_i[0] - initial_mobod_origin_location_[0];
			rr(1) = pos_0_i[1] - initial_mobod_origin_location_[1];
			rr(2) = 0.0;
			mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
			/** this is how we calculate the particle position in after transform of MBbody.
			 * const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
			 * const SimTK::Vec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
			 * const SimTK::Vec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
			 * base_particle_data_i.pos_n_ = (p_GB + r).getSubVec<2>(0);
			 */
			pos_n_[index_i] = pos.getSubVec<2>(0);
			vel_n_[index_i] = vel.getSubVec<2>(0);
			dvel_dt_[index_i] = acc.getSubVec<2>(0);
			n_[index_i] = (mobod_.getBodyRotation(*simbody_state_) 
						* upgradeToVector3D(n_0_[index_i])).getSubVec<2>(0);
		}
		//=========================================================================================//
		SimTK::SpatialVec TotalForceOnSolidBodyPartForSimBody::ReduceFunction(size_t index_i, Real dt)
		{
			Vec3 force_from_particle(0);
			force_from_particle.updSubVec<2>(0) = force_from_fluid_[index_i] + contact_force_[index_i];
			Vec3 displacement(0);
			displacement.updSubVec<2>(0) = pos_n_[index_i]
				- current_mobod_origin_location_.getSubVec<2>(0);
			Vec3 torque_from_particle = cross(displacement, force_from_particle);

			return SimTK::SpatialVec(torque_from_particle, force_from_particle);
		}
		//=================================================================================================//	
	}
}
