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
#include "polar_decomposition_3x3.h"

using namespace polar;
using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=========================================================================================//
		void UpdateElasticNormalDirection::Update(size_t index_i, Real dt)
		{
			Matd& F = F_[index_i];
			Mat3d R;
			Real Q[9], H[9], A[9];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					A[i * 3 + j] = F(i, j);

			polar::polar_decomposition(Q, H, A);
			//this decomposition has the form A = Q*H, where Q is orthogonal and H is symmetric positive semidefinite. 
			//Ref. "An algorithm to compute the polar decomposition of a 3*3 matrix, Nicholas J. Higham et al. Numer Algor(2016) "
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					R(i, j) = Q[i * 3 + j];
			n_[index_i] = R * n_0_[index_i];
		}
		//=================================================================================================//
		void ConstrainSolidBodyPartBySimBody::Update(size_t index_i, Real dt)
		{
			Vec3 rr, pos, vel, acc;
			rr = pos_0_[index_i] - initial_mobod_origin_location_;
			mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
			/** this is how we calculate the particle position in after transform of MBbody.
			 * const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
			 * const SimTK::Vec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
			 * const SimTK::Vec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
			 * base_particle_data_i.pos_n_ = (p_GB + r);
			 */
			pos_n_[index_i] = pos;
			vel_n_[index_i] = vel;
			dvel_dt_[index_i] = acc;
			n_[index_i] = (mobod_.getBodyRotation(*simbody_state_) * n_0_[index_i]);
		}
		//=================================================================================================//
		SimTK::SpatialVec TotalForceOnSolidBodyPartForSimBody
			::ReduceFunction(size_t index_i, Real dt)
		{
			Vec3 force_from_particle = force_from_fluid_[index_i] + contact_force_[index_i];
			Vec3 displacement(0);
			displacement = pos_n_[index_i] - current_mobod_origin_location_;
			Vec3 torque_from_particle = cross(displacement, force_from_particle);

			return SimTK::SpatialVec(torque_from_particle, force_from_particle);
		}
		//=================================================================================================//	
	}
}