#include "all_solid_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "neighborhood.h"
#include "base_kernel.h"
#include "base_data_package.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "cell_linked_list.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"

#include "Simbody.h"
#include "SimTKmath.h"
#include "SimTKcommon.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//=========================================================================================//
		void UpdateElasticNormalDirection::update(size_t index_i, Real dt)
		{
			Matd &F = F_[index_i];
			// deformation tensor in 2D
			Real F00 = F(0, 0);
			Real F01 = F(0, 1);
			Real F10 = F(1, 0);
			Real F11 = F(1, 1);
			// polar decomposition
			Mat2d R{
				{F00 + F11, F01 - F10},
				{F10 - F01, F00 + F11},
			};

			n_[index_i] = R * n0_[index_i] / sqrt(pow(F00 + F11, 2) + pow(F01 - F10, 2));
		}
		//=========================================================================================//
		void ConstraintBySimBody::update(size_t index_i, Real dt)
		{
			/** Change vector to Simbody. */
			SimTK::Vec3 rr, pos, vel, acc;
			rr = SimTK::Vec3(pos0_[index_i][0], pos0_[index_i][1], 0.0) - initial_mobod_origin_location_;
			mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
			/** this is how we calculate the particle position in after transform of MBbody.
			 *  const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
			 *  const SimTK::Vec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
			 *  const SimTK::Vec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
			 *  base_particle_data_i.pos_ = (p_GB + r).getSubVec<2>(0);
			 */
			SimTK::Vec3 n = (mobod_.getBodyRotation(*simbody_state_) * SimTK::Vec3(n0_[index_i][0], n0_[index_i][1], 0.0));
			/** Change vector to eigen. */
			pos_[index_i] = Vecd(pos[0], pos[1]);
			vel_[index_i] = Vecd(vel[0], vel[1]);
			n_[index_i] = Vecd(n[0], n[1]);
		}
		//=========================================================================================//
		SimTK::SpatialVec TotalForceForSimBody::reduce(size_t index_i, Real dt)
		{
			Vecd force = (acc_[index_i] + acc_prior_[index_i]) * mass_[index_i];
			SimTK::Vec3 force_from_particle(force[0], force[1], 0.0);
			SimTK::Vec3 displacement = SimTK::Vec3(pos_[index_i][0], pos_[index_i][1], 0.0) - current_mobod_origin_location_;
			SimTK::Vec3 torque_from_particle = SimTK::cross(displacement, force_from_particle);

			return SimTK::SpatialVec(torque_from_particle, force_from_particle);
		}
		//=================================================================================================//
	}
}
