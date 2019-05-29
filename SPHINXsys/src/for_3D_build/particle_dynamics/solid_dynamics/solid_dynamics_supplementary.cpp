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
#include "polar_decomposition_3x3.h"

using namespace polar;
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

			Mat3d R;
			Real Q[9], H[9], A[9];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					A[i * 3 + j] = elastic_data_i.F_(i, j);

			polar::polar_decomposition(Q, H, A);
			//this decomposition has the form A = Q*H, where Q is orthogonal and H is symmetric positive semidefinite. 
			//Ref. "An algorithm to compute the polar decomposition of a 3*3 matrix, Nicholas J. Higham et al. Numer Algor(2016) "
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					R(i, j) = Q[i * 3 + j];
			solid_data_i.n_ = R * solid_data_i.n_0_;
		}
		//=================================================================================================//
		void ConstrianSoildBodyPartBySimBody
			::ConstraintAParticle(size_t index_particle_i,
				Real dt)
		{

		}

		SpatialVec ForceOnSolidBodyPartForSimBody
			::ReduceFunction(size_t index_particle_i, Real dt)
		{
			cout << "\n This function is not done in 3D. Exit the program! \n";
			exit(0);
			return SpatialVec(Vec3(0), Vec3(0));
		}

		SpatialVec ForceOnElasticBodyPartForSimBody
			::ReduceFunction(size_t index_particle_i, Real dt)
		{
			cout << "\n This function is not done in 3D. Exit the program! \n";
			exit(0);
			return SpatialVec(Vec3(0), Vec3(0));
		}
		//=================================================================================================//	
	}
}