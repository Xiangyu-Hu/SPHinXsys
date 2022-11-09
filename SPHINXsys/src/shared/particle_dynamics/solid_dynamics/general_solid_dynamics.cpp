/**
 * @file 	general_solid_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "general_solid_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		CorrectConfiguration::
			CorrectConfiguration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), SolidDataInner(inner_relation),
			  B_(particles_->B_) {}
		//=================================================================================================//
		void CorrectConfiguration::interaction(size_t index_i, Real dt)
		{
			Matd local_configuration(Eps); // a small number added to diagonal to avoid divide zero
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration -= SimTK::outer(r_ji, gradW_ijV_j);
			}
			B_[index_i] = SimTK::inverse(local_configuration);
		}
		//=================================================================================================//
	}
}
