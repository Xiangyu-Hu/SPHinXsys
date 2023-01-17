#include "general_solid_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

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
			Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration -= r_ji* gradW_ijV_j.transpose();
			}
			B_[index_i] = local_configuration.inverse();
		}
		//=================================================================================================//
	}
}
