#include "general_dynamics_correct.h"

namespace SPH
{
	//=================================================================================================//
	GlobalCorrectConfigurationInner::
		GlobalCorrectConfigurationInner(BaseInnerRelation& inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
		Vol_(particles_->Vol_)
	{
		particles_->registerVariable(A_, "OriginalCorrectionMatrix");
		particles_->registerVariable(B_, "WeightedCorrectionMatrix");
	}
	//=================================================================================================//
	void GlobalCorrectConfigurationInner::interaction(size_t index_i, Real dt)
	{
		Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
			local_configuration -= r_ji * gradW_ij.transpose();
		}
		A_[index_i] = local_configuration;
	}
	//=================================================================================================//
	void GlobalCorrectionMatrixComplex::interaction(size_t index_i, Real dt)
	{
		GlobalCorrectConfigurationInner::interaction(index_i, dt);

		Matd local_configuration = Eps * Matd::Identity();
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				Vecd gradW_ij = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
				Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
				local_configuration -= r_ji * gradW_ij.transpose();
			}
		}
		A_[index_i] += local_configuration;
		B_[index_i] = (pow(A_[index_i].determinant(), 2) * A_[index_i].inverse() + 0.3 * Matd::Identity()) / (0.3 * 1.0 + pow(A_[index_i].determinant(), 2));
	}
}
//=================================================================================================//