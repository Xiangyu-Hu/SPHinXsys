/**
 * @file 	inelastic_dynamics.cpp
 * @author	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */

#include "inelastic_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		PlasticStressRelaxationFirstHalf::
			PlasticStressRelaxationFirstHalf(BaseBodyRelationInner* body_inner_relation) :
			StressRelaxationFirstHalf(body_inner_relation),
			plastic_solid_(dynamic_cast<PlasticSolid*>(material_))
		{
			numerical_dissipation_factor_ = 0.25;
		}
		//=================================================================================================//
		void PlasticStressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho0_ / SimTK::det(F_[index_i]);

			stress_PK1_[index_i] = plastic_solid_->PlasticConstitutiveRelation(F_[index_i], index_i, dt) * B_[index_i];
		}
		//=================================================================================================//
	}
}
