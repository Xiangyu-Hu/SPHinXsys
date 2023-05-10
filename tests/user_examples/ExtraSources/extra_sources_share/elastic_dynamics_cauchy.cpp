#include "elastic_dynamics_cauchy.h"
  
#include <numeric>

namespace SPH
{
	//=========================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		CauchyIntegration1stHalf::
			CauchyIntegration1stHalf(BaseInnerRelation &inner_relation)
			: Integration1stHalf(inner_relation) {}
		//=================================================================================================//
		void CauchyIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_[index_i] = rho0_ / F_[index_i].determinant();
			Matd inverse_F_T = F_[index_i].inverse().transpose();
			Matd almansi_strain = 0.5 * (Matd::Identity() - (F_[index_i] * F_[index_i].transpose()).inverse());
			//obtain the first Piola-Kirchhoff stress from the  Cauchy stress 
			stress_PK1_B_[index_i]  = F_[index_i].determinant() * elastic_solid_.StressCauchy(almansi_strain, F_[index_i], index_i) * 
									  inverse_F_T * B_[index_i];

		} 
		//=================================================================================================//
	}
}
