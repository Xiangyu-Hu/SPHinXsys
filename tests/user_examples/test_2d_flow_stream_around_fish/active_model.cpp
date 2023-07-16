#include "active_model.h"
#include <numeric>

namespace SPH
{
	//=========================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		ActiveIntegration1stHalf::
			ActiveIntegration1stHalf(BaseInnerRelation &inner_relation)
			: Integration1stHalfPK2(inner_relation), active_strain_(*particles_->getVariableByName<Matd>("ActiveStrain"))
			, material_id_(*particles_->getVariableByName<int>("MaterialID"))
		{
			particles_->registerVariable(F_0, "ActiveTensor");
			particles_->registerVariable(E_e, "ElasticStrain");
		}
		//=================================================================================================//
		void ActiveIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			Integration1stHalfPK2::initialization(index_i,dt);

			// Calculate the active deformation gradient
			F_0[index_i] = Matd::Identity();
			F_0[index_i] = (material_id_[index_i] == 0) ? ((2.0 * active_strain_[index_i] + Matd::Identity()).llt().matrixL()) : F_0[index_i];
			
			//  Calculate the elastic deformation 
			Matd F_e = F_[index_i] * F_0[index_i].inverse();
			Matd F_0_star = F_0[index_i].determinant() * (F_0[index_i].inverse()).transpose();

			// Calculate the elastic strain
			E_e[index_i] = (material_id_[index_i] == 0) ? (0.5 * (F_[index_i].transpose() * F_[index_i] - Matd::Identity()) - active_strain_[index_i]) : F_[index_i];

			stress_PK1_B_[index_i] = F_e * elastic_solid_.StressPK2(E_e[index_i], index_i) * F_0_star * B_[index_i];
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
