#include "elastic_solid.h"
#include "base_particles.hpp"

#ifdef max
#undef max
#endif

namespace SPH
{
    //=================================================================================================//
	OrthotropicSolid::OrthotropicSolid(Real rho_0, std::array<Vecd, 3> a, std::array<Real, 3> E,
									   std::array<Real, 3> G, std::array<Real, 3> poisson)
		// set parameters for parent class: LinearElasticSolid
		// we take the max. E and max. poisson to approximate the maximum of
		// the Bulk modulus --> for time step size calculation
		: LinearElasticSolid(rho_0, std::max({E[0], E[1], E[2]}),
							 std::max({poisson[0], poisson[1], poisson[2]})),
		  a_(a), E_(E), G_(G), poisson_(poisson)
	{
		// parameters for derived class
		material_type_name_ = "OrthotropicSolid";
		CalculateA0();
		CalculateAllMu();
		CalculateAllLambda();
	};
	//=================================================================================================//
	Matd OrthotropicSolid::StressPK2(Matd &F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (F.transpose() * F - Matd::Identity());
		Matd stress_PK2 = Matd::Zero();
		for (int i = 0; i < 3; i++)
		{
			// outer sum (a{1-3})
			Matd Summa2 = Matd::Zero();
			for (int j = 0; j < 3; j++)
			{
				// inner sum (b{1-3})
				Summa2 += Lambda_(i,j) * (CalculateDoubleDotProduct(A_[i], strain) * A_[j] +
										   CalculateDoubleDotProduct(A_[j], strain) * A_[i]);
			}
			stress_PK2 += Mu_[i] * (((A_[i] * strain) + (strain * A_[i])) + 1 / 2 * (Summa2));
		}
		return stress_PK2;
	}
	//=================================================================================================//
	Real OrthotropicSolid::VolumetricKirchhoff(Real J)
	{
		return K0_ * J * (J - 1);
	}
	//=================================================================================================//
	void OrthotropicSolid::CalculateA0()
	{
		A_[0] = a_[0] * a_[0].transpose();
		A_[1] = a_[1] * a_[1].transpose();
		A_[2] = a_[2] * a_[2].transpose();
	}
	//=================================================================================================//
	void OrthotropicSolid::CalculateAllMu()
	{

		Mu_[0] = 1 / G_[0] + 1 / G_[2] - 1 / G_[1];
		Mu_[1] = 1 / G_[1] + 1 / G_[0] - 1 / G_[2];
		Mu_[2] = 1 / G_[2] + 1 / G_[1] - 1 / G_[0];
	}
	//=================================================================================================//
	void OrthotropicSolid::CalculateAllLambda()
	{
		// first we calculate the upper left part, a 3x3 matrix of the full compliance matrix
		Matd Compliance = Matd::Zero();
		Compliance.col(0) = Vecd(1 / E_[0], -poisson_[0] / E_[1], -poisson_[1] / E_[2]);
		Compliance.col(1) = Vecd(-poisson_[0] / E_[0], 1 / E_[1], -poisson_[2] / E_[1]);
		Compliance.col(2) = Vecd(-poisson_[1] / E_[0], -poisson_[2] / E_[1], 1 / E_[2]);
		// we calculate the inverse of the Compliance matrix, and calculate the lambdas elementwise
		Matd Compliance_inv = Compliance.inverse();
		// Lambda_ is a 3x3 matrix
		Lambda_(0,0) = Compliance_inv(0,0) - 2 * Mu_[0];
		Lambda_(1,1) = Compliance_inv(1,1) - 2 * Mu_[1];
		Lambda_(2,2) = Compliance_inv(2,2) - 2 * Mu_[2];
		Lambda_(0,1) = Compliance_inv(0,1);
		Lambda_(0,2) = Compliance_inv(0,2);
		Lambda_(1,2) = Compliance_inv(1,2);
		// the matrix is symmetric
		Lambda_(1,0) = Lambda_(0,1);
		Lambda_(2,0) = Lambda_(0,2);
		Lambda_(2,1) = Lambda_(1,2);
	}
    //=================================================================================================//
}