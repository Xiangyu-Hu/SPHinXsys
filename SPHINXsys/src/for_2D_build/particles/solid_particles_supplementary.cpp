#include "solid_particles.h"
#include "solid_particles_variable.h"
#include "base_body.h"

namespace SPH
{
	//=================================================================================================//
	Real ElasticSolidParticles::getVonMisesStrain(size_t particle_i) // not tested in 2D
	{

		Matd F = F_[particle_i];
		Matd epsilon = 0.5 * (F.transpose() * F - Matd::Identity()); // calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = 0; // z-components zero for 2D measures
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = 0; // z-components zero for 2D measures
		Real epsilonyz = 0; // z-components zero for 2D measures

		return sqrt((1.0 / 3.0) * (std::pow(epsilonxx - epsilonyy, 2.0) + std::pow(epsilonyy - epsilonzz, 2.0) +
								   std::pow(epsilonzz - epsilonxx, 2.0)) +
					2.0 * (std::pow(epsilonxy, 2.0) + std::pow(epsilonyz, 2.0) + std::pow(epsilonxz, 2.0)));
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getVonMisesStrainDynamic(size_t particle_i, Real poisson) // not tested in 2D
	{
		Mat2d F = F_[particle_i];
		Mat2d epsilon = 0.5 * (F.transpose() * F - Matd::Identity());  //calculation of the Green-Lagrange strain tensor
		Vec2d principal_strains = getPrincipalValuesFromMatrix(epsilon);

		Real eps_1 = principal_strains[0];
		Real eps_2 = principal_strains[1];

		return 1.0 / (1.0 + poisson) * std::sqrt(0.5 * (powerN(eps_1 - eps_2, 2)));
	}
	//=============================================================================================//
	void VonMisesStress::update(size_t index_i, Real dt)
	{
		Real J = rho0_ / rho_[index_i];

		Mat2d F = F_[index_i];
		Mat2d stress_PK1 = F * elastic_solid_.StressPK2(F, index_i);
		Mat2d sigma = (stress_PK1 * F.transpose() ) / J;


		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		derived_variable_[index_i] =
			sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy + 3.0 * sigmaxy * sigmaxy);
	}
	//=============================================================================================//
}
