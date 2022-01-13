#include "solid_particles.h"
#include "base_body.h"

#include <iterator>

namespace SPH {
	//=============================================================================================//
	void SolidParticles::ParticleTranslationAndRotation(Transformd& transform) 
	{
			std::cout << "\n Error: the function ParticleTranslationAndRotation in 3d is not defined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_static(size_t particle_i)
	{
		Mat3d F = F_[particle_i];
		Mat3d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = epsilon(2, 2);
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = epsilon(0, 2);
		Real epsilonyz = epsilon(1, 2);

		return sqrt( (1.0 / 3.0) * (std::pow(epsilonxx - epsilonyy, 2.0) + std::pow(epsilonyy - epsilonzz, 2.0) + std::pow(epsilonzz - epsilonxx, 2.0))
		 + 2.0 * (std::pow(epsilonxy, 2.0) + std::pow(epsilonyz, 2.0) + std::pow(epsilonxz, 2.0)));
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_dynamic(size_t particle_i, Real poisson)
	{
		// calculation based on the Ansys formulation
		Mat3d F = F_[particle_i];
		Mat3d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor
		
		Vec3d principal_strains = getPrincipalValuesFromMatrix(epsilon);
		Real eps_1 = principal_strains[0];
		Real eps_2 = principal_strains[1];
		Real eps_3 = principal_strains[2];

		return 1.0/(1.0 + poisson) * std::sqrt(0.5 * (std::pow(eps_1 - eps_2, 2) + std::pow(eps_2 - eps_3, 2) + std::pow(eps_3 - eps_1, 2)));
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress_Cauchy(size_t particle_i)
	{
		Real J = rho0_ / rho_n_[particle_i];
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_PK1_[particle_i];
		Mat3d sigma = (stress * ~F) / J; // Cauchy stress

		return getVonMisesStressFromMatrix(sigma);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress_PK2(size_t particle_i)
	{
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_PK1_[particle_i];
		Mat3d sigma = SimTK::inverse(F) * stress; // Second Piola-Kirchhof stress

		return getVonMisesStressFromMatrix(sigma);
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::displacement(size_t particle_i)
	{
		Vecd disp = pos_n_[particle_i] - pos_0_[particle_i];
		return disp;
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::normal(size_t particle_i)
	{
		Vecd normal_vec = n_[particle_i];
		return normal_vec;
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::displacement(size_t particle_i)
	{
		Vecd disp = pos_n_[particle_i] - pos_0_[particle_i];
		return disp;
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::normal(size_t particle_i)
	{
		Vecd normal_vec = n_[particle_i];
		return normal_vec;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain(size_t particle_i)
	{
		
		Mat3d F = F_[particle_i];
		Mat3d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor
		

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = epsilon(2, 2);
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = epsilon(0, 2);
		Real epsilonyz = epsilon(1, 2);

		return sqrt( (1.0 / 3.0) * (std::pow(epsilonxx - epsilonyy, 2.0) + std::pow(epsilonyy - epsilonzz, 2.0) + std::pow(epsilonzz - epsilonxx, 2.0))
		 + 2.0 * (std::pow(epsilonxy, 2.0) + std::pow(epsilonyz, 2.0) + std::pow(epsilonxz, 2.0)));
	}
	//=================================================================================================//
}
