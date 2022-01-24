/**
 * @file 	solid_particles_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "solid_particles.h"
#include "base_body.h"

namespace SPH {
	//=============================================================================================//
	void SolidParticles::ParticleTranslationAndRotation(Transformd& transform) 
	{
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			pos_n_[i] = transform.imposeTransform(pos_n_[i]);
			pos_0_[i] = transform.imposeTransform(pos_0_[i]);
		}
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_GreenLagrange_strain(size_t particle_i)
	{
		Mat2d F = F_[particle_i];
		return 0.5 * (~F * F - Matd(1.0)); // calculation of the Green-Lagrange strain tensor
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::get_Principal_strains(size_t particle_i)
	{
		Mat2d epsilon = get_GreenLagrange_strain(particle_i); // calculation of the Green-Lagrange strain tensor
		return getPrincipalValuesFromMatrix(epsilon);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_static(size_t particle_i) //not tested in 2D
	{
		
		Mat2d epsilon = get_GreenLagrange_strain(particle_i); // calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = 0; 			//z-components zero for 2D measures
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = 0; 			//z-components zero for 2D measures
		Real epsilonyz = 0; 			//z-components zero for 2D measures

		return sqrt( (1.0 / 3.0) * (powerN(epsilonxx - epsilonyy, 2) + powerN(epsilonyy - epsilonzz, 2) + powerN(epsilonzz - epsilonxx, 2))
		 + 2.0 * (powerN(epsilonxy, 2) + powerN(epsilonyz, 2) + powerN(epsilonxz, 2)));
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_dynamic(size_t particle_i, Real poisson) //not tested in 2D
	{
		Mat2d F = F_[particle_i];
		Mat2d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor
		
		Vec2d principal_strains = getPrincipalValuesFromMatrix(epsilon);
		Real eps_1 = principal_strains[0];
		Real eps_2 = principal_strains[1];

		return 1.0/(1.0 + poisson) * std::sqrt(0.5 * (powerN(eps_1 - eps_2, 2)));
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_Cauchy_stress(size_t particle_i)
	{
		Real J = rho0_ / rho_n_[particle_i];
		Mat2d F = F_[particle_i];
		Mat2d stress = stress_PK1_[particle_i];

		return (stress * ~F) / J; // Cauchy stress
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_PK2_stress(size_t particle_i)
	{
		Mat2d F = F_[particle_i];
		Mat2d stress = stress_PK1_[particle_i];

		return SimTK::inverse(F) * stress; // Second Piola-Kirchhof stress
	}
	//=================================================================================================//
	Vec2d ElasticSolidParticles::get_Principal_stresses(size_t particle_i)
	{
		Mat2d sigma;
		if (stress_measure_ == "Cauchy") {
			sigma = get_Cauchy_stress(particle_i); // Cauchy stress
		} else if (stress_measure_ == "PK2") {
			sigma = get_PK2_stress(particle_i); // Second Piola-Kirchhof stress
		} else {
			throw std::runtime_error("get_Principal_stresses: wrong input");
		}

		return getPrincipalValuesFromMatrix(sigma);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::get_von_Mises_stress(size_t particle_i)
	{
		Mat2d sigma;
		if (stress_measure_ == "Cauchy") {
			sigma = get_Cauchy_stress(particle_i); // Cauchy stress
		} else if (stress_measure_ == "PK2") {
			sigma = get_PK2_stress(particle_i); // Second Piola-Kirchhof stress
		} else {
			throw std::runtime_error("get_von_Mises_stress: wrong input");
		}

		return getVonMisesStressFromMatrix(sigma);
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::displacement(size_t particle_i) //not tested in 2D
	{
		Vecd disp = pos_n_[particle_i] - pos_0_[particle_i];
		return disp;
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::normal(size_t particle_i) //not tested in 2D
	{
		Vecd normal_vec = n_[particle_i];
		return normal_vec;
	}
	//=================================================================================================//
}
