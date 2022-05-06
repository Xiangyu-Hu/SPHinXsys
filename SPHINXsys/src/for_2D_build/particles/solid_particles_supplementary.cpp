/**
 * @file 	solid_particles_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "solid_particles.h"
#include "solid_particles_variable.h"
#include "base_body.h"

namespace SPH
{
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		Real J = rho0_ / rho_n_[particle_i];
		Mat2d F = F_[particle_i];
		Mat2d stress = stress_PK1_[particle_i];
		Mat2d sigma = (stress * ~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy + 3.0 * sigmaxy * sigmaxy);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain(size_t particle_i) // not tested in 2D
	{

		Mat2d F = F_[particle_i];
		Mat2d epsilon = 0.5 * (~F * F - Mat2d(1.0)); // calculation of the Green-Lagrange strain tensor

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
	//=============================================================================================//
	void TranslationAndRotation::update(size_t index_i, Real dt)
	{
		pos_n_[index_i] = transform_.imposeTransform(pos_n_[index_i]);
		pos_0_[index_i] = transform_.imposeTransform(pos_0_[index_i]);
	}
	//=============================================================================================//
	void VonMisesStress::update(size_t index_i, Real dt)
	{
		Real J = rho0_ / rho_n_[index_i];
		Mat2d F = F_[index_i];
		Mat2d stress = stress_PK1_[index_i];
		Mat2d sigma = (stress * ~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		derived_variable_[index_i] =
			sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy + 3.0 * sigmaxy * sigmaxy);
	}
	//=============================================================================================//
	void VonMisesStrain::update(size_t index_i, Real dt)
	{
		Mat2d F = F_[index_i];
		Mat2d epsilon = 0.5 * (~F * F - Mat2d(1.0)); // calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonxy = epsilon(0, 1);

		derived_variable_[index_i] = sqrt(0.5 * ((epsilonxx - epsilonyy) * (epsilonxx - epsilonyy) +
												 epsilonyy * epsilonyy + epsilonxx * epsilonxx +
												 2.0 * epsilonxy * epsilonxy));
	}
	//=================================================================================================//

}
