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
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		Real J = rho_0_ / rho_n_[particle_i];
		Mat2d F = F_[particle_i];
		Mat2d stress = stress_PK1_[particle_i];
		Mat2d sigma = (stress * ~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy
			+ 3.0 * sigmaxy * sigmaxy);
	}
	//=================================================================================================//
}
