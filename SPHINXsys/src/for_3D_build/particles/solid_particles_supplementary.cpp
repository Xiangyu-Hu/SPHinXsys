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
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		Real J = rho0_ / rho_n_[particle_i];
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_PK1_[particle_i];
		Mat3d sigma = (stress * ~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmazz = sigma(2, 2);
		Real sigmaxy = sigma(0, 1);
		Real sigmaxz = sigma(0, 2);
		Real sigmayz = sigma(1, 2);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz
			- sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz
			+ 3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
	}
	//=================================================================================================//
}
