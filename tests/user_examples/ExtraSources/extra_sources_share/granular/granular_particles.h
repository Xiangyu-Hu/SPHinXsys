#ifndef GRANULAR_PARTICLES_H
#define GRANULAR_PARTICLES_H
#include "granular_material.h"
#include "fluid_particles.h"
#include "base_particles.h"
#include "particle_generator_lattice.h"
namespace SPH
{
	class GranularMaterial;

	class GranularMaterialParticles : public FluidParticles
	{
	public:

		StdLargeVec<Matd> strain_tensor_;
		StdLargeVec<Matd> strain_tensor_rate_;
		StdLargeVec<Vecd> acc_shear_;

		StdLargeVec<Matd> stress_tensor_;
		StdLargeVec<Matd> stress_tensor_rate_;

		StdLargeVec<Matd> shear_stress_;
		StdLargeVec<Matd> shear_stress_rate_;
		StdLargeVec<Matd> velocity_gradient_;

		StdLargeVec<Real> von_mises_stress_;

		StdLargeVec<Vecd> pos0_; /**< initial position */
		StdLargeVec<Vecd> n_;	  /**<  current normal direction */
		StdLargeVec<Vecd> n0_;	  /**<  initial normal direction */

		Real getVonMisesStress(Mat2d &stress_tensor);
		Real getVonMisesStress(Mat3d& stress_tensor);
		Real getDeviatoricPlasticStrain(Mat3d& strain_tensor);

		GranularMaterial &granular_material_;

		GranularMaterialParticles(SPHBody &sph_body, GranularMaterial *granular_material);
		virtual ~GranularMaterialParticles() {};

		virtual void initializeOtherVariables() override;
		virtual GranularMaterialParticles *ThisObjectPtr() override { return this; };
	};
}

#endif // GRANULAR_PARTICLES_H