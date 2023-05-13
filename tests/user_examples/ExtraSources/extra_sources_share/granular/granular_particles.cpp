#include "granular_particles.h"
#include "base_body.h"
#include "granular_material.h"
#include "xml_engine.h"
namespace SPH
{
	GranularMaterialParticles::
		GranularMaterialParticles(SPHBody &sph_body, GranularMaterial *granular_material)
		: FluidParticles(sph_body, granular_material), granular_material_(*granular_material) {}
	//=================================================================================================//
	void GranularMaterialParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();

		registerVariable(acc_shear_, "AccelerationByShear");
		registerVariable(stress_tensor_, "StressTensor");
		registerVariable(stress_tensor_rate_, "StressTensorRate");
		registerVariable(shear_stress_, "ShearStress");
		registerVariable(shear_stress_rate_, "ShearStressRate");
		registerVariable(von_mises_stress_, "VonMisesStress");
		registerVariable(velocity_gradient_, "VelocityGradient");
		registerVariable(strain_tensor_, "StrainTensor");
		registerVariable(strain_tensor_rate_, "StrainTensorRate");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------

		registerSortableVariable<Vecd>("AccelerationByShear");
		registerSortableVariable<Matd>("StressTensor");
		registerSortableVariable<Matd>("StressTensorRate");
		registerSortableVariable<Matd>("ShearStress");
		registerSortableVariable<Matd>("ShearStressRate");
		registerSortableVariable<Real>("VonMisesStress");
		registerSortableVariable<Matd>("VelocityGradient");
		registerSortableVariable<Matd>("StrainTensor");
		registerSortableVariable<Matd>("StrainTensorRate");
		//----------------------------------------------------------------------
		registerVariable(pos0_, "InitialPosition", [&](size_t i) -> Vecd
		{ return pos_[i]; });
		registerVariable(n_, "NormalDirection");
		registerVariable(n0_, "InitialNormalDirection",
			[&](size_t i) -> Vecd
		{ return n_[i]; });
	}

	Real GranularMaterialParticles::getVonMisesStress(Mat2d &stress_tensor)
	{
		Real sigmaxx = stress_tensor(0, 0);
		Real sigmayy = stress_tensor(1, 1);
		Real sigmaxy = stress_tensor(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy + 3.0 * sigmaxy * sigmaxy);
	}

	Real GranularMaterialParticles::getVonMisesStress(Mat3d& stress_tensor)
	{
		Real sigmaxx = stress_tensor(0, 0);
		Real sigmayy = stress_tensor(1, 1);
		Real sigmazz = stress_tensor(2, 2);
		Real sigmaxy = stress_tensor(0, 1);
		Real sigmaxz = stress_tensor(0, 2);
		Real sigmayz = stress_tensor(1, 2);

		return sqrt(0.5 * (pow(sigmaxx - sigmayy, 2) + pow(sigmaxx - sigmazz, 2) + pow(sigmayy - sigmazz, 2)) + 3 * (pow(sigmaxy, 2) + pow(sigmaxz, 2) + pow(sigmayz, 2)));
	}
}