#ifndef GRANULAR_MATERIAL_H
#define GRANULAR_MATERIAL_H

#include "base_material.h"
#include "weakly_compressible_fluid.h"
#include <fstream>

namespace SPH
{
	class GranularMaterial : public WeaklyCompressibleFluid
	{
	protected:
		Real E_;  /*< Youngs or tensile modules  */
		Real G_;  /*< shearmodules  */
		Real K_;  /*< bulkmodules  */
		Real nu_;  /*< Poisson ratio  */

	public:
		explicit GranularMaterial(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio)
			: WeaklyCompressibleFluid(rho0, c0), E_(0.0), G_(0.0), K_(0.0), nu_(0.0)
		{
			material_type_name_ = "GranularMaterial";
			E_ = youngs_modulus;
			nu_ = poisson_ratio;
			G_ = getShearModulus(youngs_modulus, poisson_ratio);
			K_ = getBulkModulus(youngs_modulus, poisson_ratio);
			lambda0_ = getLambda(youngs_modulus, poisson_ratio);
		};
		virtual ~GranularMaterial() {};

		Real lambda0_; /*< first Lame parameter */
		Real getYoungsModulus() { return E_; };
		Real getPoissonRatio() { return nu_; };
		Real getDensity() { return rho0_; };
		Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
		Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
		Real getLambda(Real youngs_modulus, Real poisson_ratio);

		virtual Matd ConstitutiveRelationShearStress(Matd& velocity_gradient, Matd& shear_stress);

		virtual GranularMaterial *ThisObjectPtr() override { return this; };
	};
}
#endif //GRANULAR_MATERIAL_H
