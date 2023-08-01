#ifndef GENERAL_CONTINUUM_H
#define GENERAL_CONTINUUM_H

#include "weakly_compressible_fluid.h"

namespace SPH
{
class GeneralContinuum : public WeaklyCompressibleFluid
{
  protected:
    Real E_;  /*< Youngs or tensile modules  */
    Real G_;  /*< shearmodules  */
    Real K_;  /*< bulkmodules  */
    Real nu_; /*< Poisson ratio  */
    Real contact_stiffness_; /**< contact-force stiffness related to bulk modulus*/

  public:
    explicit GeneralContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio)
        : WeaklyCompressibleFluid(rho0, c0), E_(0.0), G_(0.0), K_(0.0), nu_(0.0), contact_stiffness_(c0* c0)
    {
        material_type_name_ = "GeneralContinuum";
        E_ = youngs_modulus;
        nu_ = poisson_ratio;
        G_ = getShearModulus(youngs_modulus, poisson_ratio);
        K_ = getBulkModulus(youngs_modulus, poisson_ratio);
        lambda0_ = getLambda(youngs_modulus, poisson_ratio);
    };
    virtual ~GeneralContinuum(){};

    Real lambda0_; /*< first Lame parameter */
    Real getYoungsModulus() { return E_; };
    Real getPoissonRatio() { return nu_; };
    Real getDensity() { return rho0_; };
    Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
    Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
    Real getLambda(Real youngs_modulus, Real poisson_ratio);

    Real ContactStiffness() { return contact_stiffness_; };

    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress);

    virtual GeneralContinuum *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // GENERAL_CONTINUUM_H
