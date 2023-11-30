#ifndef GENERAL_CONTINUUM_H
#define GENERAL_CONTINUUM_H

#include "weakly_compressible_fluid.h"

namespace SPH
{
class GeneralContinuum : public WeaklyCompressibleFluid
{
  protected:
    Real E_;  /*< Youngs or tensile modules  */
    Real G_;  /*< shear modules  */
    Real K_;  /*< bulk modules  */
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

class PlasticContinuum : public GeneralContinuum
{
protected:
    Real c_;  /*< cohesion  */
    Real phi_;  /*< friction angle  */
    Real psi_;  /*< dilatancy angle  */
    Real alpha_phi_;  /*< Drucker¨CPrager¡¯s constants  */
    Real k_c_;  /*< Drucker¨CPrager¡¯s constants */
public:
    explicit PlasticContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real friction_angle, Real cohesion=0, Real dilatancy=0)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
        c_(cohesion), phi_(friction_angle), psi_(dilatancy), alpha_phi_(0.0), k_c_(0.0)
    {
        material_type_name_ = "PlasticContinuum";
        alpha_phi_ = getDPConstantsA(friction_angle);
        k_c_ = getDPConstantsK(cohesion, friction_angle);
    };
    virtual ~PlasticContinuum() {};

    Real getDPConstantsA(Real friction_angle);
    Real getDPConstantsK(Real cohesion, Real friction_angle);
    Real getFrictionAngle() { return phi_; };

    virtual Mat3d ConstitutiveRelation(Mat3d& velocity_gradient, Mat3d& stress_tensor);
    virtual Mat3d ReturnMapping(Mat3d& stress_tensor);

    virtual GeneralContinuum* ThisObjectPtr() override { return this; };
};

class J2Plasticity : public GeneralContinuum
{
protected:
    Real yield_stress_;  /*< cohesion  */
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
public:
    explicit J2Plasticity(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real hardening_modulus = 0.0)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
        yield_stress_(yield_stress), hardening_modulus_(hardening_modulus)
    {
        material_type_name_ = "J2Plasticity";
    };
    virtual ~J2Plasticity() {};


    Real YieldStress() { return yield_stress_; };
    Real HardeningModulus() { return hardening_modulus_; };

    virtual Matd ConstitutiveRelationShearStress(Matd& velocity_gradient, Matd& shear_stress);
    virtual Matd ReturnMappingShearStress(Matd& shear_stress);
    virtual int PlasticIndicator(Matd& shear_stress);

    virtual J2Plasticity* ThisObjectPtr() override { return this; };
};

} // namespace SPH
#endif // GENERAL_CONTINUUM_H
