#include "general_continuum.h"

namespace SPH
{
    //==============================GeneralContinuum===============================================//
    Real GeneralContinuum::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
    {
        return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
    }
    //=================================================================================================//
    Real GeneralContinuum::getShearModulus(Real youngs_modulus, Real poisson_ratio)
    {
        return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
    }
    //=================================================================================================//
    Real GeneralContinuum::getLambda(Real youngs_modulus, Real poisson_ratio)
    {
        return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
    }
    //=================================================================================================//
    Matd GeneralContinuum::ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress)
    {
        Matd strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
        Matd spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
        Matd deviatoric_strain_rate = strain_rate - (1.0 / (Real)Dimensions) * strain_rate.trace() * Matd::Identity();
        Matd stress_rate = 2.0 * G_ * deviatoric_strain_rate + shear_stress * (spin_rate.transpose()) + spin_rate * shear_stress;
        return stress_rate;
    }
} // namespace SPH