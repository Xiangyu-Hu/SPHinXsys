#include "weakly_compressible_fluid.h"

namespace SPH
{
//=================================================================================================//
Real WeaklyCompressibleFluid::getPressure(Real rho)
{
    return p0_ * (rho / rho0_ - 1.0);
}
//=================================================================================================//
Real WeaklyCompressibleFluid::DensityFromPressure(Real p)
{
    return rho0_ * (p / p0_ + 1.0);
}
//=================================================================================================//
Real WeaklyCompressibleFluid::getSoundSpeed(Real p, Real rho)
{
    return c0_;
}
//=================================================================================================//
Real SymmetricTaitFluid::getPressure(Real rho)
{
    Real rho_ratio = rho / rho0_;
    return rho_ratio > 1.0
               ? p0_ * (pow(rho_ratio, gamma_) - 1.0) / Real(gamma_)
               : -p0_ * (pow(1.0 / rho_ratio, gamma_) - 1.0) / Real(gamma_);
}
//=================================================================================================//
Real SymmetricTaitFluid::DensityFromPressure(Real p)
{
    return p > 0.0
               ? rho0_ * pow(1.0 + Real(gamma_) * p / p0_, 1.0 / Real(gamma_))
               : rho0_ / pow(1.0 - Real(gamma_) * p / p0_, 1.0 / Real(gamma_));
}
//=================================================================================================//
Real SymmetricTaitFluid::getSoundSpeed(Real p, Real rho)
{
    Real rho_ratio = rho / rho0_;
    return rho_ratio > 1.0
               ? sqrt((p0_ + Real(gamma_) * p) / rho)
               : sqrt((p0_ - Real(gamma_) * p) / rho);
}
//=================================================================================================//
Real HerschelBulkleyFluid::getViscosity(Real capped_shear_rate)
{
    return (this->yield_stress + this->consistency_index * (std::pow(capped_shear_rate, this->power_index))) / capped_shear_rate;
}
//=================================================================================================//
Real CarreauFluid::getViscosity(Real capped_shear_rate)
{
    return this->mu_infty + (this->mu_0 - this->mu_infty) * std::pow(1 + std::pow(this->characteristic_time * capped_shear_rate, 2), (this->power_index - 1) / 2);
}
//=================================================================================================//
} // namespace SPH
