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
WeaklyCompressibleFluid::EosKernel::EosKernel(WeaklyCompressibleFluid &encloser)
    : Fluid::EosKernel(encloser), p0_(encloser.p0_) {}
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
Real HerschelBulkleyFluid::getViscosity(Real shear_rate)
{

    Real effective_shear_rate = SMAX(SMIN(shear_rate, max_shear_rate_), min_shear_rate_);
    return (yield_stress_ + consistency_index_ * (std::pow(effective_shear_rate, power_index_))) /
           effective_shear_rate;
}
//=================================================================================================//
Real CarreauFluid::getViscosity(Real shear_rate)
{
    Real effective_shear_rate = SMAX(SMIN(shear_rate, max_shear_rate_), min_shear_rate_);
    return mu_infty_ + (mu0_ - mu_infty_) * std::pow(1.0 + std::pow(characteristic_time_ * effective_shear_rate, 2), 0.5 * (power_index_ - 1.0));
}

//=================================================================================================//
} // namespace SPH
