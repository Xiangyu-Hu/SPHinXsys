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
} // namespace SPH
