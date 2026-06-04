#include "weakly_compressible_fluid.h"

namespace SPH
{
//=================================================================================================//
WeaklyCompressibleFluid::WeaklyCompressibleFluid(Real rho0, Real c0)
    : Fluid(), rho0_(rho0), c0_(c0), p0_(rho0 * c0 * c0)
{
    material_type_name_ = "WeaklyCompressibleFluid";
}
//=================================================================================================//
WeaklyCompressibleFluid::WeaklyCompressibleFluid(ConstructArgs<Real, Real> args)
    : WeaklyCompressibleFluid(std::get<0>(args), std::get<1>(args)) {}
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
} // namespace SPH
