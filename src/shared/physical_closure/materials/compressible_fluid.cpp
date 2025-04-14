#include "compressible_fluid.h"

namespace SPH
{
//=============================================================================================//
CompressibleFluid::CompressibleFluid(Real rho0, Real gamma)
    : Fluid(rho0, 1.0), gamma_(gamma)
{
    material_type_name_ = "CompressibleFluid";
}
//=============================================================================================//
CompressibleFluid::CompressibleFluid(ConstructArgs<Real, Real> args)
    : CompressibleFluid(std::get<0>(args), std::get<1>(args)) {}
//=============================================================================================//
Real CompressibleFluid::getPressure(Real rho, Real rho_e)
{
    return rho_e * (gamma_ - 1.0);
}
//=============================================================================================//
Real CompressibleFluid::getSoundSpeed(Real p, Real rho)
{
    return std::sqrt(gamma_ * p / rho);
}
//=============================================================================================//
} // namespace SPH
