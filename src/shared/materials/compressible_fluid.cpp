#include "compressible_fluid.h"

namespace SPH
{
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
