#include "viscosity.h"

namespace SPH
{
//=================================================================================================//
OldroydBViscosity::OldroydBViscosity(Real mu, Real lambda, Real mu_p)
    : Viscosity(mu), lambda_(lambda), mu_p_(mu_p) {}
//=================================================================================================//
OldroydBViscosity::OldroydBViscosity(ConstructArgs<Real, Real, Real> args)
    : OldroydBViscosity(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {}
//=================================================================================================//
GeneralizedNewtonianViscosity::GeneralizedNewtonianViscosity(Real min_shear_rate, Real max_shear_rate)
    : Viscosity(max_shear_rate),
      min_shear_rate_(min_shear_rate), max_shear_rate_(max_shear_rate) {}
//=================================================================================================//
GeneralizedNewtonianViscosity::GeneralizedNewtonianViscosity(ConstructArgs<Real, Real> args)
    : GeneralizedNewtonianViscosity(std::get<0>(args), std::get<1>(args)) {}
//=================================================================================================//
HerschelBulkleyViscosity::HerschelBulkleyViscosity(
    Real min_shear_rate, Real max_shear_rate, Real consistency_index, Real power_index, Real yield_stress)
    : GeneralizedNewtonianViscosity(min_shear_rate, max_shear_rate),
      consistency_index_(consistency_index), power_index_(power_index), yield_stress_(yield_stress) {}
//=================================================================================================//
HerschelBulkleyViscosity::HerschelBulkleyViscosity(ConstructArgs<Real, Real, Real, Real, Real> args)
    : HerschelBulkleyViscosity(std::get<0>(args), std::get<1>(args),
                               std::get<2>(args), std::get<3>(args), std::get<4>(args)) {}
//=================================================================================================//
Real HerschelBulkleyViscosity::getViscosity(Real shear_rate)
{

    Real effective_shear_rate = SMAX(SMIN(shear_rate, max_shear_rate_), min_shear_rate_);
    return (yield_stress_ + consistency_index_ * (std::pow(effective_shear_rate, power_index_))) /
           effective_shear_rate;
}
//=================================================================================================//
CarreauViscosity::CarreauViscosity(Real min_shear_rate_, Real max_shear_rate_,
                                   Real characteristic_time, Real mu_infty, Real mu0, Real power_index)
    : GeneralizedNewtonianViscosity(min_shear_rate_, max_shear_rate_),
      characteristic_time_(characteristic_time), mu_infty_(mu_infty),
      mu0_(mu0), power_index_(power_index) {}
//=================================================================================================//
Real CarreauViscosity::getViscosity(Real shear_rate)
{
    Real effective_shear_rate = SMAX(SMIN(shear_rate, max_shear_rate_), min_shear_rate_);
    return mu_infty_ + (mu0_ - mu_infty_) *
                           std::pow(1.0 + std::pow(characteristic_time_ * effective_shear_rate, 2),
                                    0.5 * (power_index_ - 1.0));
}
//=================================================================================================//
} // namespace SPH
