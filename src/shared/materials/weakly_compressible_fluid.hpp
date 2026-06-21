#ifndef WEAKLY_COMPRESSIBLE_FLUID_HPP
#define WEAKLY_COMPRESSIBLE_FLUID_HPP

#include "weakly_compressible_fluid.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EnclosureType>
WeaklyCompressibleFluid::EosKernel::EosKernel(
    const ExecutionPolicy &ex_policy, EnclosureType &encloser)
    : rho0_(encloser.rho0_), c0_(encloser.c0_), p0_(encloser.p0_){};
//=================================================================================================//
inline Real WeaklyCompressibleFluid::EosKernel::PressureFromDensity(UnsignedInt, Real rho)
{
    return p0_ * (rho / rho0_ - 1.0);
}
//=================================================================================================//
inline Real WeaklyCompressibleFluid::EosKernel::DensityFromPressure(UnsignedInt, Real p)
{
    return rho0_ * (p / p0_ + 1.0);
}
//=================================================================================================//
inline Real WeaklyCompressibleFluid::EosKernel::getSoundSpeed(UnsignedInt, Real, Real)
{
    return c0_;
}
//=================================================================================================//
inline Real WeaklyCompressibleFluid::EosKernel::getReferenceDensity(UnsignedInt)
{
    return rho0_;
}
//=================================================================================================//
template <class ExecutionPolicy, class EnclosureType>
WeaklyCompressibleMixture::EosKernel::EosKernel(
    const ExecutionPolicy &ex_policy, EnclosureType &encloser)
    : c0_(encloser.c0_), c0_sq_(c0_ * c0_),
      rho0_(encloser.dv_rho0_->DelegatedDataView(ex_policy)){};
//=================================================================================================//
inline Real WeaklyCompressibleMixture::EosKernel::PressureFromDensity(UnsignedInt index_i, Real rho)
{
    return c0_sq_ * (rho - rho0_[index_i]);
}
//=================================================================================================//
inline Real WeaklyCompressibleMixture::EosKernel::DensityFromPressure(UnsignedInt index_i, Real p)
{
    return rho0_[index_i] + p / c0_sq_;
}
//=================================================================================================//
inline Real WeaklyCompressibleMixture::EosKernel::getSoundSpeed(UnsignedInt, Real, Real)
{
    return c0_;
}
//=================================================================================================//
inline Real WeaklyCompressibleMixture::EosKernel::getReferenceDensity(UnsignedInt index_i)
{
    return rho0_[index_i];
}
//=================================================================================================//
template <class ExecutionPolicy, class EnclosureType>
WeaklyCompressibleMultiSpecies::EosKernel::EosKernel(
    const ExecutionPolicy &ex_policy, EnclosureType &encloser)
    : WeaklyCompressibleMixture::EosKernel(ex_policy, encloser),
      inv_rho0_list_(encloser.ca_inv_rho0_list_->DelegatedConstantView(ex_policy)),
      Y_list_(encloser.dv_Y_list_->DelegatedMultiEntryView(ex_policy)),
      rho0_(encloser.dv_rho0_->DelegatedDataView(ex_policy)){};
//=================================================================================================//
template <typename FractionType>
void WeaklyCompressibleMultiSpecies::EosKernel::setMassFractions(
    UnsignedInt index_i, const FractionType &mass_fractions)
{
    for (size_t k = 0; k != Y_list_.Width(); ++k)
    {
        Y_list_[index_i][k] = mass_fractions[k];
    }
}
//=================================================================================================//
inline Real WeaklyCompressibleMultiSpecies::EosKernel::computeReferenceDensity(UnsignedInt index_i)
{
    Real sum = 0.0;
    for (size_t k = 0; k != Y_list_.Width(); ++k)
    {
        sum += Y_list_[index_i][k] * inv_rho0_list_[k];
    }
    return 1.0 / sum;
}
//=================================================================================================//
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_HPP