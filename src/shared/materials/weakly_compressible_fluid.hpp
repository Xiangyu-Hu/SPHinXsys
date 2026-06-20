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
template <class FirstFluidType, typename... Args>
WeaklyCompressibleMultiPhase::WeaklyCompressibleMultiPhase(Real c0, Args &&...args)
    : WeaklyCompressibleMixture(c0)
{
    addPhase<FirstFluidType>(std::forward<Args>(args)...);
}
//=================================================================================================//
template <class FluidType, typename... Args>
void WeaklyCompressibleMultiPhase::addPhase(Args &&...args)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addPhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }
    auto &fluid = fluid_ptrs_.createPtr<FluidType>(std::forward<Args>(args)...);
    phase_list_.push_back(fluid);
    phase_name_list_.push_back(fluid->Name());
}
//=================================================================================================//
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_HPP