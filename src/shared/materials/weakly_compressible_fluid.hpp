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
      inv_rho0_list_(encloser.ca_inv_rho0_list_->DelegatedArrayView(ex_policy)),
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
template <class ExecutionPolicy, class EnclosureType>
WeaklyCompressibleMultiPhase::EosKernel::EosKernel(
    const ExecutionPolicy &ex_policy, EnclosureType &encloser)
    : WeaklyCompressibleMixture::EosKernel(ex_policy, encloser),
      pure_eos_(encloser.pure_eos_kernels_->DelegatedArrayView(ex_policy)),
      multi_species_eos_(encloser.multi_species_eos_kernels_->DelegatedArrayView(ex_policy)),
      phi_list_(encloser.dv_phi_list_->DelegatedMultiEntryView(ex_policy)),
      velocity_list_(encloser.dv_velocity_list_->DelegatedMultiEntryView(ex_policy)){};
//=================================================================================================//
template <typename FractionType>
void WeaklyCompressibleMultiPhase::EosKernel::setVolumeFractions(
    UnsignedInt index_i, const FractionType &volume_fractions)
{
    for (size_t k = 0; k != phi_list_.Width(); ++k)
    {
        phi_list_[index_i][k] = volume_fractions[k];
    }
}
//=================================================================================================//
inline Real WeaklyCompressibleMultiPhase::EosKernel::computeReferenceDensity(UnsignedInt index_i)
{
    Real sum = 0.0;
    UnsignedInt num_pure_phases = pure_eos_.Size();
    for (size_t k = 0; k != phi_list_.Width(); ++k)
    {
        Real reference_density =
            k < num_pure_phases
                ? pure_eos_[k].getReferenceDensity(index_i)
                : multi_species_eos_[k - num_pure_phases].getReferenceDensity(index_i);
        sum += phi_list_[index_i][k] * reference_density;
    }
    return sum;
}
//=================================================================================================//
template <typename VelocityType>
inline Vecd WeaklyCompressibleMultiPhase::EosKernel::computeMixtureVelocity(
    UnsignedInt index_i, const VelocityType &velocities)
{
    Vecd velocity_sum = Vecd::Zero();
    Real weight_sum(0);
    UnsignedInt num_pure_phases = pure_eos_.Size();
    for (size_t k = 0; k != phi_list_.Width(); ++k)
    {
        Real reference_density =
            k < num_pure_phases
                ? pure_eos_[k].getReferenceDensity(index_i)
                : multi_species_eos_[k - num_pure_phases].getReferenceDensity(index_i);
        Real weight = phi_list_[index_i][k] * reference_density;
        velocity_sum += weight * velocities[k];
        weight_sum += weight;
    }
    return velocity_sum / weight_sum;
}
//=================================================================================================//
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_HPP