#include "weakly_compressible_fluid.h"

#include "base_particles.hpp"

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
WeaklyCompressibleMixture::WeaklyCompressibleMixture(
    StdVec<std::pair<std::string, Real>> species_data, Real c0) : Fluid(), c0_(c0)
{
    material_type_name_ = "WeaklyCompressibleMixture";
    species_name_list_.reserve(species_data.size());
    rho0_list_.reserve(species_data.size());
    for (const auto &pair : species_data)
    {
        species_name_list_.push_back(pair.first);
        rho0_list_.push_back(pair.second);
    }
    ca_inv_rho0_list_ = unique_entity_ptrs_.createPtr<ConstantArray<Real>>(
        rho0_list_.size(), [&](size_t k)
        { return 1.0 / rho0_list_[k]; });
}
//=================================================================================================//
WeaklyCompressibleMixture::~WeaklyCompressibleMixture() = default;
//=================================================================================================//
void WeaklyCompressibleMixture::initializeLocalParameters(BaseParticles *base_particles)
{
    Fluid::initializeLocalParameters(base_particles);
    dv_rho0_ = base_particles->registerStateVariable<Real>("ReferenceDensity", rho0_list_[0]);
    dv_Y_list_ = base_particles->registerStateVariable<Real>("MassFraction", species_name_list_);
    base_particles->addEvolvingVariable<Real>(dv_rho0_);
    base_particles->addEvolvingVariable<Real>(dv_Y_list_);
    base_particles->addVariableToWrite<Real>(dv_Y_list_);
}
//=================================================================================================//
Real WeaklyCompressibleMixture::getPressure(Real rho)
{
    return c0_ * c0_ * (rho - ReferenceDensity());
}
//=================================================================================================//
Real WeaklyCompressibleMixture::DensityFromPressure(Real p)
{
    return ReferenceDensity() + p / (c0_ * c0_);
}
//=================================================================================================//
Real WeaklyCompressibleMixture::getSoundSpeed(Real p, Real rho)
{
    return c0_;
}
//=================================================================================================//
} // namespace SPH
