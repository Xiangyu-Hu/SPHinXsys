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
    StdVec<std::string> species_name_list, StdVec<Real> rho0_list, Real c0)
    : Fluid(), species_name_list_(species_name_list), rho0_list_(rho0_list), c0_(c0)
{
    material_type_name_ = "WeaklyCompressibleMixture";
    ca_inv_rho0_list_ = unique_entity_ptrs_.createPtr<ConstantArray<Real>>(
        rho0_list_.size(), [&](size_t k)
        { return 1.0 / rho0_list_[k]; });
}
//=================================================================================================//
WeaklyCompressibleMixture::~WeaklyCompressibleMixture() = default;
//=================================================================================================//
void WeaklyCompressibleMixture::initializeLocalParameters(BaseParticles *base_particles)
{
    dv_rho0_ = base_particles->registerStateVariable<Real>("ReferenceDensity", rho0_list_[0]);
    dv_Y_list_ = base_particles->registerStateVariable<Real>("MassFraction", species_name_list_);
}
//=================================================================================================//
} // namespace SPH
