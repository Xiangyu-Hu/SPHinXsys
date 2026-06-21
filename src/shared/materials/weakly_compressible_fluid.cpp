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
WeaklyCompressibleMixture::WeaklyCompressibleMixture(Real c0)
    : Fluid(), c0_(c0), dv_rho0_(nullptr)
{
    material_type_name_ = "WeaklyCompressibleMixture";
}
//=================================================================================================//
WeaklyCompressibleMixture::~WeaklyCompressibleMixture() = default;
//=================================================================================================//
WeaklyCompressibleMultiSpecies::WeaklyCompressibleMultiSpecies(
    StdVec<std::pair<std::string, Real>> species_data, Real c0)
    : WeaklyCompressibleMixture(c0)
{
    material_type_name_ = "WeaklyCompressibleMultiSpecies";
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
WeaklyCompressibleMultiSpecies::~WeaklyCompressibleMultiSpecies() = default;
//=================================================================================================//
void WeaklyCompressibleMultiSpecies::initializeLocalParameters(BaseParticles *base_particles)
{
    WeaklyCompressibleMixture::initializeLocalParameters(base_particles);
    dv_rho0_ = base_particles->registerStateVariable<Real>("ReferenceDensity", rho0_list_[0]);
    dv_Y_list_ = base_particles->registerStateVariable<Real>("MassFraction", species_name_list_);
    dv_Y_list_->fill([&](UnsignedInt index) // by default first species only
                     { return Real(1); },
                     0, base_particles->TotalRealParticles());
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
WeaklyCompressibleMultiPhase::WeaklyCompressibleMultiPhase(
    std::pair<std::string, Real> first_phase, Real c0)
    : WeaklyCompressibleMixture(c0)
{
    material_type_name_ = "WeaklyCompressibleMultiPhase";
    addPurePhase(first_phase);
}
//=================================================================================================//
WeaklyCompressibleMultiPhase::~WeaklyCompressibleMultiPhase() = default;
//=================================================================================================//
Real WeaklyCompressibleMultiPhase::ReferenceDensity() const
{
    return pure_phase_list_[0]->ReferenceDensity();
}
//=================================================================================================//
void WeaklyCompressibleMultiPhase::addPurePhase(std::pair<std::string, Real> pure_phase)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addPurePhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }
    phase_name_list_.push_back(pure_phase.first);
    pure_phase_list_.push_back(
        fluid_ptrs_.createPtr<WeaklyCompressibleFluid>(pure_phase.second, c0_));
}
//=================================================================================================//
void WeaklyCompressibleMultiPhase::addMultiSpeciesPhase(
    StdVec<std::pair<std::string, Real>> species_data)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addMultiSpeciesPhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }
    phase_name_list_.push_back("MultiSpeciesPhase" + std::to_string(phase_name_list_.size()));
    multi_species_phase_list_.push_back(
        fluid_ptrs_.createPtr<WeaklyCompressibleMultiSpecies>(species_data, c0_));
}
//=================================================================================================//
void WeaklyCompressibleMultiPhase::initializeLocalParameters(BaseParticles *base_particles)
{
    if (!is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::initializeLocalParameters :"
                  << " Phases have not been set, cannot initialize local parameters ! \n ";
        exit(1);
    }
    Fluid::initializeLocalParameters(base_particles);
    dv_rho0_ = base_particles->registerStateVariable<Real>("ReferenceDensity");
    dv_phi_list_ = base_particles->registerStateVariable<Real>("VolumeFraction", phase_name_list_);
    dv_velocity_list_ = base_particles->registerStateVariable<Vecd>("Velocity", phase_name_list_);
    base_particles->addEvolvingVariable<Real>(dv_phi_list_);
    base_particles->addEvolvingVariable<Real>(dv_rho0_);
    base_particles->addEvolvingVariable<Vecd>(dv_velocity_list_);
}
//=================================================================================================//
} // namespace SPH
