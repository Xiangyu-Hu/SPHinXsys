#include "weakly_compressible_fluid.hpp"

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
    const NamesAndDensities &species_data, Real c0)
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
WeaklyCompressibleMultiPhase::WeaklyCompressibleMultiPhase(Real c0)
    : WeaklyCompressibleMixture(c0)
{
    material_type_name_ = "WeaklyCompressibleMultiPhase";
}
//=================================================================================================//
WeaklyCompressibleMultiPhase::~WeaklyCompressibleMultiPhase() = default;
//=================================================================================================//
Real WeaklyCompressibleMultiPhase::ReferenceDensity() const
{
    return pure_phase_list_[0]->ReferenceDensity();
}
//=================================================================================================//
void WeaklyCompressibleMultiPhase::addPurePhases(const NamesAndDensities &pure_phases)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addPurePhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }

    for (const auto &pure_phase : pure_phases)
    {
        phase_name_list_.push_back(pure_phase.first);
        pure_phase_list_.push_back(
            fluid_ptrs_.createPtr<WeaklyCompressibleFluid>(pure_phase.second, c0_));
    }
}
//=================================================================================================//
void WeaklyCompressibleMultiPhase::addMultiSpeciesPhases(
    const StdVec<std::pair<std::string, NamesAndDensities>> &multi_species_phases)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addMultiSpeciesPhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }

    for (const auto &multi_species_phase : multi_species_phases)
    {
        phase_name_list_.push_back(multi_species_phase.first);
        multi_species_phase_list_.push_back(
            fluid_ptrs_.createPtr<WeaklyCompressibleMultiSpecies>(multi_species_phase.second, c0_));
    }
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
    dv_rho0_ = base_particles->registerStateVariable<Real>("ReferenceDensity", ReferenceDensity());
    dv_phi_list_ = base_particles->registerStateVariable<Real>("VolumeFraction", phase_name_list_);
    dv_phi_list_->fill([&](UnsignedInt index) // by default first species only
                       { return Real(1); },
                       0, base_particles->TotalRealParticles());
    dv_velocity_list_ = base_particles->registerStateVariable<Vecd>("Velocity", phase_name_list_);
    base_particles->addEvolvingVariable<Real>(dv_phi_list_);
    base_particles->addEvolvingVariable<Real>(dv_rho0_);
    base_particles->addEvolvingVariable<Vecd>(dv_velocity_list_);

    for (size_t k = 0; k != pure_phase_list_.size(); ++k)
    {
        pure_phase_list_[k]->initializeLocalParameters(base_particles);
    }

    for (size_t k = 0; k != multi_species_phase_list_.size(); ++k)
    {
        multi_species_phase_list_[k]->initializeLocalParameters(base_particles);
    }

    pure_eos_kernels_ = unique_entity_ptrs_.createPtr<
        ComputingKernelArray<WeaklyCompressibleFluid, PureEosKernel>>(pure_phase_list_);
    multi_species_eos_kernels_ = unique_entity_ptrs_.createPtr<
        ComputingKernelArray<WeaklyCompressibleMultiSpecies, MultiSpeciesEosKernel>>(
        multi_species_phase_list_);
}
//=================================================================================================//
} // namespace SPH
