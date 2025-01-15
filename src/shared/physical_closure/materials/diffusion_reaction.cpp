#include "diffusion_reaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
BaseDiffusion::BaseDiffusion(const std::string &diffusion_species_name,
                             const std::string &gradient_species_name)
    : AbstractDiffusion(),
      diffusion_species_name_(diffusion_species_name),
      gradient_species_name_(gradient_species_name) {}
//=================================================================================================//
BaseDiffusion::BaseDiffusion(const std::string &species_name)
    : BaseDiffusion(species_name, species_name) {}
//=================================================================================================//
Real BaseDiffusion::getDiffusionTimeStepSize(Real smoothing_length)
{
    return 0.5 * smoothing_length * smoothing_length / getReferenceDiffusivity() / Real(Dimensions);
}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(const std::string &diffusion_species_name,
                                       const std::string &gradient_species_name,
                                       Real diff_cf)
    : BaseDiffusion(diffusion_species_name, gradient_species_name),
      diff_cf_(diff_cf) {}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(ConstructArgs<std::string, Real> args)
    : IsotropicDiffusion(std::get<0>(args), std::get<1>(args)) {}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(const std::string &species_name, Real diff_cf)
    : IsotropicDiffusion(species_name, species_name, diff_cf) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &diffusion_species_name,
                                                 const std::string &gradient_species_name,
                                                 Real diff_background, Real diff_max)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_background),
      diff_max_(diff_max), local_diffusivity_(nullptr) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &species_name,
                                                 Real diff_background, Real diff_max)
    : LocalIsotropicDiffusion(species_name, species_name, diff_background, diff_max) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args)
    : LocalIsotropicDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {}
//=================================================================================================//
void LocalIsotropicDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    local_diffusivity_ = base_particles->registerStateVariable<Real>(
        "ThermalConductivity", [&](size_t i) -> Real
        { return diff_cf_; });
    base_particles->addVariableToWrite<Real>("ThermalConductivity");
}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(const std::string &diffusion_species_name,
                                           const std::string &gradient_species_name,
                                           Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_cf),
      bias_direction_(bias_direction), bias_diff_cf_(bias_diff_cf),
      transformed_diffusivity_(Matd::Identity())
{
    initializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(const std::string &species_name,
                                           Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : DirectionalDiffusion(species_name, species_name, diff_cf, bias_diff_cf, bias_direction) {}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(ConstructArgs<std::string, Real, Real, Vecd> args)
    : DirectionalDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args), std::get<3>(args)) {}
//=================================================================================================//
void DirectionalDiffusion::initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
{
    bias_diff_cf_ = bias_diff_cf;
    bias_direction_ = bias_direction;
    Matd diff_i = diff_cf_ * Matd::Identity() + bias_diff_cf_ * bias_direction_ * bias_direction_.transpose();
    transformed_diffusivity_ = inverseCholeskyDecomposition(diff_i);
}
//=================================================================================================//
LocalDirectionalDiffusion::LocalDirectionalDiffusion(const std::string &diffusion_species_name,
                                                     const std::string &gradient_species_name,
                                                     Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : DirectionalDiffusion(diffusion_species_name, gradient_species_name,
                           diff_cf, bias_diff_cf, bias_direction),
      local_bias_direction_(nullptr), local_transformed_diffusivity_(nullptr) {}
//=================================================================================================//
LocalDirectionalDiffusion::LocalDirectionalDiffusion(const std::string &species_name,
                                                     Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : LocalDirectionalDiffusion(species_name, species_name, diff_cf, bias_diff_cf, bias_direction) {}
//=================================================================================================//
void LocalDirectionalDiffusion::registerLocalParameters(BaseParticles *base_particles)
{
    local_bias_direction_ = base_particles->registerStateVariable<Vecd>("Fiber");
}
//=================================================================================================//
void LocalDirectionalDiffusion::registerLocalParametersFromReload(BaseParticles *base_particles)
{
    local_bias_direction_ = base_particles->registerStateVariableFromReload<Vecd>("Fiber");
}
//=================================================================================================//
void LocalDirectionalDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    DirectionalDiffusion::initializeLocalParameters(base_particles);
    local_transformed_diffusivity_ = base_particles->registerStateVariable<Matd>(
        "LocalTransformedDiffusivity",
        [&](size_t i) -> Matd
        {
            Matd diff_i = diff_cf_ * Matd::Identity() +
                          bias_diff_cf_ * local_bias_direction_[i] * local_bias_direction_[i].transpose();
            return inverseCholeskyDecomposition(diff_i);
        });

    std::cout << "\n Local diffusion parameters setup finished " << std::endl;
};
//=================================================================================================//
} // namespace SPH
