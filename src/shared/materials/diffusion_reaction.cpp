#include "diffusion_reaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
BaseDiffusion::BaseDiffusion(const std::string &diffusion_species_name,
                             const std::string &gradient_species_name)
    : BaseMaterial(),
      diffusion_species_name_(diffusion_species_name),
      gradient_species_name_(gradient_species_name)
{
    material_type_name_ = "BaseDiffusion";
}
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
      diff_cf_(diff_cf)
{
    material_type_name_ = "IsotropicDiffusion";
}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(const std::string &species_name, Real diff_cf)
    : IsotropicDiffusion(species_name, species_name, diff_cf) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &diffusion_species_name,
                                                 const std::string &gradient_species_name,
                                                 Real diff_cf)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_cf),
      local_diffusivity_(nullptr)
{
    material_type_name_ = "LocalIsotropicDiffusion";
}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &species_name, Real diff_cf)
    : LocalIsotropicDiffusion(species_name, species_name, diff_cf) {}
//=================================================================================================//
void LocalIsotropicDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    local_diffusivity_ = base_particles->registerSharedVariable<Real>(
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
    material_type_name_ = "DirectionalDiffusion";
    initializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(const std::string &species_name,
                                           Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : DirectionalDiffusion(species_name, species_name, diff_cf, bias_diff_cf, bias_direction) {}
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
      local_bias_direction_(nullptr), local_transformed_diffusivity_(nullptr)
{
    material_type_name_ = "LocalDirectionalDiffusion";
}
//=================================================================================================//
LocalDirectionalDiffusion::LocalDirectionalDiffusion(const std::string &species_name,
                                                     Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
    : LocalDirectionalDiffusion(species_name, species_name, diff_cf, bias_diff_cf, bias_direction) {}
//=================================================================================================//
void LocalDirectionalDiffusion::registerLocalParameters(BaseParticles *base_particles)
{
    local_bias_direction_ = base_particles->registerSharedVariable<Vecd>("Fiber");
}
//=================================================================================================//
void LocalDirectionalDiffusion::registerLocalParametersFromReload(BaseParticles *base_particles)
{
    local_bias_direction_ = base_particles->registerSharedVariableFromReload<Vecd>("Fiber");
}
//=================================================================================================//
void LocalDirectionalDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    DirectionalDiffusion::initializeLocalParameters(base_particles);
    local_transformed_diffusivity_ = base_particles->registerSharedVariable<Matd>(
        "LocalTransformedDiffusivity",
        [&](size_t i) -> Matd
        {
            Matd diff_i = diff_cf_ * Matd::Identity() +
                          bias_diff_cf_ * (*local_bias_direction_)[i] * (*local_bias_direction_)[i].transpose();
            return inverseCholeskyDecomposition(diff_i);
        });

    std::cout << "\n Local diffusion parameters setup finished " << std::endl;
};
//=================================================================================================//
} // namespace SPH
