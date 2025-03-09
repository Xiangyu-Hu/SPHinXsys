#include "diffusion_reaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
BaseDiffusion::BaseDiffusion(const std::string &diffusion_species_name,
                             const std::string &gradient_species_name, Real cv)
    : AbstractDiffusion(),
      diffusion_species_name_(diffusion_species_name),
      gradient_species_name_(gradient_species_name), cv_(cv) {}
//=================================================================================================//
BaseDiffusion::BaseDiffusion(const std::string &species_name, Real cv)
    : BaseDiffusion(species_name, species_name, cv) {}
//=================================================================================================//
Real BaseDiffusion::getDiffusionTimeStepSize(Real smoothing_length)
{
    return 0.5 * smoothing_length * smoothing_length / getReferenceDiffusivity() / Real(Dimensions);
}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(const std::string &diffusion_species_name,
                                       const std::string &gradient_species_name,
                                       Real d_coeff, Real cv)
    : BaseDiffusion(diffusion_species_name, gradient_species_name, cv),
      d_coeff_(d_coeff) {}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(ConstructArgs<std::string, Real> args)
    : IsotropicDiffusion(std::get<0>(args), std::get<1>(args)) {}
//=================================================================================================//
IsotropicDiffusion::IsotropicDiffusion(const std::string &species_name, Real d_coeff, Real cv)
    : IsotropicDiffusion(species_name, species_name, d_coeff, cv) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &diffusion_species_name,
                                                 const std::string &gradient_species_name,
                                                 Real diff_background, Real diff_max, Real cv)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_background, cv),
      diff_max_(diff_max), local_diffusivity_(nullptr) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(const std::string &species_name,
                                                 Real diff_background, Real diff_max, Real cv)
    : LocalIsotropicDiffusion(species_name, species_name, diff_background, diff_max, cv) {}
//=================================================================================================//
LocalIsotropicDiffusion::LocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args)
    : LocalIsotropicDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {}
//=================================================================================================//
void LocalIsotropicDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    local_diffusivity_ = base_particles->registerStateVariable<Real>(
        "ThermalConductivity", [&](size_t i) -> Real
        { return d_coeff_; });
    base_particles->addVariableToWrite<Real>("ThermalConductivity");
}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(const std::string &diffusion_species_name,
                                           const std::string &gradient_species_name,
                                           Real d_coeff, Real bias_d_coeff_,
                                           Vecd bias_direction, Real cv)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, d_coeff, cv),
      bias_direction_(bias_direction), bias_d_coeff_(bias_d_coeff_),
      transformed_diffusivity_(Matd::Identity())
{
    initializeDirectionalDiffusivity(d_coeff, bias_d_coeff_, bias_direction);
}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(const std::string &species_name,
                                           Real d_coeff, Real bias_d_coeff_,
                                           Vecd bias_direction, Real cv)
    : DirectionalDiffusion(species_name, species_name, d_coeff, bias_d_coeff_, bias_direction, cv) {}
//=================================================================================================//
DirectionalDiffusion::DirectionalDiffusion(ConstructArgs<std::string, Real, Real, Vecd> args)
    : DirectionalDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args), std::get<3>(args)) {}
//=================================================================================================//
void DirectionalDiffusion::initializeDirectionalDiffusivity(Real d_coeff, Real bias_d_coeff_, Vecd bias_direction)
{
    bias_d_coeff_ = bias_d_coeff_;
    bias_direction_ = bias_direction;
    Matd diff_i = d_coeff_ * Matd::Identity() + bias_d_coeff_ * bias_direction_ * bias_direction_.transpose();
    transformed_diffusivity_ = inverseCholeskyDecomposition(diff_i);
}
//=================================================================================================//
LocalDirectionalDiffusion::LocalDirectionalDiffusion(const std::string &diffusion_species_name,
                                                     const std::string &gradient_species_name,
                                                     Real d_coeff, Real bias_d_coeff_,
                                                     Vecd bias_direction, Real cv)
    : DirectionalDiffusion(diffusion_species_name, gradient_species_name,
                           d_coeff, bias_d_coeff_, bias_direction, cv),
      local_bias_direction_(nullptr), local_transformed_diffusivity_(nullptr) {}
//=================================================================================================//
LocalDirectionalDiffusion::LocalDirectionalDiffusion(const std::string &species_name,
                                                     Real d_coeff, Real bias_d_coeff_,
                                                     Vecd bias_direction, Real cv)
    : LocalDirectionalDiffusion(species_name, species_name, d_coeff, bias_d_coeff_, bias_direction, cv) {}
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
            Matd diff_i = d_coeff_ * Matd::Identity() +
                          bias_d_coeff_ * local_bias_direction_[i] * local_bias_direction_[i].transpose();
            return inverseCholeskyDecomposition(diff_i);
        });

    std::cout << "\n Local diffusion parameters setup finished " << std::endl;
};
//=================================================================================================//
} // namespace SPH
