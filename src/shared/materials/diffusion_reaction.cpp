#include "diffusion_reaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void DirectionalDiffusion::initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
{
    bias_diff_cf_ = bias_diff_cf;
    bias_direction_ = bias_direction;
    Matd diff_i = diff_cf_ * Matd::Identity() + bias_diff_cf_ * bias_direction_ * bias_direction_.transpose();
    transformed_diffusivity_ = inverseCholeskyDecomposition(diff_i);
};
//=================================================================================================//
void LocalDirectionalDiffusion::registerReloadLocalParameters(BaseParticles *base_particles)
{
    base_particles->registerVariable(local_bias_direction_, "Fiber");
    base_particles->addVariableToReload<Vecd>("Fiber");
}
//=================================================================================================//
void LocalDirectionalDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    DirectionalDiffusion::initializeLocalParameters(base_particles);
    base_particles->registerVariable(
        local_transformed_diffusivity_, "LocalTransformedDiffusivity",
        [&](size_t i) -> Matd
        {
            Matd diff_i = diff_cf_ * Matd::Identity() + bias_diff_cf_ * local_bias_direction_[i] * local_bias_direction_[i].transpose();
            return inverseCholeskyDecomposition(diff_i);
        });

    std::cout << "\n Local diffusion parameters setup finished " << std::endl;
};
//=================================================================================================//
} // namespace SPH
