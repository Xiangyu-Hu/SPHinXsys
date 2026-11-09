/**
 * @file gpu_paper_local_isotropic_diffusion.cpp
 * @brief GPU-paper overlay material (A3). Does not modify src/.
 */
#include "gpu_paper_local_isotropic_diffusion.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
GpuPaperLocalIsotropicDiffusion::GpuPaperLocalIsotropicDiffusion(
    const std::string &species_name, const std::string &gradient_species_name, Real diff_background,
    Real diff_max, Real cv)
    : IsotropicDiffusion(species_name, gradient_species_name, diff_background, cv),
      diff_max_(diff_max), local_diffusivity_(nullptr), dv_thermal_conductivity_(nullptr) {}
//=================================================================================================//
GpuPaperLocalIsotropicDiffusion::GpuPaperLocalIsotropicDiffusion(const std::string &species_name,
                                                                 Real diff_background, Real diff_max,
                                                                 Real cv)
    : GpuPaperLocalIsotropicDiffusion(species_name, species_name, diff_background, diff_max, cv) {}
//=================================================================================================//
GpuPaperLocalIsotropicDiffusion::GpuPaperLocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args)
    : GpuPaperLocalIsotropicDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {}
//=================================================================================================//
void GpuPaperLocalIsotropicDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    dv_thermal_conductivity_ = base_particles->registerStateVariable<Real>(
        "ThermalConductivity", [&](size_t i) -> Real
        { return d_coeff_; });
    local_diffusivity_ = dv_thermal_conductivity_->Data();
    base_particles->addVariableToWrite<Real>("ThermalConductivity");
}
//=================================================================================================//
} // namespace SPH
