#include "diffusion_dynamics.h"

#include "adaptation.h"
#include "base_material.h"
#include "base_body.hpp"

namespace SPH
{
//=================================================================================================//
GetDiffusionTimeStepSize::GetDiffusionTimeStepSize(
    SPHBody &sph_body, AbstractDiffusion *abstract_diffusion)
    : BaseDynamics<Real>()
{
    Real smoothing_length = sph_body.getSPHAdaptation().ReferenceSmoothingLength();
    diff_time_step_ = abstract_diffusion->getDiffusionTimeStepSize(smoothing_length);
}
//=================================================================================================//
GetDiffusionTimeStepSize::GetDiffusionTimeStepSize(SPHBody &sph_body)
    : GetDiffusionTimeStepSize(sph_body, &sph_body.getMaterialProperty<AbstractDiffusion>()) {}
//=================================================================================================//
} // namespace SPH
