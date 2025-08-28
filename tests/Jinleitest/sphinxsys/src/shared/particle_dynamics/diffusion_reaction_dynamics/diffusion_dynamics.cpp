#include "diffusion_dynamics.h"

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
    : GetDiffusionTimeStepSize(
          sph_body, DynamicCast<AbstractDiffusion>(this, &sph_body.getBaseMaterial())) {}
//=================================================================================================//
} // namespace SPH
