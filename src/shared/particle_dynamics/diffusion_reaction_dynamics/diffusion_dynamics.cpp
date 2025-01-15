#include "diffusion_dynamics.h"

namespace SPH
{
//=================================================================================================//
GetDiffusionTimeStepSize::GetDiffusionTimeStepSize(SPHBody &sph_body)
    : BaseDynamics<Real>()
{
    AbstractDiffusion &diffusion = DynamicCast<AbstractDiffusion>(this, sph_body.getBaseMaterial());
    Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
    diff_time_step_ = diffusion.getDiffusionTimeStepSize(smoothing_length);
}
//=================================================================================================//
} // namespace SPH
