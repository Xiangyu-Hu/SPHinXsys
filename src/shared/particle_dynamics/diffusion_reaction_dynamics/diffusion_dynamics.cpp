#include "diffusion_dynamics.h"

namespace SPH
{
//=================================================================================================//
GetDiffusionTimeStepSize::GetDiffusionTimeStepSize(SPHBody &sph_body)
    : BaseDynamics<Real>()
{
    DiffusionType &diffusion = DynamicCast<DiffusionType>(this, sph_body.getMaterial());
    Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
    diff_time_step_ = diffusion.getDiffusionTimeStepSize(smoothing_length);
}
//=================================================================================================//
} // namespace SPH
