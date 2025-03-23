#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
DisposerOutflowDeletionCK::
    DisposerOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      remove_real_particle_method_(particles_),
      rho0_(particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")) {}
//=================================================================================================//
TagBufferParticlesCK::TagBufferParticlesCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_buffer_particle_indicator_(
        particles_->registerStateVariableOnly<int>("BufferParticleIndicator"))
{
    particles_->addEvolvingVariable<int>("BufferParticleIndicator");
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
