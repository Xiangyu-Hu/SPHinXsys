#include "bidirectional_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
TagBufferParticlesCK::TagBufferParticlesCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_buffer_particle_indicator_(
          particles_->registerStateVariableOnly<int>("BufferParticleIndicator"))
{
    particles_->addEvolvingVariable<int>("BufferParticleIndicator");
}
//=================================================================================================//
BufferOutflowDeletionCK::BufferOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(particles_->svTotalRealParticles()),
      remove_real_particle_method_(particles_),
      dv_buffer_particle_indicator_(
          particles_->registerStateVariableOnly<int>("BufferParticleIndicator")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")) {}
//=================================================================================================//
BufferOutflowDeletionCK::UpdateKernel::
    IsDeletable::IsDeletable(int part_id, AlignedBox *aligned_box,
                             Vecd *pos, int *buffer_particle_indicator)
    : part_id_(part_id), aligned_box_(aligned_box), pos_(pos),
      buffer_particle_indicator_(buffer_particle_indicator) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
