#include "bidirectional_boundary_ck.h"

#include "adaptation.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BufferIndicationCK::BufferIndicationCK(OrientedBoxByCell &oriented_box_part)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      part_id_(oriented_box_part.getPartID()),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_buffer_indicator_(particles_->registerStateVariable<int>("BufferIndicator"))
{
    particles_->addEvolvingVariable<int>("BufferIndicator");
}
//=================================================================================================//
ResetBufferCorrectionMatrixCK::ResetBufferCorrectionMatrixCK(OrientedBoxByCell &oriented_box_part)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_B_(particles_->registerStateVariable<Matd>(
          "LinearCorrectionMatrix", IdentityMatrix<Matd>::value)),
      radius_(sph_body_->getSPHAdaptation().getKernel()->CutOffRadius()) {}
//=================================================================================================//
BufferOutflowIndication::BufferOutflowIndication(OrientedBoxByCell &oriented_box_part)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      part_id_(oriented_box_part.getPartID()),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      sv_total_real_particles_(particles_->svTotalRealParticles()),
      dv_buffer_indicator_(particles_->registerStateVariable<int>("BufferIndicator")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_life_status_(particles_->registerStateVariable<int>("LifeStatus", 0)) {}
//=================================================================================================//
BufferOutflowIndication::UpdateKernel::
    IsDeletable::IsDeletable(int part_id, OrientedBox *oriented_box,
                             Vecd *pos, int *buffer_particle_indicator)
    : part_id_(part_id), oriented_box_(oriented_box), pos_(pos),
      buffer_indicator_(buffer_particle_indicator) {}
//=================================================================================================//
OutflowParticleDeletion::OutflowParticleDeletion(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      remove_real_particle_method_(particles_),
      dv_life_status_(particles_->getVariableByName<int>("LifeStatus")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
