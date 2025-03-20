#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <typename AlignedBoxPartType>
EmitterInflowInjectionCK<AlignedBoxPartType>::
    EmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      buffer_(buffer), sv_aligned_box_(aligned_box_part.svAlignedBox()),
      create_real_particle_method_(this->particles_),
      rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure"))
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <typename AlignedBoxPartType>
EmitterInflowInjectionCK<AlignedBoxPartType>::FinishDynamics::
    FinishDynamics(EmitterInflowInjectionCK<AlignedBoxPartType> &encloser)
    : particles_(encloser.particles_), buffer_(encloser.buffer_) {}
//=================================================================================================//
template <typename AlignedBoxPartType>
void EmitterInflowInjectionCK<AlignedBoxPartType>::FinishDynamics::operator()()
{
    buffer_.checkEnoughBuffer(*particles_);
}
//=================================================================================================//
DisposerOutflowDeletionCK::
    DisposerOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(particles_->svTotalRealParticles()),
      remove_real_particle_method_(particles_),
      dv_buffer_particle_indicator_(
          particles_->registerStateVariableOnly<int>("BufferParticleIndicator")),
      rho0_(particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")) {}
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
} // namespace fluid_dynamics
} // namespace SPH
