#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
EmitterInflowInjectionCK::
    EmitterInflowInjectionCK(AlignedBoxPartByParticle &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartByParticle>(aligned_box_part),
      buffer_(buffer), sv_aligned_box_(aligned_box_part.svAlignedBox()),
      create_real_particle_method_(particles_),
      rho0_(particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure"))
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
EmitterInflowInjectionCK::FinishDynamics::
    FinishDynamics(EmitterInflowInjectionCK &encloser)
    : particles_(encloser.particles_), buffer_(encloser.buffer_) {}
//=================================================================================================//
void EmitterInflowInjectionCK::FinishDynamics::operator()()
{
    buffer_.checkEnoughBuffer(*particles_);
}
//=================================================================================================//
DisposerOutflowDeletionCK::
    DisposerOutflowDeletionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      remove_real_particle_method_(particles_),
      rho0_(particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure"))
{
}
//=================================================================================================//
DisposerOutflowDeletionCK::FinishDynamics::
    FinishDynamics(DisposerOutflowDeletionCK &encloser)
    : particles_(encloser.particles_) {}
//=================================================================================================//
void DisposerOutflowDeletionCK::FinishDynamics::operator()()
{
}

//=================================================================================================//
TagBufferParticlesCK::TagBufferParticlesCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_buffer_particle_indicator_(particles_->template registerStateVariableOnly<int>("BufferParticleIndicator"))
{
}
} // namespace fluid_dynamics
} // namespace SPH
