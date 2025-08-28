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
} // namespace fluid_dynamics
} // namespace SPH
