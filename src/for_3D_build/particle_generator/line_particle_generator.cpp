#include "line_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Line>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body) {}
//=================================================================================================//
void ParticleGenerator<Line>::addLineProperties(
    const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width)
{
    line_normal_.push_back(line_normal);
    line_thickness_.push_back(thickness);
    line_binormal_.push_back(line_binormal);
    line_width_.push_back(width);
}
//=================================================================================================//
void ParticleGenerator<Line>::initializeParticleVariables()
{
    ParticleGenerator<Base>::initializeParticleVariables();
    base_particles_.registerSharedVariableFrom<Vecd>("NormalDirection", line_normal_);
    base_particles_.registerSharedVariableFrom<Vecd>("BinormalDirection", line_binormal_);
    base_particles_.registerSharedVariableFrom<Real>("Thickness", line_thickness_);
    base_particles_.registerSharedVariableFrom<Real>("Width", line_width_);
}
//=================================================================================================//
} // namespace SPH
