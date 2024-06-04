#include "line_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Line>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body),
      n_(*particles_->getVariableByName<Vecd>("NormalDirection")),
      thickness_(*particles_->getVariableByName<Real>("Thickness")),
      b_n_(*particles_->getVariableByName<Vecd>("BinormalDirection")),
      width_(*particles_->getVariableByName<Real>("Width")) {}
//=================================================================================================//
void ParticleGenerator<Line>::initializeLineProperties(
    const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width)
{
    n_.push_back(line_normal);
    thickness_.push_back(thickness);
    b_n_.push_back(line_binormal);
    width_.push_back(width);
}
//=================================================================================================//
} // namespace SPH
