

#include "beam_particle_generator.h"
#include "base_body.h"
#include "base_particles.h"
#include "io_all.h"
namespace SPH
{
//=================================================================================================//

LineParticleGenerator::LineParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator(sph_body),
      n_(*base_particles_.getVariableByName<Vecd>("NormalDirection")),
      thickness_(*base_particles_.getVariableByName<Real>("Thickness")),
      b_n_(*base_particles_.getVariableByName<Vecd>("BinormalDirection")),
      width_(*base_particles_.getVariableByName<Real>("Width")) {}
//=================================================================================================//
void LineParticleGenerator::initializeLineProperties(
    const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width)
{
    n_.push_back(line_normal);
    thickness_.push_back(thickness);
    b_n_.push_back(line_binormal);
    width_.push_back(width);
}
//=================================================================================================//

} // namespace SPH
