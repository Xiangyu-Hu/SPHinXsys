#include "particle_reserve.h"

namespace SPH
{
//=================================================================================================//
size_t ReserveSizeFactor::operator()(BaseParticles &base_particles, Real particle_spacing)
{
    return std::ceil(Real(base_particles.TotalRealParticles()) * size_factor_);
}
//=================================================================================================//
void ParticleReserve::checkParticlesReserved()
{
    if (!is_particles_reserved_)
    {
        std::cout << "\n ERROR: The reserved particles are not set yet!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        std::cout << " The reserve particles should be reserved before using!" << std::endl;
        exit(1);
    };
}
//=================================================================================================//
void ParticleBuffer<Base>::checkEnoughBuffer(BaseParticles &base_particles)
{
    if (base_particles.TotalRealParticles() >= base_particles.ParticlesBound())
    {
        std::cout << "\n ERROR: Not enough buffer particles have been reserved!" << std::endl;
        std::cout << "\n You may need to increase the particle reserve." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
void ParticleBuffer<Base>::allocateBufferParticles(BaseParticles &base_particles, size_t buffer_size)
{
    base_particles.increaseParticlesBounds(buffer_size);
}
//=================================================================================================//
void Ghost<Base>::checkWithinGhostSize(const ParticlesBound &ghost_bound)
{
    if (ghost_bound.second - ghost_bound.first > ghost_size_)
    {
        std::cout << "\n ERROR: Not enough ghost particles have been reserved!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    };
}
//=================================================================================================//
IndexRange Ghost<Base>::getGhostParticleRange(const ParticlesBound &ghost_bound)
{
    return IndexRange(ghost_bound.first, ghost_bound.second);
}
//=================================================================================================//
} // namespace SPH
