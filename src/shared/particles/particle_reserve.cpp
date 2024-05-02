#include "particle_reserve.h"

namespace SPH
{
//=================================================================================================//
size_t ReserveSizeFactor::operator()(BaseParticles &base_particles, Real particle_spacing)
{
    return std::ceil(Real(base_particles.total_real_particles_) * size_factor_);
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
    if (base_particles.total_real_particles_ >= base_particles.real_particles_bound_)
    {
        std::cout << "\n ERROR: Not enough buffer particles have been reserved!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
void ParticleBuffer<Base>::allocateBufferParticles(BaseParticles &base_particles, size_t buffer_size)
{
    size_t old_bound = base_particles.real_particles_bound_;
    base_particles.increaseAllParticlesBounds(buffer_size);
    size_t new_bound = base_particles.real_particles_bound_;

    base_particles.resize_particles_(new_bound);
    for (size_t i = old_bound; i != new_bound; ++i)
    {
        base_particles.unsorted_id_.push_back(i);
    }
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
size_t Ghost<Base>::allocateGhostParticles(BaseParticles &base_particles, size_t ghost_size)
{
    size_t ghost_lower_bound = base_particles.particles_bound_;
    base_particles.particles_bound_ += ghost_size;
    base_particles.resize_particles_(base_particles.particles_bound_);
    base_particles.unsorted_id_.resize(base_particles.particles_bound_, 0);
    return ghost_lower_bound;
}
//=================================================================================================//
} // namespace SPH
