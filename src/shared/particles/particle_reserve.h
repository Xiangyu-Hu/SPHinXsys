/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file particle_reserve.h
 * @brief There two types of particle reserve for runtime creating and deleting particles.
 * One is buffer particles, the other ghost particles.
 * Buffer particles are saved behind real particles.
 * They may be switched from real particles or switch to real particles.
 * Ghost particles whose states are updated according to
 * boundary condition if their indices are included in the neighbor particle list.
 * The ghost particles are saved behind the buffer particles in the form of one or more ghost bounds.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_RESERVE_H
#define PARTICLE_RESERVE_H

#include "base_particles.hpp"

namespace SPH
{
struct ReserveSizeFactor
{
    Real size_factor_;
    ReserveSizeFactor(Real size_factor) : size_factor_(size_factor){};
    size_t operator()(BaseParticles &base_particles, Real particle_spacing);
};

class ParticleReserve
{
  public:
    ParticleReserve(){};
    void checkParticlesReserved();
    virtual ~ParticleReserve(){};

  protected:
    bool is_particles_reserved_ = false;
};

template <>
class Buffer<Base> : public ParticleReserve
{
  public:
    Buffer() : ParticleReserve(){};
    virtual ~Buffer(){};
    void checkEnoughBuffer(BaseParticles &base_particles);
    void allocateBufferParticles(BaseParticles &base_particles, size_t buffer_size);
};

template <class BufferSizeEstimator>
class Buffer<BufferSizeEstimator> : public Buffer<Base>
{
  public:
    template <typename... Args>
    Buffer(Args &&...args) : Buffer<Base>(), buffer_size_estimator_(std::forward<Args>(args)...){};
    virtual ~Buffer(){};
    void reserveBufferParticles(BaseParticles &base_particles, Real particle_spacing)
    {
        size_t buffer_size = buffer_size_estimator_(base_particles, particle_spacing);
        allocateBufferParticles(base_particles, buffer_size);
        is_particles_reserved_ = true;
    };

  private:
    BufferSizeEstimator buffer_size_estimator_;
};

template <>
class Ghost<Base> : public ParticleReserve
{
  public:
    Ghost() : ParticleReserve(){};
    virtual ~Ghost(){};
    size_t getGhostSize() { return ghost_size_; };
    void checkWithinGhostSize(const ParticlesBound &ghost_bound);
    IndexRange getGhostParticleRange(const ParticlesBound &ghost_bound);
    size_t allocateGhostParticles(BaseParticles &base_particles, size_t ghost_size);

  protected:
    size_t ghost_size_ = 0;
};

template <class GhostSizeEstimator>
class Ghost<GhostSizeEstimator> : public Ghost<Base>
{
  public:
    template <typename... Args>
    Ghost(Args &&...args) : Ghost<Base>(), ghost_size_estimator_(std::forward<Args>(args)...){};
    virtual ~Ghost(){};
    void reserveGhostParticles(BaseParticles &base_particles, Real particle_spacing)
    {
        ghost_size_ = ghost_size_estimator_(base_particles, particle_spacing);
        ghost_bound_.first = allocateGhostParticles(base_particles, ghost_size_);
        is_particles_reserved_ = true;
    };
    ParticlesBound &GhostBound() { return ghost_bound_; };

  private:
    GhostSizeEstimator ghost_size_estimator_;
    ParticlesBound ghost_bound_ = {0, 0};
};
} // namespace SPH
#endif // PARTICLE_RESERVE_H
