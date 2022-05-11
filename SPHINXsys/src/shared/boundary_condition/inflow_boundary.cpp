
#include "inflow_boundary.h"

namespace SPH
{
    BodyRegionByParticleWithFace::BodyRegionByParticleWithFace(RealBody &real_body, SegmentFace &segment_face, Real scale)
        : BodyRegionWithFace(real_body, segment_face, scale),
          base_particles_(real_body.base_particles_)
    {
        tagParticles();
    }

    void BodyRegionByParticleWithFace::tagParticles()
    {
        for (size_t i = 0; i < base_particles_->total_real_particles_; ++i)
        {
            const auto &pos = base_particles_->pos_n_[i];
            auto d = getSignedDistance(pos);
            if (d < region_width_)
            {
                body_part_particles_.push_back(i);
            }
        }
    }

    void PartDynamicsByParticleWithFace::exec(Real dt)
    {
        setBodyUpdated();
        setupDynamics(dt);
        for (size_t i = 0; i < body_part_particles_.size(); ++i)
        {
            particle_functor_(body_part_particles_[i], dt);
        }
    }

    void PartDynamicsByParticleWithFace::parallel_exec(Real dt)
    {
        setBodyUpdated();
        setupDynamics(dt);
        parallel_for(
            blocked_range<size_t>(0, body_part_particles_.size()),
            [&](const blocked_range<size_t> &r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    particle_functor_(body_part_particles_[i], dt);
                }
            },
            ap);
    }

    PartSimpleDynamicsByParticleWithFace::PartSimpleDynamicsByParticleWithFace(RealBody &real_body, BodyRegionByParticleWithFace &body_part)
        : PartDynamicsByParticleWithFace(real_body, body_part)
    {
        particle_functor_ = std::bind(&PartSimpleDynamicsByParticleWithFace::Update, this, _1, _2);
    }

    InflowInjectingWithFace::InflowInjectingWithFace(FluidBody &fluid_body, BodyRegionByParticleWithFace &body_part, size_t body_buffer_width)
        : PartSimpleDynamicsByParticleWithFace(fluid_body, body_part),
          DataDelegateSimple<FluidBody, FluidParticles, Fluid>(fluid_body),
          body_part_(body_part),
          pos_n_(particles_->pos_n_), rho_n_(particles_->rho_n_), p_(particles_->p_),
          periodic_translation_(body_part.getRegionWidth())
    {
        size_t total_body_buffer_particles = body_part_particles_.size() * body_buffer_width;
        particles_->addBufferParticles(total_body_buffer_particles);
        sph_body_->allocateConfigurationMemoriesForBufferParticles();
    }

    void InflowInjectingWithFace::checking_bound_(size_t unsorted_index_i, Real dt)
    {
        size_t sorted_index_i = sorted_id_[unsorted_index_i];
        if (body_part_.getSignedDistance(pos_n_[sorted_index_i]) > periodic_translation_)
        {
            if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
            {
                std::cout << "InflowInjectingWithFace::checking_bound_: \n"
                          << "Not enough body buffer particles! Exit the code."
                          << "\n";
                exit(0);
            }
            /** Buffer Particle state copied from real particle. */
            particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
            /** Realize the buffer particle by increasing the number of real particle in the body.  */
            particles_->total_real_particles_ += 1;
            /** Periodic bounding. */
            pos_n_[sorted_index_i] -= periodic_translation_ * body_part_.getDirectionToFluid();
        }
    }

    InflowConditionWithFace::InflowConditionWithFace(FluidBody &fluid_body, BodyRegionByCellsWithFace &body_part)
        : PartSimpleDynamicsByCellsWithFace(fluid_body, body_part),
          DataDelegateSimple<FluidBody, FluidParticles, Fluid>(fluid_body),
          body_part_(body_part),
          vel_n_(particles_->vel_n_),
          pos_n_(particles_->pos_n_)
    {
    }

    VelocityInflowConditionWithFace::VelocityInflowConditionWithFace(FluidBody &fluid_body, BodyRegionByCellsWithFace &body_part)
        : InflowConditionWithFace(fluid_body, body_part)
    {
    }

    Vecd VelocityInflowConditionWithFace::getTargetVelocity(Vecd &position, Vecd &velocity)
    {
        if (body_part_.getSignedDistance(position) < body_part_.getRegionWidth())
        {
            return defineVelocityProfile(position, velocity);
        }

        return velocity;
    }

    void VelocityInflowConditionWithFace::Update(size_t index_i, Real dt)
    {
        vel_n_[index_i] = getTargetVelocity(pos_n_[index_i], vel_n_[index_i]);
    }

} // namespace SPH
