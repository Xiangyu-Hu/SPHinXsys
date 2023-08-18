#include "general_life_time_dynamics.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
BaseLifeTimeDynamics::BaseLifeTimeDynamics(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      particle_split_merge_(DynamicCast<ParticleSplitAndMerge>(this, *sph_body.sph_adaptation_)),
      inv_rho0_(1.0 / sph_body_.base_material_->ReferenceDensity()),
      rho_(particles_->rho_), pos_(particles_->pos_), Vol_(particles_->Vol_),
      mass_(particles_->mass_),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
RefinementInPrescribedRegion::
    RefinementInPrescribedRegion(SPHBody &sph_body, size_t body_buffer_width, Shape &refinement_region)
    : BaseSplitDynamics<Vecd>(sph_body, body_buffer_width),
      refinement_region_bounds_(refinement_region.getBounds()),
      normal_distribution_(0, 1) {}
//=================================================================================================//
void RefinementInPrescribedRegion::setupDynamics(Real dt)
{
    random_seed_ = std::mt19937(random_device_());
}
//=================================================================================================//
void RefinementInPrescribedRegion::update(size_t index_i, Real dt)
{
    if (checkSplit(index_i))
    {
        Vecd split_shift = execFirstSplit(index_i);

        mutex_split_.lock();
        execOtherSplit(index_i, split_shift);
        mutex_split_.unlock();
    }
}
//=================================================================================================//
bool RefinementInPrescribedRegion::checkSplit(size_t index_i)
{
    Real non_deformed_volume = mass_[index_i] * inv_rho0_;
    bool is_split_allowed = particle_split_merge_.isSplitAllowed(non_deformed_volume);
    bool is_split_inside = checkLocation(refinement_region_bounds_, pos_[index_i], non_deformed_volume);

    return (is_split_allowed && is_split_inside) ? true : false;
}
//=================================================================================================//
bool RefinementInPrescribedRegion::
    checkLocation(const BoundingBox &refinement_region_bounds, Vecd position, Real volume)
{
    int bound_number = 0;
    for (int axis_direction = 0; axis_direction != Dimensions; ++axis_direction)
    {
        Real particle_spacing = pow(volume, 1.0 / (Real)Dimensions);
        if (position[axis_direction] > (refinement_region_bounds.first_[axis_direction] + particle_spacing) &&
            position[axis_direction] < (refinement_region_bounds.second_[axis_direction] - particle_spacing))
            bound_number += 1;
    }
    return bound_number != Dimensions ? false : true;
}
//=================================================================================================//
Vecd RefinementInPrescribedRegion::execFirstSplit(size_t index_i)
{
    mass_[index_i] *= 0.5;
    Real split_volume = Vol_[index_i] * 0.5;
    Vol_[index_i] = split_volume;
    Real split_spacing = pow(split_volume, 1.0 / (Real)Dimensions);
    h_ratio_[index_i] = particle_split_merge_.ReferenceSpacing() / split_spacing;

    Vecd shift = Vecd::Zero();
    for (int k = 0; k < Dimensions; ++k)
    {
        shift[k] = normal_distribution_(random_seed_);
    }

    return 0.5 * split_spacing * shift / (shift.norm() + TinyReal);
}
//=================================================================================================//
void RefinementInPrescribedRegion::execOtherSplit(size_t index_i, const Vecd &split_shift)
{
    size_t total_real_particles = particles_->total_real_particles_;

    if (total_real_particles >= particles_->real_particles_bound_)
    {
        std::cout << "ParticleSplitWithPrescribedArea: \n"
                  << "Not enough body buffer particles! Exit the code."
                  << "\n";
        exit(0);
    }
    else
    {
        particles_->copyFromAnotherParticle(total_real_particles, index_i);

        pos_[index_i] += split_shift;
        pos_[total_real_particles] -= split_shift;

        particles_->total_real_particles_ += 1;
    }
}
//=================================================================================================//
} // namespace SPH
