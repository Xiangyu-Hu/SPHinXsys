#include "general_solid_dynamics.h"
#include "sph_system.hpp"
namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
DistributingPointForces::
    DistributingPointForces(SPHBody &sph_body, std::vector<Vecd> point_forces,
                            std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                            Real particle_spacing_ref, Real h_spacing_ratio)
    : LocalDynamics(sph_body),
      point_forces_(point_forces), reference_positions_(reference_positions),
      time_to_full_external_force_(time_to_full_external_force),
      particle_spacing_ref_(particle_spacing_ref), h_spacing_ratio_(h_spacing_ratio),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      thickness_(particles_->getVariableDataByName<Real>("Thickness")),
      physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
{
    weight_.resize(point_forces_.size());
    for (size_t i = 0; i < point_forces_.size(); i++)
    {
        time_dependent_point_forces_.push_back(Vecd::Zero());
        sum_of_weight_.push_back(0.0);
        weight_[i] = particles_->registerStateVariable<Real>("Weight_" + std::to_string(i));
    }

    getWeight(); // TODO: should be revised and parallelized, using SimpleDynamics
}
//=================================================================================================//
void DistributingPointForces::getWeight()
{
    Kernel *kernel_ = sph_body_.getSPHAdaptation().getKernel();
    Real reference_smoothing_length = sph_body_.getSPHAdaptation().ReferenceSmoothingLength();
    Real smoothing_length = h_spacing_ratio_ * particle_spacing_ref_;
    Real h_ratio = reference_smoothing_length / smoothing_length;
    Real cutoff_radius_sqr = pow(2.0 * smoothing_length, 2);
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        sum_of_weight_[i] = 0.0;
        for (size_t index = 0; index < particles_->TotalRealParticles(); ++index)
        {
            weight_[i][index] = 0.0;
            Vecd displacement = reference_positions_[i] - pos_[index];
            if (displacement.squaredNorm() <= cutoff_radius_sqr)
            {
                weight_[i][index] = kernel_->W(h_ratio, displacement.norm(), displacement);
                sum_of_weight_[i] += weight_[i][index];
            }
        }
    }
}
//=================================================================================================//
void DistributingPointForces::setupDynamics(Real dt)
{
    Real current_time = *physical_time_;
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        time_dependent_point_forces_[i] = current_time < time_to_full_external_force_
                                              ? current_time * point_forces_[i] / time_to_full_external_force_
                                              : point_forces_[i];
    }
}
//=================================================================================================//
void DistributingPointForces::update(size_t index_i, Real dt)
{
    force_prior_[index_i] = Vecd::Zero();
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        Vecd force = weight_[i][index_i] / (sum_of_weight_[i] + TinyReal) * time_dependent_point_forces_[i];
        force_prior_[index_i] += force;
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH