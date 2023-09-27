#include "surface_tension.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void ColorFunctionGradientInterpolationInner::interaction(size_t index_i, Real dt)
{
    Vecd grad = Vecd::Zero();
    Real weight(0);
    Real total_weight(0);
    if (indicator_[index_i] == 1 && pos_div_[index_i] > threshold_by_dimensions_)
    {
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (indicator_[index_j] == 1 && pos_div_[index_j] < threshold_by_dimensions_)
            {
                weight = inner_neighborhood.W_ij_[n] * Vol_[index_j];
                grad += weight * color_grad_[index_j];
                total_weight += weight;
            }
        }
        Vecd grad_norm = grad / (total_weight + TinyReal);
        color_grad_[index_i] = grad_norm;
        surface_norm_[index_i] = grad_norm / (grad_norm.norm() + TinyReal);
    }
}
//=================================================================================================//
ColorFunctionGradientInterpolationInner::ColorFunctionGradientInterpolationInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation), Vol_(particles_->Vol_),
      indicator_(*particles_->getVariableByName<int>("Indicator")),
      color_grad_(*particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")),
      pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
      threshold_by_dimensions_((0.75 * (Real)Dimensions))

{
    particles_->addVariableToWrite<Vecd>("SurfaceNormal");
    particles_->addVariableToWrite<Vecd>("ColorGradient");
}
//=================================================================================================//
SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation, Real gamma)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      gamma_(gamma), Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_), indicator_(*particles_->getVariableByName<int>("Indicator")),
      color_grad_(*particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")) {}
//=================================================================================================//
SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation)
    : SurfaceTensionAccelerationInner(inner_relation, 1.0) {}
//=================================================================================================//
void ColorFunctionGradientContact::interaction(size_t index_i, Real dt)
{
    Vecd gradient = Vecd::Zero();
    if (pos_div_[index_i] < threshold_by_dimensions_)
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                gradient -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }
    }
    color_grad_[index_i] += gradient;
    surface_norm_[index_i] = color_grad_[index_i] / (color_grad_[index_i].norm() + TinyReal);
}
//=================================================================================================//
SurfaceNormWithWall::SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle)
    : LocalDynamics(contact_relation.getSPHBody()), FSIContactData(contact_relation),
      contact_angle_(contact_angle),
      indicator_(*particles_->getVariableByName<int>("Indicator")),
      surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")),
      pos_div_(*particles_->getVariableByName<Real>("PositionDivergence"))
{
    particle_spacing_ = contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    smoothing_length_ = contact_relation.getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_n_.push_back(&(contact_particles_[k]->n_));
    }
}
//=================================================================================================//
void SurfaceNormWithWall::interaction(size_t index_i, Real dt)
{
    Real large_dist(1.0e6);
    Vecd n_i = surface_norm_[index_i];
    Real smoothing_factor(1.0);
    Vecd smooth_norm = Vecd::Zero();
    Vecd n_i_w = Vecd::Zero();
    /** Contact interaction. */
    if (indicator_[index_i] == 1)
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
            Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                size_t index_j = wall_neighborhood.j_[n];
                if (wall_neighborhood.r_ij_[n] < large_dist)
                {
                    Vecd n_w_t = n_i - n_i.dot(n_k[index_j]) * n_k[index_j];
                    Vecd n_t = n_w_t / (n_w_t.norm() + TinyReal);
                    n_i_w = n_t * sin(contact_angle_) + cos(contact_angle_) * n_k[index_j];
                    /** No change for multi-resolution. */
                    Real r_ij = wall_neighborhood.r_ij_[n] * n_k[index_j].dot(wall_neighborhood.e_ij_[n]);
                    if (r_ij <= smoothing_length_)
                    {
                        smoothing_factor = 0.0;
                    }
                    else
                    {
                        smoothing_factor = (r_ij - smoothing_length_) / smoothing_length_;
                    }
                    large_dist = wall_neighborhood.r_ij_[n];
                    smooth_norm = smoothing_factor * n_i + (1.0 - smoothing_factor) * n_i_w;
                    surface_norm_[index_i] = smooth_norm / (smooth_norm.norm() + TinyReal);
                }
            }
        }
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
