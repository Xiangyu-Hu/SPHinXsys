#include "fluid_surface_inner.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FreeSurfaceIndicationInner::
    FreeSurfaceIndicationInner(BaseInnerRelation &inner_relation)
    : FreeSurfaceIndication<FluidDataInner>(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
void FreeSurfaceIndicationInner::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        pos_div -= inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.r_ij_[n];
    }
    pos_div_[index_i] = pos_div;
}
//=================================================================================================//
void FreeSurfaceIndicationInner::update(size_t index_i, Real dt)
{
    indicator_[index_i] = 1;
    if (pos_div_[index_i] > threshold_by_dimensions_ && !isVeryNearFreeSurface(index_i))
        indicator_[index_i] = 0;
}
//=================================================================================================//
bool FreeSurfaceIndicationInner::isVeryNearFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        /** Two layer particles.*/
        if (pos_div_[inner_neighborhood.j_[n]] < threshold_by_dimensions_ &&
            inner_neighborhood.r_ij_[n] < smoothing_length_)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
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
} // namespace fluid_dynamics
} // namespace SPH
