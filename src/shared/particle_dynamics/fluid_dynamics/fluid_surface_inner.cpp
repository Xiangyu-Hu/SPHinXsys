#include "fluid_surface_inner.hpp"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
FreeSurfaceIndicationInner::
    FreeSurfaceIndicationInner(BaseInnerRelation &inner_relation, Real threshold)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      threshold_by_dimensions_(threshold * (Real)Dimensions),
      surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
      smoothing_length_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSmoothingLength())
{
    particles_->registerVariable(pos_div_, "PositionDivergence");
}
//=================================================================================================//
void FreeSurfaceIndicationInner::update(size_t index_i, Real dt)
{
    surface_indicator_[index_i] = 1;
    if (pos_div_[index_i] > threshold_by_dimensions_ && !isVeryNearFreeSurface(index_i))
        surface_indicator_[index_i] = 0;
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
ColorFunctionGradientInner::ColorFunctionGradientInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
      pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
      threshold_by_dimensions_((0.75 * (Real)Dimensions))
{
    particles_->registerVariable(color_grad_, "ColorGradient");
    particles_->registerVariable(surface_norm_, "SurfaceNormal");
}
//=================================================================================================//
ColorFunctionGradientInterpolationInner::ColorFunctionGradientInterpolationInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation), Vol_(particles_->Vol_),
      surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
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
      acc_prior_(particles_->acc_prior_), surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
      color_grad_(*particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")) {}
//=================================================================================================//
SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation)
    : SurfaceTensionAccelerationInner(inner_relation, 1.0) {}
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//