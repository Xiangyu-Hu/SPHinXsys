#include "general_constraint.h"

namespace SPH
{
//=================================================================================================//
ShapeSurfaceBounding::ShapeSurfaceBounding(NearShapeSurface &near_shape_surface)
    : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      constrained_distance_(0.5 * sph_body_.getSPHAdaptation().MinimumSpacing())
{
    level_set_shape_ = &near_shape_surface.getLevelSetShape();
}
//=================================================================================================//
void ShapeSurfaceBounding::update(size_t index_i, Real dt)
{
    Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

    if (phi > -constrained_distance_)
    {
        Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
        pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
    }
}
//=================================================================================================//
} // namespace SPH
