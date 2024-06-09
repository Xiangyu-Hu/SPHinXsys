#include "general_constraint.h"

namespace SPH
{
//=================================================================================================//
ShapeSurfaceBounding::ShapeSurfaceBounding(NearShapeSurface &near_shape_surface)
    : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
      GeneralDataDelegateSimple(sph_body_), pos_(particles_->pos_),
      constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
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
ComplexShapeBounding::ComplexShapeBounding(SPHBody& sph_body, ComplexShape& complex_shape)
    : BaseLocalDynamics<SPHBody>(sph_body),
    GeneralDataDelegateSimple(sph_body_), pos_(particles_->pos_),
    constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
{
    level_set_shapes_ = complex_shape.getLevelSetShapes();
};
//=================================================================================================//
void ComplexShapeBounding::update(size_t index_i, Real dt)
{
    for (size_t i = 0; i != level_set_shapes_.size(); ++i)
    {
         Real phi = level_set_shapes_[i]->findSignedDistance(pos_[index_i]);

        if (phi > -constrained_distance_)
        {
            Vecd unit_normal = level_set_shapes_[i]->findNormalDirection(pos_[index_i]);
            pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
        }
    }
}
//=================================================================================================//
} // namespace SPH
