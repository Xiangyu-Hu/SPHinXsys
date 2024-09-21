#include "base_body_part.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
void BodyPartByParticle::tagParticles(TaggingParticleMethod &tagging_particle_method)
{
    for (size_t i = 0; i < base_particles_.TotalRealParticles(); ++i)
    {
        tagging_particle_method(i);
    }
};
//=============================================================================================//
size_t BodyPartByCell::SizeOfLoopRange()
{
    size_t size_of_loop_range = 0;
    for (size_t i = 0; i != body_part_cells_.size(); ++i)
    {
        size_of_loop_range += body_part_cells_[i]->size();
    }
    return size_of_loop_range;
};
//=================================================================================================//
void BodyPartByCell::tagCells(TaggingCellMethod &tagging_cell_method)
{
    cell_linked_list_.tagBodyPartByCell(body_part_cells_, tagging_cell_method);
}
//=================================================================================================//
BodyRegionByParticle::
    BodyRegionByParticle(SPHBody &sph_body, Shape &body_part_shape)
    : BodyPartByParticle(sph_body, body_part_shape.getName()),
      body_part_shape_(body_part_shape)
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodyRegionByParticle::tagByContain, this, _1);
    tagParticles(tagging_particle_method);
}
//=================================================================================================//
BodyRegionByParticle::BodyRegionByParticle(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
    : BodyRegionByParticle(sph_body, *shape_ptr.get())
{
    shape_ptr_keeper_.assignRef(shape_ptr);
}
//==
//=================================================================================================//
void BodyRegionByParticle::tagByContain(size_t particle_index)
{
    if (body_part_shape_.checkContain(pos_[particle_index]))
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
BodySurface::BodySurface(SPHBody &sph_body, Shape &body_shape)
    : BodyPartByParticle(sph_body, "BodySurface"),
      body_shape_(body_shape),
      particle_spacing_min_(sph_body.sph_adaptation_->MinimumSpacing())
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurface::tagNearSurface, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurface::tagNearSurface(size_t particle_index)
{
    Real phi = body_shape_.findSignedDistance(pos_[particle_index]);
    if (fabs(phi) < particle_spacing_min_)
        body_part_particles_.push_back(particle_index);
}
//=================================================================================================//
BodySurfaceLayer::BodySurfaceLayer(SPHBody &sph_body, Shape &body_shape, Real layer_thickness)
    : BodyPartByParticle(sph_body, "InnerLayers"),
      body_shape_(body_shape),
      thickness_threshold_(sph_body.sph_adaptation_->ReferenceSpacing() * layer_thickness)
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurfaceLayer::tagSurfaceLayer, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurfaceLayer::tagSurfaceLayer(size_t particle_index)
{
    Real distance = fabs(body_shape_.findSignedDistance(pos_[particle_index]));
    if (distance < thickness_threshold_)
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
BodyRegionByCell::BodyRegionByCell(RealBody &real_body, Shape &body_part_shape)
    : BodyPartByCell(real_body, body_part_shape.getName()),
      body_part_shape_(body_part_shape)
{
    TaggingCellMethod tagging_cell_method = std::bind(&BodyRegionByCell::checkNotFar, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
BodyRegionByCell::BodyRegionByCell(RealBody &real_body, SharedPtr<Shape> shape_ptr)
    : BodyRegionByCell(real_body, *shape_ptr.get())
{
    shape_ptr_keeper_.assignRef(shape_ptr);
}
//=================================================================================================//
bool BodyRegionByCell::checkNotFar(Vecd cell_position, Real threshold)
{
    return body_part_shape_.checkNotFar(cell_position, threshold);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, Shape &shape)
    : BodyPartByCell(real_body, shape.getName()),
      level_set_shape_(DynamicCast<LevelSetShape>(this, shape))
{
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
bool NearShapeSurface::checkNearSurface(Vecd cell_position, Real threshold)
{
    return level_set_shape_.checkNearSurface(cell_position, threshold);
}
//=================================================================================================//
} // namespace SPH
