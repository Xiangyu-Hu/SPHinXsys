#include "base_body_part.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
BodyPart::BodyPart(SPHBody &sph_body, const std::string &body_part_name)
    : sph_body_(sph_body), part_id_(sph_body.getNewBodyPartID()),
      body_part_name_(body_part_name),
      base_particles_(sph_body.getBaseParticles()),
      dv_index_list_(nullptr), sv_range_size_(nullptr),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position")) {}
//=================================================================================================//
BodyPartByParticle::BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name)
    : BodyPart(sph_body, body_part_name),
      body_part_bounds_(Vecd::Zero(), Vecd::Zero()), body_part_bounds_set_(false)
{
    sph_body.addBodyPartByParticle(this);
}
//=================================================================================================//
void BodyPartByParticle::tagParticles(TaggingParticleMethod &tagging_particle_method)
{
    for (size_t i = 0; i < base_particles_.TotalRealParticles(); ++i)
    {
        tagging_particle_method(i);
    }
    dv_index_list_ = base_particles_.addUniqueDiscreteVariableOnly<UnsignedInt>(
        body_part_name_, body_part_particles_.size(), [&](size_t i) -> Real
        { return body_part_particles_[i]; });
    sv_range_size_ = base_particles_.addUniqueSingularVariableOnly<UnsignedInt>(
        body_part_name_ + "_Size", body_part_particles_.size());
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
    ConcurrentIndexVector cell_indexes;
    cell_linked_list_.tagBodyPartByCell(body_part_cells_, cell_indexes, tagging_cell_method);
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
//=================================================================================================//
void BodyRegionByParticle::tagByContain(size_t particle_index)
{
    if (body_part_shape_.checkContain(pos_[particle_index]))
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
BodySurface::BodySurface(SPHBody &sph_body)
    : BodyPartByParticle(sph_body, "BodySurface"),
      particle_spacing_min_(sph_body.sph_adaptation_->MinimumSpacing())
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurface::tagNearSurface, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurface::tagNearSurface(size_t particle_index)
{
    Real phi = sph_body_.getInitialShape().findSignedDistance(pos_[particle_index]);
    if (fabs(phi) < particle_spacing_min_)
        body_part_particles_.push_back(particle_index);
}
//=================================================================================================//
BodySurfaceLayer::BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness)
    : BodyPartByParticle(sph_body, "InnerLayers"),
      thickness_threshold_(sph_body.sph_adaptation_->ReferenceSpacing() * layer_thickness)
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurfaceLayer::tagSurfaceLayer, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurfaceLayer::tagSurfaceLayer(size_t particle_index)
{
    Real distance = fabs(sph_body_.getInitialShape().findSignedDistance(pos_[particle_index]));
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
NearShapeSurface::NearShapeSurface(RealBody &real_body, LevelSetShape &level_set_shape)
    : BodyPartByCell(real_body, level_set_shape.getName()), level_set_shape_(level_set_shape)
{
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr)
    : BodyPartByCell(real_body, shape_ptr->getName()),
      level_set_shape_(level_set_shape_keeper_.createRef<LevelSetShape>(real_body, *shape_ptr.get(), true))
{
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body)
    : BodyPartByCell(real_body, "NearShapeSurface"),
      level_set_shape_(DynamicCast<LevelSetShape>(this, real_body.getInitialShape()))
{
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, const std::string &sub_shape_name)
    : BodyPartByCell(real_body, sub_shape_name),
      level_set_shape_(
          DynamicCast<LevelSetShape>(this, *DynamicCast<ComplexShape>(this, real_body.getInitialShape())
                                                .getSubShapeByName(sub_shape_name)))
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
AlignedBoxPart::AlignedBoxPart(const std::string &name, const AlignedBox &aligned_box)
    : aligned_box_(*sv_aligned_box_keeper_
                        .createPtr<SingularVariable<AlignedBox>>("AlignedBox" + name, aligned_box)
                        ->Data()) {}
//=================================================================================================//
AlignedBoxPartByParticle::AlignedBoxPartByParticle(RealBody &real_body, const AlignedBox &aligned_box)
    : BodyPartByParticle(real_body, "AlignedBoxByParticle"),
      AlignedBoxPart(body_part_name_, aligned_box)
{
    TaggingParticleMethod tagging_particle_method =
        std::bind(&AlignedBoxPartByParticle::tagByContain, this, _1);
    tagParticles(tagging_particle_method);
}
//=================================================================================================//
void AlignedBoxPartByParticle::tagByContain(size_t particle_index)
{
    if (aligned_box_.checkContain(pos_[particle_index]))
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
AlignedBoxPartByCell::AlignedBoxPartByCell(RealBody &real_body, const AlignedBox &aligned_box)
    : BodyPartByCell(real_body, "AlignedBoxByCell"),
      AlignedBoxPart(body_part_name_, aligned_box)
{
    TaggingCellMethod tagging_cell_method =
        std::bind(&AlignedBoxPartByCell::checkNotFar, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
bool AlignedBoxPartByCell::checkNotFar(Vecd cell_position, Real threshold)
{
    return aligned_box_.checkNotFar(cell_position, threshold);
}
//=================================================================================================//
} // namespace SPH
