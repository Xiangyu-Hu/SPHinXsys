#include "base_body_part.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
BodyPart::BodyPart(SPHBody &sph_body, const std::string &body_part_name)
    : sph_body_(sph_body), part_id_(sph_body.getNewBodyPartID()),
      body_part_name_(body_part_name),
      base_particles_(sph_body.getBaseParticles()),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position")),
      dv_index_list_(nullptr), sv_range_size_(nullptr) {}
//=================================================================================================//
BodyPartByParticle::BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name)
    : BodyPart(sph_body, body_part_name),
      body_part_bounds_(Vecd::Zero(), Vecd::Zero()), body_part_bounds_set_(false) {}
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
}
//=============================================================================================//
BodyPartByCell::BodyPartByCell(RealBody &real_body, const std::string &body_part_name)
    : BodyPart(real_body, body_part_name),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
//=============================================================================================//
size_t BodyPartByCell::SizeOfLoopRange()
{
    size_t size_of_loop_range = 0;
    for (size_t i = 0; i != body_part_cells_.size(); ++i)
    {
        size_of_loop_range += body_part_cells_[i]->size();
    }
    return size_of_loop_range;
}
//=================================================================================================//
void BodyPartByCell::tagCells(TaggingCellMethod &tagging_cell_method)
{
    ConcurrentVec<size_t> cell_list;
    Arrayi all_cells = cell_linked_list_.AllCells();
    Real grid_spacing = cell_linked_list_.GridSpacing();
    mesh_parallel_for(
        MeshRange(Arrayi::Zero(), all_cells),
        [&](const Arrayi &cell_index)
        {
            bool is_included = false;
            mesh_for_each(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                all_cells.min(cell_index + 2 * Arrayi::Ones()),
                [&](const Arrayi &neighbor_cell_index)
                {
                    if (tagging_cell_method(
                            cell_linked_list_
                                .CellPositionFromIndex(neighbor_cell_index),
                            grid_spacing))
                    {
                        is_included = true;
                    }
                });
            if (is_included == true)
                cell_list.push_back(cell_linked_list_.LinearCellIndexFromCellIndex(cell_index));
        });

    ConcurrentIndexVector *cell_index_lists = cell_linked_list_.getCellIndexLists();
    for (size_t i = 0; i != cell_list.size(); ++i)
    {
        body_part_cells_.push_back(&cell_index_lists[cell_list[i]]);
    }

    dv_index_list_ = base_particles_.addUniqueDiscreteVariableOnly<UnsignedInt>(
        body_part_name_, cell_list.size(), [&](size_t i) -> Real
        { return cell_list[i]; });
    sv_range_size_ = base_particles_.addUniqueSingularVariableOnly<UnsignedInt>(
        body_part_name_ + "_Size", cell_list.size());
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
} // namespace SPH
