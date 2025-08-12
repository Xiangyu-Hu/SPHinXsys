#include "base_body_part.h"

#include "base_particles.hpp"
#include "cell_linked_list.hpp"
namespace SPH
{
//=================================================================================================//
BodyPart::BodyPart(SPHBody &sph_body)
    : sph_body_(sph_body), base_particles_(sph_body.getBaseParticles()),
      part_id_(base_particles_.getNewBodyPartID()),
      part_name_(sph_body.getName() + "Part" + std::to_string(part_id_)),
      sph_adaptation_(sph_body.getSPHAdaptation()),
      sv_range_size_(nullptr),
      dv_body_part_id_(base_particles_.registerStateVariable<int>(part_name_ + "ID")),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position")) {}
//=================================================================================================//
BaseCellLinkedList &BodyPart::getCellLinkedList()
{
    RealBody &real_body = DynamicCast<RealBody>(this, sph_body_);
    return real_body.getCellLinkedList();
}
//=================================================================================================//
BodyPartByID::BodyPartByID(SPHBody &sph_body) : BodyPart(sph_body)
{
    sv_range_size_ = base_particles_.svTotalRealParticles();
}
//=================================================================================================//
BodyPartByParticle::BodyPartByParticle(SPHBody &sph_body)
    : BodyPart(sph_body), body_part_bounds_(Vecd::Zero(), Vecd::Zero()),
      body_part_bounds_set_(false)
{
    base_particles_.addBodyPartByParticle(this);
    base_particles_.addEvolvingVariable<int>(dv_body_part_id_);
}
//=================================================================================================//
void BodyPartByParticle::setBodyPartBounds(BoundingBox bbox)
{
    body_part_bounds_ = bbox;
    body_part_bounds_set_ = true;
}
//=================================================================================================//
BoundingBox BodyPartByParticle::getBodyPartBounds()
{
    if (!body_part_bounds_set_)
        std::cout << "WARNING: the body part bounds are not set for BodyPartByParticle." << std::endl;
    return body_part_bounds_;
}
//=================================================================================================//
void BodyPartByParticle::tagParticles(TaggingParticleMethod &tagging_particle_method)
{
    for (size_t i = 0; i != base_particles_.TotalRealParticles(); ++i)
    {
        if (tagging_particle_method(i))
        {
            dv_body_part_id_->setValue(i, part_id_);
            body_part_particles_.push_back(i);
        }
    }

    dv_particle_list_ = unique_variable_ptrs_.createPtr<DiscreteVariable<UnsignedInt>>(
        part_name_, body_part_particles_.size(), [&](size_t i)
        { return body_part_particles_[i]; });
    sv_range_size_ = unique_variable_ptrs_.createPtr<SingularVariable<UnsignedInt>>(
        part_name_ + "_Size", body_part_particles_.size());
}
//=================================================================================================//
BodyPartByCell::BodyPartByCell(RealBody &real_body)
    : BodyPart(real_body), cell_linked_list_(real_body.getCellLinkedList()),
      dv_cell_list_(nullptr),
      dv_particle_index_(cell_linked_list_.dvParticleIndex()),
      dv_cell_offset_(cell_linked_list_.dvCellOffset()) {}
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
    ConcurrentIndexVector cell_indexes;
    cell_linked_list_.tagBodyPartByCell(body_part_cells_, cell_indexes, tagging_cell_method);

    for (size_t i = 0; i != body_part_cells_.size(); ++i)
    {
        ConcurrentIndexVector &particle_indexes = *body_part_cells_[i];
        for (size_t num = 0; num < particle_indexes.size(); ++num)
        {
            dv_body_part_id_->setValue(particle_indexes[num], part_id_);
        }
    }
    dv_cell_list_ = unique_variable_ptrs_.createPtr<DiscreteVariable<UnsignedInt>>(
        part_name_, cell_indexes.size(), [&](size_t i)
        { return cell_indexes[i]; });
    sv_range_size_ = unique_variable_ptrs_.createPtr<SingularVariable<UnsignedInt>>(
        part_name_ + "_Size", cell_indexes.size());
}
//=================================================================================================//
BodyRegionByParticle::
    BodyRegionByParticle(SPHBody &sph_body, Shape &body_part_shape)
    : BodyPartByParticle(sph_body), body_part_shape_(body_part_shape)
{
    alias_ = body_part_shape_.getName();
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
bool BodyRegionByParticle::tagByContain(size_t particle_index)
{
    return body_part_shape_.checkContain(pos_[particle_index]);
}
//=================================================================================================//
BodyRegionByCell::BodyRegionByCell(RealBody &real_body, Shape &body_part_shape)
    : BodyPartByCell(real_body), body_part_shape_(body_part_shape)
{
    alias_ = body_part_shape_.getName();
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
AlignedBoxPart::AlignedBoxPart(const std::string &part_name, const AlignedBox &aligned_box)
    : aligned_box_(*sv_aligned_box_keeper_
                        .createPtr<SingularVariable<AlignedBox>>(part_name, aligned_box)
                        ->Data()) {}
//=================================================================================================//
AlignedBoxByParticle::AlignedBoxByParticle(RealBody &real_body, const AlignedBox &aligned_box)
    : BodyPartByParticle(real_body), AlignedBoxPart(part_name_, aligned_box)
{
    TaggingParticleMethod tagging_particle_method =
        std::bind(&AlignedBoxByParticle::tagByContain, this, _1);
    tagParticles(tagging_particle_method);
}
//=================================================================================================//
bool AlignedBoxByParticle::tagByContain(size_t particle_index)
{
    return aligned_box_.checkContain(pos_[particle_index]);
}
//=================================================================================================//
AlignedBoxByCell::AlignedBoxByCell(RealBody &real_body, const AlignedBox &aligned_box)
    : BodyPartByCell(real_body), AlignedBoxPart(part_name_, aligned_box)
{
    TaggingCellMethod tagging_cell_method =
        std::bind(&AlignedBoxByCell::checkNotFar, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
bool AlignedBoxByCell::checkNotFar(Vecd cell_position, Real threshold)
{
    return aligned_box_.checkNotFar(cell_position, threshold);
}
//=================================================================================================//
} // namespace SPH
