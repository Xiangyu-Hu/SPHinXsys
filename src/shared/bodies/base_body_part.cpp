#include "base_body_part.h"

#include "base_body.h"
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
SPHSystem &BodyPart::getSPHSystem() { return sph_body_.getSPHSystem(); }
//=================================================================================================//
BodyPartByID::BodyPartByID(SPHBody &sph_body) : BodyPart(sph_body)
{
    sv_range_size_ = base_particles_.svTotalRealParticles();
}
//=================================================================================================//
BodyPartByParticle::BodyPartByParticle(SPHBody &sph_body)
    : BodyPart(sph_body)
{
    base_particles_.addBodyPartByParticle(this);
    base_particles_.addEvolvingVariable<int>(dv_body_part_id_);
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
//=============================================================================================//
template <typename DataType>
BodyPartByParticleVariableRange<DataType>::BodyPartByParticleVariableRange(
	SPHBody &sph_body, const std::string &name, const std::string &variable_name,
	DataType lower_bound, DataType upper_bound)
	: BodyPartByParticle(sph_body),
	  variable_(base_particles_.getVariableDataByName<DataType>(variable_name)),
	  variable_name_(variable_name),
	  lower_bound_(lower_bound), upper_bound_(upper_bound)
{
	if (lower_bound_ > upper_bound_)
	{
		throw std::invalid_argument(
			"BodyPartByParticleVariableRange requires upper_bound >= lower_bound for variable '" +
			variable_name_ + "'.");
	}

	alias_ = name;
	TaggingParticleMethod tagging_particle_method =
		std::bind(&BodyPartByParticleVariableRange::tagByRange, this, _1);
	tagParticles(tagging_particle_method);
}
//=============================================================================================//
template <typename DataType>
bool BodyPartByParticleVariableRange<DataType>::tagByRange(size_t index_i) const
{
	const DataType value = variable_[index_i];
	return lower_bound_ <= value && value <= upper_bound_;
}
//=============================================================================================//
template class BodyPartByParticleVariableRange<Real>;
template class BodyPartByParticleVariableRange<int>;
//=============================================================================================//
BodyPartByParticleBoolean::BodyPartByParticleBoolean(SPHBody &sph_body)
	: BodyPartByParticle(sph_body)
{
}
//=============================================================================================//
void BodyPartByParticleBoolean::initializeBooleanBodyPart(
	const std::string &name,
	TaggingParticleMethod &tagging_particle_method)
{
	alias_ = name;
	tagParticles(tagging_particle_method);
}
//=============================================================================================//
void BodyPartByParticleBoolean::validateBodyPartsInput(
	const StdVec<BodyPartByParticle *> &body_parts,
	const std::string &operation_name) const
{
	if (body_parts.empty())
	{
		throw std::invalid_argument(
			operation_name + " requires at least one input body part.");
	}

	for (size_t i = 0; i < body_parts.size(); ++i)
	{
		BodyPartByParticle *part = body_parts[i];
		if (part == nullptr)
		{
			throw std::invalid_argument(
				operation_name + " received a null body part at index " +
				std::to_string(i) + ".");
		}

		if (&part->getBaseParticles() != &base_particles_)
		{
			throw std::invalid_argument(
				operation_name + " requires all inputs to belong to the same SPHBody.");
		}
	}
}
//=============================================================================================//
void BodyPartByParticleBoolean::countParticleInclusions(
	const StdVec<BodyPartByParticle *> &body_parts,
	StdVec<int> &inclusion_mask) const
{
	const size_t total_real_particles = base_particles_.TotalRealParticles();
	if (inclusion_mask.size() != total_real_particles)
	{
		throw std::invalid_argument(
			"BodyPartByParticleBoolean::countParticleInclusions expects a pre-sized inclusion mask.");
	}
	std::fill(inclusion_mask.begin(), inclusion_mask.end(), 0);

	for (BodyPartByParticle *part : body_parts)
	{
		for (size_t index_i : part->LoopRange())
		{
			++inclusion_mask[index_i];
		}
	}
}
//=============================================================================================//
BodyPartByParticleUnion::BodyPartByParticleUnion(
	SPHBody &sph_body, const std::string &name,
	const StdVec<BodyPartByParticle *> &body_parts)
	: BodyPartByParticleBoolean(sph_body),
	  inclusion_mask_(base_particles_.TotalRealParticles())
{
	validateBodyPartsInput(body_parts, "BodyPartByParticleUnion");
	countParticleInclusions(body_parts, inclusion_mask_);
	TaggingParticleMethod tag_method = std::bind(&BodyPartByParticleUnion::tagByUnion, this, _1);
	initializeBooleanBodyPart(name, tag_method);
}
//=============================================================================================//
BodyPartByParticleUnion::BodyPartByParticleUnion(
	SPHBody &sph_body, const std::string &name,
	std::initializer_list<BodyPartByParticle *> body_parts)
	: BodyPartByParticleUnion(sph_body, name, StdVec<BodyPartByParticle *>(body_parts))
{
}
//=============================================================================================//
bool BodyPartByParticleUnion::tagByUnion(size_t particle_index) const
{
	return inclusion_mask_[particle_index] != 0;
}
//=============================================================================================//
BodyPartByParticleIntersection::BodyPartByParticleIntersection(
	SPHBody &sph_body, const std::string &name,
	const StdVec<BodyPartByParticle *> &body_parts)
	: BodyPartByParticleBoolean(sph_body),
	  inclusion_mask_(base_particles_.TotalRealParticles()),
	  required_hits_(body_parts.size())
{
	validateBodyPartsInput(body_parts, "BodyPartByParticleIntersection");
	countParticleInclusions(body_parts, inclusion_mask_);
	TaggingParticleMethod tag_method = std::bind(&BodyPartByParticleIntersection::tagByIntersection, this, _1);
	initializeBooleanBodyPart(name, tag_method);
}
//=============================================================================================//
BodyPartByParticleIntersection::BodyPartByParticleIntersection(
	SPHBody &sph_body, const std::string &name,
	std::initializer_list<BodyPartByParticle *> body_parts)
	: BodyPartByParticleIntersection(
		sph_body, name, StdVec<BodyPartByParticle *>(body_parts))
{
}
//=============================================================================================//
bool BodyPartByParticleIntersection::tagByIntersection(size_t particle_index) const
{
	return inclusion_mask_[particle_index] == static_cast<int>(required_hits_);
}
//=============================================================================================//
BodyPartByParticleDifference::BodyPartByParticleDifference(
	SPHBody &sph_body, const std::string &name,
	const StdVec<BodyPartByParticle *> &body_parts)
	: BodyPartByParticleBoolean(sph_body),
	  minuend_mask_(base_particles_.TotalRealParticles(), 0),
	  exclusion_mask_(base_particles_.TotalRealParticles(), 0)
{
	validateBodyPartsInput(body_parts, "BodyPartByParticleDifference");
	buildDifferenceMasks(body_parts);
	TaggingParticleMethod tag_method = std::bind(&BodyPartByParticleDifference::tagByDifference, this, _1);
	initializeBooleanBodyPart(name, tag_method);
}
//=============================================================================================//
void BodyPartByParticleDifference::buildDifferenceMasks(
	const StdVec<BodyPartByParticle *> &body_parts)
{
	for (size_t index_i : body_parts.front()->LoopRange())
	{
		minuend_mask_[index_i] = 1;
	}

	if (body_parts.size() == 1) return;

	StdVec<BodyPartByParticle *> exclusion_parts(body_parts.begin() + 1, body_parts.end());
	countParticleInclusions(exclusion_parts, exclusion_mask_);
}
//=============================================================================================//
BodyPartByParticleDifference::BodyPartByParticleDifference(
	SPHBody &sph_body, const std::string &name,
	std::initializer_list<BodyPartByParticle *> body_parts)
	: BodyPartByParticleDifference(
		sph_body, name, StdVec<BodyPartByParticle *>(body_parts))
{
}
//=============================================================================================//
bool BodyPartByParticleDifference::tagByDifference(size_t particle_index) const
{
	return minuend_mask_[particle_index] != 0 && exclusion_mask_[particle_index] == 0;
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
    shape_keeper_.assignRef(shape_ptr);
}
//=================================================================================================//
bool BodyRegionByParticle::tagByContain(size_t particle_index)
{
    return body_part_shape_.checkContain(pos_[particle_index]);
}
//=================================================================================================//
BodySurface::BodySurface(SPHBody &sph_body)
    : BodyPartByParticle(sph_body),
      particle_spacing_min_(sph_body.getSPHAdaptation().MinimumSpacing())
{
    alias_ = "BodySurface";
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurface::tagNearSurface, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
bool BodySurface::tagNearSurface(size_t particle_index)
{
    Real phi = sph_body_.getInitialShape().findSignedDistance(pos_[particle_index]);
    return fabs(phi) < particle_spacing_min_;
}
//=================================================================================================//
BodySurfaceLayer::BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness)
    : BodyPartByParticle(sph_body),
      thickness_threshold_(sph_body.getSPHAdaptation().ReferenceSpacing() * layer_thickness)
{
    alias_ = "InnerLayers";
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurfaceLayer::tagSurfaceLayer, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
bool BodySurfaceLayer::tagSurfaceLayer(size_t particle_index)
{
    Real distance = fabs(sph_body_.getInitialShape().findSignedDistance(pos_[particle_index]));
    return distance < thickness_threshold_;
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
    shape_keeper_.assignRef(shape_ptr);
}
//=================================================================================================//
bool BodyRegionByCell::checkNotFar(Vecd cell_position, Real threshold)
{
    return body_part_shape_.checkNotFar(cell_position, threshold);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, LevelSetShape &level_set_shape)
    : BodyPartByCell(real_body), level_set_shape_(level_set_shape)
{
    alias_ = level_set_shape.getName();
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr)
    : BodyPartByCell(real_body),
      level_set_shape_(level_set_shape_keeper_.createRef<LevelSetShape>(real_body, *shape_ptr.get(), true))
{
    alias_ = level_set_shape_.getName();
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body)
    : BodyPartByCell(real_body),
      level_set_shape_(DynamicCast<LevelSetShape>(this, real_body.getInitialShape()))
{
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
NearShapeSurface::NearShapeSurface(RealBody &real_body, const std::string &sub_shape_name)
    : BodyPartByCell(real_body),
      level_set_shape_(
          DynamicCast<LevelSetShape>(this, *DynamicCast<ComplexShape>(this, real_body.getInitialShape())
                                                .getSubShapeByName(sub_shape_name)))
{
    alias_ = sub_shape_name;
    TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
    tagCells(tagging_cell_method);
}
//=================================================================================================//
bool NearShapeSurface::checkNearSurface(Vecd cell_position, Real threshold)
{
    return level_set_shape_.checkNearSurface(cell_position, threshold);
}
//=================================================================================================//
AlignedBoxPart::AlignedBoxPart(SPHSystem &sph_system, const std::string &part_name, const AlignedBox &aligned_box)
    : sph_system_(sph_system),
      aligned_box_(*sv_aligned_box_keeper_
                        .createPtr<SingularVariable<AlignedBox>>(part_name, aligned_box)
                        ->Data())
{
    std::cout << part_name << " direction facing to fluid domain: "
              << aligned_box_.getTransform().xformFrameVecToBase(Vecd::UnitX()) << std::endl;
}
//=================================================================================================//
AlignedBoxByParticle::AlignedBoxByParticle(RealBody &real_body, const AlignedBox &aligned_box)
    : BodyPartByParticle(real_body), AlignedBoxPart(real_body.getSPHSystem(), part_name_, aligned_box)
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
    : BodyPartByCell(real_body), AlignedBoxPart(real_body.getSPHSystem(), part_name_, aligned_box)
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
