#include "base_body.h"

#include "base_body_relation.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, Shape &shape, const std::string &name)
    : sph_system_(sph_system), body_name_(name), newly_updated_(true),
      base_particles_(nullptr), is_bound_set_(false), initial_shape_(&shape),
      sph_adaptation_(sph_adaptation_ptr_keeper_.createPtr<SPHAdaptation>(sph_system.ReferenceResolution())),
      base_material_(base_material_ptr_keeper_.createPtr<BaseMaterial>())
{
    sph_system_.addSPHBody(this);
}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, Shape &shape)
    : SPHBody(sph_system, shape, shape.getName()) {}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, const std::string &name)
    : SPHBody(sph_system, makeShared<DefaultShape>(name)) {}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr, const std::string &name)
    : SPHBody(sph_system, *shape_ptr.get(), name)
{
    shape_ptr_keeper_.assignPtr(shape_ptr);
}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
    : SPHBody(sph_system, shape_ptr, shape_ptr->getName()) {}
//=================================================================================================//
BoundingBoxd SPHBody::getSPHSystemBounds()
{
    return sph_system_.getSystemDomainBounds();
}
//=================================================================================================//
SPHSystem &SPHBody::getSPHSystem()
{
    return sph_system_;
}
//=================================================================================================//
BaseParticles &SPHBody::getBaseParticles()
{
    if (base_particles_ == nullptr)
    {
        std::cout << "\n Error: BaseParticle not generated yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_particles_;
};
//=================================================================================================//
BaseMaterial &SPHBody::getBaseMaterial()
{
    if (base_material_ == nullptr)
    {
        std::cout << "\n Error: BaseMaterial not generated yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_material_;
};
//=================================================================================================//
void SPHBody::setSPHBodyBounds(const BoundingBoxd &bound)
{
    bound_ = bound;
    is_bound_set_ = true;
}
//=================================================================================================//
BoundingBoxd SPHBody::getSPHBodyBounds()
{
    return is_bound_set_ ? bound_ : initial_shape_->getBounds();
}
//=================================================================================================//
void SPHBody::defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
{
    sph_adaptation_->resetAdaptationRatios(h_spacing_ratio, new_system_refinement_ratio);
}
//=================================================================================================//
BaseCellLinkedList &RealBody::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_->createCellLinkedList(getSPHSystemBounds(), *base_particles_);
        cell_linked_list_created_ = true;
        cell_linked_list_ptr_.get()->setName(getName() + "CellLinkedList");
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
void RealBody::updateCellLinkedList()
{
    getCellLinkedList().UpdateCellLists(*base_particles_);
}
//=================================================================================================//
void RealBody::addRealBodyToSPHSystem()
{
    sph_system_.addRealBody(this);
}
//=================================================================================================//
} // namespace SPH
