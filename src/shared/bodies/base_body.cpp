#include "base_body.h"

#include "base_body_relation.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr, const std::string &body_name)
    : sph_system_(sph_system), body_name_(body_name), newly_updated_(true), base_particles_(nullptr),
      body_shape_(shape_ptr_keeper_.assignPtr(shape_ptr)),
      sph_adaptation_(sph_adaptation_ptr_keeper_.createPtr<SPHAdaptation>(*this)),
      base_material_(nullptr)
{
    sph_system_.sph_bodies_.push_back(this);
}
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
    : SPHBody(sph_system, shape_ptr, shape_ptr->getName()) {}
//=================================================================================================//
BoundingBox SPHBody::getSPHSystemBounds()
{
    return sph_system_.system_domain_bounds_;
}
//=================================================================================================//
SPHSystem &SPHBody::getSPHSystem()
{
    return sph_system_;
}
//=================================================================================================//
void SPHBody::allocateConfigurationMemoriesForBufferParticles()
{
    for (size_t i = 0; i < body_relations_.size(); i++)
    {
        body_relations_[i]->resizeConfiguration();
    }
}
//=================================================================================================//
BoundingBox SPHBody::getBodyShapeBounds()
{
    return body_shape_->getBounds();
}
//=================================================================================================//
void SPHBody::defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
{
    sph_adaptation_->resetAdaptationRatios(h_spacing_ratio, new_system_refinement_ratio);
}
//=================================================================================================//
void SPHBody::writeParticlesToVtuFile(std::ostream &output_file)
{
    base_particles_->writeParticlesToVtk(output_file);
}
//=================================================================================================//
void SPHBody::writeParticlesToVtpFile(std::ofstream &output_file)
{
    base_particles_->writeParticlesToVtk(output_file);
}
//=================================================================================================//
void SPHBody::writeSurfaceParticlesToVtuFile(std::ofstream &output_file, BodySurface &surface_particles)
{
    base_particles_->writeSurfaceParticlesToVtuFile(output_file, surface_particles);
}
//=================================================================================================//
void SPHBody::writeParticlesToPltFile(std::ofstream &output_file)
{
    base_particles_->writeParticlesToPltFile(output_file);
}
//=================================================================================================//
void SPHBody::writeParticlesToXmlForRestart(std::string &filefullpath)
{
    base_particles_->writeParticlesToXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::readParticlesFromXmlForRestart(std::string &filefullpath)
{
    base_particles_->readParticleFromXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::writeToXmlForReloadParticle(std::string &filefullpath)
{
    base_particles_->writeToXmlForReloadParticle(filefullpath);
}
//=================================================================================================//
void SPHBody::readFromXmlForReloadParticle(std::string &filefullpath)
{
    base_particles_->readFromXmlForReloadParticle(filefullpath);
}
//=================================================================================================//
BaseCellLinkedList &RealBody::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ = sph_adaptation_->createCellLinkedList(getSPHSystemBounds(), *this);
        cell_linked_list_created_ = true;
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
void RealBody::updateCellLinkedList()
{
    getCellLinkedList().UpdateCellLists(*base_particles_);
    base_particles_->total_ghost_particles_ = 0;
}
//=================================================================================================//
void RealBody::updateCellLinkedListWithParticleSort(size_t particle_sorting_period)
{
    if (iteration_count_ % particle_sorting_period == 0)
    {
        base_particles_->sortParticles(getCellLinkedList());
    }

    iteration_count_++;
    updateCellLinkedList();
}
//=================================================================================================//
} // namespace SPH
