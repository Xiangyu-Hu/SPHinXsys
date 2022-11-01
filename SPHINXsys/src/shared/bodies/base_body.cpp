#include "base_body.h"

#include "sph_system.h"
#include "base_particles.hpp"
#include "base_body_relation.h"

namespace SPH
{
	//=================================================================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
		: sph_system_(sph_system), newly_updated_(true),
		  body_shape_(shape_ptr_keeper_.assignPtr(shape_ptr)),
		  sph_adaptation_(sph_adaptation_ptr_keeper_.createPtr<SPHAdaptation>(*this)),
		  base_material_(nullptr), base_particles_(nullptr)
	{
		sph_system_.sph_bodies_.push_back(this);
	}
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
	void SPHBody::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	void SPHBody::allocateConfigurationMemoriesForBufferParticles()
	{
		for (size_t i = 0; i < body_relations_.size(); i++)
		{
			body_relations_[i]->updateConfigurationMemories();
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
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToVtpFile(std::ofstream &output_file)
	{
		base_particles_->writeParticlesToVtk(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeSurfaceParticlesToVtuFile(std::ofstream &output_file, BodySurface &surface_particles)
	{
		base_particles_->writeSurfaceParticlesToVtuFile(output_file, surface_particles);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToPltFile(std::ofstream &output_file)
	{
		if (newly_updated_)
			base_particles_->writeParticlesToPltFile(output_file);
		newly_updated_ = false;
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
	RealBody::RealBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
		: SPHBody(sph_system, shape_ptr)
		, system_domain_bounds_(this->getSPHSystem().system_domain_bounds_)
		, use_split_cell_lists_(false)
		, particle_sorting_(this)
		, iteration_count_(1)
	{
		sph_system.real_bodies_.push_back(this);
		size_t number_of_split_cell_lists = powerN(3, Dimensions);
		split_cell_lists_.resize(number_of_split_cell_lists);
		cell_linked_list_ = cell_linked_list_keeper_.movePtr( sph_adaptation_->createCellLinkedList(system_domain_bounds_, *this));
	}
	//=================================================================================================//
	void RealBody::assignBaseParticles(BaseParticles *base_particles)
	{
		SPHBody::assignBaseParticles(base_particles);
		particle_sorting_.assignBaseParticles(base_particles);
		cell_linked_list_->assignBaseParticles(base_particles);
	}
	//=================================================================================================//
	void RealBody::sortParticleWithCellLinkedList()
	{
		StdLargeVec<size_t> &sequence = base_particles_->sequence_;
		size_t size = base_particles_->total_real_particles_;
		cell_linked_list_->computingSequence(sequence);
		particle_sorting_.sortingParticleData(sequence.data(), size);
	}
	//=================================================================================================//
	void RealBody::updateCellLinkedList()
	{
		cell_linked_list_->UpdateCellLists();
		base_particles_->total_ghost_particles_ = 0;
	}
	//=================================================================================================//
	void RealBody::updateCellLinkedListWithParticleSort(size_t particle_sorting_period)
	{
		if (iteration_count_ % particle_sorting_period == 0)
			sortParticleWithCellLinkedList();
		iteration_count_++;
		updateCellLinkedList();
	}
	//=================================================================================================//
	void RealBody::defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
	{
		sph_adaptation_->resetAdaptationRatios(h_spacing_ratio, new_system_refinement_ratio);
		cell_linked_list_ = cell_linked_list_keeper_.movePtr(
			sph_adaptation_->createCellLinkedList(system_domain_bounds_, *this));
	}
	//=================================================================================================//
}
