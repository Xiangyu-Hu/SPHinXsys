/**
 * @file 	base_body.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */
#include "base_body.h"

#include "sph_system.h"
#include "base_particles.h"
#include "base_particles.hpp"
#include "body_relation.h"

namespace SPH
{
	//=================================================================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, const std::string &body_name,
					 SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: sph_system_(sph_system), body_name_(body_name), newly_updated_(true),
		  body_domain_bounds_(0, 0), is_domain_bounds_determined_(false),
		  sph_adaptation_(sph_adaptation_ptr_keeper_.assignPtr(sph_adaptation_ptr)),
		  generative_structure_(nullptr)
	{
		sph_system_.addABody(this);
		sph_adaptation_->initialize(this);
	}
	//=================================================================================================//
	BoundingBox SPHBody::getSPHSystemBounds()
	{
		return sph_system_.system_domain_bounds_;
	}
	//=================================================================================================//
	std::string SPHBody::getBodyName()
	{
		return body_name_;
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
	void SPHBody::setBodyDomainBounds(BoundingBox body_domain_bounds)
	{
		body_domain_bounds_ = body_domain_bounds;
		is_domain_bounds_determined_ = true;
	};
	//=================================================================================================//
	BoundingBox SPHBody::getBodyDomainBounds()
	{
		if (!is_domain_bounds_determined_)
		{
			body_domain_bounds_ = body_shape_.findBounds();
			is_domain_bounds_determined_ = true;
		}
		return body_domain_bounds_;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToVtuFile(std::ostream &output_file)
	{
		base_particles_->writeParticlesToVtk(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToVtpFile(std::ostream &output_file)
	{
		base_particles_->writeParticlesToVtk(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface& surface_particles)
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
	RealBody::RealBody(SPHSystem &sph_system, const std::string &body_name,
					   SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: SPHBody(sph_system, body_name, sph_adaptation_ptr),
		  particle_sorting_(this)
	{
		sph_system.addARealBody(this);
		cell_linked_list_ = cell_linked_list_keeper_.movePtr(sph_adaptation_->createCellLinkedList());
		size_t number_of_split_cell_lists = powerN(3, Vecd(0).size());
		split_cell_lists_.resize(number_of_split_cell_lists);
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
	}
	//=================================================================================================//
	FictitiousBody::
		FictitiousBody(SPHSystem &system, const std::string &body_name,
					   SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: SPHBody(system, body_name, sph_adaptation_ptr)
	{
		system.addAFictitiousBody(this);
	}
	//=================================================================================================//
	void BodyPartByParticle::tagParticles(TaggingParticleMethod &tagging_particle_method)
	{
		for (size_t i = 0; i < base_particles_->total_real_particles_; ++i)
		{
			tagging_particle_method(i);
		}
	};
	//=================================================================================================//
	void BodyPartByCell::tagCells(TaggingCellMethod &tagging_cell_method)
	{
		cell_linked_list_->tagBodyPartByCell(body_part_cells_, tagging_cell_method);
	}
	//=================================================================================================//
	BodyRegionByParticle::
		BodyRegionByParticle(SPHBody &sph_body, const std::string &body_part_name, Shape &shape)
		: BodyPartByParticle(sph_body, body_part_name), body_part_shape_(shape)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BodyRegionByParticle::tagByContain, this, _1);
		tagParticles(tagging_particle_method);
	}
	//=================================================================================================//
	void BodyRegionByParticle::tagByContain(size_t particle_index)
	{
		if (body_part_shape_.checkContain(base_particles_->pos_n_[particle_index]))
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
		Real phi = sph_body_->body_shape_.findSignedDistance(base_particles_->pos_n_[particle_index]);
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
		Real distance = fabs(sph_body_->body_shape_.findSignedDistance(base_particles_->pos_n_[particle_index]));
		if (distance < thickness_threshold_)
		{
			body_part_particles_.push_back(particle_index);
		}
	}
	//=================================================================================================//
	BodyRegionByCell::BodyRegionByCell(RealBody &real_body, const std::string &body_part_name, Shape &shape)
		: BodyPartByCell(real_body, body_part_name), body_part_shape_(shape)
	{
		TaggingCellMethod tagging_cell_method = std::bind(&BodyRegionByCell::checkNotFar, this, _1, _2);
		tagCells(tagging_cell_method);
	};
	//=================================================================================================//
	bool BodyRegionByCell::checkNotFar(Vecd cell_position, Real threshold)
	{
		return body_part_shape_.checkNotFar(cell_position, threshold);
	}
	//=================================================================================================//
	NearShapeSurface::
		NearShapeSurface(RealBody &real_body, const std::string &body_part_name, Shape &shape)
		: BodyPartByCell(real_body, body_part_name),
		  level_set_shape_(
			  level_set_shape_keeper_.createRef<LevelSetShape>(&real_body, shape, true))
	{
		TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurface::checkNearSurface, this, _1, _2);
		tagCells(tagging_cell_method);
	}
	//=================================================================================================//
	NearShapeSurface::NearShapeSurface(RealBody &real_body)
		: BodyPartByCell(real_body, "NearShapeSurface"),
		  level_set_shape_(
			  DynamicCast<LevelSetShape>(this, *real_body.body_shape_.getShapeByName(real_body.getBodyName())))
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
	TreeTerminates::TreeTerminates(SPHBody &sph_body)
		: BodyPartByParticle(sph_body, "Leaves"),
		  tree_(*DynamicCast<GenerativeTree>(this, sph_body.generative_structure_))
	{
		for (const auto *branch : tree_.branches_)
		{
			if (branch->is_terminated_)
			{
				size_t particle_index = branch->inner_particles_.back();
				body_part_particles_.push_back(particle_index);
			}
		}
	}
	//=================================================================================================//
}
