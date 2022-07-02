/**
 * @file 	base_particle_generator.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.h"
#include "in_output.h"

namespace SPH
{
	//=================================================================================================//
	BaseParticleGenerator::BaseParticleGenerator(SPHBody &sph_body)
		: base_particles_(sph_body.base_particles_),
		pos_(base_particles_->pos_), unsorted_id_(base_particles_->unsorted_id_)
	{
		if (sph_body.base_particles_ == nullptr || sph_body.base_material_ == nullptr)
		{
			std::cout << "\n Error: Particles or Materials have not been defined yet!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void BaseParticleGenerator::initializePosition(const Vecd &position)
	{
		pos_.push_back(position);
		unsorted_id_.push_back(base_particles_->total_real_particles_);
		base_particles_->total_real_particles_ ++;
	}
	//=================================================================================================//
	ParticleGenerator::ParticleGenerator(SPHBody &sph_body)
		: BaseParticleGenerator(sph_body), Vol_(base_particles_->Vol_) {}
	//=================================================================================================//
	void ParticleGenerator::initializePositionAndVolume(const Vecd &position, Real volume)
	{
		initializePosition(position);
		Vol_.push_back(volume);
	}
	//=================================================================================================//
	SurfaceParticleGenerator::SurfaceParticleGenerator(SPHBody &sph_body)
		: ParticleGenerator(sph_body),
		  n_(*base_particles_->getVariableByName<Vecd>("NormalDirection")),
		  thickness_(*base_particles_->getVariableByName<Real>("Thickness"))
	{
		sph_body.sph_adaptation_->getKernel()->reduceOnce();
	}
	//=================================================================================================//
	void SurfaceParticleGenerator::initializeSurfaceProperties(const Vecd &surface_normal, Real thickness)
	{
		n_.push_back(surface_normal);
		thickness_.push_back(thickness);
	}
	//=================================================================================================//
	void ObserverParticleGenerator::initializeGeometricVariables()
	{
		for (size_t i = 0; i < positions_.size(); ++i)
		{
			initializePositionAndVolume(positions_[i], 0.0);
		}
	}
	//=================================================================================================//
	ParticleGeneratorReload::
		ParticleGeneratorReload(SPHBody &sph_body, InOutput &in_output, const std::string &reload_body_name)
		: ParticleGenerator(sph_body)
	{
		if (!fs::exists(in_output.reload_folder_))
		{
			std::cout << "\n Error: the particle reload folder:" << in_output.reload_folder_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		file_path_ = in_output.reload_folder_ + "/" + reload_body_name + "_rld.xml";
	}
	//=================================================================================================//
	void ParticleGeneratorReload::initializeGeometricVariables()
	{
		base_particles_->readFromXmlForReloadParticle(file_path_);
	}
	//=================================================================================================//
}
