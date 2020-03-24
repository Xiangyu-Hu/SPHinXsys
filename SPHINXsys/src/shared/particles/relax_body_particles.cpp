/**
 * @file relax_body_particles.cpp
 * @brief Definition of funcitons declared in relax_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "relax_body_particles.h"
#include "base_material.h"
//=============================================================================================//
namespace SPH 
{
	//=============================================================================================//
	RelaxBodyParticleData::RelaxBodyParticleData(Vecd position)
	:pos_0_(position),f_(0.0), s_(0.0)
	{
		/** Nothing to do rightnow. */
	}
	//=============================================================================================//
	RelaxBodyParticles::RelaxBodyParticles(SPHBody *body)
		: BaseParticles(body, new BaseMaterial())
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i) 
		{
			Point pnt = base_particle_data_[i].pos_n_;
			relax_body_data_.push_back(RelaxBodyParticleData(pnt));
		}
	}
	//=============================================================================================//
	RelaxBodyParticles::RelaxBodyParticles(SPHBody* body, BaseMaterial* base_material)
		: BaseParticles(body, base_material)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i) 
		{
			Point pnt = base_particle_data_[i].pos_n_;
			relax_body_data_.push_back(RelaxBodyParticleData(pnt));
		}
	}
	//=============================================================================================//
	void RelaxBodyParticles
		::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		BaseParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		relax_body_data_[this_particle_index] = relax_body_data_[another_particle_index];
	}
	//=============================================================================================//
	void RelaxBodyParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		BaseParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(relax_body_data_[this_particle_index], relax_body_data_[that_particle_index]);
	}
	//=============================================================================================//
	RelaxBodyParticles* RelaxBodyParticles::PointToThisObject()
	{
		return this;
	}
	//=============================================================================================//
	void RelaxBodyParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		BaseParticles::WriteParticlesToVtuFile(output_file);
	}
	//=============================================================================================//
	void RelaxBodyParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function RelaxBodyParticles::WriteParticlesToXmlForRestart is not done. Exit the program! \n";
		exit(0);

	}
	//=============================================================================================//
	void RelaxBodyParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function RelaxBodyParticles::WriteParticlesToXmlForRestart is not done. Exit the program! \n";
		exit(0);

	}
	//=============================================================================================//
}
//=============================================================================================//