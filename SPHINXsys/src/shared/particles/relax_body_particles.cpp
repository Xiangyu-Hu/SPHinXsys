/**
 * @file relax_body_particles.cpp
 * @brief Definition of funcitons declared in relax_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "relax_body_particles.h"

namespace SPH {
	//===============================================================//
	RelaxBodyParticleData::RelaxBodyParticleData(Vecd position)
	:pos_0_(position)
	{

	}
	//===============================================================//
	RelaxBodyParticles::RelaxBodyParticles(SPHBody *body)
		: Particles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			Point pnt = base_particle_data_[i].pos_n_;
			relax_body_data_.push_back(RelaxBodyParticleData(pnt));
		}
	}
	//===============================================================//
	RelaxBodyParticles* RelaxBodyParticles::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	void RelaxBodyParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function RelaxBodyParticles::WriteParticlesToXmlForRestart is not done. Exit the program! \n";
		exit(0);

	}
	//===========================================================//
	void RelaxBodyParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function RelaxBodyParticles::WriteParticlesToXmlForRestart is not done. Exit the program! \n";
		exit(0);

	}
	//===============================================================//
}
