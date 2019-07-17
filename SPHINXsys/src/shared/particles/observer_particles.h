/**
 * @file 	observer_particles.h
 * @brief 	This is the calss for observer particle,  derived class of base particle.
 * @author	Chi Zhang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_particles.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @class ObserverParticleData
	 * @brief Data for observer body particles.
	 */
	class ObserverParticleData 
	{
	public:
		ObserverParticleData();
		virtual ~ObserverParticleData() {};
	};
	/**
	 * @class ObserverParticles
	 * @brief A group of particles with observer body particle data.
	 */
	class ObserverParticles : public Particles
	{
	public:
		explicit ObserverParticles(string body_name);
		virtual ~ObserverParticles() {};

		/** Vector of observer body data. */
		StdLargeVec<ObserverParticleData> observer_body_data_; 
		/** Initialize a prticle by input a postion and volume. */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) override;
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override {};
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override {};
		/** Write particle data in XML format. */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override {};
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override {};
		/** Pointer to this object. */
		virtual ObserverParticles* PointToThisObject();
	};
}
