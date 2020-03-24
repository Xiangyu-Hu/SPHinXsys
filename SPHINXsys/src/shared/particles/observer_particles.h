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

		// curently nothing here yet
	};
	/**
	 * @class ObserverParticles
	 * @brief A group of particles for observing.
	 */
	class ObserverParticles : public BaseParticles
	{
	public:
		explicit ObserverParticles(SPHBody *body);
		virtual ~ObserverParticles() {};

		/** Vector of observer body data. */
		StdLargeVec<ObserverParticleData> observer_body_data_; 

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override {};
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override {};

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override {};
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override {};

		/** Pointer to this object. */
		virtual ObserverParticles* PointToThisObject();
	};
}
