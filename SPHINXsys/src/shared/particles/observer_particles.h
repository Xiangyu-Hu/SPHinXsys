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
		/**
		 * @brief Default constrcutor.
		 */
		ObserverParticleData();
		/**
		 * @brief Default destructor.
		 */
		virtual ~ObserverParticleData() {};
	};
	/**
	 * @class ObserverParticles
	 * @brief A group of particles with observer body particle data.
	 */
	class ObserverParticles : public Particles
	{
	public:
		/**
		 * @brief Default constrcutor.
		 */
		explicit ObserverParticles(string body_name);
		/**
		 * @brief Default destructor.
		 */
		virtual ~ObserverParticles() {};

		StdLargeVec<ObserverParticleData> observer_body_data_; /**< Vector of observer body data. */
		/**
		 * @brief Initialize a prticle by input a postion and volume. 
		 * @param[in] pnt Vecotor of particle position.
		 * @param[in] particle_volume Volume of particle.
		 */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) override;
		/**
		 * @brief Write particle data in VTU format for Paraview.
		 * @param[out] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override {};
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[out] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override {};
		/**
		 * @brief Write particle data in XML format.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override{};
		/**
		 * @brief Initialize particle data from restart xml file.
		 */
		virtual void InitialParticleFromRestartXmlFile(std::string &filefullpath) override {};
		/**
		 * @brief Pointer to this object. 
		 */
		virtual ObserverParticles* PointToThisObject();
	};
}
