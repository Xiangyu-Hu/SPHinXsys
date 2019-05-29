/**
 * @file 	relax_body_particle.h
 * @brief 	This is the derived class of base particle.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particles.h"
#include "xml_engine.h"

#include <fstream>
#include<atomic>

using namespace std;

namespace SPH {
	/**
	 * @class RelaxBodyParticleData
	 * @brief Data for relax body particles.
	 */
	class RelaxBodyParticleData
	{

	public:
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 * @param[in] positioni Particle positiion.
		 */
		RelaxBodyParticleData(Vecd position);
		/**
		 * @brief Default destructor.
		 */
		virtual ~RelaxBodyParticleData() {};

		Vecd pos_0_;	/**< Paticle postion after mapping particles to surface. */
	};
	/**
	 * @class RelaxBodyParticles
	 * @brief A group of particles with relax body particle data.
	 */
	class RelaxBodyParticles : public Particles
	{
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		RelaxBodyParticles(string body_name);
		/**
		 * @brief Destructor.
		 */
		virtual ~RelaxBodyParticles() {};

		StdLargeVec<RelaxBodyParticleData> relax_body_data_;	/**< Vector of particle data. */
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
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[out] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in XML format.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override;
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
		virtual RelaxBodyParticles* PointToThisObject() override;
	};
}