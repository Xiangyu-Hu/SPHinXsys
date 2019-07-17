/**
 * @file 	solid_particle.h
 * @brief 	This is the derived class of base particle.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particles.h"
#include "xml_engine.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @class SolidParticleData 
	 * @brief Data for solid body particles.
	 */
	class SolidParticleData 
	{
	public:
		SolidParticleData(Vecd position);
		virtual ~SolidParticleData() {};

		/** Inital position, and inital and current normal direction. */
		Vecd pos_0_, n_0_, n_;			
		/** Linear reproducing configuration correction. */
		Matd B_;
		/** fluid time-step averaged particle velocity and acceleration,
		  * or applying fluid structure interaction. */
		Vecd vel_ave_, dvel_dt_ave_;	
		/** Forces from fluid. */
		Vecd force_from_fluid_, viscous_force_from_fluid_;	

		/** Temporary data for intermediate usage. */
		Real temp_real_;
		/** Temporary vector data for intermediate usage. */
		Vecd temp_vec_;
		/** Temporary matrix data for intermediate usage. */
		Matd temp_matrix_;
	};

	/**
	 * @class ElasticSolidParticleData 
	 * @brief Data for elastic solid body particles.
	 */	
	class ElasticSolidParticleData 
	{
	public:
		ElasticSolidParticleData();
		virtual ~ElasticSolidParticleData() {};

		/** mass, reference density and current density. */
		Real mass_, rho_0_, rho_n_;	
		/** elastic body strain varaibles, deformation tensor. */
		Matd F_, dF_dt_, stress_;

		/** temporally particle position for computing average velocity. */
		Vecd pos_temp_;				
	};
	/**
	 * @class SolidParticles
	 * @brief A group of particles with solid body particle data.
	 */
	class SolidParticles : public Particles
	{
	public:
		SolidParticles(string body_name);
		virtual ~SolidParticles() {};

		/** Vector of solid body data. */
		StdLargeVec<SolidParticleData> solid_body_data_; 
		/** Initialize a prticle by input a postion and volume. */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) override;
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;
		/** Write particle data in XML format. */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/* Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override ;
		/** Reload particle position and volume from XML files. */
		virtual void ReadFromXmlForReloadParticle(std::string &filefullpath) override;
		/** Pointer to this object. */
		virtual SolidParticles* PointToThisObject() override;
	};
	
	/**
	 * @class ElasticSoildParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticSolidParticles : public SolidParticles
	{
	protected:
		/** Computing von_Mises_stress. */
		Real von_Mises_stress(size_t particle_i);

	public:
		ElasticSolidParticles(string body_name);
		virtual ~ElasticSolidParticles() {};

		/** Vector of elastic solid particle data. */
		StdLargeVec<ElasticSolidParticleData> elastic_body_data_;
		/** Initialize a prticle by input a postion and volume. */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) override;
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;
		/** Write particle data in XML format. */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override ;
		/** Pointer to this object.  */
		virtual ElasticSolidParticles* PointToThisObject() override;

	};
}