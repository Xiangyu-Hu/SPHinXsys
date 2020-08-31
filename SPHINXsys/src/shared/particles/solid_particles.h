/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	solid_particles.h
 * @brief 	This is the derived class of base particle.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 * @version 0.2.1
 * 			add muscle particles and muscle data.
 */
#pragma once

#include "base_particles.h"
#include "xml_engine.h"
#include <fstream>

using namespace std;
namespace SPH {

	/** preclaimed classes. */
	class Solid;
	class ElasticSolid;
	class ActiveMuscle;
	class BaseMaterial;

	/**
	 * @class SolidParticleData 
	 * @brief Data for solid body particles.
	 */
	class SolidParticleData 
	{
	public:
		/** default constructor */
		SolidParticleData();
		/** in constructor, set the particle at rest*/
		SolidParticleData(BaseParticleData& base_particle_data, Solid* solid);
		virtual ~SolidParticleData() {};

		/** mass, reference density and current density. */
		Real rho_0_, rho_n_, mass_;
		/** Inital position, and inital and current normal direction. */
		Vecd n_0_, n_;
		/** Linear reproducing configuration correction. */
		Matd B_;
	
		/** fluid time-step averaged particle velocity and acceleration,
		  * or applying fluid structure interaction. */
		Vecd vel_ave_, dvel_dt_ave_;	
		/** Forces from fluid. */
		Vecd force_from_fluid_, viscous_force_from_fluid_;
	};

	/**
	 * @class ElasticSolidParticleData 
	 * @brief Data for elastic solid body particles.
	 */	
	class ElasticSolidParticleData 
	{
	public:
		ElasticSolidParticleData(BaseParticleData &base_particle_data,
			ElasticSolid *elastic_solid);
		virtual ~ElasticSolidParticleData() {};

		/** elastic body strain variables, deformation tensor. */
		Matd F_, dF_dt_, stress_;

		/** temporally particle position for computing average velocity. */
		Vecd pos_temp_;				
	};
	/**
	 * @class SolidParticles
	 * @brief A group of particles with solid body particle data.
	 */
	class SolidParticles : public BaseParticles
	{
	public:
		/** Constructor as the most derived object. */
		SolidParticles(SPHBody* body);
		/** Constructor */
		SolidParticles(SPHBody* body, Solid* solid);
		virtual ~SolidParticles() {};

		/** Vector of solid body data. */
		StdLargeVec<SolidParticleData> solid_body_data_; 
	
		/** Set initial condition for a solid body with different material. */
		virtual void OffsetInitialParticlePosition(Vecd offset);
		/** add buffer particles which latter may be realized for particle dynamics*/
		virtual void AddABufferParticle() override;
		/** Copy state from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override;

		/**
		 * @brief Write particle data in VTU format for Paraview.
		 * @param[inout] output_file Ofstream of particle data.
		 */

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in XML format.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) {};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[inout] filefullpath Full path to file being write.
		 */

		/* Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/**
		 * @brief Initialize particle data from restart xml file.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override ;
		/**
		 * @brief Reload particle position and volume from XML files.
		 * @param[inout] filefullpath Full path to file being write.
		 */	
		virtual void ReadFromXmlForReloadParticle(std::string &filefullpath) override;
		/** Pointer to this object. */
		virtual SolidParticles* PointToThisObject() override;
		/** Normalize a gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& gradient) override;
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) override;
	};
	
	/**
	 * @class ElasticSolidParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticSolidParticles : public SolidParticles
	{
	protected:
		/** Computing von_Mises_stress. */
		Real von_Mises_stress(size_t particle_i);
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		ElasticSolidParticles(SPHBody* body, ElasticSolid* elastic_solid);
		/**
		 * @brief Destructor.
		 */
		virtual ~ElasticSolidParticles() {};

		/** Vector of elastic solid particle data. */
		StdLargeVec<ElasticSolidParticleData> elastic_body_data_;

		/** Copy state from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override;

		/**
		 * @brief Write particle data in VTU format for Paraview.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;
		/**
		 * @brief Write particle data in XML format.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/**
		 * @brief Initialize particle data from restart xml file.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override ;

		/** Pointer to this object.  */
		virtual ElasticSolidParticles* PointToThisObject() override;
	};

	/**
	 * @class ActiveMuscleParticleData
	 * @brief Data for active muscle.
	 */
	class ActiveMuscleParticleData
	{
	public:
		/** default constructor. */
		ActiveMuscleParticleData() 
			: active_contraction_stress_(0.0), active_stress_(0.0) {};
		virtual ~ActiveMuscleParticleData() {};

		/** Active contraction stress. */
		Real active_contraction_stress_;
		Matd active_stress_;
	};

	/**
	 * @class ActiveMuscleParticles
	 * @brief A group of particles with active muscle particle data.
	 */
	class ActiveMuscleParticles : public ElasticSolidParticles
	{
	public:
		/** including of electrophysiology data_. */
		StdLargeVec<ActiveMuscleParticleData> active_muscle_data_;

		/** Constructor. */
		ActiveMuscleParticles(SPHBody* body, ActiveMuscle* active_muscle);
		/** Default destructor. */
		virtual ~ActiveMuscleParticles() {};

		/** Copy state from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override;

		/**
		 * @brief Write particle data in VTU format for Paraview.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToVtuFile(ofstream& output_file) override;
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream& output_file) override;
		/**
		 * @brief Write particle data in XML format.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string& filefullpath) override {};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string& filefullpath) override;
		/**
		 * @brief Initialize particle data from restart xml file.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void ReadParticleFromXmlForRestart(std::string& filefullpath) override;

		/** Pointer to this object.  */
		virtual ActiveMuscleParticles* PointToThisObject() override;

		/** Access a real data*/
		Real getActiveContractionStress(size_t particle_index) {
			return  active_muscle_data_[particle_index].active_contraction_stress_;
		}
		/** Access a matrix data*/
		Matd getActiveStress(size_t particel_index) {
			return active_muscle_data_[particel_index].active_stress_;
		};
	};
}