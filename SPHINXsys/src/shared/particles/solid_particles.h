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
 * @brief 	This is the derived class of base particles.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particles.h"
#include "xml_engine.h"
#include <fstream>

using namespace std;
namespace SPH {

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class Solid;
	class ElasticSolid;
	class ActiveMuscle;
	class BaseMaterial;

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

		StdLargeVec<Vecd>	pos_0_;	/**< initial position */
		StdLargeVec<Vecd>	n_;		/**<  current normal direction */
		StdLargeVec<Vecd>	n_0_;	/**<  inital normal direction */  //seems to be moved to method class
		StdLargeVec<Matd>	B_;		/**<  configuration correction for linear reproducing */

		//----------------------------------------------------------------------
		//		for fluid-structure interaction (FSI) 
		//----------------------------------------------------------------------
		StdLargeVec<Vecd>	vel_ave_;	/**<  fluid time-step averaged particle velocity */
		StdLargeVec<Vecd>	dvel_dt_ave_;	/**<  fluid time-step averaged particle acceleration */
		StdLargeVec<Vecd>	force_from_fluid_;	/**<  forces (including pressure and viscous) from fluid */
		StdLargeVec<Vecd>	viscous_force_from_fluid_;	/**<  viscous forces from fluid */

		//----------------------------------------------------------------------
		//		for soild-soild contact dynmaics 
		//----------------------------------------------------------------------
		StdLargeVec<Real> 	contact_density_;		/**< density due to contact of solid-solid. */
		StdLargeVec<Vecd>	contact_force_;			/**< contact force from other solid body or bodies */

		/** shift the initial position of the solid particles. */
		void OffsetInitialParticlePosition(Vecd offset);
		/** initialize normal direction along solid body shape. */
		void initializeNormalDirectionFromGeometry();
		/** particle translation and rotation */
		void ParticleTranslationAndRotation(Transformd& transform);
		/** Write particle data in PLT format for Tecplot */
		virtual void writeParticlesToPltFile(ofstream &output_file) override;
		/** Write particle data in XML format for restart */
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file */
		virtual void readParticleFromXmlForRestart(std::string &filefullpath) override ;
		/** Reload particle position and volume from XML files */	
		virtual void readFromXmlForReloadParticle(std::string &filefullpath) override;

		/** Pointer to this object. */
		virtual SolidParticles* pointToThisObject() override;

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

		StdLargeVec<Matd>	F_;			/**<  deformation tensor */
		StdLargeVec<Matd>	dF_dt_;		/**<  deformation tensor change rate */
		StdLargeVec<Matd>	stress_;	/**<  stress tensor */

		/** Write particle data in VTU format for Paraview */
		virtual void writeParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot */
		virtual void writeParticlesToPltFile(ofstream &output_file) override;
		/** Write particle data in XML format for restart */
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file */
		virtual void readParticleFromXmlForRestart(std::string &filefullpath) override ;

		/** Pointer to this object.  */
		virtual ElasticSolidParticles* pointToThisObject() override;
	};

	/**
	 * @class ActiveMuscleParticles
	 * @brief A group of particles with active muscle particle data.
	 */
	class ActiveMuscleParticles : public ElasticSolidParticles
	{
	public:

		StdLargeVec<Real>	active_contraction_stress_;			/**<  active contraction stress */
		StdLargeVec<Matd>	active_stress_;		/**<  active stress */ //seems to be moved to method class

		/** Constructor. */
		ActiveMuscleParticles(SPHBody* body, ActiveMuscle* active_muscle);
		/** Default destructor. */
		virtual ~ActiveMuscleParticles() {};

		/** Write particle data in VTU format for Paraview */
		virtual void writeParticlesToVtuFile(ofstream& output_file) override;
		/** Write particle data in PLT format for Tecplot */
		virtual void writeParticlesToPltFile(ofstream& output_file) override;
		/** Write particle data in XML format for restart */
		virtual void writeParticlesToXmlForRestart(std::string& filefullpath) override;
		/** Initialize particle data from restart xml file */
		virtual void readParticleFromXmlForRestart(std::string& filefullpath) override;

		/** Pointer to this object.  */
		virtual ActiveMuscleParticles* pointToThisObject() override;
	};
}
