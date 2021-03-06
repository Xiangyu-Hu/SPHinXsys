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
 */
#pragma once

#include "base_particles.h"

using namespace std;
namespace SPH {

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class Solid;
	class ElasticSolid;
	class ActiveMuscle;

	/**
	 * @class SolidParticles
	 * @brief A group of particles with solid body particle data.
	 */
	class SolidParticles : public BaseParticles
	{
	public:
		SolidParticles(SPHBody* body);
		SolidParticles(SPHBody* body, Solid* solid);
		virtual ~SolidParticles() {};

		StdLargeVec<Vecd>	pos_0_;	/**< initial position */
		StdLargeVec<Vecd>	n_;		/**<  current normal direction */
		StdLargeVec<Vecd>	n_0_;	/**<  inital normal direction */
		StdLargeVec<Matd>	B_;		/**<  configuration correction for linear reproducing */
		//----------------------------------------------------------------------
		//		for fluid-structure interaction (FSI) 
		//----------------------------------------------------------------------
		StdLargeVec<Vecd>	vel_ave_;	/**<  fluid time-step averaged particle velocity */
		StdLargeVec<Vecd>	dvel_dt_ave_;	/**<  fluid time-step averaged particle acceleration */
		StdLargeVec<Vecd>	force_from_fluid_;	/**<  forces (including pressure and viscous) from fluid */
		StdLargeVec<Vecd>	viscous_force_from_fluid_;	/**<  viscous forces from fluid */
		//----------------------------------------------------------------------
		//		for solid-solid contact dynamics 
		//----------------------------------------------------------------------
		StdLargeVec<Real> 	contact_density_;		/**< density due to contact of solid-solid. */
		StdLargeVec<Vecd>	contact_force_;			/**< contact force from other solid body or bodies */

		void offsetInitialParticlePosition(Vecd offset);
		void initializeNormalDirectionFromGeometry();
		void ParticleTranslationAndRotation(Transformd& transform);
		virtual void writeParticlesToPltFile(ofstream &output_file) override;
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath) override;
		virtual void readParticleFromXmlForRestart(std::string &filefullpath) override;
		virtual void readFromXmlForReloadParticle(std::string &filefullpath) override;

		/** Normalize a gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& gradient) override;
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) override;

		virtual SolidParticles* pointToThisObject() override {return this;};
	};
	
	/**
	 * @class ElasticSolidParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticSolidParticles : public SolidParticles
	{
	protected:
		Real von_Mises_stress(size_t particle_i); /**< Computing von_Mises_stress. */
	public:
		ElasticSolidParticles(SPHBody* body, ElasticSolid* elastic_solid);
		virtual ~ElasticSolidParticles() {};

		StdLargeVec<Matd>	F_;			/**<  deformation tensor */
		StdLargeVec<Matd>	dF_dt_;		/**<  deformation tensor change rate */
		StdLargeVec<Matd>	stress_PK1_;	/**<  first Piola-Kirchhoff stress tensor */
		StdLargeVec<Matd>	corrected_stress_;	/**<  corrected stress for pair interaction */

		virtual void writeParticlesToVtuFile(ofstream &output_file) override;
		virtual void writeParticlesToPltFile(ofstream &output_file) override;
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath) override;
		virtual void readParticleFromXmlForRestart(std::string &filefullpath) override;
		virtual ElasticSolidParticles* pointToThisObject() override {return this;};
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

		ActiveMuscleParticles(SPHBody* body, ActiveMuscle* active_muscle);
		virtual ~ActiveMuscleParticles() {};

		virtual void writeParticlesToPltFile(ofstream& output_file) override;
		virtual void writeParticlesToXmlForRestart(std::string& filefullpath) override;
		virtual void readParticleFromXmlForRestart(std::string& filefullpath) override;
		virtual ActiveMuscleParticles* pointToThisObject() override {return this;};
	};

	/**
	 * @class ShellParticles
	 * @brief A group of particles with shell particle data.
	 */
	class ShellParticles : public ElasticSolidParticles
	{
	public:
		ShellParticles(SPHBody* body, ElasticSolid* elastic_solid, Real thickness);
		virtual ~ShellParticles() {};

		StdLargeVec<Real> shell_thickness_;		 /**< the thickness of the shell should be passed in by the main function */
		//----------------------------------------------------------------------
		//	extra generalized coordinates compared to elastic solid particles
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> pseudo_n_;		     /**< current pseudo-normal vector */
		StdLargeVec<Vecd> dpseudo_n_dt_;		 /**< pseudo-normal vector change rate */
		StdLargeVec<Vecd> dpseudo_n_d2t_;		 /**< pseudo-normal vector second order time derivation */
		//----------------------------------------------------------------------
		//	extra generalized velocity compared to elastic solid particles
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> rotation_;		    /**< rotation angle of the initial normal respective to each axis */
		StdLargeVec<Vecd> angular_vel_;	        /**< angular velocity respective to each axis */
		StdLargeVec<Vecd> dangular_vel_dt_;	    /**< angular accelration of respective to each axis*/
		//----------------------------------------------------------------------
		//	extra deformation and deformation rate compared to elastic solid particles
		//----------------------------------------------------------------------
		StdLargeVec<Matd> F_bending_;		    /**< bending deformation gradient	*/
		StdLargeVec<Matd> dF_bending_dt_;  	    /**< bending deformation gradient change rate	*/
		//----------------------------------------------------------------------
		//	extra stress compared to elastic solid particles
		//----------------------------------------------------------------------
		StdLargeVec<Matd> resultant_stress_;	/**< membrane stress */
		StdLargeVec<Matd> resultant_moment_;	/**< bending moment */

		StdLargeVec<Matd> transformation_matrix_;   /**< transformation matrix from initial global to initial local coordinates */


		virtual void writeParticlesToPltFile(ofstream &output_file) override;
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath) override;
		virtual void readParticleFromXmlForRestart(std::string &filefullpath) override;
		virtual ShellParticles* pointToThisObject() override {return this;};
	};
}
