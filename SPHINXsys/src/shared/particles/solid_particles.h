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

#ifndef SOLID_PARTICLES_H
#define SOLID_PARTICLES_H


#include "base_particles.h"
	
namespace SPH {

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class Solid;
	class ElasticSolid;
	class PlasticSolid;
	template<class MuscleType> class ActiveMuscle;
	class ShapeSurface;

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
		//----------------------------------------------------------------------
		//		for solid-solid contact dynamics 
		//----------------------------------------------------------------------
		StdLargeVec<Real> 	contact_density_;		/**< density due to contact of solid-solid. */
		StdLargeVec<Vecd>	contact_force_;			/**< contact force from other solid body or bodies */

		void offsetInitialParticlePosition(Vecd offset);
		void initializeNormalDirectionFromGeometry();
		void ParticleTranslationAndRotation(Transformd& transform);

		/** Normalize a gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& gradient) override;
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) override;

		virtual SolidParticles* ThisObjectPtr() override {return this;};
	};
	
	/**
	 * @class ElasticSolidParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticSolidParticles : public SolidParticles
	{
	protected:
		virtual void writePltFileHeader(std::ofstream& output_file);
		virtual void writePltFileParticleData(std::ofstream& output_file, size_t index_i);

	public:
		ElasticSolidParticles(SPHBody* body, ElasticSolid* elastic_solid);
		virtual ~ElasticSolidParticles() {};

		StdLargeVec<Matd>	F_;			/**<  deformation tensor */
		StdLargeVec<Matd>	dF_dt_;		/**<  deformation tensor change rate */
		StdLargeVec<Matd>	stress_PK1_;	/**<  first Piola-Kirchhoff stress tensor */

		/**< Computing von_Mises_stress. */
		Real von_Mises_stress(size_t particle_i);
		StdLargeVec<Real> getVonMisesStress();
		Real getMaxVonMisesStress();

		/**< Computing displacemnt. */
		Vecd displacement(size_t particle_i);
		StdLargeVec<Vecd> getDisplacement();

		/**< Computing normal vector. */
		Vecd normal (size_t particle_i);
		StdLargeVec<Vecd> getNormal();

		/**< Computing von Mises equivalent stress. */
		Real von_Mises_strain (size_t particle_i);
		StdLargeVec<Real> getVonMisesStrain();
		Real getMaxVonMisesStrain();

		virtual void writeParticlesToVtuFile(std::ostream &output_file) override;
		/** Write only surface particle data in VTU format for Paraview. */
		virtual void writeSurfaceParticlesToVtuFile(std::ofstream& output_file, ShapeSurface& surface_particles);
		virtual ElasticSolidParticles* ThisObjectPtr() override {return this;};
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

		template<class MuscleType>
		ActiveMuscleParticles(SPHBody* body, ActiveMuscle<MuscleType>* active_muscle) : 
			ElasticSolidParticles(body, active_muscle)
		{
			active_muscle->assignActiveMuscleParticles(this);
			initializeActiveMuscleParticleData();
		};
		virtual ~ActiveMuscleParticles() {};

		virtual ActiveMuscleParticles* ThisObjectPtr() override { return this; };
	private:
		void initializeActiveMuscleParticleData();
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

		StdLargeVec<Matd> transformation_matrix_;	/**< initial transformation matrix from global to local coordinates */
		StdLargeVec<Real> shell_thickness_;			/**< shell thickness */
		//----------------------------------------------------------------------
		//	extra generalized coordinates in global coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> pseudo_n_;		     /**< current pseudo-normal vector */
		StdLargeVec<Vecd> dpseudo_n_dt_;		 /**< pseudo-normal vector change rate */
		StdLargeVec<Vecd> dpseudo_n_d2t_;		 /**< pseudo-normal vector second order time derivation */
		//----------------------------------------------------------------------
		//	extra generalized coordinate and velocity in local coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> rotation_;		    /**< rotation angle of the initial normal respective to each axis */
		StdLargeVec<Vecd> angular_vel_;	        /**< angular velocity respective to each axis */
		StdLargeVec<Vecd> dangular_vel_dt_;	    /**< angular accelration of respective to each axis*/
		//----------------------------------------------------------------------
		//	extra deformation and deformation rate in local coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Matd> F_bending_;		    /**< bending deformation gradient	*/
		StdLargeVec<Matd> dF_bending_dt_;  	    /**< bending deformation gradient change rate	*/
		//----------------------------------------------------------------------
		//	extra stress for pair interaction in global coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> global_shear_stress_;	    /**< global shear stress */
		StdLargeVec<Matd> global_stress_;	/**<  global stress for pair interaction */
		StdLargeVec<Matd> global_moment_;	/**<  global bending moment for pair interaction */

		virtual ShellParticles* ThisObjectPtr() override {return this;};
	};
}
#endif //SOLID_PARTICLES_H
