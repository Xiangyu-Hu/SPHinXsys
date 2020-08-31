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
 * @file 	fluid_particle.h
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

	class Fluid;
	class Oldroyd_B_Fluid;
	/**
	 * @class FluidParticleData 
	 * @brief Data for newtonian fluid particles.
	 */
	class FluidParticleData 
	{
	public:
		/** default constructor */
		FluidParticleData();
		/** in the constructor, particles is set at rest */
		FluidParticleData(BaseParticleData &base_particle_data, Fluid *fluid);
		virtual ~FluidParticleData() {};

		Real p_; /**< Particle pressure. */
		/** Particle mass, initial number density, initial density and current density. */
		Real rho_0_, rho_n_, mass_;
		Real drho_dt_; /**< Particle density change rate */
		Vec3d vorticity_;	/**< Vorticity of fluid in 3D. */				
	};

	/**
	 * @class FluidParticles
	 * @brief newtonian flud particles.
	 */
	class FluidParticles : public BaseParticles
	{
	public:
		explicit FluidParticles(SPHBody *body, Fluid *fluid);
		virtual ~FluidParticles() {};

		/** vector of fluid particle data. */
		StdLargeVec<FluidParticleData> fluid_particle_data_; 	

		//----------------------------------------------------------------------
		//Global data
		//----------------------------------------------------------------------
		/** add buffer particles which latter may be realized for particle dynamics*/
		virtual void AddABufferParticle() override;
		/** copy particle data from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Update the state of a particle from another particle */
		virtual void UpdateFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override;
		/** Get mirror a particle along an axis direaction. */
		virtual void mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction) override;

		/**
		 * @brief Write particle data in XML format.
		 * @param[inout] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath){};
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override;

		/** Pointer to this object. */
		virtual FluidParticles* PointToThisObject() override;
	};

	/**
	 * @class ViscoelasticFluidParticleData
	 * @brief Data for viscoelastic non-Newtanian flud particles.
	 */
	class ViscoelasticFluidParticleData
	{
	public: 
		/** in constructor, set the particle at rest*/
		ViscoelasticFluidParticleData();
		virtual ~ViscoelasticFluidParticleData() {};
		/** Particle elastic stress. */
		Matd tau_, dtau_dt_;	
	};

	/**
	 * @class ViscoelasticFluidParticles
	 * @brief Viscoelastic fluid particles.
	 */	
	class ViscoelasticFluidParticles : public FluidParticles
	{
	public:
		//constructor
		explicit ViscoelasticFluidParticles(SPHBody *body, Oldroyd_B_Fluid* oldroyd_b_fluid);
		virtual ~ViscoelasticFluidParticles() {};
		
		/** Vector of Oldroyd b particle data. */
		StdLargeVec<ViscoelasticFluidParticleData> viscoelastic_particle_data_;	

		/** add buffer particles which latter may be realized for particle dynamics*/
		virtual void AddABufferParticle() override;
		/** copy particle data from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index) override;
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override;

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override;

		/** Pointer to this object. */
		virtual ViscoelasticFluidParticles* PointToThisObject() override;
	};
}
