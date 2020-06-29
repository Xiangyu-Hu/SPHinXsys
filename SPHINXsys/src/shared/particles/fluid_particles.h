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
		/** Particle mass, initial number desity, initial density and current density. */
		Real rho_0_, rho_n_, mass_;
		Real drho_dt_; /**< Paticle desity change rate */
		/** Paticle transport acceleration and velocity. */
		Vecd dvel_dt_trans_, vel_trans_;
		Vec3d vorticity_;	/**< Vorticcity of fluid in 3D. */				
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
		/** Maximum signal speed.*/
		Real signal_speed_max_;

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
	 * @brief Viscoelastic flud particles.
	 */	
	class ViscoelasticFluidParticles : public FluidParticles
	{
	public:
		//constructor
		explicit ViscoelasticFluidParticles(SPHBody *body, Oldroyd_B_Fluid* oldroyd_b_fluid);
		virtual ~ViscoelasticFluidParticles() {};
		
		/** Vector of oldroyd b particle data. */
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
