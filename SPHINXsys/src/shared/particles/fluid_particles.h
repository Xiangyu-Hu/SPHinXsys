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
	/**
	 * @class FluidParticleData 
	 * @brief Data for newtonian fluid particles.
	 */
	class FluidParticleData 
	{
	public:
		FluidParticleData();
		virtual ~FluidParticleData() {};

		/** Particle mass, initial number desity, initial density and current density. */
		Real mass_, sigma_0_, rho_0_, rho_n_;	
		/** Particle denisty from summation, pressure and sound speed. */
		Real rho_sum_, p_, c_;
		/** Paticle desity change rate and divergence correction. */
		Real drho_dt_, div_correction_;
		/** Paticle transport acceleration and velocity. */
		Vecd dvel_dt_trans_, vel_trans_;		
		/** Vorticcity of fluid in 3D. */
		Vec3d vorticity_;					
		/** Vorticcity of fluid in 2D. */
		Real vort_2d_; 						

		/** Temporary data for initermediate usage. */
		Real temp_real_;		
	};

	/**
	 * @class FluidParticles
	 * @brief newtonian flud particles.
	 */
	class FluidParticles : public Particles
	{
	public:
		explicit FluidParticles(SPHBody *body);
		virtual ~FluidParticles() {};

		/** vector of fluid particle data. */
		StdLargeVec<FluidParticleData> fluid_particle_data_; 	

		//----------------------------------------------------------------------
		//Global data
		//----------------------------------------------------------------------
		/** Maximum signal speed.*/
		Real signal_speed_max_;

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
		virtual FluidParticles* PointToThisObject();
	};

	/**
	 * @class ViscoelasticFluidParticleData
	 * @brief Data for viscoelastic non-Newtanian flud particles.
	 */
	class ViscoelasticFluidParticleData
	{
	public: 
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
		explicit ViscoelasticFluidParticles(SPHBody *body);
		virtual ~ViscoelasticFluidParticles() {};
		
		/** Vector of oldroyd b particle data. */
		StdLargeVec<ViscoelasticFluidParticleData> viscoelastic_particle_data_;	

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override;

		/** Pointer to this object. */
		virtual ViscoelasticFluidParticles* PointToThisObject();
	};
}
