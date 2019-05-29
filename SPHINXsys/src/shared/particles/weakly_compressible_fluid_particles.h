/**
 * @file 	weakly_compressible_fluid_particle.h
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
	 * @brief Related material.
	 */
	class WeaklyCompressibleFluid;

	/**
	 * @class WeaklyCompressibleFluidParticleData 
	 * @brief Data for weakly compressible flud body particles.
	 */
	class WeaklyCompressibleFluidParticleData 
	{
	public:
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 */
		WeaklyCompressibleFluidParticleData();
		/**
		 * @brief Default destructor.
		 */
		virtual ~WeaklyCompressibleFluidParticleData() {};

		Real mass_, sigma_0_, rho_0_, rho_n_;	/**< Particle mass, initial number desity, initial density and current density. */
		Real rho_sum_, p_, c_;					/**< Particle denisty from summation, pressure. */
		Real drho_dt_, div_correction_;			/**< Paticle desity change rate and divergence correction. */
		Vecd dvel_dt_trans_, vel_trans_;		/**< Paticle transport acceleration and velocity. */

		Vec3d vorticity_;					/**< Vorticcity of fluid in 3D. */
		Real vort_2d_; 						/**< Vorticcity of fluid in 2D. */

		Real temp_real_, temp_real_1_;		/**< Temporary data for initermediate usage. */
	protected:

	};

	/**
	 * @class WeaklyCompressibleFluidParticleData
	 * @brief Weakly compressible flud body particles.
	 */
	class WeaklyCompressibleFluidParticles : public Particles
	{
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		explicit WeaklyCompressibleFluidParticles(string body_name);
		/**
		 * @brief Default destructor.
		 */
		virtual ~WeaklyCompressibleFluidParticles() {};

		StdLargeVec<WeaklyCompressibleFluidParticleData> fluid_data_; 	/**< vector of fluid particle data. */
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
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/**
		 * @brief Initialize particle data from restart xml file.
		 */
		virtual void InitialParticleFromRestartXmlFile(std::string &filefullpath) override ;
		/**
		 * @brief Pointer to this object. 
		 */
		virtual WeaklyCompressibleFluidParticles* PointToThisObject();

	protected:

	};

	/**
	 * @class Oldroyd_B_FluidParticleData
	 * @brief Data for Oldroyd B type non-Newtanian flud body particles.
	 */
	class Oldroyd_B_FluidParticleData
	{
	public: 
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 */
		Oldroyd_B_FluidParticleData();
		/**
		 * @brief Default destructor.
		 */
		virtual ~Oldroyd_B_FluidParticleData() {};
		Matd tau_, dtau_dt_;	/**< Particle elastic stress. */
	};

	/**
	 * @class Oldroyd_B_FluidParticles
	 * @brief Oldroyd_B flud body particles.
	 */	class Oldroyd_B_FluidParticles : public WeaklyCompressibleFluidParticles
	{
	public:
		//constructor
		explicit Oldroyd_B_FluidParticles(string body_name);
		/**
		 * @brief Default destructor.
		 */
		virtual ~Oldroyd_B_FluidParticles() {};

		StdLargeVec<Oldroyd_B_FluidParticleData> oldroyd_b_data_;	/**< Vector of oldroyd b particle data. */
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
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) override{};
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/**
		 * @brief Initialize particle data from restart xml file.
		 */
		virtual void InitialParticleFromRestartXmlFile(std::string &filefullpath) override ;
		/**
		 * @brief Pointer to this object. 
		 */
		virtual Oldroyd_B_FluidParticles* PointToThisObject();

	};
}
