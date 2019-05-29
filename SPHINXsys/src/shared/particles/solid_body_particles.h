/**
 * @file 	solid_body_particle.h
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
	 * @class SolidBodyParticleData 
	 * @brief Data for solid body particles.
	 */
	class SolidBodyParticleData 
	{
	public:
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 * @param[in] positioni Particle positiion.
		 */
		SolidBodyParticleData(Vecd position);
		/**
		 * @brief Default destructor.
		 */
		virtual ~SolidBodyParticleData() {};

		Vecd pos_0_, n_0_, n_;			/**< Inital position, and inital and current normal direction. */
		/** wealkly compressible fluid time-step averaged particle velocity and acceleration,
		 * or applying fluid structure interaction. */
		Vecd vel_ave_, dvel_dt_ave_;	
		Vecd force_from_fluid_, viscous_force_from_fluid_;	/**< Force from fluid. */
		Real temp_real_;				/**< Temporary data for initermediate usage. */
		Vecd temp_vec_;					/**< Temporary velocity data for initermediate usage. */
	};

	/**
	 * @class ElasticBodyParticleData 
	 * @brief Data for elastic solid body particles.
	 */	
	class ElasticBodyParticleData 
	{
	public:
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 */
		ElasticBodyParticleData();
		/**
		 * @brief Default destructor.
		 */
		virtual ~ElasticBodyParticleData() {};

		Real mass_, rho_0_, rho_n_;		/**< mass, reference density and current density. */
		Matd B_, F_, dF_dt_, stress_;	/**< elastic body strain varaibles, deformation tensor. */
		Matd temp_matrix_;				/**< Temporary deformation tensor. */
		Real local_G_, local_lambda_, local_eta_, local_c_;	/**< Local material properties. */
		Vecd pos_temp_;				/**< temporally particle position for computing average velocity. */
	};
	/**
	 * @class MuscleBodyData
	 * @brief Data for muscl body particles.
	 */	
	class MuscleBodyData
	{
	public:
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 */
		MuscleBodyData();
		/**
		 * @brief Default destructor.
		 */		
		virtual ~MuscleBodyData() {};
		Vecd local_f0_, local_s0_;					/**< Fiber and sheet direction and related matrices. */
		Matd local_f0f0_, local_s0s0_, local_f0s0_; /**< Local fiber and sheet direction and related matrices. */
		Real local_a_[4], local_b_[4];				/**< Local elastic properties. */
	};
	/**
	 * @class SolidBodyParticles
	 * @brief A group of particles with solid body particle data.
	 */
	class SolidBodyParticles : public Particles
	{
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		SolidBodyParticles(string body_name);
		/**
		 * @brief Destructor.
		 */
		virtual ~SolidBodyParticles() {};

		StdLargeVec<SolidBodyParticleData> solid_body_data_; /**< Vector of solid body data. */
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
		virtual SolidBodyParticles* PointToThisObject() override;
	};
	
	/**
	 * @class ElasticBodyParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticBodyParticles : public SolidBodyParticles
	{
	protected:
		Real von_Mises_stress(size_t particle_i);

	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		ElasticBodyParticles(string body_name);
		/**
		 * @brief Destructor.
		 */
		virtual ~ElasticBodyParticles() {};

		StdLargeVec<ElasticBodyParticleData> elastic_body_data_;
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
		virtual ElasticBodyParticles* PointToThisObject() override;

	};
	
	/**
	 * @class MuscleBodyParticles
	 * @brief A group of particles with muscle body particle data.
	 */
	class MuscleBodyParticles : public ElasticBodyParticles
	{
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		MuscleBodyParticles(string body_name);
		/**
		 * @brief Destructor.
		 */
		virtual ~MuscleBodyParticles() {};

		StdLargeVec<MuscleBodyData> muscle_body_data_;
		/**
		 * @brief Initialize a prticle by input a postion and volume. 
		 * @param[in] pnt Vecotor of particle position.
		 * @param[in] particle_volume Volume of particle.
		 */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) override;
		/**
		 * @brief Pointer to this object. 
		 */
		virtual MuscleBodyParticles* PointToThisObject() override;

	};
}