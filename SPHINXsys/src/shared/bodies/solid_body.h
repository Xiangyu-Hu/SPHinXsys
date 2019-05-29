/**
 * @file    solid_body.h
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#pragma once

#include "base_body.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief Friend Class.
	 */
	class SPHSystem;
	class ElasticSolid;
	class ExternalForce;
	class MeshCellLinkedList;
	class SolidBodyParticles;
	class ElasticBodyParticles;
	class Muscle;
	class MuscleBodyParticles;

	/**
	 * @class SolidBody
	 * @brief Declaration of solidbody which is used for Solid BCs and derived from RealBody.
	 */
	class SolidBody : public RealBody
	{
	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] solid_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		SolidBody(SPHSystem &system, string body_name, 
			SolidBodyParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default destructor.
		 */
		virtual ~SolidBody() {};

		SolidBodyParticles &solid_particles_;	/**< Particles in this body. */

		/**
		 * @brief Build inner configuration.
		 */
		virtual void BuildInnerConfiguration() override;
		/**
		 * @brief Build contact configuration.
		 */
		virtual void BuildContactConfiguration() override;
		/**
		 * @brief Initial condition.
		 */
		virtual void InitialCondition() = 0;
		/**
		 * @brief Set a particle at rest for easy initial condition.
		 */
		virtual void SetAllParticleAtRest() override;
		/**
		 * @brief offset initial particle position.
		 */
		virtual void OffsetInitialParticlePosition(Vecd offset) override;
		/**
		 * @brief Override the virtual function to output global basic parameters in SolidBody.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	private:

	};

	/**
	 * @class SolidBodyPart
	 * @brief A auxillariy class for SoildBody to
	 * indicate a part of the body moving together with particles.
	 */
	class SolidBodyPart : public LagrangianBodyPart
	{
	protected:
		SolidBody *solid_body_;
		Region soild_body_part_region_;

		virtual void TagBodyPartParticles() override;
	public:
		SolidBodyPart(SolidBody *solid_body, string soild_body_part_name);
		virtual~SolidBodyPart() {};
	};
	
	/**
	 * @class SolidBodyPartForSimbody
	 * @brief A SolidBodyPart for coupling with Simbody.
	 * The mass, origin, and unit inertial matrix are computed
	 * Note: Simbody spatial vectors are three dimensional
	 */
	class SolidBodyPartForSimbody : public SolidBodyPart
	{
	protected:
		Real solid_body_density_;


		virtual void TagBodyPartParticles() override;
	public:
		SolidBodyPartForSimbody(SolidBody *solid_body, 
			string soild_body_part_name, Real solid_body_density);
		virtual~SolidBodyPartForSimbody() {};

		SimTK::MassProperties *body_part_mass_properties_;
		Vec3 initial_mass_center_;
	};

	/**
	 * @class ElasticBody
	 * @brief Declaration of elasticbody which is used for Elastic structure, 
	 * and derived from SolidBody.
	 */
	class ElasticBody : public SolidBody
	{
	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] material Material property of this body.
		 * @param[in] elastic_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		ElasticBody(SPHSystem &system, string body_name, 
			ElasticSolid* material,
			ElasticBodyParticles &elastic_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default destructor.
		 */
		virtual ~ElasticBody() {};
 
		ElasticSolid* material_;					/**< Functions on material properties. */
		ElasticBodyParticles &elastic_particles_;	/**< Particles in this body.*/

		/**
		 * @brief Set a particle at rest for easy initial condition.
		 */
		virtual void SetAllParticleAtRest() override;
		/**
		 * @brief Interface for particle dependent material properties.
		 */
		virtual void InitializeLocalMaterialProperties();
		/**
		 * @brief Calculate the elastic stree.
		 * @param[in] deform_grad Deform gradient.
		 * @param[in] index_particle_i Index of particle.
		 */		
		virtual Matd GetElasticStress(Matd &deform_grad, size_t index_particle_i);
		/**
		 * @brief Override the virtual function to output global basic parameters in ElasticBody.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	};
	/**
	 * @class MuscleBody
	 * @brief Declaration of musclbody which is used for muscl structure, 
	 * and derived from ElasticBody.
	 */
	class MuscleBody : public ElasticBody
	{
	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] material Material property of this body.
		 * @param[in] muscle_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		MuscleBody(SPHSystem &system, string body_name, Muscle* material, 
			MuscleBodyParticles &muscle_particles, int refinement_level, ParticlesGeneratorOps op);
		virtual ~MuscleBody() {};

		Muscle* material_;						/**< Functions on material properties. */
		MuscleBodyParticles &muscle_particles_;	/**< Particles in this body. */
		/**
		 * @brief Interface for particle dependent material properties.
		 */
		virtual void InitializeLocalMaterialProperties();
		/**
		 * @brief Calculate the elastic stree.
		 * @param[in] deform_grad Deform gradient.
		 * @param[in] index_particle_i Index of particle.
		 */	
		virtual Matd GetElasticStress(Matd &deform_grad, size_t index_particle_i) override;
		/**
		 * @brief Override the virtual function to output global basic parameters in MusclBody.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	};
}