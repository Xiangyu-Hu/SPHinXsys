/**
 * @file 	base_body.h
 * @brief 	This is the base classes of SPH bodies. The real body is for 
 *			that with cell linked list and the fivtitious one doesnot.     
 * 			Before the defination of the SPH bodies, the shapes with complex 
 *			geometries, i.e. those are produced by adavnced binary operation, 
 * 			such as intersection, should be produced first.
 * 			Then, all shapes used in body definition should be either contain 
 * 			or not contain each other. 
 *			Partial overlap between them are not premitted.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "geometry.h"

#include <string>
using namespace std;

namespace SPH {
	/**
	 * @brief Friend Class.
	 */
	class SPHSystem;
	class Particles;
	class Kernel;
	class MeshCellLinkedList;
	class MeshBackground;
	class Region;

	/**
	 * @class ParticlesGeneratorOps
	 * @brief Serval manners are provied for particles generator.
	 * @details lattice : Generate partice from lattcie grid.
	 *			direct  : Input particle position and volume directly.
	 *			relax   : Read relaxed particles from specific files.
	 *			restart : Read particles info from restart files.
	 */
	enum class ParticlesGeneratorOps {lattice, direct, relax, restart};
	/**
	 * @class SPHBody
	 * @brief SPHBody is a base body with basic data and functions.
	 *		  Its derived class can be a real fluid body, a real deformable solid body,
	 *        a static or moving solid body or a fictitious body.
	 * 		  Note that only real bodies have cell linked list.
	 */
	class SPHBody
	{
	protected:
		SPHSystem &sph_system_; 	/**< SPHSystme. */
		string body_name_; 			/**< name of this body */
		/**
		 * @brief Computing particle spacing from refinement level.
		 * @returns current particle sapcing calculated 
		 * as particle_spacing_ref_ * powern(2.0, refinement_level_).
		 */
		Real RefinementLevelToParticleSpacing();

	public:
		/**
		 * @brief Defaut constructor of SPHBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] base_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		explicit SPHBody(SPHSystem &sph_system, string body_name,
				Particles &base_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default distructor.
		 */
		virtual ~SPHBody() {};

		Kernel *kernel_; 		/**< sph kernel function specific to a SPHBody */
		int refinement_level_;	/**< refinement level of this body */
		int rst_step_; 			/**< Time step for restart simulaiton
		/**
		 * @brief Allocate memory for cell linked list.
		 */
		virtual void AllocateMeoemryCellLinkedList() {};
		/**
		 * @brief Allocate memory for back ground mesh.
		 */
		virtual void AllocateMeoemryBackgroundMesh() {};
		/**
		 * @brief Update cell linked list.
		 */
		virtual void UpdateCellLinkedList() = 0;
		/**
		 * @brief Build inner configuration.
		 */
		virtual void BuildInnerConfiguration() = 0;
		/**
		 * @brief Build contact configuration.
		 */
		virtual void BuildContactConfiguration() = 0;
		/**
		 * @brief Update inner configuration.
		 */
		virtual void UpdateInnerConfiguration() = 0;
		/**
		 * @brief Update contact configuration.
		 */
		virtual void UpdateContactConfiguration() = 0;
		/**
		 * @brief Update interactiong configuration, e.g., both inner and contact.
		 * @param[in] interacting_bodies Interacting bodies of this body. 
		 */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) = 0;

		//----------------------------------------------------------------------
		//Global variables
		//----------------------------------------------------------------------
		Real speed_max_;							/**< Maxium particle speed. */
		Real particle_spacing_;						/**< Particle spacing of the body. */
		size_t number_of_particles_;				/**< Number of particles of the body. */
		Particles &base_particles_;					/**< Base particles of this body. */
		MeshCellLinkedList *mesh_cell_linked_list_; /**< Cell linked mesh of this body. */
		/**
 		 * @brief particle lists by cells for parallel splitting algorithm
		 * particles in each cell are collected together, 
		 *if two partiles each belongs two different entries,
		 * they have no interaction because they are too far.
		 */
		ByCellLists by_cell_lists_particle_indexes_;
		

		/** 
		  * @brief Inner configurations for totoal Lagrangian formulation. 
		  */
		ReferenceNeighborList reference_inner_configuration_;	
		/**
		  * @brief Inner configurations for updated Lagrangian formulation.
		  */
		NeighborList current_inner_configuration_;
		
		/**
		 * @brief Contact configurations
		 * @details Note that contact configuration only gives all topological relation to this body.
		 * the specific physical interaction, which may not involving all contact bodies,
		 * will be defined in the specific particle dynmaics
		 */
		 /** Contact map: pointing to contacted bodies. **/
		SPHBodyContactMap contact_map_;
		/** Pointing to Base particles of all contatced bodies. **/
		StdVec<Particles *>  base_particles_contact_bodies_;
		/** Lists of particles has a ocnfiguration with particles in contaced bodies. **/
		ContactParticleList indexes_contact_particles_;
		/** Configurations for total Lagrangian formulation.**/
		ReferenceContactNeighborList reference_contact_configuration_;
		/** Configurations for updated Lagrangian formulation. **/
		ContactNeighborList current_contact_configuration_;				

		/**
		 * @brief Get the name of this body for out file name.
		 * @returns Name of the body.
		 */
		string GetBodyName(); 
		/**
		 * @brief Set up the contact map.
		 */
		void SetContactMap(SPHBodyContactMap &contact_map);
		/**
		 * @brief Allocate memories for configuration
		 */
		void AllocateMemoriesForConfiguration();

		/**
		 * @brief the reagion describe the geometry of the body
		 * static member, so the geoemtry head file is included.
		 */
		Region body_region_;
		/**
		 * @brief Check wether a point within the geometry of this body.
		 * @returns TRUE if a point within body's region other wise FALSE.
		 */
		bool BodyContain(Vecd pnt); 
		/**
		 * @brief Find closest point from a given point to body surface.
		 * @param[in] input_pnt The given point.
		 * @param[in,out] closest_pnt The found point.
		 * @param[in,out] phi The distance from given point to closest point. 
		 */
		void ClosestPointOnBodySurface(Vecd input_pnt, Vecd& closest_pnt, Real& phi);
		/**
		 * @brief Find the lower and upper bounds of the body.
		 * @param[in,out] lower_bound Lower bound of this body.
		 * @param[in,out] upper_bound Upper bound of this body.
		 */
		void BodyBounds(Vecd &lower_bound, Vecd &upper_bound);
		/**
		 * @brief The positions and volumes for genearting particles directly.
		 */		
		PositionsAndVolumes body_input_points_volumes_;
		ParticlesGeneratorOps particle_generator_op_;	/**< Particle generator manner */
		/**
		 * @brief Generate particles for the body in specific manner.
		 * @details Several manners are defind in the Body constructor
		 *			lattice Generate particles on lattice point;
		 *			direct  Input particles with pos and vol;
		 *    		relaxed Read particles info from relax scheme;
		 *			restart Read particles info from restart files.
		 */		
		virtual void CreateParticelsInSpecificManner();
		/**
		 * @brief Generate a particle for the body.
		 * @param[in] pnt The position the particle.
		 * @param[in] particle_volume The volume of the particle.
		 */
		virtual void GenerateAParticle(Vecd pnt, Real particle_volume = 0.0);
		/**
		 * @brief Output particle data in VTU file for visuallization in Paraview.
		 * @param[in,out] output_file Ofstream for output.
		 */
		virtual void WriteParticlesToVtuFile(ofstream &output_file);
		/**
		 * @brief Output particle data in PLT file for visuallization in Tecplot.
		 * @param[in,out] output_file Ofstream for output.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file);
		/**
		 * @brief Output particle data in XML file for visuallization.
		 * @param[in,out] output_file Ofstream for output.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath);
		/**
		 * @brief Output particle data in XML file for restart simulation.
		 * @param[in,out] output_file Ofstream for writing.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath);
		/**
		 * @brief The pointer to derived class object.
		 */
		virtual SPHBody* PointToThisObject();
		/**
		 * @brief Output global data in dat file for later reference.
		 * @param[in,out] output_file Ofstream for output.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) = 0;
		/**
		 * @brief The initial condition is readed from restart files.
		 */
		virtual void InitialConditionFromRestartFile();
		/**
		 * @brief offset initial particle position.
		 */
		virtual void OffsetInitialParticlePosition(Vecd offset);
	};
	/**
	 * @class RealBody
	 * @brief Derived class from SPHBody. 
	 * With inner particle configuration or inner interactions.
	 */
	class RealBody : public SPHBody
	{
	protected:

	public:
		/**
		 * @brief Defaut constructor of RealBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] base_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		RealBody(SPHSystem &sph_system, string body_name, 
				Particles &base_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default distructor.
		 */
		virtual ~RealBody() {};
		/**
		 * @brief Allocate memory for cell linked list.
		 */
		virtual void AllocateMeoemryCellLinkedList() override;
		/**
		 * @brief Update cell linked list.
		 */
		virtual void UpdateCellLinkedList() override;
		/**
		 * @brief Update inner configuration.
		 */
		virtual void UpdateInnerConfiguration() override;
		/**
		 * @brief Update contact configuration.
		 */
		virtual void UpdateContactConfiguration() override;
		/**
		 * @brief Update interaction configuration, e.g., both ineer and contact.
		 * @param[in] interacting_bodies Interacting bodies of this body. 
		 */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) override;
		/**
		 * @brief Initilized local material properties.
		 */
		virtual void InitializeLocalMaterialProperties() {};
		/**
		 * @brief Set initial condition for the body.
		 */
		virtual void InitialCondition() = 0;
		/**
		 * @brief Set all particle at rest for easy initial condition.
		 */
		virtual void SetAllParticleAtRest();
		/**
		 * @brief The pointer to derived class object.
		 */
		virtual RealBody* PointToThisObject() override;
	};

	/**
	 * @class FictitiousBody.
	 * @brief Derived class from SPHBody. 
	 * Without innner configuration or inner interaction.
	 */
	class FictitiousBody : public SPHBody
	{
	protected:

	public:
		/**
		 * @brief Defaut constructor of FictitiousBodyBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] base_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		FictitiousBody(SPHSystem &system, string body_name, 
					Particles &base_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default distructor.
		 */
		virtual ~FictitiousBody() {};

		/**
		 * @brief Build inner configuration.
		*/
		virtual void BuildInnerConfiguration() override;
		/**
		 * @brief Update cell linked list.
		 */
		virtual void UpdateCellLinkedList() override;
		/**
		 * @brief Update inner configuration.
		 * doing nothing because there is no inner configuration
		 */
		virtual void UpdateInnerConfiguration() override;
		/**
		 * @brief Update contact configuration.
		 */
		virtual void UpdateContactConfiguration() override;
		/**
		 * @brief Update interaction configuration, e.g., both ineer and contact.
		 * @param[in] interacting_bodies Interacting bodies of this body. 
		 */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) override;

		/**
		 * @brief The pointer to derived class object.
		 */
		virtual FictitiousBody* PointToThisObject() override;
		/**
		 * @brief Output global data in dat file for visuallization.
		 * @param[in,out] output_file Ofstream for output.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	};

	/**
	 * @class BodyPart
	 * @brief A auxillariy class for SPHBody to indicate a part of the body.
	 */
	class BodyPart
	{
	protected:
		SPHBody *body_;
		string body_part_name_;
	public:
		BodyPart(SPHBody *body, string body_part_name)
			: body_(body), body_part_name_(body_part_name) {};
		virtual ~BodyPart() {};
	};

	/**
	 * @class LagrangianBodyPart
	 * @brief A auxillariy class for SPHBody to 
	 * indicate a part of the body moving together with particles.
	 */
	class LagrangianBodyPart : public BodyPart
	{
	protected:
		virtual void TagBodyPartParticles() = 0;
	public:
		IndexVector body_part_particles_;

		LagrangianBodyPart(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~LagrangianBodyPart() {};
	};

	/**
	 * @class EulerianBodyPart
	 * @brief A auxillariy class for SPHBody to
	 * indicate a part of the body fixed in space.
	 */
	class EulerianBodyPart : public BodyPart
	{
	protected:
		virtual void TagBodyPartCells() = 0;
	public:
		CellVector body_part_cells_;

		EulerianBodyPart(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~EulerianBodyPart() {};
	};
}
