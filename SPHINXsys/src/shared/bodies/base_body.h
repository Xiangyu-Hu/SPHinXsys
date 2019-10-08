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
	 * @brief preclaimed class.
	 */
	class SPHSystem;
	class Output;
	class Particles;
	class Material;
	class Reaction;
	class Kernel;
	class MeshCellLinkedList;
	class MeshBackground;

	/**
	 * @class ParticlesGeneratorOps
	 * @brief Serval manners are provied for particles generator.
	 * @details lattice : Generate partice from lattcie grid.
	 *			direct  : Input particle position and volume directly.
	 */
	enum class ParticlesGeneratorOps {lattice, direct};
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
		SPHSystem &sph_system_; 	/**< SPHSystem. */
		string body_name_; 			/**< name of this body */
		/** Ratio between smoothing length to particle spacing. */
		Real smoothinglength_ratio_;		
		/** Computing particle spacing from refinement level */
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
			int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op);
		virtual ~SPHBody() {};

		Kernel *kernel_; 		/**< sph kernel function specific to a SPHBody */
		int refinement_level_;	/**< refinement level of this body */
		/** Allocate memory for cell linked list. */
		virtual void AllocateMeoemryCellLinkedList() {};
		/** Allocate memory for back ground mesh. */
		virtual void AllocateMeoemryBackgroundMesh() {};
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() = 0;
		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() = 0;
		/** Build contact configuration. */
		virtual void BuildContactConfiguration() = 0;
		/** Update inner configuration. */
		virtual void UpdateInnerConfiguration() = 0;
		/** Update contact configuration. */
		virtual void UpdateContactConfiguration() = 0;
		/** Update interactiong configuration. */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) = 0;

		//----------------------------------------------------------------------
		//Global variables
		//----------------------------------------------------------------------
		Real particle_spacing_;						/**< Particle spacing of the body. */
		size_t number_of_particles_;				/**< Number of real particles of the body. */
		Particles *base_particles_;					/**< Base particles of this body. */
		Material *base_material_;					/**< Base material of this body. */
		Reaction *base_reaction_;					/**< Base reaction model of this body */
		MeshCellLinkedList *mesh_cell_linked_list_; /**< Cell linked mesh of this body. */
		MeshBackground *mesh_background_;			/**< Background mesh.*/

		/**
 		 * @brief particle by cells lists is for parallel splitting algorithm
		 * particles in each cell are collected together, 
		 * if two partiles each belongs two different cell entries,
		 * they have no interaction because they are too far.
		 */
		size_t number_of_by_cell_lists_;
		ByCellLists by_cell_lists_particle_indexes_;
		
		/** Reference inner configurations for totoal Lagrangian formulation. */
		ReferenceNeighborList reference_inner_configuration_;	
		/** current inner configurations for updated Lagrangian formulation. */
		NeighborList current_inner_configuration_;
		
		/**
		 * @brief Contact configurations
		 * @details Note that contact configuration only gives all topological relation to this body.
		 * The specific physical interaction, which may not involving all contact bodies,
		 * will be defined in the specific particle dynamics
		 */
		 /** Contact map: pointing to toplogically contacted bodies. **/
		SPHBodyContactMap contact_map_;
		/** Lists of particles has a ocnfiguration with particles in contaced bodies. **/
		ContactParticleList indexes_contact_particles_;
		/** Configurations for total Lagrangian formulation.**/
		ReferenceContactNeighborList reference_contact_configuration_;
		/** Configurations for updated Lagrangian formulation. **/
		ContactNeighborList current_contact_configuration_;				

		/** Get the name of this body for out file name. */
		string GetBodyName(); 
		/** Set up the contact map. */
		void SetContactMap(SPHBodyContactMap &contact_map);
		/** Allocate memories for configuration. */
		void AllocateMemoriesForConfiguration();

		/** the reagion describe the geometry of the body.
		 * static member, so the geoemtry head file is included. */
		Region body_region_;
		/** Check wether a point within the geometry of this body.
		 * @returns TRUE if a point within body's region otherwise FALSE. 
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
		/** Output particle data in VTU file for visuallization in Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file);
		/** Output particle data in PLT file for visuallization in Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file);

		/** Output particle data in XML file for restart simulation. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath);
		/** Read particle data in XML file for restart simulation. */
		virtual void ReadParticlesFromXmlForRestart(std::string &filefullpath);

		/** Output particle position and volume in XML file for reloading particles. */
		virtual void WriteToXmlForReloadParticle(std::string &filefullpath);
		/** Reload particle position and volume from XML files. */
		virtual void ReadFromXmlForReloadParticle(std::string &filefullpath);
		
		/** The pointer to derived class object. */
		virtual SPHBody* PointToThisObject();
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
			int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op);
		virtual ~RealBody() {};
		/** Allocate memory for cell linked list. */
		virtual void AllocateMeoemryCellLinkedList() override;
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() override;
		/** Update inner configuration. */
		virtual void UpdateInnerConfiguration() override;
		/** Update contact configuration. */
		virtual void UpdateContactConfiguration() override;
		/** Update interaction configuration, e.g., both ineer and contact. */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) override;

		/** The pointer to derived class object. */
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
			int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op);
		virtual ~FictitiousBody() {};

		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() override;
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() override;
		/** Update inner configuration. */
		virtual void UpdateInnerConfiguration() override;
		/** Update contact configuration. */
		virtual void UpdateContactConfiguration() override;
		/** Update interaction configuration, e.g., both ineer and contact. */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) override;

		/** The pointer to derived class object. */
		virtual FictitiousBody* PointToThisObject() override;
	};

	/**
	 * @class BodyPart
	 * @brief An auxillariy class for SPHBody to indicate a part of the body.
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
	 * @brief An auxillariy class for SPHBody to 
	 * indicate a part of the body moving together with particles.
	 */
	class LagrangianBodyPart : public BodyPart
	{
	protected:
		virtual void TagBodyPartParticles() = 0;
	public:
		/** Collection particle in this body part. */
		IndexVector body_part_particles_;

		LagrangianBodyPart(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~LagrangianBodyPart() {};
	};

	/**
	 * @class EulerianBodyPart
	 * @brief An auxillariy class for SPHBody to
	 * indicate a part of the body fixed in space.
	 */
	class EulerianBodyPart : public BodyPart
	{
	protected:
		virtual void TagBodyPartCells() = 0;
	public:
		/** Collection of cells to indicate the body part. */
		CellVector body_part_cells_;

		EulerianBodyPart(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~EulerianBodyPart() {};
	};
}
