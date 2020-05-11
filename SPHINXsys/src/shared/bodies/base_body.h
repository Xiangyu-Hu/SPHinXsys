/**
 * @file 	base_body.h
 * @brief 	This is the base classes of SPH bodies. The real body is for 
 *			that with cell linked list and the fictitious one does not.     
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
#include "neighbor_relation.h"
#include "geometry.h"
#include <string>
using namespace std;

namespace SPH 
{
	/**
	 * @brief preclaimed classes.
	 */
	class SPHSystem;
	class BaseParticles;
	class Kernel;
	class BaseMeshCellLinkedList;
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
		/** the reagion describe the geometry of the body.
		 * static member, so the geoemtry head file is included. */
		Region body_region_;
		/** smoothing length. */
		Real smoothinglength_;
		/** Computational domain bounds of the body for boundry conditions. */
		Vecd body_lower_bound_, body_upper_bound_;
		/** Whether the computational domain bound for this body is prescribed. */
		bool prescribed_body_bounds_;
		/** Computing particle spacing from refinement level. */
		Real RefinementLevelToParticleSpacing();

		/** Generate a kernel. */
		Kernel* GenerateAKernel(Real smoothing_lenght);
		/** Change kernel function specific for this body. */
		void ReplaceKernelFunction(Kernel* kernel);
	public:
		//----------------------------------------------------------------------
		//Global variables
		//----------------------------------------------------------------------
		int refinement_level_;	/**< refinement level of this body */
		Kernel* kernel_; 		/**< sph kernel function specific to a SPHBody */
		Real particle_spacing_;						/**< Particle spacing of the body. */
		size_t number_of_particles_;				/**< Number of real particles of the body. */
		BaseParticles* base_particles_;				/**< Base particles of this body. */
		BaseMeshCellLinkedList* base_mesh_cell_linked_list_; /**< Cell linked mesh of this body. */
		MeshBackground* mesh_background_;			/**< Background mesh.*/
		ParticlesGeneratorOps particle_generator_op_;	/**< Particle generator manner */
		PositionsAndVolumes body_input_points_volumes_; /**< For direct generate particles. /

		/**
		 * @brief particle by cells lists is for parallel splitting algorithm.
		 * All particles in each cell are collected together.
		 * If two partiles each belongs two different cell entries,
		 * they have no interaction because they are too far.
		 */
		SplitCellLists split_cell_lists_;

		/** inner configuration for the neighbor relations. */
		ParticleConfiguration inner_configuration_;

		/**
		 * @brief Contact configurations
		 * @details Note that contact configuration only gives all topological relation to this body.
		 * The specific physical interaction, which may not involving all contact bodies,
		 * will be defined in the specific particle dynamics
		 */
		 /** Contact map: pointing to toplogically contacted bodies. **/
		SPHBodyContactMap contact_map_;
		/** Lists of particles has a ocnfiguration with particles in contaced bodies. **/
		ContactParticles indexes_contact_particles_;
		/** Configurations for updated Lagrangian formulation. **/
		ContatcParticleConfiguration contact_configuration_;

		/**
		 * @brief Constructor of SPHBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] smoothinglength_ratio The ratio between smoothinglength to particle spacing.
		 * @param[in] op Partciel generator manner.
		 */
		explicit SPHBody(SPHSystem &sph_system, string body_name, 
			int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op);
		virtual ~SPHBody() {};

		/** Get the name of this body for out file name. */
		string GetBodyName();
		/** Get the name of this body for out file name. */
		Region& getBodyReagion() { return body_region_; };
		/** Set up the contact map. */
		void SetContactMap(SPHBodyContactMap& contact_map);
		/** Set up the contact map. */
		SPHBodyContactMap& getContactMap() { return contact_map_; };
		/** Allocate memory for cell linked list. */
		virtual void AllocateMeoemryCellLinkedList() {};
		/** add the back ground mesh particle mesh interaction. */
		virtual void addBackgroundMesh(Real mesh_size_ratio = 0.5);
		/** Allocate memories for configuration. */
		void AllocateMemoriesForConfiguration();
		/** Allocate extra configuration memories for body buffer particles. */
		void AllocateConfigurationMemoriesForBodyBuffer(size_t body_buffer_particles);
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() = 0;
		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() = 0;
		/** Build contact configuration. */
		virtual void BuildContactConfiguration();
		/** Update inner configuration. */
		virtual void UpdateInnerConfiguration() = 0;
		/** Update contact configuration. */
		virtual void UpdateContactConfiguration() = 0;
		/** Update interactiong configuration. */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) = 0;

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
		/** Constructor of RealBody. */
		RealBody(SPHSystem &sph_system, string body_name, 
			int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op);
		virtual ~RealBody() {};

		/** Allocate memory for cell linked list. */
		virtual void AllocateMeoemryCellLinkedList() override;
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
		/** Constructor of FictitiousBodyBody. */
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
		Region body_part_region_;

	public:
		BodyPart(SPHBody *body, string body_part_name)
			: body_(body), body_part_name_(body_part_name),
			body_part_region_(body_part_name) {};
		virtual ~BodyPart() {};

		Region* GetRegion() { return &body_part_region_; };
	};

	/**
	 * @class BodyPartByParticle
	 * @brief An auxillariy class for SPHBody to 
	 * indicate a part of the body moving together with particles.
	 */
	class BodyPartByParticle : public BodyPart
	{
	protected:
		void tagAParticle(size_t particle_index);
		virtual void TagBodyPartParticles();
	public:
		/** Collection particle in this body part. */
		IndexVector body_part_particles_;

		BodyPartByParticle(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~BodyPartByParticle() {};
	};

	/**
	 * @class BodySurface
	 * @brief A auxillariy class for Body to
	 * indicate the surface particles from background mesh
	 */
	class BodySurface : public BodyPartByParticle
	{
	protected:
		virtual void TagBodyPartParticles() override;
	public:
		BodySurface(SPHBody* body);
		virtual~BodySurface() {};
	};

	/**
	 * @class BodyPartByCell
	 * @brief An auxillariy class for SPHBody to
	 * indicate a part of the body fixed in space.
	 */
	class BodyPartByCell : public BodyPart
	{
	protected:
		virtual void TagBodyPartCells();
	public:
		/** Collection of cells to indicate the body part. */
		CellLists body_part_cells_;

		BodyPartByCell(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~BodyPartByCell() {};
	};

	/**
	 * @class NearBodySurface
	 * @brief An auxillariy class for SPHBody to
	 * indicate region close the body surface.
	 */
	class NearBodySurface : public BodyPartByCell
	{
	protected:
		virtual void TagBodyPartCells();
	public:
		NearBodySurface(SPHBody* body);
		virtual ~NearBodySurface() {};
	};

	/**
	 * @class SolidBodyPartForSimbody
	 * @brief A SolidBodyPart for coupling with Simbody.
	 * The mass, origin, and unit inertial matrix are computed.
	 * Note: In Simbody, all spatial vectors are three dimensional.
	 */
	class SolidBodyPartForSimbody : public BodyPartByParticle
	{
	protected:
		Real solid_body_density_;

		virtual void TagBodyPartParticles() override;
	public:
		SolidBodyPartForSimbody(SPHBody* body, string soild_body_part_name);
		virtual~SolidBodyPartForSimbody() {};

		Vec3 initial_mass_center_;
		SimTK::MassProperties* body_part_mass_properties_;
	};

}
