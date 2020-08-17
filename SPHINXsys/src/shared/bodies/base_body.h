/**
 * @file 	base_body.h
 * @brief 	This is the base classes of SPH bodies. The real body is for 
 *			that with cell linked list and the fictitious one does not.     
 * 			Before the definition of the SPH bodies, the shapes with complex 
 *			geometries, i.e. those are produced by advanced binary operation, 
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
	class SPHBodyBaseRelation;
	class BaseLevelSet;

	/**
	 * @class ParticlesGeneratorOps
	 * @brief Serval manners are provied for particles generator.
	 * @details lattice : Generate partice from lattice grid.
	 *			direct  : Input particle position and volume directly.
	 *			regularized : geometry will be regularized with level set technique.
	 */
	enum class ParticlesGeneratorOps {lattice, direct, regularized};
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
		string body_name_; 		/**< name of this body */
		bool newly_updated_;		/**< whether this body is in a newly updated state */
		/** describe the geometry of the body, static member, so the geometry head file is included. */
		ComplexShape  body_shape_;
		/** smoothing length. */
		Real smoothing_length_;
		/** Computational domain bounds of the body for boundary conditions. */
		Vecd body_lower_bound_, body_upper_bound_;
		/** Whether the computational domain bound for this body is prescribed. */
		bool prescribed_body_bounds_;
		/** Computing particle spacing from refinement level. */
		Real RefinementLevelToParticleSpacing();

		/** Generate a kernel. */
		Kernel* GenerateAKernel(Real smoothing_length);
		/** Change kernel function specific for this body. */
		void ReplaceKernelFunction(Kernel* kernel);
	public:
		int refinement_level_;	/**< refinement level of this body */
		Kernel* kernel_; 		/**< sph kernel function specific to a SPHBody */
		Real particle_spacing_;						/**< Particle spacing of the body. */
		size_t number_of_particles_;				/**< Number of real particles of the body. */
		BaseParticles* base_particles_;				/**< Base particles of this body. */
		BaseMeshCellLinkedList* base_mesh_cell_linked_list_; /**< Cell linked mesh of this body. */
		ParticlesGeneratorOps particle_generator_op_;	/**< Particle generator manner */
		PositionsAndVolumes body_input_points_volumes_; /**< For direct generate particles. */
		BaseLevelSet* levelset_mesh_;					/**< narrow bounded levelset mesh. */
		/**
		 * @brief particle by cells lists is for parallel splitting algorithm.
		 * All particles in each cell are collected together.
		 * If two partiles each belongs two different cell entries,
		 * they have no interaction because they are too far.
		 */
		SplitCellLists split_cell_lists_;

		/** all contact relations centered from this body **/
		StdVec<SPHBodyBaseRelation*> body_relations_;

		/**
		 * @brief Constructor of SPHBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] smoothing_length_ratio The ratio between smoothinglength to particle spacing.
		 * @param[in] op Particle generator manner.
		 */
		explicit SPHBody(SPHSystem &sph_system, string body_name, 
			int refinement_level, Real smoothing_length_ratio, ParticlesGeneratorOps op);
		virtual ~SPHBody() {};

		/** Get the name of this body for out file name. */
		string GetBodyName();
		void setNewlyUpdated() { newly_updated_ = true; };
		bool checkNewlyUpdated() { return newly_updated_; };
		void setNotNewlyUpdated() { newly_updated_ = false; };

		/** Get the name of this body for out file name. */
		ComplexShape& getBodyShape() { return body_shape_; };
		void setBodyLowerBound(Vecd lower_bound) { body_lower_bound_ = lower_bound; };
		void setBodyUpperBound(Vecd upper_bound) { body_upper_bound_ = upper_bound; };
		Vecd getBodyLowerBound() { return body_lower_bound_; };
		Vecd getBodyUpperBound() { return body_upper_bound_; };
		void getSPHSystemBound(Vecd& system_lower_bound, Vecd& system_uppwer_bound);

		/** add levelset mesh */
		virtual void addLevelsetMesh(Real mesh_size_ratio = 4.0);
		/** Allocate memory for cell linked list. */
		virtual void AllocateMemoryCellLinkedList() = 0;
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() = 0;
		/** Allocate extra configuration memories for body buffer particles. */
		void AllocateConfigurationMemoriesForBodyBuffer();

		/** Check wether a point within the geometry of this body.
		 * @returns TRUE if a point within body's region otherwise FALSE. 
		 */
		bool checkBodyShapeContain(Vecd pnt); 
		/** Find closest point from a given point to body surface. */
		Vecd ClosestPointOnBodySurface(Vecd input_pnt);
		/**
		 * @brief Find the lower and upper bounds of the body.
		 * @param[in,out] lower_bound Lower bound of this body.
		 * @param[in,out] upper_bound Upper bound of this body.
		 */
		void findBodyShapeBounds(Vecd &lower_bound, Vecd &upper_bound);

		/** Output particle data in VTU file for visualization in Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file);
		/** Output particle data in PLT file for visualization in Tecplot. */
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
			int refinement_level, Real smoothing_length_ratio, ParticlesGeneratorOps op);
		virtual ~RealBody() {};

		/** Allocate memory for cell linked list. */
		virtual void AllocateMemoryCellLinkedList() override;
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() override;
		/** The pointer to derived class object. */
		virtual RealBody* PointToThisObject() override;
	};

	/**
	 * @class FictitiousBody.
	 * @brief Derived class from SPHBody. 
	 * Without inner configuration or inner interaction.
	 */
	class FictitiousBody : public SPHBody
	{
	protected:

	public:
		/** Constructor of FictitiousBodyBody. */
		FictitiousBody(SPHSystem &system, string body_name, 
			int refinement_level, Real smoothing_length_ratio, ParticlesGeneratorOps op);
		virtual ~FictitiousBody() {};

		/** Allocate memory for cell linked list. */
		virtual void AllocateMemoryCellLinkedList() override;
		/** Update cell linked list. */
		virtual void UpdateCellLinkedList() override;
		/** The pointer to derived class object. */
		virtual FictitiousBody* PointToThisObject() override;
	};

	/**
	 * @class BodyPart
	 * @brief An auxillary class for SPHBody to indicate a part of the body.
	 */
	class BodyPart
	{
	public:
		BodyPart(SPHBody *body, string body_part_name)
			: body_(body), body_part_name_(body_part_name),
			body_part_shape_(body_part_name) {};
		virtual ~BodyPart() {};

		ComplexShape& getBodyPartShape() { return body_part_shape_; };
		SPHBody* getBody() { return body_; };

	protected:
		SPHBody* body_;
		string body_part_name_;
		ComplexShape body_part_shape_;

		virtual void TagBodyPart() = 0;
	};

	/**
	 * @class BodyPartByParticle
	 * @brief An auxillary class for SPHBody to 
	 * indicate a part of the body moving together with particles.
	 */
	class BodyPartByParticle : public BodyPart
	{
	public:
		/** Collection particle in this body part. */
		IndexVector body_part_particles_;

		BodyPartByParticle(SPHBody* body, string body_part_name)
			: BodyPart(body, body_part_name) {};
	virtual ~BodyPartByParticle() {};

	protected:
		void tagAParticle(size_t particle_index);
		virtual void TagBodyPart() override;

	};

	/**
	 * @class BodySurface
	 * @brief A auxillary class for Body to
	 * indicate the surface particles from background mesh
	 */
	class BodySurface : public BodyPartByParticle
	{
	public:
		BodySurface(SPHBody* body);
		virtual~BodySurface() {};

	protected:
		virtual void TagBodyPart() override;
	};

	/**
	 * @class BodySurfaceLayer
	 * @brief A auxillary class for Body to
	 * indicate the particles within the inner layers
	 */
	class BodySurfaceLayer : public BodyPartByParticle
	{
	public:
		BodySurfaceLayer(SPHBody* body, Real layer_thickness = 3.0);
		virtual~BodySurfaceLayer() {};

	protected:
		Real layer_thickness_;

		virtual void TagBodyPart() override;
	};

	/**
	 * @class BodyPartByCell
	 * @brief An auxillary class for SPHBody to
	 * indicate a part of the body fixed in space.
	 */
	class BodyPartByCell : public BodyPart
	{
	public:
		/** Collection of cells to indicate the body part. */
		CellLists body_part_cells_;

		BodyPartByCell(SPHBody *body, string body_part_name)
			: BodyPart(body, body_part_name) {};
		virtual ~BodyPartByCell() {};

	protected:
		virtual void TagBodyPart() override;
	};

	/**
	 * @class NearBodySurface
	 * @brief An auxillary class for SPHBody to
	 * indicate region close the body surface.
	 */
	class NearBodySurface : public BodyPartByCell
	{
	public:
		NearBodySurface(SPHBody* body);
		virtual ~NearBodySurface() {};

	protected:
		virtual void TagBodyPart() override;
	};

	/**
	 * @class SolidBodyPartForSimbody
	 * @brief A SolidBodyPart for coupling with Simbody.
	 * The mass, origin, and unit inertial matrix are computed.
	 * Note: In Simbody, all spatial vectors are three dimensional.
	 */
	class SolidBodyPartForSimbody : public BodyPartByParticle
	{
	public:
		Vec3d initial_mass_center_;
		SimTK::MassProperties* body_part_mass_properties_;
		
		SolidBodyPartForSimbody(SPHBody* body, string solid_body_part_name);
		virtual~SolidBodyPartForSimbody() {};
	protected:
		Real solid_body_density_;

		virtual void TagBodyPart() override;
	};
}
