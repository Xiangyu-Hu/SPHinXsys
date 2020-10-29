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
#include "all_particle_generators.h"
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
		/** Computational domain bounds of the body for boundary conditions. */
		Vecd body_lower_bound_, body_upper_bound_;
		/** Whether the computational domain bound for this body is prescribed. */
		bool prescribed_body_bounds_;
		/** smoothing length. */
		Real smoothing_length_;
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
		BaseMeshCellLinkedList* mesh_cell_linked_list_; /**< Cell linked mesh of this body. */
		ParticleGenerator* particle_generator_;	/**< Particle generator manner */
		PositionsAndVolumes body_input_points_volumes_; /**< For direct generate particles. */
		ComplexShape*  body_shape_;		/** describe the geometry of the body*/
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
		 * @param[in] particle_generator Particle generator.
		 */
		explicit SPHBody(SPHSystem &sph_system, string body_name, int refinement_level, Real smoothing_length_ratio, 
			ParticleGenerator* particle_generator = new ParticleGeneratorLattice());
		virtual ~SPHBody() {};

		/** Get the name of this body for out file name. */
		string GetBodyName();
		void setNewlyUpdated() { newly_updated_ = true; };
		bool checkNewlyUpdated() { return newly_updated_; };
		void setNotNewlyUpdated() { newly_updated_ = false; };
		SPHSystem& getSPHSystem();

		/** Get the name of this body for out file name. */
		void setBodyLowerBound(Vecd lower_bound) { body_lower_bound_ = lower_bound; };
		void setBodyUpperBound(Vecd upper_bound) { body_upper_bound_ = upper_bound; };
		Vecd getBodyLowerBound() { return body_lower_bound_; };
		Vecd getBodyUpperBound() { return body_upper_bound_; };
		void getSPHSystemBound(Vecd& system_lower_bound, Vecd& system_uppwer_bound);

		/** assign base particle to the body and cell linked list. */
		void assignBaseParticle(BaseParticles* base_particles);
		/** Compute reference number density*/
		virtual Real computeReferenceNumberDensity();
		/** Update cell linked list. */
		virtual void updateCellLinkedList() = 0;
		/** Allocate extra configuration memories for body buffer particles. */
		void allocateConfigurationMemoriesForBodyBuffer();

		/**
		 * @brief Find the lower and upper bounds of the body.
		 * @param[in,out] lower_bound Lower bound of this body.
		 * @param[in,out] upper_bound Upper bound of this body.
		 */
		void findBodyDomainBounds(Vecd &lower_bound, Vecd &upper_bound);

		/** Output particle data in VTU file for visualization in Paraview. */
		virtual void writeParticlesToVtuFile(ofstream &output_file);
		/** Output particle data in PLT file for visualization in Tecplot. */
		virtual void writeParticlesToPltFile(ofstream &output_file);

		/** Output particle data in XML file for restart simulation. */
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath);
		/** Read particle data in XML file for restart simulation. */
		virtual void readParticlesFromXmlForRestart(std::string &filefullpath);

		/** Output particle position and volume in XML file for reloading particles. */
		virtual void writeToXmlForReloadParticle(std::string &filefullpath);
		/** Reload particle position and volume from XML files. */
		virtual void readFromXmlForReloadParticle(std::string &filefullpath);
		
		/** The pointer to derived class object. */
		virtual SPHBody* pointToThisObject();
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
		RealBody(SPHSystem &sph_system, string body_name, int refinement_level, Real smoothing_length_ratio, 
			ParticleGenerator* particle_generator = new ParticleGeneratorLattice());
		virtual ~RealBody() {};

		/** Update cell linked list. */
		virtual void updateCellLinkedList() override;
		/** The pointer to derived class object. */
		virtual RealBody* pointToThisObject() override;
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
		FictitiousBody(SPHSystem &system, string body_name, int refinement_level, Real smoothing_length_ratio,
			ParticleGenerator* particle_generator = new ParticleGeneratorDirect());
		virtual ~FictitiousBody() {};

		/** Update cell linked list. */
		virtual void updateCellLinkedList() override;
		/** The pointer to derived class object. */
		virtual FictitiousBody* pointToThisObject() override;
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
			body_part_shape_(NULL) {};
		virtual ~BodyPart() {};

		ComplexShape* getBodyPartShape() { return body_part_shape_; };
		SPHBody* getBody() { return body_; };
		string BodyPartName() { return body_part_name_; };
		/**
		 * @brief Find the lower and upper bounds of the body part.
		 * @param[in,out] lower_bound Lower bound of this body part.
		 * @param[in,out] upper_bound Upper bound of this body part.
		 */
		void BodyPartBounds(Vecd &lower_bound, Vecd &upper_bound)
		{
			body_part_shape_->findBounds(lower_bound, upper_bound);
		};
	protected:
		SPHBody* body_;
		string body_part_name_;
		ComplexShape* body_part_shape_;

		virtual void tagBodyPart() = 0;
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
		virtual void tagBodyPart() override;

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
		virtual void tagBodyPart() override;
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

		virtual void tagBodyPart() override;
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
		virtual void tagBodyPart() override;
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
		virtual void tagBodyPart() override;
	};
}
