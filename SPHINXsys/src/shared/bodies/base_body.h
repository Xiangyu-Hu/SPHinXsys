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
 */

#ifndef BASE_BODY_H
#define BASE_BODY_H

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "adaptation.h"
#include "cell_linked_list.h"
#include "particle_sorting.h"
#include "all_geometries.h"
#include "generative_structures.h"

#include <string>

namespace SPH
{
	class SPHSystem;
	class BaseParticles;
	class SPHBodyRelation;
	class ComplexShape;
	class BodySurface;

	/**
	 * @class SPHBody
	 * @brief SPHBody is a base body with basic data and functions.
	 *		  Its derived class can be a real fluid body, a real deformable solid body,
	 *        a static or moving solid body or a fictitious body.
	 * 		  Note that only real bodies have cell linked list.
	 */
	class SPHBody
	{
	private:
		SharedPtrKeeper<SPHAdaptation> sph_adaptation_ptr_keeper_;
		UniquePtrKeeper<GenerativeStructure> generative_structure_ptr_keeper_;

	protected:
		SPHSystem &sph_system_;
		std::string body_name_;
		bool newly_updated_; /**< whether this body is in a newly updated state */
		/**< Computational domain bounds for boundary conditions. 
		 * Note that domain bounds may be different from those of the initial body geometry. */
		BoundingBox body_domain_bounds_;
		bool is_domain_bounds_determined_;

	public:
		ComplexShape body_shape_;					/**< describe the volumetric geometry of the body */
		SPHAdaptation *sph_adaptation_;				/**< Particle adapation policy. */
		BaseParticles *base_particles_;				/**< Base particles of this body. */
		GenerativeStructure *generative_structure_; /**< structure which can be used to generate particles or/and configurations directly*/
		/**
		 * @brief particle by cells lists is for parallel splitting algorithm.
		 * All particles in each cell are collected together.
		 * If two partiles each belongs two different cell entries,
		 * they have no interaction because they are too far.
		 */
		SplitCellLists split_cell_lists_;

		StdVec<SPHBodyRelation *> body_relations_; /**< all contact relations centered from this body **/

		explicit SPHBody(SPHSystem &sph_system, const std::string &body_name,
						 SharedPtr<SPHAdaptation> sph_adaptation_ptr);
		virtual ~SPHBody(){};

		std::string getBodyName();
		SPHSystem &getSPHSystem();
		Real getSPHBodyResolutionRef() { return sph_adaptation_->ReferenceSpacing(); };
		void setNewlyUpdated() { newly_updated_ = true; };
		void setNotNewlyUpdated() { newly_updated_ = false; };
		bool checkNewlyUpdated() { return newly_updated_; };
		void setBodyDomainBounds(BoundingBox body_domain_bounds);
		BoundingBox getBodyDomainBounds();
		BoundingBox getSPHSystemBounds();
		/** create the generative structure from outside */
		template <class GenerativeStructureType, typename... ConstructorArgs>
		GenerativeStructureType *createGenerativeStructure(ConstructorArgs &&...args)
		{
			generative_structure_ = generative_structure_ptr_keeper_.createPtr<GenerativeStructureType>(std::forward<ConstructorArgs>(args)...);
			return DynamicCast<GenerativeStructureType>(this, generative_structure_);
		};

		/** This will be called in BaseParticle constructor
		 * and is important because particles are not defined in SPHBody constructor.  */
		virtual void assignBaseParticles(BaseParticles *base_particles);
		void allocateConfigurationMemoriesForBufferParticles();

		virtual void writeParticlesToVtuFile(std::ostream &output_file);	
		virtual void writeParticlesToVtpFile(std::ostream &output_file);;
		virtual void writeParticlesToPltFile(std::ofstream &output_file);
		virtual void writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface& surface_particles);
		virtual void writeParticlesToXmlForRestart(std::string &filefullpath);
		virtual void readParticlesFromXmlForRestart(std::string &filefullpath);
		virtual void writeToXmlForReloadParticle(std::string &filefullpath);
		virtual void readFromXmlForReloadParticle(std::string &filefullpath);
		virtual SPHBody *ThisObjectPtr() { return this; };
	};

	/**
	 * @class RealBody
	 * @brief Derived class from SPHBody. 
	 * With inner particle configuration or inner interactions.
	 */
	class RealBody : public SPHBody
	{
	private:
		UniquePtrKeeper<BaseCellLinkedList> cell_linked_list_keeper_;

	public:
		ParticleSorting particle_sorting_;
		BaseCellLinkedList *cell_linked_list_; /**< Cell linked mesh of this body. */

		RealBody(SPHSystem &sph_system, const std::string &body_name,
				 SharedPtr<SPHAdaptation> sph_adaptation_ptr);
		virtual ~RealBody(){};

		/** This will be called in BaseParticle constructor
		 * and is important because particles are not defined in FluidBody constructor.  */
		virtual void assignBaseParticles(BaseParticles *base_particles) override;
		virtual void sortParticleWithCellLinkedList();
		virtual void updateCellLinkedList();
	};

	/**
	 * @class FictitiousBody
	 * @brief Derived class from SPHBody. 
	 * Without inner configuration or inner interaction.
	 */
	class FictitiousBody : public SPHBody
	{
	public:
		FictitiousBody(SPHSystem &system, const std::string &body_name,
					   SharedPtr<SPHAdaptation> sph_adaptation_ptr);
		virtual ~FictitiousBody(){};
	};

	/**
	 * @class BodyPart
	 * @brief An auxillary class for SPHBody to indicate a part of the body.
	 */
	using namespace std::placeholders;
	class BodyPart
	{
	public:
		BodyPart(SPHBody &sph_body, const std::string &body_part_name)
			: sph_body_(&sph_body), body_part_name_(body_part_name), body_part_bounds_(Vecd(0), Vecd(0)), body_part_bounds_set_(false)
			{};
		virtual ~BodyPart(){};

		SPHBody *getSPHBody() { return sph_body_; };
		std::string BodyPartName() { return body_part_name_; };

		void setBodyPartBounds(BoundingBox bbox){
			body_part_bounds_ = bbox;
			body_part_bounds_set_ = true;
		};

		BoundingBox getBodyPartBounds(){
			if (!body_part_bounds_set_) std::cout << "WARNING: the body part bounds are not set for BodyPart." << std::endl;
			return body_part_bounds_;
		}

	protected:
		SPHBody *sph_body_;
		std::string body_part_name_;

		BoundingBox body_part_bounds_;
		bool body_part_bounds_set_;
	};

	/**
	 * @class BodyPartByParticle
	 * @brief A body part with a collection of particles.
	 */
	class BodyPartByParticle : public BodyPart
	{
	public:
		IndexVector body_part_particles_; /**< Collection particle in this body part. */

		BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name)
			: BodyPart(sph_body, body_part_name), base_particles_(sph_body.base_particles_)
			  {};
		virtual ~BodyPartByParticle(){};

	protected:
		BaseParticles *base_particles_;

		typedef std::function<void(size_t)> TaggingParticleMethod;
		void tagParticles(TaggingParticleMethod &tagging_particle_method);
	};

	/**
	 * @class BodyPartByCell
	 * @brief A body part with a collection of cell lists.
	 */
	class BodyPartByCell : public BodyPart
	{
	public:
		CellLists body_part_cells_; /**< Collection of cells to indicate the body part. */

		BodyPartByCell(RealBody &real_body, const std::string &body_part_name)
			: BodyPart(real_body, body_part_name), cell_linked_list_(real_body.cell_linked_list_){};
		virtual ~BodyPartByCell(){};

	protected:
		BaseCellLinkedList *cell_linked_list_;
		typedef std::function<bool(Vecd, Real)> TaggingCellMethod;
		void tagCells(TaggingCellMethod &tagging_cell_method);
	};

	/**
	 * @class BodyRegionByParticle
	 * @brief A  body part with the collection of particles within by a presribed shape.
	 */
	class BodyRegionByParticle : public BodyPartByParticle
	{
	public:
		Shape &body_part_shape_;

		BodyRegionByParticle(SPHBody &sph_body, const std::string &body_part_name, Shape &shape);
		virtual ~BodyRegionByParticle(){};

	private:
		void tagByContain(size_t particle_index);
	};

	/**
	 * @class BodySurface
	 * @brief A  body part with the collection of particles at surface of a body
	 */
	class BodySurface : public BodyPartByParticle
	{
	public:
		explicit BodySurface(SPHBody &sph_body);
		virtual ~BodySurface(){};

	private:
		Real particle_spacing_min_;
		void tagNearSurface(size_t particle_index);
	};

	/**
	 * @class BodySurfaceLayer
	 * @brief A  body part with the collection of particles within the surface layers of a body.
	 */
	class BodySurfaceLayer : public BodyPartByParticle
	{
	public:
		explicit BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness = 3.0);
		virtual ~BodySurfaceLayer(){};

	private:
		Real thickness_threshold_;
		void tagSurfaceLayer(size_t particle_index);
	};

	/**
	 * @class BodyRegionByCell
	 * @brief A body part with the cell lists within a prescribed shape.
	 */
	class BodyRegionByCell : public BodyPartByCell
	{
	public:
		Shape &body_part_shape_;

		BodyRegionByCell(RealBody &real_body, const std::string &body_part_name, Shape &shape);
		virtual ~BodyRegionByCell(){};

	private:
		bool checkNotFar(Vecd cell_position, Real threshold);
	};

	/**
	 * @class NearShapeSurface
	 * @brief A body part with the cell lists near the surface of a prescribed shape.
	 */
	class NearShapeSurface : public BodyPartByCell
	{
	private:
		UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;

	public:
		LevelSetShape &level_set_shape_;

		/** for the case that the body part shape is not that of the body */
		NearShapeSurface(RealBody &real_body, const std::string &body_part_name, Shape &shape);
		/** for the case that the body part is the surface of the body shape */
		explicit NearShapeSurface(RealBody &real_body);
		virtual ~NearShapeSurface(){};

	private:
		/** only cells near the surface of the body part shape are included */
		bool checkNearSurface(Vecd cell_position, Real threshold);
	};

	/**
	 * @class TreeTerminates
	 * @brief A  body part with the collection of particles as the terminates of the tree. 
	 */
	class TreeTerminates : public BodyPartByParticle
	{
	public:
		GenerativeTree &tree_;

		explicit TreeTerminates(SPHBody &sph_body);
		virtual ~TreeTerminates(){};
	};
}
#endif //BASE_BODY_H