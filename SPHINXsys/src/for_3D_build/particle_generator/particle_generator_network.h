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
 * @file 	particle_generator_network.h
 * @brief 	This is a class of particle generator, which generates particles
 * 			with in network or tree form. 
 * @author	Chi ZHang and Xiangyu Hu
 */
#pragma once

#include "sph_data_conainers.h"
#include "base_particle_generator.h"
#include "base_geometry.h"

namespace SPH 
{
	class Tree;
	class BaseLevelSet;
	class BaseMeshCellLinkedList;
	class ComplexShape;

	/**
	 * @class ParticleGeneratorNetwork
	 * @brief Generate a tree-shape network for the conduction system of a heart with particles.
	 */
	class ParticleGeneratorNetwork : public ParticleGenerator
	{
	public:
		ParticleGeneratorNetwork(Vecd starting_pnt, Vecd second_pnt);
		virtual ~ParticleGeneratorNetwork() {};

		virtual void initialize(SPHBody* sph_body) override;
		virtual void createBaseParticles(BaseParticles* base_particles) override;
	protected:
		Vecd starting_pnt_;	/**< Starting point for net work. */
		Vecd second_pnt_;		/**< Second point, approximate the growing direction. */
		size_t n_it_; 			/**< Number of iterations (generations of branch. */
		bool fascicles_;		/**< Create fascicles? */
		size_t	segments_in_branch_;	/**< approximated number of segments in a branch. */
		Real segment_length_;			/**< segment length of the branch. */
		Real angle_ = 0.2;		/**< angle with respect to the direction of the previous edge and the new edge. */
		Real repulsivity_ = 0.1; 			/**< repulsivity parameter. */
		std::vector<Real> fascicle_angles_ = {-1.5, 0.2}; 	/**< angles with respect to the initial edge of the fascicles.*/
		Real fascicle_ratio_ = 5.0; 	/**< ratio of length  of the fascicles. Include one per fascicle to include.*/
		ComplexShape* body_shape_;
		RealBody* real_body_;

		Vecd getGradientFromNearestPoints(Vecd pt, Real delta, BaseMeshCellLinkedList* mesh_cell_linked_list);
		bool createABranchIfValid(SPHBody* sph_body, size_t parent_id, Real angle,
			Real repulsivity, size_t number_segments, Tree* tree);
		/** Creates a new tentative point at the body surface */
		Vecd creatATentativeNewBranchVecd(Vecd init_point, Vecd dir);
		bool isCollision(Vecd& new_point, ListData& nearest_neighbor, size_t parent_id, Tree* tree);
	};
}