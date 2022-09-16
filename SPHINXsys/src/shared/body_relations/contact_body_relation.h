/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	contact_body_relation.h
 * @brief 	The topological relations between bodies,
 * 			and the corresponding local topologies (particle configurations)
 * 			are constructed in these classes.
 * @author	Xiangyu Hu
 */

#ifndef CONTACT_BODY_RELATION_H
#define CONTACT_BODY_RELATION_H

#include "base_body_relation.h"
#include "inner_body_relation.h"

namespace SPH
{
	/**
	 * @class BodyRelationContact
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class BodyRelationContact : public BaseBodyRelationContact
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;

		void initialization();

	public:
		BodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies);
		BodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts);
		virtual ~BodyRelationContact(){};
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SolidBodyRelationContact
	 * @brief The relation between a solid body and its contact solid bodies
	 * TODO: better called BodySurfaceContact
	 */
	class SolidBodyRelationContact : public BaseBodyRelationContact
	{
	private:
		UniquePtrKeeper<BodySurfaceLayer> shape_surface_ptr_keeper_;

	public:
		BodySurfaceLayer *body_surface_layer_;

		SolidBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies);
		SolidBodyRelationContact(SolidBodyRelationSelfContact &solid_body_relation_self_contact,
								 RealBodyVector contact_bodies);
		virtual ~SolidBodyRelationContact(){};
		BodyPartByParticle &getDynamicsRange() { return *body_surface_layer_; };

		virtual void updateConfiguration() override;

	protected:
		IndexVector &body_part_particles_;
		BodyPartParticlesIndex get_body_part_particle_index_;

		void initialization();
		virtual void resetNeighborhoodCurrentSize() override;
	};

	/**
	 * @class BodyRelationContactToBodyPart
	 * @brief The relation between a SPH body and a vector of body parts.
	 */
	class BodyRelationContactToBodyPart : public BodyRelationContact
	{
	protected:
		UniquePtrKeepers<NeighborBuilderContactBodyPart> neighbor_relation_contact_body_part_ptr_vector_keeper_;

	public:
		BodyPartVector contact_body_parts_;
		StdVec<NeighborBuilderContactBodyPart *> get_part_contact_neighbors_;

		BodyRelationContactToBodyPart(RealBody &real_body, BodyPartVector contact_body_parts);
		virtual ~BodyRelationContactToBodyPart(){};

		virtual void updateConfiguration() override;
	};
}
#endif // CONTACT_BODY_RELATION_H