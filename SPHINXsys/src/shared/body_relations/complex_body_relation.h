/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and Hu1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	complex_body_relation.h
 * @brief 	The topological relations within one body and to other bodies.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
 */
 
#ifndef COMPLEX_BODY_RELATION_H
#define COMPLEX_BODY_RELATION_H

#include "base_body_relation.h"

namespace SPH
{
	/**
	 * @class ComplexBodyRelation
	 * @brief The relation combined an inner and a contact body relation.
	 * The interaction is in a inner-boundary-condition fashion. Here inner interaction is
	 * different from contact interaction.
	 */
	class ComplexBodyRelation : public SPHBodyRelation
	{
	private:
		UniquePtrKeeper<BaseBodyRelationInner> base_body_relation_inner_ptr_keeper_;
		UniquePtrKeeper<BaseBodyRelationContact> base_body_relation_contact_ptr_keeper_;

	public:
		BaseBodyRelationInner &inner_relation_;
		BaseBodyRelationContact &contact_relation_;
		RealBodyVector contact_bodies_;
		ParticleConfiguration &inner_configuration_;
		ContactParticleConfiguration &contact_configuration_;

		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation);
		ComplexBodyRelation(RealBody &real_body, RealBodyVector contact_bodies);
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, RealBodyVector contact_bodies);
		ComplexBodyRelation(RealBody &real_body, BodyPartVector contact_body_parts);
		virtual ~ComplexBodyRelation(){};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};
}
#endif //COMPLEX_BODY_RELATION_H