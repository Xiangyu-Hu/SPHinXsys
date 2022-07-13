/**
 * @file 	complex_body_relation.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "complex_body_relation.h"

#include "inner_body_relation.h"
#include "contact_body_relation.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	ComplexBodyRelation::
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation)
		: SPHBodyRelation(*inner_relation.sph_body_),
		  inner_relation_(inner_relation),
		  contact_relation_(contact_relation),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(RealBody &real_body, RealBodyVector contact_bodies)
		: SPHBodyRelation(real_body),
		  inner_relation_(base_body_relation_inner_ptr_keeper_.createRef<BodyRelationInner>(real_body)),
		  contact_relation_(base_body_relation_contact_ptr_keeper_
								.createRef<BodyRelationContact>(real_body, contact_bodies)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, RealBodyVector contact_bodies)
		: SPHBodyRelation(*inner_relation.sph_body_),
		  inner_relation_(inner_relation),
		  contact_relation_(base_body_relation_contact_ptr_keeper_.createRef<BodyRelationContact>(
			  DynamicCast<RealBody>(this, *sph_body_), contact_bodies)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(RealBody &real_body, BodyPartVector contact_body_parts)
		: SPHBodyRelation(real_body),
		  inner_relation_(base_body_relation_inner_ptr_keeper_.createRef<BodyRelationInner>(real_body)),
		  contact_relation_(base_body_relation_contact_ptr_keeper_
								.createRef<BodyRelationContactToBodyPart>(real_body, contact_body_parts)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfigurationMemories()
	{
		inner_relation_.updateConfigurationMemories();
		contact_relation_.updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfiguration()
	{
		inner_relation_.updateConfiguration();
		contact_relation_.updateConfiguration();
	}
	//=================================================================================================//
}
