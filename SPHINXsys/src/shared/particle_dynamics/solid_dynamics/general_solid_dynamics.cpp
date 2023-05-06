#include "general_solid_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		CorrectConfiguration::
			CorrectConfiguration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), SolidDataInner(inner_relation),
			  B_(particles_->B_) {}
		//=================================================================================================//
	}
}
