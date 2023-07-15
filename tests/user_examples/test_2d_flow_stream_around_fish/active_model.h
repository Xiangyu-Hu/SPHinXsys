#ifndef ACTIVE_MODEL_H
#define ACTIVE_MODEL_H

#include "elastic_dynamics.h"
#include "composite_material.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		 * @class ActiveIntegration1stHalf
		 */
		class ActiveIntegration1stHalf : public Integration1stHalfPK2
		{
		public:
			explicit ActiveIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~ActiveIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> &active_strain_,F_0, E_e;
			StdLargeVec<int>& material_id_;
		};
	}
}
#endif // ACTIVE_MODEL_H