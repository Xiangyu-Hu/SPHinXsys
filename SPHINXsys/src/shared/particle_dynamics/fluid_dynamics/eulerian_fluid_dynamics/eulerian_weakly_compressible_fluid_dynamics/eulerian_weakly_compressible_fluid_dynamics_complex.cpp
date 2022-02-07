/**
 * @file 	eulerian_weakly_compressible_fluid_dynamics_complex.cpp
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
 */

#include "eulerian_weakly_compressible_fluid_dynamics_complex.h"
#include "eulerian_weakly_compressible_fluid_dynamics_complex.hpp"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.hpp"
//=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation, Real thereshold)
			: FreeSurfaceIndicationInner(inner_relation, thereshold), WCFluidContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				Real rho0_k = contact_particles_[k]->rho0_;
				contact_inv_rho0_.push_back(1.0 / rho0_k);
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(ComplexBodyRelation &complex_relation, Real thereshold)
			: FreeSurfaceIndicationComplex(complex_relation.inner_relation_,
				complex_relation.contact_relation_, thereshold) {}
		//=================================================================================================//
		void FreeSurfaceIndicationComplex::Interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIndicationInner::Interaction(index_i, dt);

			Real pos_div = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
				Real contact_inv_rho0_k = contact_inv_rho0_[k];
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					pos_div -= contact_neighborhood.dW_ij_[n] * contact_neighborhood.r_ij_[n] *
						contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//