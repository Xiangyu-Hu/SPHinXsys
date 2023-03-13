#include "fluid_dynamics_multi_phase.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		ViscousAccelerationMultiPhase::ViscousAccelerationMultiPhase(BaseInnerRelation &inner_relation,
																	 BaseContactRelation &contact_relation)
			: ViscousAccelerationInner(inner_relation), MultiPhaseContactData(contact_relation)
		{
			if (&inner_relation.getSPHBody() != &contact_relation.getSPHBody())
			{
				std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_fluids_.push_back(&contact_particles_[k]->fluid_);
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
			}
		}
		//=================================================================================================//
		ViscousAccelerationMultiPhase::
			ViscousAccelerationMultiPhase(ComplexRelation &complex_relation)
			: ViscousAccelerationMultiPhase(complex_relation.getInnerRelation(),
											complex_relation.getContactRelation()) {}
		//=================================================================================================//
		MultiPhaseColorFunctionGradient::
			MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation)
			: LocalDynamics(contact_relation.getSPHBody()), MultiPhaseData(contact_relation),
			  rho0_(sph_body_.base_material_->ReferenceDensity()), Vol_(particles_->Vol_),
			  pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
			  surface_indicator_(particles_->surface_indicator_)
		{
			particles_->registerVariable(color_grad_, "ColorGradient");
			particles_->registerVariable(surface_norm_, "SurfaceNormal");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
				contact_rho0_.push_back(rho0_k);
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//