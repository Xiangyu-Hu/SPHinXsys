/**
 * @file 	fluid_dynamics_complex.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(BaseBodyRelationInner &inner_relation,
											   BaseBodyRelationContact &contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(inner_relation, contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation)
			: TransportVelocityCorrectionComplex(complex_relation.inner_relation_,
												 complex_relation.contact_relation_) {}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation,
											   BaseBodyRelationContact &extra_contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(complex_relation, extra_contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::prepareContactData()
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::Interaction(size_t index_i, Real dt)
		{
			TransportVelocityCorrectionInner::Interaction(index_i, dt);

			Real rho_i = rho_[index_i];

			Vecd acceleration_trans(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * Vol_k[index_j] * nablaW_ij / rho_i;
				}
			}

			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
				pos_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxation<PressureRelaxationInnerOldroyd_B>::Interaction(index_i, dt);

			Real rho_i = rho_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ij * Vol_k[index_j] / rho_i;
				}
			}

			acc_[index_i] += acceleration;
		}
		//=================================================================================================//
		void DensityRelaxationWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxation<DensityRelaxationInnerOldroyd_B>::Interaction(index_i, dt);

			Vecd vel_i = vel_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];

					Matd velocity_gradient = -SimTK::outer((vel_i - vel_ave_k[index_j]), nablaW_ij) * Vol_k[index_j] * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
								   (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}
			dtau_dt_[index_i] += stress_rate;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//