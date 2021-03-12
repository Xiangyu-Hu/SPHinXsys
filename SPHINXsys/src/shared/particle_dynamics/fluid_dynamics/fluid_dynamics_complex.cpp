/**
 * @file 	fluid_dynamics_complex.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_complex.h"
#include "in_output.h"
#include "geometry_level_set.h"
//=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(ComplexBodyRelation* body_complex_relation, Real thereshold) :
			FreeSurfaceIndicationInner(body_complex_relation->inner_relation_, thereshold),
			FluidContactData(body_complex_relation->contact_relation_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k) {
				Real rho_0_k = contact_particles_[k]->rho_0_;
				contact_inv_rho_0_.push_back(1.0 / rho_0_k);
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		void FreeSurfaceIndicationComplex::Interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIndicationInner::Interaction(index_i, dt);

			Real pos_div = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& contact_mass_k = *(contact_mass_[k]);
				Real contact_inv_rho_0_k = contact_inv_rho_0_[k];
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					pos_div -= contact_neighborhood.dW_ij_[n] * contact_neighborhood.r_ij_[n]
						* contact_inv_rho_0_k * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
		ViscousAccelerationWithWall::ViscousAccelerationWithWall(ComplexBodyRelation* body_complex_relation) :
			ViscousAccelerationInner(body_complex_relation->inner_relation_),
			FluidWallData(body_complex_relation->contact_relation_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_ave_.push_back(&(contact_particles_[k]->vel_ave_));
			}
		}		
		//=================================================================================================//
		void ViscousAccelerationWithWall::Interaction(size_t index_i, Real dt)
		{
			ViscousAccelerationInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0*(vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * smoothing_length_);
					acceleration += 2.0 * mu_ * vel_derivative 
								  * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		AngularConservativeViscousAccelerationWithWall::
			AngularConservativeViscousAccelerationWithWall(ComplexBodyRelation* body_complex_relation) :
			AngularConservativeViscousAccelerationInner(body_complex_relation->inner_relation_),
			FluidWallData(body_complex_relation->contact_relation_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_ave_.push_back(&(contact_particles_[k]->vel_ave_));
			}
		}
		//=================================================================================================//
		void AngularConservativeViscousAccelerationWithWall::Interaction(size_t index_i, Real dt)
		{
			AngularConservativeViscousAccelerationInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			Vecd acceleration(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
					 * is formulation is more accurate thant the previous one for Taylor-Green-Vortex flow. */
					Real v_r_ij = 2.0 * dot(vel_i - vel_ave_k[index_j], r_ij * e_ij);
					Real vel_difference = 0.0 * (vel_i - vel_ave_k[index_j]).norm() * r_ij;
					Real eta_ij = 8.0 * SMAX(mu_, rho_i * vel_difference) * v_r_ij /
						(r_ij * r_ij + 0.01 * smoothing_length_);
					acceleration += eta_ij * Vol_k[index_j] * contact_neighborhood.dW_ij_[n] * e_ij / rho_i;
				}
			}

			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(BaseInnerBodyRelation* inner_relation,
			BaseContactBodyRelation* conatct_relation) :
			ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(
				inner_relation, conatct_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation* body_complex_relation) :
			TransportVelocityCorrectionComplex(body_complex_relation->inner_relation_,
				body_complex_relation->contact_relation_) {}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation* complex_relation,
			BaseContactBodyRelation* extra_conatct_relation) :
			ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(
				complex_relation, extra_conatct_relation)
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

			Real rho_i = rho_n_[index_i];

			Vecd acceleration_trans(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * Vol_k[index_j] * nablaW_ij / rho_i;
				}
			}

			/** correcting particle position */
			if (!is_free_surface_[index_i]) pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationRiemannWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxation<PressureRelaxationRiemannInnerOldroyd_B>::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(wall_Vol_[k]);
				Neighborhood& wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ij * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		void DensityRelaxationRiemannWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxation<DensityRelaxationRiemannInnerOldroyd_B>::Interaction(index_i, dt);

			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(wall_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood& wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];

					Matd velocity_gradient = -SimTK::outer((vel_i - vel_ave_k[index_j]), nablaW_ij) * Vol_k[index_j] * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient
						- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}

			dtau_dt_[index_i] += stress_rate;
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//