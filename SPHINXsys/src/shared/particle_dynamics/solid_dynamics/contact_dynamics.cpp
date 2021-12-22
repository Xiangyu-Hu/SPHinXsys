/**
 * @file 	contact_dynamics.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "contact_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		SelfContactDensitySummation::
			SelfContactDensitySummation(SolidBodyRelationSelfContact* self_contact_relation)
			: PartInteractionDynamicsByParticle(self_contact_relation->sph_body_,
												&self_contact_relation->body_surface_layer_),
			  SolidDataInner(self_contact_relation),
			  mass_(particles_->mass_), contact_density_(particles_->contact_density_) {}
		//=================================================================================================//
		void SelfContactDensitySummation::Interaction(size_t index_i, Real dt)
		{
			Real sigma = 0.0;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				sigma += inner_neighborhood.W_ij_[n] * mass_[inner_neighborhood.j_[n]];
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ContactDensitySummation::
			ContactDensitySummation(SolidBodyRelationContact *solid_body_contact_relation)
			: PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
												&solid_body_contact_relation->body_surface_layer_),
			  ContactDynamicsData(solid_body_contact_relation),
			  mass_(particles_->mass_), contact_density_(particles_->contact_density_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		void ContactDensitySummation::Interaction(size_t index_i, Real dt)
		{
			/** Contact interaction. */
			Real sigma = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					sigma += contact_neighborhood.W_ij_[n] * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ShellContactDensity::
			ShellContactDensity(SolidBodyRelationContact *solid_body_contact_relation)
			: PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
												&solid_body_contact_relation->body_surface_layer_),
			  ContactDynamicsData(solid_body_contact_relation), pos_n_(particles_->pos_n_),
			  contact_density_(particles_->contact_density_),
			  kernel_(solid_body_contact_relation->sph_body_->particle_adaptation_->getKernel()), 
			  spacing_ref_(solid_body_contact_relation->sph_body_->particle_adaptation_->ReferenceSpacing())
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_pos_.push_back(&(contact_particles_[k]->pos_0_));
			}
		}
		//=================================================================================================//
		// template<int intOpt>
		void ShellContactDensity::Interaction(size_t index_i, Real dt)
		{
			/** Shell contact interaction. */
			Real sigma = 0.0;
			/** Abscissas and weights for Gauss-Legendre quadrature for n=3 nodes */
			Real x_0 = 0.774596669241483377035853079956;
			Real x_1 = 0.000000000000000000000000000000;
			Real x_2 = -x_0;

			Real w_0 = 0.555555555555555555555555555556;
			Real w_1 = 0.888888888888888888888888888889;
			Real w_2 = w_0;

			boundary_factor_ = material_->ReferenceDensity() / (1.3*5.2*spacing_ref_*spacing_ref_*kernel_->W(0., Vecd(0., 0.)));

			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &contact_pos_k = *(contact_pos_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					
					const Real dp_2 = spacing_ref_/2.;
					const Vecd contact_pos_j = contact_pos_k[contact_neighborhood.j_[n]];
					
					const Vecd dp_2_t_0 = pos_n_[index_i] - Vecd(dp_2*x_0, 0.) - contact_pos_j;
					const Vecd dp_2_t_1 = pos_n_[index_i] - Vecd(dp_2*x_1, 0.) - contact_pos_j;
					const Vecd dp_2_t_2 = pos_n_[index_i] - Vecd(dp_2*x_2, 0.) - contact_pos_j;

					const Real W_rij_t_0 = kernel_->W(dp_2_t_0.norm(), dp_2_t_0);
					const Real W_rij_t_1 = kernel_->W(dp_2_t_1.norm(), dp_2_t_1);
					const Real W_rij_t_2 = kernel_->W(dp_2_t_2.norm(), dp_2_t_2);

					sigma  += (w_0 * W_rij_t_0 + w_1 * W_rij_t_1 + w_2 * W_rij_t_2) * dp_2 * boundary_factor_ * 1.3 * spacing_ref_;
				}				
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		SelfContactForce::
			SelfContactForce(SolidBodyRelationSelfContact* self_contact_relation)
			: PartInteractionDynamicsByParticle(self_contact_relation->sph_body_,
												&self_contact_relation->body_surface_layer_),
			  SolidDataInner(self_contact_relation),
			  mass_(particles_->mass_), contact_density_(particles_->contact_density_), Vol_(particles_->Vol_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_), contact_force_(particles_->contact_force_) {}
		//=================================================================================================//
		void SelfContactForce::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real p_i = contact_density_[index_i] * material_->ContactStiffness();

			/** Inner interaction. */
			Vecd force(0.0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real p_star = 0.5 * (p_i + contact_density_[index_j] * material_->ContactStiffness());
				//force to mimic pressure
				force -= 2.0 * p_star * e_ij * Vol_i * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}
			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		DynamicSelfContactForce::
			DynamicSelfContactForce(SolidBodyRelationSelfContact* self_contact_relation, Real penalty_strength)
			: PartInteractionDynamicsByParticle(self_contact_relation->sph_body_,
				&self_contact_relation->body_surface_layer_),
			SolidDataInner(self_contact_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_),
			vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			contact_force_(particles_->contact_force_), penalty_strength_(penalty_strength),
			contact_impedance_(material_->ReferenceDensity()* sqrt(material_->ContactStiffness())),
			contact_reference_pressure_(material_->ReferenceDensity()* material_->ContactStiffness())
		{
			particle_spacing_j1_ = 1.0 / body_->particle_adaptation_->ReferenceSpacing();
			particle_spacing_ratio2_ =
				1.0 / (body_->particle_adaptation_->ReferenceSpacing() * particle_spacing_j1_);
			particle_spacing_ratio2_ *= 0.1 * particle_spacing_ratio2_;
		}
		//=================================================================================================//
		void DynamicSelfContactForce::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd vel_i = vel_n_[index_i];

			/** Contact interaction. */
			Vecd force(0.0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];

				Real impedance_p = contact_impedance_ * (SimTK::dot(vel_i - vel_n_[index_j], -e_ij));
				Real overlap = inner_neighborhood.r_ij_[n];
				Real delta = 2.0 * overlap * particle_spacing_j1_;
				Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2_ : 0.0;
				Real penalty_p = penalty_strength_* beta* overlap* contact_reference_pressure_;

				//force due to pressure
				force -= 2.0 * (impedance_p + penalty_p) * e_ij *
					Vol_i * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}
			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForce::ContactForce(SolidBodyRelationContact *solid_body_contact_relation)
			: PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
												&solid_body_contact_relation->body_surface_layer_),
			  ContactDynamicsData(solid_body_contact_relation),
			  contact_density_(particles_->contact_density_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  contact_force_(particles_->contact_force_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_contact_density_.push_back(&(contact_particles_[k]->contact_density_));
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
			
		}
		//=================================================================================================//
		void ContactForce::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real p_i = contact_density_[index_i] * material_->ContactStiffness();
			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_density_k = *(contact_contact_density_[k]);
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real> &mass_k = *(contact_mass_[k]);
				Solid *solid_k = contact_material_[k];

				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					Real p_star = 0.5 * (p_i + contact_density_k[index_j] * solid_k->ContactStiffness());
					//force due to pressure
					force -= 2.0 * p_star * e_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}
			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		DynamicContactForce::
			DynamicContactForce(SolidBodyRelationContact *solid_body_contact_relation, Real penalty_strength)
			: PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
												&solid_body_contact_relation->body_surface_layer_),
			  ContactDynamicsData(solid_body_contact_relation),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			  contact_force_(particles_->contact_force_), penalty_strength_(penalty_strength)
		{
			Real impedance = material_->ReferenceDensity() * sqrt(material_->ContactStiffness());
			Real reference_pressure = material_->ReferenceDensity() * material_->ContactStiffness();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
				Real contact_impedance =
					contact_material_[k]->ReferenceDensity() * sqrt(contact_material_[k]->ContactStiffness());
				contact_impedance_.push_back(2.0 * impedance * contact_impedance / (impedance + contact_impedance));
				Real contact_reference_pressure =
					contact_material_[k]->ReferenceDensity() * contact_material_[k]->ContactStiffness();
				contact_reference_pressure_.push_back(2.0 * reference_pressure * contact_reference_pressure /
													  (reference_pressure + contact_reference_pressure));
			}
		}
		//=================================================================================================//
		void DynamicContactForce::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd vel_i = vel_n_[index_i];

			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->particle_adaptation_->ReferenceSpacing();
				Real particle_spacing_ratio2 =
					1.0 / (this->body_->particle_adaptation_->ReferenceSpacing() * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);

				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					Real impedance_p = 0.5 * contact_impedance_[k] * (SimTK::dot(vel_i - vel_n_k[index_j], -e_ij));
					Real overlap = contact_neighborhood.r_ij_[n];
					Real delta = 2.0 * overlap * particle_spacing_j1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty_p = penalty_strength_ * beta * overlap * contact_reference_pressure_[k];

					//force due to pressure
					force -= 2.0 * (impedance_p + penalty_p) * e_ij *
							 Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForceWithWall::
			ContactForceWithWall(SolidBodyRelationContact *solid_body_contact_relation, Real penalty_strength)
			: PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
												&solid_body_contact_relation->body_surface_layer_),
			  ContactDynamicsData(solid_body_contact_relation),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			  contact_force_(particles_->contact_force_), penalty_strength_(penalty_strength)
		{
			impedance_ = material_->ReferenceDensity() * sqrt(material_->ContactStiffness());
			reference_pressure_ = material_->ReferenceDensity() * material_->ContactStiffness();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
				contact_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void ContactForceWithWall::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd vel_i = vel_n_[index_i];

			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->particle_adaptation_->ReferenceSpacing();
				Real particle_spacing_ratio2 =
					1.0 / (body_->particle_adaptation_->ReferenceSpacing() * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd> &n_k = *(contact_n_[k]);
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);

				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];
					Vecd n_k_j = n_k[index_j];

					Real impedance_p = 0.5 * impedance_ * (SimTK::dot(vel_i - vel_n_k[index_j], -n_k_j));
					Real overlap = contact_neighborhood.r_ij_[n] * SimTK::dot(n_k_j, e_ij);
					Real delta = 2.0 * overlap * particle_spacing_j1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty_p = penalty_strength_ * beta * overlap * reference_pressure_;

					//force due to pressure
					force -= 2.0 * (impedance_p + penalty_p) * dot(e_ij, n_k_j) *
							 n_k_j * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
	}
}
