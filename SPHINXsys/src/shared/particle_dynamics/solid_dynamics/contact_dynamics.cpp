/**
 * @file 	contact_dynamics.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "contact_dynamics.h"

#ifdef max
#undef max
#endif

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		SelfContactDensitySummation::
			SelfContactDensitySummation(SelfSurfaceContactRelation &self_contact_relation)
			: LocalDynamics(self_contact_relation.sph_body_), SolidDataInner(self_contact_relation),
			  mass_(particles_->mass_)
		{
			particles_->registerVariable(self_contact_density_, "SelfContactDensity");
			Real dp_1 = self_contact_relation.sph_body_.sph_adaptation_->ReferenceSpacing();
			offset_W_ij_ = self_contact_relation.sph_body_.sph_adaptation_->getKernel()->W(dp_1, Vecd(0.0));
		}
		//=================================================================================================//
		void SelfContactDensitySummation::interaction(size_t index_i, Real dt)
		{
			Real sigma = 0.0;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				Real corrected_W_ij = std::max(inner_neighborhood.W_ij_[n] - offset_W_ij_, 0.0);
				sigma += corrected_W_ij * mass_[inner_neighborhood.j_[n]];
			}
			self_contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ContactDensitySummation::
			ContactDensitySummation(SurfaceContactRelation &solid_body_contact_relation)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactDynamicsData(solid_body_contact_relation), mass_(particles_->mass_),
			  offset_W_ij_(StdVec<Real>(contact_configuration_.size(), 0.0))
		{
			particles_->registerVariable(contact_density_, "ContactDensity");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}

			// we modify the default formulation by an offset, so that exactly touching bodies produce 0 initial force
			// subtract summation of the kernel function of 2 particles at 1 particle distance, and if the result is negative, we take 0
			// different resolution: distance = 0.5 * dp1 + 0.5 * dp2
			// dp1, dp2 half reference spacing
			Real dp_1 = solid_body_contact_relation.sph_body_.sph_adaptation_->ReferenceSpacing();
			// different resolution: distance = 0.5 * dp1 + 0.5 * dp2
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real dp_2 = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
				Real distance = 0.5 * dp_1 + 0.5 * dp_2;
				offset_W_ij_[k] = solid_body_contact_relation.sph_body_.sph_adaptation_->getKernel()->W(distance, Vecd(0.0));
			}
		}
		//=================================================================================================//
		void ContactDensitySummation::interaction(size_t index_i, Real dt)
		{
			/** Contact interaction. */
			Real sigma = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], 0.0);
					sigma += corrected_W_ij * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ShellContactDensity::ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactDynamicsData(solid_body_contact_relation),
			  solid_(particles_->solid_), pos_(particles_->pos_),
			  kernel_(solid_body_contact_relation.sph_body_.sph_adaptation_->getKernel()),
			  spacing_ref_(solid_body_contact_relation.sph_body_.sph_adaptation_->ReferenceSpacing())
		{
			particles_->registerVariable(contact_density_, "ContactDensity");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_pos_.push_back(&(contact_particles_[k]->pos_));
			}
		}
		//=================================================================================================//
		void ShellContactDensity::interaction(size_t index_i, Real dt)
		{
			/** shell contact interaction. */
			Real sigma = 0.0;
			const int dimension = Vecd(0).size();
			/** a calibration factor to avoid particle penetration into shell structure */
			boundary_factor_ = solid_.ReferenceDensity() /
							   (kernel_->SmoothingLength() * kernel_->W0(Vecd(0.)) * Pi * std::pow(kernel_->CutOffRadius(), dimension - 1));

			const Real dp_2 = 0.5 * spacing_ref_;

			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &contact_pos_k = *(contact_pos_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					const Vecd contact_pos_j = contact_pos_k[contact_neighborhood.j_[n]];

					const Vecd dp_2_t_0 = pos_[index_i] - Vecd(dp_2 * x_0, dp_2 * x_0) - contact_pos_j;
					const Vecd dp_2_t_1 = pos_[index_i] - Vecd(dp_2 * x_1, dp_2 * x_1) - contact_pos_j;
					const Vecd dp_2_t_2 = pos_[index_i] - Vecd(dp_2 * x_2, dp_2 * x_2) - contact_pos_j;

					const Real W_rij_t_0 = kernel_->W(dp_2_t_0.norm(), dp_2_t_0);
					const Real W_rij_t_1 = kernel_->W(dp_2_t_1.norm(), dp_2_t_1);
					const Real W_rij_t_2 = kernel_->W(dp_2_t_2.norm(), dp_2_t_2);

					sigma += (w_0 * W_rij_t_0 + w_1 * W_rij_t_1 + w_2 * W_rij_t_2) * dp_2;
				}
			}
			contact_density_[index_i] = sigma * boundary_factor_ * kernel_->SmoothingLength();
		}
		//=================================================================================================//
		SelfContactForce::
			SelfContactForce(SelfSurfaceContactRelation &self_contact_relation)
			: LocalDynamics(self_contact_relation.sph_body_),
			  SolidDataInner(self_contact_relation),
			  solid_(particles_->solid_), mass_(particles_->mass_),
			  self_contact_density_(*particles_->getVariableByName<Real>("SelfContactDensity")),
			  Vol_(particles_->Vol_), acc_prior_(particles_->acc_prior_),
			  contact_impedance_(solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness())),
			  vel_(particles_->vel_) {}
		//=================================================================================================//
		void SelfContactForce::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd vel_i = vel_[index_i];
			Real p_i = self_contact_density_[index_i] * solid_.ContactStiffness();

			/** Inner interaction. */
			Vecd force(0.0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real p_star = 0.5 * (p_i + self_contact_density_[index_j] * solid_.ContactStiffness());
				Real impedance_p = 0.5 * contact_impedance_ * (SimTK::dot(vel_i - vel_[index_j], -e_ij));
				// force to mimic pressure
				force -= 2.0 * (p_star + impedance_p) * e_ij * Vol_i * inner_neighborhood.dW_ijV_j_[n];
			}
			acc_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForce::ContactForce(SurfaceContactRelation &solid_body_contact_relation)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactDynamicsData(solid_body_contact_relation),
			  solid_(particles_->solid_),
			  contact_density_(*particles_->getVariableByName<Real>("ContactDensity")),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  acc_prior_(particles_->acc_prior_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_solids_.push_back(&contact_particles_[k]->solid_);
				contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("ContactDensity"));
			}
		}
		//=================================================================================================//
		void ContactForce::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real p_i = contact_density_[index_i] * solid_.ContactStiffness();
			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_density_k = *(contact_contact_density_[k]);
				Solid *solid_k = contact_solids_[k];

				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					Real p_star = 0.5 * (p_i + contact_density_k[index_j] * solid_k->ContactStiffness());
					// force due to pressure
					force -= 2.0 * p_star * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
				}
			}
			acc_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForceFromWall::ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactWithWallData(solid_body_contact_relation), solid_(particles_->solid_),
			  contact_density_(*particles_->getVariableByName<Real>("ContactDensity")),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  acc_prior_(particles_->acc_prior_) {}
		//=================================================================================================//
		void ContactForceFromWall::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real p_i = contact_density_[index_i] * solid_.ContactStiffness();
			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					// force due to pressure
					force -= 2.0 * p_i * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
				}
			}
			acc_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForceToWall::ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactDynamicsData(solid_body_contact_relation),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  acc_prior_(particles_->acc_prior_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_solids_.push_back(&contact_particles_[k]->solid_);
				contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("ContactDensity"));
			}
		}
		//=================================================================================================//
		void ContactForceToWall::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_density_k = *(contact_contact_density_[k]);
				Solid *solid_k = contact_solids_[k];

				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					Real p_star = contact_density_k[index_j] * solid_k->ContactStiffness();
					// force due to pressure
					force -= 2.0 * p_star * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
				}
			}
			acc_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		PairwiseFrictionFromWall::
			PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
			: LocalDynamics(contact_relation.sph_body_), ContactWithWallData(contact_relation),
			  eta_(eta), Vol_(particles_->Vol_), mass_(particles_->mass_),
			  vel_(particles_->vel_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				wall_vel_n_.push_back(&contact_particles_[k]->vel_);
				wall_n_.push_back(&contact_particles_[k]->n_);
			}
		}
		//=================================================================================================//
		void PairwiseFrictionFromWall::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real mass_i = mass_[index_i];
			Vecd &vel_i = vel_[index_i];

			std::array<Real, MaximumNeighborhoodSize> parameter_b;

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &vel_k = *(wall_vel_n_[k]);
				StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				// forward sweep
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd &e_ij = contact_neighborhood.e_ij_[n];

					parameter_b[n] = eta_ * contact_neighborhood.dW_ijV_j_[n] * Vol_i * dt / contact_neighborhood.r_ij_[n];

					// only update particle i
					Vecd vel_derivative = (vel_i - vel_k[index_j]);
					Vecd n_j = SimTK::dot(e_ij, n_k[index_j]) > 0.0 ? n_k[index_j] : -1.0 * n_k[index_j];
					vel_derivative -= SMAX(0.0, SimTK::dot(vel_derivative, n_j)) * n_j;
					vel_i += parameter_b[n] * vel_derivative / (mass_i - 2.0 * parameter_b[n]);
				}
				// backward sweep
				for (size_t n = contact_neighborhood.current_size_; n != 0; --n)
				{
					size_t index_j = contact_neighborhood.j_[n - 1];
					Vecd &e_ij = contact_neighborhood.e_ij_[n];

					// only update particle i
					Vecd vel_derivative = (vel_i - vel_k[index_j]);
					Vecd n_j = SimTK::dot(e_ij, n_k[index_j]) > 0.0 ? n_k[index_j] : -1.0 * n_k[index_j];
					vel_derivative -= SMAX(0.0, SimTK::dot(vel_derivative, n_j)) * n_j;
					vel_i += parameter_b[n - 1] * vel_derivative / (mass_i - 2.0 * parameter_b[n - 1]);
				}
			}
		}
		//=================================================================================================//
		DynamicContactForceWithWall::
			DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength)
			: LocalDynamics(solid_body_contact_relation.sph_body_),
			  ContactDynamicsData(solid_body_contact_relation),
			  solid_(particles_->solid_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_),
			  vel_(particles_->vel_), acc_prior_(particles_->acc_prior_),
			  penalty_strength_(penalty_strength)
		{
			impedance_ = solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness());
			reference_pressure_ = solid_.ReferenceDensity() * solid_.ContactStiffness();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
				contact_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void DynamicContactForceWithWall::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd vel_i = vel_[index_i];

			/** Contact interaction. */
			Vecd force(0.0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
				Real particle_spacing_ratio2 =
					1.0 / (sph_body_.sph_adaptation_->ReferenceSpacing() * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

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
					Real penalty_p = penalty_strength_ * beta * fabs(overlap) * reference_pressure_;

					// force due to pressure
					force -= 2.0 * (impedance_p + penalty_p) * dot(e_ij, n_k_j) *
							 n_k_j * Vol_i * contact_neighborhood.dW_ijV_j_[n];
				}
			}

			acc_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
	}
}
