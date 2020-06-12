/**
 * @file 	fluid_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_dynamics.h"
//=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		void DensityBySummation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real Vol_0_i = base_particle_data_i.Vol_0_;

			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];

				sigma += neighboring_particle->W_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					sigma += neighboring_particle->W_ij_*base_particle_data_j.Vol_0_ / Vol_0_i;
				}
			}

			/** Update Particle. */
			UpdateDensity(index_particle_i, sigma);
		}
		//=================================================================================================//
		void DensityBySummation::UpdateDensity(size_t index_particle_i, Real sigma)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			fluid_data_i.rho_n_ = sigma * fluid_data_i.rho_0_ / base_particle_data_i.sigma_0_;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
		}
		//=================================================================================================//
		void DensityBySummationFreeSurface::UpdateDensity(size_t index_particle_i, Real sigma)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_sum = sigma * fluid_data_i.rho_0_ / base_particle_data_i.sigma_0_;
			fluid_data_i.rho_n_ = rho_sum + SMAX(0.0, (fluid_data_i.rho_n_ - rho_sum)) * fluid_data_i.rho_0_ / fluid_data_i.rho_n_;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
		}
		//=================================================================================================//
		void DivergenceCorrection::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			/** Inner interaction. */
			Real div_correction = 0.0;
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				div_correction -= neighboring_particle->dW_ij_ * neighboring_particle->r_ij_ * base_particle_data_j.Vol_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					div_correction -= neighboring_particle->dW_ij_ * neighboring_particle->r_ij_ * base_particle_data_j.Vol_;
				}
			}

			/** Particle summation. */
			Real div_correction_1 = div_correction / dimension_;
			fluid_data_i.div_correction_
				= 1.0 / (div_correction_1 + 0.1*(1.0 - div_correction_1)*(1.0 - div_correction_1));
		}
		//=================================================================================================//
		void ComputingViscousAcceleration::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Vecd& vel_i = base_particle_data_i.vel_n_;
			Real p_i = fluid_data_i.p_;

			/** Inner interaction. */
			Vecd acceleration(0);
			Vecd vel_derivative(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				Real dW_ij = neighboring_particle->dW_ij_;
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				Real Vol_j = base_particle_data_j.Vol_;

				//viscous force
				vel_derivative = (vel_i - base_particle_data_j.vel_n_)
					/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * Vol_j * dW_ij / rho_i;				/** acceleration in strong form */
			}


			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				/** computing the accelerations of near wall particles without considering wall. */
				Vecd acceleration_inner = acceleration;
				Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
				NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
				for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = inner_neighors[n];
					Vecd nablaW_ij = neighboring_particle->getNablaWij();
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
					FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];
					Real Vol_j = base_particle_data_j.Vol_;

					acceleration_inner += (p_i - fluid_data_j.p_) * Vol_j * nablaW_ij / rho_i;
					acceleration_inner += SimTK::outer(vel_i - base_particle_data_j.vel_n_, nablaW_ij) * vel_i * Vol_j;
				}

				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//viscous force with a simple wall model for high-Reynolds number flow
					vel_derivative = 2.0*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_)
						/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_);
					Real vel_difference =  0.0 * (base_particle_data_i.vel_n_ - solid_data_j.vel_ave_).norm()
						* neighboring_particle->r_ij_;
					acceleration += 2.0*SMAX(mu_, rho_i * vel_difference) * vel_derivative 
						* neighboring_particle->dW_ij_ * base_particle_data_j.Vol_ / rho_i;
				}
			}

			/** Particle summation. */
			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//=================================================================================================//
		void ComputingAngularConservativeViscousAcceleration::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real r_ij = neighboring_particle->r_ij_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
				Real v_r_ij = dot(base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_, r_ij * e_ij);
				Real eta_ij = 8.0 * mu_  * v_r_ij /	(r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle->dW_ij_ * e_ij;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					Vecd& e_ij = neighboring_particle->e_ij_;
					Real r_ij = neighboring_particle->r_ij_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
				Real v_r_ij = 2.0 * dot(base_particle_data_i.vel_n_ -  solid_data_j.vel_ave_, r_ij * e_ij);
				Real vel_difference = 0.0*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_).norm() * r_ij;
				Real eta_ij = 8.0 * SMAX(mu_, fluid_data_i.rho_n_*vel_difference) * v_r_ij / 
					(r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle->dW_ij_ * e_ij;
				}
			}

			/** Particle summation. */
			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//=================================================================================================//
		void TransportVelocityStress::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			/** Inner interaction. */
			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Real dW_ij = neighboring_particle->dW_ij_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real r_ij = neighboring_particle->r_ij_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//exra stress
				acceleration += 0.5*dt*((fluid_data_i.rho_n_ * base_particle_data_i.vel_n_
					* dW_ij * dot(fluid_data_i.dvel_dt_trans_, e_ij))
					+ (fluid_data_j.rho_n_*base_particle_data_j.vel_n_
					* dW_ij * dot(fluid_data_j.dvel_dt_trans_, e_ij)))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//exra stress
					acceleration += 0.5 * dt * (fluid_data_i.rho_n_ * base_particle_data_i.vel_n_
						* neighboring_particle->dW_ij_ * dot(fluid_data_i.dvel_dt_trans_, neighboring_particle->e_ij_))
						* base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
				}
			}

			/** Particle summation. */
			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//=================================================================================================//
		void TransportVelocityCorrection::SetupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->GetReferenceDensity();
			p_background_ =  10.0 * density * speed_max * speed_max;
		}
		//=================================================================================================//
		void TransportVelocityCorrection::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Vecd acceleration_trans(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_*base_particle_data_j.Vol_ 
					* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_ / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * base_particle_data_j.Vol_
						* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_ / rho_i;
				}
			}

			/** Particle summation. */
			fluid_data_i.dvel_dt_trans_ = acceleration_trans;
			base_particle_data_i.pos_n_ += acceleration_trans * dt*dt*0.5;
		}
		//=================================================================================================//
		TotalMechanicalEnergy::TotalMechanicalEnergy(FluidBody* body, Gravity* gravity)
			: WeaklyCompressibleFluidDynamicsSum<Real>(body), gravity_(gravity)
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real TotalMechanicalEnergy::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			return 0.5 * fluid_data_i.mass_* base_particle_data_i.vel_n_.normSqr()
				+ fluid_data_i.mass_* gravity_->getPotential(base_particle_data_i.pos_n_);
		}
		//=================================================================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(FluidBody* body)
			: WeaklyCompressibleFluidDynamicsMaximum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time setep size due to linear viscosity
			Real rho_0 = material_->GetReferenceDensity();
			Real mu = material_->getReferenceViscosity();
			initial_reference_ = mu / rho_0 / smoothing_length_;
		}
		//=================================================================================================//
		Real GetAcousticTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i	= particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			return material_->GetSoundSpeed(fluid_data_i.p_, fluid_data_i.rho_n_)
				+ base_particle_data_i.vel_n_.norm();
		}
		//=================================================================================================//
		Real GetAcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			particles_->signal_speed_max_ = reduced_value;
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + 1.0e-15);
		}
		//=================================================================================================//
		GetAdvectionTimeStepSize::GetAdvectionTimeStepSize(FluidBody* body, Real U_f)
			: GetAcousticTimeStepSize(body)
		{
			Real u_max = SMAX(initial_reference_, U_f);
			initial_reference_ = u_max * u_max;
		}
		//=================================================================================================//
		Real GetAdvectionTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			return base_particle_data_i.vel_n_.normSqr();
		}
//=================================================================================================//
		Real GetAdvectionTimeStepSize::OutputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			particles_->speed_max_ = speed_max;
			return 0.25 * smoothing_length_ / (speed_max + 1.0e-15);
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_ * dt * 0.5;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
			fluid_data_i.p_ = material_->GetPressure(fluid_data_i.rho_n_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd vel_i = base_particle_data_i.vel_n_;

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Real dW_ij = neighboring_particle->dW_ij_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real r_ij = neighboring_particle->r_ij_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				/** Solving Riemann problem or not. */
				Real p_star = getPStar(e_ij, vel_i, p_i, rho_i,
					base_particle_data_j.vel_n_, fluid_data_j.p_, fluid_data_j.rho_n_);
			
				acceleration -= 2.0 * p_star * base_particle_data_j.Vol_* dW_ij * e_ij / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Vecd dvel_dt_others_i = base_particle_data_i.dvel_dt_others_;
				Real particle_spacing_j1 = 1.0 / interacting_bodies_[k]->particle_spacing_;
				Real particle_spacing_ratio2 = 1.0 / (body_->particle_spacing_ * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					Vecd& n_j = solid_data_j.n_;
					Vecd& e_ij = neighboring_particle->e_ij_;
					Real dW_ij = neighboring_particle->dW_ij_;
					Real r_ij = neighboring_particle->r_ij_;
					Real face_wall_external_acceleration
						= dot((dvel_dt_others_i - solid_data_j.dvel_dt_ave_), -e_ij);
					Real p_star = p_i + 0.5 * rho_i * r_ij * SMAX(0.0, face_wall_external_acceleration);

					/** penalty method to prevent particle running into boundary */
					Real projection = dot(e_ij, n_j);
					Real delta = 2.0 * projection * r_ij * particle_spacing_j1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty = beta * projection * fabs(p_star);

					//pressure force
					acceleration -= 2.0 * (p_star * e_ij + penalty * n_j)
						* base_particle_data_j.Vol_ * dW_ij / rho_i;
				}
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//=================================================================================================//
		Real PressureRelaxationFirstHalfRiemann::getPStar(Vecd& e_ij,
			Vecd& vel_i, Real p_i, Real rho_i, Vecd& vel_j, Real p_j, Real rho_j)
		{
			//low dissipation Riemann problem
			Real ul = dot(-e_ij, vel_i);
			Real ur = dot(-e_ij, vel_j);
			return material_->RiemannSolverForPressure(rho_i, rho_j, p_i, p_j, ul, ur);
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_ * dt;
		}
		//=================================================================================================//
		Real PressureRelaxationFirstHalf::getPStar(Vecd& e_ij,
			Vecd& vel_i, Real p_i, Real rho_i, Vecd& vel_j, Real p_j, Real rho_j)
		{
			return (p_i * rho_j + p_j * rho_i) 	/ (rho_i + rho_j);;
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd vel_i = base_particle_data_i.vel_n_;

			Real density_change_rate = 0.0;
			Vecd vel_star(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real dW_ij = neighboring_particle->dW_ij_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				/** Solving Riemann problem or not. */
				vel_star = getVStar(e_ij, vel_i, p_i, rho_i,
					base_particle_data_j.vel_n_, fluid_data_j.p_, fluid_data_j.rho_n_);

				density_change_rate += 2.0 * rho_i * base_particle_data_j.Vol_
					* dot(vel_i - vel_star, e_ij) * dW_ij;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Vecd dvel_dt_others_i = base_particle_data_i.dvel_dt_others_;

				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					Vecd& e_ij = neighboring_particle->e_ij_;
					Real r_ij = neighboring_particle->r_ij_;
					Vecd vel_in_wall = 2.0 * solid_data_j.vel_ave_ - vel_i;
					Real dW_ij = neighboring_particle->dW_ij_;
					Real face_wall_external_acceleration
						= dot((dvel_dt_others_i - solid_data_j.dvel_dt_ave_), e_ij);
					Real p_in_wall = p_i + rho_i * r_ij * SMAX(0.0, face_wall_external_acceleration);
					Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

					//soliving Riemann or not
					vel_star = getVStar(solid_data_j.n_, vel_i, p_i, rho_i, vel_in_wall, p_in_wall, rho_in_wall);

					density_change_rate += 2.0 * rho_i * base_particle_data_j.Vol_
						* dot(vel_i - vel_star, e_ij) * dW_ij;
				}
			}

			fluid_data_i.drho_dt_ = fluid_data_i.div_correction_ * density_change_rate;
		}
		//=================================================================================================//
		Vecd PressureRelaxationSecondHalfRiemann::getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
			Vecd& vel_j, Real p_j, Real rho_j)
		{
			//low dissipation Riemann problem
			Real ul = dot(-e_ij, vel_i);
			Real ur = dot(-e_ij, vel_j);
			Real u_star = material_->RiemannSolverForVelocity(rho_i, rho_j,	p_i, p_j, ul, ur);
			return 0.5 * (vel_i +vel_j) - e_ij * (u_star - 0.5 * (ul + ur));
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::Update(size_t index_particle_i, Real dt)
		{
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_ * dt * 0.5;
		}
		//=================================================================================================//
		Vecd PressureRelaxationSecondHalf::getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
			Vecd& vel_j, Real p_j, Real rho_j)
		{
			return 0.5 * (vel_i + vel_j);
		}
		//=================================================================================================//
		PressureRelaxationFirstHalfOldroyd_B
			::PressureRelaxationFirstHalfOldroyd_B(FluidBody* body, StdVec<SolidBody*> interacting_bodies)
			: PressureRelaxationFirstHalfRiemann(body, interacting_bodies) {
			viscoelastic_fluid_particles_
				= dynamic_cast<ViscoelasticFluidParticles*>(body->base_particles_->PointToThisObject());
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::Initialization(size_t index_particle_i, Real dt)
		{
			PressureRelaxationFirstHalfRiemann::Initialization(index_particle_i, dt);

			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i 
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			non_newtonian_fluid_data_i.tau_ += non_newtonian_fluid_data_i.dtau_dt_ * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			PressureRelaxationFirstHalfRiemann::ComplexInteraction(index_particle_i, dt);

			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				ViscoelasticFluidParticleData& non_newtonian_fluid_data_j
					= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_j];

				//elastic force
				acceleration += (tau_i + non_newtonian_fluid_data_j.tau_)
					* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_
					* base_particle_data_j.Vol_ / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * neighboring_particle->dW_ij_ * neighboring_particle->e_ij_
						* base_particle_data_j.Vol_ / rho_i;
				}
			}

			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//=================================================================================================//
		PressureRelaxationSecondHalfOldroyd_B
			::PressureRelaxationSecondHalfOldroyd_B(FluidBody* body, StdVec<SolidBody*> interacting_bodies)
			: PressureRelaxationSecondHalfRiemann(body, interacting_bodies) {
			viscoelastic_fluid_particles_
				= dynamic_cast<ViscoelasticFluidParticles*>(body->base_particles_->PointToThisObject());
			Oldroyd_B_Fluid *oldroy_b_fluid 
				= dynamic_cast<Oldroyd_B_Fluid*>(body->base_particles_->base_material_->PointToThisObject());
			mu_p_ = oldroy_b_fluid->getReferencePloymericViscosity();
			lambda_ = oldroy_b_fluid->getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			PressureRelaxationSecondHalfRiemann::ComplexInteraction(index_particle_i, dt);
			
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Vecd vel_trans_i = fluid_data_i.vel_trans_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Matd stress_rate(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				Matd velocity_gradient = - SimTK::outer((vel_trans_i - fluid_data_j.vel_trans_),
					neighboring_particle->dW_ij_ * neighboring_particle->e_ij_) * base_particle_data_j.Vol_;
				stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient 
					- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					Matd velocity_gradient = - SimTK::outer((vel_trans_i - solid_data_j.vel_ave_),
						neighboring_particle->dW_ij_ * neighboring_particle->e_ij_) * base_particle_data_j.Vol_ * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient
						- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}

			non_newtonian_fluid_data_i.dtau_dt_ = stress_rate;
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::Update(size_t index_particle_i, Real dt)
		{
			PressureRelaxationSecondHalfRiemann::Update(index_particle_i, dt);

			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			non_newtonian_fluid_data_i.tau_ += non_newtonian_fluid_data_i.dtau_dt_ * dt * 0.5;
		}
		//=================================================================================================//
		void InflowBoundaryCondition
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.vel_n_ = base_particle_data_i.vel_n_*(1.0 - constrain_strength_)
				+ constrain_strength_*GetInflowVelocity(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_);
		}
		//=================================================================================================//
		void EmitterInflowCondition
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			/** constrain the state*/
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_particle_data_i
				= particles_->fluid_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_
				= GetInflowVelocity(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_);
			fluid_particle_data_i.rho_n_ = fluid_particle_data_i.rho_0_;
			fluid_particle_data_i.p_ = material_->GetPressure(fluid_particle_data_i.rho_n_);
		}
		//=================================================================================================//
		EmitterInflowInjecting
			::EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_size, int axis_direction, bool positive)
			: WeaklyCompressibleFluidConstraintByParticle(body, body_part), 
			axis_(axis_direction), periodic_translation_(0), body_buffer_size_(body_buffer_size) 
		{
			body_part->GetRegion()->regionbound(body_part_lower_bound_, body_part_upper_bound_);
			periodic_translation_[axis_] = body_part_upper_bound_[axis_] - body_part_lower_bound_[axis_];
			size_t total_body_buffer_particles = constrained_particles_.size() * body_buffer_size_;
			for (size_t i = 0; i < total_body_buffer_particles; ++i)
			{
				particles_->AddABufferParticle();
			}
			particles_->real_particles_bound_ += total_body_buffer_particles;
			body_->AllocateConfigurationMemoriesForBodyBuffer(total_body_buffer_particles);

			checking_bound_ = positive ?
				std::bind(&EmitterInflowInjecting::CheckUpperBound, this, _1, _2)
				: std::bind(&EmitterInflowInjecting::CheckLowerBound, this, _1, _2);
		}
		//=================================================================================================//
		void EmitterInflowInjecting::CheckUpperBound(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			if (base_particle_data_i.pos_n_[axis_] > body_part_upper_bound_[axis_]) {
				if (body_->number_of_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->CopyFromAnotherParticle(body_->number_of_particles_, index_particle_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				/** Periodic bounding. */
				base_particle_data_i.pos_n_[axis_] -= periodic_translation_[axis_];

			}
		}
		//=================================================================================================//
		void EmitterInflowInjecting::CheckLowerBound(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			if (base_particle_data_i.pos_n_[axis_] < body_part_lower_bound_[axis_]) {
				if (body_->number_of_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->CopyFromAnotherParticle(body_->number_of_particles_, index_particle_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				base_particle_data_i.pos_n_[axis_] += periodic_translation_[axis_];
			}
		}
		//=================================================================================================//
	    void ImplicitComputingViscousAcceleration::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.dvel_dt_ = 0.0;
		}
		//=================================================================================================//
	    void ImplicitComputingViscousAcceleration::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Real A_sum(0.0);
			Real A_square_sum(0.0);
			Vecd A_v_j_sum(0.0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				/** Viscous force. */
				Real A_ij = dt * 2.0 * mu_ * base_particle_data_j.Vol_ * neighboring_particle->dW_ij_ 
					/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_) / rho_i;
				A_sum += A_ij;
				A_square_sum += A_ij * A_ij;
				A_v_j_sum += A_ij * base_particle_data_j.vel_n_;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];
					/** Viscous force. */
					Real A_ij = dt * 2.0 * mu_ * base_particle_data_j.Vol_ * neighboring_particle->dW_ij_ 
						/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_) / rho_i;
					A_sum += A_ij;
					A_square_sum += A_ij * A_ij;
					A_v_j_sum += A_ij * solid_data_j.vel_ave_;
				}
			}
			/** Compute the lambda. */
			Vecd lambda = (base_particle_data_i.vel_n_ * A_sum - A_v_j_sum) / ((1.0 - A_sum) * (1.0 - A_sum) + A_square_sum) ;
			/** Force calculation. */
			/** Inner interaction. */
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				/** Viscous force. */
				Real A_ij = dt * 2.0 * mu_ * base_particle_data_j.Vol_ * neighboring_particle->dW_ij_ 
					/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_) / rho_i;
				base_particle_data_j.dvel_dt_ += A_ij * lambda;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];
					/** Viscous force. */
					Real A_ij = 2.0 * mu_ * base_particle_data_j.Vol_ * neighboring_particle->dW_ij_ 
						/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_) / rho_i;
					base_particle_data_j.dvel_dt_ += A_ij * lambda;
				}
			}

			base_particle_data_i.dvel_dt_ += (1.0 - A_sum) * lambda;
		}
		//=================================================================================================//
	    void ImplicitComputingViscousAcceleration::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_ * dt;
		}
		//=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//