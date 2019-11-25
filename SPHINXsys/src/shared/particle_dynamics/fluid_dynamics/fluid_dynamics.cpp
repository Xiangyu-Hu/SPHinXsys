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
		void DensityBySummation::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real Vol_0_i = base_particle_data_i.Vol_0_;

			/** Inner interaction. */
			Real sigma = W0_;
			StdVec<NeighboringParticle>  &inner_neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = inner_neighors[n];

				sigma += neighboring_particle.W_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>  &contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle &neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					sigma += neighboring_particle.W_ij_*base_particle_data_j.Vol_0_ / Vol_0_i;
				}
			}

			/** Particle summation. */
			base_particle_data_i.sigma_ = sigma;

			/** Particle updating. */
			ParticleUpdate(index_particle_i, dt);
		}
//=================================================================================================//
		void DensityBySummation::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			fluid_data_i.rho_n_ = base_particle_data_i.sigma_ * fluid_data_i.rho_0_ / base_particle_data_i.sigma_0_;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
		}
//=================================================================================================//
		void DensityBySummationFreeSurface::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_sum = base_particle_data_i.sigma_ * fluid_data_i.rho_0_ / base_particle_data_i.sigma_0_;
			fluid_data_i.rho_n_ = SMAX(rho_sum, fluid_data_i.rho_0_);
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
		}
//=================================================================================================//
		void DivergenceCorrection::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			/** Inner interaction. */
			Real div_correction = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				div_correction -= neighboring_particle.dW_ij_ * neighboring_particle.r_ij_ * base_particle_data_j.Vol_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>  &contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle &neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					div_correction -= neighboring_particle.dW_ij_ * neighboring_particle.r_ij_ * base_particle_data_j.Vol_;
				}
			}

			/** Particle summation. */
			Real div_correction_1 = div_correction / dimension_;
			fluid_data_i.div_correction_
				= 1.0 / (div_correction_1 + 0.1*(1.0 - div_correction_1)*(1.0 - div_correction_1));
		}
//=================================================================================================//
		void ComputingViscousAcceleration::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				acceleration += 2.0*mu_*vel_detivative*base_particle_data_j.Vol_ / rho_i
					*neighboring_particle.dW_ij_;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>  &contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle &neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//viscous force with a simple wall model for high-Reynolds number flow
					Vecd vel_detivative = 2.0*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_)
						/ (neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
					Real vel_difference = 0.03*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_).norm()
						* neighboring_particle.r_ij_;
					acceleration += 2.0*SMAX(mu_, rho_i * vel_difference)
						*vel_detivative*base_particle_data_j.Vol_ / rho_i
						*neighboring_particle.dW_ij_;
				}
			}

			/** Particle summation. */
			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
//=================================================================================================//
		void ComputingAngularConservitiveViscousAcceleration::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
				Real v_r_ij = dot(base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_, 
					neighboring_particle.r_ij_ * neighboring_particle.e_ij_);
				Real eta_ij = 16.0 * (mu_ * mu_) * v_r_ij / (mu_ + mu_) / 
					(neighboring_particle.r_ij_ * neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				acceleration += eta_ij * base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>  &contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle &neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData &base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData &solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
				Real v_r_ij = dot(base_particle_data_i.vel_n_ -  solid_data_j.vel_ave_, 
					neighboring_particle.r_ij_ * neighboring_particle.e_ij_);
				Real vel_difference = 0.03*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_).norm()
					* neighboring_particle.r_ij_;
				Real eta_ij = 8.0 * SMAX(mu_, fluid_data_i.rho_n_*vel_difference) * v_r_ij / 
					(neighboring_particle.r_ij_ * neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				acceleration += eta_ij * base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				}
			}

			/** Particle summation. */
			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
//=================================================================================================//
		void TransportVelocityStress::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			/** Inner interaction. */
			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//exra stress
				acceleration += 0.5*dt*((fluid_data_i.rho_n_ * base_particle_data_i.vel_n_
					* neighboring_particle.dW_ij_ * dot(fluid_data_i.dvel_dt_trans_, neighboring_particle.e_ij_))
					+ (fluid_data_j.rho_n_*base_particle_data_j.vel_n_
					* neighboring_particle.dW_ij_ * dot(fluid_data_j.dvel_dt_trans_, neighboring_particle.e_ij_)))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>& contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle& neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//exra stress
					acceleration += 0.5 * dt * (fluid_data_i.rho_n_ * base_particle_data_i.vel_n_
						* neighboring_particle.dW_ij_ * dot(fluid_data_i.dvel_dt_trans_, neighboring_particle.e_ij_))
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
			p_background_ =  8.0 * density * speed_max * speed_max;
		}
//=================================================================================================//
		void TransportVelocityCorrection::ParticleInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;

			/** Inner interaction. */
			Vecd acceleration_trans(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_*base_particle_data_j.Vol_ 
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_ / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				StdVec<NeighboringParticle>& contact_neighors
					= (*current_interacting_configuration_[k])[index_particle_i];
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					NeighboringParticle& neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle.j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					SolidParticleData& solid_data_j
						= (*interacting_particles_[k]).solid_body_data_[index_particle_j];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * base_particle_data_j.Vol_
						* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_ / rho_i;
				}
			}

			/** Particle summation. */
			fluid_data_i.dvel_dt_trans_ = acceleration_trans;
			base_particle_data_i.pos_n_ += acceleration_trans * dt*dt*0.5;
		}
//=================================================================================================//
		TotalMechanicalEnergy::TotalMechanicalEnergy(FluidBody* body, ExternalForce *external_force)
			: WeaklyCompressibleFluidDynamicsSum<Real>(body)
		{
			initial_reference_ = 0.0;
			average_farctor_ = 1.0;// / Real(body_->number_of_real_particles_);
			potential_ = external_force->InducedAcceleration().norm();
		}
//=================================================================================================//
		Real TotalMechanicalEnergy::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			return average_farctor_ * (
				0.5 * fluid_data_i.mass_* base_particle_data_i.vel_n_.normSqr()
				+ fluid_data_i.mass_*potential_*base_particle_data_i.pos_n_[1]);
		}
//=================================================================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(FluidBody* body)
			: WeaklyCompressibleFluidDynamicsMaximum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time setep size due to linear viscosity
			Real rho_0 = material_->rho_0_;
			Real mu = material_->mu_;
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
		void PressureRelaxationFirstHalfRiemann::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd vel_i = base_particle_data_i.vel_n_;

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			Vecd acceleration_strong_form = base_particle_data_i.dvel_dt_others_;

			StdVec<NeighboringParticle>& neighors = (*current_inner_configuration_)[index_particle_i];

			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//low dissipation Riemann problem
				Real ul = dot(-neighboring_particle.e_ij_, vel_i);
				Real ur = dot(-neighboring_particle.e_ij_, base_particle_data_j.vel_n_);
				Real p_star = material_->RiemannSolverForPressure(rho_i, fluid_data_j.rho_n_,
					p_i, fluid_data_j.p_, ul, ur);

				acceleration -= 2.0 * p_star * base_particle_data_j.Vol_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_ / rho_i;
				/** acceleration in strong form */
				acceleration_strong_form += (p_i - p_star) * base_particle_data_j.Vol_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_ / rho_i;
			}

			fluid_data_i.dvel_dt_inner_ = acceleration_strong_form;
			base_particle_data_i.dvel_dt_ = acceleration;
		}
//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd dvel_dt_others_i = base_particle_data_i.dvel_dt_others_;
			Real particle_spacing_j1 = 1.0 / interacting_bodies_[interacting_body_index]->particle_spacing_;
			Real particle_spacing_ratio2 = 1.0 / (body_->particle_spacing_ * particle_spacing_j1);
			particle_spacing_ratio2 *= 0.1*particle_spacing_ratio2;

			Vecd acceleration(0);
			StdVec<NeighboringParticle>& neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData& solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Vecd n_j = solid_data_j.n_;
				Vecd e_ij = neighboring_particle.e_ij_;
				Real face_wall_external_acceleration
					= dot((dvel_dt_others_i - solid_data_j.dvel_dt_ave_), - e_ij);
				Real p_star = p_i + 0.5 * rho_i * neighboring_particle.r_ij_ * SMAX(0.0, face_wall_external_acceleration);

				/** penalty method to prevent particle running into boundary */
				Real projection = dot(e_ij, n_j);
				Real delta = 2.0 * projection * neighboring_particle.r_ij_ * particle_spacing_j1;
				Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
				Real penalty = beta * projection * fabs(p_star);

				//pressure force
				acceleration -= 2.0 * (p_star * e_ij + penalty * n_j)
								* base_particle_data_j.Vol_ * neighboring_particle.dW_ij_  / rho_i;
			}

			base_particle_data_i.dvel_dt_ += acceleration;
		}
//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_ * dt;
		}
//=================================================================================================//
		void PressureRelaxationFirstHalf::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				Real p_star = (p_i*fluid_data_j.rho_n_ + fluid_data_j.p_*rho_i)
					/ (rho_i + fluid_data_j.rho_n_);

				//pressure force
				acceleration -= 2.0*p_star*base_particle_data_j.Vol_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_ / rho_i;
			}

			base_particle_data_i.dvel_dt_ = acceleration;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd vel_i = base_particle_data_i.vel_n_;

			Real density_change_rate = 0.0;
			StdVec<NeighboringParticle>& neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//low dissipation Riemann problem
				Real ul = dot(-neighboring_particle.e_ij_, vel_i);
				Real ur = dot(-neighboring_particle.e_ij_, base_particle_data_j.vel_n_);
				Real u_star = material_->RiemannSolverForVelocity(rho_i, fluid_data_j.rho_n_,
					p_i, fluid_data_j.p_, ul, ur);
				Vecd vel_star = 0.5 * (vel_i + base_particle_data_j.vel_n_) 
							  - neighboring_particle.e_ij_ * (u_star - 0.5 * (ul + ur));
				//simplify the equation to cancle vector operations
				density_change_rate += 2.0 * rho_i * base_particle_data_j.Vol_
					* dot(vel_i - vel_star, neighboring_particle.e_ij_) * neighboring_particle.dW_ij_;
			}
			fluid_data_i.drho_dt_ = fluid_data_i.div_correction_ * density_change_rate;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Real p_i = fluid_data_i.p_;
			Vecd vel_i = base_particle_data_i.vel_n_;
			Vecd pos_i = base_particle_data_i.pos_n_;
			Vecd dvel_dt_others_i = base_particle_data_i.dvel_dt_others_;

			Real density_change_rate = 0.0;
			StdVec<NeighboringParticle>& neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData& solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Vecd e_ij = neighboring_particle.e_ij_;
				Vecd vel_in_wall = 2.0 * solid_data_j.vel_ave_ - vel_i;
				Real face_wall_external_acceleration
					= dot((dvel_dt_others_i - solid_data_j.dvel_dt_ave_), e_ij);
				Real p_in_wall = p_i + rho_i * neighboring_particle.r_ij_ * SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

				//soliving Riemann problem along the wall normal direction
				Vecd n_j = solid_data_j.n_;
				Real ul = dot(-n_j, vel_i);
				Real ur = dot(-n_j, vel_in_wall);
				Real u_star = material_->RiemannSolverForVelocity(rho_i, rho_in_wall, p_i, p_in_wall, ul, ur);
				Vecd vel_star = 0.5 * (vel_i + vel_in_wall) - n_j * (u_star - 0.5 * (ul + ur));

				density_change_rate += 2.0 * rho_i * base_particle_data_j.Vol_
					* dot(vel_i - vel_star, neighboring_particle.e_ij_)
					* neighboring_particle.dW_ij_;
			}
			fluid_data_i.drho_dt_ += fluid_data_i.div_correction_ * density_change_rate;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::Update(size_t index_particle_i, Real dt)
		{
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_ * dt * 0.5;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalf::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Vecd vel_i = base_particle_data_i.vel_n_;

			Real density_change_rate = 0.0;
			StdVec<NeighboringParticle>& neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				density_change_rate += rho_i * base_particle_data_j.Vol_
					* dot(vel_i - base_particle_data_j.vel_n_, neighboring_particle.e_ij_)
					* neighboring_particle.dW_ij_;
			}
			fluid_data_i.drho_dt_ = fluid_data_i.div_correction_ * density_change_rate;
		}
//=================================================================================================//
		PressureRelaxationFirstHalfOldroyd_B
			::PressureRelaxationFirstHalfOldroyd_B(FluidBody* body, StdVec<SolidBody*> interacting_bodies)
			: PressureRelaxationFirstHalf(body, interacting_bodies) {
			viscoelastic_fluid_particles_
				= dynamic_cast<ViscoelasticFluidParticles*>(body->base_particles_->PointToThisObject());
		}
//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::Initialization(size_t index_particle_i, Real dt)
		{
			PressureRelaxationFirstHalf::Initialization(index_particle_i, dt);

			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i 
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			non_newtonian_fluid_data_i.tau_ += non_newtonian_fluid_data_i.dtau_dt_ * dt * 0.5;
		}
//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::InnerInteraction(size_t index_particle_i, Real dt)
		{
			PressureRelaxationFirstHalf::InnerInteraction(index_particle_i, dt);

			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Vecd acceleration(0);
			StdVec<NeighboringParticle>& neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				ViscoelasticFluidParticleData& non_newtonian_fluid_data_j
					= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_j];

				//elastic force
				acceleration -= (tau_i + non_newtonian_fluid_data_j.tau_)
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					* base_particle_data_j.Vol_ / rho_i;
			}
			base_particle_data_i.dvel_dt_ += acceleration;
		}
//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			PressureRelaxationFirstHalf::ContactInteraction(index_particle_i, interacting_body_index, dt);

			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Real rho_i = fluid_data_i.rho_n_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Vecd acceleration(0);
			StdVec<NeighboringParticle>& neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				//stress boundary condition 
				acceleration -= 2.0 * tau_i * neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					* base_particle_data_j.Vol_ / rho_i;
			}
			base_particle_data_i.dvel_dt_ += acceleration;
		}
//=================================================================================================//
		PressureRelaxationSecondHalfOldroyd_B
			::PressureRelaxationSecondHalfOldroyd_B(FluidBody* body, StdVec<SolidBody*> interacting_bodies)
			: PressureRelaxationSecondHalf(body, interacting_bodies) {
			viscoelastic_fluid_particles_
				= dynamic_cast<ViscoelasticFluidParticles*>(body->base_particles_->PointToThisObject());
			Oldroyd_B_Fluid *oldroy_b_fluid 
				= dynamic_cast<Oldroyd_B_Fluid*>(body->base_material_->PointToThisObject());
			mu_p_ = oldroy_b_fluid->mu_p_;
			lambda_ = oldroy_b_fluid->lambda_;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::InnerInteraction(size_t index_particle_i, Real dt)
		{
			PressureRelaxationSecondHalf::InnerInteraction(index_particle_i, dt);
			
			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Vecd vel_trans_i = fluid_data_i.vel_trans_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Matd stress_rate(0);
			StdVec<NeighboringParticle>& neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData& fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				Matd velocity_gradient = SimTK::outer((vel_trans_i - fluid_data_j.vel_trans_),
					neighboring_particle.dW_ij_ * neighboring_particle.e_ij_) * base_particle_data_j.Vol_;
				stress_rate -= ~velocity_gradient * tau_i + tau_i * velocity_gradient 
					+ tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}
			non_newtonian_fluid_data_i.dtau_dt_ = stress_rate;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			PressureRelaxationSecondHalf::ContactInteraction(index_particle_i, interacting_body_index, dt);

			FluidParticleData& fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData& non_newtonian_fluid_data_i
				= viscoelastic_fluid_particles_->viscoelastic_particle_data_[index_particle_i];
			Vecd vel_trans_i = fluid_data_i.vel_trans_;
			Matd tau_i = non_newtonian_fluid_data_i.tau_;

			Matd stress_rate(0);
			StdVec<NeighboringParticle>& neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle& neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData& base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData& solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Matd velocity_gradient = SimTK::outer((vel_trans_i - solid_data_j.vel_ave_),
					neighboring_particle.dW_ij_ * neighboring_particle.e_ij_) * base_particle_data_j.Vol_ * 2.0;
				stress_rate -= ~velocity_gradient * tau_i + tau_i * velocity_gradient 
					+ tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}
			non_newtonian_fluid_data_i.dtau_dt_ += stress_rate;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::Update(size_t index_particle_i, Real dt)
		{
			PressureRelaxationSecondHalf::Update(index_particle_i, dt);

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
				/** increase the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				if (body_->number_of_particles_ > particles_->base_particle_data_.size())
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** realize a particle */
				particles_->RealizeABufferParticle(body_->number_of_particles_ - 1, index_particle_i);
				base_particle_data_i.pos_n_[axis_] -= periodic_translation_[axis_];

			}
		}
//=================================================================================================//
		void EmitterInflowInjecting::CheckLowerBound(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			if (base_particle_data_i.pos_n_[axis_] < body_part_lower_bound_[axis_]) {
				/** increase the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				if (body_->number_of_particles_ > particles_->base_particle_data_.size())
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** realize a particle */
				particles_->RealizeABufferParticle(body_->number_of_particles_ - 1, index_particle_i);
				base_particle_data_i.pos_n_[axis_] += periodic_translation_[axis_];
			}
		}
//=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//