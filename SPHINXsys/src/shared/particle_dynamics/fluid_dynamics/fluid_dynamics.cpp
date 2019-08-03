#include "fluid_dynamics.h"

using namespace std;

namespace SPH
{
	namespace fluid_dynamics
	{
		//===========================================================//
		void WeaklyCompressibleFluidInitialCondition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_(0);
			base_particle_data_i.dvel_dt_(0);

			fluid_data_i.p_ = 0.0;
			fluid_data_i.vel_trans_(0);
			fluid_data_i.rho_0_
				= material_->ReinitializeRho(fluid_data_i.p_);
			fluid_data_i.rho_n_ = fluid_data_i.rho_0_;
			fluid_data_i.mass_
				= fluid_data_i.rho_0_*base_particle_data_i.Vol_;
		}
		//===========================================================//
		void InitialNumberDensity::InnerInteraction(size_t index_particle_i, Real dt)
		{
			Real sigma = W0_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;

				sigma += neighors[n].W_ij_;
			}

			particles_->fluid_particle_data_[index_particle_i].sigma_0_ = sigma;
		}
		//===========================================================//
		void InitialNumberDensity::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real sigma = 0.0;
			StdVec<NeighboringParticle>  &neighors 
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				sigma += neighors[n].W_ij_*base_particle_data_j.Vol_0_/ base_particle_data_i.Vol_0_;
			}

			particles_->fluid_particle_data_[index_particle_i].sigma_0_ += sigma;
		}
		//===========================================================//
		void DensityBySummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			//initial value for summation
			Real sigma = W0_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				sigma += neighboring_particle.W_ij_;
			}

			fluid_data_i.rho_n_ = sigma * fluid_data_i.rho_0_ / fluid_data_i.sigma_0_;
		}
		//===========================================================//
		void DensityBySummation::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real sigma = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				sigma += neighboring_particle.W_ij_*base_particle_data_j.Vol_0_ / base_particle_data_i.Vol_0_;
			}

			fluid_data_i.rho_n_ += sigma * fluid_data_i.rho_0_ / fluid_data_i.sigma_0_;
		}
		//===========================================================//
		void  DensityBySummation::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
		}
		//===========================================================//
		void DensityBySummationFreeSurface::InnerInteraction(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.temp_real_ = fluid_data_i.rho_n_;

			//initial value for summation
			Real sigma = W0_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
	
				sigma += neighboring_particle.W_ij_;
			}

			fluid_data_i.rho_sum_ = sigma * fluid_data_i.rho_0_ / fluid_data_i.sigma_0_;
			fluid_data_i.rho_n_ = SMAX(fluid_data_i.temp_real_, fluid_data_i.rho_sum_);
		}
		//===========================================================//
		void DensityBySummationFreeSurface::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real sigma = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				sigma += neighboring_particle.W_ij_*base_particle_data_j.Vol_0_ / base_particle_data_i.Vol_0_;
			}

			fluid_data_i.rho_sum_ += sigma * fluid_data_i.rho_0_ / fluid_data_i.sigma_0_;
			fluid_data_i.rho_n_ = SMAX(fluid_data_i.temp_real_, fluid_data_i.rho_sum_);
		}
		//===========================================================//
		void DivergenceCorrection::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real div_correction = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				div_correction += neighboring_particle.dW_ij_ * neighboring_particle.r_ij_ *
					dot(-neighboring_particle.e_ij_, neighboring_particle.e_ij_) * base_particle_data_j.Vol_;
			}

			fluid_data_i.div_correction_ = div_correction;
		}
		//===========================================================//
		void DivergenceCorrection::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real div_correction = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				div_correction += neighboring_particle.dW_ij_ * neighboring_particle.r_ij_ *
					dot(-neighboring_particle.e_ij_, neighboring_particle.e_ij_)*base_particle_data_j.Vol_;
			}

			fluid_data_i.div_correction_ += div_correction;
		}
		//===========================================================//
		void  DivergenceCorrection::Update(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			Real div_correction_1 = fluid_data_i.div_correction_ / dimension_;
			fluid_data_i.div_correction_
				= 1.0 / (div_correction_1 + 0.1*(1.0 - div_correction_1)*(1.0 - div_correction_1));
		}
		//===========================================================//
		void ComputingViscousAcceleration::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				acceleration += 2.0*mu_*vel_detivative*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					*neighboring_particle.dW_ij_;
			}

			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//===========================================================//
		void ComputingViscousAcceleration
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				//viscous force with a simple wall model for high-Reynolds number flow
				Vecd vel_detivative = 2.0*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_)
					/ (neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				Real vel_difference = 0.03*(base_particle_data_i.vel_n_ - solid_data_j.vel_ave_).norm()
					* neighboring_particle.r_ij_;
				acceleration += 2.0*SMAX(mu_, fluid_data_i.rho_n_*vel_difference)
					*vel_detivative*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					*neighboring_particle.dW_ij_;
			}

			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//===========================================================//
		void TransportVelocityStress::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

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
					* neighboring_particle.dW_ij_ * dot(fluid_data_i.dvel_dt_trans_, -neighboring_particle.e_ij_))
					+ (fluid_data_j.rho_n_*base_particle_data_j.vel_n_
					* neighboring_particle.dW_ij_ * dot(fluid_data_j.dvel_dt_trans_, -neighboring_particle.e_ij_)))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//===========================================================//
		void TransportVelocityStress
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				//exra stress
				acceleration += 0.5*dt*(fluid_data_i.rho_n_*base_particle_data_i.vel_n_
					* neighboring_particle.dW_ij_ * dot(fluid_data_i.dvel_dt_trans_, -neighboring_particle.e_ij_))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			base_particle_data_i.dvel_dt_others_ += acceleration;
		}
		//===========================================================//
		void TransportVelocityCorrection::SetupDynamics(Real dt)
		{
			Real speed_max = body_->speed_max_;
			Real density = material_->GetReferenceDensity();
			p_background_ =  10.0 * density * speed_max * speed_max;
		}
		//===========================================================//
		void TransportVelocityCorrection::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration_trans(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}
			fluid_data_i.dvel_dt_trans_ = acceleration_trans;
			base_particle_data_i.pos_n_ += acceleration_trans * dt*dt*0.5;
		}
		//===========================================================//
		void TransportVelocityCorrection::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration_trans(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				//acceleration for transport velocity
				acceleration_trans -= 2.0*p_background_*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}

			fluid_data_i.dvel_dt_trans_ += acceleration_trans;
			base_particle_data_i.pos_n_ += acceleration_trans * dt * dt * 0.5;
		}
		//===========================================================//
		TotalMechanicalEnergy::TotalMechanicalEnergy(FluidBody* body, ExternalForce *external_force)
			: WeaklyCompressibleFluidDynamicsSum<Real>(body)
		{
			initial_reference_ = 0.0;
			average_farctor_ = 1.0;// / Real(body_->number_of_real_particles_);
			potential_ = external_force->InducedAcceleration().norm();
		}
		//===========================================================//
		Real TotalMechanicalEnergy::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			return average_farctor_ * (
				0.5 * fluid_data_i.mass_* base_particle_data_i.vel_n_.normSqr()
				+ fluid_data_i.mass_*potential_*base_particle_data_i.pos_n_[1]);
		}
		//===========================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(FluidBody* body)
			: WeaklyCompressibleFluidDynamicsMaximum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time setep size due to linear viscosity
			Real rho_0 = material_->rho_0_;
			Real mu = material_->mu_;
			initial_reference_ = mu / rho_0 / smoothing_length_;
		}
		//===========================================================//
		Real GetAcousticTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i	= particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return SMAX(smoothing_length_ / sqrt(smoothing_length_
				/ (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				material_->GetSoundSpeed(fluid_data_i.p_, fluid_data_i.rho_n_)
				+ base_particle_data_i.vel_n_.norm());
		}
		//===========================================================//
		Real GetAcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			body_->signal_speed_max_ = reduced_value;
			return 0.6 * smoothing_length_ / (reduced_value + 1.0e-15);
		}
		//===========================================================//
		GetAdvectionTimeStepSize::GetAdvectionTimeStepSize(FluidBody* body, Real U_f)
			: GetAcousticTimeStepSize(body)
		{
			initial_reference_ = SMAX(initial_reference_, U_f);
		}
		//===========================================================//
		Real GetAdvectionTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			return SMAX(smoothing_length_ / sqrt(smoothing_length_
				/ (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				base_particle_data_i.vel_n_.norm() + 1.0e-15);
		}
		//===========================================================//
		Real GetAdvectionTimeStepSize::OutputResult(Real reduced_value)
		{
			body_->speed_max_ = reduced_value;
			return 0.25 * smoothing_length_ / (reduced_value + 1.0e-15);
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_* dt * 0.5;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
			fluid_data_i.p_ = material_->GetPressure(fluid_data_i.rho_n_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];

			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//low dissipation Riemann problem
				Real ul = dot(neighboring_particle.e_ij_, base_particle_data_i.vel_n_);
				Real ur = dot(neighboring_particle.e_ij_, base_particle_data_j.vel_n_);
				Real p_star = material_->RiemannSolverForPressure(fluid_data_i.rho_n_, fluid_data_j.rho_n_,
					fluid_data_i.p_, fluid_data_j.p_, ul, ur);
				acceleration -= 2.0*p_star*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}

			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Real face_wall_external_acceleration 
					= dot((external_force_->InducedAcceleration(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_, dt)
					- solid_data_j.dvel_dt_ave_), -solid_data_j.n_);
				Real p_in_wall = fluid_data_i.p_ + fluid_data_i.rho_n_* neighboring_particle.r_ij_ *
				 		dot(neighboring_particle.e_ij_, -solid_data_j.n_)
						*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

				Real p_star = (fluid_data_i.p_*rho_in_wall + p_in_wall * fluid_data_i.rho_n_)
					/ (fluid_data_i.rho_n_ + rho_in_wall);

				//pressure force
				acceleration -= 2.0*p_star*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}

			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface::Intermediate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface::InnerInteraction2nd(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real density_change_rate = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//low dissipation Riemann problem
				Real ul = dot(neighboring_particle.e_ij_, base_particle_data_i.vel_n_);
				Real ur = dot(neighboring_particle.e_ij_, base_particle_data_j.vel_n_);
				Real u_star = material_->RiemannSolverForVelocity(fluid_data_i.rho_n_, fluid_data_j.rho_n_,
					fluid_data_i.p_, fluid_data_j.p_, ul, ur);

				//simplify the equation to cancle vector operations
				density_change_rate += 2.0*fluid_data_i.rho_n_*base_particle_data_j.Vol_
					*(u_star - ul)*neighboring_particle.dW_ij_;
			}
			fluid_data_i.drho_dt_ = fluid_data_i.div_correction_*density_change_rate;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface
			::ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Real density_change_rate = 0.0;
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				//seems not using rimann problem solution is better
				Vecd vel_in_wall = 2.0*solid_data_j.vel_ave_ - base_particle_data_i.vel_n_;
				density_change_rate += fluid_data_i.rho_n_*base_particle_data_j.Vol_
						* dot(base_particle_data_i.vel_n_ - vel_in_wall, -neighboring_particle.e_ij_) 
						* neighboring_particle.dW_ij_;
			}

			fluid_data_i.drho_dt_ += fluid_data_i.div_correction_*density_change_rate;
		}
		//===========================================================//
		void PressureRelaxationVerletFreeSurface::Update(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_ * dt * 0.5;
		}
		//===========================================================//
		void PressureRelaxationVerlet::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				Real p_star = (fluid_data_i.p_*fluid_data_j.rho_n_ + fluid_data_j.p_*fluid_data_i.rho_n_)
						/ (fluid_data_i.rho_n_ + fluid_data_j.rho_n_);

				//pressure force
				acceleration -= 2.0*p_star*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}

			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void PressureRelaxationVerlet
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Real face_wall_external_acceleration 
					= dot((external_force_->InducedAcceleration(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_, dt)
					- solid_data_j.dvel_dt_ave_), -solid_data_j.n_);
				Real p_in_wall = fluid_data_i.p_ + fluid_data_i.rho_n_ * neighboring_particle.r_ij_ * 
					dot(neighboring_particle.e_ij_, -solid_data_j.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

				Real p_star = (fluid_data_i.p_*rho_in_wall + p_in_wall *fluid_data_i.rho_n_)
					/ (fluid_data_i.rho_n_ + rho_in_wall);

				//pressure force
				acceleration -= 2.0*p_star*base_particle_data_j.Vol_ / fluid_data_i.rho_n_
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_);
			}

			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void Oldroyd_B_FluidInitialCondition::Update(size_t index_particle_i, Real dt)
		{

			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_(0);
			base_particle_data_i.dvel_dt_(0);

			fluid_data_i.p_ = 0.0;
			fluid_data_i.vel_trans_(0);
			fluid_data_i.rho_0_
				= material_->ReinitializeRho(fluid_data_i.p_);
			fluid_data_i.rho_n_ = fluid_data_i.rho_0_;
			fluid_data_i.mass_
				= fluid_data_i.rho_0_*base_particle_data_i.Vol_;

			non_newtonian_fluid_data_i.tau_(0);
		}
		//===========================================================//
		VerletOldroyd_B_Fluid::VerletOldroyd_B_Fluid(FluidBody *body,
			StdVec<SolidBody*> interacting_bodies, ExternalForce *external_force)
			:Oldroyd_B_FluidDynamicsComplex2Levels(body, interacting_bodies),
			external_force_(external_force)
		{
			mu_ = material_->mu_;
			mu_p_ = material_->mu_p_;
			lambda_ = material_->lambda_;
			Real sound_speed = material_->GetReferenceSoundSpeed();
			Real density = material_->GetReferenceDensity();
			p_background_ = 0.25 * density * sound_speed * sound_speed;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_* dt * 0.5;
			non_newtonian_fluid_data_i.tau_ += non_newtonian_fluid_data_i.dtau_dt_ * dt * 0.5;
			base_particle_data_i.Vol_ = fluid_data_i.mass_ / fluid_data_i.rho_n_;
			fluid_data_i.p_ = material_->GetPressure(fluid_data_i.rho_n_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];
				ViscoelasticFluidParticleData &non_newtonian_fluid_data_j = particles_->viscoelastic_particle_data_[index_particle_j];

				//simple average for pressure
				Real p_star = (fluid_data_i.p_*fluid_data_j.rho_n_ + fluid_data_j.p_*fluid_data_i.rho_n_)
					/ (fluid_data_i.rho_n_ + fluid_data_j.rho_n_);

				//pressure and elastic force
				acceleration -= (2.0 * p_star * neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_)
					+ (non_newtonian_fluid_data_i.tau_ + non_newtonian_fluid_data_j.tau_)
					* neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Real face_wall_external_acceleration
					= dot((external_force_->InducedAcceleration(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_, dt)
						- solid_data_j.dvel_dt_ave_), -solid_data_j.n_);
				Real p_in_wall = fluid_data_i.p_ + fluid_data_i.rho_n_ * neighboring_particle.r_ij_ * 
					dot(neighboring_particle.e_ij_, -solid_data_j.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

				Real p_star = (fluid_data_i.p_*rho_in_wall + p_in_wall * fluid_data_i.rho_n_)
					/ (fluid_data_i.rho_n_ + rho_in_wall);

				//stress boundary condition 
				acceleration -= (2.0 * p_star * neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_)
					+ 2.0 * non_newtonian_fluid_data_i.tau_ * neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_))
					*base_particle_data_j.Vol_ / fluid_data_i.rho_n_;
			}

			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid::Intermediate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid::InnerInteraction2nd(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			Real density_change_rate = 0.0;
			Matd stress_rate(0);

			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				//low dissipation Riemann problem
				Real ul = dot(neighboring_particle.e_ij_, fluid_data_i.vel_trans_);
				Real ur = dot(neighboring_particle.e_ij_, fluid_data_j.vel_trans_);
				Real u_star = material_->RiemannSolverForVelocity(fluid_data_i.rho_n_, fluid_data_j.rho_n_,
					fluid_data_i.p_, fluid_data_j.p_, ul, ur);;
				density_change_rate += 2.0*fluid_data_i.rho_n_*base_particle_data_j.Vol_
					*(u_star - ul)*neighboring_particle.dW_ij_;
				//transport velocity
				Matd velocity_gradient = SimTK::outer((fluid_data_i.vel_trans_ - fluid_data_j.vel_trans_),
					neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_))*base_particle_data_j.Vol_;
				stress_rate -= ~velocity_gradient*non_newtonian_fluid_data_i.tau_ 
					+ non_newtonian_fluid_data_i.tau_*velocity_gradient + non_newtonian_fluid_data_i.tau_/lambda_
					+ (~velocity_gradient + velocity_gradient)*mu_p_/lambda_;
			}
			fluid_data_i.drho_dt_ = density_change_rate;
			non_newtonian_fluid_data_i.dtau_dt_ = stress_rate;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid
			::ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			Real density_change_rate = 0.0;
			Matd stress_rate(0);

			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				density_change_rate += 2.0*fluid_data_i.rho_n_
					*dot(fluid_data_i.vel_trans_ - solid_data_j.vel_ave_, 
							neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_)) * base_particle_data_j.Vol_;
				Matd velocity_gradient = SimTK::outer((fluid_data_i.vel_trans_ - solid_data_j.vel_ave_),
					neighboring_particle.dW_ij_ * (-neighboring_particle.e_ij_)) * base_particle_data_j.Vol_*2.0;
				stress_rate -= ~velocity_gradient*non_newtonian_fluid_data_i.tau_
					+ non_newtonian_fluid_data_i.tau_*velocity_gradient + non_newtonian_fluid_data_i.tau_ / lambda_
					+ (~velocity_gradient + velocity_gradient)*mu_p_ / lambda_;
			}

			fluid_data_i.drho_dt_ += density_change_rate;
			non_newtonian_fluid_data_i.dtau_dt_ += stress_rate;
		}
		//===========================================================//
		void VerletOldroyd_B_Fluid::Update(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			ViscoelasticFluidParticleData &non_newtonian_fluid_data_i = particles_->viscoelastic_particle_data_[index_particle_i];

			fluid_data_i.rho_n_ += fluid_data_i.drho_dt_ * dt * 0.5;
			non_newtonian_fluid_data_i.tau_ += non_newtonian_fluid_data_i.dtau_dt_ * dt * 0.5;
		}
		//===========================================================//
		void InflowBoundaryCondition
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.vel_n_ = base_particle_data_i.vel_n_*(1.0 - constrain_strength_)
				+ constrain_strength_*GetInflowVelocity(base_particle_data_i.pos_n_, base_particle_data_i.vel_n_);
		}
	}
}
