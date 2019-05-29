#include "solid_dynamics.h"
#include "solid_body.h"
#include "solid_body_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "mesh_cell_linked_list.h"
#include "weakly_compressible_fluid_particles.h"
#include "weakly_compressible_fluid.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//===========================================================//
		void NormalDirectionSummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd gradient(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				gradient += neighors[n].gradW_ij_;
			}
			solid_data_i.temp_vec_ = gradient;
			solid_data_i.n_0_ = -solid_data_i.temp_vec_ / (solid_data_i.temp_vec_.norm() + 1.0e-2);;
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//===========================================================//
		void NormalDirectionSummation::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd gradient(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				gradient += neighors[n].gradW_ij_;
			}
			solid_data_i.temp_vec_ += gradient;
			solid_data_i.n_0_ = -solid_data_i.temp_vec_ / (solid_data_i.temp_vec_.norm() + 1.0e-2);;
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//===========================================================//
		void NormalDirectionReNormalization::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				local_configuration += base_particle_data_j.Vol_
					*SimTK::outer(neighboring_particle.r_ij_, neighboring_particle.gradW_ij_);
				gradient += neighboring_particle.gradW_ij_ * base_particle_data_j.Vol_;
			}

			solid_data_i.temp_vec_ = gradient;
			elastic_data_i.temp_matrix_ = local_configuration;
		}
		//===========================================================//
		void NormalDirectionReNormalization::ContactInteraction(size_t index_particle_i,
			size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				local_configuration += base_particle_data_j.Vol_
					*SimTK::outer(neighboring_particle.r_ij_, neighboring_particle.gradW_ij_);
				gradient += neighboring_particle.gradW_ij_ * base_particle_data_j.Vol_;
			}

			solid_data_i.temp_vec_ += gradient;
			elastic_data_i.temp_matrix_ += local_configuration;
		}
		//===========================================================//
		void  NormalDirectionReNormalization::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd correction_matrix = GeneralizedInverse(elastic_data_i.temp_matrix_);

			Vecd n_temp_ = ~correction_matrix * solid_data_i.temp_vec_;
			if (n_temp_.norm() <= 0.75)
			{
				solid_data_i.n_0_ = Vecd(0.0);
			}else{
				solid_data_i.n_0_ =  - n_temp_ / (n_temp_.norm() + 1.0e-2);
			}
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//===========================================================//
		void InitializeDisplacement::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i
				= particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.pos_temp_ = base_particle_data_i.pos_n_;
		}
		//===========================================================//
		void UpdateAverageVelocity::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			solid_data_i.vel_ave_ = (base_particle_data_i.pos_n_ - elastic_data_i.pos_temp_) / (dt + 1.0e-15);
		}
		//===========================================================//
		void FluidViscousForceOnSolid::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.viscous_force_from_fluid_ = Vecd(0);
		}
		//===========================================================//
		void FluidViscousForceOnSolid
			::ContactInteraction(size_t index_particle_i,
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				WeaklyCompressibleFluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_data_[index_particle_j];

				//froce due to viscousity
				//viscous force with a simple wall model for high-Reynolds number flow
				Vecd vel_detivative = 2.0*(solid_data_i.vel_ave_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle.r_ij_.norm() + 0.01 * smoothing_length_);
				Real vel_difference = 0.03*(solid_data_i.vel_ave_ - base_particle_data_j.vel_n_).norm()
					*neighboring_particle.r_ij_.norm();

				force += 2.0*SMAX(mu_, fluid_data_j.rho_n_*vel_difference)
					*vel_detivative*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					*neighboring_particle.dW_ij_;

			}

			solid_data_i.viscous_force_from_fluid_ += force;
		}
		//===========================================================//
		TotalViscousForceOnSolid
			::TotalViscousForceOnSolid(SolidBody *body)
			: ParticleDynamicsSum<Vecd, SolidBody, SolidBodyParticles>(body)
		{
			initial_reference_(0);
		}
	//===========================================================//
		Vecd TotalViscousForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.viscous_force_from_fluid_;
		}
		//===========================================================//
		Vecd TotalForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.force_from_fluid_;
		}
		//===========================================================//
		void FluidPressureForceOnSolid::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.force_from_fluid_ = solid_data_i.viscous_force_from_fluid_;
		}
		//===========================================================//
		void FluidPressureForceOnSolid
			::ContactInteraction(size_t index_particle_i, 
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors 
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				WeaklyCompressibleFluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_data_[index_particle_j];

				Real face_wall_external_acceleration
					= dot((external_force_->InducedAcceleration(base_particle_data_j.pos_n_, base_particle_data_j.vel_n_, dt)
						- solid_data_i.dvel_dt_ave_), -solid_data_i.n_);
				Real p_in_wall = fluid_data_j.p_ + fluid_data_j.rho_n_*dot(neighboring_particle.r_ij_, -solid_data_i.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);

				Real p_star = (fluid_data_j.p_*rho_in_wall + p_in_wall * fluid_data_j.rho_n_)
					/ (fluid_data_j.rho_n_ + rho_in_wall);

				//force due to pressure
				force -= 2.0*p_star*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					*neighboring_particle.gradW_ij_;

			}

			solid_data_i.force_from_fluid_ += force;
		}
		//===========================================================//
		void FluidPressureForceOnSolidFreeSurface
			::ContactInteraction(size_t index_particle_i,
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				WeaklyCompressibleFluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_data_[index_particle_j];

				Real face_wall_external_acceleration
					= dot((external_force_->InducedAcceleration(base_particle_data_j.pos_n_, base_particle_data_j.vel_n_, dt)
						- solid_data_i.dvel_dt_ave_), -solid_data_i.n_);
				Real p_in_wall = fluid_data_j.p_ + fluid_data_j.rho_n_*dot(neighboring_particle.r_ij_, -solid_data_i.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = material_->ReinitializeRho(p_in_wall);
				Vecd vel_in_wall = 2.0*solid_data_i.vel_ave_ - base_particle_data_j.vel_n_;

				//low dissipation Riemann problem
				Real ul = dot(-solid_data_i.n_, base_particle_data_j.vel_n_);
				Real ur = dot(-solid_data_i.n_, vel_in_wall);
				Real p_star = material_->RiemannSolverForPressure(fluid_data_j.rho_n_, rho_in_wall,
					fluid_data_j.p_, p_in_wall, ul, ur);

				//force due to pressure
				force -= 2.0*p_star*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					*neighboring_particle.gradW_ij_;

			}

			solid_data_i.force_from_fluid_ += force;
		}
		//===========================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(ElasticBody* body)
			: ParticleDynamicsMinimum<ElasticBody, ElasticBodyParticles>(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			ElasticSolid *material = body_->material_;
			Real rho_0 = material->rho_0_;
			Real eta_0 = material->eta_0_;
			//time setep size due to linear viscosity
			initial_reference_ = 0.5 * smoothing_length_*smoothing_length_ * rho_0 / (eta_0 + 1.0e-15);
		}
		Real GetAcousticTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			Real sound_speed = elastic_data_i.local_c_;
			return 0.6 * SMIN(sqrt(smoothing_length_ / (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				smoothing_length_ / (sound_speed + base_particle_data_i.vel_n_.norm()));
		}
		//===========================================================//
		void CorrectConfiguration::InnerInteraction(size_t index_particle_i, Real dt)
		{
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				local_configuration += base_particle_data_j.Vol_
					*SimTK::outer(neighboring_particle.r_ij_, neighboring_particle.gradW_ij_);
			}
			elastic_data_i.temp_matrix_ = local_configuration;
		}
		//===========================================================//
		void CorrectConfiguration::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				local_configuration += base_particle_data_j.Vol_
					*SimTK::outer(neighboring_particle.r_ij_, neighboring_particle.gradW_ij_);
			}
			elastic_data_i.temp_matrix_ += local_configuration;
		}
		//===========================================================//
		void CorrectConfiguration::Update(size_t index_particle_i, Real dt)
		{
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.B_ = GeneralizedInverse(elastic_data_i.temp_matrix_);
		}
		//===========================================================//
		void DeformationGradientTensorBySummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//unint matrix
			Matd deformation(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				deformation -= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.pos_n_ - base_particle_data_j.pos_n_), 
						neighboring_particle.gradW_ij_);
			}
			elastic_data_i.F_ = deformation * elastic_data_i.B_;
		}
		//===========================================================//
		void DeformationGradientTensorBySummation
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				deformation -= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.pos_n_ - base_particle_data_j.pos_n_),
						neighboring_particle.gradW_ij_);
			}
			elastic_data_i.F_ += deformation * elastic_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_ 
				= body_->GetElasticStress(elastic_data_i.F_, index_particle_i)
				+ material_->DampingStress(elastic_data_i.F_, 
					elastic_data_i.dF_dt_, elastic_data_i.local_eta_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_*dt * 0.5;
		}
		//===========================================================//
		void StressRelaxation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_
				+ solid_data_i.force_from_fluid_ / elastic_data_i.mass_;

			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_*elastic_data_i.B_
					+ elastic_data_j.stress_ *elastic_data_j.B_)
					*neighboring_particle.gradW_ij_
					*base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void StressRelaxation::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ *elastic_data_i.B_
					+ elastic_data_j.stress_*elastic_data_j.B_)
					*neighboring_particle.gradW_ij_
					*base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void StressRelaxation::Intermediate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ 
				+= base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;

		}
		//===========================================================//
		void StressRelaxation::InnerInteraction2nd(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_),
						neighboring_particle.gradW_ij_);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate * elastic_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation
			::ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - solid_data_j.vel_ave_),
						neighboring_particle.gradW_ij_);
			}
			elastic_data_i.dF_dt_ += deformation_gradient_change_rate * elastic_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation::Update(size_t index_particle_i, Real dt)
		{
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void StressRelaxationFirstStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_ = body_->GetElasticStress(elastic_data_i.F_, index_particle_i)
				+ material_->DampingStress(elastic_data_i.F_, elastic_data_i.dF_dt_, elastic_data_i.local_eta_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_*dt * 0.5;

		}
		//===========================================================//
		void StressRelaxationFirstStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//including gravity and force from fluid
			Vecd acceleration = base_particle_data_i.dvel_dt_others_ 
				+ solid_data_i.force_from_fluid_/ elastic_data_i.mass_;

			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ *elastic_data_i.B_
					+ elastic_data_j.stress_*elastic_data_j.B_)
					*neighboring_particle.gradW_ij_
					*base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void StressRelaxationFirstStep::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
		}
		//===========================================================//
		void StressRelaxationSecondStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i 	= particles_->base_particle_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

		base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;

		}
		//===========================================================//
		void StressRelaxationSecondStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				//deformtion
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_),
						neighboring_particle.gradW_ij_);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate* elastic_data_i.B_;
		}
		//===========================================================//
		void StressRelaxationSecondStep::Update(size_t index_particle_i, Real dt)
		{
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void StressInConstrinedElasticBodyFirstHalf::ParticleUpdate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_
				= body_->GetElasticStress(elastic_data_i.F_, index_particle_i)
				+ material_->DampingStress(elastic_data_i.F_, elastic_data_i.dF_dt_, elastic_data_i.local_eta_);
		}
		//===========================================================//
		void StressInConstrinedElasticBodySecondHalf
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidBodyParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((solid_data_i.vel_ave_ - base_particle_data_j.vel_n_), neighboring_particle.gradW_ij_);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate * elastic_data_i.B_;
			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void ConstrainSolidBodyRegion::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i
				= particles_->solid_body_data_[index_particle_i];

			Vecd pos_old = base_particle_data_i.pos_n_;
			base_particle_data_i.pos_n_ = solid_data_i.pos_0_ + GetDisplacement(pos_old);
			base_particle_data_i.vel_n_ = GetVelocity(pos_old);
			base_particle_data_i.dvel_dt_ = GetAcceleration(pos_old);
			/** the average values are prescirbed also. */
			solid_data_i.vel_ave_ = base_particle_data_i.vel_n_;
			solid_data_i.dvel_dt_ave_ = base_particle_data_i.dvel_dt_;
		}
		//===========================================================//
		void ImposeExternalForce
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd induced_acceleration = GetAcceleration(base_particle_data_i.pos_n_);
			base_particle_data_i.vel_n_ += induced_acceleration * dt;
			solid_data_i.vel_ave_ = base_particle_data_i.vel_n_;
		}
		//=================================================================================================//
		ConstrianSoildBodyPartBySimBody
			::ConstrianSoildBodyPartBySimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ)
			: LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>(body, body_part),
			MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
		void  ConstrianSoildBodyPartBySimBody
			::PrepareConstraint()
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
		}
		//=================================================================================================//
		ForceOnSolidBodyPartForSimBody
			::ForceOnSolidBodyPartForSimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ)
			: LagrangianConstraintSum<SpatialVec, SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>(body, body_part),
			MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			initial_reference_ = SpatialVec(Vec3(0), Vec3(0));
		}
		//=================================================================================================//
		void ForceOnSolidBodyPartForSimBody::SetupReduce()
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
		ForceOnElasticBodyPartForSimBody
			::ForceOnElasticBodyPartForSimBody(ElasticBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ)
			: LagrangianConstraintSum<SpatialVec, ElasticBody, ElasticBodyParticles, SolidBodyPartForSimbody>(body, body_part),
			MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			initial_reference_ = SpatialVec(Vec3(0), Vec3(0));
		}
		//=================================================================================================//
		void ForceOnElasticBodyPartForSimBody::SetupReduce()
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
	//===========================================================//
		DampingBySplittingAlgorithm
			::DampingBySplittingAlgorithm(ElasticBody *elastic_body, Real eta)
			: ParticleDynamicsInnerSplitting<ElasticBody, ElasticBodyParticles>(elastic_body)
		{
			eta_ = eta;
		}
		//===========================================================//
		void DampingBySplittingAlgorithm
			::InnerInteraction(size_t index_particle_i, Real dt)
		{

			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidBodyParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticBodyParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= particles_->base_particle_data_[index_particle_j];
				ElasticBodyParticleData &elastic_data_j
					= particles_->elastic_body_data_[index_particle_j];

				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ neighboring_particle.r_ij_.norm();
				acceleration += 0.5*eta_*vel_detivative
					*neighboring_particle.dW_ij_*elastic_data_j.mass_
					/ elastic_data_i.rho_0_ / elastic_data_j.rho_0_;
			}
			
			base_particle_data_i.dvel_dt_ += acceleration * 0.5;
			base_particle_data_i.vel_n_ += acceleration * dt* 0.5;
		}
		//===========================================================//
	}
}
