#include "solid_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "mesh_cell_linked_list.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//===========================================================//
		void SolidDynamicsInitialCondition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_body_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
			base_particle_data_i.dvel_dt_ = zero;

			solid_body_data_i.vel_ave_ = zero;
			solid_body_data_i.dvel_dt_ave_ = zero;
		}
		//===========================================================//
		void OffsetInitialParticlePosition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_particle_data_i = particles_->solid_body_data_[index_particle_i];

			base_particle_data_i.pos_n_ += offset_;
			solid_particle_data_i.pos_0_ += offset_;
		}
		//===========================================================//
		void TwistingInitialCondition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_body_data_i = particles_->solid_body_data_[index_particle_i];

			Real x = solid_body_data_i.pos_0_[0];
			Real y = solid_body_data_i.pos_0_[1];
			Real z = solid_body_data_i.pos_0_[2];
			Real angular_velocity = dt * sin((pi * x) / (2.0 * 6.0));
			Real local_radius = sqrt(pow(y, 2.0) + pow(z, 2.0));
			Real angular = atan2(y, z);
			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;

			if (x > 0.0) {
				base_particle_data_i.vel_n_[1] = angular_velocity * local_radius * cos(angular);
				base_particle_data_i.vel_n_[2] = -angular_velocity * local_radius * sin(angular);
			}
		}
		//===========================================================//
		void NormalDirectionSummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd gradient(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				gradient += neighors[n].dW_ij_ * neighors[n].e_ij_;
			}
			solid_data_i.temp_vec_ = gradient;
			solid_data_i.n_0_ = -solid_data_i.temp_vec_ / (solid_data_i.temp_vec_.norm() + 1.0e-2);;
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//===========================================================//
		void NormalDirectionSummation::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd gradient(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				gradient += neighors[n].dW_ij_ * neighors[n].e_ij_;
			}
			solid_data_i.temp_vec_ += gradient;
			solid_data_i.n_0_ = -solid_data_i.temp_vec_ / (solid_data_i.temp_vec_.norm() + 1.0e-2);;
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//===========================================================//
		void NormalDirectionReNormalization::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				Vecd r_ij = -neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
				gradient += gradw_ij * base_particle_data_j.Vol_;
			}

			solid_data_i.temp_vec_ = gradient;
			solid_data_i.temp_matrix_ = local_configuration;
		}
		//===========================================================//
		void NormalDirectionReNormalization::ContactInteraction(size_t index_particle_i,
			size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				Vecd r_ij = -neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
				gradient += gradw_ij * base_particle_data_j.Vol_;
			}

			solid_data_i.temp_vec_ += gradient;
			solid_data_i.temp_matrix_ += local_configuration;
		}
		//===========================================================//
		void  NormalDirectionReNormalization::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd correction_matrix = inverse(solid_data_i.temp_matrix_);

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
		void ElasticSolidDynamicsInitialCondition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_body_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_particle_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
			base_particle_data_i.dvel_dt_ = zero;

			solid_body_data_i.vel_ave_ = zero;
			solid_body_data_i.dvel_dt_ave_ = zero;

			elastic_particle_data_i.rho_0_ = material_->rho_0_;
			elastic_particle_data_i.rho_n_ = material_->rho_0_;
			elastic_particle_data_i.mass_ = material_->rho_0_*base_particle_data_i.Vol_;
		}
		//===========================================================//
		void InitializeDisplacement::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i
				= particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.pos_temp_ = base_particle_data_i.pos_n_;
		}
		//===========================================================//
		void UpdateAverageVelocity::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			solid_data_i.vel_ave_ = (base_particle_data_i.pos_n_ - elastic_data_i.pos_temp_) / (dt + 1.0e-15);
		}
		//===========================================================//
		void FluidViscousForceOnSolid::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.viscous_force_from_fluid_ = Vecd(0);
		}
		//===========================================================//
		void FluidViscousForceOnSolid
			::ContactInteraction(size_t index_particle_i,
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_particle_data_[index_particle_j];

				//froce due to viscousity
				//viscous force with a simple wall model for high-Reynolds number flow
				Vecd vel_detivative = 2.0*(solid_data_i.vel_ave_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle.r_ij_ + 0.01 * smoothing_length_);
				Real vel_difference = 0.03*(solid_data_i.vel_ave_ - base_particle_data_j.vel_n_).norm()
					*neighboring_particle.r_ij_ ;

				force += 2.0*SMAX(mu_, fluid_data_j.rho_n_*vel_difference)
					*vel_detivative*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					*neighboring_particle.dW_ij_;

			}

			solid_data_i.viscous_force_from_fluid_ += force;
		}
		//===========================================================//
		TotalViscousForceOnSolid
			::TotalViscousForceOnSolid(SolidBody *body) : SolidDynamicsSum<Vecd>(body)
		{
			initial_reference_(0);
		}
	//===========================================================//
		Vecd TotalViscousForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.viscous_force_from_fluid_;
		}
		//===========================================================//
		Vecd TotalForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.force_from_fluid_;
		}
		//===========================================================//
		void FluidPressureForceOnSolid::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.force_from_fluid_ = solid_data_i.viscous_force_from_fluid_;
		}
		//===========================================================//
		void FluidPressureForceOnSolid
			::ContactInteraction(size_t index_particle_i, 
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors 
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_particle_data_[index_particle_j];

				Real face_wall_external_acceleration
					= dot((external_force_->InducedAcceleration(base_particle_data_j.pos_n_, base_particle_data_j.vel_n_, dt)
						- solid_data_i.dvel_dt_ave_), -solid_data_i.n_);
				Real p_in_wall = fluid_data_j.p_ + fluid_data_j.rho_n_ * neighboring_particle.r_ij_ *
					dot(neighboring_particle.e_ij_, -solid_data_i.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = interacting_material_[interacting_body_index]->ReinitializeRho(p_in_wall);

				Real p_star = (fluid_data_j.p_*rho_in_wall + p_in_wall * fluid_data_j.rho_n_)
					/ (fluid_data_j.rho_n_ + rho_in_wall);

				//force due to pressure
				force += 2.0*p_star*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;

			}

			solid_data_i.force_from_fluid_ += force;
		}
		//===========================================================//
		void FluidPressureForceOnSolidFreeSurface
			::ContactInteraction(size_t index_particle_i,
				size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd force(0);
			StdVec<NeighboringParticle>  &neighors
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				WeaklyCompressibleFluid *interacting_material = interacting_material_[interacting_body_index];
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index])
					.base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = (*interacting_particles_[interacting_body_index])
					.fluid_particle_data_[index_particle_j];

				Real face_wall_external_acceleration
					= dot((external_force_->InducedAcceleration(base_particle_data_j.pos_n_, base_particle_data_j.vel_n_, dt)
						- solid_data_i.dvel_dt_ave_), -solid_data_i.n_);
				Real p_in_wall = fluid_data_j.p_ + fluid_data_j.rho_n_ * neighboring_particle.r_ij_ *
					dot(neighboring_particle.e_ij_, -solid_data_i.n_)
					*SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = interacting_material->ReinitializeRho(p_in_wall);
				Vecd vel_in_wall = 2.0*solid_data_i.vel_ave_ - base_particle_data_j.vel_n_;

				//low dissipation Riemann problem
				Real ul = dot(-solid_data_i.n_, base_particle_data_j.vel_n_);
				Real ur = dot(-solid_data_i.n_, vel_in_wall);
				Real p_star = interacting_material->RiemannSolverForPressure(fluid_data_j.rho_n_, rho_in_wall,
					fluid_data_j.p_, p_in_wall, ul, ur);

				//force due to pressure
				force += 2.0*p_star*base_particle_data_i.Vol_ * base_particle_data_j.Vol_
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;

			}

			solid_data_i.force_from_fluid_ += force;
		}
		//===========================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(SolidBody* body)
			: ElasticSolidDynamicsMinimum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			Real rho_0 = material_->rho_0_;
			Real eta_0 = material_->eta_0_;
			//time setep size due to linear viscosity
			initial_reference_ = 0.5 * smoothing_length_*smoothing_length_ * rho_0 / (eta_0 + 1.0e-15);
		}
		//===========================================================//
		Real GetAcousticTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			Real sound_speed = material_->GetSoundSpeed(index_particle_i);
			return 0.6 * SMIN(sqrt(smoothing_length_ / (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				smoothing_length_ / (sound_speed + base_particle_data_i.vel_n_.norm()));
		}
		//===========================================================//
		void CorrectConfiguration::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				Vecd r_ij = -neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
			}
			solid_data_i.temp_matrix_ = local_configuration;
		}
		//===========================================================//
		void CorrectConfiguration::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				Vecd r_ij = -neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
			}
			solid_data_i.temp_matrix_ += local_configuration;
		}
		//===========================================================//
		void CorrectConfiguration::Update(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			/** note that the generalized inverse only works here*/
			solid_data_i.B_ = GeneralizedInverse(solid_data_i.temp_matrix_);
		}
		//===========================================================//
		void DeformationGradientTensorBySummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//unint matrix
			Matd deformation(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation -= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.pos_n_ - base_particle_data_j.pos_n_), gradw_ij);
			}
			elastic_data_i.F_ = deformation * solid_data_i.B_;
		}
		//===========================================================//
		void DeformationGradientTensorBySummation
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation -= base_particle_data_j.Vol_
					* SimTK::outer((base_particle_data_i.pos_n_ - base_particle_data_j.pos_n_), gradw_ij);
			}
			elastic_data_i.F_ += deformation * solid_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_ 
				= material_->ConstitutiveRelation(elastic_data_i.F_, index_particle_i)
				+ material_->NumericalDampingStress(elastic_data_i.F_, 
					elastic_data_i.dF_dt_, numerical_viscosity_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_*dt * 0.5;
		}
		//===========================================================//
		void StressRelaxation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_
				+ solid_data_i.force_from_fluid_ / elastic_data_i.mass_;

			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_*solid_data_i.B_
					+ elastic_data_j.stress_ *solid_data_j.B_)
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					* base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void StressRelaxation::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ *solid_data_i.B_
					+ elastic_data_j.stress_*solid_data_j.B_)
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					*base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void StressRelaxation::Intermediate(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ 
				+= base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;

		}
		//===========================================================//
		void StressRelaxation::InnerInteraction2nd(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_), gradw_ij);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate * solid_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation
			::ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - solid_data_j.vel_ave_), gradw_ij);
			}
			elastic_data_i.dF_dt_ += deformation_gradient_change_rate * solid_data_i.B_;
		}
		//===========================================================//
		void StressRelaxation::Update(size_t index_particle_i, Real dt)
		{
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void StressRelaxationFirstStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_ = material_->ConstitutiveRelation(elastic_data_i.F_, index_particle_i)
				+ material_->NumericalDampingStress(elastic_data_i.F_, elastic_data_i.dF_dt_, numerical_viscosity_);
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_*dt * 0.5;

		}
		//===========================================================//
		void StressRelaxationFirstStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

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
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ *solid_data_i.B_
					+ elastic_data_j.stress_*solid_data_j.B_)
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					* base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void StressRelaxationFirstStep::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
		}
		//===========================================================//
		void StressRelaxationSecondStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i 	= particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

		base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;

		}
		//===========================================================//
		void StressRelaxationSecondStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				//deformtion
				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_), gradw_ij);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate* solid_data_i.B_;
		}
		//===========================================================//
		void StressRelaxationSecondStep::Update(size_t index_particle_i, Real dt)
		{
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void StressInConstrinedElasticBodyFirstHalf::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_
				= material_->ConstitutiveRelation(elastic_data_i.F_, index_particle_i)
				+ material_->NumericalDampingStress(elastic_data_i.F_, elastic_data_i.dF_dt_, numerical_viscosity_);
		}
		//===========================================================//
		void StressInConstrinedElasticBodySecondHalf
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j
					= (*interacting_particles_[interacting_body_index]).solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j
					= (*interacting_particles_[interacting_body_index]).elastic_body_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((solid_data_i.vel_ave_ - base_particle_data_j.vel_n_), gradw_ij);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate * solid_data_i.B_;
			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//===========================================================//
		void ConstrainSolidBodyRegion::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
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
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

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
			: LagrangianConstraint<SolidBody, SolidParticles, SolidBodyPartForSimbody>(body, body_part),
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
			: SolidDynamicsConstraintForSimbodySum<SpatialVec>(body, body_part),
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
			::ForceOnElasticBodyPartForSimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ)
			: ElasticSolidDynamicsConstraintForSimbodySum<SpatialVec>(body, body_part),
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
		DampingBySplittingAlgorithm
			::DampingBySplittingAlgorithm(SolidBody *elastic_body, Real eta)
			: ElasticSolidDynamicsInnerSplitting(elastic_body)
		{
			eta_ = eta;
		}
		//===========================================================//
		void DampingBySplittingAlgorithm
			::InnerInteraction(size_t index_particle_i, Real dt)
		{

			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= particles_->base_particle_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j
					= particles_->elastic_body_data_[index_particle_j];

				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ neighboring_particle.r_ij_ ;
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
