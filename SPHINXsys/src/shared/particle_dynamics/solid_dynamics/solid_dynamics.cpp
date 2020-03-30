/**
 * @file 	solid_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_dynamics.h"

using namespace SimTK;
//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		void NormalDirectionSummation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd gradient(0.0);
			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				gradient += inner_neighors[n]->dW_ij_ * inner_neighors[n]->e_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					gradient += contact_neighors[n]->dW_ij_ * contact_neighors[n]->e_ij_;
				}
			}

			solid_data_i.n_0_ = - gradient / (gradient.norm() + 1.0e-2);;
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//=================================================================================================//
		void NormalDirectionReNormalization::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				Vecd r_ij = -neighboring_particle->r_ij_ * neighboring_particle->e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
				gradient += gradw_ij * base_particle_data_j.Vol_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
					Vecd r_ij = -neighboring_particle->r_ij_ * neighboring_particle->e_ij_;
					local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
					gradient += gradw_ij * base_particle_data_j.Vol_;
				}
			}

			Matd correction_matrix = inverse(local_configuration);
			Vecd n_temp_ = ~correction_matrix * gradient;
			if (n_temp_.norm() <= 0.75)	{
				solid_data_i.n_0_ = Vecd(0.0);
			}
			else {
				solid_data_i.n_0_ = -n_temp_ / (n_temp_.norm() + 1.0e-2);
			}
			solid_data_i.n_ = solid_data_i.n_0_;
		}
		//=================================================================================================//
		void InitializeDisplacement::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i
				= particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.pos_temp_ = base_particle_data_i.pos_n_;
		}
		//=================================================================================================//
		void UpdateAverageVelocity::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			solid_data_i.vel_ave_ = (base_particle_data_i.pos_n_ - elastic_data_i.pos_temp_) / (dt + 1.0e-15);
		}
		//=================================================================================================//
		void FluidViscousForceOnSolid::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.viscous_force_from_fluid_ = Vecd(0);

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					FluidParticleData& fluid_data_j = (*interacting_particles_[k])
						.fluid_particle_data_[index_particle_j];

					//froce due to viscousity
					//viscous force with a simple wall model for high-Reynolds number flow
					Vecd vel_detivative = 2.0 * (solid_data_i.vel_ave_ - base_particle_data_j.vel_n_)
						/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_);
					Real vel_difference = 0.0 * (solid_data_i.vel_ave_ - base_particle_data_j.vel_n_).norm()
						* neighboring_particle->r_ij_;

					force += 2.0 * SMAX(mu_, fluid_data_j.rho_n_ * vel_difference)
						* vel_detivative * base_particle_data_i.Vol_ * base_particle_data_j.Vol_
						* neighboring_particle->dW_ij_;
				}
			}

			solid_data_i.viscous_force_from_fluid_ += force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnSolid::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.viscous_force_from_fluid_ = Vecd(0);

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					FluidParticleData& fluid_data_j = (*interacting_particles_[k])
						.fluid_particle_data_[index_particle_j];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
					 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
					Real v_r_ij = dot(solid_data_i.vel_ave_ - base_particle_data_j.vel_n_,
						neighboring_particle->r_ij_ * neighboring_particle->e_ij_);
					Real vel_difference = 0.0 * (solid_data_i.vel_ave_ - base_particle_data_j.vel_n_).norm()
						* neighboring_particle->r_ij_;
					Real eta_ij = 8.0 * SMAX(mu_, fluid_data_j.rho_n_ * vel_difference) * v_r_ij /
						(neighboring_particle->r_ij_ * neighboring_particle->r_ij_ + 0.01 * smoothing_length_);
					force += eta_ij * base_particle_data_i.Vol_ * base_particle_data_j.Vol_
						* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				}
			}

			solid_data_i.viscous_force_from_fluid_ += force;
		}
		//=================================================================================================//
		TotalViscousForceOnSolid
			::TotalViscousForceOnSolid(SolidBody *body) : SolidDynamicsSum<Vecd>(body)
		{
			initial_reference_(0);
		}
		//=================================================================================================//
		Vecd TotalViscousForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.viscous_force_from_fluid_;
		}
		//=================================================================================================//
		Vecd TotalForceOnSolid::ReduceFunction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			return solid_data_i.force_from_fluid_;
		}
		//=================================================================================================//
		void FluidPressureForceOnSolid::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			solid_data_i.force_from_fluid_ = solid_data_i.viscous_force_from_fluid_;

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Real particle_spacing_i1 = 1.0 / body_->particle_spacing_;
				Vecd n_i = solid_data_i.n_;
				Real particle_spacing_ratio2 = 1.0 / (interacting_bodies_[k]->particle_spacing_
					* particle_spacing_i1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					FluidParticleData& fluid_data_j = (*interacting_particles_[k])
						.fluid_particle_data_[index_particle_j];

					Vecd e_ij = neighboring_particle->e_ij_;
					Real face_wall_external_acceleration
						= dot((base_particle_data_j.dvel_dt_others_ - solid_data_i.dvel_dt_ave_), e_ij);
					Real p_star = fluid_data_j.p_ + 0.5 * fluid_data_j.rho_n_
						* neighboring_particle->r_ij_ * SMAX(0.0, face_wall_external_acceleration);

					/** penalty correction to prevent particle running into boundary */
					Real projection = dot(-e_ij, n_i);
					Real delta = 2.0 * projection * neighboring_particle->r_ij_ * particle_spacing_i1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty = beta * projection * fabs(p_star);

					//force due to pressure
					force -= 2.0 * (p_star * e_ij - penalty * n_i)
						* base_particle_data_i.Vol_ * base_particle_data_j.Vol_ * neighboring_particle->dW_ij_;
				}
			}
			
			solid_data_i.force_from_fluid_ += force;
		}
		//=================================================================================================//
		GetAcousticTimeStepSize::GetAcousticTimeStepSize(SolidBody* body)
			: ElasticSolidDynamicsMinimum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time setep size due to linear viscosity
			initial_reference_ = material_->getViscousTimeStepSize(smoothing_length_);
		}
		//=================================================================================================//
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
		//=================================================================================================//
		void CorrectConfiguration::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				Vecd r_ij = -neighboring_particle->r_ij_ * neighboring_particle->e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
			}

			/** note that the generalized inverse only works here*/
			solid_data_i.B_ = GeneralizedInverse(local_configuration);
		}
		//=================================================================================================//
		void DeformationGradientTensorBySummation::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation(0.0);
			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				deformation -= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.pos_n_ - base_particle_data_j.pos_n_), gradw_ij);
			}

			elastic_data_i.F_ = deformation * solid_data_i.B_;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Initialization(size_t index_particle_i, Real dt)
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
		//=================================================================================================//
		void StressRelaxationFirstHalf::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			//including gravity and force from fluid
			Vecd acceleration = base_particle_data_i.dvel_dt_others_ 
				+ solid_data_i.force_from_fluid_/ elastic_data_i.mass_;

			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ *solid_data_i.B_
					+ elastic_data_j.stress_*solid_data_j.B_)
					* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_
					* base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i 	= particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_), gradw_ij);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate* solid_data_i.B_;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Update(size_t index_particle_i, Real dt)
		{
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
		//=================================================================================================//
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
		//=================================================================================================//
		void ConstrainSolidBodyRegionSinusoidalMotion::ConstraintAParticle(size_t index_particle_i, Real dt)
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
		//=================================================================================================//
		Vecd ConstrainSolidBodyRegionSinusoidalMotion::GetDisplacement(Vecd &pos)
		{	
			Vecd disp(0.0);
			disp[1] = h_m_ * sin(2.0 * pi * f_ * GlobalStaticVariables::physical_time_ + Real((id_ -1)) * phi_);
			return disp;
		}
		//=================================================================================================//
		Vecd ConstrainSolidBodyRegionSinusoidalMotion::GetVelocity(Vecd &pos)
		{
			Vecd disp(0.0);
			disp[1] = h_m_ * 2.0 * pi * f_ * cos(2.0 * pi * f_ * GlobalStaticVariables::physical_time_ + Real((id_ -1)) * phi_);
			return disp;
		}
		//=================================================================================================//
		Vecd ConstrainSolidBodyRegionSinusoidalMotion::GetAcceleration(Vecd &pos)
		{
			Vecd disp(0.0);
			disp[1] = -h_m_ * 2.0 * pi * f_ * 2.0 * pi * f_ * cos(2.0 * pi * f_ * GlobalStaticVariables::physical_time_ + Real((id_ -1)) * phi_);
			return disp;
		}
		//=================================================================================================//
		void constrainNormDirichletBoundary::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= particles_->solid_body_data_[index_particle_i];

			Vecd pos_old = base_particle_data_i.pos_n_;
			base_particle_data_i.pos_n_[axis_id_] = solid_data_i.pos_0_[axis_id_];
			base_particle_data_i.vel_n_[axis_id_] = 0.0;
			base_particle_data_i.dvel_dt_[axis_id_] = 0.0;
			/** the average values are prescirbed also. */
			solid_data_i.vel_ave_ = base_particle_data_i.vel_n_;
			solid_data_i.dvel_dt_ave_ = base_particle_data_i.dvel_dt_;
		}
		//=================================================================================================//
		void ImposeExternalForce
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Vecd induced_acceleration = GetAcceleration(solid_data_i.pos_0_);
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
			: ConstraintByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>(body, body_part),
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
		//=================================================================================================//
		void DampingBySplittingAlgorithm
			::InnerInteraction(size_t index_particle_i, Real dt)
		{

			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			Vecd acceleration(0);
			NeighborList& inner_neighors
				= getNeighborList(reference_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j
					= particles_->base_particle_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j
					= particles_->elastic_body_data_[index_particle_j];

				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ neighboring_particle->r_ij_ ;
				acceleration += 0.5*eta_*vel_detivative
					*neighboring_particle->dW_ij_*elastic_data_j.mass_
					/ elastic_data_i.rho_0_ / elastic_data_j.rho_0_;
			}
			
			base_particle_data_i.dvel_dt_ += acceleration * 0.5;
			base_particle_data_i.vel_n_ += acceleration * dt* 0.5;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
