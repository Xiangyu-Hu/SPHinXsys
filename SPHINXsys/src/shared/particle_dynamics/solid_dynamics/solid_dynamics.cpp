/**
 * @file 	solid_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "solid_dynamics.h"
#include "general_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		ContactDensitySummation::
			ContactDensitySummation(SolidContactBodyRelation* solid_body_contact_relation) :
			PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
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
				StdLargeVec<Real>& contact_mass_k = *(contact_mass_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					sigma += contact_neighborhood.W_ij_[n] * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ContactForce::ContactForce(SolidContactBodyRelation* solid_body_contact_relation) :
			PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_, 
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
				StdLargeVec<Real>& contact_density_k = *(contact_contact_density_[k]);
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Solid* solid_k = contact_material_[k];

				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
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
			DynamicContactForce(SolidContactBodyRelation* solid_body_contact_relation, Real penalty_strength) :
			PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
				&solid_body_contact_relation->body_surface_layer_),
			ContactDynamicsData(solid_body_contact_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_),
			vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			contact_force_(particles_->contact_force_), penalty_strength_(penalty_strength)
		{
			Real impedence = material_->ReferenceDensity() * sqrt(material_->ContactStiffness());
			Real reference_pressure = material_->ReferenceDensity() * material_->ContactStiffness();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
				Real contact_impedence = 
					contact_material_[k]->ReferenceDensity()*sqrt(contact_material_[k]->ContactStiffness());
				contact_impedence_.push_back(2.0 * impedence * contact_impedence
					/ (impedence + contact_impedence));
				Real contact_reference_pressure =
					contact_material_[k]->ReferenceDensity() * contact_material_[k]->ContactStiffness();
				contact_reference_pressure_.push_back(2.0 * reference_pressure * contact_reference_pressure
					/ (reference_pressure + contact_reference_pressure));
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
				Real particle_spacing_ratio2 = 1.0 / (this->body_->particle_adaptation_->ReferenceSpacing() * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Solid* solid_k = contact_material_[k];

				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];

					Real impedence_p = 0.5 * contact_impedence_[k] * (SimTK::dot(vel_i - vel_n_k[index_j], -e_ij));
					Real overlap = contact_neighborhood.r_ij_[n];
					Real delta = 2.0 * overlap * particle_spacing_j1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty_p =penalty_strength_* beta* overlap* contact_reference_pressure_[k];
					
					//force due to pressure
					force -= 2.0 * (impedence_p + penalty_p) * e_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}
	
			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForceWithWall::
			ContactForceWithWall(SolidContactBodyRelation* solid_body_contact_relation, Real penalty_strength) :
			PartInteractionDynamicsByParticle(solid_body_contact_relation->sph_body_,
				&solid_body_contact_relation->body_surface_layer_),
			ContactDynamicsData(solid_body_contact_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_),
			vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			contact_force_(particles_->contact_force_), penalty_strength_(penalty_strength)
		{
			impedence_ = material_->ReferenceDensity() * sqrt(material_->ContactStiffness());
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
				Real particle_spacing_ratio2 = 1.0 / (body_->particle_adaptation_->ReferenceSpacing() * particle_spacing_j1);
				particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Solid* solid_k = contact_material_[k];

				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];
					Vecd n_k_j = n_k[index_j];

					Real impedence_p = 0.5 * impedence_ * (SimTK::dot(vel_i - vel_n_k[index_j], -n_k_j));
					Real overlap = contact_neighborhood.r_ij_[n] * SimTK::dot(n_k_j, e_ij);
					Real delta = 2.0 * overlap * particle_spacing_j1;
					Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
					Real penalty_p = penalty_strength_ * beta * overlap * reference_pressure_;

					//force due to pressure
					force -= 2.0 * (impedence_p + penalty_p) * dot(e_ij, n_k_j) *
						n_k_j * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			contact_force_[index_i] = force;
			dvel_dt_prior_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SolidBody* body, Real CFL) :
			ParticleDynamicsReduce<Real, ReduceMin>(body),
			ElasticSolidDataSimple(body), CFL_(CFL),
			vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_)
		{
			smoothing_length_ = particle_adaptation_->ReferenceSmoothingLength();
			initial_reference_ = DBL_MAX;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			Real sound_speed = material_->SoundWaveSpeed();
			return CFL_ * SMIN(sqrt(smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal)),
				smoothing_length_ / (sound_speed + vel_n_[index_i].norm()));
		}
		//=================================================================================================//
		CorrectConfiguration::
			CorrectConfiguration(BaseInnerBodyRelation* body_inner_relation) :
			InteractionDynamics(body_inner_relation->sph_body_),
			SolidDataInner(body_inner_relation),
			Vol_(particles_->Vol_), B_(particles_->B_)
		{
		}
		//=================================================================================================//
		void CorrectConfiguration::Interaction(size_t index_i, Real dt)
		{
			Matd local_configuration(Eps); // a small number added to diagonal to avoid divide zero
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration -= Vol_[index_j] * SimTK::outer(r_ji, gradw_ij);
			}
			B_[index_i] = SimTK::inverse(local_configuration);
		}
		//=================================================================================================//
		ConstrainSolidBodyRegion::
			ConstrainSolidBodyRegion(SPHBody* body, BodyPartByParticle* body_part) :
			PartSimpleDynamicsByParticle(body, body_part), SolidDataSimple(body),
			pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			n_(particles_->n_), n_0_(particles_->n_0_),
			vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
		}
		//=================================================================================================//
		void ConstrainSolidBodyRegion::Update(size_t index_i, Real dt)
		{
			Vecd pos_0 = pos_0_[index_i];
			Vecd pos_n = pos_n_[index_i];
			Vecd vel_n = vel_n_[index_i];
			Vecd dvel_dt = dvel_dt_[index_i];

			pos_n_[index_i] = getDisplacement(pos_0, pos_n);
			vel_n_[index_i] = getVelocity(pos_0, pos_n, vel_n);
			dvel_dt_[index_i] = getAcceleration(pos_0, pos_n, dvel_dt);
			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		SoftConstrainSolidBodyRegion::
			SoftConstrainSolidBodyRegion(BaseInnerBodyRelation* body_inner_relation, BodyPartByParticle* body_part) :
			PartInteractionDynamicsByParticleWithUpdate(body_inner_relation->sph_body_, body_part),
			SolidDataInner(body_inner_relation),
			Vol_(particles_->Vol_),	
			vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_),
			vel_temp_(*particles_->createAVariable<indexVector, Vecd>("TemporaryVelocity")),
			dvel_dt_temp_(*particles_->createAVariable<indexVector, Vecd>("TemporaryAcceleration")) {}
		//=================================================================================================//
		void SoftConstrainSolidBodyRegion::Interaction(size_t index_i, Real dt)
		{
			Real ttl_weight(Eps);			
			Vecd vel_i = vel_n_[index_i];
			Vecd dvel_dt_i = dvel_dt_[index_i];

			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real weight_j = inner_neighborhood.W_ij_[n] * Vol_[index_j];

				ttl_weight += weight_j;
				vel_i += vel_n_[index_j] * weight_j;
				dvel_dt_i += dvel_dt_[index_j] * weight_j;
			}

			vel_temp_[index_i] = vel_i / ttl_weight;
			dvel_dt_temp_[index_i] = dvel_dt_i / ttl_weight;
		}
		//=================================================================================================//
		void SoftConstrainSolidBodyRegion::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] = vel_temp_[index_i];
			dvel_dt_[index_i] = dvel_dt_temp_[index_i];
			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		ClampConstrainSolidBodyRegion::
			ClampConstrainSolidBodyRegion(BaseInnerBodyRelation* body_inner_relation, BodyPartByParticle* body_part) :
			ParticleDynamics<void>(body_inner_relation->sph_body_),
			constrianing_(new ConstrainSolidBodyRegion(body_inner_relation->sph_body_, body_part)),
			softing_(new SoftConstrainSolidBodyRegion(body_inner_relation, body_part)) {}
		//=================================================================================================//
		void ClampConstrainSolidBodyRegion::exec(Real dt)
		{
			constrianing_->exec();
			softing_->exec();
		}
		//=================================================================================================//
		void ClampConstrainSolidBodyRegion::parallel_exec(Real dt)
		{
			constrianing_->parallel_exec();
			softing_->parallel_exec();
		}
		//=================================================================================================//
		ConstrainSolidBodyMassCenter::
			ConstrainSolidBodyMassCenter(SPHBody* body, Vecd constrain_direction) :
			ParticleDynamicsSimple(body), SolidDataSimple(body),
			correction_matrix_(Matd(1.0)), vel_n_(particles_->vel_n_)
		{
			for (int i = 0; i != Dimensions; ++i) correction_matrix_[i][i] = constrain_direction[i];
			BodySummation<indexScalar, Real> compute_total_mass_(body, "Mass");
			total_mass_ = compute_total_mass_.parallel_exec();
			compute_total_momentum_ = new BodyMoment<indexVector, Vecd>(body, "Velocity");
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::setupDynamics(Real dt)
		{
			velocity_correction_ = 
				correction_matrix_ * compute_total_momentum_->parallel_exec(dt) / total_mass_;
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] -= velocity_correction_;
		}
		//=================================================================================================//
		ImposeExternalForce::
			ImposeExternalForce(SolidBody* body, SolidBodyPartForSimbody* body_part) :
			PartSimpleDynamicsByParticle(body, body_part), SolidDataSimple(body),
			pos_0_(particles_->pos_0_), vel_n_(particles_->vel_n_),
			vel_ave_(particles_->vel_ave_) {}
		//=================================================================================================//
		void ImposeExternalForce::Update(size_t index_i, Real dt)
		{
			Vecd induced_acceleration = getAcceleration(pos_0_[index_i]);
			vel_n_[index_i] += induced_acceleration * dt;
			vel_ave_[index_i] = vel_n_[index_i];
		}
		//=================================================================================================//
		ParticleWiseAcceleration
			::ParticleWiseAcceleration(SolidBody* body, Vecd stiffness)
			: ParticleDynamicsSimple(body), SolidDataSimple(body),
			mass_(particles_->mass_),
			pos_n_(particles_->pos_n_),
			pos_0_(particles_->pos_0_),
			dvel_dt_prior_(particles_->dvel_dt_prior_) {}
		//=================================================================================================//
		Vecd ParticleWiseAcceleration::getAcceleration(Vecd& disp, Real mass)
		{
			Vecd spring_force(0);
			for(int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i] / mass;
			}
			return spring_force;
		}
		//=================================================================================================//
		void ParticleWiseAcceleration::Update(size_t index_i, Real dt)
		{	
			Vecd disp_from_0 = pos_n_[index_i] - pos_0_[index_i];
			dvel_dt_prior_[index_i] += getAcceleration(disp_from_0, mass_[index_i]);
		}
		//=================================================================================================//	
		ElasticDynamicsInitialCondition::
			ElasticDynamicsInitialCondition(SolidBody* body) :
			ParticleDynamicsSimple(body),
			ElasticSolidDataSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		UpdateElasticNormalDirection::
			UpdateElasticNormalDirection(SolidBody* elastic_body) :
			ParticleDynamicsSimple(elastic_body),
			ElasticSolidDataSimple(elastic_body),
			n_(particles_->n_), n_0_(particles_->n_0_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		DeformationGradientTensorBySummation::
			DeformationGradientTensorBySummation(BaseInnerBodyRelation* body_inner_relation) :
			InteractionDynamics(body_inner_relation->sph_body_),
			ElasticSolidDataInner(body_inner_relation),
			Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			B_(particles_->B_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		void DeformationGradientTensorBySummation::Interaction(size_t index_i, Real dt)
		{
			Vecd& pos_n_i = pos_n_[index_i];

			Matd deformation(0.0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation -= Vol_[index_j] * SimTK::outer((pos_n_i - pos_n_[index_j]), gradw_ij);
			}

			F_[index_i] = B_[index_i] * deformation;
		}
		//=================================================================================================//
		BaseElasticRelaxation::
			BaseElasticRelaxation(BaseInnerBodyRelation* body_inner_relation) :
			ParticleDynamics1Level(body_inner_relation->sph_body_),
			ElasticSolidDataInner(body_inner_relation), Vol_(particles_->Vol_),
			rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) {}
		//=================================================================================================//
		StressRelaxationFirstHalf::
			StressRelaxationFirstHalf(BaseInnerBodyRelation* body_inner_relation) :
			BaseElasticRelaxation(body_inner_relation), 
			dvel_dt_prior_(particles_->dvel_dt_prior_), force_from_fluid_(particles_->force_from_fluid_),
			stress_PK1_(particles_->stress_PK1_)
		{
			rho0_ = material_->ReferenceDensity();
			inv_rho0_ = 1.0 / rho0_;
			smoothing_length_ = particle_adaptation_->ReferenceSmoothingLength();
			numerical_dissipation_factor_ = 0.25;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho0_ / det(F_[index_i]);
			//obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
			stress_PK1_[index_i] = F_[index_i] * material_->ConstitutiveRelation(F_[index_i], index_i);
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Interaction(size_t index_i, Real dt)
		{
			//including gravity and force from fluid
			Vecd acceleration = dvel_dt_prior_[index_i]
				+ force_from_fluid_[index_i] / mass_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Real dim_r_ij_1 = Dimensions / r_ij;
				Vecd pos_jump = pos_n_[index_i] - pos_n_[index_j];
				Vecd vel_jump = vel_n_[index_i] - vel_n_[index_j];
				Real strain_rate = SimTK::dot(pos_jump, vel_jump) * dim_r_ij_1 * dim_r_ij_1;
				Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
				Matd numerical_stress_ij = 0.5 * (F_[index_i] + F_[index_j]) * material_->NumericalDamping(strain_rate, smoothing_length_);
				acceleration += (stress_PK1_[index_i] + stress_PK1_[index_j] + numerical_dissipation_factor_ * weight * numerical_stress_ij)
					* inner_neighborhood.dW_ij_[n] * e_ij * Vol_[index_j] * inv_rho0_;
			}

			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Interaction(size_t index_i, Real dt)
		{
			Vecd& vel_n_i = vel_n_[index_i];

			Matd deformation_gradient_change_rate(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation_gradient_change_rate -= Vol_[index_j] * SimTK::outer((vel_n_i - vel_n_[index_j]), gradw_ij);
			}

			dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		ConstrainSolidBodyPartBySimBody::ConstrainSolidBodyPartBySimBody(SolidBody* body,
			SolidBodyPartForSimbody* body_part,
			SimTK::MultibodySystem& MBsystem,
			SimTK::MobilizedBody& mobod,
			SimTK::Force::DiscreteForces& force_on_bodies,
			SimTK::RungeKuttaMersonIntegrator& integ)
			: ConstrainSolidBodyRegion(body, body_part),
			MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
		void  ConstrainSolidBodyPartBySimBody::setupDynamics(Real dt)
		{
			body_->setNewlyUpdated();
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
		}
		//=================================================================================================//
		TotalForceOnSolidBodyPartForSimBody
			::TotalForceOnSolidBodyPartForSimBody(SolidBody* body,
				SolidBodyPartForSimbody* body_part,
				SimTK::MultibodySystem& MBsystem,
				SimTK::MobilizedBody& mobod,
				SimTK::Force::DiscreteForces& force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator& integ)
			: PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>(body, body_part),
			SolidDataSimple(body),
			force_from_fluid_(particles_->force_from_fluid_), contact_force_(particles_->contact_force_),
			pos_n_(particles_->pos_n_),
			MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			initial_reference_ = SpatialVec(Vec3(0), Vec3(0));
		}
		//=================================================================================================//
		void TotalForceOnSolidBodyPartForSimBody::SetupReduce()
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
	}
}
