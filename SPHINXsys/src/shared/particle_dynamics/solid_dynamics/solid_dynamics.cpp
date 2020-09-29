/**
 * @file 	solid_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		NormalDirectionSummation::
			NormalDirectionSummation(SPHBodyComplexRelation* body_complex_relation) :
			ParticleDynamicsComplex(body_complex_relation),
			SolidDataDelegateComplex(body_complex_relation),
			n_(particles_->n_), n_0_(particles_->n_0_) {}
		//=================================================================================================//
		void NormalDirectionSummation::ComplexInteraction(size_t index_i, Real dt)
		{
			Vecd gradient(0.0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				gradient += inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					gradient += contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			n_0_[index_i] = -gradient / (gradient.norm() + Eps);
			n_[index_i] = n_0_[index_i];
		}
		//=================================================================================================//
		NormalDirectionReNormalization::
			NormalDirectionReNormalization(SPHBodyComplexRelation* body_complex_relation)
			: NormalDirectionSummation(body_complex_relation),
			Vol_0_(particles_->Vol_0_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
				contact_Vol_0_.push_back(&(contact_particles_[k]->Vol_0_));
		}
		//=================================================================================================//
		void NormalDirectionReNormalization::ComplexInteraction(size_t index_i, Real dt)
		{
			Matd local_configuration(0.0);
			Vecd gradient(0.0);

			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration += Vol_0_[index_j] * SimTK::outer(r_ij, gradw_ij);
				gradient += gradw_ij * Vol_0_[index_j];
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_0_k = *(contact_Vol_0_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd gradw_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
					Vecd r_ij = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
					local_configuration += Vol_0_k[index_j] * SimTK::outer(r_ij, gradw_ij);
					gradient += gradw_ij * Vol_0_k[index_j];
				}
			}

			Matd correction_matrix = inverse(local_configuration);
			Vecd n_temp_ = ~correction_matrix * gradient;
			if (n_temp_.norm() <= 0.75) {
				n_0_[index_i] = Vecd(0.0);
			}
			else {
				n_0_[index_i] = -n_temp_ / (n_temp_.norm() + Eps);
			}
			n_[index_i] = n_0_[index_i];
		}
		//=================================================================================================//
		ElasticSolidDynamicsInitialCondition::
			ElasticSolidDynamicsInitialCondition(SolidBody* body) :
			ParticleDynamicsSimple(body),
			ElasticSolidDataDelegateSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		UpdateElasticNormalDirection::
			UpdateElasticNormalDirection(SolidBody* elastic_body) :
			ParticleDynamicsSimple(elastic_body),
			ElasticSolidDataDelegateSimple(elastic_body),
			n_(particles_->n_), n_0_(particles_->n_0_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		InitializeDisplacement::
			InitializeDisplacement(SolidBody* body, StdLargeVec<Vecd>& pos_temp) :
			ParticleDynamicsSimple(body), ElasticSolidDataDelegateSimple(body),
			pos_temp_(pos_temp), pos_n_(particles_->pos_n_),
			vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
		}
		//=================================================================================================//
		void InitializeDisplacement::Update(size_t index_i, Real dt)
		{
			pos_temp_[index_i] = pos_n_[index_i];
		}
		//=================================================================================================//
		void UpdateAverageVelocityAndAcceleration::Update(size_t index_i, Real dt)
		{
			Vecd updated_vel_ave = (pos_n_[index_i] - pos_temp_[index_i]) / (dt + TinyReal);
			dvel_dt_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + TinyReal);
			vel_ave_[index_i] = updated_vel_ave;
		}
		//=================================================================================================//
		AverageVelocityAndAcceleration::
			AverageVelocityAndAcceleration(SolidBody* body) :
			initialize_displacement_(body, pos_temp_),
			update_averages_(body, pos_temp_)
		{
			BaseParticles* base_particles = body->base_particles_;
			//register particle varibales defined in this class
			base_particles->registerAVariable(pos_temp_, base_particles->registered_vectors_,
				base_particles->vectors_map_, base_particles->vectors_to_write_, "TemporaryPosition", false);
		}
		//=================================================================================================//
		FluidViscousForceOnSolid::
			FluidViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation) :
			ParticleDynamicsContact(body_contact_relation),
			FSIDataDelegateContact(body_contact_relation),
			Vol_(particles_->Vol_), vel_ave_(particles_->vel_ave_),
			viscous_force_from_fluid_(particles_->viscous_force_from_fluid_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
			}

			//more work should be done for more general cases with multiple resolutions
			//and for fluids with different viscosities
			mu_ = contact_material_[0]->ReferenceViscosity();
			/** the smoothing length should be discuss more. */
			smoothing_length_ = powern(2.0, body_->refinement_level_) * body_->kernel_->GetSmoothingLength();
		}
		//=================================================================================================//
		void FluidViscousForceOnSolid::ContactInteraction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					//froce due to viscousity
					//viscous force with a simple wall model for high-Reynolds number flow
					Vecd vel_detivative = 2.0 * (vel_ave_i - vel_n_k[index_j])
						/ (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);

					force += 2.0 * mu_ * vel_detivative * Vol_i * Vol_k[index_j]
						   * contact_neighborhood.dW_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnSolid::ContactInteraction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real>& rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
					 * is formulation is more accurate thant the previsou one for Taygree-Vortex flow. */
					Real v_r_ij = dot(vel_ave_i - vel_n_k[index_j],
						contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n]);
					Real vel_difference = 0.0 * (vel_ave_i - vel_n_k[index_j]).norm()
						* contact_neighborhood.r_ij_[n];
					Real eta_ij = 8.0 * SMAX(mu_, rho_n_k[index_j] * vel_difference) * v_r_ij /
						(contact_neighborhood.r_ij_[n] * contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
					force += eta_ij * Vol_i * Vol_k[index_j]
						* contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		FluidPressureForceOnSolid::
			FluidPressureForceOnSolid(SPHBodyContactRelation* body_contact_relation) :
			FluidViscousForceOnSolid(body_contact_relation),
			force_from_fluid_(particles_->force_from_fluid_),
			dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_p_.push_back(&(contact_particles_[k]->p_));
				contact_dvel_dt_others_.push_back(&(contact_particles_[k]->dvel_dt_others_));
			}

		}
		//=================================================================================================//
		void FluidPressureForceOnSolid::ContactInteraction(size_t index_i, Real dt)
		{
			force_from_fluid_[index_i] = viscous_force_from_fluid_[index_i];
			Vecd& dvel_dt_ave_i = dvel_dt_ave_[index_i];
			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real>& rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Real>& p_k = *(contact_p_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				StdLargeVec<Vecd>& dvel_dt_others_k = *(contact_dvel_dt_others_[k]);
				Fluid* fluid_k = contact_material_[k];
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd e_ij = contact_neighborhood.e_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Real face_wall_external_acceleration
						= dot((dvel_dt_others_k[index_j] - dvel_dt_ave_i), e_ij);
					Real p_in_wall = p_k[index_j] + rho_n_k[index_j] * r_ij * SMAX(0.0, face_wall_external_acceleration);
					Real rho_in_wall = fluid_k->DensityFromPressure(p_in_wall);

					//solving Riemann problem
					Real p_star = fluid_k->RiemannSolverForPressure(rho_n_k[index_j], rho_in_wall, p_k[index_j], p_in_wall,
						dot(e_ij, vel_n_k[index_j]), dot(e_ij, vel_ave_i));;

					//force due to pressure
					force -= 2.0 * p_star * e_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			force_from_fluid_[index_i] += force;
		}
		//=================================================================================================//
		TotalViscousForceOnSolid
			::TotalViscousForceOnSolid(SolidBody* body) :
			ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>(body),
			SolidDataDelegateSimple(body),
			viscous_force_from_fluid_(particles_->viscous_force_from_fluid_)
		{
			initial_reference_ = Vecd(0);
		}
		//=================================================================================================//
		Vecd TotalViscousForceOnSolid::ReduceFunction(size_t index_i, Real dt)
		{
			return viscous_force_from_fluid_[index_i];
		}
		//=================================================================================================//
		TotalForceOnSolid::TotalForceOnSolid(SolidBody* body) :
			ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>(body),
			SolidDataDelegateSimple(body),
			force_from_fluid_(particles_->force_from_fluid_)
		{
			initial_reference_ = Vecd(0);
		}
		//=================================================================================================//
		Vecd TotalForceOnSolid::ReduceFunction(size_t index_i, Real dt)
		{
			return force_from_fluid_[index_i];
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SolidBody* body) :
			ParticleDynamicsReduce<Real, ReduceMin>(body),
			ElasticSolidDataDelegateSimple(body),
			vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			initial_reference_ = DBL_MAX;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			Real sound_speed = material_->ReferenceSoundSpeed();
			return 0.6 * SMIN(sqrt(smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal)),
				smoothing_length_ / (sound_speed + vel_n_[index_i].norm()));
		}
		//=================================================================================================//
		CorrectConfiguration::
			CorrectConfiguration(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamicsInner(body_inner_relation),
			SolidDataDelegateInner(body_inner_relation),
			Vol_0_(particles_->Vol_0_), B_(particles_->B_)
		{
		}
		//=================================================================================================//
		void CorrectConfiguration::InnerInteraction(size_t index_i, Real dt)
		{
			/** a small number added to diagnal to avoid divide zero */
			Matd local_configuration(Eps);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration += Vol_0_[index_j] * SimTK::outer(r_ji, gradw_ij);
			}
			/** note that I have changed to use stadrad linear solver here*/
			B_[index_i] = SimTK::inverse(local_configuration);
		}
		//=================================================================================================//
		DeformationGradientTensorBySummation::
			DeformationGradientTensorBySummation(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamicsInner(body_inner_relation),
			ElasticSolidDataDelegateInner(body_inner_relation),
			Vol_0_(particles_->Vol_0_), pos_n_(particles_->pos_n_),
			B_(particles_->B_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		void DeformationGradientTensorBySummation::InnerInteraction(size_t index_i, Real dt)
		{
			Vecd& pos_n_i = pos_n_[index_i];

			Matd deformation(0.0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation -= Vol_0_[index_j] * SimTK::outer((pos_n_i - pos_n_[index_j]), gradw_ij);
			}

			F_[index_i] = deformation * B_[index_i];
		}
		//=================================================================================================//
		StressRelaxationFirstHalf::
			StressRelaxationFirstHalf(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamicsInner1Level(body_inner_relation),
			ElasticSolidDataDelegateInner(body_inner_relation), Vol_0_(particles_->Vol_0_),
			rho_n_(particles_->rho_n_), rho_0_(particles_->rho_0_), mass_(particles_->mass_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			dvel_dt_others_(particles_->dvel_dt_others_), force_from_fluid_(particles_->force_from_fluid_),
			B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
			stress_(particles_->stress_)
		{
			numerical_viscosity_
				= material_->getNumericalViscosity(body_->kernel_->GetSmoothingLength());
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho_0_[index_i] / det(F_[index_i]);
			stress_[index_i] = material_->ConstitutiveRelation(F_[index_i], index_i)
				+ material_->NumericalDampingStress(F_[index_i], dF_dt_[index_i], numerical_viscosity_, index_i);
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;

		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::InnerInteraction(size_t index_i, Real dt)
		{
			Real rho_0_i = rho_0_[index_i];
			Matd& stress_i = stress_[index_i];
			Matd& B_i = B_[index_i];

			//including gravity and force from fluid
			Vecd acceleration = dvel_dt_others_[index_i]
				+ force_from_fluid_[index_i] / mass_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				acceleration += (stress_i * B_i + stress_[index_j] * B_[index_j])
					* inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_0_[index_j] / rho_0_i;
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
		void StressRelaxationSecondHalf::InnerInteraction(size_t index_i, Real dt)
		{
			Vecd& vel_n_i = vel_n_[index_i];

			Matd deformation_gradient_change_rate(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation_gradient_change_rate
					-= Vol_0_[index_j] * SimTK::outer((vel_n_i - vel_n_[index_j]), gradw_ij);
			}

			dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		ConstrainSolidBodyRegion::
			ConstrainSolidBodyRegion(SolidBody* body, BodyPartByParticle* body_part) :
			PartDynamicsByParticle(body, body_part), SolidDataDelegateSimple(body),
			pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			n_(particles_->n_), n_0_(particles_->n_0_),
			vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
		}
		//=================================================================================================//
		void ConstrainSolidBodyRegion::Update(size_t index_i, Real dt)
		{
			Point pos_0 = pos_0_[index_i];
			Point pos_n = pos_n_[index_i];
			Vecd vel_n = vel_n_[index_i];
			Vecd dvel_dt = dvel_dt_[index_i];

			pos_n_[index_i] = GetDisplacement(pos_0, pos_n);
			vel_n_[index_i] = GetVelocity(pos_0, pos_n, vel_n);
			dvel_dt_[index_i] = GetAcceleration(pos_0, pos_n, dvel_dt);
			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		ImposeExternalForce::
			ImposeExternalForce(SolidBody* body, SolidBodyPartForSimbody* body_part) :
			PartDynamicsByParticle(body, body_part), SolidDataDelegateSimple(body),
			pos_0_(particles_->pos_0_), vel_n_(particles_->vel_n_),
			vel_ave_(particles_->vel_ave_)
		{
		}
		//=================================================================================================//
		void ImposeExternalForce::Update(size_t index_i, Real dt)
		{
			Vecd induced_acceleration = GetAcceleration(pos_0_[index_i]);
			vel_n_[index_i] += induced_acceleration * dt;
			vel_ave_[index_i] = vel_n_[index_i];
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
			SolidDataDelegateSimple(body),
			force_from_fluid_(particles_->force_from_fluid_), pos_n_(particles_->pos_n_),
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
		DampingBySplittingAlgorithm
			::DampingBySplittingAlgorithm(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamicsInnerSplitting(body_inner_relation),
			ElasticSolidDataDelegateInner(body_inner_relation),
			Vol_0_(particles_->Vol_0_), mass_(particles_->mass_), vel_n_(particles_->vel_n_),
			eta_(material_->getPhysicalViscosity())
		{

		}
		//=================================================================================================//
		void DampingBySplittingAlgorithm
			::InnerInteraction(size_t index_i, Real dt)
		{

			Real Vol_0_i = Vol_0_[index_i];
			Real mass_i = mass_[index_i];
			Vecd& vel_n_i = vel_n_[index_i];

			Vecd error(0);
			Real parameter_a(0);
			Real parameter_c(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//linear projection 
				Vecd vel_detivative = (vel_n_i - vel_n_[index_j]);
				Real parameter_b = 2.0 * eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_0_i * Vol_0_[index_j] * dt / inner_neighborhood.r_ij_[n];

				error -= vel_detivative * parameter_b;
				parameter_a += parameter_b;
				parameter_c += parameter_b * parameter_b;
			}

			parameter_a -= mass_i;
			Real parameter_l = parameter_a * parameter_a + parameter_c;
			Vecd parameter_k = error / (parameter_l + TinyReal);
			vel_n_[index_i] += parameter_k * parameter_a;

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Real parameter_b = 2.0 * eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_0_i * Vol_0_[index_j] * dt / inner_neighborhood.r_ij_[n];

				//predicted velocity at particle j
				Vecd vel_j = vel_n_[index_j] - parameter_k * parameter_b;
				Vecd vel_detivative = (vel_n_[index_j] - vel_j);

				//viscous force in conservation form
				vel_n_[index_i] -= vel_detivative * parameter_b / mass_[index_i];
			}
		}
		//=================================================================================================//
		DampingBySplittingPairwise
			::DampingBySplittingPairwise(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamicsInnerSplitting(body_inner_relation),
			ElasticSolidDataDelegateInner(body_inner_relation),
			Vol_0_(particles_->Vol_0_), mass_(particles_->mass_), vel_n_(particles_->vel_n_),
			eta_(material_->getPhysicalViscosity())
		{

		}
		//=================================================================================================//
		void DampingBySplittingPairwise
			::InnerInteraction(size_t index_i, Real dt)
		{

			Real Vol_0_i = Vol_0_[index_i];
			Real mass_i = mass_[index_i];
			Vecd& vel_n_i = vel_n_[index_i];

			StdVec<Real> parameter_b(50);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			size_t number_of_neighbors = inner_neighborhood.current_size_;
			parameter_b.resize(number_of_neighbors);
			//forward sweep
			for (size_t n = 0; n != number_of_neighbors; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real mass_j = mass_[index_j];

				Vecd vel_detivative = (vel_n_i - vel_n_[index_j]);
				parameter_b[n] = eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_0_i * Vol_0_[index_j] * dt / inner_neighborhood.r_ij_[n];

				Vecd increment = parameter_b[n] * vel_detivative
					/ (mass_i * mass_j - parameter_b[n] * (mass_i + mass_j));
				vel_n_[index_i] += increment * mass_j;
				vel_n_[index_j] -= increment * mass_i;
			}

			//backward sweep
			for (size_t n = number_of_neighbors; n != 0; --n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real mass_j = mass_[index_j];

				Vecd vel_detivative = (vel_n_i - vel_n_[index_j]);
				Vecd increment = parameter_b[n - 1] * vel_detivative
					/ (mass_i * mass_j - parameter_b[n - 1] * (mass_i + mass_j));

				vel_n_[index_i] += increment * mass_j;
				vel_n_[index_j] -= increment * mass_i;
			}
		}
		//=================================================================================================//
		DampingBySplittingWithRandomChoice
			::DampingBySplittingWithRandomChoice(SPHBodyInnerRelation* body_inner_relation, Real random_ratio)
			: DampingBySplittingAlgorithm(body_inner_relation), random_ratio_(random_ratio)
		{
			eta_ = eta_ / random_ratio_;
		}
		//=================================================================================================//
		bool DampingBySplittingWithRandomChoice::RandomChoice()
		{
			return ((double)rand() / (RAND_MAX)) < random_ratio_ ? true : false;
		}
		//=================================================================================================//
		void DampingBySplittingWithRandomChoice::exec(Real dt)
		{
			if (RandomChoice()) DampingBySplittingAlgorithm::exec(dt);
		}
		//=================================================================================================//
		void DampingBySplittingWithRandomChoice::parallel_exec(Real dt)
		{
			if (RandomChoice()) DampingBySplittingAlgorithm::parallel_exec(dt);
		}
		//=================================================================================================//
		FluidViscousForceOnSolidWallModel::
			FluidViscousForceOnSolidWallModel(SPHBodyContactRelation* body_contact_relation,
				fluid_dynamics::ViscousAccelerationWallModel* viscous_acceleration_wall_modeling)
			: FluidViscousForceOnSolid(body_contact_relation),
			n_(particles_->n_),
			gradient_p_(viscous_acceleration_wall_modeling->gradient_p_),
			gradient_vel_(viscous_acceleration_wall_modeling->gradient_vel_) {};
		//=================================================================================================//
		void FluidViscousForceOnSolidWallModel::ContactInteraction(size_t index_i, Real dt)
		{
			FluidViscousForceOnSolid::ContactInteraction(index_i, dt);

			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real>& rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					//obtain the unit vector in the tangential direction
					Real height = 0.5 * contact_neighborhood.r_ij_[n];
					Vecd vel_to_wall = (vel_n_k[index_j] - vel_ave_i);
					Real vel_to_wall_n = SimTK::dot(vel_to_wall, n_[index_i]);
					Vecd vel_to_wall_t = vel_to_wall - vel_to_wall_n * n_[index_i];
					Vecd unit_t = vel_to_wall_t / (vel_to_wall_t.norm() + TinyReal);

					Real coefficient_A = SimTK::dot(gradient_p_[index_j], unit_t);
					Real coefficient_B = SimTK::dot(gradient_vel_[index_j] * vel_to_wall, unit_t);
					Vecd v_derivative_add = - height * (coefficient_A + rho_n_k[index_j] * coefficient_B) * unit_t;

					force -= 2.0 * v_derivative_add * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] += force;
		}
		//=================================================================================================//
	}
}
