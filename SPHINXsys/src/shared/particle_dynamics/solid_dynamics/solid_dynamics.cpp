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
		SummationContactDensity::
			SummationContactDensity(SolidBodyContactRelation* soild_body_contact_relation) :
			PartInteractionDynamicsByParticle(soild_body_contact_relation->sph_body_,
				&soild_body_contact_relation->body_surface_layer_),
			ContactDynamicsDataDelegate(soild_body_contact_relation),
			mass_(particles_->mass_), contact_density_(particles_->contact_density_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		void SummationContactDensity::Interaction(size_t index_i, Real dt)
		{
			/** Contact interaction. */
			Real sigma = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& contact_mass_k = *(contact_mass_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					sigma += contact_neighborhood.W_ij_[n] * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			contact_density_[index_i] = sigma;
		}
		//=================================================================================================//
		ContactForce::ContactForce(SolidBodyContactRelation* soild_body_contact_relation) :
			PartInteractionDynamicsByParticle(soild_body_contact_relation->sph_body_, 
				&soild_body_contact_relation->body_surface_layer_),
			ContactDynamicsDataDelegate(soild_body_contact_relation),
			contact_density_(particles_->contact_density_), 
			Vol_(particles_->Vol_), mass_(particles_->mass_),
			dvel_dt_others_(particles_->dvel_dt_others_),
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

				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
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
			dvel_dt_others_[index_i] += force / mass_[index_i];
		}
		//=================================================================================================//
		ContactForceFromFriction::ContactForceFromFriction(SPHBodyContactRelation* body_contact_relation,
			StdLargeVec<Vecd>& vel_n, Real eta) :
			InteractionDynamicsSplitting(body_contact_relation->sph_body_),
			ContactDynamicsDataDelegate(body_contact_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_),
			contact_force_(particles_->contact_force_), vel_n_(vel_n), eta_(eta)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
				contact_contact_force_.push_back(&(contact_particles_[k]->contact_force_));
			}
		}
		//=================================================================================================//
		void ContactForceFromFriction::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Real mass_i = mass_[index_i];
			Vecd& vel_n_i = vel_n_[index_i];

			StdVec<Real> parameter_b(MaximumNeighborhoodSize);
			size_t neighbors = 0;
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real>& mass_k = *(contact_mass_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				StdLargeVec<Vecd>& contact_force_k = *(contact_contact_force_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				//forward sweep
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real mass_j = mass_k[index_j];

					Vecd velocity_derivative = (vel_n_i - vel_n_k[index_j]);
					parameter_b[neighbors] = eta_ * contact_neighborhood.dW_ij_[n]
						* Vol_i * Vol_k[index_j] * dt / contact_neighborhood.r_ij_[n];
					Vecd increment = parameter_b[neighbors] * velocity_derivative
						/ (mass_i * mass_j - parameter_b[neighbors] * (mass_i + mass_j));
					neighbors++;

					vel_n_[index_i] += increment * mass_j;
					vel_n_k[index_j] -= increment * mass_i;
					contact_force_[index_i] += increment * mass_i * mass_j / (dt + TinyReal);
					contact_force_k[index_j] -= increment * mass_i * mass_j / (dt + TinyReal);
				}
			}

			/** Contact interaction. */
			for (size_t k = contact_configuration_.size(); k != 0; --k)
			{
				StdLargeVec<Real>& mass_k = *(contact_mass_[k - 1]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k - 1]);
				StdLargeVec<Vecd>& contact_force_k = *(contact_contact_force_[k - 1]);
				Neighborhood& contact_neighborhood = contact_configuration_[k - 1][index_i];
				//backward sweep
				for (size_t n = contact_neighborhood.current_size_; n != 0; --n)
				{
					size_t index_j = contact_neighborhood.j_[n - 1];
					Real mass_j = mass_k[index_j];

					Vecd velocity_derivative = (vel_n_i - vel_n_k[index_j]);
					neighbors--;
					Vecd increment = parameter_b[neighbors] * velocity_derivative
						/ (mass_i * mass_j - parameter_b[neighbors] * (mass_i + mass_j));

					vel_n_[index_i] += increment * mass_j;
					vel_n_k[index_j] -= increment * mass_i;
					contact_force_[index_i] += increment * mass_i * mass_j / (dt + TinyReal);
					contact_force_k[index_j] -= increment * mass_i * mass_j / (dt + TinyReal);
				}
			}
		}
		//=================================================================================================//
		FluidViscousForceOnSolid::
			FluidViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation) :
			InteractionDynamics(body_contact_relation->sph_body_),
			FSIDataDelegateContact(body_contact_relation),
			Vol_(particles_->Vol_), vel_ave_(particles_->vel_ave_),
			viscous_force_from_fluid_(particles_->viscous_force_from_fluid_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));

				mu_.push_back(contact_material_[k]->ReferenceViscosity());
				smoothing_length_.push_back(contact_bodies_[k]->kernel_->GetSmoothingLength());
			}
		}
		//=================================================================================================//
		void FluidViscousForceOnSolid::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					//force due to viscosity
					//viscous force with a simple wall model for high-Reynolds number flow
					Vecd vel_derivative = 2.0 * (vel_ave_i - vel_n_k[index_j])
						/ (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);

					force += 2.0 * mu_k * vel_derivative * Vol_i * Vol_k[index_j]
						   * contact_neighborhood.dW_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnSolid::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd& vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real>& rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
					 * is formulation is more accurate thant the previous one for Taylor-Green-Vortex flow. */
					Real v_r_ij = dot(vel_ave_i - vel_n_k[index_j],
						contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n]);
					Real vel_difference = 0.0 * (vel_ave_i - vel_n_k[index_j]).norm()
						* contact_neighborhood.r_ij_[n];
					Real eta_ij = 8.0 * SMAX(mu_k, rho_n_k[index_j] * vel_difference) * v_r_ij /
						(contact_neighborhood.r_ij_[n] * contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);
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
		void FluidPressureForceOnSolid::Interaction(size_t index_i, Real dt)
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
			InteractionDynamics(body_inner_relation->sph_body_),
			SolidDataDelegateInner(body_inner_relation),
			Vol_(particles_->Vol_), B_(particles_->B_)
		{
		}
		//=================================================================================================//
		void CorrectConfiguration::Interaction(size_t index_i, Real dt)
		{
			/** a small number added to diagonal to avoid divide zero */
			Matd local_configuration(Eps);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration -= Vol_[index_j] * SimTK::outer(r_ji, gradw_ij);
			}
			/** note that I have changed to use standard linear solver here*/
			B_[index_i] = SimTK::inverse(local_configuration);
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
		FluidViscousForceOnSolidWallModel::
			FluidViscousForceOnSolidWallModel(SPHBodyContactRelation* body_contact_relation,
				fluid_dynamics::ViscousAccelerationWallModel* viscous_acceleration_wall_modeling)
			: FluidViscousForceOnSolid(body_contact_relation),
			n_(particles_->n_),
			gradient_p_(viscous_acceleration_wall_modeling->gradient_p_),
			gradient_vel_(viscous_acceleration_wall_modeling->gradient_vel_) {};
		//=================================================================================================//
		void FluidViscousForceOnSolidWallModel::Interaction(size_t index_i, Real dt)
		{
			FluidViscousForceOnSolid::Interaction(index_i, dt);

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
		DeformationGradientTensorBySummation::
			DeformationGradientTensorBySummation(SPHBodyInnerRelation* body_inner_relation) :
			InteractionDynamics(body_inner_relation->sph_body_),
			ElasticSolidDataDelegateInner(body_inner_relation),
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

			F_[index_i] = deformation * B_[index_i];
		}
		//=================================================================================================//
		StressRelaxationFirstHalf::
			StressRelaxationFirstHalf(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamics1Level(body_inner_relation->sph_body_),
			ElasticSolidDataDelegateInner(body_inner_relation), Vol_(particles_->Vol_),
			rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			dvel_dt_others_(particles_->dvel_dt_others_), force_from_fluid_(particles_->force_from_fluid_),
			B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
			stress_(particles_->stress_)
		{
			rho_0_ = material_->ReferenceDensity();
			inv_rho_0_ = 1.0 / rho_0_;
			numerical_viscosity_
				= material_->getNumericalViscosity(body_->kernel_->GetSmoothingLength());
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho_0_ / det(F_[index_i]);
			//obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress,
			// including numerical disspation stress  
			stress_[index_i] = F_[index_i] * (material_->ConstitutiveRelation(F_[index_i], index_i)
				+ material_->NumericalDampingStress(F_[index_i], dF_dt_[index_i], numerical_viscosity_, index_i));
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;

		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Interaction(size_t index_i, Real dt)
		{
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
					* inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j] * inv_rho_0_;
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
				deformation_gradient_change_rate
					-= Vol_[index_j] * SimTK::outer((vel_n_i - vel_n_[index_j]), gradw_ij);
			}

			dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
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
			//register particle variables defined in this class
			base_particles->registerAVariable(pos_temp_, base_particles->registered_vectors_,
				base_particles->vectors_map_, base_particles->vectors_to_write_, "TemporaryPosition", false);
		}
		//=================================================================================================//
	}
}
