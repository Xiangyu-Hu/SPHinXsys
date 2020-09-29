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
		FluidInitialCondition::
			FluidInitialCondition(FluidBody* body)
			: ParticleDynamicsSimple(body), FluidDataDelegateSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		DensityBySummation::DensityBySummation(SPHBodyComplexRelation* body_complex_relation) :
			ParticleDynamicsComplex(body_complex_relation), FluidDataDelegateComplex(body_complex_relation),
			Vol_(particles_->Vol_), Vol_0_(particles_->Vol_0_), sigma_0_(particles_->sigma_0_), 
			rho_n_(particles_->rho_n_), rho_0_(particles_->rho_0_), mass_(particles_->mass_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
				contact_Vol_0_.push_back(&(contact_particles_[k]->Vol_0_));
			W0_ = body_->kernel_->W(Vecd(0));
		}
		//=================================================================================================//
		void DensityBySummation::ComplexInteraction(size_t index_i, Real dt)
		{
			Real Vol_0_i = Vol_0_[index_i];

			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_0_k = *(contact_Vol_0_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					sigma += contact_neighborhood.W_ij_[n] * Vol_0_k[contact_neighborhood.j_[n]] / Vol_0_i;
				}
			}

			Real rho_sum = sigma * rho_0_[index_i] / sigma_0_[index_i];
			rho_n_[index_i] = ReinitializedDensity(rho_sum, rho_0_[index_i], rho_n_[index_i]);
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		ViscousAcceleration::ViscousAcceleration(SPHBodyComplexRelation* body_complex_relation) :
			ParticleDynamicsComplex(body_complex_relation), FluidDataDelegateComplex(body_complex_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->mass_), 
			vel_n_(particles_->vel_n_), dvel_dt_others_(particles_->dvel_dt_others_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_ave_.push_back(&(contact_particles_[k]->vel_ave_));
			}
			mu_ = material_->ReferenceViscosity();
			smoothing_length_ = body_->kernel_->GetSmoothingLength();
		}		
		//=================================================================================================//
		void ViscousAcceleration::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			/** Inner interaction. */
			Vecd acceleration(0);
			Vecd vel_derivative(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//viscous force
				vel_derivative = (vel_i - vel_n_[index_j])
								/ (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative 
								* Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}


			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0*(vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * smoothing_length_);
					acceleration += 2.0 * mu_ * vel_derivative 
								  * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
				}
			}

			/** Particle summation. */
			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		void AngularConservativeViscousAcceleration::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			/** Inner interaction. */
			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
				Real v_r_ij = dot(vel_i - vel_n_[index_j], r_ij * e_ij);
				Real eta_ij = 8.0 * mu_  * v_r_ij /	(r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * Vol_[index_j] / rho_i
								* inner_neighborhood.dW_ij_[n] * e_ij;
			}
			
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
					 * is formulation is more accurate thant the previous one for Taylor-Green-Vortex flow. */
					Real v_r_ij = 2.0 * dot(vel_i - vel_ave_k[index_j], r_ij * e_ij);
					Real vel_difference = 0.0 * (vel_i - vel_ave_k[index_j]).norm() * r_ij;
					Real eta_ij = 8.0 * SMAX(mu_, rho_i * vel_difference) * v_r_ij /
						(r_ij * r_ij + 0.01 * smoothing_length_);
					acceleration += eta_ij * Vol_k[index_j] / rho_i
						* inner_neighborhood.dW_ij_[n] * e_ij;
				}
			}

			/** Particle summation. */
			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrection
			::TransportVelocityCorrection(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans)
			: ParticleDynamicsComplex(body_complex_relation), FluidDataDelegateComplex(body_complex_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), pos_n_(particles_->pos_n_),
			dvel_dt_trans_(dvel_dt_trans), p_background_(0)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void TransportVelocityCorrection::setupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->ReferenceDensity();
			p_background_ =  10.0 * density * speed_max * speed_max;
		}
		//=================================================================================================//
		void TransportVelocityCorrection::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];

			/** Inner interaction. */
			Vecd acceleration_trans(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_*Vol_[index_j] * nablaW_ij / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * Vol_k[index_j] * nablaW_ij / rho_i;
				}
			}

			dvel_dt_trans_[index_i] = acceleration_trans;
			/** correcting particle position */
			pos_n_[index_i] += acceleration_trans * dt*dt*0.5;
		}
		//=================================================================================================//
		TransportVelocityStress::
			TransportVelocityStress(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans) : 
			TransportVelocityCorrection(body_complex_relation, dvel_dt_trans),
			vel_n_(particles_->vel_n_), dvel_dt_others_(particles_->dvel_dt_others_)
		{
		}
		//=================================================================================================//
		void TransportVelocityStress::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];
			Vecd& dvel_dt_trans_i = dvel_dt_trans_[index_i];

			/** Inner interaction. */
			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				//extra stress
				acceleration += 0.5 * dt * dW_ij  * (rho_i * SimTK::dot(e_ij, vel_i) * dvel_dt_trans_i
					+ rho_n_[index_j] * SimTK::dot(e_ij, vel_n_[index_j]) * dvel_dt_trans_[index_j]) * Vol_[index_j] / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real dW_ij = inner_neighborhood.dW_ij_[n];
					Vecd& e_ij = inner_neighborhood.e_ij_[n];

					//extra stress
					acceleration += 0.5 * dt * dW_ij * SimTK::dot(e_ij, vel_i) *  dvel_dt_trans_i * Vol_k[index_j];
				}
			}

			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityFormulation::
			TransportVelocityFormulation(SPHBodyComplexRelation* body_complex_relation) : 
			correction_(body_complex_relation, dvel_dt_trans_), stress_(body_complex_relation, dvel_dt_trans_) 
		{
			SPHBody* sph_body = body_complex_relation->sph_body_;
			BaseParticles* base_particles = sph_body->base_particles_;
			//register particle variable defined in this class
			base_particles->registerAVariable(dvel_dt_trans_, base_particles->registered_vectors_, 
				base_particles->vectors_map_, base_particles->vectors_to_write_, "TransportAccesleration", false);
		}
		//=================================================================================================//
		TotalMechanicalEnergy::TotalMechanicalEnergy(FluidBody* body, Gravity* gravity)
			: ParticleDynamicsReduce<Real, ReduceSum<Real>>(body), 
			FluidDataDelegateSimple(body), mass_(particles_->mass_), 
			vel_n_(particles_->vel_n_), pos_n_(particles_->pos_n_), gravity_(gravity)
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real TotalMechanicalEnergy::ReduceFunction(size_t index_i, Real dt)
		{
			return 0.5 * mass_[index_i] * vel_n_[index_i].normSqr()
				+ mass_[index_i] * gravity_->getPotential(pos_n_[index_i]);
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(FluidBody* body)
			: ParticleDynamicsReduce<Real, ReduceMax>(body),
			FluidDataDelegateSimple(body), rho_n_(particles_->rho_n_),
			p_(particles_->p_), vel_n_(particles_->vel_n_)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time step size due to linear viscosity
			Real rho_0 = material_->ReferenceDensity();
			Real mu = material_->ReferenceViscosity();
			initial_reference_ = mu / rho_0 / smoothing_length_;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return material_->GetSoundSpeed(p_[index_i], rho_n_[index_i]) + vel_n_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			particles_->signal_speed_max_ = reduced_value;
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(FluidBody* body, Real U_max)
			: AcousticTimeStepSize(body)
		{
			Real u_max = SMAX(initial_reference_, U_max);
			initial_reference_ = u_max * u_max;
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return vel_n_[index_i].normSqr();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::OutputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			particles_->speed_max_ = speed_max;
			return 0.25 * smoothing_length_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		VorticityInFluidField::
			VorticityInFluidField(SPHBodyInnerRelation* body_inner_relation) : 
			ParticleDynamicsInner(body_inner_relation), FluidDataDelegateInner(body_inner_relation),
			Vol_(particles_->Vol_), vel_n_(particles_->vel_n_), vorticity_(particles_->vorticity_) {};
		//=================================================================================================//
		void VorticityInFluidField::InnerInteraction(size_t index_i, Real dt)
		{
			Vecd& vel_i = vel_n_[index_i];

			Vecd vorticity(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd r_ij = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

				Vecd vel_diff = vel_i - vel_n_[index_j];
				vorticity += upgradeVector<Vecd>(SimTK::cross(vel_diff, r_ij)) * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}

			vorticity_[index_i] = vorticity;
		}
		//=================================================================================================//
		PressureRelaxationFirstHalfRiemann::
			PressureRelaxationFirstHalfRiemann(SPHBodyComplexRelation* body_complex_relation) : 
			ParticleDynamicsComplex1Level(body_complex_relation), 
			FluidDataDelegateComplex(body_complex_relation),  
			Vol_(particles_->Vol_), mass_(particles_->mass_), rho_n_(particles_->rho_n_), 
			p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			dvel_dt_(particles_->dvel_dt_), dvel_dt_others_(particles_->dvel_dt_others_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_ave_.push_back(&(contact_particles_[k]->vel_ave_));
				contact_dvel_dt_ave_.push_back(&(contact_particles_[k]->dvel_dt_ave_));
				contact_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::Initialization(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
			p_[index_i] = material_->GetPressure(rho_n_[index_i]);
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Real p_i = p_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			Vecd acceleration = dvel_dt_others_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				/** Solving Riemann problem or not. */
				Real p_star = getPStar(e_ij, vel_i, p_i, rho_i,	vel_n_[index_j], p_[index_j], rho_n_[index_j]);
			
				acceleration -= 2.0 * p_star * Vol_[index_j] * dW_ij * e_ij / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Vecd& dvel_dt_others_i = dvel_dt_others_[index_i];

				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				StdLargeVec<Vecd>& dvel_dt_ave_k = *(contact_dvel_dt_ave_[k]);
				StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					Real face_wall_external_acceleration
						= dot((dvel_dt_others_i - dvel_dt_ave_k[index_j]), -e_ij);
					Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_i;
					Real p_in_wall = p_i + rho_i * r_ij * SMAX(0.0, face_wall_external_acceleration);
					Real rho_in_wall = material_->DensityFromPressure(p_in_wall);

					//solving Riemann problem or not
					Real p_star = getPStar(n_k[index_j], vel_i, p_i, rho_i, vel_in_wall, p_in_wall, rho_in_wall);

					//pressure force
					acceleration -= 2.0 * p_star * e_ij * Vol_k[index_j] * dW_ij / rho_i;
				}
			}
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		Real PressureRelaxationFirstHalfRiemann::getPStar(Vecd& e_ij,
			Vecd& vel_i, Real p_i, Real rho_i, Vecd& vel_j, Real p_j, Real rho_j)
		{
			//low dissipation Riemann problem
			return material_->RiemannSolverForPressure(rho_i, rho_j, p_i, p_j, dot(-e_ij, vel_i), dot(-e_ij, vel_j));
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfRiemann::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		Real PressureRelaxationFirstHalf::getPStar(Vecd& e_ij,
			Vecd& vel_i, Real p_i, Real rho_i, Vecd& vel_j, Real p_j, Real rho_j)
		{
			return (p_i * rho_j + p_j * rho_i) 	/ (rho_i + rho_j);;
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
//=================================================================================================//
		void PressureRelaxationSecondHalfRiemann::ComplexInteraction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Real p_i = p_[index_i];
			Vecd vel_i = vel_n_[index_i];

			Real density_change_rate = 0.0;
			Vecd vel_star(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];

				/** Solving Riemann problem or not. */
				vel_star = getVStar(e_ij, vel_i, p_i, rho_i, vel_n_[index_j], p_[index_j], rho_n_[index_j]);

				density_change_rate += 2.0 * rho_i * Vol_[index_j] * dot(vel_i - vel_star, e_ij) * dW_ij;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Vecd& dvel_dt_others_i = dvel_dt_others_[index_i];

				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				StdLargeVec<Vecd>& dvel_dt_ave_k = *(contact_dvel_dt_ave_[k]);
				StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];

					Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_i;
					Real face_wall_external_acceleration
						= dot((dvel_dt_others_i - dvel_dt_ave_k[index_j]), e_ij);
					Real p_in_wall = p_i + rho_i * r_ij * SMAX(0.0, face_wall_external_acceleration);
					Real rho_in_wall = material_->DensityFromPressure(p_in_wall);

					//solving Riemann problem or not
					vel_star = getVStar(n_k[index_j], vel_i, p_i, rho_i, vel_in_wall, p_in_wall, rho_in_wall);

					density_change_rate += 2.0 * rho_i * Vol_k[index_j]	* dot(vel_i - vel_star, e_ij) * dW_ij;
				}
			}

			drho_dt_[index_i] = density_change_rate;
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
		void PressureRelaxationSecondHalfRiemann::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		Vecd PressureRelaxationSecondHalf::getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
			Vecd& vel_j, Real p_j, Real rho_j)
		{
			return 0.5 * (vel_i + vel_j);
		}
		//=================================================================================================//
		PressureRelaxationFirstHalfOldroyd_B
			::PressureRelaxationFirstHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation)
			: PressureRelaxationFirstHalfRiemann(body_complex_relation),
			tau_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->tau_),
			dtau_dt_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->dtau_dt_)
		{
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::Initialization(size_t index_i, Real dt)
		{
			PressureRelaxationFirstHalfRiemann::Initialization(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationFirstHalfOldroyd_B::ComplexInteraction(size_t index_i, Real dt)
		{
			PressureRelaxationFirstHalfRiemann::ComplexInteraction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//elastic force
				acceleration += (tau_i + tau_[index_j]) * nablaW_ij * Vol_[index_j] / rho_i;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ij * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		PressureRelaxationSecondHalfOldroyd_B
			::PressureRelaxationSecondHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation)
			: PressureRelaxationSecondHalfRiemann(body_complex_relation),
			tau_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->tau_),
			dtau_dt_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->dtau_dt_) 
		{
			Oldroyd_B_Fluid *oldroy_b_fluid 
				= dynamic_cast<Oldroyd_B_Fluid*>(body_->base_particles_->base_material_);
			mu_p_ = oldroy_b_fluid->ReferencePolymericViscosity();
			lambda_ = oldroy_b_fluid->getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::ComplexInteraction(size_t index_i, Real dt)
		{
			PressureRelaxationSecondHalfRiemann::ComplexInteraction(index_i, dt);
			
			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				Matd velocity_gradient = - SimTK::outer((vel_i - vel_n_[index_j]), nablaW_ij) * Vol_[index_j];
				stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient 
					- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					Matd velocity_gradient = - SimTK::outer((vel_i - vel_ave_k[index_j]), nablaW_ij) * Vol_k[index_j] * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient
						- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}

			dtau_dt_[index_i] = stress_rate;
		}
		//=================================================================================================//
		void PressureRelaxationSecondHalfOldroyd_B::Update(size_t index_i, Real dt)
		{
			PressureRelaxationSecondHalfRiemann::Update(index_i, dt);

			tau_[index_i] +=  dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		FlowRelaxationBuffer::
			FlowRelaxationBuffer(FluidBody* body, BodyPartByCell* body_part) :
			PartDynamicsByCell(body, body_part), FluidDataDelegateSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), relaxation_rate_(0.1)
		{
		};
		//=================================================================================================//
		void FlowRelaxationBuffer
			::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += 
				relaxation_rate_ * ( getTargetVelocity(pos_n_[index_i], vel_n_[index_i]) - vel_n_[index_i]);
		}
		//=================================================================================================//
		DampingBoundaryCondition::
			DampingBoundaryCondition(FluidBody* body, BodyPartByCell* body_part) :
			PartDynamicsByCell(body, body_part), FluidDataDelegateSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), strength_(5.0)
		{
			body_part->BodyPartBounds(damping_zone_lower_bound_, damping_zone_upper_bound_);
		};
		//=================================================================================================//
		void DampingBoundaryCondition::Update(size_t index_i, Real dt)
		{
			Real damping_factor = (pos_n_[index_i][0] - damping_zone_lower_bound_[0]) /
								  (damping_zone_upper_bound_[0]-damping_zone_lower_bound_[0]);
			vel_n_[index_i] *=  (1.0 - dt * strength_ * damping_factor * damping_factor);
		}
		//=================================================================================================//
		EmitterInflowCondition::
			EmitterInflowCondition(FluidBody* body, BodyPartByParticle* body_part) :
			PartDynamicsByParticle(body, body_part), FluidDataDelegateSimple(body),
			rho_n_(particles_->rho_n_), rho_0_(particles_->rho_0_), p_(particles_->p_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), inflow_pressure_(0)
		{

		}
		//=================================================================================================//
		void EmitterInflowCondition
			::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] = getTargetVelocity(pos_n_[index_i], vel_n_[index_i]);
			rho_n_[index_i] = rho_0_[index_i];
			p_[index_i] = material_->GetPressure(rho_n_[index_i]);
		}
		//=================================================================================================//
		EmitterInflowInjecting
			::EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_width, int axis_direction, bool positive)
			: PartDynamicsByParticle(body, body_part), FluidDataDelegateSimple(body), pos_n_(particles_->pos_n_),
			axis_(axis_direction), periodic_translation_(0), body_buffer_width_(body_buffer_width) 
		{
			body_part->getBodyPartShape()->findBounds(body_part_lower_bound_, body_part_upper_bound_);
			periodic_translation_[axis_] = body_part_upper_bound_[axis_] - body_part_lower_bound_[axis_];
			size_t total_body_buffer_particles = constrained_particles_.size() * body_buffer_width_;
			for (size_t i = 0; i < total_body_buffer_particles; ++i)
			{
				particles_->addABufferParticle();
			}
			particles_->real_particles_bound_ += total_body_buffer_particles;
			body_->allocateConfigurationMemoriesForBodyBuffer();

			checking_bound_ = positive ?
				std::bind(&EmitterInflowInjecting::CheckUpperBound, this, _1, _2)
				: std::bind(&EmitterInflowInjecting::CheckLowerBound, this, _1, _2);
		}
		//=================================================================================================//
		void EmitterInflowInjecting::CheckUpperBound(size_t index_i, Real dt)
		{
			if (pos_n_[index_i][axis_] > body_part_upper_bound_[axis_]) {
				if (body_->number_of_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(body_->number_of_particles_, index_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				/** Periodic bounding. */
				pos_n_[index_i][axis_] -= periodic_translation_[axis_];

			}
		}
		//=================================================================================================//
		void EmitterInflowInjecting::CheckLowerBound(size_t index_i, Real dt)
		{
			if (pos_n_[index_i][axis_] < body_part_lower_bound_[axis_]) {
				if (body_->number_of_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(body_->number_of_particles_, index_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				body_->number_of_particles_ += 1;
				pos_n_[index_i][axis_] += periodic_translation_[axis_];
			}
		}
		//=================================================================================================//
		ViscousAccelerationWallModel::ViscousAccelerationWallModel(SPHBodyComplexRelation* body_complex_relation)
			: ViscousAcceleration(body_complex_relation)
		{
			SPHBody* sph_body = body_complex_relation->sph_body_;
			BaseParticles* base_particles = sph_body->base_particles_;
			//register particle variables defined in this class
			base_particles->registerAVariable(gradient_p_, base_particles->registered_vectors_,
				base_particles->vectors_map_, base_particles->vectors_to_write_, "NearWallPressureGradient", false);
			base_particles->registerAVariable(gradient_vel_, base_particles->registered_matrices_,
				base_particles->matrices_map_, base_particles->matrices_to_write_, "NearWallVelocityGradient", false);

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void ViscousAccelerationWallModel::ComplexInteraction(size_t index_i, Real dt)
		{
			ViscousAcceleration::ComplexInteraction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];
			Real p_i = p_[index_i];

			Vecd acceleration(0), gradient_p(0);
			Matd gradient_vel(0);
			// computing the  outer region values for near wall particles
			if (contact_configuration_.size() != 0)
			{
				Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
					size_t index_j = inner_neighborhood.j_[n];
					Real Vol_j = Vol_[index_j];

					gradient_p += (p_i - p_[index_j]) * nablaW_ij * Vol_j;
					gradient_vel += SimTK::outer(vel_i - vel_n_[index_j], nablaW_ij) * Vol_j;
				}

				/** apply boundary condition for the  outer region values */
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
					Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
						size_t index_j = contact_neighborhood.j_[n];

						gradient_vel += 2.0 * SimTK::outer(vel_i - vel_ave_k[index_j], nablaW_ij) * Vol_k[index_j];
					}
				}
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				StdLargeVec<Vecd>& n_k = *(contact_n_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				// solving inner region momentum balance equation in the tangential direction
				// using outer region values as upper boundary conditions
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					Real dW_ij = contact_neighborhood.dW_ij_[n];
					size_t index_j = contact_neighborhood.j_[n];

					//obtain the unit vector in the tangential direction
					Real height = 0.5 * contact_neighborhood.r_ij_[n];
					Vecd vel_to_wall = (vel_i - vel_ave_k[index_j]);
					Real vel_to_wall_n = SimTK::dot(vel_to_wall, n_k[index_j]);
					Vecd vel_to_wall_t = vel_to_wall - vel_to_wall_n * n_k[index_j];
					Vecd unit_t = vel_to_wall_t / (vel_to_wall_t.norm() + TinyReal);

					Real coefficient_A = SimTK::dot(gradient_p, unit_t);
					Real coefficient_B = SimTK::dot(gradient_vel * vel_to_wall, unit_t);
					Vecd v_derivative_add = -height * (coefficient_A + rho_i * coefficient_B) * unit_t;

					acceleration += 2.0 * v_derivative_add * Vol_k[index_j] * dW_ij / rho_i;
				}
			}

			// save the pressure and velocity gradient
			gradient_p_[index_i] = gradient_p;
			gradient_vel_[index_i] = gradient_vel;

			/** Particle summation. */
			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//