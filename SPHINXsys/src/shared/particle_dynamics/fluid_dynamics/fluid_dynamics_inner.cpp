#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FluidInitialCondition::
			FluidInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body), FluidDataSimple(sph_body),
			  pos_(particles_->pos_), vel_(particles_->vel_) {}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseBodyRelationInner &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), mass_(particles_->mass_),
			  rho_sum_(particles_->rho_sum_),
			  W0_(sph_body_.sph_adaptation_->getKernel()->W0(zero_vec)),
			  rho0_(particles_->rho0_), inv_sigma0_(1.0 / particles_->sigma0_) {}
		//=================================================================================================//
		void DensitySummationInner::interaction(size_t index_i, Real dt)
		{
			/** Inner interaction. */
			Real sigma = W0_;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
		}
		//=================================================================================================//
		void DensitySummationInner::update(size_t index_i, Real dt)
		{
			rho_[index_i] = ReinitializedDensity(rho_sum_[index_i], rho0_, rho_[index_i]);
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		BaseViscousAccelerationInner::BaseViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			  vel_(particles_->vel_), acc_prior_(particles_->acc_prior_),
			  mu_(material_->ReferenceViscosity()),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd &vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				// viscous force
				vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}

			acc_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		void AngularConservativeViscousAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd &vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
				 * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
				Real v_r_ij = (vel_i - vel_[index_j]).dot(r_ij * e_ij);
				Real eta_ij = 8.0 * mu_ * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * Vol_[index_j] / rho_i * inner_neighborhood.dW_ij_[n] * e_ij;
			}

			acc_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseBodyRelationInner &inner_relation, Real coefficient)
			: LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), pos_(particles_->pos_),
			  surface_indicator_(particles_->surface_indicator_), p_background_(0),
			  coefficient_(coefficient) {}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::setupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->ReferenceDensity();
			p_background_ = coefficient_ * density * speed_max * speed_max;
		}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];

			Vecd acceleration_trans = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				// acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_ * Vol_[index_j] * nablaW_ij / rho_i;
			}

			if (surface_indicator_[index_i] == 0)
				pos_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  FluidDataSimple(sph_body), rho_(particles_->rho_),
			  p_(particles_->p_), vel_(particles_->vel_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return material_->getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSizeForImplicitViscosity::
			AdvectionTimeStepSizeForImplicitViscosity(SPHBody &sph_body, Real U_max)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max * U_max),
			  FluidDataSimple(sph_body), vel_(particles_->vel_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::reduce(size_t index_i, Real dt)
		{
			return vel_[index_i].squaredNorm();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::outputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			particles_->speed_max_ = speed_max;
			return 0.25 * smoothing_length_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(SPHBody &sph_body, Real U_max)
			: AdvectionTimeStepSizeForImplicitViscosity(sph_body, U_max)
		{
			Real viscous_speed = material_->ReferenceViscosity() / material_->ReferenceDensity() / smoothing_length_;
			reference_ = SMAX(viscous_speed * viscous_speed, reference_);
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return AdvectionTimeStepSizeForImplicitViscosity::reduce(index_i, dt);
		}
		//=================================================================================================//
		VorticityInner::VorticityInner(BaseBodyRelationInner &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), vel_(particles_->vel_)
		{
			particles_->registerVariable(vorticity_, "VorticityInner");
			particles_->addVariableToWrite<AngularVecd>("VorticityInner");
		}
		//=================================================================================================//
		void VorticityInner::interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_i = vel_[index_i];

			AngularVecd vorticity(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd vel_diff = vel_i - vel_[index_j];
				vorticity += getCrossProduct(vel_diff, inner_neighborhood.e_ij_[n]) * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}

			vorticity_[index_i] = vorticity;
		}
		//=================================================================================================//
		BaseRelaxation::BaseRelaxation(BaseBodyRelationInner &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), mass_(particles_->mass_), rho_(particles_->rho_),
			  p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  acc_(particles_->acc_),
			  acc_prior_(particles_->acc_prior_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BasePressureRelaxation::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
			p_[index_i] = material_->getPressure(rho_[index_i]);
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BasePressureRelaxation::update(size_t index_i, Real dt)
		{
			vel_[index_i] += acc_[index_i] * dt;
		}
		//=================================================================================================//
		Vecd BasePressureRelaxation::computeNonConservativeAcceleration(size_t index_i)
		{
			Real rho_i = rho_[index_i];
			Real p_i = p_[index_i];
			Vecd acceleration = acc_prior_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				Real rho_j = rho_[index_j];
				Real p_j = p_[index_j];

				Real p_star = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
				acceleration += (p_i - p_star) * Vol_[index_j] * dW_ij * e_ij / rho_i;
			}
			return acceleration;
		}
		//=================================================================================================//
		BaseDensityRelaxation::
			BaseDensityRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BaseDensityRelaxation::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BaseDensityRelaxation::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		PressureRelaxationInnerOldroyd_B ::
			PressureRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation)
			: PressureRelaxationDissipativeRiemannInner(inner_relation),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_.base_particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_.base_particles_)->dtau_dt_) {}
		//=================================================================================================//
		void PressureRelaxationInnerOldroyd_B::initialization(size_t index_i, Real dt)
		{
			PressureRelaxationDissipativeRiemannInner::initialization(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationInnerOldroyd_B::interaction(size_t index_i, Real dt)
		{
			PressureRelaxationDissipativeRiemannInner::interaction(index_i, dt);

			Real rho_i = rho_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration = Vecd::Zero();
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				// elastic force
				acceleration += (tau_i + tau_[index_j]) * nablaW_ij * Vol_[index_j] / rho_i;
			}

			acc_[index_i] += acceleration;
		}
		//=================================================================================================//
		DensityRelaxationInnerOldroyd_B::
			DensityRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation)
			: DensityRelaxationDissipativeRiemannInner(inner_relation),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_.base_particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_.base_particles_)->dtau_dt_)
		{
			Oldroyd_B_Fluid *oldroy_b_fluid = DynamicCast<Oldroyd_B_Fluid>(this, sph_body_.base_material_);
			mu_p_ = oldroy_b_fluid->ReferencePolymericViscosity();
			lambda_ = oldroy_b_fluid->getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void DensityRelaxationInnerOldroyd_B::interaction(size_t index_i, Real dt)
		{
			DensityRelaxationDissipativeRiemannInner::interaction(index_i, dt);

			Vecd vel_i = vel_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate = Matd::Zero();
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				Matd velocity_gradient = - Vol_[index_j] * (vel_i - vel_[index_j]) * nablaW_ij.transpose();
				stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
							   (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
			}

			dtau_dt_[index_i] = stress_rate;
		}
		//=================================================================================================//
		void DensityRelaxationInnerOldroyd_B::update(size_t index_i, Real dt)
		{
			DensityRelaxationDissipativeRiemannInner::update(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
	}
	//=====================================================================================================//
}
//=========================================================================================================//