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
		BaseDensitySummationInner::BaseDensitySummationInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  rho_(particles_->rho_), rho_sum_(particles_->rho_sum_), mass_(particles_->mass_),
			  rho0_(sph_body_.base_material_->ReferenceDensity()) {}
		//=================================================================================================//
		void BaseDensitySummationInner::update(size_t index_i, Real dt)
		{
			rho_[index_i] = rho_sum_[index_i];
		}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseInnerRelation &inner_relation)
			: BaseDensitySummationInner(inner_relation),
			  W0_(sph_body_.sph_adaptation_->getKernel()->W0( ZeroVecd )),
			  inv_sigma0_(1.0 / sph_body_.sph_adaptation_->ReferenceNumberDensity()) {}
		//=================================================================================================//
		DensitySummationInnerAdaptive::
			DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation)
			: BaseDensitySummationInner(inner_relation),
			  sph_adaptation_(*sph_body_.sph_adaptation_),
			  kernel_(*sph_adaptation_.getKernel()),
			  h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
		//=================================================================================================//
		BaseViscousAccelerationInner::BaseViscousAccelerationInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  rho_(particles_->rho_), vel_(particles_->vel_), acc_prior_(particles_->acc_prior_),
			  mu_(particles_->fluid_.ReferenceViscosity()),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseInnerRelation &inner_relation, Real coefficient)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			  smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			  coefficient_(coefficient) {}
		//=================================================================================================//
		TransportVelocityCorrectionInnerAdaptive::
			TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation &inner_relation, Real coefficient)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  sph_adaptation_(*sph_body_.sph_adaptation_),
			  pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			  smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			  coefficient_(coefficient) {}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  FluidDataSimple(sph_body), fluid_(particles_->fluid_), rho_(particles_->rho_),
			  p_(particles_->p_), vel_(particles_->vel_),
			  smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			  acousticCFL_(acousticCFL) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return acousticCFL_ * smoothing_length_min_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSizeForImplicitViscosity::
			AdvectionTimeStepSizeForImplicitViscosity(SPHBody &sph_body, Real U_max, Real advectionCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max * U_max),
			  FluidDataSimple(sph_body), vel_(particles_->vel_),
			  smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			  advectionCFL_(advectionCFL) {}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::reduce(size_t index_i, Real dt)
		{
			return vel_[index_i].squaredNorm();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::outputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			return advectionCFL_ * smoothing_length_min_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL)
			: AdvectionTimeStepSizeForImplicitViscosity(sph_body, U_max, advectionCFL),
			  fluid_(particles_->fluid_)
		{
			Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
			reference_ = SMAX(viscous_speed * viscous_speed, reference_);
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return AdvectionTimeStepSizeForImplicitViscosity::reduce(index_i, dt);
		}
		//=================================================================================================//
		VorticityInner::VorticityInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  vel_(particles_->vel_)
		{
			particles_->registerVariable(vorticity_, "VorticityInner");
			particles_->addVariableToWrite<AngularVecd>("VorticityInner");
		}
		//=================================================================================================//
		BaseIntegration::BaseIntegration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			  fluid_(particles_->fluid_), rho_(particles_->rho_),
			  p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  acc_(particles_->acc_), acc_prior_(particles_->acc_prior_) {}
		//=================================================================================================//
		Oldroyd_BIntegration1stHalf ::
			Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation)
			: Integration1stHalfDissipativeRiemann(inner_relation),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, particles_)->dtau_dt_) {}
		//=================================================================================================//
		void Oldroyd_BIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			Integration1stHalfDissipativeRiemann::initialization(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		Oldroyd_BIntegration2ndHalf::
			Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation)
			: Integration2ndHalfDissipativeRiemann(inner_relation),
			  oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->fluid_)),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, particles_)->dtau_dt_)
		{
			mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
			lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void Oldroyd_BIntegration2ndHalf::update(size_t index_i, Real dt)
		{
			Integration2ndHalfDissipativeRiemann::update(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
	}
	//=====================================================================================================//
}
//=========================================================================================================//