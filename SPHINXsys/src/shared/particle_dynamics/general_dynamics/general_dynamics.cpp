/**
 * @file 	general_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "general_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	TimeStepInitialization::TimeStepInitialization(SPHBody &sph_body)
		: ParticleDynamicsSimple(sph_body), GeneralDataDelegateSimple(sph_body),
		  pos_n_(particles_->pos_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
		  gravity_(gravity_ptr_keeper_.createPtr<Gravity>(Vecd(0))) {}
	//=================================================================================================//
	TimeStepInitialization ::TimeStepInitialization(SPHBody &sph_body, Gravity &gravity)
		: ParticleDynamicsSimple(sph_body), GeneralDataDelegateSimple(sph_body),
		  pos_n_(particles_->pos_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
		  gravity_(&gravity) {}
	//=================================================================================================//
	void TimeStepInitialization::setupDynamics(Real dt)
	{
		particles_->total_ghost_particles_ = 0;
	}
	//=================================================================================================//
	void TimeStepInitialization::Update(size_t index_i, Real dt)
	{
		dvel_dt_prior_[index_i] = gravity_->InducedAcceleration(pos_n_[index_i]);
	}
	//=================================================================================================//
	RandomizePartilePosition::RandomizePartilePosition(SPHBody &sph_body)
		: ParticleDynamicsSimple(sph_body), DataDelegateSimple<SPHBody, BaseParticles>(sph_body),
		  pos_n_(particles_->pos_n_), randomize_scale_(sph_body.sph_adaptation_->MinimumSpacing()) {}
	//=================================================================================================//
	void RandomizePartilePosition::Update(size_t index_i, Real dt)
	{
		Vecd &pos_n_i = pos_n_[index_i];
		for (int k = 0; k < pos_n_i.size(); ++k)
		{
			pos_n_i[k] += dt * (((double)rand() / (RAND_MAX)) - 0.5) * 2.0 * randomize_scale_;
		}
	}
	//=================================================================================================//
	VelocityBoundCheck::
		VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound)
		: ParticleDynamicsReduce<bool, ReduceOR>(sph_body),
		  GeneralDataDelegateSimple(sph_body),
		  vel_n_(particles_->vel_n_), velocity_bound_(velocity_bound)
	{
		initial_reference_ = false;
	}
	//=================================================================================================//
	bool VelocityBoundCheck::ReduceFunction(size_t index_i, Real dt)
	{
		return vel_n_[index_i].norm() > velocity_bound_;
	}
	//=================================================================================================//
	UpperFrontInXDirection::
		UpperFrontInXDirection(SPHBody &sph_body) : ParticleDynamicsReduce<Real, ReduceMax>(sph_body),
													GeneralDataDelegateSimple(sph_body),
													pos_n_(particles_->pos_n_)
	{
		quantity_name_ = "UpperFrontInXDirection";
		initial_reference_ = 0.0;
	}
	//=================================================================================================//
	Real UpperFrontInXDirection::ReduceFunction(size_t index_i, Real dt)
	{
		return pos_n_[index_i][0];
	}
	//=================================================================================================//
	MaximumSpeed::
		MaximumSpeed(SPHBody &sph_body) : ParticleDynamicsReduce<Real, ReduceMax>(sph_body),
										  GeneralDataDelegateSimple(sph_body),
										  vel_n_(particles_->vel_n_)
	{
		quantity_name_ = "MaximumSpeed";
		initial_reference_ = 0.0;
	}
	//=================================================================================================//
	Real MaximumSpeed::ReduceFunction(size_t index_i, Real dt)
	{
		return vel_n_[index_i].norm();
	}
	//=================================================================================================//
	BodyLowerBound::BodyLowerBound(SPHBody &sph_body)
		: ParticleDynamicsReduce<Vecd, ReduceLowerBound>(sph_body),
		  GeneralDataDelegateSimple(sph_body),
		  pos_n_(particles_->pos_n_)
	{
		constexpr double max_real_number = (std::numeric_limits<double>::max)();
		initial_reference_ = Vecd(max_real_number);
	}
	//=================================================================================================//
	Vecd BodyLowerBound::ReduceFunction(size_t index_i, Real dt)
	{
		return pos_n_[index_i];
	}
	//=================================================================================================//
	BodyUpperBound::
		BodyUpperBound(SPHBody &sph_body) : ParticleDynamicsReduce<Vecd, ReduceUpperBound>(sph_body),
											GeneralDataDelegateSimple(sph_body),
											pos_n_(particles_->pos_n_)
	{
		constexpr double min_real_number = (std::numeric_limits<double>::min)();
		initial_reference_ = Vecd(min_real_number);
	}
	//=================================================================================================//
	Vecd BodyUpperBound::ReduceFunction(size_t index_i, Real dt)
	{
		return pos_n_[index_i];
	}
	//=================================================================================================//
	TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody &sph_body)
		: ParticleDynamicsReduce<Real, ReduceSum<Real>>(sph_body),
		  GeneralDataDelegateSimple(sph_body), mass_(particles_->mass_),
		  vel_n_(particles_->vel_n_), pos_n_(particles_->pos_n_),
		  gravity_(gravity_ptr_keeper_.createPtr<Gravity>(Vecd(0)))
	{
		quantity_name_ = "TotalMechanicalEnergy"; // TODO: this need to ben changed as "TotalKineticEnergy"
		initial_reference_ = 0.0;
	}
	//=================================================================================================//
	TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity)
		: ParticleDynamicsReduce<Real, ReduceSum<Real>>(sph_body),
		  GeneralDataDelegateSimple(sph_body), mass_(particles_->mass_),
		  vel_n_(particles_->vel_n_), pos_n_(particles_->pos_n_),
		  gravity_(&gravity)
	{
		quantity_name_ = "TotalMechanicalEnergy";
		initial_reference_ = 0.0;
	}
	//=================================================================================================//
	Real TotalMechanicalEnergy::ReduceFunction(size_t index_i, Real dt)
	{
		return 0.5 * mass_[index_i] * vel_n_[index_i].normSqr() + mass_[index_i] * gravity_->getPotential(pos_n_[index_i]);
	}
	//=================================================================================================//
}
//=================================================================================================//