#include "general_dynamics.h"
namespace SPH {
	//===================================================================//
	InitializeOtherAccelerations
		::InitializeOtherAccelerations(SPHBody* body)
		: ParticleDynamicsSimple<SPHBody, Particles>(body)
	{
		initial_value_ = Vecd(0);
	}
	//===================================================================//
	InitializeOtherAccelerations
		::InitializeOtherAccelerations(SPHBody* body, ExternalForce *external_force)
		: ParticleDynamicsSimple<SPHBody, Particles>(body)
	{
		initial_value_ = external_force->InducedAcceleration();
	}
	//===================================================================//
	void InitializeOtherAccelerations::Update(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		base_particle_data_i.dvel_dt_others_ = initial_value_;
	}
	//===================================================================//
	RandomizePartilePosition::RandomizePartilePosition(SPHBody* body)
		: ParticleDynamicsSimple<SPHBody, Particles>(body)
	{
		particle_spacing_ = body->particle_spacing_;
	}
	//===================================================================//
	void RandomizePartilePosition::Update(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];

		for (int i = 0; i < base_particle_data_i.pos_n_.size(); ++i)
		{
			base_particle_data_i.pos_n_[i] += 0.25 * (((double)rand() / (RAND_MAX)) - 1.0) * particle_spacing_;
		}
	}
	//===================================================================//
	BoundingBodyDomain
		::BoundingBodyDomain(SPHBody* body)
		: ParticleDynamicsByCells(body)
	{
		body_->BodyBounds(body_lower_bound_, body_upper_bound_);
		SetCellBounds();
	}
	//===================================================================//
	void BoundingBodyDomain::SetCellBounds()
	{
		Vecd rltpos = body_lower_bound_ - mesh_lower_bound_;
		for (size_t i = 0; i < rltpos.size(); ++i) {
			body_lower_bound_cell_[i] = floor(rltpos[i] / cell_spacing_);
			// floor in c++ Rounds x downward, returning the largest integral value that is not greater than x.
		}

		rltpos = body_upper_bound_ - mesh_lower_bound_;
		for (size_t i = 0; i < rltpos.size(); ++i) {
			body_upper_bound_cell_[i] = floor(rltpos[i] / cell_spacing_);
		}
	}
	//===================================================================//
	PeriodicBoundingInXDirection
		::PeriodicBoundingInXDirection(SPHBody* body)
		: BoundingInXDirection(body), periodic_translation_(0)
	{
		periodic_translation_[0] = body_upper_bound_[0] - body_lower_bound_[0];
		if (periodic_translation_.norm() < body->particle_spacing_) {
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			std::cout << "\n Periodic bounding failure: bounds not defined!" << std::endl;
			exit(1);
		}
	}
	//===================================================================//
	void PeriodicBoundingInXDirection
		::CheckLowerBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[0] < body_lower_bound_[0])
			base_particle_data_i.pos_n_[0] += periodic_translation_[0];
	}
	//===================================================================//
	void PeriodicBoundingInXDirection
		::CheckUpperBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[0] > body_upper_bound_[0])
			base_particle_data_i.pos_n_[0] -= periodic_translation_[0];
	}
	//===================================================================//
	PeriodicConditionInXDirection::
		PeriodicConditionInXDirection(SPHBody* body)
		: PeriodicBoundingInXDirection(body)
	{

	}
	//===================================================================//
	VelocityBoundCheck::
		VelocityBoundCheck(SPHBody* body, Real velocity_bound)
		: ParticleDynamicsReduce<bool, ReduceOR, SPHBody, Particles>(body),
		velocity_bound_(velocity_bound)
	{
		initial_reference_ = false;
	}
	//===================================================================//
	bool VelocityBoundCheck::ReduceFunction(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];

		return base_particle_data_i.vel_n_.norm() > velocity_bound_;
	}
	//===================================================================//
	Real UpperBoundInXDirection::ReduceFunction(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

		return base_particle_data_i.pos_n_[0];
	}
	//===================================================================//
}