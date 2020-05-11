/**
 * @file 	general_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "general_dynamics.h"

namespace SPH {
//=================================================================================================//
	InitializeATimeStep
		::InitializeATimeStep(SPHBody* body)
		: ParticleDynamicsSimple<SPHBody, BaseParticles>(body)
	{
		initial_value_ = Vecd(0);
	}
//=================================================================================================//
	InitializeATimeStep
		::InitializeATimeStep(SPHBody* body, ExternalForce *external_force)
		: ParticleDynamicsSimple<SPHBody, BaseParticles>(body)
	{
		initial_value_ = external_force->InducedAcceleration();
	}
//=================================================================================================//
	void InitializeATimeStep::SetupDynamics(Real dt)
	{
		particles_->number_of_ghost_particles_ = 0;
	}
//=================================================================================================//
	void InitializeATimeStep::Update(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		base_particle_data_i.dvel_dt_others_ = initial_value_;
	}
//=================================================================================================//
	RandomizePartilePosition::RandomizePartilePosition(SPHBody* body)
		: ParticleDynamicsSimple<SPHBody, BaseParticles>(body)
	{
		particle_spacing_ = body->particle_spacing_;
	}
//=================================================================================================//
	void RandomizePartilePosition::Update(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];

		for (int i = 0; i < base_particle_data_i.pos_n_.size(); ++i)
		{
			base_particle_data_i.pos_n_[i] += dt * (((double)rand() / (RAND_MAX)) - 0.5) * 2.0 * particle_spacing_;
		}
	}
//=================================================================================================//
	BoundingBodyDomain
		::BoundingBodyDomain(SPHBody* body)
		: ParticleDynamicsByCells(body)
	{
		body_->BodyBounds(body_lower_bound_, body_upper_bound_);
		SetCellBounds();
	}
//=================================================================================================//
	void BoundingBodyDomain::SetCellBounds()
	{
		Vecd rltpos = body_lower_bound_ - mesh_lower_bound_;
		for (int i = 0; i < rltpos.size(); ++i) {
			body_lower_bound_cell_[i] = (size_t)floor(rltpos[i] / cell_spacing_);
			/**< floor in c++ Rounds x downward, returning the largest integral value 
			 		that is not greater than x. */
		}

		rltpos = body_upper_bound_ - mesh_lower_bound_;
		for (int i = 0; i < rltpos.size(); ++i) {
			body_upper_bound_cell_[i] = (size_t)floor(rltpos[i] / cell_spacing_);
		}
	}
	//=================================================================================================//
	void PeriodicBoundingInAxisDirection::setPeriodicTranslation()
	{
		periodic_translation_[axis_] = body_upper_bound_[axis_] - body_lower_bound_[axis_];
		if (periodic_translation_.norm() < body_->particle_spacing_) {
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			std::cout << "\n Periodic bounding failure: bounds not defined!" << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void PeriodicBoundingInAxisDirection::CheckLowerBound(size_t index_particle_i, Vecd& pnt, Real dt)
	{
		BaseParticleData& base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[axis_] < body_lower_bound_[axis_])
			base_particle_data_i.pos_n_[axis_] += periodic_translation_[axis_];
	}
//=================================================================================================//
	void PeriodicBoundingInAxisDirection::CheckUpperBound(size_t index_particle_i, Vecd& pnt, Real dt)
	{
		BaseParticleData& base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[axis_] > body_upper_bound_[axis_])
			base_particle_data_i.pos_n_[axis_] -= periodic_translation_[axis_];
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::CheckLowerBound(size_t index_particle_i, Vecd& pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[axis_] > body_lower_bound_[axis_]
			&& particle_position[axis_] < (body_lower_bound_[axis_] + cell_spacing_))
		{
			Vecd translated_position = particle_position + periodic_translation_;
			/** insert ghost particle to cell linked list */
			mesh_cell_linked_list_->InsertACellLinkedListEntry(index_particle_i, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::CheckUpperBound(size_t index_particle_i, Vecd& pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[axis_] < body_upper_bound_[axis_]
			&& particle_position[axis_] > (body_upper_bound_[axis_] - cell_spacing_))
		{
			Vecd translated_position = particle_position - periodic_translation_;
			/** insert ghost particle to cell linked list */
			mesh_cell_linked_list_->InsertACellLinkedListEntry(index_particle_i, translated_position);
		}
	}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection::Bounding
		::Bounding(CellVector& bound_cells, SPHBody* body, int axis_direction, bool positive)
		: BoundingInAxisDirection(body, axis_direction),
		bound_cells_(bound_cells)
	{
		checking_bound_ = positive ? 
			std::bind(&MirrorBoundaryConditionInAxisDirection::Bounding::checkUpperBound, this, _1, _2) 
			: std::bind(&MirrorBoundaryConditionInAxisDirection::Bounding::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::CreatingGhostParticles(IndexVector& ghost_particles, 
			CellVector& bound_cells, SPHBody* body, int axis_direction, bool positive)
		: Bounding(bound_cells, body, axis_direction, positive), ghost_particles_(ghost_particles) {}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::UpdatingGhostStates(IndexVector& ghost_particles,
			SPHBody* body, int axis_direction, bool positive)
		: BoundingInAxisDirection(body, axis_direction), ghost_particles_(ghost_particles) 
	{
		checking_bound_update_ = positive ?
			std::bind(&MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates::updateForUpperBound, this, _1, _2)
			: std::bind(&MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates::updateForLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::Bounding::checkLowerBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[axis_] < body_lower_bound_[axis_]) {
			particles_->mirrorInAxisDirection(index_particle_i, body_lower_bound_, axis_);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::Bounding
		::checkUpperBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[axis_] > body_upper_bound_[axis_]) {
			particles_->mirrorInAxisDirection(index_particle_i, body_upper_bound_, axis_);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::checkLowerBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

		Vecd particle_position = base_particle_data_i.pos_n_;
		if (particle_position[axis_] > body_lower_bound_[axis_]
			&& particle_position[axis_] < (body_lower_bound_[axis_] + cell_spacing_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_particle_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			particles_->mirrorInAxisDirection(expected_particle_index, body_lower_bound_, axis_);
			Vecd translated_position = particles_->base_particle_data_[expected_particle_index].pos_n_;
			/** insert ghost particle to cell linked list */
			mesh_cell_linked_list_->InsertACellLinkedListEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::checkUpperBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

		Vecd particle_position = base_particle_data_i.pos_n_;
		if (particle_position[axis_] < body_upper_bound_[axis_]
			&& particle_position[axis_] > (body_upper_bound_[axis_] - cell_spacing_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_particle_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			particles_->mirrorInAxisDirection(expected_particle_index, body_upper_bound_, axis_);
			Vecd translated_position = particles_->base_particle_data_[expected_particle_index].pos_n_;
			/** insert ghost particle to cell linked list */
			mesh_cell_linked_list_->InsertACellLinkedListEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::updateForLowerBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		particles_->UpdateFromAnotherParticle(index_particle_i, base_particle_data_i.particle_id_);
		particles_->mirrorInAxisDirection(index_particle_i, body_lower_bound_, axis_);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::updateForUpperBound(size_t index_particle_i, Real dt)
	{
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		particles_->UpdateFromAnotherParticle(index_particle_i, base_particle_data_i.particle_id_);
		particles_->mirrorInAxisDirection(index_particle_i, body_upper_bound_, axis_);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::exec(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_.size(); ++i) {
				checking_bound_update_(ghost_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, ghost_particles_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					checking_bound_update_(ghost_particles_[i], dt);
				}
			}, ap);
	}
	//=================================================================================================//
	VelocityBoundCheck::
		VelocityBoundCheck(SPHBody* body, Real velocity_bound)
		: ParticleDynamicsReduce<bool, ReduceOR, SPHBody, BaseParticles>(body),
		velocity_bound_(velocity_bound)
	{
		initial_reference_ = false;
	}
//=================================================================================================//
	bool VelocityBoundCheck::ReduceFunction(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i
			= particles_->base_particle_data_[index_particle_i];

		return base_particle_data_i.vel_n_.norm() > velocity_bound_;
	}
//=================================================================================================//
	Real UpperFrontInXDirection::ReduceFunction(size_t index_particle_i, Real dt)
	{
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

		return base_particle_data_i.pos_n_[0];
	}
//=================================================================================================//
}
//=================================================================================================//