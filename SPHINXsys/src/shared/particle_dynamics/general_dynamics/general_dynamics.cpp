/**
 * @file 	general_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "general_dynamics.h"

namespace SPH {
	//=================================================================================================//
	TimeStepInitialization
		::TimeStepInitialization(SPHBody* body, Gravity* gravity)
		: ParticleDynamicsSimple(body), GeneralDataDelegateSimple(body),
		pos_n_(particles_->pos_n_), dvel_dt_prior_(particles_->dvel_dt_prior_),
		gravity_(gravity)
	{
	}
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
	RandomizePartilePosition::RandomizePartilePosition(SPHBody* body)
		: ParticleDynamicsSimple(body), DataDelegateSimple<SPHBody, BaseParticles>(body),
		pos_n_(particles_->pos_n_)
	{
		randomize_scale_ = body->particle_adaptation_->MinimumSpacing();
	}
	//=================================================================================================//
	void RandomizePartilePosition::Update(size_t index_i, Real dt)
	{
		Vecd& pos_n_i = pos_n_[index_i];
		for (int k = 0; k < pos_n_i.size(); ++k)
		{
			pos_n_i[k] += dt * (((double)rand() / (RAND_MAX)) - 0.5) * 2.0 * randomize_scale_;
		}
	}
	//=================================================================================================//
	BoundingInAxisDirection::BoundingInAxisDirection(RealBody* real_body, int axis_direction) :
		ParticleDynamics<void>(real_body), DataDelegateSimple<SPHBody, BaseParticles>(real_body),
		axis_(axis_direction), body_domain_bounds_(real_body->getBodyDomainBounds()),
		pos_n_(particles_->pos_n_),
		cell_linked_list_(real_body->cell_linked_list_)
	{
		Real h_ratio_min = 1.0 / particle_adaptation_->MaximumSpacingRatio();
		cut_off_radius_max_ = particle_adaptation_->getKernel()->CutOffRadius(h_ratio_min);
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::
		setPeriodicTranslation(BoundingBox& body_domain_bounds, int axis_direction)
	{
		periodic_translation_[axis_direction] = 
			body_domain_bounds.second[axis_direction] - body_domain_bounds.first[axis_direction];
	}
	//=================================================================================================//
	PeriodicConditionInAxisDirection::
		PeriodicConditionInAxisDirection(RealBody* real_body, int axis_direction) :
		periodic_translation_(0.0)
	{
		BoundingBox body_domain_bounds = real_body->getBodyDomainBounds();
		setPeriodicTranslation(body_domain_bounds, axis_direction);
		bound_cells_.resize(2);
		BaseCellLinkedList* cell_linked_list = real_body->cell_linked_list_;
		cell_linked_list->tagBodyDomainBoundingCells(bound_cells_, body_domain_bounds, axis_direction);
		if (periodic_translation_.norm() < real_body->particle_adaptation_->ReferenceSpacing())
		{
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			std::cout << "\n Periodic bounding failure: bounds not defined!" << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::checkLowerBound(size_t index_i, Real dt)
	{
		if (pos_n_[index_i][axis_] < body_domain_bounds_.first[axis_])
			pos_n_[index_i][axis_] += periodic_translation_[axis_];
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_n_[index_i][axis_] > body_domain_bounds_.second[axis_])
			pos_n_[index_i][axis_] -= periodic_translation_[axis_];
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::exec(Real dt)
	{
		setupDynamics(dt);

		//check lower bound
		CellLists& lower_bound_cells = bound_cells_[0];
		for (size_t i = 0; i != lower_bound_cells.size(); ++i) 
		{
			IndexVector& particle_indexes = lower_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num) checkLowerBound(particle_indexes[num], dt);
		}

		//check upper bound
		CellLists& upper_bound_cells = bound_cells_[1];
		for (size_t i = 0; i != upper_bound_cells.size(); ++i) 
		{
			IndexVector& particle_indexes = upper_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num) checkUpperBound(particle_indexes[num], dt);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		//check lower bound
		CellLists& lower_bound_cells = bound_cells_[0];
		parallel_for(blocked_range<size_t>(0, lower_bound_cells.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					IndexVector& particle_indexes = lower_bound_cells[i]->real_particle_indexes_;
					for (size_t num = 0; num < particle_indexes.size(); ++num) checkLowerBound(particle_indexes[num], dt);
				}
			}, ap);

		//check upper bound
		CellLists& upper_bound_cells = bound_cells_[1];
		parallel_for(blocked_range<size_t>(0, upper_bound_cells.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					IndexVector& particle_indexes = upper_bound_cells[i]->real_particle_indexes_;
					for (size_t num = 0; num < particle_indexes.size(); ++num) checkUpperBound(particle_indexes[num], dt);
				}
			}, ap);
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicCondition::exec(Real dt)
	{
		setupDynamics(dt);

		//check lower bound
		CellLists& lower_bound_cells = bound_cells_[0];
		for (size_t i = 0; i != lower_bound_cells.size(); ++i) {
			ListDataVector& cell_list_data = lower_bound_cells[i]->cell_list_data_;
			for (size_t num = 0; num < cell_list_data.size(); ++num) checkLowerBound(cell_list_data[num], dt);
		}

		//check upper bound
		CellLists& upper_bound_cells = bound_cells_[1];
		for (size_t i = 0; i != upper_bound_cells.size(); ++i) {
			ListDataVector& cell_list_data = upper_bound_cells[i]->cell_list_data_;
			for (size_t num = 0; num < cell_list_data.size(); ++num) checkUpperBound(cell_list_data[num], dt);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingCellLinkedList::
		PeriodicCellLinkedList::checkUpperBound(ListData& list_data, Real dt)
	{
		Vecd particle_position = list_data.second;
		if (particle_position[axis_] < body_domain_bounds_.second[axis_]
			&& particle_position[axis_] > (body_domain_bounds_.second[axis_] - cut_off_radius_max_))
		{
			Vecd translated_position = particle_position - periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(list_data.first, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingCellLinkedList::
		PeriodicCellLinkedList::checkLowerBound(ListData& list_data, Real dt)
	{
		Vecd particle_position = list_data.second;
		if (particle_position[axis_] > body_domain_bounds_.first[axis_]
			&& particle_position[axis_] < (body_domain_bounds_.first[axis_] + cut_off_radius_max_))
		{
			Vecd translated_position = particle_position + periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(list_data.first, translated_position);
		}
	}
	//=================================================================================================//
	OpenBoundaryConditionInAxisDirection::
		OpenBoundaryConditionInAxisDirection(RealBody* real_body, int axis_direction, bool positive):
		particle_type_transfer(this->bound_cells_, real_body, axis_direction, positive)
	{
		BoundingBox body_domain_bounds = real_body->getBodyDomainBounds();
		bound_cells_.resize(2);
		BaseCellLinkedList* cell_linked_list = real_body->cell_linked_list_;
		cell_linked_list->tagBodyDomainBoundingCells(bound_cells_, body_domain_bounds, axis_direction);
	}
	//=================================================================================================//
	void OpenBoundaryConditionInAxisDirection ::
		ParticleTypeTransfer::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] < body_domain_bounds_.first[axis_])
			particles_->switchToBufferParticle(index_i);
	}
	//=================================================================================================//
	void OpenBoundaryConditionInAxisDirection ::
		ParticleTypeTransfer::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] > body_domain_bounds_.second[axis_])
			particles_->switchToBufferParticle(index_i);
	}
	//=================================================================================================//
	void OpenBoundaryConditionInAxisDirection::ParticleTypeTransfer::exec(Real dt)
	{
		setupDynamics(dt);

		//check lower bound
		CellLists& lower_bound_cells = bound_cells_[0];
		for (size_t i = 0; i != lower_bound_cells.size(); ++i)
		{
			IndexVector& particle_indexes = lower_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num) checking_bound_(particle_indexes[num], dt);
		}

		//check upper bound
		CellLists& upper_bound_cells = bound_cells_[1];
		for (size_t i = 0; i != upper_bound_cells.size(); ++i)
		{
			IndexVector& particle_indexes = upper_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num) checking_bound_(particle_indexes[num], dt);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		CreatPeriodicGhostParticles::setupDynamics(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_.size(); ++i) ghost_particles_[i].clear();
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		CreatPeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] > body_domain_bounds_.first[axis_]
			&& particle_position[axis_] < (body_domain_bounds_.first[axis_] + cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[0].push_back(expected_particle_index);
			Vecd translated_position = particle_position + periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		CreatPeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] < body_domain_bounds_.second[axis_]
			&& particle_position[axis_] > (body_domain_bounds_.second[axis_] - cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[1].push_back(expected_particle_index);
			Vecd translated_position = particle_position - periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		UpdatePeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_n_[index_i] += periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		UpdatePeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_n_[index_i] -= periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		UpdatePeriodicGhostParticles::exec(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_[0].size(); ++i) 
		{
			checkLowerBound(ghost_particles_[0][i], dt);
		}
		for (size_t i = 0; i != ghost_particles_[1].size(); ++i)
		{
			checkUpperBound(ghost_particles_[1][i], dt);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirectionUsingGhostParticles::
		UpdatePeriodicGhostParticles::parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, ghost_particles_[0].size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					checkLowerBound(ghost_particles_[0][i], dt);
				}
			}, ap);
	
		parallel_for(blocked_range<size_t>(0, ghost_particles_[1].size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					checkUpperBound(ghost_particles_[1][i], dt);
				}
			}, ap);
	}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection::MirrorBounding
		::MirrorBounding(CellLists& bound_cells, RealBody* real_body, int axis_direction, bool positive)
		: BoundingInAxisDirection(real_body, axis_direction),
		bound_cells_(bound_cells), vel_n_(particles_->vel_n_)
	{
		checking_bound_ = positive ?
			std::bind(&MirrorBoundaryConditionInAxisDirection::MirrorBounding::checkUpperBound, this, _1, _2)
			: std::bind(&MirrorBoundaryConditionInAxisDirection::MirrorBounding::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::CreatingGhostParticles(IndexVector& ghost_particles,
			CellLists& bound_cells, RealBody* real_body, int axis_direction, bool positive)
		: MirrorBounding(bound_cells, real_body, axis_direction, positive), ghost_particles_(ghost_particles) {}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::UpdatingGhostStates(IndexVector& ghost_particles, CellLists& bound_cells,
			RealBody* real_body, int axis_direction, bool positive)
		: MirrorBounding(bound_cells, real_body, axis_direction, positive), ghost_particles_(ghost_particles)
	{
		checking_bound_update_ = positive ?
			std::bind(&MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates::checkUpperBound, this, _1, _2)
			: std::bind(&MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::MirrorBounding::checkLowerBound(size_t index_i, Real dt)
	{
		if (pos_n_[index_i][axis_] < body_domain_bounds_.first[axis_]) {
			mirrorInAxisDirection(index_i, body_domain_bounds_.first, axis_);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::MirrorBounding
		::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_n_[index_i][axis_] > body_domain_bounds_.second[axis_]) {
			mirrorInAxisDirection(index_i, body_domain_bounds_.second, axis_);
		}
	}
	//=================================================================================================//
	void  MirrorBoundaryConditionInAxisDirection::MirrorBounding::
		mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction)
	{
		pos_n_[particle_index_i][axis_direction]
			= 2.0 * body_bound[axis_direction] - pos_n_[particle_index_i][axis_direction];
		vel_n_[particle_index_i][axis_direction] *= -1.0;
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::MirrorBounding::exec(Real dt)
	{
		setupDynamics(dt);
		for (size_t i = 0; i != bound_cells_.size(); ++i) {
			ListDataVector& list_data = bound_cells_[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num) checking_bound_(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::MirrorBounding ::parallel_exec(Real dt)
	{
		setupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					ListDataVector& list_data = bound_cells_[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num) checking_bound_(list_data[num].first, dt);
				}
			}, ap);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] > body_domain_bounds_.first[axis_]
			&& particle_position[axis_] < (body_domain_bounds_.first[axis_] + cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorInAxisDirection(expected_particle_index, body_domain_bounds_.first, axis_);
			Vecd translated_position = particles_->pos_n_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection
		::CreatingGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_n_[index_i];
		if (particle_position[axis_] < body_domain_bounds_.second[axis_]
			&& particle_position[axis_] > (body_domain_bounds_.second[axis_] - cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorInAxisDirection(expected_particle_index, body_domain_bounds_.second, axis_);
			Vecd translated_position = particles_->pos_n_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorInAxisDirection(index_i, body_domain_bounds_.first, axis_);
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::UpdatingGhostStates
		::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorInAxisDirection(index_i, body_domain_bounds_.second, axis_);
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
		: ParticleDynamicsReduce<bool, ReduceOR>(body),
		GeneralDataDelegateSimple(body),
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
		UpperFrontInXDirection(SPHBody* body) :
		ParticleDynamicsReduce<Real, ReduceMax>(body),
		GeneralDataDelegateSimple(body),
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
		MaximumSpeed(SPHBody* body) :
		ParticleDynamicsReduce<Real, ReduceMax>(body),
		GeneralDataDelegateSimple(body),
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
	BodyLowerBound::BodyLowerBound(SPHBody* body)
		: ParticleDynamicsReduce<Vecd, ReduceLowerBound>(body),
		GeneralDataDelegateSimple(body),
		pos_n_(particles_->pos_n_)
	{
		constexpr  double max_real_number = (std::numeric_limits<double>::max)();
		initial_reference_ = Vecd(max_real_number);
	}
	//=================================================================================================//
	Vecd BodyLowerBound::ReduceFunction(size_t index_i, Real dt)
	{
		return pos_n_[index_i];
	}
	//=================================================================================================//
	BodyUpperBound::
		BodyUpperBound(SPHBody* body) :
		ParticleDynamicsReduce<Vecd, ReduceUpperBound>(body),
		GeneralDataDelegateSimple(body),
		pos_n_(particles_->pos_n_) 
	{
		constexpr  double min_real_number = (std::numeric_limits<double>::min)();
		initial_reference_ = Vecd(min_real_number);
	}
	//=================================================================================================//
	Vecd BodyUpperBound::ReduceFunction(size_t index_i, Real dt)
	{
		return pos_n_[index_i];
	}
	//=================================================================================================//
	TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody* body, Gravity* gravity)
		: ParticleDynamicsReduce<Real, ReduceSum<Real>>(body), 
		GeneralDataDelegateSimple(body), mass_(particles_->mass_), 
		vel_n_(particles_->vel_n_), pos_n_(particles_->pos_n_), gravity_(gravity)
	{
		quantity_name_ = "TotalMechanicalEnergy";
		initial_reference_ = 0.0;
	}
	//=================================================================================================//
	Real TotalMechanicalEnergy::ReduceFunction(size_t index_i, Real dt)
	{
		return 0.5 * mass_[index_i] * vel_n_[index_i].normSqr()
			+ mass_[index_i] * gravity_->getPotential(pos_n_[index_i]);
	}
	//=================================================================================================//
}
//=================================================================================================//