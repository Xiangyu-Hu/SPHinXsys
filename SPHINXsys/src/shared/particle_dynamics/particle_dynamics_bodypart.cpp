/**
* @file 	particle_dynamics_bodypart.cpp
* @brief 	This is the implementation of the template class
* @author	Chi ZHang and Xiangyu Hu
*/

#include "particle_dynamics_bodypart.h"

namespace SPH
{
	//=================================================================================================//
	void PartDynamicsByParticle::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			particle_functor_(body_part_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void PartDynamicsByParticle::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					particle_functor_(body_part_particles_[i], dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	PartSimpleDynamicsByParticle::
		PartSimpleDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part) 
		: PartDynamicsByParticle(sph_body, body_part)
	{
		particle_functor_ = std::bind(&PartSimpleDynamicsByParticle::Update, this, _1, _2);
	};
	//=================================================================================================//
	PartInteractionDynamicsByParticle::
		PartInteractionDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part) 
		: PartDynamicsByParticle(sph_body, body_part)
	{
		particle_functor_ = std::bind(&PartInteractionDynamicsByParticle::Interaction, this, _1, _2);
	}
	//=================================================================================================//
	PartInteractionDynamicsByParticleWithUpdate::
		PartInteractionDynamicsByParticleWithUpdate(SPHBody &sph_body, BodyPartByParticle &body_part) 
		: PartInteractionDynamicsByParticle(sph_body, body_part)
	{
		functor_update_ = std::bind(&PartInteractionDynamicsByParticleWithUpdate::Update, this, _1, _2);
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticleWithUpdate::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);

		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			particle_functor_(body_part_particles_[i], dt);
		}

		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			functor_update_(body_part_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticleWithUpdate::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);

		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					particle_functor_(body_part_particles_[i], dt);
				}
			},
			ap);

		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					functor_update_(body_part_particles_[i], dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);

		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			functor_initialization_(body_part_particles_[i], dt);
		}

		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			particle_functor_(body_part_particles_[i], dt);
		}

		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			functor_update_(body_part_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);

		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					functor_initialization_(body_part_particles_[i], dt);
				}
			},
			ap);

		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					particle_functor_(body_part_particles_[i], dt);
				}
			},
			ap);

		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					functor_update_(body_part_particles_[i], dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	void PartDynamicsByCell::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t i = 0; i != body_part_cells_.size(); ++i)
		{
			ListDataVector &list_data = body_part_cells_[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num)
				Update(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void PartDynamicsByCell::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		parallel_for(
			blocked_range<size_t>(0, body_part_cells_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ListDataVector &list_data = body_part_cells_[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num)
						Update(list_data[num].first, dt);
				}
			},
			ap);
	}
	//=================================================================================================//
}
