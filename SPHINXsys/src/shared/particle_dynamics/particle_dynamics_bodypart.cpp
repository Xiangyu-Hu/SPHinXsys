/**
 * @file 	particle_dynamics_bodypart.cpp
 * @brief 	This is the implementation of the template class
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "particle_dynamics_bodypart.h"

namespace SPH
{
	//=================================================================================================//
	void PartIteratorByParticle(const IndexVector &body_part_particles, const ParticleFunctor &particle_functor, Real dt)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
		{
			particle_functor(body_part_particles[i], dt);
		}
	}
	//=================================================================================================//
	void PartIteratorByParticle_parallel(const IndexVector &body_part_particles, const ParticleFunctor &particle_functor, Real dt)
	{
		parallel_for(
			blocked_range<size_t>(0, body_part_particles.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					particle_functor(body_part_particles[i], dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	void PartIteratorByCell(const CellLists &body_part_cells, const ParticleFunctor &particle_functor, Real dt)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ListDataVector &list_data = body_part_cells[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num)
				particle_functor(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void PartIteratorByCell_parallel(const CellLists &body_part_cells, const ParticleFunctor &particle_functor, Real dt)
	{
		parallel_for(
			blocked_range<size_t>(0, body_part_cells.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ListDataVector &list_data = body_part_cells[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num)
						particle_functor(list_data[num].first, dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	PartDynamicsByParticle::PartDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part)
		: ParticleDynamics<void>(sph_body),
		  body_part_particles_(body_part.body_part_particles_) {}
	//=================================================================================================//
	PartSimpleDynamicsByParticle::
		PartSimpleDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part)
		: PartDynamicsByParticle(sph_body, body_part),
		  functor_update_(std::bind(&PartSimpleDynamicsByParticle::Update, this, _1, _2)) {}
	//=================================================================================================//
	void PartSimpleDynamicsByParticle::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	void PartSimpleDynamicsByParticle::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	PartInteractionDynamicsByParticle::
		PartInteractionDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part)
		: PartDynamicsByParticle(sph_body, body_part),
		  functor_interaction_(std::bind(&PartInteractionDynamicsByParticle::Interaction, this, _1, _2)) {}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle(body_part_particles_, functor_interaction_, dt);
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_interaction_, dt);
	}
	//=================================================================================================//
	PartInteractionDynamicsByParticleWithUpdate::
		PartInteractionDynamicsByParticleWithUpdate(SPHBody &sph_body, BodyPartByParticle &body_part)
		: PartInteractionDynamicsByParticle(sph_body, body_part),
		  functor_update_(std::bind(&PartInteractionDynamicsByParticleWithUpdate::Update, this, _1, _2)) {}
	//=================================================================================================//
	void PartInteractionDynamicsByParticleWithUpdate::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle(body_part_particles_, functor_interaction_, dt);
		PartIteratorByParticle(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticleWithUpdate::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_interaction_, dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	PartInteractionDynamicsByParticle1Level::
		PartInteractionDynamicsByParticle1Level(SPHBody &sph_body, BodyPartByParticle &body_part)
		: PartInteractionDynamicsByParticleWithUpdate(sph_body, body_part),
		  functor_initialization_(
			  std::bind(&PartInteractionDynamicsByParticle1Level::Initialization, this, _1, _2)) {}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle(body_part_particles_, functor_initialization_, dt);
		PartIteratorByParticle(body_part_particles_, functor_interaction_, dt);
		PartIteratorByParticle(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	void PartInteractionDynamicsByParticle1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_initialization_, dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_interaction_, dt);
		PartIteratorByParticle_parallel(body_part_particles_, functor_update_, dt);
	}
	//=================================================================================================//
	PartDynamicsByCell::PartDynamicsByCell(SPHBody &sph_body, BodyPartByCell &body_part)
		: ParticleDynamics<void>(sph_body),
		  body_part_cells_(body_part.body_part_cells_),
		  functor_update_(std::bind(&PartDynamicsByCell::Update, this, _1, _2)){};
	//=================================================================================================//
	void PartDynamicsByCell::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByCell(body_part_cells_, functor_update_, dt);
	}
	//=================================================================================================//
	void PartDynamicsByCell::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		PartIteratorByCell_parallel(body_part_cells_, functor_update_, dt);
	}
	//=================================================================================================//
}
