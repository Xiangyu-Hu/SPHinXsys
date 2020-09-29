/**
* @file 	particle_dynamics_algorithms.cpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#include "particle_dynamics_algorithms.h"

//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	void ParticleDynamicsSimple::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsSimple::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInner::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInner::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInnerWithUpdate::exec(Real dt)
	{
		ParticleDynamicsInner::exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInnerWithUpdate::parallel_exec(Real dt)
	{
		ParticleDynamicsInner::parallel_exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInner1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_initialization_, dt);
		InnerIterator(number_of_particles, functor_inner_interaction_, dt);
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInner1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles, functor_inner_interaction_, dt);
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsContact::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsContact::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplex::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplex::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplexWithUpdate::exec(Real dt)
	{
		ParticleDynamicsComplex::exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplexWithUpdate::parallel_exec(Real dt)
	{
		ParticleDynamicsComplex::parallel_exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplex1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_initialization_, dt);
		InnerIterator(number_of_particles, functor_complex_interaction_, dt);
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsComplex1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles, functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//===============================================================//
	void ParticleDynamicsComplexSplit::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_initialization_, dt);
		InnerIteratorSplitting(split_cell_lists_, functor_complex_interaction_, dt);
		InnerIterator(number_of_particles, this->functor_update_, dt);
	}
	//===============================================================//
	void ParticleDynamicsComplexSplit::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIteratorSplitting_parallel(split_cell_lists_, functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	ParticleDynamicsCellListSplitting
		::ParticleDynamicsCellListSplitting(SPHBody* body)
		: ParticleDynamics<void>(body),
		functor_cell_list_(std::bind(&ParticleDynamicsCellListSplitting::CellListInteraction, this, _1, _2))
	{
		cutoff_radius_ = mesh_cell_linked_list_->CellSpacing();
		kernel_ = body->kernel_;
	};
	//=================================================================================================//
	void ParticleDynamicsCellListSplitting::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		CellListIteratorSplitting(split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsCellListSplitting::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		CellListIteratorSplitting_parallel(split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInnerSplitting::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		InnerIteratorSplittingSweeping(split_cell_lists_, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsInnerSplitting::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		InnerIteratorSplittingSweeping_parallel(split_cell_lists_, functor_inner_interaction_, dt);
	}
	//=============================================================================================//
	void ParticleDynamicsComplexSplitting::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		InnerIteratorSplittingSweeping(split_cell_lists_, functor_particle_interaction_, dt);
	}
	//=============================================================================================//
	void ParticleDynamicsComplexSplitting::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		InnerIteratorSplittingSweeping_parallel(split_cell_lists_, functor_particle_interaction_, dt);
	}
	//=============================================================================================//
}
//=================================================================================================//