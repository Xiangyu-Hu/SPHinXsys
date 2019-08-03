/**
* @file 	particle_dynamics_algorithms.hpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "particle_dynamics_algorithms.h"

namespace SPH {
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsInner1Level(BodyType* body)
		: ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body),
		functor_initialization_(std::bind(&ParticleDynamicsInner1Level::Initialization, this, _1, _2)),
		functor_inner_interaction_(std::bind(&ParticleDynamicsInner1Level::InnerInteraction, this, _1, _2)),
		functor_update_(std::bind(&ParticleDynamicsInner1Level::Update, this, _1, _2)) {
		number_of_particles_ = body->number_of_particles_;
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles_, functor_initialization_, dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_, dt);
		InnerIterator(number_of_particles_, functor_update_, dt);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles_, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_, dt);
		InnerIterator_parallel(number_of_particles_, functor_update_, dt);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplex(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_inner_interaction_(std::bind(&ParticleDynamicsComplex::InnerInteraction, this, _1, _2)),
		functor_contact_interaction_(std::bind(&ParticleDynamicsComplex::ContactInteraction, this, _1, _2, _3)) 
	{
		number_of_particles_ = body->number_of_particles_;
	}
	//===============================================================//
template <class BodyType, class ParticlesType, class MaterialType,
	class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplexWithUpdate(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_inner_interaction_(std::bind(&ParticleDynamicsComplexWithUpdate::InnerInteraction, this, _1, _2)),
		functor_contact_interaction_(std::bind(&ParticleDynamicsComplexWithUpdate::ContactInteraction, this, _1, _2, _3)),
		functor_update_(std::bind(&ParticleDynamicsComplexWithUpdate::Update, this, _1, _2)) {
		number_of_particles_ = body->number_of_particles_;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator_parallel(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplex1Level(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies), 
		functor_initialization_(std::bind(&ParticleDynamicsComplex1Level::Initialization, this, _1, _2)),
		functor_inner_interaction_(std::bind(&ParticleDynamicsComplex1Level::InnerInteraction, this, _1, _2)),
		functor_contact_interaction_(std::bind(&ParticleDynamicsComplex1Level::ContactInteraction, this, _1, _2, _3)),
		functor_update_(std::bind(&ParticleDynamicsComplex1Level::Update, this, _1, _2)) {
		number_of_particles_ = body->number_of_particles_;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles_, functor_initialization_, dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles_, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator_parallel(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplex2Levels<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplex2Levels(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_initialization_(std::bind(&ParticleDynamicsComplex2Levels::Initialization, this, _1, _2)),
		functor_inner_interaction_(std::bind(&ParticleDynamicsComplex2Levels::InnerInteraction, this, _1, _2)),
		functor_contact_interaction_(std::bind(&ParticleDynamicsComplex2Levels::ContactInteraction, this, _1, _2, _3)),
		functor_intermediate_(std::bind(&ParticleDynamicsComplex2Levels::Intermediate, this, _1, _2)),
		functor_inner_interaction_2nd_(std::bind(&ParticleDynamicsComplex2Levels::InnerInteraction2nd, this, _1, _2)),
		functor_contact_interaction_2nd_(std::bind(&ParticleDynamicsComplex2Levels::ContactInteraction2nd, this, _1, _2, _3)),
		functor_update_(std::bind(&ParticleDynamicsComplex2Levels::Update, this, _1, _2)) {
		number_of_particles_ = body->number_of_particles_;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles_, functor_initialization_, dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator(number_of_particles_, functor_intermediate_, dt);
		InnerIterator(number_of_particles_, functor_inner_interaction_2nd_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_2nd_, dt);
		InnerIterator(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles_, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
		InnerIterator_parallel(number_of_particles_, functor_intermediate_, dt);
		InnerIterator_parallel(number_of_particles_, functor_inner_interaction_2nd_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_2nd_, dt);
		InnerIterator_parallel(number_of_particles_, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplexSplitting(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies), 
		by_cell_lists_particle_indexes_(body->by_cell_lists_particle_indexes_),
		functor_inner_interaction_(std::bind(&ParticleDynamicsComplexSplitting::InnerInteraction, this, _1, _2)),
		functor_contact_interaction_(std::bind(&ParticleDynamicsComplexSplitting::ContactInteraction, this, _1, _2, _3)) 
	{
	
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting(by_cell_lists_particle_indexes_, functor_inner_interaction_, dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting_parallel(by_cell_lists_particle_indexes_, functor_inner_interaction_, dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
}
