/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file base_particle_dynamics.h
* @brief This is for the base classes of particle dynamics, which describe the
* interaction between particles. These interactions are used to define  
* differential operators for surface forces or fluxes in continuum mechanics
* @author  Xiangyu Hu, Luhui Han and Chi Zhang
*/

#ifndef BASE_PARTICLE_DYNAMICS_H
#define BASE_PARTICLE_DYNAMICS_H


#include "base_data_package.h"
#include "sph_data_containers.h"
#include "all_particles.h"
#include "all_materials.h"
#include "neighbor_relation.h"
#include "all_bodies.h"
#include "cell_linked_list.h"
#include "external_force.h"
#include "body_relation.h"
#include <functional>

using namespace std::placeholders;

namespace SPH 
{
	/** Functor for operation on particles. */
	typedef std::function<void(size_t, Real)> ParticleFunctor;
	/** Functors for reducing operation on particles. */
	template <class ReturnType>
	using ReduceFunctor = std::function<ReturnType(size_t, Real)>;

	/** Iterators for particle functors. sequential computing. */
	void ParticleIterator(size_t total_real_particles, ParticleFunctor &particle_functor, Real dt = 0.0);
	/** Iterators for particle functors. parallel computing. */
	void ParticleIterator_parallel(size_t total_real_particles, ParticleFunctor &particle_functor, Real dt = 0.0);

	/** Iterators for reduce functors. sequential computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t total_real_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt = 0.0);
	/** Iterators for reduce functors. parallel computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t total_real_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt = 0.0);

	/** Iterators for particle functors with splitting. sequential computing. */
	void ParticleIteratorSplittingSweep(SplitCellLists& split_cell_lists,
		ParticleFunctor& particle_functor, Real dt = 0.0);
	/** Iterators for particle functors with splitting. parallel computing. */
	void ParticleIteratorSplittingSweep_parallel(SplitCellLists& split_cell_lists,
		ParticleFunctor& particle_functor, Real dt = 0.0);


	/** A Functor for Summation */
	template <class ReturnType>
	struct ReduceSum { ReturnType operator () (const ReturnType& x, const ReturnType& y) const { return x + y; }; };
	/** A Functor for Maximum */
	struct ReduceMax { Real operator () (Real x, Real y) const { return SMAX(x, y); }; };
	/** A Functor for Minimum */
	struct ReduceMin { Real operator () (Real x, Real y) const { return SMIN(x, y); }; };
	/** A Functor for OR operator */
	struct ReduceOR { bool operator () (bool x, bool y) const { return x || y; }; };
	/** A Functor for lower bound */
	struct ReduceLowerBound {
		Vecd operator () (const Vecd& x, const Vecd& y) const {
			Vecd lower_bound;
			for (int i = 0; i < lower_bound.size(); ++i) lower_bound[i] = SMIN(x[i], y[i]);
			return lower_bound;
		}; 
	};
	/** A Functor for upper bound */
	struct ReduceUpperBound {
		Vecd operator () (const Vecd& x, const Vecd& y) const {
			Vecd upper_bound;
			for (int i = 0; i < upper_bound.size(); ++i) upper_bound[i] = SMAX(x[i], y[i]);
			return upper_bound;
		};
	};

	/**
	 * @class GlobalStaticVariables
	 * @brief A place to put all global variables
	 */
	class GlobalStaticVariables
	{
	public:
		explicit GlobalStaticVariables() {};
		virtual ~GlobalStaticVariables() {};

		/** the physical time is global value for all dynamics */
		static Real physical_time_;
	};

	/**
	* @class ParticleDynamics
	* @brief The base class for all particle dynamics
	* This class contains the only two interface functions available
	* for particle dynamics. An specific implementation should be realized.
	*/
	template <class ReturnType = void>
	class ParticleDynamics : public GlobalStaticVariables
	{
	public:
		explicit ParticleDynamics(SPHBody* sph_body)
			: GlobalStaticVariables(), sph_body_(sph_body),
			particle_adaptation_(sph_body->particle_adaptation_),
			base_particles_(sph_body->base_particles_) {};
		virtual ~ParticleDynamics() {};

		SPHBody* getSPHBody() { return sph_body_; };
		/** The only two functions can be called from outside
		  * One is for sequential execution, the other is for parallel. */
		virtual ReturnType exec(Real dt = 0.0) = 0;
		virtual ReturnType parallel_exec(Real dt = 0.0) = 0;
	protected:
		SPHBody* sph_body_;
		ParticleAdaptation* particle_adaptation_;
		BaseParticles* base_particles_;

		void setBodyUpdated() { sph_body_->setNewlyUpdated(); };
		/** the function for set global parameters for the particle dynamics */
		virtual void setupDynamics(Real dt = 0.0) {};
	};
	
	/**
	* @class DataDelegateBase
	* @brief empty base class mixin template.
	*/
	class DataDelegateEmptyBase
	{
	public:
		explicit DataDelegateEmptyBase(SPHBody* sph_body) {};
		virtual ~DataDelegateEmptyBase() {};
	};

	/**
	* @class DataDelegateSimple
	* @brief prepare data for simple particle dynamics.
	*/
	template <class BodyType = SPHBody, 
			  class ParticlesType = BaseParticles, 
			  class MaterialType = BaseMaterial>
	class DataDelegateSimple
	{
	public:
		explicit DataDelegateSimple(SPHBody* body) :
			body_(dynamic_cast<BodyType*>(body)),
			particles_(dynamic_cast<ParticlesType*>(body->base_particles_)),
			material_(dynamic_cast<MaterialType*>(body->base_particles_->base_material_)),
			sorted_id_(body_->base_particles_->sorted_id_),
			unsorted_id_(body_->base_particles_->unsorted_id_) {};
		virtual ~DataDelegateSimple() {};
	protected:
		BodyType* body_;
		ParticlesType* particles_;
		MaterialType* material_;
		StdLargeVec<size_t>& sorted_id_;
		StdLargeVec<size_t>& unsorted_id_;
	};

	/**
	* @class DataDelegateInner
	* @brief prepare data for inner particle dynamics
	*/
	template <class BodyType = SPHBody,
			  class ParticlesType = BaseParticles,
			  class MaterialType = BaseMaterial,
			  class BaseDataDelegateType = DataDelegateSimple<BodyType, ParticlesType, MaterialType>>
	class DataDelegateInner : public BaseDataDelegateType
	{
	public:
		explicit DataDelegateInner(BaseBodyRelationInner* body_inner_relation) : 
			BaseDataDelegateType(body_inner_relation->sph_body_),
			inner_configuration_(body_inner_relation->inner_configuration_) {};
		virtual ~DataDelegateInner() {};
	protected:
		/** inner configuration of the designated body */
		ParticleConfiguration& inner_configuration_;
	};

	/**
	* @class DataDelegateContact
	* @brief prepare data for contact particle dynamics
	*/
	template <class BodyType = SPHBody,
			  class ParticlesType = BaseParticles,
			  class MaterialType = BaseMaterial,
			  class ContactBodyType = SPHBody,
			  class ContactParticlesType = BaseParticles,
			  class ContactMaterialType = BaseMaterial,
			  class BaseDataDelegateType = DataDelegateSimple<BodyType, ParticlesType, MaterialType>>
	class DataDelegateContact : public BaseDataDelegateType
	{
	public:
		explicit DataDelegateContact(BaseBodyRelationContact* body_contact_relation);
		virtual ~DataDelegateContact() {};
	protected:
		StdVec<ContactBodyType*>  contact_bodies_;
		StdVec<ContactParticlesType*>  contact_particles_;
		StdVec<ContactMaterialType*>  contact_material_;
		/** Configurations for particle interaction between bodies. */
		StdVec<ParticleConfiguration*> contact_configuration_;
	};

	/**
	* @class DataDelegateComplex
	* @brief prepare data for complex particle dynamics
	*/
	template <class BodyType = SPHBody,
		class ParticlesType = BaseParticles,
		class MaterialType = BaseMaterial,
		class ContactBodyType = SPHBody,
		class ContactParticlesType = BaseParticles,
		class ContactMaterialType = BaseMaterial>
	class DataDelegateComplex : 
		public DataDelegateInner<BodyType, ParticlesType, MaterialType>,
		public DataDelegateContact<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType, DataDelegateEmptyBase>
	{
	public:
		explicit DataDelegateComplex(ComplexBodyRelation* body_complex_relation) :
			DataDelegateInner<BodyType, ParticlesType, MaterialType>(body_complex_relation->inner_relation_),
			DataDelegateContact<BodyType, ParticlesType, MaterialType, ContactBodyType, ContactParticlesType, 
			ContactMaterialType, DataDelegateEmptyBase>(body_complex_relation->contact_relation_) {};
		virtual ~DataDelegateComplex() {};
	};

	/**
	* @class ParticleDynamicsComplex
	* @brief particle dynamics by considering  contribution from extra contact bodies
	*/
	template<class ParticleDynamicsInnerType, class ContactDataType>
	class ParticleDynamicsComplex : public ParticleDynamicsInnerType, public ContactDataType
	{
	public:
		ParticleDynamicsComplex(BaseBodyRelationInner* inner_relation, 
			BaseBodyRelationContact* contact_relation) :
			ParticleDynamicsInnerType(inner_relation), ContactDataType(contact_relation) {};

		ParticleDynamicsComplex(ComplexBodyRelation* complex_relation, 
			BaseBodyRelationContact* extra_contact_relation);

		virtual ~ParticleDynamicsComplex() {};
	protected:
		virtual void prepareContactData() = 0;
	};
}
#endif //BASE_PARTICLE_DYNAMICS_H