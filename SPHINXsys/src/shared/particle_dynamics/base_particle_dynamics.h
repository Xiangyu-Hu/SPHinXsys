/**
* @file base_particle_dynamics.h
* @brief This is for the base classes of particle dynamics, which describe the
* intection between particles. These interactions are used to define  
* differential operators for surface forces or fluxes in continuum mechanics
* @author  Xiangyu Hu, Luhui Han and Chi Zhang
*/

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_particles.h"
#include "all_materials.h"
#include "neighboring_particle.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
#include "external_force.h"

#include <functional>
using namespace std::placeholders;

namespace SPH 
{
	/** Functor for operation of inner particles. */
	typedef std::function<void(size_t, Real)> InnerFunctor;
	/** Functor for operation of particles bewteen different bodies. */
	typedef std::function<void(size_t, size_t, Real)> ContactFunctor;
	/** Functors for reducing operation of inner particles. */
	template <class ReturnType>
	using ReduceFunctor = std::function<ReturnType(size_t, Real)>;

	/** Iterators for inner functors. sequential computing. */
	void InnerIterator(size_t number_of_particles, InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors. parallel computing. */
	void InnerIterator_parallel(size_t number_of_particles, InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for contact functors. sequential computing. */
	void ContactIterator(StdVec<ListIndexVector*> indexes_interacting_particles, 
		ContactFunctor &contact_functor, Real dt = 0.0);
	/** Iterators for contact functors. parallel computing. */
	void ContactIterator_parallel(StdVec<ListIndexVector*> indexes_interacting_particles, 
		ContactFunctor &contact_functor, Real dt = 0.0);

	/** Iterators for reduce functors. sequential computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt = 0.0);
	/** Iterators for reduce functors. parallel computing. */
	template <class ReturnType, class ReduceFunction, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt = 0.0);

	/** Iterators for inner functors with splitting. sequential computing. */
	void InnerIteratorSplitting(ByCellLists by_cell_lists_particle_indexes, 
		InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. parallel computing. */
	void InnerIteratorSplitting_parallel(ByCellLists by_cell_lists_particle_indexes, 
		InnerFunctor &inner_functor, Real dt = 0.0);

	/** A Functor for Summation */
	template <class ReturnType>
	struct ReduceSum { ReturnType operator () (ReturnType x, ReturnType y) const { return x + y; }; };
	/** A Functor for Maximum */
	struct ReduceMax { Real operator () (Real x, Real y) const { return SMAX(x, y); }; };
	/** A Functor for Minimum */
	struct ReduceMin { Real operator () (Real x, Real y) const { return SMIN(x, y); }; };
	/** A Functor for OR operator */
	struct ReduceOR { bool operator () (bool x, bool y) const { return x || y; }; };
	/** A Functor for lower bound */
	struct ReduceLowerBound {
		Vecd operator () (Vecd x, Vecd y) const {
			Vecd lower_bound;
			for (int i = 0; i < lower_bound.size(); ++i) 
				lower_bound[i] = SMIN(x[i], y[i]);
			return lower_bound;
		}; };
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
	* @class Dynamics
	* @brief The base class for all dynamics
	* This class contains the only two interface functions available
	* for particle dynamics. An specific implementation should be realized.
	*/
	template <class ReturnType>
	class Dynamics : public GlobalStaticVariables
	{
	protected:
		virtual void SetupDynamics(Real dt = 0.0) = 0;
	public:
		/** Constructor */
		explicit Dynamics() : GlobalStaticVariables() {};
		virtual ~Dynamics() {};

		/** The only two functions can be called from outside
		  * One is for sequential excution, the other is for parallel. */
		virtual ReturnType exec(Real dt = 0.0) = 0;
		virtual ReturnType parallel_exec(Real dt = 0.0) = 0;
	};

	/**
	* @class ParticleDynamics
	* @brief Particle dynamics base class
	* Bodies are involved the dynamics of a designated body
	* and the bodies interacting with this body.
	*/
	template <class ReturnType, class BodyType, 
		class ParticlesType = Particles, class MaterialType = Material>
	class ParticleDynamics : public Dynamics<ReturnType>
	{
	protected:
		BodyType * body_;
		ParticlesType *particles_;
		MaterialType *material_;

		virtual void SetupDynamics(Real dt = 0.0) override {};
	public:
		/** Constructor */
		explicit ParticleDynamics(BodyType* body) : Dynamics<ReturnType>(), body_(body), 
			particles_(dynamic_cast<ParticlesType*>(body->base_particles_->PointToThisObject())),
			material_(dynamic_cast<MaterialType*>(body->base_material_->PointToThisObject())) {};
		virtual ~ParticleDynamics() {};
	};

	/**
	* @class ParticleDynamicsWithInnerConfigurations
	* @brief Particle dynamics base class for the case 
	* with the inner configuration
	*/
	template <class BodyType, class ParticlesType = Particles, class MaterialType = Material>	
	class ParticleDynamicsWithInnerConfigurations
		: public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		/** current inner confifuration of the designated body */
		NeighborList *current_inner_configuration_;
		/** reference inner confifuration of the designated body */
		ReferenceNeighborList *reference_inner_configuration_;
	public:
		explicit ParticleDynamicsWithInnerConfigurations(BodyType* body);
		virtual ~ParticleDynamicsWithInnerConfigurations() {};
	};

	/**
	* @class ParticleDynamicsWithContactConfigurations
	* @brief This is the bas class for the case with contact configurations
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
		class ParticleDynamicsWithContactConfigurations
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		StdVec<InteractingBodyType *>  interacting_bodies_;
		StdVec<InteractingParticlesType *>  interacting_particles_;
		StdVec<InteractingMaterialType *>  interacting_material_;

		/** lists of the original body particles contacted to interacting bodies*/
		StdVec<ListIndexVector*> indexes_interacting_particles_;
		/** current interaction confifuration of the interacting bodies */
		StdVec<NeighborList *> current_interacting_configuration_;
		/** reference interaction confifuration of the designated body */
		StdVec<ReferenceNeighborList *> reference_interacting_configuration_;
	public:
		explicit ParticleDynamicsWithContactConfigurations(BodyType *body, 
			StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsWithContactConfigurations() {};
	};

	/**
	* @class ParticleDynamicsByCells
	* @brief Simple particle dynamics base class
	* The particles are iterated according a list of cells.
	* This is mainly used for imposing Eulerian boundary conditions.
	*/
	template <class BodyType, class ParticlesType = Particles, class MaterialType = Material>
	class ParticleDynamicsByCells : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		MeshCellLinkedList *mesh_cell_linked_list_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Real cell_spacing_;
		Vecd mesh_lower_bound_, mesh_upper_bound_;

	public:
		/** Constructor */
		explicit ParticleDynamicsByCells(BodyType* body);
		virtual ~ParticleDynamicsByCells() {};
	};

	/**
	* @class ParticleDynamicsSimple
	* @brief Simple particle dynamics base class
	*/
	template <class BodyType, class ParticlesType = Particles, class MaterialType = Material>
	class ParticleDynamicsSimple : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;
		public:
		explicit ParticleDynamicsSimple(BodyType* body)
			: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body), 
			functor_update_(std::bind(&ParticleDynamicsSimple::Update, this, _1, _2)) {};
		virtual ~ParticleDynamicsSimple() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/**
	* @class ParticleDynamicsReduce
	* @brief Base abstract class for reduce
	*/
	template <class ReturnType, typename ReduceOperation,
		class BodyType, class ParticlesType = Particles, class MaterialType = Material>
	class ParticleDynamicsReduce : public ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>
	{
	protected:
		ReduceOperation reduce_operation_;

		//inital or refence value
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
		ReduceFunctor<ReturnType> functor_reduce_function_;
	public:
		explicit ParticleDynamicsReduce(BodyType* body) : ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>(body), 
			functor_reduce_function_(std::bind(&ParticleDynamicsReduce::ReduceFunction, this, _1, _2)) {};
		virtual ~ParticleDynamicsReduce() {};
	
		virtual ReturnType exec(Real dt = 0.0) override;
		virtual ReturnType parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsInner
	* @brief This is the class for inner interactions
	* in which one the particles from the same body
	* interact with each other
	*/
	template <class BodyType, class ParticlesType = Particles, class MaterialType = Material>
	class ParticleDynamicsInner : public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
	public:
		explicit ParticleDynamicsInner(BodyType* body) : 
			ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body),
			functor_inner_interaction_(std::bind(&ParticleDynamicsInner::InnerInteraction, this, _1, _2)) {};
		virtual ~ParticleDynamicsInner() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class ParticleDynamicsInnerSplitting
	 * @brief This is for the splitting algorihm
	 */
	template <class BodyType, class ParticlesType, class MaterialType = Material>
	class ParticleDynamicsInnerSplitting 
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		ByCellLists by_cell_lists_particle_indexes_;

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
	public:
		explicit ParticleDynamicsInnerSplitting(BodyType* body)
			: ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body),
			by_cell_lists_particle_indexes_(body->by_cell_lists_particle_indexes_),
			functor_inner_interaction_(std::bind(&ParticleDynamicsInnerSplitting::InnerInteraction, this, _1, _2)) {};
		virtual ~ParticleDynamicsInnerSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class ParticleDynamicsContact
	 * @brief This is the class for contact interactions
	 */
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsContact 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType, 
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		ContactFunctor functor_contact_interaction_;
	public:
		explicit ParticleDynamicsContact(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
			: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType, 
				InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies), 
			functor_contact_interaction_(std::bind(&ParticleDynamicsContact::ContactInteraction, this, _1, _2, _3)) {};
		virtual ~ParticleDynamicsContact() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}
