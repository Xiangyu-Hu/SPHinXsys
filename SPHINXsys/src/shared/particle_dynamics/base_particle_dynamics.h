/**
* @file base_particle_dynamics.h
* @brief This is for the base classes of particle dynamics, which describe the
* interaction between particles. These interactions are used to define  
* differential operators for surface forces or fluxes in continuum mechanics
* @author  Xiangyu Hu, Luhui Han and Chi Zhang
*/
#pragma once
#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_particles.h"
#include "all_materials.h"
#include "neighbor_relation.h"
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
	void ContactIterator(InteractingParticles& indexes_interacting_particles,
		ContactFunctor &contact_functor, Real dt = 0.0);
	/** Iterators for contact functors. parallel computing. */
	void ContactIterator_parallel(InteractingParticles& indexes_interacting_particles,
		ContactFunctor &contact_functor, Real dt = 0.0);

	/** Iterators for reduce functors. sequential computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt = 0.0);
	/** Iterators for reduce functors. parallel computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt = 0.0);

	/** Functor for cofiguration operation. */
	typedef std::function<void(CellList*, Real)> CellListFunctor;
	/** Iterators for inner functors with splitting for configuration dynamics. sequential computing. */
	void CellListIteratorSplitting(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting for configuration dynamics. parallel computing. */
	void CellListIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt = 0.0);

	/** Iterators for inner functors with splitting. sequential computing. */
	void InnerIteratorSplitting(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. parallel computing. */
	void InnerIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. sequential computing. */
	void InnerIteratorSplittingSweeping(SplitCellLists& split_cell_lists,
		InnerFunctor& inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. parallel computing. */
	void InnerIteratorSplittingSweeping_parallel(SplitCellLists& split_cell_lists,
		InnerFunctor& inner_functor, Real dt = 0.0);


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
			for (int i = 0; i < lower_bound.size(); ++i) lower_bound[i] = SMIN(x[i], y[i]);
			return lower_bound;
		}; 
	};
	/** A Functor for upper bound */
	struct ReduceUpperBound {
		Vecd operator () (Vecd x, Vecd y) const {
			Vecd upper_bound;
			for (int i = 0; i < upper_bound.size(); ++i) upper_bound[i] = SMIN(x[i], y[i]);
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
		class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
	class ParticleDynamics : public Dynamics<ReturnType>
	{
	protected:
		/** the body involving the particle dynamics */
		BodyType * body_;
		/** the particles involving the particle dynamics */
		ParticlesType *particles_;
		/** the material involving the particle dynamics */
		MaterialType *material_;
		/** Split cell lists*/
		SplitCellLists& split_cell_lists_;


		/** the function for set global parameters for the particle dynamics */
		virtual void SetupDynamics(Real dt = 0.0) override {};
	public:
		/** Constructor */
		explicit ParticleDynamics(BodyType* body) : Dynamics<ReturnType>(), body_(body), 
			particles_(dynamic_cast<ParticlesType*>(body->base_particles_->PointToThisObject())),
			material_(dynamic_cast<MaterialType*>(body->base_particles_->base_material_->PointToThisObject())),
			split_cell_lists_(body->split_cell_lists_){};
		virtual ~ParticleDynamics() {};
	};

	/**
	* @class ParticleDynamicsWithInnerConfigurations
	* @brief Particle dynamics base class for the case 
	* with the inner configuration
	*/
	template <class BodyType, class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
	class ParticleDynamicsWithInnerConfigurations
		: public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		/** inner confifuration of the designated body */
		ParticleConfiguration* inner_configuration_;
		/** Get neighbor list for particle interaction. */
		NeighborList& getNeighborList(ParticleConfiguration* particle_configuration, 
			size_t index_particle_i) {
			Neighborhood& neighborhood = (*particle_configuration)[index_particle_i];
			return std::get<0>(neighborhood);
		}
	public:
		explicit ParticleDynamicsWithInnerConfigurations(BodyType* body);
		virtual ~ParticleDynamicsWithInnerConfigurations() {};
	};

	/**
	* @class ParticleDynamicsWithContactConfigurations
	* @brief This is the bas class for the case with contact configurations
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsWithContactConfigurations
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		StdVec<InteractingBodyType*>  interacting_bodies_;
		StdVec<InteractingParticlesType*>  interacting_particles_;
		StdVec<InteractingMaterialType*>  interacting_material_;

		/** lists of the original body particles contacted to interacting bodies*/
		InteractingParticles indexes_interacting_particles_;
		/** current interaction confifuration of the interacting bodies */
		InteractingParticleConfiguration current_interacting_configuration_;
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
	template <class BodyType, class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
	class ParticleDynamicsByCells : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		BaseMeshCellLinkedList *mesh_cell_linked_list_;
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
	  * @class ParticleDynamicsCellListSplitting
	  * @brief This is for using splitting algorihm 
	  * which does not use particle conmfiguration data for 
	  * particle interaction
	  */
	template <class BodyType, class ParticlesType, class MaterialType = BaseMaterial>
	class ParticleDynamicsCellListSplitting
		: public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		BaseMeshCellLinkedList* mesh_cell_linked_list_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Kernel* kernel_;
		Real cutoff_radius_;
		Real cell_spacing_;
		Vecd mesh_lower_bound_, mesh_upper_bound_;

		virtual void CellListInteraction(CellList* cell_list, Real dt = 0.0) = 0;
		CellListFunctor functor_cell_list_;
	public:
		explicit ParticleDynamicsCellListSplitting(BodyType* body);
		virtual ~ParticleDynamicsCellListSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/**
	 * @class ParticleDynamicsInnerSplitting
	 * @brief This is for the splitting algorihm
	 */
	template <class BodyType, class ParticlesType, class MaterialType = BaseMaterial>
	class ParticleDynamicsInnerSplitting 
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
	public:
		explicit ParticleDynamicsInnerSplitting(BodyType* body)
			: ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body),
			functor_inner_interaction_(std::bind(&ParticleDynamicsInnerSplitting::InnerInteraction, this, _1, _2)) {};
		virtual ~ParticleDynamicsInnerSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class ParticleDynamicsComplexSplitting
	 * @brief This is for the splitting algorihm
	 * which taking account wall boundary conditions
	 */
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
	class ParticleDynamicsComplexSplitting
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		/** the particle interaction also taking account wall boundary conditions,
		  * but the function form is the same as the inner interaction. */
		virtual void ParticleInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_particle_interaction_;
	public:
		explicit ParticleDynamicsComplexSplitting(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
			: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
			functor_particle_interaction_(std::bind(&ParticleDynamicsComplexSplitting::ParticleInteraction, this, _1, _2)) {};
		virtual ~ParticleDynamicsComplexSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}
