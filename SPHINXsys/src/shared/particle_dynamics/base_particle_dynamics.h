/**
* @file base_particle_dynamics.h
* @brief This is the base classes of particle dynamics, which describe the
* intection between particles. These interactions are used to define  
*differential operators for surface forces or fluxes in continuum mechanics
* @author  Xiangyu Hu, Luhui Han and Chi Zhang
*/

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_particles.h"
#include "neighboring_particle.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
#include "external_force.h"

namespace SPH 
{
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
	*/
	template <class ReturnType>
	class Dynamics : public GlobalStaticVariables
	{
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
	* This class contains all the interface functions available
	* for particle dynamics. An specific implementation should be realized.
	* Bodies are involved the dynamics of a designated body
	* and the bodies interacting with this body.
	*/
	template <class ReturnType, class BodyType, class ParticlesType>
	class ParticleDynamics : public Dynamics<ReturnType>
	{
	protected:
		BodyType * body_;
		ParticlesType *particles_;
	public:
		/** Constructor */
		explicit ParticleDynamics(BodyType* body) : Dynamics(), body_(body),
			particles_(dynamic_cast<ParticlesType*>(body->base_particles_.PointToThisObject())) {};
		virtual ~ParticleDynamics() {};
	};

	//------------------------------------------------------------------------//
	/**
	* @class ParticleDynamicsWithInnerConfigurations
	* @brief Particle dynamics base class for the case 
	* with the inner configuration
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsWithInnerConfigurations : public ParticleDynamics<void, BodyType, ParticlesType>
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
	* @brief This is the bas class for the case contact configurations
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsWithContactConfigurations 
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType>
	{
	protected:
		StdVec<InteractingBodytype*>  interacting_bodies_;
		StdVec<InteractingParticlesType*>  interacting_particles_;
		StdVec<ListIndexVector*> indexes_interacting_particles_;
		/** current interaction confifuration of the designated body */
		StdVec<NeighborList *> current_interacting_configuration_;
		/** reference interaction confifuration of the designated body */
		StdVec<ReferenceNeighborList *> reference_interacting_configuration_;
	public:
		explicit ParticleDynamicsWithContactConfigurations(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies);
		virtual ~ParticleDynamicsWithContactConfigurations() {};
	};

	/**
	* @class ParticleDynamicsByCells
	* @brief Simple particle dynamics base class
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsByCells : public ParticleDynamics<void, BodyType, ParticlesType>
	{
	protected:
		/** information of mesh cell linked list */
		MeshCellLinkedList & mesh_cell_linked_list_;
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
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsSimple : public ParticleDynamics<void, BodyType, ParticlesType>
	{
	protected:
		size_t number_of_particles_;
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void ParticleUpdate(size_t index_particle_i, Real dt = 0.0) = 0;
	public:
		/** Constructor */
		explicit ParticleDynamicsSimple(BodyType* body) 
			: ParticleDynamics<void, BodyType, ParticlesType>(body) {
			number_of_particles_ = body->number_of_particles_;
		};
		virtual ~ParticleDynamicsSimple() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsReduce
	* @brief Base abstract class for reduce
	*/
	template <class ReturnType, class BodyType, class ParticlesType>
	class ParticleDynamicsReduce : public ParticleDynamics<ReturnType, BodyType, ParticlesType>
	{
	protected:
		size_t number_of_particles_;
		//inital or refence value
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
		virtual ReturnType ReduceOperation(ReturnType x, ReturnType y) = 0;
	public:
		/** Constructor */
		explicit ParticleDynamicsReduce(BodyType* body)
			: ParticleDynamics<ReturnType, BodyType, ParticlesType>(body) {
			number_of_particles_ = body->number_of_particles_;
		};
		virtual ~ParticleDynamicsReduce() {};

		virtual ReturnType exec(Real dt = 0.0) override;
		virtual ReturnType parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsMinimum
	* @brief For computing minimum value
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsMinimum : public ParticleDynamicsReduce<Real, BodyType, ParticlesType>
	{
	protected:
		virtual Real ReduceOperation(Real x, Real y) override { return SMIN(x, y); };
	public:
		ParticleDynamicsMinimum(BodyType* body)
			: ParticleDynamicsReduce<Real, BodyType, ParticlesType>(body) {};
		virtual ~ParticleDynamicsMinimum() {};
	};

	/**
	* @class ParticleDynamicsMaximum
	* @brief For computing maximum value
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsMaximum : public ParticleDynamicsReduce<Real, BodyType, ParticlesType>
	{
	protected:
		virtual Real ReduceOperation(Real x, Real y) override { return SMAX(x, y); };
	public:
		ParticleDynamicsMaximum(BodyType* body)
			: ParticleDynamicsReduce<Real, BodyType, ParticlesType>(body) {};
		virtual ~ParticleDynamicsMaximum() {};
	};

	/**
	* @class ParticleDynamicsSum
	* @brief For computing sum
	*/
	template <class ReturnType, class BodyType, class ParticlesType>
	class ParticleDynamicsSum : public ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>
	{
	protected:

		virtual ReturnType ReduceOperation(ReturnType x, ReturnType y) override { return x + y; };
	public:
		ParticleDynamicsSum(BodyType* body)
			: ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>(body) {};
		virtual ~ParticleDynamicsSum() {};
	};

	/**
	* @class ParticleDynamicsOR
	* @brief This for asserting the solution of the simulation
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsOR : public ParticleDynamicsReduce<bool, BodyType, ParticlesType>
	{
	protected:

		virtual bool ReduceOperation(bool x, bool y) override { return x || y; };
	public:
		ParticleDynamicsOR(BodyType* body)
			: ParticleDynamicsReduce<bool, BodyType, ParticlesType>(body) {};
		virtual ~ParticleDynamicsOR() {};
	};

	/**
	* @class ParticleDynamicsInner
	* @brief This is the class for inner interactions
	* in which one the particles from the same body
	* interact with each other
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsInner
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType>
	{
	protected:
		size_t number_of_particles_;
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;

	public:
		ParticleDynamicsInner(BodyType* body) : ParticleDynamicsWithInnerConfigurations(body) {
			number_of_particles_ = body->number_of_particles_;
		};
		virtual ~ParticleDynamicsInner() {};

		virtual void exec(Real dt = 0.0);
		virtual void parallel_exec(Real dt = 0.0);
	};

	/**
	 * @class ParticleDynamicsInnerSplitting
	 * @brief This is for the splitting algorihm
	 */
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsInnerSplitting
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType>
	{
	protected:
		int numer_of_lists_;
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;

	public:
		ParticleDynamicsInnerSplitting(BodyType* body) : ParticleDynamicsWithInnerConfigurations(body) {
			numer_of_lists_ = powern(3, Vecd(0).size());
		};
		virtual ~ParticleDynamicsInnerSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsContact
	 * @brief This is the class for contact interactions
	 */
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsContact
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;

	public:
		explicit ParticleDynamicsContact(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsWithContactConfigurations(body, interacting_bodies) {};
		virtual ~ParticleDynamicsContact() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}
