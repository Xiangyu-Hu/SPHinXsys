/**
 * @file sph_system.h
 * @brief The SPH_System managing objects in the system level.
 * @details Note that the system only provided functions and data, and
 * 			execute these functions are determined in main file of the examples.
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */
#pragma once

#include "base_data_package.h"
#include "output.h"

namespace SPH 
{
	/**
	 * @brief Friend classes.
	 */
	class Kernel;
	class SPHBody;
	class ParticleGeneratorLattice;
	class ParticleGeneratorDirect;

	/**
	 * @class SPHSystem
	 * @brief The SPHsystem managing objects in the system level.
	 */
	class SPHSystem
	{

	public:
		/**
		 * @brief Default constructor.
		 * @param[in] lower_bound Lower bound of the computational domain.
		 * @param[in] upper_bound Upper bound of the computational domain.
		 * @param[in] particle_spacing_ref Reference particle spacing.
		 * @param[in] smoothinglength_ratio The Referncen ratio of smoothing length to particle spacing.
		 */
 		SPHSystem(Vecd lower_bound, Vecd upper_bound,
			Real particle_spacing_ref, Real smoothinglength_ratio = 1.3);
 		/**
 		 * @brief Default descructor.
 		 */
		virtual ~SPHSystem();

		Vecd lower_bound_, upper_bound_;	/**< Lower and Upper domain bound. */
		Real particle_spacing_ref_;			/**< Refernce initial particle spacing. */
		Real smoothinglength_ratio_;		/**< Ratio between smoothing length to particle spacing. */
		int restart_step_;

		task_scheduler_init tbb_init_;		/**< TBB library. */

		StdVec<SPHBody*> bodies_;			/**< All sph bodies. */
		StdVec<SPHBody*> fictitious_bodies_;/**< The bodies without inner particle configuration. */
		StdVec<SPHBody*> real_bodies_;		/**< The bodies with inner particle configuration. */
		SPHBodyTopology* body_topology_;	/**< SPH body topology. */

		/**
		 * @brief Generate a kernel.
		 * @param[in] particle_spacing The reference initial particle spacing.
		 */
		Kernel* GenerateAKernel(Real particle_spacing);
		/**
		 * @brief Add a new body to the SPH system.
		 * @param[in] Pointer to a sph body.
		 */
		void AddBody(SPHBody* body);
		/**
		 * @brief Add a new body to the SPH real bodies.
		 * @param[in] body Pointer to a sph body.
		 */
		void AddRealBody(SPHBody* body);
		/**
		 * @brief Add a new body to the SPH fictitious bodies.
		 * @param[in] body Pointer to a sph body.
		 */
		void AddFictitiousBody(SPHBody* body);
		/**
		 * @brief Set up the body topology.
		 * @param[in] body_topology Topology of bobies.
		 */
		void SetBodyTopology(SPHBodyTopology* body_topology);
		/**
		 * @brief Create particles for the system, and particles are created and saved in the bodies.
		 */
		void CreateParticelsForAllBodies();
		/**
		 * @brief Set up the initial conditions for all bodies.
		 */
		void InitializeAllRealBodies();
		/**
		 * @brief Initialize cell linked lists system. 
		 */
		void InitializeSystemCellLinkedLists();
		/**
		 * @brief Initialize particle interacting configurations.
		 */
		void InitializeSystemConfigurations();
		/**
		 * @brief Set up, create particle, initializaiton, cell-linked list, configuration, to start a simulation.
		 */
		void SetupSPHSimulation();
		/**
		 * @brief Set up the initial conditions from restart files for all bodies.
		 */
		void ReInitializeAllRealBodiesFromRestart();
		/**
		 * @brief Set up, create particle, initializaiton, cell-linked list, configuration, to "Restart" a simulation.
		 */
		void ResetSPHSimulationFromRestart();
	};
}