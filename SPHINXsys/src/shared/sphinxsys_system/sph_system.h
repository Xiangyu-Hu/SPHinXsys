/**
 * @file sph_system.h
 * @brief The SPH_System managing objects in the system level.
 * @details Note that the system operation prefer these are application independent.
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */

#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H


#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>
#ifdef BOOST_AVAILABLE
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#endif

#include "base_data_package.h"
#include "sph_data_conainers.h"

#include <thread>
#include <fstream>
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH 
{
	/**
	 * @brief Preclaimed classes.
	 */
	class SPHBody;

	/**
	 * @class SPHSystem
	 * @brief The SPH system managing objects in the system level.
	 */
	class SPHSystem
	{
	public:
 		SPHSystem(BoundingBox system_domain_bounds, Real resolution_ref,
			size_t number_of_threads = std::thread::hardware_concurrency());
		virtual ~SPHSystem() {};

		BoundingBox system_domain_bounds_;			/**< Lower and Upper domain bounds. */
		Real resolution_ref_;						/**< refernce resolution of the SPH system */
		tbb::global_control tbb_global_control_;	/**< global controling on the total number parallel threads */
	
		size_t restart_step_;				/**< restart step */
		bool run_particle_relaxation_;		/**< run particle relaxation for body fitted particle distribution */
		bool reload_particles_;				/**< start the simulation with relaxed particles. */
		std::string output_folder_;			/**< Folder for saving output files. */
		std::string restart_folder_;		/**< Folder for saving restart files. */
		std::string reload_folder_;			/**< Folder for saving particle reload files. */

		SPHBodyVector bodies_;				/**< All sph bodies. */
		SPHBodyVector fictitious_bodies_;	/**< The bodies without inner particle configuration. */
		SPHBodyVector real_bodies_;			/**< The bodies with inner particle configuration. */

		void addABody(SPHBody* sph_body);
		void addARealBody(RealBody* Real_body);
		void addAFictitiousBody(FictitiousBody* fictitious_body);
		void initializeSystemCellLinkedLists();
		void initializeSystemConfigurations();
		#ifdef BOOST_AVAILABLE
		void handleCommandlineOptions(int ac, char* av[]);
		#endif
	};
}
#endif //SPH_SYSTEM_H