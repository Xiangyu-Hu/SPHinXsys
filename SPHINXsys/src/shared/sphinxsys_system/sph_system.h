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
#include "sph_data_containers.h"

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
	class In_Output;

	/**
	 * @class SPHSystem
	 * @brief The SPH system managing objects in the system level.
	 */
	class SPHSystem
	{
	public:
		SPHSystem(BoundingBox system_domain_bounds, Real resolution_ref,
				  size_t number_of_threads = std::thread::hardware_concurrency());
		virtual ~SPHSystem(){};

		BoundingBox system_domain_bounds_;		 /**< Lower and Upper domain bounds. */
		Real resolution_ref_;					 /**< reference resolution of the SPH system */
		tbb::global_control tbb_global_control_; /**< global controlling on the total number parallel threads */

		In_Output *in_output_;			/**< in_output setup */
		size_t restart_step_;			/**< restart step */
		bool run_particle_relaxation_;	/**< run particle relaxation for body fitted particle distribution */
		bool reload_particles_;			/**< start the simulation with relaxed particles. */
		bool generate_regression_data_; /**< run and generate or enhancethe regression test data set. */

		SPHBodyVector bodies_;			  /**< All sph bodies. */
		SPHBodyVector fictitious_bodies_; /**< The bodies without inner particle configuration. */
		SPHBodyVector real_bodies_;		  /**< The bodies with inner particle configuration. */
		SolidBodyVector solid_bodies_;	  /**< The bodies with inner particle configuration and acoustic time steps . */

		void addABody(SPHBody *sph_body);
		void addARealBody(RealBody *real_body);
		void addASolidBody(SolidBody *solid_body);
		void addAFictitiousBody(FictitiousBody *fictitious_body);
		void initializeSystemCellLinkedLists();
		void initializeSystemConfigurations();
		Real getSmallestTimeStepAmongSolidBodies(Real CFL = 0.6);
#ifdef BOOST_AVAILABLE
		void handleCommandlineOptions(int ac, char *av[]);
#endif
	};
}
#endif //SPH_SYSTEM_H