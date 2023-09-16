/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file sph_system.h
 * @brief The SPH_System managing objects in the system level.
 * @details Note that the system operation prefer these are application independent.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>
#ifdef BOOST_AVAILABLE
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#endif

#include "base_data_package.h"
#include "sph_data_containers.h"

#include <filesystem>
#include <fstream>
#include <thread>
namespace fs = std::filesystem;

namespace SPH
{
/**
 * Pre-claimed classes.
 */
class SPHBody;
class IOEnvironment;
class ComplexShape;

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

    void setRunParticleRelaxation(bool run_particle_relaxation) { run_particle_relaxation_ = run_particle_relaxation; };
    bool RunParticleRelaxation() { return run_particle_relaxation_; };
    void setReloadParticles(bool reload_particles) { reload_particles_ = reload_particles; };
    bool ReloadParticles() { return reload_particles_; };
    void setRestartStep(size_t restart_step) { restart_step_ = restart_step; };
    size_t RestartStep() { return restart_step_; };
    BoundingBox system_domain_bounds_;       /**< Lower and Upper domain bounds. */
    Real resolution_ref_;                    /**< reference resolution of the SPH system */
    tbb::global_control tbb_global_control_; /**< global controlling on the total number parallel threads */

    IOEnvironment *io_environment_; /**< io_environment setup */
    bool generate_regression_data_; /**< run and generate or enhance the regression test data set. */

    SPHBodyVector sph_bodies_;         /**< All sph bodies. */
    SPHBodyVector observation_bodies_; /**< The bodies without inner particle configuration. */
    SPHBodyVector real_bodies_;        /**< The bodies with inner particle configuration. */
    SolidBodyVector solid_bodies_;     /**< The bodies with inner particle configuration and acoustic time steps . */
    /** Initialize cell linked list for the SPH system. */
    void initializeSystemCellLinkedLists();
    /** Initialize particle configuration for the SPH system. */
    void initializeSystemConfigurations();
    /** get the min time step from all bodies. */
    Real getSmallestTimeStepAmongSolidBodies(Real CFL = 0.6);
    /** Command line handle for Ctest. */
#ifdef BOOST_AVAILABLE
    void handleCommandlineOptions(int ac, char *av[]);
#endif
  protected:
    bool run_particle_relaxation_; /**< run particle relaxation for body fitted particle distribution */
    bool reload_particles_;        /**< start the simulation with relaxed particles. */
    size_t restart_step_;          /**< restart step */
};
} // namespace SPH
#endif // SPH_SYSTEM_H