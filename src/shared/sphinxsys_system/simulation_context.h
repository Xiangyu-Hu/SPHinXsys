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
 * @file simulation_context.h
 * @brief Handling global simulation context and variables.
 * @details TBD
 * @author	Xiangyu Hu
 */

#ifndef SIMULATION_CONTEXT_H
#define SIMULATION_CONTEXT_H

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>
#ifdef BOOST_AVAILABLE
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#endif

#include "base_data_package.h"
#include "execution_policy.h"
#include "io_environment.h"

namespace SPH
{
using namespace execution;

class SimulationContext
{
    UniquePtrKeeper<IOEnvironment> io_ptr_keeper_;
    DataContainerUniquePtrAssemble<SingularVariable> all_context_variable_ptrs_;
    UniquePtrsKeeper<Entity> unique_context_variable_ptrs_;

  public:
    SimulationContext(BoundingBox system_domain_bounds, Real resolution_ref,
                      size_t number_of_threads = std::thread::hardware_concurrency());
    virtual ~SimulationContext() {};

#ifdef BOOST_AVAILABLE
    SimulationContext *handleCommandlineOptions(int ac, char *av[]);
#endif

    SimulationContext *setIOEnvironment(bool delete_output = true);
    IOEnvironment &getIOEnvironment();
    Real ReferenceResolution() { return resolution_ref_; };
    BoundingBox SystemDomainBounds() { return system_domain_bounds_; };
    void setRunParticleRelaxation(bool run_particle_relaxation) { run_particle_relaxation_ = run_particle_relaxation; };
    bool RunParticleRelaxation() { return run_particle_relaxation_; };
    void setReloadParticles(bool reload_particles) { reload_particles_ = reload_particles; };
    bool ReloadParticles() { return reload_particles_; };
    bool GenerateRegressionData() { return generate_regression_data_; };
    void setGenerateRegressionData(bool generate_regression_data) { generate_regression_data_ = generate_regression_data; };
    bool StateRecording() { return state_recording_; };
    void setStateRecording(bool state_recording) { state_recording_ = state_recording; };
    void setRestartStep(size_t restart_step) { restart_step_ = restart_step; };
    size_t RestartStep() { return restart_step_; };

    template <typename DataType>
    SingularVariable<DataType> *registerContextVariable(
        const std::string &name, DataType initial_value = ZeroData<DataType>::value)
    {
        SingularVariable<DataType> *variable =
            findVariableByName<DataType>(all_system_variables_, name);

        return variable != nullptr
                   ? variable
                   : addVariableToAssemble<DataType>(
                         all_system_variables_, all_system_variable_ptrs_, name, initial_value);
    };

    template <typename DataType>
    SingularVariable<DataType> *getContextVariableByName(const std::string &name)
    {
        SingularVariable<DataType> *variable =
            findVariableByName<DataType>(all_system_variables_, name);

        if (variable == nullptr)
        {
            std::cout << "\nError: the system variable '" << name << "' is not registered!\n";
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        }

        return variable;
    };

    template <typename DataType>
    DataType *getContextVariableDataByName(const std::string &name)
    {
        SingularVariable<DataType> *variable =
            findVariableByName<DataType>(all_system_variables_, name);

        if (variable == nullptr)
        {
            std::cout << "\nError: the system variable '" << name << "' is not registered!\n";
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        }

        return variable->Data();
    };

  protected:
    BoundingBox system_domain_bounds_;       /**< Lower and Upper domain bounds. */
    Real resolution_ref_;                    /**< reference resolution of the SPH system */
    tbb::global_control tbb_global_control_; /**< global controlling on the total number parallel threads */
    bool run_particle_relaxation_;           /**< run particle relaxation for body fitted particle distribution */
    bool reload_particles_;                  /**< start the simulation with relaxed particles. */
    size_t restart_step_;                    /**< restart step */
    bool generate_regression_data_;          /**< run and generate or enhance the regression test data set. */
    bool state_recording_;                   /**< Record state in output folder. */
};
} // namespace SPH
#endif // SIMULATION_CONTEXT_H
