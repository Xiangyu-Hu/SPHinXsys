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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	io_base_ck.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef IO_BASE_CK_H
#define IO_BASE_CK_H

#include "execution_policy.h"
#include "io_base.h"
#include "io_vtk.h"
#include "simple_algorithms_ck.h"

namespace SPH
{

template <class ExecutionPolicy>
class BodyStatesRecordingToVtpCK : public BodyStatesRecordingToVtp
{
  protected:
    OperationOnDataAssemble<DiscreteVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_write_;

    void prepareToWrite();

  public:
    template <typename... Args>
    BodyStatesRecordingToVtpCK(Args &&...args) : BodyStatesRecordingToVtp(std::forward<Args>(args)...){};
    virtual ~BodyStatesRecordingToVtpCK() {};
    virtual void writeToFile() override;
    virtual void writeToFile(size_t iteration_step) override;

    template <typename DerivedVariableMethod, typename DynamicsIdentifier, typename... Args>
    BodyStatesRecording &addDerivedVariableToWrite(DynamicsIdentifier &identifier, Args &&...args);
};

template <class ExecutionPolicy>
class RestartIOCK : public RestartIO
{
    UniquePtrsKeeper<AbstractDynamics> particle_dynamics_keeper_;

  public:
    template <typename... Args>
    RestartIOCK(Args &&...args);
    virtual ~RestartIOCK() {};
    virtual void writeToFile(size_t iteration_step) override;
    virtual void readFromFile(size_t iteration_step) override;
    void reportEvolvingVariablesBounds(size_t restart_step);

  protected:
    StdVec<StdVec<BaseDynamics<std::pair<Real, UnsignedInt>> *>> output_evolving_variables_bounds_[3];
    StdVec<StdVec<std::string>> evolving_variables_names_[3];
    OperationOnDataAssemble<DiscreteVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_write_;
    OperationOnDataAssemble<DiscreteVariables, FinalizeVariablesAfterRead<DiscreteVariable>> finalize_variables_after_read_;
};

template <class ExecutionPolicy>
class ReloadParticleIOCK : public ReloadParticleIO
{
  public:
    template <typename... Args>
    ReloadParticleIOCK(Args &&...args) : ReloadParticleIO(std::forward<Args>(args)...){};
    virtual ~ReloadParticleIOCK() {};
    virtual void writeToFile(size_t iteration_step = 0) override;

  protected:
    OperationOnDataAssemble<DiscreteVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_reload_;
};
} // namespace SPH
#endif // IO_BASE_CK_H
